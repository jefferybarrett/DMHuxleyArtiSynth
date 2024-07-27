package artisynth.models.pendulum;


import java.awt.Color;
import java.io.IOException;

import artisynth.core.workspace.RootModel;
import artisynth.core.materials.*;
import artisynth.core.mechmodels.*;
import artisynth.core.modelbase.*; 
import artisynth.core.mechmodels.MechSystemSolver.Integrator;
import maspack.geometry.PolygonalMesh;
import maspack.matrix.AxisAlignedRotation;
import maspack.matrix.*;
import maspack.render.RenderProps;
import maspack.util.PathFinder;
import maspack.spatialmotion.*; 
import maspack.util.*;
import maspack.interpolation.*;
import artisynth.core.probes.*; 
import maspack.interpolation.*;



public class InvertedPendulum extends RootModel {


    String geodir = PathFinder.getSourceRelativePath (this, "geometry/");

    // some other physical constants
    double g = 9.81 * (Units.m / (Units.s*Units.s));

    // include some geometry here
    double a = 0.05 * Units.m;
    double b = 14.2 * Units.cm;//42.6 * Units.cm; //(0.3052) * Units.m;

    // commented here is from the old geometry of the inverted pendulum
    double L = b; //Math.sqrt(a*a + b*b);    // length in upright position
    double r = a; // a*b / L;                 // moment arm in upright position

    // cervical spine properties
    double pendulum_length = 0.3052 * Units.m; // this is set in stone by the geometry
    double pendulum_mass = 30.0 * Units.kg; // kg (5 kg for cervical)

    // this will be whether or not we use springs, or muscles
    boolean doSprings = false;
    boolean doDampedOscillator = false;
    boolean doSLS = true;
    boolean doHuxley = false;

    // spring properties
    // >> xstar is the resting length
    // according to my original theory, the stiffness needs to satisfy:
    // r2 k1 + r2 k2 > mgl
    // so then, assuming symmetry:
    // >> k > mgl/(2 r^2)
    // new theory (this is actually the correct theory... annoyingly):
    // >> k > mgl/(2 r^2) x [L/L0]
    // with the new symmetry (more like Bergmark) the resting length's cancel out
    // (I don't know why they cancel out, but they do.)
    // >> k > mgl / (2 a^2) [critical stiffness]
    double stiffness_multiplier = 0.8; // like a factor of safety.
    double critical_stiffness = (pendulum_mass * g * pendulum_length) / (2 * r * r);
    double L0 = L;
    //double common_stiffness = stiffness_multiplier * (L/L0) * (pendulum_mass * g * pendulum_length) / (2 * r * r);
    double common_stiffness = stiffness_multiplier * critical_stiffness;
    //double common_stiffness = 2565711.26 * (Units.N / Units.m);// 1214436.66; // (max stiffness for each muscle)
    double flexor_xstar = L0;
    double extensor_xstar = L0;
    double flexor_k = common_stiffness;// * (Units.N / Units.m);
    double extensor_k = common_stiffness;// * (Units.N / Units.m);

    double SRS_stiffness_multiplier = 1000.0; // make this 1.0 for relatively nice simulations. 10x for high short-range stiffness
    // Here are the properties for the SLS
    double SLS_common_stiffness = common_stiffness;
    double SLS_short_range_stiffness = SRS_stiffness_multiplier * (10.0/9.0) * SLS_common_stiffness;

    // default muscle properties
    double relaxation_time = 10.0 * Units.ms;
    double muscle_fibre_length = 14.2 * Units.cm; //9.6 * Units.cm;
    //double muscle_damping = relaxation_time * SLS_short_range_stiffness;//1000.0 * (Units.N * Units.s / Units.m);
    double muscle_damping = 10000.0 * (Units.N * Units.s / Units.m);
    double muscle_pcsa = 24.0 * (Units.cm2);
    double[] activations = {0.005, 0.01, 0.10, 0.5, 1.0};
    int act_index = 4;


    // perturbation properties
    double perturbation_magnitude = -15.0 * Units.Nm; // Nm
    double[] perturbation_magnitudes = {-0.5 * Units.Nm, -5.0 * Units.Nm, -50.0 * Units.Nm};
    int pert_index = 1;

    // collection parameters
    double Trange = 60.0 * Units.s;

    public void build(String[] args) throws IOException {
        MechModel mech = new MechModel("mech");
        addModel(mech);
        mech.setGravity(0.0, -g, 0.0); // y-up

        // stuff for saving the file
        Interpolation interp = new Interpolation();
        interp.setOrder(Interpolation.Order.Cubic);
        String angle_filename = "TEST_ANGLE.txt";
        String force_filename = "TEST_FORCE.txt";
        String flexor_filename = "TEST_FLEXOR.txt";
        String extensor_filename = "TEST_EXTENSOR.txt";
        
        // adding in the rigid bodies
        PolygonalMesh floor_mesh = new PolygonalMesh(geodir + "floor.obj");
        floor_mesh.triangulate();
        RigidBody floor = new RigidBody("floor");
        floor.setSurfaceMesh(floor_mesh);
        floor.setDynamic(false);
        mech.addRigidBody(floor);


        PolygonalMesh pendulum_mesh = new PolygonalMesh(geodir + "pendulum2.obj");
        pendulum_mesh.triangulate();
        RigidBody pendulum = new RigidBody("pendulum");
        SpatialInertia pendulum_inertia = new SpatialInertia(0.0, 0.0, 0.0, 0.0);
        pendulum_inertia.addPointMass(pendulum_mass, new Vector3d(0.0, pendulum_length, 0.0));
        pendulum.setSurfaceMesh(pendulum_mesh);
        pendulum.setInertia(pendulum_inertia);
        mech.addRigidBody(pendulum);

        // now we just need to add the hinge-joint at O
        Point3d ORIGIN = new Point3d(0.0, 0.0, 0.0);
        HingeJoint hinge = new HingeJoint(floor, pendulum, ORIGIN, new Vector3d(0.0, 0.0, 1.0));
        hinge.setThetaRange(-45.0, 45.0);
        mech.addBodyConnector(hinge);

        // muscle geometry
        FrameMarker MASS_POSITION = new FrameMarker(pendulum, new Point3d(0.0, pendulum_length, 0.0));
        FrameMarker M1_ORIGIN = new FrameMarker(floor, new Point3d(a, 0.0, 0.0));
        FrameMarker M2_ORIGIN = new FrameMarker(floor, new Point3d(-a, 0.0, 0.0));
        FrameMarker M1_INSERTION = new FrameMarker(pendulum, new Point3d(a, b, 0.0)); 
        FrameMarker M2_INSERTION = new FrameMarker(pendulum, new Point3d(-a, b, 0.0));
        mech.add(M1_ORIGIN);
        mech.add(M2_ORIGIN);
        mech.add(M1_INSERTION);
        mech.add(M2_INSERTION);
        mech.add(MASS_POSITION);

        if (doSprings){
            // here we set up the logic to attach springs to the model
            AxialSpring flexor = new AxialSpring("Muscle1", flexor_xstar);
            AxialSpring extensor = new AxialSpring("Muscle2", extensor_xstar);
            flexor.setPoints(M1_ORIGIN, M1_INSERTION);
            extensor.setPoints(M2_ORIGIN, M2_INSERTION);

            flexor.setMaterial(
                new LinearAxialMaterial(flexor_k, 0.0)
            );

            extensor.setMaterial(
                new LinearAxialMaterial(extensor_k, 0.0)
            );

            mech.addAxialSpring(flexor);
            mech.addAxialSpring(extensor);
            RenderProps.setSpindleLines(flexor, 0.02, Color.BLUE);
            RenderProps.setSpindleLines(extensor, 0.02, Color.BLUE);

            // saving
            angle_filename = "angle_springs.csv";
            force_filename = "force_springs.csv";
            flexor_filename = "flexor_springs.csv";
            extensor_filename = "extensor_springs.csv";

            String[] muscleProperties = {"length"};
            NumericOutputProbe flexorProbe = new NumericOutputProbe(
                flexor,
                muscleProperties,
                PathFinder.getSourceRelativePath(this, "/outputs/" + flexor_filename),
                0.001
            );
            flexorProbe.setFormat(",%g");
            flexorProbe.setShowHeader(false);
            flexorProbe.setInterpolation(interp);
            flexorProbe.setName("Flexor");
            flexorProbe.setStopTime(Trange);
            addOutputProbe(flexorProbe);

            NumericOutputProbe extensorProbe = new NumericOutputProbe(
                extensor,
                muscleProperties,
                PathFinder.getSourceRelativePath(this, "/outputs/" + extensor_filename),
                0.001
            );
            extensorProbe.setFormat(",%g");
            extensorProbe.setShowHeader(false);
            extensorProbe.setInterpolation(interp);
            extensorProbe.setName("Extensor");
            extensorProbe.setStopTime(Trange);
            addOutputProbe(extensorProbe);
        }
        else if (doDampedOscillator){
            // here we set up the logic to attach springs to the model
            AxialSpring flexor = new AxialSpring("Muscle1", flexor_xstar);
            AxialSpring extensor = new AxialSpring("Muscle2", extensor_xstar);
            flexor.setPoints(M1_ORIGIN, M1_INSERTION);
            extensor.setPoints(M2_ORIGIN, M2_INSERTION);

            flexor.setMaterial(
                new LinearAxialMaterial(flexor_k, muscle_damping)
            );

            extensor.setMaterial(
                new LinearAxialMaterial(extensor_k, muscle_damping)
            );

            mech.addAxialSpring(flexor);
            mech.addAxialSpring(extensor);
            RenderProps.setSpindleLines(flexor, 0.02, Color.BLUE);
            RenderProps.setSpindleLines(extensor, 0.02, Color.BLUE);

            // saving
            angle_filename = "angle_damped.csv";
            force_filename = "force_damped.csv";
            flexor_filename = "flexor_damped.csv";
            extensor_filename = "extensor_damped.csv";

        }
        else if (doSLS){
            // set up the stuff to have the Maxwell Elements here
            // stuff to try:
            // >>> see about making these AxialSprings themselves implement HasNumericState
            AxialSpring flexor = new AxialSpring("Muscle1", flexor_xstar);
            AxialSpring extensor = new AxialSpring("Muscle2", extensor_xstar);
            flexor.setPoints(M1_ORIGIN, M1_INSERTION);
            extensor.setPoints(M2_ORIGIN, M2_INSERTION);

            System.out.println("Common Stiffness = " + SLS_common_stiffness / (Units.N/Units.cm) + " N/cm");

            flexor.setMaterial(
                new StandardLinearSolidAxialMaterial(SLS_common_stiffness, SLS_short_range_stiffness, muscle_damping)
            );
            extensor.setMaterial(
                new StandardLinearSolidAxialMaterial(SLS_common_stiffness, SLS_short_range_stiffness, muscle_damping)
            );
            
            mech.add(flexor);
            mech.add(extensor);
            RenderProps.setSpindleLines(flexor, 0.02, Color.GREEN);
            RenderProps.setSpindleLines(extensor, 0.02, Color.GREEN);

            angle_filename = "angle_SLS.csv";
            force_filename = "force_SLS.csv";
        }
        else if (doHuxley){
            // set up the logic to do the Huxley muscle model
            Muscle flexor = new Muscle("Muscle1", flexor_xstar);
            Muscle extensor = new Muscle("Muscle2", extensor_xstar);
            flexor.setPoints(M1_ORIGIN, M1_INSERTION);
            extensor.setPoints(M2_ORIGIN, M2_INSERTION);

            UniformDMMuscle flexorMat = new UniformDMMuscle(muscle_fibre_length, b - muscle_fibre_length, muscle_pcsa);
            UniformDMMuscle extensorMat = new UniformDMMuscle(muscle_fibre_length, b - muscle_fibre_length, muscle_pcsa);
            double critical_activation = critical_stiffness / flexorMat.maxIsoStiffness;
            double activation = activations[act_index];
            flexorMat.isometricEquilibriate(muscle_fibre_length, 0.0, activation);
            extensorMat.isometricEquilibriate(muscle_fibre_length, 0.0, activation);

            System.out.println("Crtical stiffness = " + (critical_stiffness) + " N/m");
            System.out.println("Critical Activation = " + (critical_activation*100) + "%");
            System.out.println("Activation set to " + 100*activation + "%");
            System.out.println("Muscle Stiffness = " + flexorMat.computeDFdl(b, 0.0, 0.0, activation) + " N/m");
            
            flexor.setMaterial(flexorMat);
            extensor.setMaterial(extensorMat);

            //UniformDMMuscle matTest = new UniformDMMuscle(b, 0.0, 16.0 * (Units.cm2));
            //double activation = common_stiffness / matTest.maxIsoStiffness;
            //matTest.isometricEquilibriate(flexor.getLength(), 0.0, activation);
            // now let's calculate the force in muscle 1 at 0% activation
            //System.out.println("Force = " + matTest.computeF(flexor.getLength(), 0.0, 0.0, 0.0));
            
            flexor.setExcitation(activation);
            extensor.setExcitation(activation);
            //flexor.getMaterial().isometricEquilibriate(b, 0.0, activation);
            //extensor.getMaterial().isometricEquilibriate(b, 0.0, activation);


            mech.add(flexor);
            mech.add(extensor);
            RenderProps.setSpindleLines(flexor, 0.02, Color.RED);
            RenderProps.setSpindleLines(extensor, 0.02, Color.RED);

            angle_filename = "angle_DM_" + act_index + ".csv";
            force_filename = "force_DM_" + act_index + ".csv";
            flexor_filename = "flexor_DM_" + act_index + ".csv";
            extensor_filename = "extensor_DM_" + act_index + ".csv";

            String[] muscleProperties = {"forceNorm", "length"};
            NumericOutputProbe flexorProbe = new NumericOutputProbe(
                flexor,
                muscleProperties,
                PathFinder.getSourceRelativePath(this, "/outputs/" + flexor_filename),
                0.001
            );
            flexorProbe.setFormat(",%g");
            flexorProbe.setShowHeader(false);
            flexorProbe.setInterpolation(interp);
            flexorProbe.setName("Flexor");
            flexorProbe.setStopTime(Trange);
            addOutputProbe(flexorProbe);

            NumericOutputProbe extensorProbe = new NumericOutputProbe(
                extensor,
                muscleProperties,
                PathFinder.getSourceRelativePath(this, "/outputs/" + extensor_filename),
                0.001
            );
            extensorProbe.setFormat(",%g");
            extensorProbe.setShowHeader(false);
            extensorProbe.setInterpolation(interp);
            extensorProbe.setName("Extensor");
            extensorProbe.setStopTime(Trange);
            addOutputProbe(extensorProbe);

        } else {
            // LOL
            System.out.println("Not supported by anything!");

            angle_filename = "angle_unsupported.csv";
            force_filename = "force_unsupported.csv";
        }

        //hinge.setTheta(0.1);
        //PointForce perturbForce = new PointForce(new Vector3d(0.0, 0.0, 0.0), MASS_POSITION);
        //PerturbationController controller = new PerturbationController(perturbForce);
        PerturbationController controller = new PerturbationController(pendulum);
        controller.momentZ = new CubicHermiteSpline1d(
            controller.CreateKnotPointsForPerturbation(
                1.0 * Units.s,                          // start time
                50.0 * Units.ms,                        // duration
                0.0 * Units.Nm,                         // starting force
                perturbation_magnitudes[pert_index]     // perturbation force
            )
        );

        //pendulum.addExternalForce(new Wrench(1.0, 0.0, 0.0, 0.0, 0.0, 0.0));

        //mech.add(perturbForce);
        addController(controller);
        //perturbForce.setForceLengthRatio(0.005);
        //perturbForce.setAxisRadiusRatio(0.0500);
        //perturbForce.setRenderOutward(true);
        //RenderProps.setLineColor(perturbForce, Color.BLUE);
        mech.setIntegrator(Integrator.Trapezoidal);
        mech.setMaxStepSize(0.0001);
        // make sure we're dealing with y-up
        setDefaultViewOrientation(AxisAlignedRotation.X_Y);


        // ----------------------------------------
        // STUFF FOR SAVING DATA
        // ----------------------------------------
        NumericOutputProbe angleProbe = new NumericOutputProbe(
            hinge, 
            "theta", 
            PathFinder.getSourceRelativePath(this, "/outputs/" + angle_filename),
            0.0001
        );
        angleProbe.setFormat(",%g");
        angleProbe.setName("Angle");
        angleProbe.setShowHeader(false);
        angleProbe.setStopTime(Trange);
        angleProbe.setInterpolation(interp);
        addOutputProbe(angleProbe);

        NumericOutputProbe forceProbe = new NumericOutputProbe(
           pendulum,
            "externalForce",
            PathFinder.getSourceRelativePath(this, "/outputs/" + force_filename),
            0.0001
        );
        forceProbe.setFormat(",%g");
        forceProbe.setShowHeader(false);
        forceProbe.setInterpolation(interp);
        forceProbe.setName("Force");
        forceProbe.setStopTime(Trange);
        addOutputProbe(forceProbe);
    }


    
    private class PerturbationController extends ControllerBase {

        // this is the direction of the perturbation, as a wrench
        RigidBody myBody;
        Wrench myWrench = new Wrench();

        public double[] defaultZero = {
            0.0 * Units.ms, 0.0, 0.0,
            0.0 * Units.ms, 0.0, 0.0,
            10.0 * Units.ms, 0.0, 0.0,
        };

        public CubicHermiteSpline1d forceX = new CubicHermiteSpline1d(defaultZero);
        public CubicHermiteSpline1d forceY = new CubicHermiteSpline1d(defaultZero);
        public CubicHermiteSpline1d forceZ = new CubicHermiteSpline1d(defaultZero);
        public CubicHermiteSpline1d momentX = new CubicHermiteSpline1d(defaultZero);
        public CubicHermiteSpline1d momentY = new CubicHermiteSpline1d(defaultZero);
        public CubicHermiteSpline1d momentZ = new CubicHermiteSpline1d(defaultZero);


        public PerturbationController(RigidBody body){
            myBody = body;
            updateWrench(0.0);
        }

        @Override
        public void apply(double t0, double t1){
            //Vector3d newForce = new Vector3d(forceX.evalY(t0), forceY.evalY(t0), forceZ.evalY(t0));
            //myForce.setForce(newForce);
            updateWrench(t1);
            myBody.setExternalForce(myWrench);
        }

        public void updateWrench(double t){
            myWrench.set(
                forceX.evalY(t),
                forceY.evalY(t),
                forceZ.evalY(t),
                momentX.evalY(t),
                momentY.evalY(t),
                momentZ.evalY(t)
            );
        }

        public double[] CreateKnotPointsForPerturbation(double start_time, 
        double duration, double startForce, double perturbationForce){
            double ramp_time = 1.0 * Units.ms;

            double[] points = {
                0.0 * Units.s, startForce * Units.N, 0.0,
                (start_time - ramp_time) * Units.s, (startForce) * Units.N, 0.0,
                (start_time + ramp_time) * Units.s, (startForce + perturbationForce) * Units.N, 0.0,
                (start_time + duration - ramp_time) * Units.s - ramp_time, (startForce + perturbationForce) * Units.N, 0.0,
                (start_time + duration + ramp_time) * Units.s + ramp_time, (startForce) * Units.N, 0.0,
                (start_time + duration + 2*ramp_time) * Units.s, (startForce) * Units.N, 0.0,
            };
            return points;
        }
    }


}



