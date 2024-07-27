package artisynth.models.pendulum;

import artisynth.core.modelbase.HasNumericState;
import artisynth.core.materials.*;
import maspack.util.DataBuffer;
import maspack.properties.PropertyList;


public class UniformDMMuscle 
    extends AxialMuscleMaterialBase implements HasNumericState {

    /*  This is an implementation of a Distribution Moment Approximation
     *  with a rigid tendon. The arrangement of this model is:
     * 
     *      |
     *      |  tendon (rigid)              X    
     *      |----/\/\/\/\------------<==========>----> Force
     *      |       Y             Contractile Element          
     *      |                 (Distribution Moment Approx)
     * 
     * The total length of the element is L = X + Y. In a refernce position
     * this is L0 = X0 + Y0. Because the tendon is assumed to be rigid, Y = Y0
     * Regardless of the length.
     * 
     * The contractile element is a simplification of Huxley's Equation
     *          dn/dt - v dn/dx = f(x) - (f(x) + g(x))n(x,t)        (1)
     * 
     * We renormalize to xi = x/h, where h is the "characteristic bond length".
     * This gives:
     *          dn/dt - u dn/dxi = f(xi) - (f(xi) + g(xi))n(xi,t)
     * Where u = v/h. We further add activation and force-length characteristics:
     *          dn/dt - u dn/dxi = a(t) F(X) f(xi) - (f(xi) + g(xi))n(xi,t)
     * 
     * Next, we recast this to involve the distribution moment approximation.
     * 
     * 
     * By default, this uses the uniform distribution since it is computationally
     * the simplest.
     * 
     * By: Jeff M. Barrett
    */

    // macroscopic property default values
    protected static double DEFAULT_CONTRACTILE_LENGTH = 0.20 * Units.m;                     // X_0 in Zahalak
    protected static double DEFAULT_CROSS_SECTIONAL_AREA = 1500 * Units.mm2; 
    protected static double DEFAULT_DAMPING = 40.0 * Units.N/(Units.m / Units.s);
    protected static double DEFAULT_TENDON_LENGTH = 0.0 * Units.m;                          // Y_0 in Zahalak

    // setting these up as normal properties
    protected double X_0 = DEFAULT_CONTRACTILE_LENGTH;
    protected double Y_0 = DEFAULT_TENDON_LENGTH;
    protected double L_0 = X_0 + Y_0;
    protected double A_0 = DEFAULT_CROSS_SECTIONAL_AREA;
    protected double damping  = DEFAULT_DAMPING;

    // binding rates
    protected double f1 = 15.0 * Units.Hz; // Binding rate between 0 and h

    // detachment rates
    protected double g0 = 170.0 * Units.Hz; // Dissociation rate constant for compression
    protected double g1 = 8.0 * Units.Hz;   // Dissociation rate for short displacements
    protected double g2 = 25.0 * Units.Hz;  // Dissociation rate for large displacements

    // these are currently not editable, but maybe they'll be available later on
    // these parameters give a theoretical maximum stress of 34.5 N/cm2, which I think is reasonable
    protected double s_sarcomere = 2.2 * Units.um;          // [sarcomere length] (s_0 in Zahalak)
    protected double h_myosin = 5.0 * Units.nm;             // [characteristic bond length] (h in Zahalak and Huxley)
    protected double ell_actin = 36.0 * Units.nm;           // [distance between attachment sites] (this should be larger than h)
    protected double k_myosin = 5.0 * Units.pN / Units.nm;  // [myosin cross-bridge stiffness] 5 pN/nm (0.005 N/m)
    // this could be as low as 1-1.2 pN/nm for skeletal muscle (Pinzauti et al., 2018; Wang et al., 2020)
    // Huxley and Simmons had this at 2.5E-4 N/m (0.25 pN/nm)
    // the slow-twitch stuff is around 0.2, fast twitch is around 5. Which is interesting.
    protected double m_myosin = 2.3*2.0E-4 * Units.mol/Units.L; // [myosin concentration] mols/L
    // this concentration of myosin (which I think is slightly too high) makes up
    // the bulk of the correction

    // derived variables
    protected double Gamma0 = (A_0 * m_myosin * s_sarcomere * s_sarcomere * k_myosin * h_myosin) / (4 * ell_actin * X_0);
    protected double Gamma1 = (A_0 * m_myosin * s_sarcomere * k_myosin * h_myosin * h_myosin) / (2 * ell_actin);
    protected double gamma = (2.0 * X_0 * h_myosin) / (s_sarcomere); // I think Cholewicki has an extra L0 here. It does not seem necessary (for rigid tendons)
    protected double varphi = f1 / (f1 + g1);
    public double maxIsoForce = Gamma1 * (varphi/2);
    public double maxIsoStress = maxIsoForce / A_0;
    public double maxIsoStiffness = Gamma0 * varphi;

    // derived from state-variables
    protected double N = 0.0;
    protected double xi0 = 0.5;
    protected double Delta = Math.sqrt(0.5);

    
    // state variables
    double Q0 = varphi;
    double Q1 = varphi/2;
    double Q2 = varphi/3;
    protected double activation = 0.0;  // hack for storing activaiton
    protected double u_normalized = 0.0;         // hack for storing the velocity
    

    // numeric parameters
    private double Q0_MIN = 1E-6; // minimum Q0

    // force-length relationship
    //public CubicHermiteSpline1d forceLengthRelationship = new CubicHermiteSpline1d();

    /* =======================================
     * Property List, getters and setters
     * =======================================
    */
    // default values
    public static double sarcomere_length_default = 2.2 * Units.um;
    
    public static PropertyList myProps = new PropertyList(UniformDMMuscle.class, AxialMuscleMaterialBase.class);
    
    public PropertyList getAllPropertyInfor(){
        return myProps;
    }

    static {
        myProps.add("s0", "Sarcomere length at optimum", sarcomere_length_default);
    }

    public double getS0(){
        return s_sarcomere;
    }

    public void setS0(double _s0){
        s_sarcomere = _s0;
    }


    /* =======================================
     * Creation Methods
     * =======================================
    */

    public UniformDMMuscle(){
        // uses the default values (not recommended?)
        updateDerivedGeometricValues();
        updateDerivedStateValues();
    }

    public UniformDMMuscle(
        double fibreLength, 
        double tendonLength,
        double physiologicalCrossSectionalArea){
        
        // set geometric properties
        // note: I will need to do this with setters though
        X_0 = fibreLength;
        Y_0 = tendonLength;
        A_0 = physiologicalCrossSectionalArea;

        updateDerivedGeometricValues();
        updateDerivedStateValues();
    }

    // updateDerivedValues : void -> void
    // Purpose: updates the derived values
    private void updateDerivedStateValues(){
        if (Q0 > Q0_MIN){
            N = N_approx(Q0, Q1, Q2);
            xi0 = xi0_approx(Q0, Q1, Q2);
            Delta = Delta_approx(Q0, Q1, Q2);
        } else{
            System.out.println("Warning: Minimum Q0 not achieved!");
            N = 0.0;
            xi0 = 1/2.0;
            Delta = 1/Math.sqrt(2.0);
        }
    }

    // setQstate : double double double double -> void
    // Purpose: Explicitly set Q0, Q1, Q2 and ldot to the specified magnitudes
    //          This is mostly for debugging purposes.
    public void setQstate(double _Q0, double _Q1, double _Q2, double ldot){
        Q0 = _Q0;
        Q1 = _Q1;
        Q2 = _Q2;
        u_normalized = -ldot/gamma;
        updateDerivedStateValues();
    }
    
    // updateDerivedGeometricValues : void -> void
    // Purpose: Updates variables derived from the material properties and geometry of the muscle
    //      Gamma0 : constant of proportionality with Q0 and stiffness (i.e. K = Gamma0 * Q0)
    //      Gamma1 : constant of proportionality with Q1 and force (i.e. P = Gamma1 * Q1)
    //      gamma  : constant of proportionality between macroscopic velocity and microscopic.
    //                          (i.e. Xdot = gamma * u)
    //      varphi : steady-state ratio of binding to unbinding rates.
    private void updateDerivedGeometricValues(){
        Gamma0 = (A_0 * m_myosin * s_sarcomere * s_sarcomere * k_myosin * h_myosin) / (4 * ell_actin * X_0);
        Gamma1 = (A_0 * m_myosin * s_sarcomere * k_myosin * h_myosin * h_myosin) / (2 * ell_actin);
        gamma = (2.0 * X_0 * h_myosin) / (s_sarcomere); // extra factor of L0 here?
        varphi = f1 / (f1 + g1);
        maxIsoForce = Gamma1 * (varphi/2);
        maxIsoStiffness = Gamma0 * varphi;
        maxIsoStress = maxIsoForce / A_0;
    }

    /* =======================================
     * Methods that implement the DM Model
     * =======================================
     * 
     * dQ_lambda/dt + lambda u Q_(lambda-1) = Beta_lambda - (Phi_lambda + Psi_lambda)
     * 
     * Where:
     *      Beta_lambda : the partial moment associated with f(xi)
     *      Phi_lambda  : the partial moment associated with f(xi)n(xi,t)
     *      Psi_lambda  : the partial moment associated with g(xi)n(xi,t)
     * 
    */

    // Beta : int -> double
    // Purpose: Evaluates the integral:
    //              integral xi^lambda f(xi) dxi from -oo to oo
    //          For the usual rate function, we get f1/(lambda + 2) for this.
    public double Beta(int lambda){
        return f1/(lambda + 2);
    }

    // J_lambda : int, double, double -> double
    // Purpose: Calculates the J_lambda terms. These are the partial moments
    //          of the distribution.
    // Inputs: 
    // Outputs:
    public double J_lambda(int lambda, double a, double b){
        double upper = Math.max(Math.min(xi0 + Delta, b), xi0 - Delta);
        double lower = Math.min(Math.max(xi0 - Delta, a), xi0 + Delta);

        return (N/(lambda + 1)) * (Math.pow(upper, lambda + 1) - Math.pow(lower, lambda + 1));
    }

    // Phi : int -> double
    // Purpose: calculates the partial moment weighted by the binding rate
    //         integral xi^lambda f(xi) n(xi, t) dxi from 0 to 1
    public double Phi(int lambda){
        return f1 * J_lambda(lambda + 1, 0.0, 1.0);
    }

    // Psi : int -> double
    // Purpose: Caluclates the weighted partial moments for the unbinding rates
    //    integral xi^lambda g(xi) n(xi, t) dxi from -oo to oo
    public double Psi(int lambda){
        double first_term = g0 * J_lambda(lambda, Double.NEGATIVE_INFINITY, 0.0);
        double second_term = g1 * J_lambda(lambda + 1, 0.0, 1.0);
        double third_term = g2 * J_lambda(lambda + 1, 1.0, Double.POSITIVE_INFINITY);
        double fourth_term = (g1 - g2) * J_lambda(lambda, 1.0, Double.POSITIVE_INFINITY);

        return first_term + second_term + third_term + fourth_term;
    }

    // isometricEquilibriate : double double double -> void
    // Purpose: Sets the initial condition for the 
    public void isometricEquilibriate(double l, double l0, double excitation){
        //double fa = myActiveForceLengthCurve.evalY(l);
        u_normalized = 0.0;
        activation = excitation;
        
        // for the uniform distribution, these are the steady state values
        Q0 = activation * varphi; // times force-length
        Q1 = activation * varphi / 2.0;
        Q2 = activation * varphi / 3.0;

        // update the N, Delta, and xi0.
        updateDerivedStateValues();

        //System.out.println("Max Stiffness = " + maxIsoStiffness/(Units.N/Units.cm) + " N/cm");
        //System.out.println("Current Stiffness = " + activation * maxIsoStiffness/(Units.N/Units.cm) + " N/cm");
    }


    // todo: make these splines
    private double ActiveForceLengthCurve(double stretch){
        return 1.0;
    }
    
    private double PassiveForceLengthCurve(double stretch){
        return 0.0;
    }


    /* =======================================
     * LinearAxialMaterial Methods
     * =======================================
    */

    public double computeF(double l, double ldot, double l0, double excitation){
        // note: if tendon and pennation angle are added, this will need to change
        u_normalized = -ldot/gamma;
        activation = excitation;
        return Gamma1*Q1 + damping*ldot; // plus passive force
    }

    public double computeDFdl(double l, double ldot, double l0, double excitation){
        return Gamma0*Q0; // plus passive stiffness
    }

    public double computeDFdldot(double l, double ldot, double l0, double excitation) {
        return damping;
    }

    public boolean isDFdldotZero() {
        return false;
    }

    
    /* =======================================
     * HasNumericState Functions
     * =======================================
    */

    public boolean hasState(){
        return true;
    }

    public boolean requiresAdvance(){
        return true;
    }

    public void getState (DataBuffer data){
        data.dput(Q0);
        data.dput(Q1);
        data.dput(Q2);
    }

    public void setState (DataBuffer data){
        Q0 = data.dget();
        Q1 = data.dget();
        Q2 = data.dget();
    }

    public void advanceState (double t0, double t1){
        //System.out.println("=> Advancin' the state!");
    }

    public int getStateVersion(){
        return 0;
    }

    public int numAuxVars(){
        return 3;
    }

    public int getAuxVarState(double[] buf, int idx){
        buf[idx]   = Q0;
        buf[idx+1] = Q1;
        buf[idx+2] = Q2;
        return idx + numAuxVars();
    }

    public int setAuxVarState(double[] buf, int idx){
        Q0 = buf[idx];
        Q1 = buf[idx+1];
        Q2 = buf[idx+2];

        // possible add-on for stability
        // make sure these are all above their minimum values
        updateDerivedStateValues();
        //System.out.println("Setting state vars");

        return idx + numAuxVars();
    }

    public double dQ0_func(){
        return activation * Beta(0) - (Phi(0) + Psi(0));
    }

    public double dQ1_func(){
        return activation * Beta(1) - (Phi(1) + Psi(1)) - u_normalized*Q0;
    }

    public double dQ2_func(){
        return activation * Beta(2) - (Phi(2) + Psi(2)) - 2*u_normalized*Q1;
    }

    public int getAuxVarDerivative(double[] buf, int idx){
        // these are the workhorse of this model.
        //double dQ0 = activation * Beta(0) - (Phi(0) + Psi(0));
        //double dQ1 = activation * Beta(1) - (Phi(1) + Psi(1)) - u_normalized*Q0;
        //double dQ2 = activation * Beta(2) - (Phi(2) + Psi(2)) - 2*u_normalized*Q1;

        buf[idx] = dQ0_func();
        buf[idx+1] = dQ1_func();
        buf[idx+2] = dQ2_func();

        //System.out.println("========================");
        //System.out.println("activation = " + activation);
        //System.out.println("dQ0/dt = " + dQ0);
        //System.out.println("dQ1/dt = " + dQ1);
        //System.out.println("dQ2/dt = " + dQ2);
        //System.out.println("========================");

        return idx + numAuxVars();
    }



    /* =======================================
     * Updating Value methods
     * =======================================
    */
    
    private double N_approx(double q0, double q1, double q2){
        double numerator = q0*q0;
        double denominator = 2 * Math.sqrt(3 * (q0*q2 - q1*q1));

        // possible that denominator = 0?
        return numerator / denominator;
    }

    private double xi0_approx(double q0, double q1, double q2){
        if (q0 > Q0_MIN){
            return q1/q0;
        } else {
            System.out.println("Warning: minimum Q0-magnitude not achieved!");
            return 0.0;
        }
    }

    private double Delta_approx(double q0, double q1, double q2){
        return Math.sqrt(3 * (q0*q2 - q1*q1)) / q0;
    }

    // StateFromUniformDistribution : void -> void
    // This is essentially the inverse of updating the derived state variables.
    // Essentially, we can move from N, Delta and xi0 to Q0, Q1, and Q2.
    private void StateFromUniformDistribution(){
        Q0 = 2.0 * Delta * N;
        Q1 = 2.0 * Delta * N * xi0;
        Q2 = (2.0 * Delta * N * (Delta*Delta + 3*xi0*xi0)) / 3.0;
    }

    public void StretchMuscle(double deltaX){
        // nudge xi0 by deltaX, but adjusted for sarcomere lengths
        xi0 += deltaX/gamma;
        //System.out.println("Sarcomere's stretched by " + deltaX/gamma);

        // update Q0, Q1, and Q2 so that force can be calculated accurately
        StateFromUniformDistribution();
    }


}
