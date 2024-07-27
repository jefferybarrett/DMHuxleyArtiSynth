package artisynth.models.pendulum;

import artisynth.core.modelbase.HasNumericState;
import artisynth.core.materials.*;
import maspack.util.DataBuffer;


public class StandardLinearSolidAxialMaterial 
    extends LinearAxialMaterial implements HasNumericState {

    /*  This is an implementation of the Standard Linear Viscoelastic Solid
     *  in the Maxwell Arrangment:
     * 
     *             k'     eta
     * |     |---\/\/\/---=]--|
     * |-----|                |-------> F
     * |     |---/\/\/\/\/\/--|
     *               k
     * 
     * We formulate this material using the state-space method (Onat)
     * 
     *      dx0/dt = dxi/dt                             (1)
     *      dx1/dt = dxi/dt - (k'/eta) * x1             (2)
     * 
     * Where xi(t) is the length of the overall material. Further, the force being:
     * 
     *           F = k x0 + k' x1                       (3)
     * 
     * So the material's state is defined by the deflections in both springs (x0, x1),
     * and whose state evolution is governed by Equations 1 and 2.
     * 
     * 
     * By: Jeff M. Barrett
    */

    protected double k = 0.0;
    protected double kp = 0.0;
    protected double eta = 1.0;

    // state variables
    protected double x0 = 0.0;
    protected double x1 = 0.0;
    protected double dxi = 0.0; // hack for storing the velocity
    
    public StandardLinearSolidAxialMaterial(double _k, double _kprime, double _eta){
        k = _k;
        kp = _kprime;
        eta = _eta;

        System.out.println("k   = " + k/(Units.N/Units.cm) + " N/cm");
        System.out.println("k'  = " + kp/(Units.N/Units.cm) + " N/cm");
        System.out.println("eta = " + eta/(Units.N*Units.s/Units.cm) + " Ns/cm");
        System.out.println("tau = " + (eta/kp)/(Units.s) + " seconds");
    }

    public double computeF(double l, double ldot, double l0, double excitation){
        dxi = ldot;
        //double force = k * x0 + kp * x1;
        return k * x0 + kp * x1;
    }

    public double computeDFdl(double l, double ldot, double l0, double excitation){
        return k + kp;
    }

    public double computeDFdldot(double l, double ldot, double l0, double excitation) {
        return 0.0;
    }

    public boolean isDFdldotZero() {
        return true;
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
        data.dput(x0);
        data.dput(x1);
    }

    public void setState (DataBuffer data){
        x0 = data.dget();
        x1 = data.dget();
    }

    public void advanceState (double t0, double t1){
        //System.out.println("=> Advancin' the state!");
    }

    public int getStateVersion(){
        return 0;
    }

    public int numAuxVars(){
        return 2;
    }

    public int getAuxVarState(double[] buf, int idx){
        buf[idx] = x0;
        buf[idx+1] = x1;
        return idx + numAuxVars();
    }

    public int setAuxVarState(double[] buf, int idx){
        x0 = buf[idx];
        x1 = buf[idx+1];
        return idx + numAuxVars();
    }

    public int getAuxVarDerivative(double[] buf, int idx){
        double dx0 = dxi;
        double dx1 = dxi - (kp/eta) * x1;

        buf[idx] = dx0;
        buf[idx+1] = dx1;

        return idx + numAuxVars();
    }

}
