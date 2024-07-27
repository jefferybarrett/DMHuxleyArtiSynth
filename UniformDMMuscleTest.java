package artisynth.models.pendulum;

import maspack.util.*;

public class UniformDMMuscleTest extends UnitTest {

    // test Stiffness
    // In theory, stretching the muscle by deltaX should result in a stiffness
    // given by (s0/h) F/X0. This should be a way to predict that
    // Here we test 5 different deltaX's and 5 different levels of activation
    public void testStiffness(UniformDMMuscle mat){

        int nsamps = 2;

        // musclelength
        double fibreLength = mat.X_0;
        double s0 = mat.s_sarcomere;
        double h = mat.h_myosin;
        double q = (s0 / h);

        // min and max activation
        double act_min = 0.0;
        double act_max = 1.0;

        // min and max displacements
        double dx_min = 0.01;
        double dx_max = 0.5;

        for (double dx = dx_min; dx <= dx_max; dx += (dx_max - dx_min)/(nsamps+1)){
            for (double act = act_min; act <= act_max; act += (act_max - act_min)/(nsamps+1)){

                // equilibriate at the level of activation
                mat.isometricEquilibriate(fibreLength, 0.0, act);
                double preForce = mat.computeF(fibreLength, 0.0, fibreLength, act);

                // then nudge shift the distribution according to
                mat.StretchMuscle(dx);
                double postForce = mat.computeF(fibreLength + dx, 0.0, fibreLength, act);

                double approximateStiffness = (postForce - preForce) / dx;
                double analyticStiffness = q * preForce / fibreLength;

                double err = (analyticStiffness - approximateStiffness) / approximateStiffness;

                if (err > 1e-6){
                    throw new TestException(
                        "Error exceeded tolerance for activation = " + act + " and dx = " + dx
                    );
                }
                //System.out.println("Approximate Stiffness = " + approximateStiffness);
                //System.out.println("Analytic Stiffness = " + analyticStiffness);
            }
        }
    }

    private double beta_numeric(int lambda, double f1){
        return f1/(lambda + 2);
    }

    private double Phi_numeric(double N, double Delta, double xi0){
        return 0.0;
    }

    private double Psi_numeric(double N, double Delta, double xi0){
        return 0.0;
    }

    public void testDerivatives(UniformDMMuscle mat){

        double f1 = mat.f1;
        double g0 = mat.g0;
        double g1 = mat.g1;
        double g2 = mat.g2;
        double h = mat.h_myosin;

        UniformNumericApproximator approximator = new UniformNumericApproximator(
            50000,
            -5.0,
            15.0,
            f1, g0, g1, g2
        );

        double[] x = {0.0, 0.25, 0.5, 0.75, 1.0, 1.25};
        double[] y = {0.0, 0.25*0.25, 0.25, 0.75*0.75, 1.0, 1.25*1.25};
        // TrapezoidIntegrator

        System.out.println("Q_0 ~= " + approximator.Q_numeric(0));
        System.out.println("Q_1 ~= " + approximator.Q_numeric(1));
        System.out.println("Q_2 ~= " + approximator.Q_numeric(2));

        System.out.println("Q_0 = " + approximator.Q_analytic(0));
        System.out.println("Q_1 = " + approximator.Q_analytic(1));
        System.out.println("Q_2 = " + approximator.Q_analytic(2));

        System.out.println("Beta_0 ~= " + approximator.Beta_numeric(0));
        System.out.println("Beta_1 ~= " + approximator.Beta_numeric(1));
        System.out.println("Beta_2 ~= " + approximator.Beta_numeric(2));

        System.out.println("Beta_0 = " + approximator.Beta_analytic(0));
        System.out.println("Beta_1 = " + approximator.Beta_analytic(1));
        System.out.println("Beta_2 = " + approximator.Beta_analytic(2));

        System.out.println("Phi_0 ~= " + approximator.Phi_numeric(0));
        System.out.println("Phi_1 ~= " + approximator.Phi_numeric(1));
        System.out.println("Phi_2 ~= " + approximator.Phi_numeric(2));

        System.out.println("Phi_0 = " + approximator.Phi_analytic(0));
        System.out.println("Phi_1 = " + approximator.Phi_analytic(1));
        System.out.println("Phi_2 = " + approximator.Phi_analytic(2));
        
        System.out.println("Psi_0 ~= " + approximator.Psi_numeric(0));
        System.out.println("Psi_1 ~= " + approximator.Psi_numeric(1));
        System.out.println("Psi_2 ~= " + approximator.Psi_numeric(2));

        System.out.println("Psi_0 = " + approximator.Psi_analytic(0));
        System.out.println("Psi_1 = " + approximator.Psi_analytic(1));
        System.out.println("Psi_2 = " + approximator.Psi_analytic(2));
        

        int nsteps = 5;

        double N_min = 0.01;
        double N_max = 1.0;

        double Delta_min = 0.01;
        double Delta_max = 1.0;

        double xi0_min = -5.0;
        double xi0_max = 5.0;

        double ldot_min = -0.1;
        double ldot_max = 0.1;

        double q0_numeric = 0.0;
        double q1_numeric = 0.0;
        double q2_numeric = 0.0;
        double q0_model = 0.0;
        double q1_model = 0.0;
        double q2_model = 0.0;

        double q_tolerance = 1e-2;
        double err = 0.0;

        for(double N = N_min; N <= N_max; N += (N_max - N_min)/nsteps){
            for(double Delta = Delta_min; Delta <= Delta_max; Delta += (Delta_max - Delta_min)/nsteps){
                for(double xi0 = xi0_min; xi0 <= xi0_max; xi0 += (xi0_max - xi0_min)/nsteps){
                    double q0 = 2.0 * Delta * N;
                    double q1 = 2.0 * Delta * N * xi0;
                    double q2 = (2.0 * Delta * N * (Delta*Delta + 3*xi0*xi0)) / 3.0;
                    approximator.setDistribution(N, Delta, xi0);
                    mat.setQstate(q0, q1, q2, 0.0);

                    //System.out.println("N     = " + N);
                    //System.out.println("Delta = " + Delta);
                    //System.out.println("xi0   = " + xi0);

                    //System.out.println("Checking Qs");
                    checkError(q_tolerance, mat.Q0, approximator.Q_numeric(0));
                    checkError(q_tolerance, mat.Q1, approximator.Q_numeric(1));
                    checkError(q_tolerance, mat.Q2, approximator.Q_numeric(2));
                    
                    //System.out.println("Checking Betas");
                    checkError(q_tolerance, mat.Beta(0), approximator.Beta_numeric(0));
                    checkError(q_tolerance, mat.Beta(1), approximator.Beta_numeric(1));
                    checkError(q_tolerance, mat.Beta(2), approximator.Beta_numeric(2));

                    //System.out.println("Checking Phis");
                    checkError(q_tolerance, mat.Phi(0), approximator.Phi_numeric(0));
                    checkError(q_tolerance, mat.Phi(1), approximator.Phi_numeric(1));
                    checkError(q_tolerance, mat.Phi(2), approximator.Phi_numeric(2));

                    //System.out.println("Checking Psis");
                    checkError(q_tolerance, mat.Psi(0), approximator.Psi_numeric(0));
                    checkError(q_tolerance, mat.Psi(1), approximator.Psi_numeric(1));
                    checkError(q_tolerance, mat.Psi(2), approximator.Psi_numeric(2));


                    //checkError(q_tolerance, mat.dQ0_func(), approximator.dQ_numeric(0, 0.0));
                    //checkError(q_tolerance, mat.dQ1_func(), approximator.dQ_numeric(1, 0.0));
                    //checkError(q_tolerance, mat.dQ1_func(), approximator.dQ_numeric(2, 0.0));


                    //for (double ldot = ldot_min; ldot <= ldot_max; ldot += (ldot_max - ldot_min)/nsteps){
                    //    mat.setQstate(q0, q1, q2, ldot);
                    //}
                }
            }
        }

    }
    
    public int checkError(double tolerance, double val0, double val1){
        double err = (val1 - val0)/val0;
        if (err > tolerance){
            System.out.println("Problem!");
            System.out.println("val_0 = " + val0);
            System.out.println("val_1 = " + val1);
            return -1;
        } else{
            return 0;
        }
    }

    public void test(){
        testStiffness(new UniformDMMuscle());
        testDerivatives(new UniformDMMuscle());
    }


    public static void main (String[] args) {
        UniformDMMuscleTest tester = new UniformDMMuscleTest();
        tester.runtest();
        
        //System.out.println("Hello World!");
     }

     
    public class UniformNumericApproximator{

        // binding rates
        protected double f1 = 15.0 * Units.Hz; // Binding rate between 0 and h

        // detachment rates
        protected double g0 = 170.0 * Units.Hz; // Dissociation rate constant for compression
        protected double g1 = 8.0 * Units.Hz;   // Dissociation rate for short displacements
        protected double g2 = 25.0 * Units.Hz;  // Dissociation rate for large displacements
        protected double h_myosin = 1.0;//5.0 * Units.nm;    // characteristic bond-length

        // distribution properties
        protected double Height = 1.0;   // Height
        protected double Delta = 1.0;   // spread
        protected double mu = 1.0;      // centre

        protected int N = 1000;
        protected double xmin = -5.0;
        protected double xmax = 5.0;
        double dx = (10.0)/1001;

        protected double[] x; // this is the x-discretization
        protected double[] n; // bond distribution

        // the rate functions here
        protected double[] G;
        protected double[] F;

        public UniformNumericApproximator(int _N, double _xmin, double _xmax, double _f1, double _g0, double _g1, double _g2) {
            N = _N;
            xmin = _xmin;
            xmax = _xmax;
            dx = (xmax - xmin)/(N-1);

            x = new double[N];
            for (int i = 0; i < N; ++i){
                x[i] = xmin + i * dx;
            }
            //System.out.println("x[0] = " + x[0]);
            //System.out.println("x[1] = " + x[1]);
            //System.out.println("dx = " + dx);
            //System.out.println("x[N-1] = " + x[N-1]);

            f1 = _f1;
            g0 = _g0;
            g1 = _g1;
            g2 = _g2;
            evalFunctions();
        }

        public void setDistribution(double _Height, double _Delta, double _mu){
            Height = _Height;
            Delta = _Delta;
            mu = _mu;
            updateXCoordinates();
            evalFunctions();
        }

        public void updateXCoordinates(){
            xmin = mu - Delta;
            xmax = mu + Delta;
            dx = (xmax - xmin)/(N-1);
            x = new double[N];
            for (int i = 0; i < N; ++i){
                x[i] = xmin + i * dx;
            }
        }


        public double evalDistribution(double x){
            if ((x > mu - Delta) & (x < mu + Delta)){
                return Height;
            } else{
                return 0.0;
            }
        }

        public double evalf(double x){
            if ((x > 0) & (x <= h_myosin)){
                return f1 * x / h_myosin;
            } else{
                return 0.0;
            }
        }

        public double evalg(double x){
            if (x < 0){
                return g0;
            } else if ((x >= 0) & (x <= h_myosin)){
                return g1 * x / h_myosin;
            } else{
                return g2 * (x/h_myosin - 1) + g1;
            }
        }

        public void evalFunctions(){
            G = new double[N];
            F = new double[N];
            n = new double[N];
            for (int i = 0; i < N; ++i){
                G[i] = evalg(x[i]);
                F[i] = evalf(x[i]);
                n[i] = evalDistribution(x[i]);
            }
        }

        public double ForwardEulerIntegrate(double[] x, double[] y){
            assert(x.length == y.length);
            int M = x.length;
            
            double retVal = y[0] * (x[1] - x[0]);

            for (int i = 1; i < M; ++i){
                retVal += y[i] * (x[i] - x[i-1]);
            }

            return retVal;
        }

        public double TrapezoidIntegrator(double[] x, double[] y){
            assert(x.length == y.length);
            int M = x.length;
            double retVal = 0.0;

            for (int i = 1; i < M; ++i){
                retVal += 0.5 * (y[i] + y[i-1]) * (x[i] - x[i-1]);
            }

            return retVal;
        }

        public double SimpsonsRule(double[] x, double[] y){
            assert((x.length == y.length) & (x.length > 1));
            int M = x.length;
            double dx = (x[M - 1] - x[0]) / (M - 1);
            double retVal = y[0] + y[M-1]; // start with the boudnarys

            for(int i = 1; i < M-1; i+=2){
                retVal += 4*y[i];
            }
            for(int i = 2; i < M-1; i+=2){
                retVal += 2*y[i];
            }
            return (dx/3.0) * retVal;
        }

        public double Q_numeric(int lambda){
            if (lambda >= 0){
                double[] integrand = new double[N];
                for (int i = 0; i < N; ++i){
                    integrand[i] = Math.pow(x[i], lambda) * n[i];
                }
                return SimpsonsRule(x, integrand);
            } else{
                return 0.0;
            }
        }

        public double Q_analytic(int lambda){
            if (lambda >= 0){
                return Height/(lambda + 1) * (Math.pow(mu + Delta, lambda+1) - Math.pow(mu - Delta, lambda+1));
            } else{
                return 0.0;
            }
        }

        public double Beta_numeric(int lambda){
            double[] integrand = new double[N];
            for (int i = 0; i < N; ++i){
                integrand[i] = Math.pow(x[i], lambda) * F[i];
            }
            return SimpsonsRule(x, integrand);
        }

        public double Beta_analytic(int lambda){
            return h_myosin*h_myosin*f1/(lambda + 2);
        }

        public double Phi_numeric(int lambda){
            double[] integrand = new double[N];
            for (int i = 0; i < N; ++i){
                integrand[i] = Math.pow(x[i], lambda) * F[i] * n[i];
            }
            return SimpsonsRule(x, integrand);
        }

        public double Psi_numeric(int lambda){
            double[] integrand = new double[N];
            for (int i = 0; i < N; ++i){
                integrand[i] = Math.pow(x[i], lambda) * G[i] * n[i];
            }
            return SimpsonsRule(x, integrand);
        }

        private double J_lambda(int lambda, double a, double b){
            double upper = Math.max(Math.min(mu + Delta, b), mu - Delta);
            double lower = Math.min(Math.max(mu - Delta, a), mu + Delta);
    
            return (Height/(lambda + 1)) * (Math.pow(upper, lambda + 1) - Math.pow(lower, lambda + 1));
        }

        public double Phi_analytic(int lambda){
            return f1 * J_lambda(lambda + 1, 0.0, 1.0);
        }

        public double Psi_analytic(int lambda){
            double term0 = g0 * J_lambda(lambda, Double.NEGATIVE_INFINITY, 0.0);
            double term1 = g1 * J_lambda(lambda + 1, 0.0, 1.0);
            double term2 = g2 * J_lambda(lambda + 1, 1.0, Double.POSITIVE_INFINITY);
            double term3 = (g1 - g2) * J_lambda(lambda, 1.0, Double.POSITIVE_INFINITY);
            return term0 + term1 + term2 + term3;
        }

        public double dQ_numeric(int lambda, double udot){
            return Beta_numeric(lambda) - (Phi_numeric(lambda) + Psi_numeric(lambda)) - lambda * udot * Q_numeric(lambda - 1);
        }

        public double dQ_analytic(int lambda, double udot){
            return Beta_analytic(lambda) - (Phi_analytic(lambda) + Psi_analytic(lambda)) - lambda * udot * Q_analytic(lambda -1 );
        }
    }

}

