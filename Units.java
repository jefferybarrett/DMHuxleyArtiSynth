package artisynth.models.pendulum;


// this is just a unit manager
public class Units{

    private static double tera = 1E12;
    private static double giga = 1E9;
    private static double mega = 1E6;
    private static double kilo = 1E3;
    private static double milli = 1E-3;
    private static double micro = 1E-6;
    private static double nano = 1E-9;
    private static double pico = 1E-12;
    private static double fempto = 1E-15;

    // Base Units
    private static double K = 1.0; // kelvin
    public static double m = 1.0;  // metre
    public static double N = 1.0;  // newton
    public static double s = 1.0;  // second
    public static double Hz = 1.0; // the Hz
    public static double g = 1.0E-3; // gram (base is kg)

    // masses
    public static double kg = kilo * g;
    public static double mg = milli * g;
    public static double ug = micro * g;
    public static double ng = nano * g;
    public static double pg = pico * g;

    // lengths
    public static double km = kilo * m;
    public static double cm = 1.0E-2 * m; // not an SI unit
    public static double mm = milli * m;
    public static double um = micro * m;
    public static double nm = nano * m;
    public static double pm = pico * m;

    
    // areas
    public static double m2 = m * m;
    public static double cm2 = cm * cm;
    public static double mm2 = mm * mm;
    public static double um2 = um * um;
    public static double nm2 = nm * nm;
    public static double pm2 = pm * pm;


    // stresses
    public static double Pa = N / m2;
    public static double kPa = kilo * Pa;
    public static double MPa = mega * Pa;
    public static double GPa = giga * Pa;

    // forces
    public static double kN = kilo * N;
    public static double mN = milli * N;
    public static double uN = micro * N;
    public static double nN = nano * N;
    public static double pN = pico * N;

    // energy
    public static double J = N * m;
    public static double kJ = kilo * J;
    public static double mJ = milli * J;

    // time
    public static double ms = milli * s;
    public static double us = micro * s;
    public static double ns = nano * s;
    public static double ps = pico * s;
    public static double min = 60 * s;
    public static double hr = 60 * min;
    public static double day = 24 * hr;

    // constants
    public static double AvogadrosNumber = 6.022E23;
    public static double BoltzmansConstant = 1.380649E-23 * J / K;

    // amounts
    public static double mol = AvogadrosNumber;
    public static double mmol = milli * mol;
    public static double umol = micro * mol;
    public static double nmol = nano * mol;
    public static double pmol = pico * mol;

    // volumes
    public static double mL = cm*cm*cm;
    public static double L = kilo * mL;
    public static double uL = micro * L;
    public static double nL = nano * L;
    public static double pL = pico * L;

    // torques
    public static double Nm = N * m;
    public static double Nmm = N * mm;
    public static double Num = N * um;

    // frequency
    public static double MHz = mega * Hz;
    public static double kHz = kilo * Hz;
    public static double mHz = milli * Hz;    

    // temperature
    public static double CtoK(double TinC){
        return TinC + 273.5;
    }
}