import static java.lang.Math.pow;

public class Function {
    public static double f(double x) {
        return 3 * pow(x, 3) - 3 * pow(x, 2) * x + 6 * x / x;
    }

    public static double trueU(double x) {
        return pow(x, 3);
    }
    public static double trueW(double x) {
        return 3 * pow(x, 2);
    }
    public static double trueZ(double x) {
        return 6 * x;
    }
    public static double fz(double x, double u, double w, double z) {
        return 3 * u - w * x + z / x;
    }
}
