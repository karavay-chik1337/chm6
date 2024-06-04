import static java.lang.Math.pow;

public class Function {
    //u'''(x)=3u-u'-xu'+2u''+1/2 xu''
    public static double exactU(double x) {
        return pow(x, 3) - 6 * x;
    }

    public static double exactW(double x) {
        return 3 * pow(x, 2) - 6;
    }

    public static double exactZ(double x) {
        return 6 * x;
    }

    public static double f(double x, double u, double w, double z) {
        return 3 * u - w - w * x + 2 * z + 0.5 * x * z;
    }
}
