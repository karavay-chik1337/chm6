import java.util.ArrayList;
import java.util.List;

import static java.lang.Math.abs;

public class ShootingMethod {
    double a;
    double b;
    int n;
    double eps;
    int K;
    double alpha0;
    double A;
    double B;
    double C;

    double[] x;
    double[] u;
    double[] w;
    double[] z;
    double h;

    public ShootingMethod() {}

    public ShootingMethod(double a, double b, int N, double eps, int K, double alpha0) {
        this.a = a;
        this.b = b;
        this.n = N;
        this.eps = eps;
        this.K = K;
        this.alpha0 = alpha0;
        this.A = Function.exactU(a);
        this.B = Function.exactZ(a);
        this.C = Function.exactW(b);

        this.h = (b - a) / n;

        x = new double[n + 1];
        x[0] = a;
        for (int i = 1; i <= n; i++)//равномерное разбиение
            x[i] = a + (h * i);
    }

    public void bisectionMethod() {
        long L = 0;
        int ier = 0;
        double delta = 2 * alpha0;
        double alpha1 = alpha0 + delta;
        double phi1 = RungeKutt(alpha0);
        double phi2 = RungeKutt(alpha1);
        if (phi1 * phi2 > 0) {
            alpha1 = alpha0 - delta;
            phi2 = RungeKutt(alpha1);
        }
        double alpha = 0;
        while (L < K && Math.abs(alpha1 - alpha0) > 2 * eps) {
            alpha = (alpha0 + alpha1) / 2;
            double phi3 = RungeKutt(alpha);
            if (phi1 * phi3 < 0) {
                alpha1 = alpha;
            } else if (phi2 * phi3 < 0) {
                alpha0 = alpha;
            } else if (abs(phi3) < eps) {
                break;
            }
            alpha = (alpha0 + alpha1) / 2;
            L++;
        }
        if (L == K) {
            ier = 1;
        }
        System.out.println("L: " + L + "\nalpha: " + alpha + "\nIER: " + ier);
//        int L = 0;
//        List<Double> alpha = new ArrayList<>();
//        alpha.add(alpha0);
//        alpha.add(alpha1);
//        int alphaPos = 0;
//
//        List<Double> target = new ArrayList<>();
//        target.add(RungeKutt(alpha.get(alphaPos)));
//        alphaPos++;
//        target.add(RungeKutt(alpha.get(alphaPos)));
//
//        System.out.println("\nSecant method:");
//        while (L < K && target.getLast() > eps) {
//            alpha.add(alpha.get(alphaPos) - target.get(alphaPos) * (alpha.get(alphaPos) - alpha.get(alphaPos - 1)) / (target.get(alphaPos) - target.get(alphaPos - 1)));
//            target.add(RungeKutt(alpha.getLast()));
//            alphaPos++;
//            L++;
//        }
//
//        if (L >= K) {
//            System.out.println("IER: " + 1);
//        }
//        else {
//            System.out.println("IER: " + 0);
//            System.out.println("Iterations: " + L);
//            System.out.println("Final alpha: " + alpha.getLast());
//        }

    }

    public double RungeKutt(double curAlpha) {
        u = new double[n + 1];
        w = new double[n + 1];
        z = new double[n + 1];
        u[0] = A;
        w[0] = curAlpha;
        z[0] = B;

        for (int i = 1; i <= n; i++) {
            if (i == 1) {
                System.out.printf("%-10s \t %-10s \t %-10s \t %-10s \t %-10s \t %-10s \t %-10s \n", "x", "u", "du", "w", "dw", "z", "dz");
                System.out.printf("%-10f \t %-10f \t %-10f \t %-10f \t %-10f \t %-10f \t %-10f \n", x[0], u[0],
                        abs(u[0] - Function.exactU(x[0])), w[0], abs(w[0] - Function.exactW(x[0])), z[0], abs(z[0] - Function.exactZ(x[0])));
            }
            double k11 = h * w[i-1];
            double k12 = h * z[i-1];
            double k13 = h * Function.f(x[i-1], u[i-1], w[i-1], z[i-1]);

            double k21 = h * (w[i-1] + k12 / 2);
            double k22 = h * (z[i-1] + k13 / 2);
            double k23 = h * Function.f(x[i-1] + h / 2, u[i-1] + k11 / 2, w[i-1] + k12 / 2, z[i-1] + k13 / 2);

            double k31 = h * (w[i-1] + k22 / 2);
            double k32 = h * (z[i-1] + k23 / 2);
            double k33 = h * Function.f(x[i-1] + h / 2, u[i-1] + k21 / 2, w[i-1] + k22 / 2, z[i-1] + k23 / 2);

            double k41 = h * (w[i-1] + k32);
            double k42 = h * (z[i-1] + k33);
            double k43 = h * Function.f(x[i-1] + h, u[i-1] + k31, w[i-1] + k32, z[i-1] + k33);

            u[i] = u[i-1] + (k11 + 2 * k21 + 2 * k31 + k41) / 6;
            w[i] = w[i-1] + (k12 + 2 * k22 + 2 * k32 + k42) / 6;
            z[i] = z[i-1] + (k13 + 2 * k23 + 2 * k33 + k43) / 6;
            System.out.printf("%-10f \t %-10f \t %-10f \t %-10f \t %-10f \t %-10f \t %-10f \n", x[i], u[i], abs(u[i] - Function.exactU(x[i])), w[i], abs(w[i] - Function.exactW(x[i])), z[i], abs(z[i] - Function.exactZ(x[i])));
        }

        return w[n] - C;
    }
}
