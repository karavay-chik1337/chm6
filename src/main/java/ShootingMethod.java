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
    double alpha1;
    double A;
    double B;
    double C;

    double[] x;
    double[] u;
    double[] w;
    double[] z;
    double h;

    public ShootingMethod() {}

    public ShootingMethod(double a, double b, int N, double eps, int K, double alpha0, double alpha1, double A, double B, double C) {
        this.a = a;
        this.b = b;
        this.n = N;
        this.eps = eps;
        this.K = K;
        this.alpha0 = alpha0;
        this.alpha1 = alpha1;
        this.A = A;
        this.B = B;
        this.C = C;

        this.h = (b - a) / n;

        x = new double[n + 1];
        x[0] = a;
        for (int i = 1; i <= n; i++)
            x[i] = a + (h * i);

    }

    public void method() {
        int L = 0;
        List<Double> alpha = new ArrayList<>();
        alpha.add(alpha0);
        alpha.add(alpha1);
        int alphaPos = 0;

        List<Double> target = new ArrayList<>();
        target.add(RungeKutt(alpha.get(alphaPos)));
        alphaPos++;
        target.add(RungeKutt(alpha.get(alphaPos)));

        System.out.println("\nSecant method:");
        while (L < K && target.getLast() > eps) {
            alpha.add(alpha.get(alphaPos) - target.get(alphaPos) * (alpha.get(alphaPos) - alpha.get(alphaPos - 1)) / (target.get(alphaPos) - target.get(alphaPos - 1)));
            target.add(RungeKutt(alpha.getLast()));
            alphaPos++;
            L++;
        }

        if (L >= K) {
            System.out.println("IER: " + 1);
        }
        else {
            System.out.println("IER: " + 0);
            System.out.println("Iterations: " + L);
            System.out.println("Final alpha: " + alpha.getLast());
        }

    }

    public double RungeKutt(double curAlpha) {
        u = new double[n + 1];
        w = new double[n + 1];
        z = new double[n + 1];
        u[0] = A;
        w[0] = B;
        z[0] = curAlpha;

        for (int i = 1; i <= n; i++) {
            if (i == 1) {
                System.out.printf("%-10s \t %-10s \t %-10s \t %-10s \t %-10s \t %-10s \t %-10s \n", "x", "u", "du", "w", "dw", "z", "dz");
                System.out.printf("%-10f \t %-10f \t %-10f \t %-10f \t %-10f \t %-10f \t %-10f \n", x[0], u[0],
                        abs(u[0] - Function.trueU(x[0])), w[0], abs(w[0] - Function.trueW(x[0])), z[0], abs(z[0] - Function.trueZ(x[0])));
            }
            double k11 = h * w[i-1];
            double k21 = h * z[i-1];
            double k31 = h * Function.fz(x[i-1], u[i-1], w[i-1], z[i-1]);

            double k12 = h * (w[i-1] + k21 / 2);
            double k22 = h * (z[i-1] + k31 / 2);
            double k32 = h * Function.fz(x[i-1] + h / 2, u[i-1] + k11 / 2, w[i-1] + k21 / 2, z[i-1] + k31 / 2);

            double k13 = h * (w[i-1] + k22 / 2);
            double k23 = h * (z[i-1] + k32 / 2);
            double k33 = h * Function.fz(x[i-1] + h / 2, u[i-1] + k12 / 2, w[i-1] + k22 / 2, z[i-1] + k32 / 2);

            double k14 = h * (w[i-1] + k23);
            double k24 = h * (z[i-1] + k33);
            double k34 = h * Function.fz(x[i-1] + h, u[i-1] + k13, w[i-1] + k23, z[i-1] + k33);

            u[i] = u[i-1] + (k11 + 2 * k12 + 2 * k13 + k14) / 6;
            w[i] = w[i-1] + (k21 + 2 * k22 + 2 * k23 + k24) / 6;
            z[i] = z[i-1] + (k31 + 2 * k32 + 2 * k33 + k34) / 6;
            System.out.printf("%-10f \t %-10f \t %-10f \t %-10f \t %-10f \t %-10f \t %-10f \n", x[i], u[i], abs(u[i] - Function.trueU(x[i])), w[i], abs(w[i] - Function.trueW(x[i])), z[i], abs(z[i] - Function.trueZ(x[i])));
        }

        return abs(u[n] - C);
    }
}
