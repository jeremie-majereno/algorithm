package cubic;

public class ClampedCubicSpline {

    public static void main(String[] args) {
        // Example input data (hardcoded for testing)
        int n = 3;
        double[] x = {1,2,3,4};
        double[] a = {1,8,27,64};
        double FPO = 3;
        double FPN = 48;
        double x0 = 2.5;

        double[] h = new double[n];
        for (int i = 0; i < n; i++) {
            h[i] = x[i + 1] - x[i];
        }

        // Calculate alpha values
        double[] alpha = new double[n + 1];
        alpha[0] = 3.0 * ((a[1] - a[0]) / h[0] - FPO);
        alpha[n] = 3.0 * (FPN - (a[n] - a[n - 1]) / h[n - 1]);

        for (int i = 1; i < n; i++) {
            alpha[i] = (3.0 / h[i]) * (a[i + 1] - a[i]) - (3.0 / h[i - 1]) * (a[i] - a[i - 1]);
        }

        // Solve tridiagonal system
        double[] l = new double[n + 1];
        double[] mu = new double[n + 1];
        double[] z = new double[n + 1];
        double[] c = new double[n + 1];

        // Forward elimination
        l[0] = 2 * h[0];
        mu[0] = 0.5;
        z[0] = alpha[0] / l[0];

        for (int i = 1; i < n; i++) {
            l[i] = 2 * (h[i] + h[i - 1]) - h[i - 1] * mu[i - 1];
            mu[i] = h[i] / l[i];
            z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
        }

        l[n] = h[n - 1] * (2 - mu[n - 1]);
        z[n] = (alpha[n] - h[n - 1] * z[n - 1]) / l[n];
        c[n] = z[n];

        // Back substitution
        for (int j = n - 1; j >= 0; j--) {
            c[j] = z[j] - mu[j] * c[j + 1];
        }

        // Calculate b and d coefficients
        double[] b = new double[n];
        double[] d = new double[n];

        for (int j = 0; j < n; j++) {
            b[j] = (a[j + 1] - a[j]) / h[j] - h[j] * (c[j + 1] + 2 * c[j]) / 3.0;
            d[j] = (c[j + 1] - c[j]) / (3.0 * h[j]);
        }

        // Print results
        System.out.println("\nClamped Cubic Spline Coefficients:");
        for (int j = 0; j < n; j++) {
            System.out.printf("Interval [%.1f, %.1f]:%n", x[j], x[j+1]);
            System.out.printf("a: %.6f%n", a[j]);
            System.out.printf("b: %.6f%n", b[j]);
            System.out.printf("c: %.6f%n", c[j]);
            System.out.printf("d: %.6f%n%n", d[j]);
        }

        // Evaluate at x0 = 2.5
        double splineValue = evaluateSpline(x, a, b, c, d, n, x0);
        System.out.printf("The value of the spline at x = %.1f is %.6f%n", x0, splineValue);
    }

    private static double evaluateSpline(double[] x, double[] a, double[] b, double[] c, double[] d, int n, double x0) {
        // Find the interval containing x0
        int j = 0;
        while (j < n && x[j + 1] < x0) {
            j++;
        }
        if (j >= n) {
            j = n - 1;
        }

        double dx = x0 - x[j];
        return a[j] + b[j] * dx + c[j] * dx * dx + d[j] * dx * dx * dx;
    }
}