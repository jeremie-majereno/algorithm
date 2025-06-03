package cubic;

import java.util.Scanner;

public class CubicSpline {

    public static void main(String[] args) {
        Scanner scanner = new Scanner(System.in);

        System.out.print("Enter the number of intervals (n): ");
        int n = scanner.nextInt();

        double[] x = new double[n + 1];
        System.out.println("Enter x values (x0 < x1 < ... < xn):");
        for (int i = 0; i <= n; i++) {
            x[i] = scanner.nextDouble();
        }

        double[] a = new double[n + 1];
        System.out.println("Enter function values (a0, a1, ..., an):");
        for (int i = 0; i <= n; i++) {
            a[i] = scanner.nextDouble();
        }

        // Calculate step sizes h[i] = x[i+1] - x[i]
        double[] h = new double[n];
        for (int i = 0; i < n; i++) {
            h[i] = x[i + 1] - x[i];
        }

        // Calculate alpha coefficients
        double[] alpha = new double[n + 1];
        for (int i = 1; i < n; i++) {
            alpha[i] = (3.0 / h[i]) * (a[i + 1] - a[i])
                    - (3.0 / h[i - 1]) * (a[i] - a[i - 1]);
        }

        // Declare arrays for solving the tridiagonal system
        double[] l = new double[n + 1];
        double[] mu = new double[n + 1];
        double[] z = new double[n + 1];
        double[] c = new double[n + 1];

        // Boundary conditions (example: f'(x0) = 2 and f'(xn) = 2)
        l[0] = -1; // Represents a boundary condition at the start
        mu[0] = 0;
        z[0] = 0;

        // Forward elimination to set up the system
        for (int i = 1; i < n; i++) {
            l[i] = 2 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1];
            mu[i] = h[i] / l[i];
            z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
        }

        l[n] = 1; // Represents a boundary condition at the end
        z[n] = 0.0;
        c[n] = 0.0;

        // Back substitution to compute coefficients c
        for (int j = n - 1; j >= 0; j--) {
            c[j] = z[j] - mu[j] * c[j + 1];
        }

        // Calculate b and d coefficients for each interval
        double[] b = new double[n];
        double[] d = new double[n];
        for (int j = 0; j < n; j++) {
            b[j] = (a[j + 1] - a[j]) / h[j] - h[j] * (c[j + 1] + 2 * c[j]) / 3.0;
            d[j] = (c[j + 1] - c[j]) / (3.0 * h[j]);
        }

        System.out.println("\nCubic Spline Coefficients:");
        for (int j = 0; j < n; j++) {
            System.out.printf("Interval [x%d, x%d]:\n", j, j + 1);
            System.out.printf("a%d = %.6f\n", j, a[j]);
            System.out.printf("b%d = %.6f\n", j, b[j]);
            System.out.printf("c%d = %.6f\n", j, c[j]);
            System.out.printf("d%d = %.6f\n\n", j, d[j]);
        }

        // Evaluate the spline at a specific point x0
        double x0 = 2.5; // You can adjust this point as needed
        double splineValue = evaluateSpline(x, a, b, c, d, n, x0);
        System.out.printf("The value of the spline at x = %.6f is %.6f\n", x0, splineValue);

        scanner.close();
    }

    // Function to evaluate the spline at x0
    private static double evaluateSpline(double[] x, double[] a, double[] b, double[] c, double[] d, int n, double x0) {
        // Check if x0 is within the interval [x[0], x[n]]
        if (x0 < x[0] || x0 > x[n]) {
            System.out.println("x0 is out of bounds.");
            return Double.NaN;
        }

        // Find the interval [x[i], x[i+1]] containing x0
        int i = 0;
        while (i < n - 1 && x[i + 1] < x0) {
            i++;
        }

        double dx = x0 - x[i];
        return a[i] + b[i] * dx + c[i] * dx * dx + d[i] * dx * dx * dx;
    }
}
