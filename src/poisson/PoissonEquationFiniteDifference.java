package poisson;


import java.util.Arrays;

public class PoissonEquationFiniteDifference {

    public static void main(String[] args) {
        // Input parameters:
        // x-domain: [a, b] and y-domain: [c, d]
        double a = 0.0, b = 0.5;
        double c = 0.0, d = 0.5;
        int n = 10, m = 10;       // Number of intervals along x and y
        double TOL = 1e-10;       // Convergence tolerance
        int MAX_ITER = 10000;     // Maximum number of iterations

        // Solve the Poisson (or Laplace) equation using finite difference method
        solvePoisson(a, b, c, d, n, m, TOL, MAX_ITER);
    }

    // Source function f(x, y). Here, for the Laplace equation, f is zero.
    public static double f(double x, double y) {
        return 0.0;
    }

    // Boundary condition function g(x, y, a, b, c, d)
    // Left (x=a) and bottom (y=c) boundaries set to 0
    // Right (x=b): u(x,y) = 200*y; Top (y=d): u(x,y) = 200*x
    public static double g(double x, double y, double a, double b, double c, double d) {
        if (x == a || y == c) return 0.0;
        if (x == b) return 200 * y;
        if (y == d) return 200 * x;
        return 0.0;
    }

    public static void solvePoisson(double a, double b, double c, double d,
                                    int n, int m, double TOL, int MAX_ITER) {
        // Step 1: Calculate step sizes in x (h) and y (k) directions
        double h = (b - a) / n;
        double k = (d - c) / m;
        double lambda = (h * h) / (k * k);
        double mu = 2 * (1 + lambda);

        // Steps 2-3: Define mesh points in x- and y-directions
        double[] x = new double[n + 1];
        double[] y = new double[m + 1];
        for (int i = 0; i <= n; i++) {
            x[i] = a + i * h;
        }
        for (int j = 0; j <= m; j++) {
            y[j] = c + j * k;
        }

        // Step 4: Initialize the solution grid w (including boundaries) to zero.
        double[][] w = new double[n + 1][m + 1];
        for (int i = 0; i <= n; i++) {
            Arrays.fill(w[i], 0.0);
        }

        // Apply boundary conditions. Here we update grid points on the boundary.
        for (int i = 0; i <= n; i++) {
            for (int j = 0; j <= m; j++) {
                if (i == 0 || i == n || j == 0 || j == m) {
                    w[i][j] = g(x[i], y[j], a, b, c, d);
                }
            }
        }

        // Steps 5-6: Gauss-Seidel iteration
        int iter = 0;
        double NORM;
        do {
            NORM = 0.0;

            // Step 7: Update first corner point in the top row (at (1, m-1))
            double z = (-h * h * f(x[1], y[m - 1])
                    + g(a, y[m - 1], a, b, c, d)
                    + lambda * g(x[1], d, a, b, c, d)
                    + lambda * w[1][m - 2]
                    + w[2][m - 1]) / mu;
            NORM = Math.abs(z - w[1][m - 1]);
            w[1][m - 1] = z;

            // Step 8: Update middle points in the top row (i = 2 to n-2, with j = m-1)
            for (int i = 2; i <= n - 2; i++) {
                z = (-h * h * f(x[i], y[m - 1])
                        + lambda * g(x[i], d, a, b, c, d)
                        + w[i - 1][m - 1]
                        + w[i + 1][m - 1]
                        + lambda * w[i][m - 2]) / mu;
                NORM = Math.max(NORM, Math.abs(w[i][m - 1] - z));
                w[i][m - 1] = z;
            }

            // Step 9: Update last corner point in the top row (at (n-1, m-1))
            z = (-h * h * f(x[n - 1], y[m - 1])
                    + g(b, y[m - 1], a, b, c, d)
                    + lambda * g(x[n - 1], d, a, b, c, d)
                    + w[n - 2][m - 1]
                    + lambda * w[n - 1][m - 2]) / mu;
            NORM = Math.max(NORM, Math.abs(w[n - 1][m - 1] - z));
            w[n - 1][m - 1] = z;

            // Steps 10-13: Update the middle rows, for j = m-2 down to 2
            for (int j = m - 2; j >= 2; j--) {
                // Step 11: Update first column of the current row (i = 1)
                z = (-h * h * f(x[1], y[j])
                        + g(a, y[j], a, b, c, d)
                        + lambda * w[1][j + 1]
                        + lambda * w[1][j - 1]
                        + w[2][j]) / mu;
                NORM = Math.max(NORM, Math.abs(w[1][j] - z));
                w[1][j] = z;

                // Step 12: Update middle columns (i = 2 to n-2)
                for (int i = 2; i <= n - 2; i++) {
                    z = (-h * h * f(x[i], y[j])
                            + w[i - 1][j]
                            + lambda * w[i][j + 1]
                            + w[i + 1][j]
                            + lambda * w[i][j - 1]) / mu;
                    NORM = Math.max(NORM, Math.abs(w[i][j] - z));
                    w[i][j] = z;
                }

                // Step 13: Update last column (i = n - 1)
                z = (-h * h * f(x[n - 1], y[j])
                        + g(b, y[j], a, b, c, d)
                        + w[n - 2][j]
                        + lambda * w[n - 1][j + 1]
                        + lambda * w[n - 1][j - 1]) / mu;
                NORM = Math.max(NORM, Math.abs(w[n - 1][j] - z));
                w[n - 1][j] = z;
            }

            // Step 14: Update first corner point in the bottom row (at (1, 1))
            z = (-h * h * f(x[1], y[1])
                    + g(a, y[1], a, b, c, d)
                    + lambda * g(x[1], c, a, b, c, d)
                    + lambda * w[1][2]
                    + w[2][1]) / mu;
            NORM = Math.max(NORM, Math.abs(w[1][1] - z));
            w[1][1] = z;

            // Step 15: Update middle points in the bottom row (i = 2 to n-2, with j = 1)
            for (int i = 2; i <= n - 2; i++) {
                z = (-h * h * f(x[i], y[1])
                        + lambda * g(x[i], c, a, b, c, d)
                        + w[i - 1][1]
                        + lambda * w[i][2]
                        + w[i + 1][1]) / mu;
                NORM = Math.max(NORM, Math.abs(w[i][1] - z));
                w[i][1] = z;
            }

            // Step 16: Update last corner point in the bottom row (at (n - 1, 1))
            z = (-h * h * f(x[n - 1], y[1])
                    + g(b, y[1], a, b, c, d)
                    + lambda * g(x[n - 1], c, a, b, c, d)
                    + w[n - 2][1]
                    + lambda * w[n - 1][2]) / mu;
            NORM = Math.max(NORM, Math.abs(w[n - 1][1] - z));
            w[n - 1][1] = z;

            iter++;
            if (iter % 100 == 0) {
                System.out.printf("Iteration %d, Norm: %.8f%n", iter, NORM);
            }

            // Step 17: Check for convergence
            if (NORM <= TOL) {
                System.out.println("\nConverged after " + iter + " iterations");
                // Step 18: Output the solution for interior points
                System.out.println("\nSolution (interior points):");
                for (int i = 1; i < n; i++) {
                    for (int j = 1; j < m; j++) {
                        System.out.printf("u(%.3f, %.3f) = %.6f%n", x[i], y[j], w[i][j]);
                    }
                }
                return;
            }

        } while (iter < MAX_ITER);

        // Step 21: If maximum iterations are exceeded, notify the user.
        System.out.println("Maximum number of iterations exceeded");
    }
}
