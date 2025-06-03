package poisson;

import java.util.Arrays;

public class PoissonEquation {

    // Function for f(x, y) in the Poisson equation
    public static double f(double x, double y) {
        return 0.0;
    }

    // Boundary condition function: defines u(x,y) on the boundaries
    public static double g(double x, double y, double a, double b, double c, double d) {
        if (x == a || y == c)
            return 0;          // u(0, y) = 0, u(x, 0) = 0
        if (x == b)
            return 200 * y;    // u(1, y) = 200y (if b equals 1)
        if (y == d)
            return 200 * x;    // u(x, 1) = 200x (if d equals 1)
        return 0;              // Default condition
    }

    // Solve the Poisson equation using Gauss-Seidel iteration
    public static void solvePoissonEquation(double a, double b, double c, double d,
                                            int n, int m, double TOL, int N) {
        double h = (b - a) / n;  // Grid spacing in the x-direction
        double k = (d - c) / m;  // Grid spacing in the y-direction
        double lambda = (h * h) / (k * k);
        double mu = 2 * (1 + lambda);

        // Initialize mesh grid with zeros
        double[][] w = new double[n + 1][m + 1];
        for (double[] row : w) {
            Arrays.fill(row, 0.0);
        }

        // Apply boundary conditions from function g(x,y,...)
        for (int i = 0; i <= n; i++) {
            for (int j = 0; j <= m; j++) {
                double x = a + i * h;
                double y = c + j * k;
                w[i][j] = g(x, y, a, b, c, d);
            }
        }

        // Gauss-Seidel iteration loop
        int iteration = 0;
        while (iteration < N) {
            double norm = 0.0; // Maximum error across grid points in one iteration

            // Update each interior grid point (excluding the boundaries)
            for (int i = 1; i < n; i++) {
                for (int j = 1; j < m; j++) {
                    double x = a + i * h;
                    double y = c + j * k;

                    // Compute new value at grid point using finite differences and known neighboring points
                    double z = (h * h * f(x, y) + w[i - 1][j] + w[i + 1][j]
                            + lambda * (w[i][j - 1] + w[i][j + 1])) / mu;

                    // Determine the maximum difference (norm) for convergence test
                    norm = Math.max(norm, Math.abs(w[i][j] - z));

                    // Update the grid point with the new value
                    w[i][j] = z;
                }
            }

            // Check if the maximum change is within the tolerance level
            if (norm <= TOL) {
                System.out.printf("Solution grid (w[i][j]):%n");
                for (int i = 1; i < n; i++) {
                    for (int j = 1; j < m; j++) {
                        double x = a + i * h;
                        double y = c + j * k;
                        System.out.printf("(%.2f, %.2f) = %.6f%n", x, y, w[i][j]);
                    }
                }
                return; // Convergence achieved: exit the method
            }

            iteration++;
        }

        System.out.println("Maximum number of iterations exceeded.");
    }

    public static void main(String[] args) {
        // Domain boundaries (you can adjust these as needed)
        double a = 0.0, b = 0.5, c = 0.0, d = 0.5;

        // Grid dimensions: number of intervals in x and y directions
        int n = 10, m = 10;

        // Convergence tolerance and maximum number of iterations
        double TOL = 1e-10;
        int N = 10000;

        // Solve the Poisson equation using the finite difference method and Gauss-Seidel iteration
        solvePoissonEquation(a, b, c, d, n, m, TOL, N);
    }
}
