package heatEquation;



public class HeatEquationBackwardDifference {

    public static void main(String[] args) {
        // Input parameters
        double l = 1.0;       // Length of spatial domain (endpoint l)
        double T = 1.0;       // Maximum time T
        double alpha = 1.0;   // Thermal diffusivity α
        int m = 10;           // Number of spatial steps (m ≥ 3)
        int N = 100;          // Number of time steps (N ≥ 1)

        // Solve the heat equation using the Backward Difference Method
        double[][] solution = solveHeatEquationBackwardDifference(l, T, alpha, m, N);

        // Output the solution: for each time step, output the interior spatial points
        for (int j = 0; j < N; j++) {
            double t = (j + 1) * T / N;
            System.out.printf("Time t = %.4f:%n", t);
            for (int i = 0; i < m - 1; i++) {
                double x = (i + 1) * l / m;
                System.out.printf("   x = %.4f, u(x,t) ≈ %.6f%n", x, solution[i][j]);
            }
        }
    }

    public static double[][] solveHeatEquationBackwardDifference(double l, double T, double alpha, int m, int N) {
        // Step 1: Initialize parameters for spatial and temporal discretization
        double h = l / m;         // Spatial step size
        double k = T / N;         // Time step size
        double lambda = Math.pow(alpha, 2) * k / Math.pow(h, 2);

        // Allocate solution array (for interior spatial points only: m-1 points, for N time steps)
        double[][] w = new double[m - 1][N];

        // Step 2: Set the initial values (using f(x) as the initial condition)
        for (int i = 0; i < m - 1; i++) {
            double x = (i + 1) * h;
            w[i][0] = f(x);
        }

        // Steps 3-5: Set up the coefficients of the tridiagonal system for the implicit method.
        // The system arises from the backward difference discretization of
        // (u^{j} - u^{j-1})/k = α² * uₓₓ^{j}.
        // This leads to, for each interior point i:
        //   -λ * u_{i-1}^{j} + (1 + 2λ) * u_{i}^{j} - λ * u_{i+1}^{j} = u_{i}^{j-1}.
        // We set:
        //   a[i] = -λ, b[i] = 1 + 2λ, c[i] = -λ.
        // Here we precompute new diagonal coefficients (lCoeff) and upper diagonal coefficients (uCoeff)
        // for the Thomas algorithm.
        double[] lCoeff = new double[m - 1];  // Modified main diagonal coefficients
        double[] uCoeff = new double[m - 2];  // Modified upper diagonal coefficients

        // For the first interior point:
        lCoeff[0] = 1 + 2 * lambda;
        uCoeff[0] = -lambda / lCoeff[0];

        // For interior points 1 to m-2 (indices 1 to m-3)
        for (int i = 1; i < m - 2; i++) {
            lCoeff[i] = 1 + 2 * lambda + lambda * uCoeff[i - 1];
            uCoeff[i] = -lambda / lCoeff[i];
        }
        // Last interior point (index m-2)
        lCoeff[m - 2] = 1 + 2 * lambda + lambda * uCoeff[m - 3];

        // Steps 6-11: Time-stepping using the implicit backward difference method.
        for (int j = 1; j < N; j++) {
            // Step 7: Forward substitution to solve Lz = w(:, j-1)
            double[] z = new double[m - 1];
            z[0] = w[0][j - 1] / lCoeff[0];
            // Step 8: Proceed with forward substitution for the remaining points.
            for (int i = 1; i < m - 1; i++) {
                z[i] = (w[i][j - 1] + lambda * z[i - 1]) / lCoeff[i];
            }

            // Steps 9-10: Backward substitution to solve U * w^{j} = z.
            w[m - 2][j] = z[m - 2];
            for (int i = m - 3; i >= 0; i--) {
                w[i][j] = z[i] - uCoeff[i] * w[i + 1][j];
            }
        }

        return w;
    }

    // Example initial condition function f(x).
    // Here u(x, 0) = sin(πx), which satisfies u(0,0) = 0 and u(l,0) = 0 when l = 1.
    private static double f(double x) {
        return Math.sin(Math.PI * x);
    }
}

