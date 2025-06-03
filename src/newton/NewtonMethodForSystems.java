package newton;

import java.util.Arrays;

public class NewtonMethodForSystems {

    // Example: Define the nonlinear system of equations F(x)
    // F(x) = [f1(x), f2(x), ..., fn(x)]
    public static double[] F(double[] x) {
        // Example for a 2x2 system:
        // f1(x1, x2) = x1^2 + x2^2 - 4
        // f2(x1, x2) = x1 - x2^2
        return new double[]{
                x[0] * x[0] + x[1] * x[1] - 4,  // f1(x)
                x[0] - x[1] * x[1]              // f2(x)
        };
    }

    // Define the Jacobian matrix J(x) for the system F(x)
    // J(x)_ij = ∂fi(x) / ∂xj
    public static double[][] J(double[] x) {
        // Example for a 2x2 system:
        return new double[][]{
                {2 * x[0], 2 * x[1]},  // Partial derivatives of f1
                {1, -2 * x[1]}         // Partial derivatives of f2
        };
    }

    // Solve the linear system J(x)y = -F(x) using Gaussian elimination
    // This method returns the solution vector y
    public static double[] solveLinearSystem(double[][] J, double[] F) {
        int n = F.length;
        double[] y = new double[n];
        double[][] A = new double[n][n + 1];

        // Build augmented matrix A = [J | -F]
        for (int i = 0; i < n; i++) {
            System.arraycopy(J[i], 0, A[i], 0, n);
            A[i][n] = -F[i];  // Negative of F
        }

        // Gaussian elimination
        for (int i = 0; i < n; i++) {
            // Pivoting
            int maxRow = i;
            for (int k = i + 1; k < n; k++) {
                if (Math.abs(A[k][i]) > Math.abs(A[maxRow][i])) {
                    maxRow = k;
                }
            }
            double[] temp = A[i];
            A[i] = A[maxRow];
            A[maxRow] = temp;

            // Make the diagonal element 1 and eliminate below
            for (int k = i + 1; k < n; k++) {
                double factor = A[k][i] / A[i][i];
                for (int j = i; j <= n; j++) {
                    A[k][j] -= factor * A[i][j];
                }
            }
        }

        // Back substitution
        for (int i = n - 1; i >= 0; i--) {
            y[i] = A[i][n] / A[i][i];
            for (int k = i - 1; k >= 0; k--) {
                A[k][n] -= A[k][i] * y[i];
            }
        }

        return y;
    }

    // Newton's Method for Systems
    public static void newtonMethod(double[] x0, double TOL, int N) {
        int n = x0.length;
        double[] x = Arrays.copyOf(x0, n);
        int k = 1;

        while (k <= N) {
            // Step 3: Calculate F(x) and J(x)
            double[] Fx = F(x);
            double[][] Jx = J(x);

            // Step 4: Solve J(x)y = -F(x)
            double[] y = solveLinearSystem(Jx, Fx);

            // Step 5: Update x = x + y
            for (int i = 0; i < n; i++) {
                x[i] += y[i];
            }

            // Step 6: Check if ||y|| < TOL
            double normY = 0;
            for (double yi : y) {
                normY += yi * yi;
            }
            normY = Math.sqrt(normY);

            if (normY < TOL) {
                System.out.println("The procedure was successful.");
                System.out.println("Approximate solution: " + Arrays.toString(x));
                return; // STOP
            }

            // Step 7: Increment iteration counter
            k++;
        }

        // Step 8: If maximum iterations exceeded
        System.out.println("Maximum number of iterations exceeded.");
    }

    public static void main(String[] args) {
        // Example inputs
        double[] initialGuess = {1.0, 1.0}; // Initial approximation x0
        double tolerance = 1e-6;           // Tolerance TOL
        int maxIterations = 100;           // Maximum number of iterations N

        // Run Newton's Method for Systems
        newtonMethod(initialGuess, tolerance, maxIterations);
    }
}