package hyperBolic;


class CrankNicolsonSolver {
    public static double[][] solveCrankNicolson(double l, double T, double alpha, int m, int N) {
        // Compute spatial step (h) and time step (k)
        double h = l / m;
        double k = T / N;
        // For the heat equation, lambda appears as: lambda = (alpha^2 * k) / (h^2)
        double lambda = alpha * alpha * k / (h * h);

        // Create the solution array for the interior points only (m-1 grid points) at N time levels.
        double[][] w = new double[m - 1][N];

        // Set the initial condition at t = 0 for all interior points.
        for (int i = 0; i < m - 1; i++) {
            w[i][0] = f((i + 1) * h);
        }

        // Set up the tridiagonal system coefficients.
        // The Crank-Nicolson discretization leads to:
        //     a u_{i-1}^j + b u_i^j + c u_{i+1}^j = RHS,
        // where a, b, and c are defined as below.
        double[] a = new double[m - 1]; // sub-diagonal coefficients
        double[] b = new double[m - 1]; // main diagonal coefficients
        double[] c = new double[m - 1]; // super-diagonal coefficients
        double[] rhs = new double[m - 1]; // right-hand side vector

        for (int i = 0; i < m - 1; i++) {
            a[i] = -lambda / 2;
            b[i] = 1 + lambda;
            c[i] = -lambda / 2;
        }

        // Time stepping loop (from time level 1 up to N-1)
        for (int j = 1; j < N; j++) {
            // Construct rhs using values from previous time step (j-1)
            // The formulation is:
            //   rhs = (λ/2)[u_{i-1}^{j-1} + u_{i+1}^{j-1}] + (1 - λ)u_i^{j-1}
            for (int i = 0; i < m - 1; i++) {
                double leftVal  = (i > 0 ? w[i - 1][j - 1] : 0);
                double rightVal = (i < m - 2 ? w[i + 1][j - 1] : 0);
                rhs[i] = lambda / 2 * (leftVal + rightVal) + (1 - lambda) * w[i][j - 1];
            }

            // Solve the tridiagonal system for the new time level.
            w = solveTridiagonal(a, b, c, rhs, w, j);
        }

        return w;
    }

    // This helper method solves the tridiagonal system using the Thomas algorithm.
    // It performs forward elimination and then back substitution to update the solution at time level j.
    private static double[][] solveTridiagonal(double[] a, double[] b, double[] c,
                                               double[] rhs, double[][] w, int j) {
        int n = rhs.length;
        double[] newC = new double[n];
        double[] newRhs = new double[n];

        // Forward elimination.
        newC[0] = c[0] / b[0];
        newRhs[0] = rhs[0] / b[0];

        for (int i = 1; i < n; i++) {
            double denom = b[i] - a[i] * newC[i - 1];
            newC[i] = c[i] / denom;
            newRhs[i] = (rhs[i] - a[i] * newRhs[i - 1]) / denom;
        }

        // Back substitution: compute the solution at time level j.
        w[n - 1][j] = newRhs[n - 1];
        for (int i = n - 2; i >= 0; i--) {
            w[i][j] = newRhs[i] - newC[i] * w[i + 1][j];
        }

        return w;
    }

    // Initial condition function u(x,0) = sin(πx)
    private static double f(double x) {
        return Math.sin(Math.PI * x);
    }
}
