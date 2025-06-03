package heatEquation;


public class CrankNicolsonSolver {

    public static double[][] solveCrankNicolson(double l, double T, double alpha, int m, int N) {
        double h = l / m;       // spatial step size
        double k = T / N;       // time step size
        double lambda = alpha * alpha * k / (h * h);

        // Only interior grid points are computed: there are (m - 1) of them.
        double[][] w = new double[m - 1][N];

        // Set initial condition at t = 0 for interior nodes,
        // where x = (i + 1) * h for i = 0, 1, ..., m - 2.
        for (int i = 0; i < m - 1; i++) {
            w[i][0] = f((i + 1) * h);
        }

        // Set up the tridiagonal system coefficients for the implicit Crank-Nicolson method
        double[] a = new double[m - 1]; // sub-diagonal
        double[] b = new double[m - 1]; // main diagonal
        double[] c = new double[m - 1]; // super-diagonal
        double[] rhs = new double[m - 1];

        for (int i = 0; i < m - 1; i++) {
            a[i] = -lambda / 2;
            b[i] = 1 + lambda;
            c[i] = -lambda / 2;
        }

        // Time stepping loop: for each new time level j
        for (int j = 1; j < N; j++) {
            // Construct the right-hand side (rhs) using values from the previous time level (j-1).
            // The formulation assumes u(0,t) = u(l,t) = 0 for all time t.
            for (int i = 0; i < m - 1; i++) {
                double leftNeighbor = (i > 0) ? w[i - 1][j - 1] : 0;       // boundary at x = 0
                double rightNeighbor = (i < m - 2) ? w[i + 1][j - 1] : 0;     // boundary at x = l
                rhs[i] = lambda / 2 * (leftNeighbor + rightNeighbor) + (1 - lambda) * w[i][j - 1];
            }

            // Solve the resulting tridiagonal system to update w(:, j)
            w = solveTridiagonal(a, b, c, rhs, w, j);
        }

        return w;
    }

    // Solving the tridiagonal system using the Thomas algorithm.
    // This method performs forward elimination and then back substitution
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

        // Back substitution.
        w[n - 1][j] = newRhs[n - 1];
        for (int i = n - 2; i >= 0; i--) {
            w[i][j] = newRhs[i] - newC[i] * w[i + 1][j];
        }

        return w;
    }

    // Initial condition: u(x,0) = sin(πx)
    private static double f(double x) {
        return Math.sin(Math.PI * x);
    }

    // Utility method to print the solution matrix.
    // Each column corresponds to a time level; each row corresponds to an interior spatial node.
    public static void printSolution(double[][] w) {
        int nTimeSteps = w[0].length;
        System.out.println("Computed solution at interior points:");
        for (int j = 0; j < nTimeSteps; j++) {
            System.out.printf("Time level %d: ", j);
            for (int i = 0; i < w.length; i++) {
                System.out.printf("%.4f ", w[i][j]);
            }
            System.out.println();
        }
    }

    // Main method to test the Crank-Nicolson solver.
    public static void main(String[] args) {
        double l = 1.0;
        double T = 0.1;
        double alpha = 1.0;
        int m = 10;
        int N = 10;
        double[][] solution = solveCrankNicolson(l, T, alpha, m, N);
        printSolution(solution);
    }
}


