package hyperBolic;



public class WaveEquationComparison {
    public static void main(String[] args) {
        // Input parameters
        double l = 1.0;      // endpoint (spatial domain [0,1])
        double T = 1.0;      // maximum time
        double alpha = 2.0;  // wave speed (adjusted to match exact solution)
        int m = 100;         // space divisions (m ≥ 2; gives m+1 grid points)
        int N = 2000;        // time divisions (N ≥ 2)

        // Solve the wave equation numerically using our finite difference method
        double[][] numericalSolution = solveWaveEquation(l, T, alpha, m, N);

        // Compare with the exact solution
        System.out.println("Comparison of Numerical and Exact Solutions");
        System.out.println("==========================================");
        System.out.printf("Parameters: l=%.1f, T=%.1f, α=%.1f, m=%d, N=%d%n", l, T, alpha, m, N);
        System.out.println("Note: Exact solution is u(x,t) = sin(πx)cos(2πt)");
        System.out.println();

        // Compare at selected time points
        int[] timePoints = {0, N/4, N/2, 3*N/4, N};

        for (int j : timePoints) {
            double t = j * (T / N);
            System.out.printf("%nTime t = %.4f:%n", t);
            System.out.printf("%8s %15s %15s %15s%n", "x", "Numerical", "Exact", "Error");

            // Compare at all spatial points
            for (int i = 0; i <= m; i++) {
                double x = i * (l / m);
                double numSol = numericalSolution[i][j];
                double exactSol = exactSolution(x, t);
                double error = Math.abs(numSol - exactSol);

                System.out.printf("%8.4f %15.8f %15.8f %15.8f%n",
                        x, numSol, exactSol, error);
            }
        }

        // Calculate and display maximum error over the whole domain
        double maxError = 0;
        for (int j = 0; j <= N; j++) {
            for (int i = 0; i <= m; i++) {
                double x = i * (l / m);
                double error = Math.abs(numericalSolution[i][j] - exactSolution(x, j * (T / N)));
                if (error > maxError) {
                    maxError = error;
                }
            }
        }
        System.out.printf("%nMaximum error: %.8f%n", maxError);
    }

    // Exact solution for the wave equation with these parameters:
    // u(x,t) = sin(πx)cos(2πt)
    private static double exactSolution(double x, double t) {
        return Math.sin(Math.PI * x) * Math.cos(2 * Math.PI * t);
    }

    // Wave equation solver using a finite difference scheme.
    // The solver assumes spatial domain [0,l] and time interval [0,T].
    // It uses:
    //   - Dirichlet boundary conditions: u(0,t) = u(l,t) = 0  for all t
    //   - Initial displacement: u(x, 0) = sin(πx)
    //   - Initial velocity: u_t(x, 0) = 0
    public static double[][] solveWaveEquation(double l, double T, double alpha, int m, int N) {
        double h = l / m;   // Spatial step-size
        double k = T / N;   // Time step-size
        double lambda = k * alpha / h;
        double lambda2 = lambda * lambda;

        // w[i][j] will hold the numerical solution at spatial point i and time level j.
        // There are m+1 spatial points and N+1 time levels.
        double[][] w = new double[m+1][N+1];

        // Boundary conditions for all time levels:
        for (int j = 1; j <= N; j++) {
            w[0][j] = 0;   // u(0,t) = 0
            w[m][j] = 0;   // u(l,t) = 0
        }

        // Initial conditions at t = 0:
        w[0][0] = f(0);    // At x = 0
        w[m][0] = f(l);    // At x = l
        for (int i = 1; i < m; i++) {
            double x = i * h;
            w[i][0] = f(x); // Initial displacement u(x,0) = sin(πx)

            // For the first time level (t = k), a Taylor expansion can be used:
            // u(x,k) ≈ u(x,0) + k*u_t(x,0) + (k^2/2)*u_tt(x,0)
            // With u_t(x,0) = 0 and using a centered spatial difference for u_tt, we obtain:
            w[i][1] = (1 - lambda2) * f(x) +
                    (lambda2 / 2) * (f(x + h) + f(x - h)) +
                    k * g(x);
        }

        // Time evolution:
        // For each time level j, compute u(x, t_{j+1}) using the finite difference scheme:
        // u_i^{j+1} = 2*(1 - λ²)*u_i^j + λ²*(u_{i+1}^j + u_{i-1}^j) - u_i^{j-1}
        for (int j = 1; j < N; j++) {
            for (int i = 1; i < m; i++) {
                w[i][j+1] = 2 * (1 - lambda2) * w[i][j] +
                        lambda2 * (w[i+1][j] + w[i-1][j]) -
                        w[i][j-1];
            }
        }

        return w;
    }

    // Initial displacement function: f(x) = sin(πx)
    private static double f(double x) {
        return Math.sin(Math.PI * x);
    }

    // Initial velocity function: g(x) = 0 for t = 0 (zero initial velocity)
    private static double g(double x) {
        return 0;
    }
}
