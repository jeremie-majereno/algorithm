package hyperBolic;


public class HyperbolicPartialDifferentialEquations {
    public static void main(String[] args) {
        // Input parameters
        double l = 1.0;     // Right endpoint of the spatial domain [0, l]
        double T = 1.0;     // Maximum time
        double alpha = 1.0; // Constant α (affects wave speed)
        int m = 10;         // Number of spatial divisions (m ≥ 2 implies m+1 grid points)
        int N = 100;        // Number of time divisions (N ≥ 2)

        // Solve the wave equation using a finite difference scheme
        double[][] solution = solveWaveEquation(l, T, alpha, m, N);

        // Output the results: Loop over time and space grid points and print u(x,t).
        for (int j = 0; j <= N; j++) {
            double t = j * (T / N);
            for (int i = 0; i <= m; i++) {
                double x = i * (l / m);
                System.out.printf("x=%.4f, t=%.4f, u=%.6f\n", x, t, solution[i][j]);
            }
        }
    }

    public static double[][] solveWaveEquation(double l, double T, double alpha, int m, int N) {
        // Step 1: Calculate step sizes and the stability parameter λ
        double h = l / m;         // Spatial step size
        double k = T / N;         // Time step size
        double lambda = k * alpha / h;
        double lambda2 = lambda * lambda; // Square of lambda

        // Initialize the solution grid:
        // There are m+1 points in space (including boundaries) and N+1 time levels.
        double[][] w = new double[m+1][N+1];

        // Step 2: Set boundary conditions for all time levels.
        // These are Dirichlet boundaries: u(0,t) = 0 and u(l,t) = 0.
        for (int j = 1; j <= N; j++) {
            w[0][j] = 0;  // At x = 0
            w[m][j] = 0;  // At x = l
        }

        // Step 3: Set initial conditions at t = 0 for the boundaries.
        w[0][0] = f(0);   // u(0,0) = f(0)
        w[m][0] = f(l);   // u(l,0) = f(l)

        // Step 4: Set initial conditions for interior points at t=0 and compute value at t=k.
        // For interior points: 1 <= i <= m-1.
        for (int i = 1; i < m; i++) {
            double x = i * h;
            w[i][0] = f(x); // Initial displacement: u(x,0) = f(x)
            // Compute u(x, k) using a finite difference approximation that includes initial velocity.
            w[i][1] = (1 - lambda2) * f(x) +
                    (lambda2 / 2) * (f(x + h) + f(x - h)) +
                    k * g(x);
        }

        // Step 5: Time evolution using the explicit finite difference scheme for the wave equation.
        // The recurrence relation is:
        // u_i^(j+1) = 2*(1 - λ^2)*u_i^j + λ^2*(u_(i+1)^j + u_(i-1)^j) - u_i^(j-1)
        for (int j = 1; j < N; j++) {
            for (int i = 1; i < m; i++) {
                w[i][j+1] = 2 * (1 - lambda2) * w[i][j] +
                        lambda2 * (w[i+1][j] + w[i-1][j]) -
                        w[i][j-1];
            }
        }

        return w;
    }

    // Example initial displacement function: f(x) = sin(πx)
    private static double f(double x) {
        return Math.sin(Math.PI * x);
    }

    // Example initial velocity function: g(x) = 0 (zero initial velocity)
    private static double g(double x) {
        return 0;
    }
}

