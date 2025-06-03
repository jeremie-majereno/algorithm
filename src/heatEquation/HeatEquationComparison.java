package heatEquation;



public class HeatEquationComparison {
    public static void main(String[] args) {
        // Common parameters for both methods
        double l = 1.0;     // endpoint l
        double T = 1.0;     // maximum time T
        double alpha = 1.0; // thermal diffusivity α
        int m = 10;         // number of spatial steps (m ≥ 3), so there are m-1 interior points
        int N = 100;        // number of time steps (N ≥ 1)

        // Solve with both methods.
        // (Ensure that HeatEquationBackwardDifference and CrankNicolsonSolver are compiled and available)
        double[][] backwardSolution = HeatEquationBackwardDifference.solveHeatEquationBackwardDifference(l, T, alpha, m, N);
        double[][] crankNicolsonSolution = CrankNicolsonSolver.solveCrankNicolson(l, T, alpha, m, N);

        // List of time indices to compare (e.g., beginning, quarter, halfway, three-quarters, and final time)
        int[] timeIndicesToCompare = {0, N / 4, N / 2, 3 * N / 4, N - 1};

        System.out.println("Comparison of Heat Equation Solvers");
        System.out.println("==================================");
        System.out.printf("Parameters: l=%.1f, T=%.1f, α=%.1f, m=%d, N=%d%n", l, T, alpha, m, N);
        System.out.println();

        // Compare at selected time instants
        for (int j : timeIndicesToCompare) {
            double t = (j) * T / N;
            System.out.printf("Time t = %.4f:%n", t);
            System.out.println("  x\t\tBackward Diff\tCrank-Nicolson\tDifference");

            // Loop over all interior spatial points for this time level.
            for (int i = 0; i < m - 1; i++) {
                double x = (i + 1) * l / m;
                double backwardVal = backwardSolution[i][j];
                double cnVal = crankNicolsonSolution[i][j];
                double difference = Math.abs(backwardVal - cnVal);

                System.out.printf("  %.3f\t\t%.6f\t\t%.6f\t\t%.6f%n",
                        x, backwardVal, cnVal, difference);
            }
            System.out.println();
        }

        // Calculate and display maximum difference between the two methods overall
        double maxDifference = 0;
        for (int j = 0; j < N; j++) {
            for (int i = 0; i < m - 1; i++) {
                double diff = Math.abs(backwardSolution[i][j] - crankNicolsonSolution[i][j]);
                if (diff > maxDifference) {
                    maxDifference = diff;
                }
            }
        }
        System.out.printf("Maximum difference between methods: %.6f%n", maxDifference);
    }
}
