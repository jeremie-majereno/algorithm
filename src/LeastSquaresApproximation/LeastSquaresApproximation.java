package LeastSquaresApproximation;



public class LeastSquaresApproximation {

    interface MyFunction {
        double apply(double x);
    }

    // Method to generate least squares polynomial coefficients
    public static double[] leastSquaresPolynomial(int degree, double a, double b, int samples, MyFunction f) {
        double[] x = new double[samples];
        double[] y = new double[samples];

// Sample points from function
        for (int i = 0; i < samples; i++) {
            x[i] = a + (b - a) * i / (samples - 1);
            y[i] = f.apply(x[i]);
        }

// Construct normal equations
        int n = degree + 1;
        double[][] X = new double[n][n]; // Coefficient matrix
        double[] Y = new double[n]; // Right-hand side

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                X[i][j] = sumPower(x, i + j);
            }
            Y[i] = sumProduct(x, y, i);
        }

        return solveLinearSystem(X, Y);
    }

    // Helper function to compute sum of powers: Σ(x^p)
    private static double sumPower(double[] x, int p) {
        double sum = 0;
        for (double xi : x) {
            sum += Math.pow(xi, p);
        }
        return sum;
    }

    // Helper function to compute Σ(y * x^p)
    private static double sumProduct(double[] x, double[] y, int p) {
        double sum = 0;
        for (int i = 0; i < x.length; i++) {
            sum += y[i] * Math.pow(x[i], p);
        }
        return sum;
    }

    // Solve linear system Ax = B using Gaussian elimination
    private static double[] solveLinearSystem(double[][] A, double[] B) {
        int n = A.length;

        for (int i = 0; i < n; i++) {
// Partial pivoting for numerical stability
            int maxRow = i;
            for (int j = i + 1; j < n; j++) {
                if (Math.abs(A[j][i]) > Math.abs(A[maxRow][i])) {
                    maxRow = j;
                }
            }
            swapRows(A, B, i, maxRow);

// Eliminate below current row
            for (int j = i + 1; j < n; j++) {
                double factor = A[j][i] / A[i][i];
                for (int k = i; k < n; k++) {
                    A[j][k] -= factor * A[i][k];
                }
                B[j] -= factor * B[i];
            }
        }

// Back-substitution
        double[] x = new double[n];
        for (int i = n - 1; i >= 0; i--) {
            x[i] = B[i];
            for (int j = i + 1; j < n; j++) {
                x[i] -= A[i][j] * x[j];
            }
            x[i] /= A[i][i];
        }
        return x;
    }

    // Swap rows for partial pivoting
    private static void swapRows(double[][] A, double[] B, int row1, int row2) {
        double[] tempA = A[row1];
        A[row1] = A[row2];
        A[row2] = tempA;

        double tempB = B[row1];
        B[row1] = B[row2];
        B[row2] = tempB;
    }

    public static void main(String[] args) {
        int degree = 15;
        double a = 0, b = Math.PI;
        int n = 20;

        double[] coefficients = leastSquaresPolynomial(degree, a, b, n, x->Math.sin(x));
        printPolynomial(coefficients);
        double error = calculateError(coefficients, x -> Math.sin(x), a, b, n);
        System.out.println("Mean Squared Error: " + String.format("%.6f", error));

    }
    public static void printPolynomial(double[] coefficients) {
        System.out.print("P(x) = ");
        for (int i = coefficients.length - 1; i >= 0; i--) {
            System.out.printf("%.6f * x^%d ", coefficients[i], i);
            if (i > 0) System.out.print("+ ");
        }
        System.out.println();
    }
    public static double calculateError(double[] coefficients, MyFunction f, double a, double b, int samples) {
        double error = 0;
        double step = (b - a) / samples;

        for (int i = 0; i < samples; i++) {
            double x = a + i * step;
            double approx = evaluatePolynomial(coefficients, x);
            error += Math.pow(f.apply(x) - approx, 2);
        }

        return error / samples;
    }

    public static double evaluatePolynomial(double[] coefficients, double x) {
        double sum = 0;
        for (int i = 0; i < coefficients.length; i++) {
            sum += coefficients[i] * Math.pow(x, i);
        }
        return sum;
    }

}
