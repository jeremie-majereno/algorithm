package LeastSquaresApproximation;

import java.util.function.DoubleUnaryOperator;

public class LeastSquaresApproximation00 {

    public static void main(String[] args) {
        int degree = 5; // Degree of the polynomial basis
        double a = -5, b = 5; // Interval bounds
        DoubleUnaryOperator w = x -> 1.0 / Math.sqrt(100 - x * x); // Weight function
        DoubleUnaryOperator f = x -> Math.sin(x) - Math.cos(x); // Target function

        // Compute orthogonal polynomials basis
        double[][] basis = computeOrthogonalBasis(degree, a, b, w);

        // Display the orthogonal basis
        displayOrthogonalBasis(basis, degree);

        // Perform least squares polynomial approximation
        double[] coefficients = leastSquaresApproximation(degree, basis, a, b, f, w);

        // Print the resulting polynomial
        System.out.println("Polynomial approximation:");
        printPolynomial(coefficients);

        // Calculate and display the mean squared error
        double mse = calculateMeanSquaredError(coefficients, f, a, b, 100);
        System.out.println("Mean squared error: " + mse);
    }

    // ----- Orthogonal Basis Computation -----
    public static double[][] computeOrthogonalBasis(int degree, double a, double b, DoubleUnaryOperator w) {
        double[][] phi = new double[degree + 1][100]; // Store φₖ(x)
        double[] x = generateSamples(a, b, 100); // Generate x sample points


        for (int i = 0; i < x.length; i++) {
            phi[0][i] = 1.0;
        }

        // Recursively calculate φ₁(x) and φₖ(x) for k ≥ 2
        for (int k = 1; k <= degree; k++) {
            double Bk = computeBk(k, x, phi, w);
            double Ck = (k >= 2) ? computeCk(k, x, phi, w) : 0;

            for (int i = 0; i < x.length; i++) {
                phi[k][i] = (x[i] - Bk) * phi[k - 1][i] - (k >= 2 ? Ck * phi[k - 2][i] : 0);
            }
        }

        return phi;
    }

    // ----- Display Orthogonal Basis -----
    public static void displayOrthogonalBasis(double[][] basis, int degree) {
        System.out.println("Base des polynômes orthogonaux :");
        for (int k = 0; k <= degree; k++) {
            System.out.print("φ" + k + "(x) = ");
            for (int i = 0; i < basis[k].length; i++) {
                System.out.printf("%.5f", basis[k][i]);
                if (i < basis[k].length - 1) {
                    System.out.print(", ");
                }
            }
            System.out.println();
        }
    }

    // ----- Least Squares Approximation -----
    public static double[] leastSquaresApproximation(int degree, double[][] basis, double a, double b,
                                                     DoubleUnaryOperator f, DoubleUnaryOperator w) {
        double[] coefficients = new double[degree + 1];
        double[] x = generateSamples(a, b, 100);

        for (int k = 0; k <= degree; k++) {
            double numerator = 0;
            double denominator = 0;

            for (int i = 0; i < x.length; i++) {
                double weight = w.applyAsDouble(x[i]);
                numerator += weight * f.applyAsDouble(x[i]) * basis[k][i];
                denominator += weight * Math.pow(basis[k][i], 2);
            }
            coefficients[k] = numerator / denominator;
        }

        return coefficients;
    }

    // ----- Error Calculation -----
    public static double calculateMeanSquaredError(double[] coefficients, DoubleUnaryOperator f, double a, double b, int samples) {
        double[] x = generateSamples(a, b, samples);
        double mse = 0;

        for (double xi : x) {
            double approximation = evaluatePolynomial(coefficients, xi);
            double error = f.applyAsDouble(xi) - approximation;
            mse += Math.pow(error, 2);
        }

        return mse / samples;
    }

    // ----- Helper Methods -----
    private static double computeBk(int k, double[] x, double[][] phi, DoubleUnaryOperator w) {
        double numerator = 0;
        double denominator = 0;

        for (int i = 0; i < x.length; i++) {
            double weight = w.applyAsDouble(x[i]);
            numerator += x[i] * weight * Math.pow(phi[k - 1][i], 2);
            denominator += weight * Math.pow(phi[k - 1][i], 2);
        }

        return numerator / denominator;
    }

    private static double computeCk(int k, double[] x, double[][] phi, DoubleUnaryOperator w) {
        double numerator = 0;
        double denominator = 0;

        for (int i = 0; i < x.length; i++) {
            double weight = w.applyAsDouble(x[i]);
            numerator += weight * phi[k - 1][i] * phi[k - 2][i];
            denominator += weight * Math.pow(phi[k - 2][i], 2);
        }

        return numerator / denominator;
    }

    private static double[] generateSamples(double a, double b, int samples) {
        double[] x = new double[samples];
        for (int i = 0; i < samples; i++) {
            x[i] = a + i * (b - a) / (samples - 1);
        }
        return x;
    }

    private static double evaluatePolynomial(double[] coefficients, double x) {
        double sum = 0;
        for (int i = 0; i < coefficients.length; i++) {
            sum += coefficients[i] * Math.pow(x, i);
        }
        return sum;
    }

    private static void printPolynomial(double[] coefficients) {
        StringBuilder polynomial = new StringBuilder("P(x) = ");
        for (int i = 0; i < coefficients.length; i++) {
            polynomial.append(coefficients[i]).append(" * x^").append(i);
            if (i < coefficients.length - 1) {
                polynomial.append(" + ");
            }
        }
        System.out.println(polynomial);
    }
}
