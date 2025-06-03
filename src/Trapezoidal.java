import java.util.function.DoubleUnaryOperator;

public class Trapezoidal {

    public static double trapezoidalRule(double a, double b, int n, DoubleUnaryOperator f) {
        if (n <= 0) {
            throw new IllegalArgumentException("Number of subintervals (n) must be greater than 0.");
        }

        double h = (b - a) / n; // Step size
        double sum = 0.5 * (f.applyAsDouble(a) + f.applyAsDouble(b)); // First and last terms

        // Loop through intermediate points
        for (int i = 1; i < n; i++) {
            double x = a + i * h;
            sum += f.applyAsDouble(x); // Add each intermediate term
        }

        return h * sum; // Multiply the sum by the step size
    }

    public static void main(String[] args) {
        try {

            double a1 = 0;
            double b1 = 10;
            int n1 = 5;
            DoubleUnaryOperator f1 = x -> Math.abs(x);
            double result1 = trapezoidalRule(a1, b1, n1, f1);
            System.out.printf("Approximation of the integral (|x|): %.6f%n", result1);

            // Example 2: f(x) = sin(x) - cos(x), a = -π, b = π
            double a2 = -Math.PI;
            double b2 = Math.PI;
            int n2 = 999;
            DoubleUnaryOperator f2 = x -> Math.sin(x) - Math.cos(x);
            double result2 = trapezoidalRule(a2, b2, n2, f2);
            System.out.printf("Approximation of the integral (sin(x) - cos(x)): %.6f%n", result2);

            // Example 3: f(x) = 4x^3, a = 0, b = 4
            double a3 = 0.0;
            double b3 = 4.0;
            int n3 = 999;
            DoubleUnaryOperator f3 = x -> 4 * Math.pow(x, 3);
            double result3 = trapezoidalRule(a3, b3, n3, f3);
            System.out.printf("Approximation of the integral (4x^3): %.6f%n", result3);

        } catch (IllegalArgumentException e) {
            System.err.println("Error: " + e.getMessage());
        }
    }
}
