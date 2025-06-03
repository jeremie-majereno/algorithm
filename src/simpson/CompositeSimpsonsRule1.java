package simpson;

import java.util.function.Function;

public class CompositeSimpsonsRule1 {

    // Composite Simpson's Rule implementation
    public static double compositeSimpsonRule(double a, double b, int n, Function<Double, Double> f) {
        if (n % 2 != 0) {
            throw new IllegalArgumentException("Number of subintervals (n) must be even.");
        }

        double h = (b - a) / n; // Step size
        double sum = f.apply(a) + f.apply(b); // Initialize with first and last terms

        // Loop through odd indexed points
        for (int i = 1; i < n; i += 2) {
            double x = a + i * h;
            sum += 4 * f.apply(x); // Odd terms are multiplied by 4
        }

        // Loop through even indexed points
        for (int i = 2; i < n; i += 2) {
            double x = a + i * h;
            sum += 2 * f.apply(x); // Even terms are multiplied by 2
        }

        return (h / 3) * sum; // Multiply the sum by h/3
    }

    public static void main(String[] args) {
        int n = 6;

        // Example 1: f(x) = |x|
        double a1 = 0;
        double b1 = 10;
        Function<Double, Double> f1 = x -> Math.abs(x);
        double result1 = compositeSimpsonRule(a1, b1, n, f1);
        System.out.printf("Approximation of the integral (|x|): %.6f%n", result1);


    }
}
