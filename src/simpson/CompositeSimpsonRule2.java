package simpson;

import java.util.Scanner;
import java.util.function.Function;

public class CompositeSimpsonRule2 {

    // Function to calculate the integral using Composite Simpson's Rule
    public static double compositeSimpsonsRule(Function<Double, Double> f, double a, double b, int n) {
        double h = (b - a) / n; // Step size
        double XI0 = f.apply(a) + f.apply(b); // Initial sum with endpoints
        double XI1 = 0.0; // Sum of odd-indexed terms
        double XI2 = 0.0; // Sum of even-indexed terms

        // Loop through intervals
        for (int i = 1; i < n; i++) {
            double X = a + i * h; // Compute the current point
            if (i % 2 == 0) {
                XI2 += f.apply(X); // Add to even-indexed sum
            } else {
                XI1 += f.apply(X); // Add to odd-indexed sum
            }
        }

        // Final calculation
        return h * (XI0 + 2 * XI2 + 4 * XI1) / 3.0;
    }

    public static void main(String[] args) {
        Scanner scanner = new Scanner(System.in);

        // Input limits and subintervals
        System.out.print("Enter the lower limit (a): ");
        double a = scanner.nextDouble();
        System.out.print("Enter the upper limit (b): ");
        double b = scanner.nextDouble();
        System.out.print("Enter the number of subintervals (even positive integer, n): ");
        int n = scanner.nextInt();

        if (n <= 0 || n % 2 != 0) {
            System.err.println("Error: n must be an even positive integer.");
            System.exit(1);
        }

        // Define the function f(x) = 64^(5/6)
        Function<Double, Double> f = x -> Math.pow(64, 5.0 / 6.0);

        // Calculate the integral approximation
        double result = compositeSimpsonsRule(f, a, b, n);

        // Output the result
        System.out.printf("The approximate integral of f(x) = 64^(5/6) from %.2f to %.2f using %d subintervals is: %.6f%n", a, b, n, result);

        scanner.close();
    }
}
