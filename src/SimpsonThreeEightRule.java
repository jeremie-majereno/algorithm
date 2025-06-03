import java.util.Scanner;
import java.util.function.DoubleUnaryOperator;

public class SimpsonThreeEightRule {

    // Function to calculate integral using Simpson's 3/8 Rule
    public static double simpsonThreeEightRule(double a, double b, int n, DoubleUnaryOperator f) {
        if (n % 3 != 0) {
            throw new IllegalArgumentException("Number of subintervals (n) must be divisible by 3.");
        }

        double h = (b - a) / n;  // step size
        double sum = f.applyAsDouble(a) + f.applyAsDouble(b);  // first and last terms

        // Loop through intermediate points and apply weighting
        for (int i = 1; i < n; i++) {
            double x = a + i * h;
            if (i % 3 == 0) {
                sum += 2 * f.applyAsDouble(x);
            } else {
                sum += 3 * f.applyAsDouble(x);
            }
        }

        return (3 * h / 8) * sum;
    }

    public static void main(String[] args) {
        try {
            // Example 1: f(x) = |x|, a = 0, b = 10, n = 6

            System.out.println("Exemple 1");
            Scanner sc =new Scanner(System.in);
            System.out.println("enter a");
            double a1 = sc.nextDouble();
            System.out.println("enter b");
            double b1 = sc.nextDouble();
            System.out.println("enter n");
            int n1 = sc.nextInt();  // Ensure n is divisible by 3
            DoubleUnaryOperator f1 = x -> Math.pow(64,5.0/6.0);
            double result1 = simpsonThreeEightRule(a1, b1, n1, f1);
            System.out.printf("Approximation of the integral (64^(5/6): %.20f%n", result1);

            // Example 2: f(x) = sin(x) - cos(x), a = -π, b = π, n = 6
            double a2 = -Math.PI;
            double b2 = Math.PI;
            int n2 = 6;  // Ensure n is divisible by 3
            DoubleUnaryOperator f2 = x -> Math.sin(x) - Math.cos(x);
            double result2 = simpsonThreeEightRule(a2, b2, n2, f2);
            System.out.println("\nExemple 2");
            System.out.println("a="+ a2);
            System.out.println("b="+ b2);
            System.out.println("n="+ n2);
            System.out.printf("Approximation of the integral (sin(x) - cos(x)): %.20f%n", result2);

            // Example 3: f(x) = 4x^3, a = 0, b = 4, n = 6
       //     Scanner sc = new Scanner(System.in);
            double a3 = 0.0;
            double b3 = 4.0;
            int n3 = 6;  // Ensure n is divisible by 3
            DoubleUnaryOperator f3 = x -> 4 * Math.pow(x, 3);
            double result3 = simpsonThreeEightRule(a3, b3, n3, f3);
            System.out.println("\nExemple 3");
            System.out.println("a="+ a3);
            System.out.println("b="+ b3);
            System.out.println("n="+ n3);
            System.out.printf("Approximation of the integral (4x^3): %.20f%n", result3);
        } catch (IllegalArgumentException e) {
            System.err.println("Error: " + e.getMessage());
        }
    }
}
