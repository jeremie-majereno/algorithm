import java.util.Scanner;

public class Main {

    // Example function: f(x) = x^2 - 2
    public static double f(double x) {
        return x * x - 2;
    }

    // Derivative of the example function: f'(x) = 2x
    public static double fPrime(double x) {
        return 2 * x;
    }

    public static void newtonRaphson(double p0, double TOL, int N0) {
        int i = 1;
        while (i <= N0) {
            double fp0 = f(p0);
            double fprimep0 = fPrime(p0);

            // Check if derivative is zero to avoid division by zero
            if (Math.abs(fprimep0) < 1e-12) {
                System.out.println("Derivative is zero. Method aborted.");
                return;
            }

            double p = p0 - fp0 / fprimep0;
            System.out.printf("Iteration %d: p = %.6f%n", i, p);

            if (Math.abs(p - p0) < TOL) {
                System.out.printf("Approximate solution: %.6f%n", p);
                return;
            }

            p0 = p;
            i++;
        }
        System.out.println("Method failed after " + N0 + " iterations.");
    }

    public static void main(String[] args) {
        Scanner scanner = new Scanner(System.in);

        System.out.print("Enter initial approximation p0: ");
        double p0 = scanner.nextDouble();

        System.out.print("Enter tolerance TOL: ");
        double TOL = scanner.nextDouble();

        System.out.print("Enter maximum number of iterations N0: ");
        int N0 = scanner.nextInt();

        newtonRaphson(p0, TOL, N0);

        scanner.close();
    }
}