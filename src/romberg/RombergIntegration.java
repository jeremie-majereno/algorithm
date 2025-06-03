package romberg;

import java.util.Scanner;

public class RombergIntegration {

    // Define the function to integrate
    public static double f(double x) {
        // Example function: f(x) = x^2; modify this as needed
        return x * x; // Change this to your desired function
    }

    public static void main(String[] args) {
        Scanner scanner = new Scanner(System.in);

        // Input endpoints a and b, and the integer n
        System.out.print("Enter the value of a: ");
        double a = scanner.nextDouble();

        System.out.print("Enter the value of b: ");
        double b = scanner.nextDouble();

        System.out.print("Enter the integer (n > 0): ");
        int n = scanner.nextInt();

        // Validate input for n
        if (n <= 0) {
            System.out.println("Error: n must be greater than 0.");
            return;
        }

        // Array to store results, only 2 rows needed
        double[][] R = new double[2][n + 1];

        // Step 1: Compute initial h and R[1][1]
        double h = b - a;
        R[0][1] = h / 2 * (f(a) + f(b));

        // Step 2: Output R[1][1]
        System.out.printf("R[%d][%d] = %.6f%n", 1, 1, R[0][1]);

        // Step 3: Start iteration from i = 2 to n
        for (int i = 2; i <= n; i++) {
            // Step 4: Compute R[2][1] (Trapezoidal approximation)
            double midSum = 0;
            for (int k = 1; k <= (1 << (i - 2)); k++) { // 2^(i-2) new points
                midSum += f(a + (k - 0.5) * h / (1 << (i - 1))); // Points at (k-0.5)*h/2^(i-1)
            }
            R[1][1] = 0.5 * (R[0][1] + h / (1 << (i - 1)) * midSum);

            // Step 5: Update R[2][j] for j = 2 to i (Extrapolation)
            for (int j = 2; j <= i; j++) {
                R[1][j] = R[1][j - 1] + (R[1][j - 1] - R[0][j - 1]) / (4 * j - 1);
            }

            // Step 6: Output R[2][j] for j = 1, 2, ..., i
            for (int j = 1; j <= i; j++) {
                System.out.printf("R[%d][%d] = %.6f%n", 2, j, R[1][j]);
            }

            // Step 7: Update h
            h /= 2;

            // Step 8: Update row 1 of R
            for (int j = 1; j <= i; j++) {
                R[0][j] = R[1][j];
            }
        }

        // Step 9: End of program
        System.out.println("Computation completed.");
        scanner.close();
    }
}