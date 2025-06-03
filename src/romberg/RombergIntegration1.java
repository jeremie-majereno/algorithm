package romberg;

import java.util.Scanner;

public class RombergIntegration1 {

    public static void main(String[] args) {
        Scanner scanner = new Scanner(System.in);

        System.out.print("Enter lower limit a: ");
        double a = scanner.nextDouble();
        System.out.print("Enter upper limit b: ");
        double b = scanner.nextDouble();
        System.out.print("Enter the number of rows n (positive integer): ");
        int n = scanner.nextInt();

        if (n <= 0) {
            System.out.println("Error: n must be a positive integer.");
            System.exit(1);
        }

        double h = b - a;

        // Initialize the first row R1 with R[1,1]
        double[] R1 = new double[1];
        R1[0] = (h / 2.0) * (f(a) + f(b));
        System.out.printf("R1: %.10f%n", R1[0]);

        for (int i = 2; i <= n; i++) {
            int m = (int) Math.pow(2, i - 2);
            double sum = 0.0;

            // Calculate the sum of midpoints
            for (int k = 1; k <= m; k++) {
                double x = a + (k - 0.5) * h;
                sum += f(x);
            }

            // Create new row R2
            double[] R2 = new double[i];
            R2[0] = 0.5 * (R1[0] + h * sum);

            // Extrapolation for each j from 2 to i
            for (int j = 2; j <= i; j++) {
                int idx = j - 1;
                double denominator = Math.pow(4, j - 1) - 1;
                R2[idx] = R2[idx - 1] + (R2[idx - 1] - R1[idx - 1]) / denominator;
            }

            // Output the current row R2
            System.out.print("R" + i + ": ");
            for (double val : R2) {
                System.out.printf("%.10f ", val);
            }
            System.out.println();

            // Update h and R1 for next iteration
            h /= 2.0;
            R1 = R2;
        }

        scanner.close();
    }


    public static double f(double x) {
        return Math.exp(x);
    }
}