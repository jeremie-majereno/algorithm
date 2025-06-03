package newton;

public class NewtonMethod {

    // Function f(x) - Define your specific function here
    public static double f(double x) {
        return x * x - 2; // Example: f(x) = x^2 - 2 (finding the square root of 2)
    }

    // Derivative of f(x) - Define the derivative of your specific function here
    public static double fPrime(double x) {
        return 2 * x; // Example: f'(x) = 2x
    }

    public static void newton(double p0, double TOL, int N0) {
        int i = 1; // Step 1: Initialize iteration counter

        // Step 2: Start the iteration loop
        while (i <= N0) {
            // Step 3: Compute the next approximation p
            double p = p0 - f(p0) / fPrime(p0);

            // Step 4: Check if the tolerance condition is met
            if (Math.abs(p - p0) < TOL) {
                System.out.println("The procedure was successful.");
                System.out.println("Approximate solution: " + p);
                return; // STOP
            }

            // Step 5: Increment the iteration counter
            i++;

            // Step 6: Update p0 for the next iteration
            p0 = p;
        }

        // Step 7: If the loop exits, the method failed
        System.out.println("The method failed after " + N0 + " iterations.");
    }

    public static void main(String[] args) {
        // Example inputs
        double initialGuess = 1.0; // Initial approximation p0
        double tolerance = 1e-6;  // Tolerance TOL
        int maxIterations = 100;  // Maximum number of iterations N0

        // Run Newton's method
        newton(initialGuess, tolerance, maxIterations);
    }
}