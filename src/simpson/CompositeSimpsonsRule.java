package simpson;

import java.util.function.DoubleUnaryOperator;

public class CompositeSimpsonsRule {
    public static void main(String[] args) {
        DoubleUnaryOperator f = x -> Math.pow(64, 5/6);
        double a = 0;
        double b = 10;
        int n = 6;
        double result = compositeSimpson(a, b, n, f);
        System.out.printf("Approximated integral of f(x) = x^(5/6) from %f to %f: %f%n", a, b, result);
    }

    public static double compositeSimpson(double a, double b, int n, DoubleUnaryOperator f) {
        if (n <= 0 || n % 2 != 0) {
            throw new IllegalArgumentException("n must be a positive even integer");
        }
        double h = (b - a) / n;
        double XI0 = f.applyAsDouble(a) + f.applyAsDouble(b);
        double XI1 = 0.0;
        double XI2 = 0.0;
        for (int i = 1; i < n; i++) {
            double x = a + i * h;
            if (i % 2 == 0) {
                XI2 += f.applyAsDouble(x);
            } else {
                XI1 += f.applyAsDouble(x);
            }
        }
        double XI = ((h) * (XI0 + 2 * XI2 + 4 * XI1))/3;
        return XI;
    }
}
