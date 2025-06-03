package pointMidPoint;

public class ThreePointMidPoint {
    public class ThreePointMidpointDerivative {
        public static void main(String[] args) {
            double x = Math.PI;
            double h = 0.0001;
            double derivative = threePointMidpointDerivative(x, h);
            double exactValue = 3 * Math.cos(3 * x);
            System.out.println("The result: " + derivative);
            System.out.println("theoretical result: " + exactValue);
            System.out.println("the error: " + Math.abs(derivative - exactValue));
        }

        public static double threePointMidpointDerivative(double x, double h) {
            // f(x) = sin(3x)
            double fx_plus_h = Math.sin(3 * (x + h));
            double fx_minus_h = Math.sin(3 * (x - h));
            // [f(x+h) - f(x-h)] / (2h)
            return (fx_plus_h - fx_minus_h) / (2 * h);
        }
    }
}
