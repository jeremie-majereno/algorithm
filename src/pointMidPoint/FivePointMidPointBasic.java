package pointMidPoint;

public class FivePointMidPointBasic {
    public static void main(String[] args) {
        double x = Math.PI;
        double h = 0.01;
        double derivative = fivePointMidpointBasic(x, h);
        double exactValue = 3 * Math.cos(3 * x);//f=cos(3x)
        System.out.println("The result: " + derivative);
        System.out.println("theoretical result: " + exactValue);
        System.out.println("the error: " + Math.abs(derivative - exactValue));
    }
    public static double fivePointMidpointBasic(double x, double h) {
        double f_x_minus_2h = Math.sin(3 * (x - 2 * h));
        double f_x_minus_h = Math.sin(3 * (x - h));
        double f_x_plus_h = Math.sin(3 * (x + h));
        double f_x_plus_2h = Math.sin(3 * (x + 2 * h));
        return (f_x_minus_2h - 8 * f_x_minus_h + 8 * f_x_plus_h - f_x_plus_2h) / (12 * h);
    }
}

