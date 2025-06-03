package runge_kunta;

public class RungeKuttaMethod {

    // Fonction représentant l'EDO : dy/dt = 2*t*y
    public static double f(double t, double y) {
        return Math.cos(t)*y;
    }

    // Solution analytique de l'équation différentielle : y(t) = exp(t²) y(t) = exp(sin(t))
    public static double analyticalSolution(double t) {
        return Math.exp(Math.sin(t));
    }

    // Méthode Runge-Kutta d'ordre 4
    public static double[] rungeKutta(double alpha, double h, double t0, double tEnd) {
        // Nombre d'étapes
        int n = (int) ((tEnd - t0) / h);
        double[] w = new double[n + 1];
        double[] tArr = new double[n + 1];

        // Conditions initiales
        w[0] = alpha;
        tArr[0] = t0;

        for (int i = 0; i < n; i++) {
            tArr[i + 1] = tArr[i] + h;

            // Calcul des coefficients Runge-Kutta
            double k1 = h * f(tArr[i], w[i]);
            double k2 = h * f(tArr[i] + h / 2.0, w[i] + k1 / 2.0);
            double k3 = h * f(tArr[i] + h / 2.0, w[i] + k2 / 2.0);
            double k4 = h * f(tArr[i] + h, w[i] + k3);

            // Mise à jour de w[i+1]
            w[i + 1] = w[i] + (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;
        }

        return w; // Renvoie le tableau des valeurs approchées
    }

    public static void main(String[] args) {
        // Conditions initiales
        double alpha = 1.0;  // y(0) = 1
        double h = 0.01;     // Pas de temps
        double t0 = 0.0;     // Temps initial
        double tEnd = 2.00*Math.PI;   // Temps final

        double[] numericalResults = rungeKutta(alpha, h, t0, tEnd);


        System.out.println("t      Numerical(w)   Analytical(y)   Error");
        for (int i = 0; i < numericalResults.length; i++) {
            double t = t0 + i * h;
            double analytical = analyticalSolution(t);
            double error = Math.abs(analytical - numericalResults[i]);
            System.out.printf("%6.2f    %10.5f    %10.5f    %e%n", t, numericalResults[i], analytical, error);
        }
    }
}

