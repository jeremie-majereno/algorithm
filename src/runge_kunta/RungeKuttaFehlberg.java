package runge_kunta;


public class RungeKuttaFehlberg {

    public static double f(double t, double y) {
        return Math.cos(t) * y;
    }


    // Méthode Runge-Kutta-Fehlberg (RKF45)
    public static void RKF45(double a, double b, double alpha, double TOL, double hmax, double hmin) {
        double t = a;
        double w = alpha;
        double h = hmax;
        boolean flag = true;

        // Affichage des valeurs initiales
        System.out.println("t\t\tw (Approx. y(t))\tStep size (h)");
        System.out.println("--------------------------------------------------------");
        System.out.printf("%.6f\t\t%.6f\t\t%.6f%n", t, w, h);

        while (flag) {
            // Calcul des coefficients de Runge-Kutta
            double K1 = h * f(t, w);
            double K2 = h * f(t + (1.0 / 4.0) * h, w + (1.0 / 4.0) * K1);
            double K3 = h * f(t + (3.0 / 8.0) * h, w + (3.0 / 32.0) * K1 + (9.0 / 32.0) * K2);
            double K4 = h * f(t + (12.0 / 13.0) * h, w + (1932.0 / 2197.0) * K1 - (7200.0 / 2197.0) * K2 + (7296.0 / 2197.0) * K3);
            double K5 = h * f(t + h, w + (439.0 / 216.0) * K1 - 8.0 * K2 + (3680.0 / 513.0) * K3 - (845.0 / 4104.0) * K4);
            double K6 = h * f(t + (1.0 / 2.0) * h, w - (8.0 / 27.0) * K1 + 2.0 * K2 - (3544.0 / 2565.0) * K3
                    + (1859.0 / 4104.0) * K4 - (11.0 / 40.0) * K5);

            // Calcul de l'estimation de l'erreur
            double R = (1.0 / h) * Math.abs((1.0 / 360.0) * K1 - (128.0 / 4275.0) * K3
                    - (2197.0 / 75240.0) * K4 + (1.0 / 50.0) * K5 + (2.0 / 55.0) * K6);

            if (R <= TOL) { // Si l'erreur est acceptable, on met à jour t et w
                t += h;
                w += (25.0 / 216.0) * K1 + (1408.0 / 2565.0) * K3 + (2197.0 / 4104.0) * K4 - (1.0 / 5.0) * K5;
                System.out.printf("%.6f\t\t%.6f\t\t%.6f%n", t, w, h);
            }

            // Calcul de la nouvelle taille de pas
            double delta = 0.84 * Math.pow(TOL / R, 0.25);
            if (delta <= 0.1) {
                h = 0.1 * h;
            } else if (delta >= 4) {
                h = 4 * h;
            } else {
                h = delta * h;
            }

            // S'assurer que la taille de pas reste dans les limites
            if (h > hmax) {
                h = hmax;
            }

            if (t >= b) {
                flag = false; // Fin de l'intégration si le bout de l'intervalle est atteint
            } else if (t + h > b) {
                h = b - t; // Ajuster le pas final pour atteindre exactement b
            } else if (h < hmin) {
                flag = false;
                System.out.println("Minimum step size exceeded.");
            }
        }
    }

    public static void main(String[] args) {
        // Entrées

        double a = 0.0, b = 2.0 * Math.PI;   // Intervalle [0,2]
        double alpha = 1.0;        // Condition initiale y(0) = 1.0
        double TOL = 1e-6;         // Tolérance pour le contrôle de l'erreur
        double hmax = 0.1;         // Taille de pas maximale
        double hmin = 1e-5;        // Taille de pas minimale

        // Résoudre en utilisant la méthode Runge-Kutta-Fehlberg
        RKF45(a, b, alpha, TOL, hmax, hmin);
    }
}
