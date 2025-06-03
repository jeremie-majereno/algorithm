package runge_kunta;

public class RungeKuttaMidpoint {

    // Fonction définissant l'EDO : y' = y - t² + 1
    public static double f(double t, double y) {
        return Math.cos(t)*y;
    }

    // Méthode Runge-Kutta Midpoint
    public static void RK_midpoint(double a, double b, double alpha, int N) {
        double h = (b - a) / N;  // Taille du pas
        double t = a;
        double w = alpha;

        // Affichage des en-têtes
        System.out.println("t\t\tw (Approx. y(t))");
        System.out.println("--------------------------------------------------------");
        System.out.printf("%.6f\t\t%.6f%n", t, w);

        // Itération de la méthode RK Midpoint
        for (int i = 0; i < N; i++) {
            double K1 = h * f(t, w);
            double K2 = h * f(t + h / 2.0, w + K1 / 2.0);

            // Mise à jour de w selon la formule du point milieu
            w = w + K2;
            t = t + h;

            // Affichage du résultat courant
            System.out.printf("%.6f\t\t%.6f%n", t, w);
        }
    }

    public static void main(String[] args) {
        // Entrées
        double a = 0.0, b = 2.0*Math.PI; // Intervalle [0,2]
        double alpha = 1.0;      // Condition initiale y(0) = 5
        int N = 100;             // Nombre de pas

        // Résolution de l'EDO avec la méthode Runge-Kutta Midpoint
        RK_midpoint(a, b, alpha, N);
    }
}

