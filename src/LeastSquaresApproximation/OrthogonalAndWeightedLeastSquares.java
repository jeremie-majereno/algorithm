package LeastSquaresApproximation;

public class OrthogonalAndWeightedLeastSquares {

    // Interface pour la fonction de poids.
    interface WeightFunction {
        double apply(double x);
    }

    // Méthode qui calcule la base orthogonale (à la Gram–Schmidt) sur un ensemble de points x.
    public static double[][] computeOrthogonalBasis(double[] x, int numBasis, WeightFunction w) {
        int samples = x.length;
        double[][] phi = new double[numBasis][samples];

        // Définir φ₀(x) = 1 pour tout x.
        for (int k = 0; k < samples; k++) {
            phi[0][k] = 1.0;
        }

        // Calcul de φ₁ par orthogonalisation (pour toutes les valeurs x).
        if (numBasis > 1) {
            double numerator = 0, denominator = 0;
            for (int k = 0; k < samples; k++) {
                double weight = w.apply(x[k]);
                numerator += x[k] * weight;
                denominator += weight;
            }
            double B1 = numerator / denominator;
            for (int k = 0; k < samples; k++) {
                phi[1][k] = x[k] - B1;
            }
        }

        // Pour j >= 2, on part du monôme x^j et on orthogonalise par rapport aux fonctions précédentes.
        for (int j = 2; j < numBasis; j++) {
            // Initialisation avec le monôme x^j.
            for (int k = 0; k < samples; k++) {
                phi[j][k] = Math.pow(x[k], j);
            }
            // Orthogonalisation par rapport à φ₀, φ₁, ..., φ_{j-1}.
            for (int i = 0; i < j; i++) {
                double projNumerator = 0, projDenom = 0;
                for (int k = 0; k < samples; k++) {
                    double weight = w.apply(x[k]);
                    projNumerator += weight * phi[i][k] * phi[j][k];
                    projDenom += weight * phi[i][k] * phi[i][k];
                }
                double projCoeff = 0.0;
                if (Math.abs(projDenom) > 1e-10) {
                    projCoeff = projNumerator / projDenom;
                }
                for (int k = 0; k < samples; k++) {
                    phi[j][k] -= projCoeff * phi[i][k];
                }
            }
        }
        return phi;
    }

    /**
     * Pour une série de points dans [a, b] (numSamples points),
     * cette méthode calcule la base orthogonale (6 fonctions : φ₀ à φ₅)
     * puis pour chaque point affiche la valeur de chaque φ_i(x)
     * et le terme correspondant dans l'expression :
     *
     * P(x) = -0,002433*φ₅(x) + 0,018929*φ₄(x) + 0,037309*φ₃(x)
     *        - 0,135616*φ₂(x) - 0,080596*φ₁(x) - 0,197534*φ₀(x)
     */
    public static void printCustomPolynomialValues(double a, double b, int numSamples, WeightFunction w) {
        // Création du tableau d'échantillons.
        double[] xSamples = new double[numSamples];
        for (int i = 0; i < numSamples; i++) {
            xSamples[i] = a + (b - a) * i / (numSamples - 1);
        }
        // Calcul de la base orthogonale sur l'ensemble des points.
        double[][] phi = computeOrthogonalBasis(xSamples, 6, w);
        // Les coefficients, donnés dans l'ordre naturel :
        // coeff[0] pour φ₀, coeff[1] pour φ₁, …, coeff[5] pour φ₅.
        double[] coeff = new double[]{-0.197534, -0.080596, -0.135616, 0.037309, 0.018929, -0.002433};

        // Pour chaque valeur de x, on affiche :
        // - La valeur de chaque φ_i(x)
        // - Le calcul du terme coeff[i] * φ_i(x)
        // - La somme finale P(x)
        for (int i = 0; i < numSamples; i++) {
            System.out.printf("Pour x = %.4f:%n", xSamples[i]);
            double total = 0.0;
            // On affiche dans l'ordre de φ₅ jusqu'à φ₀, pour coller à l'expression donnée.
            for (int j = 5; j >= 0; j--) {
                double term = coeff[j] * phi[j][i];
                total += term;
                System.out.printf("  φ_%d(x) = %.6f   =>   %+.6f * φ_%d(x) = %+.6f%n",
                        j, phi[j][i], coeff[j], j, term);
            }
            System.out.printf("==> P(%.4f) = %.6f%n%n", xSamples[i], total);
        }
    }

    public static void main(String[] args) {
        // Choix de l'intervalle d'évaluation.
        double a = -3, b = 4;
        int numSamples = 20;

        // Définition de la fonction de poids.
        // Pour que w(x) = 1/sqrt(25 - x²) soit bien définie, il faut |x| < 5.
        WeightFunction w = x -> 1.0 / Math.sqrt(25 - x * x);

        // Affichage détaillé du calcul du polynôme personnalisé pour chaque x.
        printCustomPolynomialValues(a, b, numSamples, w);
    }
}
