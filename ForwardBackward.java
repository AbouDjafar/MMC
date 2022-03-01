public class ForwardBackward {

    protected double[][] getAlphas(String O, MMC model) {
        // Fonction de calcul et remplissage de la matrice des variables Forward
        int[] O_indice = model.get_Tab_indices(O, 0);
        double[][] alpha = new double[model.getEtats().length][O_indice.length];
        if (O_indice.length > 0) {
            for (int j = 0; j < alpha.length; j++) // Calcul de alpha[1][j]
                alpha[j][0] = model.getPI()[j] * model.getB()[j][O_indice[0]];
            for (int t = 1; t < O_indice.length; t++) { // Calcul de alpha[t+1:T][j]
                for (int j = 0; j < model.getA().length; j++) {
                    double tmp = 0;
                    for (int i = 0; i < model.getA().length; i++) {
                        tmp += alpha[i][t - 1] * model.getA()[i][j];
                    }
                    alpha[j][t] = model.getB()[j][O_indice[t]] * tmp;
                }
            }
        }
        return alpha;
    }

    //*************************** variables backward ************************************************
    protected double[][] getBetas(String O, MMC model) {
        // Fonction de calcul et remplissage de la matrice des variables Backward
        int[] O_indice = model.get_Tab_indices(O,0);
        double[][] beta = new double[model.getEtats().length][O_indice.length];
        if (O_indice.length > 0) {
            for (int j = 0; j < beta.length; j++) // Calcul de beta[T][i]
                beta[j][O_indice.length-1] = 1;
            for(int t = O_indice.length-1; t > 0; t--){ // Calcul de beta[T-1:1][i]
                for (int i = 0; i < model.getA().length; i++){
                    double tmp = 0;
                    for (int j = 0; j < model.getA().length; j++)
                        tmp += model.getB()[j][O_indice[t]]*beta[j][t]*model.getA()[i][j];
                    beta[i][t-1] = tmp;
                }
            }
        }
        return beta;
    }

    //*********************** evaluation d'un MMC selon forward-backward ******************************
    public double evaluer(double[][] alpha, double[][] beta, MMC model){
        double ev = 0.0;
        // Choix aléatoire de t barre
        int t_bar = (int) (Math.random() * alpha[0].length);
        //Calcul de Pr(O|Lambda)
        for (int i = 0; i < model.getA().length; i++)
            ev += alpha[i][t_bar]*beta[i][t_bar];

        return ev;
    }
    public double evaluer(String O, MMC model){
        double ev = 0.0;
        double[][] alpha = getAlphas(O, model);
        double[][] beta = getBetas(O, model);
        return evaluer(alpha, beta, model);
    }
}
