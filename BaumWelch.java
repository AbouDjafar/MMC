import java.util.ArrayList;

public class BaumWelch {
    private MMC model;
    private ForwardBackward FB;

    public BaumWelch(MMC model){
        this.model = new MMC(model);
        FB = new ForwardBackward();
    }

//************************ Optimisation des paramètres ***************************************************
    //---------------------- Création de la matrice des Xi -----------------------------------------
    public double[][][] getXiTab(String O, double[][] alpha, double[][] beta, MMC model){
        int[] O_i = model.get_Tab_indices(O,0);
        double ev = FB.evaluer(alpha, beta, model);
        double[][][] Xi = new double[model.getA().length][model.getA().length][O_i.length-1];
        for (int i = 0; i < Xi.length; i++){
            for (int j = 0; j < Xi[0].length; j++){
                for (int t = 0; t < Xi[0][0].length; t++){
                    Xi[i][j][t] = (alpha[i][t] * model.getA()[i][j] * model.getB()[j][O_i[t+1]] * beta[j][t+1]) / ev;
                }
            }
        }
        return Xi;
    }
    protected double[][][] getXiTab(String O){
        double[][] alpha = FB.getAlphas(O, model);
        double[][] beta = FB.getBetas(O, model);
        return getXiTab(O, alpha, beta, model);
    }

    //----------------- Optimisation de Pi --------------------------------------------
    public double[] getOptmizedPI(double[][] iGammas){
        double[] PI_barre = new double[model.getPI().length];
        for (int i = 0; i < PI_barre.length; i++) {
            PI_barre[i] = iGammas[i][0];
        }
        return PI_barre;
    }

    protected double[][] getIGammas2(double[][][] Xi){ // Comme dans le cours
        // iGamma = somme de j allant de 1 à N de Xi_t(i,j)
        double[][] GTab = new double[Xi[0].length][Xi[0][0].length];
        for (int i = 0; i < GTab.length; i++){
            for (int t = 0; t < GTab[0].length; t++){
                double tmp = 0;
                for (int j = 0; j < model.getA().length; j++)
                    tmp += Xi[i][j][t];
                GTab[i][t] = tmp;
            }
        }
        return GTab;
    }

    protected double[][] getIGammas(double[][] alpha, double[][] beta){
        // iGamma = beta_t(i) * alpha_t(i) / Pr(O|lambda)
        double[][] GTab = new double[alpha.length][alpha[0].length];
        double ev = FB.evaluer(alpha, beta, model);
        for (int t = 0; t < GTab[0].length; t++){
            for (int i = 0; i < GTab.length; i++){
                GTab[i][t] = alpha[i][t] * beta[i][t] / ev;
            }
        }
        return GTab;
    }

    // ---------------- Optimisation de la matrice des transitions d'états ------------------------
    protected double[][] getOptimizedA(double[][][] Xi, double[][] iGammas){
        // A_bar = somme t=1:T-1 de Xi_t(i,j) / somme t=1:T-1 de gamma_t(i)
        double[][] A_bar = new double[model.getA().length][model.getA().length];
        for (int i = 0; i < model.getA().length; i++){
            for (int j = 0; j < model.getA()[0].length; j++) {
                double numerateur = 0, denominateur = 0;
                for (int t = 0; t < Xi[0][0].length; t++)
                    numerateur += Xi[i][j][t];
                for (int t = 0; t < Xi[0][0].length; t++)
                    denominateur += iGammas[i][t];
                A_bar[i][j] = numerateur / denominateur;
            }
        }
        return A_bar;
    }

    //-------------------- Optimisation de la matrice des observations --------------------------------------
    protected double[][] getOptimizedB(double[][] iGammas, String O){
        double[][] B_bar = new double[model.getB().length][model.getB()[0].length];
        final int T = O.split(" ").length;
        for (int i = 0; i < B_bar.length; i++){
            for (int j = 0; j < B_bar[0].length; j++){
                ArrayList<Integer> Ul = getLIndex(model.getSymboles()[j].trim(), O);
                if (Ul.isEmpty()) { // Si le symbole n'est pas contenu dans la séquence
                    B_bar[i][j] = 0;
                }else {
                    double numerateur = 0;
                    double denominateur = 0;
                    for (int t : Ul){  // Somme_t inclu dans Ul de gamma_t(i)
                        numerateur += iGammas[i][t];
                    }
                    for (int t = 0; t < T; t++){ // Somme_t=1:T de gamma_t[i]
                        denominateur += iGammas[i][t];
                    }
                    B_bar[i][j] = numerateur / denominateur;
                }
            }
        }
        return  B_bar;
    }

    protected ArrayList<Integer> getLIndex(String l, String O){
        /*
         * Méthode pour déterminer les positions auxquelles un symbole du modèle apparait dans la séquence O
         */
        ArrayList<Integer> lIndex = new ArrayList<>();
        String[] Os = O.split(" ");
        int i = 0;
        do{
            if (l.equals(Os[i])){
                lIndex.add(i);
            }
            ++i;
        }while (i < Os.length);
        return (lIndex);
    }

//********************** Algorithmes de Baum-Welch *******************************************************
    //------------------- Baum-Welch mono séquence -----------------------------------------
    public MMC MonoSequence(String O, MMC model, int MaxIter, double epsilon, boolean impression){
        MMC model_bar = new MMC(model);
        double[][] A_bar, B_bar;
        double[] Pi_bar;
        // Calculs des variables forward et backward
        double[][] alpha, beta, iGammas, alpha_bar, beta_bar;
        double[][][] Xi;
        boolean fin = false;
        double ev, ev_bar;
        int iter = 0;
        do{
            alpha = FB.getAlphas(O, model);
            beta = FB.getBetas(O, model);
            Xi = getXiTab(O, alpha, beta, model);
            iGammas = getIGammas(alpha, beta);
            // Calcul des optimisés des éléments du modèle
            A_bar = getOptimizedA(Xi, iGammas);
            B_bar = getOptimizedB(iGammas, O);
            Pi_bar = getOptmizedPI(iGammas);
            // Recalcul des variables forward et backward mais avec les données optimisées
            model_bar.setA(A_bar);
            model_bar.setB(B_bar);
            model_bar.setPI(Pi_bar);
            alpha_bar = FB.getAlphas(O, model_bar);
            beta_bar = FB.getBetas(O, model_bar);
            // Calcul de la valeur de l'évaluation du modèle initial et du modèle optimisé
            ev = FB.evaluer(alpha, beta, model);
            ev_bar = FB.evaluer(alpha_bar, beta_bar, model_bar);
            iter++;
            if ((ev_bar - ev <= epsilon) || (iter > MaxIter))
                fin = true;
            else {
                // Lambda = Lambda_bar
                model.setA(model_bar.getA());
                model.setB(model_bar.getB());
                model.setPI(model_bar.getPI());
            }
        }while (fin);
        // test
        System.out.println("Iterations: "+iter);
        //------------ écriture sur fichier du nouveau model -----------------------------
        if (impression)
            model.modelOUT(model.getEspaceNom()+"_optimized.txt", model_bar);

        return model_bar;
    }

    //--------- Baum-Welch multi séquences --------------------------
    public MMC MultiSequence(String[] O, MMC model, int MaxIter, double epsilon, boolean impression){
        final int K = O.length;
        MMC model_bar = new MMC(model);
        double[][] A_bar, B_bar;
        double[][][] alphas = new double[K][][], betas = new double[K][][], iGammas = new double[K][][];
        double[][][][] Xis = new double[K][][][];
        double numerateur, denominateur;
        for (int k = 0; k < K; k++){
            alphas[k] = FB.getAlphas(O[k], model);
            betas[k] = FB.getBetas(O[k], model);
            iGammas[k] = getIGammas(alphas[k], betas[k]);
            Xis[k] = getXiTab(O[k], alphas[k], betas[k], model);
        }
        //----------- Optimisation du vecteur PI -----------------------
        double[] Pi_bar = model_bar.getPI();
        denominateur = K;
        for (int i = 0; i < Pi_bar.length; i++){
            numerateur = 0;
            for (int k = 0; k < K; k++){
                numerateur += iGammas[k][i][0];
            }
            Pi_bar[i] = numerateur / denominateur;
        }
        //---------- Optimisation de la matrice A -----------------------
        A_bar = model_bar.getA();
        for (int i = 0; i < A_bar.length; i++){
            for (int j = 0; j < A_bar[0].length; j++){
                numerateur = 0;
                denominateur = 0;
                for (int k = 0; k < K; k++){
                    for (int t = 0; t < Xis[i][j].length; t++) {
                        numerateur += Xis[k][i][j][t];
                        denominateur += iGammas[k][i][t];
                    }
                }
                A_bar[i][j] = numerateur / denominateur;
            }
        }
        //--------------- Optimisation de la matrice B -----------------------
        B_bar = model_bar.getB();
        for (int i = 0; i < B_bar.length; i++){
            for (int j = 0; j < B_bar[0].length; j++){
                numerateur = 0;
                denominateur = 0;
                for (int k = 0; k < K; k++){
                    ArrayList<Integer> Ul = getLIndex(model.getSymboles()[j].trim(), O[k]);
                    if (Ul.isEmpty()) { // Si le symbole n'est pas contenu dans la séquence
                        B_bar[i][j] = 0;
                    }else {
                        for (int t = 0; t < Ul.size(); t++)
                            numerateur += iGammas[k][i][t];
                    }
                    for (int t = 0; t < iGammas[k][i].length; t++)
                        denominateur += iGammas[k][i][t];
                }
                B_bar[i][j] = numerateur / denominateur;
            }
        }
        model_bar.setPI(Pi_bar);
        model_bar.setA(A_bar);
        model_bar.setB(B_bar);
        if (impression)
            model.modelOUT(model.getEspaceNom()+"_optimized.txt", model_bar);

        return model_bar;
    }

    public MMC MultiSequence2(String[] O, MMC model, int MaxIter, double epsilon, boolean impression){ // Comme dans le cours
        final int K = O.length;
        MMC model_bar = new MMC(model);
        double[][] A_bar, B_bar;
        double[][][] alphas = new double[K][][], betas = new double[K][][], iGammas = new double[K][][];
        double[][][][] Xis = new double[K][][][];
        double Y1, Y2, Y3, Y4;
        for (int k = 0; k < K; k++){
            alphas[k] = FB.getAlphas(O[k], model);
            betas[k] = FB.getBetas(O[k], model);
            iGammas[k] = getIGammas(alphas[k], betas[k]);
            Xis[k] = getXiTab(O[k], alphas[k], betas[k], model);
        }
        //----------- Optimisation du vecteur PI -----------------------
        double[] Pi_bar = model_bar.getPI();
        for (int i = 0; i < Pi_bar.length; i++){
            Y1 = 0;
            for (int k = 0; k < K; k++){
                Y1 += iGammas[k][i][0];
            }
            Pi_bar[i] = Y1 / K;
        }
        //---------- Optimisation de la matrice A -----------------------
        A_bar = model_bar.getA();
        for (int i = 0; i < A_bar.length; i++){
            for (int j = 0; j < A_bar[0].length; j++){
                Y2 = 0;
                Y3 = 0;
                for (int k = 0; k < K; k++){
                    for (int t = 0; t < Xis[i][j].length; t++) {
                        Y2 += Xis[k][i][j][t];
                        Y3 += iGammas[k][i][t];
                    }
                }
                A_bar[i][j] = Y2 / Y3;
            }
        }
        //--------------- Optimisation de la matrice B -----------------------
        B_bar = model_bar.getB();
        for (int i = 0; i < B_bar.length; i++){
            for (int j = 0; j < B_bar[0].length; j++){
                Y4 = 0;
                Y3 = 0;
                for (int k = 0; k < K; k++){
                    ArrayList<Integer> Ul = getLIndex(model.getSymboles()[j].trim(), O[k]);
                    if (Ul.isEmpty()) { // Si le symbole n'est pas contenu dans la séquence
                        B_bar[i][j] = 0;
                    }else {
                        for (int t = 0; t < Ul.size(); t++)
                            Y4 += iGammas[k][i][t];
                    }
                    for (int t = 0; t < iGammas[k][i].length; t++)
                        Y3 += iGammas[k][i][t];
                }
                B_bar[i][j] = Y4 / Y3;
            }
        }
        model_bar.setPI(Pi_bar);
        model_bar.setA(A_bar);
        model_bar.setB(B_bar);
        if (impression)
            model.modelOUT(model.getEspaceNom()+"_optimized.txt", model_bar);

        return model_bar;
    }
}
