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

    protected double[][] getIGammas(double[][] alpha, double[][] beta, MMC model){
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
        for (int j = 0; j < model.getA().length; j++){
            for (int i = 0; i < model.getA()[0].length; i++) {
                double numerateur = 0, denominateur = 0;
                for (int t = 0; t < Xi[0][0].length; t++) {
                    numerateur += Xi[i][j][t];
                    denominateur += iGammas[i][t];
                }
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
        boolean iterer = true;
        double ev, ev_bar;
        int iter = 0;
        do{
            alpha = FB.getAlphas(O, model);
            beta = FB.getBetas(O, model);
            Xi = getXiTab(O, alpha, beta, model);
            iGammas = getIGammas(alpha, beta, model);
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
                iterer = false;
            else {
                // Lambda := Lambda_bar
                model.setA(model_bar.getA());
                model.setB(model_bar.getB());
                model.setPI(model_bar.getPI());
            }
        }while (iterer);
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
        double ev, ev_bar;
        int iter = 0, k = 0;
        boolean fin = false;
        do {
            alphas[k] = FB.getAlphas(O[k], model);
            betas[k] = FB.getBetas(O[k], model);
            iGammas[k] = getIGammas(alphas[k], betas[k], model);
            Xis[k] = getXiTab(O[k], alphas[k], betas[k], model);
            //----------- Optimisation du vecteur PI -----------------------
            double[] Pi_bar = model_bar.getPI();
            for (int i = 0; i < Pi_bar.length; i++) {
                Pi_bar[i] = iGammas[k][i][0]/K;
            }
            //---------- Optimisation de la matrice A -----------------------
            A_bar = model_bar.getA();
            for (int i = 0; i < A_bar.length; i++) {
                for (int j = 0; j < A_bar[0].length; j++) {
                    numerateur = 0;
                    denominateur = 0;
                    for (int t = 0; t < Xis[k][i][j].length; t++) {
                        numerateur += Xis[k][i][j][t];
                        denominateur += iGammas[k][i][t];
                    }
                    A_bar[i][j] = numerateur / denominateur;
                }
            }
            //--------------- Optimisation de la matrice B -----------------------
            B_bar = model_bar.getB();
            for (int i = 0; i < B_bar.length; i++) {
                for (int j = 0; j < B_bar[0].length; j++) {
                    numerateur = 0;
                    denominateur = 0;
                    ArrayList<Integer> Ul = getLIndex(model.getSymboles()[j].trim(), O[k]);
                    if (Ul.isEmpty()) { // Si le symbole n'est pas contenu dans la séquence
                        B_bar[i][j] = 0;
                    } else {
                        for (int t = 0; t < Ul.size(); t++)
                            numerateur += iGammas[k][i][t];
                    }
                    for (int t = 0; t < iGammas[k][i].length; t++)
                        denominateur += iGammas[k][i][t];
                    B_bar[i][j] = numerateur / denominateur;
                }
            }
            ev = FB.evaluer(alphas[k], betas[k], model);
            ev_bar = FB.evaluer(alphas[k], betas[k], model_bar);
            if (ev_bar - ev <= epsilon) {
                ++k;
                iter = 0;
            }else{
                model_bar.setPI(Pi_bar);
                model_bar.setA(A_bar);
                model_bar.setB(B_bar);
                ++iter;
            }
        }while (iter < MaxIter);
        if (impression)
            model.modelOUT(model.getEspaceNom()+"_optimized.txt", model_bar);

        return model_bar;
    }

    private double[] getY1(String[] O, MMC model, double[][][] iGammas){
        final int K = O.length;
        double[] Y1 = new double[model.getPI().length];
        for (int i = 0; i < Y1.length; i++){
            double tmp = 0;
            for (int k = 0; k < K; k++){
                tmp += iGammas[k][i][0];
            }
            Y1[i] = tmp / K;
        }

        return Y1;
    }

    private double[][] getY2(String[] O, MMC model, double[][][][] Xis){
        final int K = O.length;
        double[][] Y2 = new double[model.getA().length][model.getA()[0].length];
        for (int i = 0; i < Y2.length; i++){
            for (int j = 0; j < Y2[0].length; j++){
                double tmp = 0;
                for (int k = 0; k < K; k++){
                    for (int t = 0; t < Xis[i][j].length; t++) {
                        tmp += Xis[k][i][j][t];
                    }
                }
                Y2[i][j] = tmp;
            }
        }

        return Y2;
    }

    private double[][] getY3(String[] O, MMC model, double[][][] iGammas, int parametre){
        final int K = O.length;
        double[][] Y3;
        parametre = (parametre == 1) ? parametre : 0;
        if (parametre == 0)
            Y3 = new double[model.getA().length][model.getA()[0].length];
        else
            Y3 = new double[model.getB().length][model.getB()[0].length];
        for (int i = 0; i < Y3.length; i++){
            for (int j = 0; j < Y3[0].length; j++){
                double tmp = 0;
                for (int k = 0; k < K; k++){
                    for (int t = 0; t < iGammas[k][0].length; t++)
                        tmp += iGammas[k][i][t];
                }
                Y3[i][j] = tmp;
            }
        }

        return Y3;
    }

    private double[][] getY4(String[] O, MMC model, double[][][] iGammas){
        final int K = O.length;
        double[][] Y4 = new double[model.getB().length][model.getB()[0].length];
        for (int i = 0; i < Y4.length; i++){
            for (int j = 0; j < Y4[0].length; j++){
                double tmp = 0;
                for (int k = 0; k < K; k++){
                    ArrayList<Integer> Ul = getLIndex(model.getSymboles()[j].trim(), O[k]);
                    if (Ul.isEmpty()) { // Si le symbole n'est pas contenu dans la séquence
                        Y4[i][j] = 0;
                    }else {
                        for (int t = 0; t < Ul.size(); t++)
                            tmp += iGammas[k][i][t];
                    }
                }
                Y4[i][j] = tmp;
            }
        }

        return Y4;
    }

    public MMC MultiSequence2(String[] O, MMC model, int MaxIter, double epsilon, boolean impression){ // Comme dans le cours
        final int K = O.length;
        MMC model_bar = new MMC(model);
        double[] Pi_bar;
        double[][] A_bar, B_bar, Y2, Y3, Y4;
        double[][][] alphas = new double[K][][], betas = new double[K][][], alphas_bar = new double[K][][], betas_bar = new double[K][][], iGammas = new double[K][][];
        double[][][][] Xis = new double[K][][][];
        double ev = 0, ev_bar = 0;
        int iter = 0;
        boolean iterer = true;
        do {
            for (int k = 0; k < K; k++) {
                alphas[k] = FB.getAlphas(O[k], model);
                betas[k] = FB.getBetas(O[k], model);
                iGammas[k] = getIGammas(alphas[k], betas[k], model);
                Xis[k] = getXiTab(O[k], alphas[k], betas[k], model);
            }
            //----------- Optimisation du vecteur PI -----------------------
            Pi_bar = getY1(O, model, iGammas);
            //---------- Optimisation de la matrice A -----------------------
            Y2 = getY2(O, model, Xis);
            Y3 = getY3(O, model, iGammas, 0);
            A_bar = model_bar.getA();
            for (int i = 0; i < A_bar.length; i++) {
                for (int j = 0; j < A_bar[0].length; j++)
                    A_bar[i][j] = Y2[i][j] / Y3[i][j];
            }
            //--------------- Optimisation de la matrice B -----------------------
            Y3 = getY3(O, model, iGammas, 1);
            Y4 = getY4(O, model, iGammas);
            B_bar = model_bar.getB();
            for (int i = 0; i < B_bar.length; i++) {
                for (int j = 0; j < B_bar[0].length; j++)
                    B_bar[i][j] = Y4[i][j] / Y3[i][j];
            }
            //------------- Modèle optimisé --------------------------------
            model_bar.setPI(Pi_bar);
            model_bar.setA(A_bar);
            model_bar.setB(B_bar);
            //------------ calcul des des alphas et betas du modèle optimisé -------
            for (int k = 0; k < K; k++) {
                alphas_bar[k] = FB.getAlphas(O[k], model_bar);
                betas_bar[k] = FB.getBetas(O[k], model_bar);
            }
            //--------------- Evaluation des différents modèles ------------------------
            for (int k = 0; k < K; k++){
                ev = FB.evaluer(alphas[k], betas[k], model);
                ev_bar = FB.evaluer(alphas_bar[k], betas_bar[k], model_bar);
            }
            ++iter;
            if ((ev_bar - ev <= epsilon) || iter > MaxIter)
                iterer = false;
            else {
                // Mise à jour du modèle initial
                model.setPI(model_bar.getPI());
                model.setA(model_bar.getA());
                model.setB(model_bar.getB());
            }
        }while (iterer);

        if (impression)
            model.modelOUT(model.getEspaceNom()+"_optimized.txt", model_bar);

        System.out.println(iter+" itération(s)"); //test

        return model_bar;
    }

}
