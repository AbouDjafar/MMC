// @author: Djafar
/*
    /!\ Formatage des fichier d'entrée des MMC /!\
    ligne 1: nombre total des états (N)
    ligne 2: nombre total des symboles (T)
    ligne 3: lister les états, séparés par un espace (e1 e2 e3 e4 ... eN)
    ligne 4: lister les symboles, séparés par un espace (s1 s2 s3 s4 ... sT)
    ligne 5: dérouler la matrice de transition d'états (A11 A12 .. A1N A21 A22 ... A2N A31 A32 ... A3N ... AN1 AN2 ... ANN)
    ligne 6: dérouler la matrice d'observations (B11 B12 ...B1T B21 B22 ... B2T ... BN1 BN2 ... BNT)
    ligne 7: lister le vecteur des états initiaux (pi1 pi2 pi3 ... piN)

    /!\ Formatage des séquences d'entrée /!\
    sur une seule ligne: s1 s2 s3 ... sT
*/
import java.io.*;
import java.util.ArrayList;

public class MMC {
    private String[] Etats;
    private String[] Symboles;
    private double[][] A;
    private double[][] B;
    private double[] PI;

    public double[][] getA() {
        return A;
    }

    public void setA(double[][] a) {
        A = a;
    }

    public double[][] getB() {
        return B;
    }

    public void setB(double[][] b) {
        B = b;
    }

    public double[] getPI() {
        return PI;
    }

    public void setPI(double[] PI) {
        this.PI = PI;
    }

    public MMC(String fileURL){
        File model = new File(fileURL);
        FileReader fr;
        StringBuilder s = new StringBuilder();
        int k;
//************** Lecture du fichier source **************************************
        try {
            fr = new FileReader(model);
            while((k = fr.read()) != -1){
                s.append((char) k);
            }
            fr.close();
        } catch (IOException e){
            e.printStackTrace();
        }
//*************** Allocation des dimensions aux différents tableaux ******************
        String[] sTab = s.toString().split("\n");
        Etats = new String[Integer.parseInt(sTab[0].trim())];
        Symboles = new String[Integer.parseInt(sTab[1].trim())];
        A = new double[Etats.length][Etats.length];
        B = new double[Etats.length][Symboles.length];
        PI = new double[Etats.length];
//**************** Remplissage des tableaux alloués ************************************
        Etats = sTab[2].split(" ");
        Symboles = sTab[3].split(" ");
        String[] tmp = sTab[4].split(" ");
        k = 0;
        for (int i = 0; i < A.length; i++){
            for (int j = 0; j < A[i].length; j++) {
                A[i][j] = Double.parseDouble(tmp[k]);
                ++k;
            }
        }
        tmp = sTab[5].split(" ");
        k = 0;
        for (int i = 0; i < B.length; i++){
            for (int j = 0; j < B[i].length; j++) {
                B[i][j] = Double.parseDouble(tmp[k]);
                ++k;
            }
        }
        tmp = sTab[6].split(" ");
        for (int i = 0; i < PI.length; i++)
            PI[i] = Double.parseDouble(tmp[i]);

        /* test
        System.out.println("Etats: ");
        affichage(Etats);
        System.out.println("Symboles: ");
        affichage(Symboles);
        System.out.println("Matrice de transition d'états: ");
        affichage2(A);
        System.out.println("Matrice d'observations: ");
        affichage2(B);
        System.out.println("Matrice des états initiaux: ");
        affichage1(PI);*/
    }

    private int[] getO_indices(String O){
    // Vecteur des indices de la suite d'observations. Au lieu d'avoir O = (o1 02...0n) on aura O_indice = (3 0 2 ...) où
    // chaque chiffre correspond à l'indice du symbole correspondant dans le vecteur des Symboles S
        String[] O_Tab = O.split(" ");
        int[] O_indice = new int[O_Tab.length];
        int k = 0, compteur;
        boolean trouve;
        for (String o : O_Tab) {
            trouve = false;
            compteur = 0;
            while (!trouve && compteur < Symboles.length) {
                if (o.equals(Symboles[compteur].trim())) {
                    O_indice[k] = compteur;
                    ++k;
                    trouve = true;
                } else {
                    ++compteur;
                }
            }
        }
        return O_indice;
    }

    protected double[][] getAlphas(String O) {
        // Fonction de calcul et remplissage de la matrice des variables Forward
        int[] O_indice = getO_indices(O);
        double[][] alpha = new double[Etats.length][O_indice.length];
        if (O_indice.length > 0) {
            for (int j = 0; j < alpha.length; j++) // Calcul de alpha[1][j]
                alpha[j][0] = PI[j] * B[j][O_indice[0]];
            for (int t = 1; t < O_indice.length; t++) { // Calcul de alpha[t+1:T][j]
                for (int j = 0; j < A.length; j++) {
                    double tmp = 0;
                    for (int i = 0; i < A.length; i++) {
                        tmp += alpha[i][t - 1] * A[i][j];
                    }
                    alpha[j][t] = B[j][O_indice[t]] * tmp;
                }
            }
        }
        return alpha;
    }

    protected double[][] getAlphas(double[][] A, double[][] B, double[] PI, String O) {
        // Fonction de calcul et remplissage de la matrice des variables Forward
        int[] O_indice = getO_indices(O);
        double[][] alpha = new double[Etats.length][O_indice.length];
        if (O_indice.length > 0) {
            for (int j = 0; j < alpha.length; j++) // Calcul de alpha[1][j]
                alpha[j][0] = PI[j] * B[j][O_indice[0]];
            for (int t = 1; t < O_indice.length; t++) { // Calcul de alpha[t+1:T][j]
                for (int j = 0; j < A.length; j++) {
                    double tmp = 0;
                    for (int i = 0; i < A.length; i++) {
                        tmp += alpha[i][t - 1] * A[i][j];
                    }
                    alpha[j][t] = B[j][O_indice[t]] * tmp;
                }
            }
        }
        return alpha;
    }

    protected double[][] getAlphas(String O, int t_bar) {
        // Fonction de calcul et remplissage de la matrice des variables Forward
        int[] O_indice = getO_indices(O);
        double[][] alpha = new double[Etats.length][t_bar];
        if (O_indice.length > 0) {
            for (int j = 0; j < alpha.length; j++) // Calcul de alpha[1][j]
                alpha[j][0] = PI[j] * B[j][O_indice[0]];
            for (int t = 1; t < t_bar; t++) { // Calcul de alpha[t+1:t_bar][j]
                for (int j = 0; j < A.length; j++) {
                    double tmp = 0;
                    for (int i = 0; i < A.length; i++) {
                        tmp += alpha[i][t - 1] * A[i][j];
                    }
                    alpha[j][t] = B[j][O_indice[t]] * tmp;
                }
            }
        }
        return alpha;
    }

    protected double[][] getBetas(double[][] A, double[][] B, String O) {
    // Fonction de calcul et remplissage de la matrice des variables Backward
        int[] O_indice = getO_indices(O);
        double[][] beta = new double[Etats.length][O_indice.length];
        if (O_indice.length > 0) {
            for (int j = 0; j < beta.length; j++) // Calcul de beta[T][i]
                beta[j][O_indice.length-1] = 1;
            for(int t = O_indice.length-1; t > 0; t--){ // Calcul de beta[T-1:1][i]
                for (int i = 0; i < A.length; i++){
                    double tmp = 0;
                    for (int j = 0; j < A.length; j++)
                        tmp += B[j][O_indice[t]]*beta[j][t]*A[i][j];
                    beta[i][t-1] = tmp;
                }
            }
        }
        return beta;
    }

    protected double[][] getBetas(String O) {
        // Fonction de calcul et remplissage de la matrice des variables Backward
        int[] O_indice = getO_indices(O);
        double[][] beta = new double[Etats.length][O_indice.length];
        if (O_indice.length > 0) {
            for (int j = 0; j < beta.length; j++) // Calcul de beta[T][i]
                beta[j][O_indice.length-1] = 1;
            for(int t = O_indice.length-1; t > 0; t--){ // Calcul de beta[T-1:1][i]
                for (int i = 0; i < A.length; i++){
                    double tmp = 0;
                    for (int j = 0; j < A.length; j++)
                        tmp += B[j][O_indice[t]]*beta[j][t]*A[i][j];
                    beta[i][t-1] = tmp;
                }
            }
        }
        return beta;
    }

    protected double[][] getBetas(String O, int t_bar) {
        // Fonction de calcul et remplissage de la matrice des variables Backward
        int[] O_indice = getO_indices(O);
        double[][] beta = new double[Etats.length][t_bar];
        if (O_indice.length > 0) {
            for (int j = 0; j < beta.length; j++) // Calcul de beta[t_bar][i]
                beta[j][t_bar-1] = 1;
            for(int t = t_bar-1; t > 0; t--){ // Calcul de beta[t_bar-1:1][i]
                for (int i = 0; i < A.length; i++){
                    double tmp = 0;
                    for (int j = 0; j < A.length; j++)
                        tmp += B[j][O_indice[t]]*beta[j][t]*A[i][j];
                    beta[i][t-1] = tmp;
                }
            }
        }
        return beta;
    }

    public double evaluer(String O){ // Evaluation selon l'algo de Forward-Backward
        // Choix aléatoire de t barre
        int t_bar = 1 + (int) (Math.random() * O.split(" ").length - 1);
        double ev = 0.0;
        double[][] alpha = getAlphas(O, t_bar);
        double[][] beta = getBetas(O, t_bar);
        //Calcul de Pr(O|Lambda)
        for (int i = 0; i < A.length; i++)
            ev += alpha[i][t_bar - 1]*beta[i][t_bar - 1];

        return ev;
    }
    // Cette méthode est juste la surcharge de la première afin d'avoir un gain du temps d'exécution
    public double evaluer(double[][] alpha, double[][] beta){
        double ev = 0.0;
        // Choix aléatoire de t barre
        int t_bar = (int) (Math.random() * alpha[0].length);
        //Calcul de Pr(O|Lambda)
        for (int i = 0; i < A.length; i++)
            ev += alpha[i][t_bar]*beta[i][t_bar];

        return ev;
    }

// ******************* Optimisation des paramètres du MMC **************************
    //Création de la matrice des Xi
    public double[][][] getXiTab(String O){
        double[][] alpha = getAlphas(O);
        double[][] beta = getBetas(O);
        int[] O_i = getO_indices(O);
        double ev = evaluer(alpha, beta);
        double[][][] Xi = new double[A.length][A.length][O_i.length-1];
        for (int i = 0; i < Xi.length; i++){
            for (int j = 0; j < Xi[0].length; j++){
                for (int t = 0; t < Xi[0][0].length; t++){
                   Xi[i][j][t] = (alpha[i][t] * A[i][j] * B[j][O_i[t+1]] * beta[j][t+1]) / ev;
                }
            }
        }
        return Xi;
    }

    // Optimisation de Pi
    public double[] getOptmizedPI(double[][] iGammas){
        double[] PI_barre = new double[PI.length];
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
                for (int j = 0; j < A.length; j++)
                    tmp += Xi[i][j][t];
                GTab[i][t] = tmp;
            }
        }
        return GTab;
    }

    protected double[][] getIGammas(double[][] alpha, double[][] beta){
        // iGamma = beta_t(i) * alpha_t(i) / Pr(O|lambda)
        double[][] GTab = new double[alpha.length][alpha[0].length];
        double ev = evaluer(alpha, beta);
        for (int t = 0; t < GTab[0].length; t++){
            for (int i = 0; i < GTab.length; i++){
                GTab[i][t] = alpha[i][t] * beta[i][t] / ev;
            }
        }
        return GTab;
    }

    protected double[][] getOptimizedA(double[][][] Xi, double[][] iGammas){
        // A_bar = somme t=1:T-1 de Xi_t(i,j) / somme t=1:T-1 de gamma_t(i)
        double[][] A_bar = new double[A.length][A.length];
        for (int i = 0; i < A.length; i++){
            for (int j = 0; j < A[0].length; j++) {
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

    protected double[][] getOptimizedB(double[][] iGammas, String O){
        double[][] B_bar = new double[B.length][B[0].length];
        final int T = O.split(" ").length;
        for (int i = 0; i < B_bar.length; i++){
            for (int j = 0; j < B_bar[0].length; j++){
                ArrayList<Integer> Ul = getLIndex(Symboles[j].trim(), O);
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

//*********************** Algorithmes de Baum-Welch ***********************************************

    public void BaumWelchMono(String O, int MaxIter, double epsilon){
        double[][] A_bar, B_bar;
        double[] Pi_bar;
        // Calculs des variables forward et backward
        double[][] alpha, beta, iGammas, alpha_bar, beta_bar;
        double[][][] Xi;
        alpha = getAlphas(O);
        beta = getBetas(O);
        boolean fin = false;
        double ev, ev_bar;
        int iter = 0;
        do{
            Xi = getXiTab(O);
            iGammas = getIGammas(alpha, beta);
            // Calcul des optimisés des éléments du modèle
            A_bar = getOptimizedA(Xi, iGammas);
            B_bar = getOptimizedB(iGammas, O);
            Pi_bar = getOptmizedPI(iGammas);
            // Recalcul des variables forward et backward mais avec les données optimisées
            alpha_bar = getAlphas(A_bar, B_bar, Pi_bar, O);
            beta_bar = getBetas(A_bar, B_bar, O);
            // Calcul de la valeur de l'évaluation du modèle initial et du modèle optimisé
            ev = evaluer(alpha, beta);
            ev_bar = evaluer(alpha_bar, beta_bar);
            iter++;
            if ((ev_bar - ev <= epsilon) || (iter > MaxIter))
                fin = true;
            else {
                // Lambda = Lambda_bar
                A = A_bar;
                B = B_bar;
                PI = Pi_bar;
            }
        }while (fin);
        // test
        System.out.println("Iterations: "+iter+"\nNouvelle matrice A");
        affichage2(A_bar);
        System.out.println("Nouvelle matrice B");
        affichage2(B_bar);
        System.out.println("Nouveau vecteur Pi");
        affichage1(Pi_bar);
    }

    protected void affichage3(double[][][] o){
        for (double[][] o1: o) {
            for (double[] o2 : o1){
                for (double o3 : o2)
                    System.out.print(o3+" ");
                System.out.println();
            }
        }
    }

    protected void affichage2(double[][] o){
        for (double[] oT: o) {
            for (Double oi : oT)
                System.out.print(oi+" ");
            System.out.println();
        }
    }
    protected void affichage(Object[] o){
        for (Object oi:o) {
            System.out.print(oi+" ");
        }
        System.out.println();
    }
    protected void affichage1(double[] o){
        for (double oi:o) {
            System.out.print(oi+" ");
        }
        System.out.println();
    }

}
