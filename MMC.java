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
    private String espaceNom;
    private String[] Etats;
    private String[] Symboles;
    private double[][] A;
    private double[][] B;
    private double[] PI;

//************ Getters & Setters ***********************************
    public String getEspaceNom(){
        return espaceNom;
    }

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

    public String[] getEtats() {
        return Etats;
    }

    public void setEtats(String[] etats) {
        Etats = etats;
    }

    public String[] getSymboles() {
        return Symboles;
    }

    public void setSymboles(String[] symboles) {
        Symboles = symboles;
    }

//************ Constructeurs ********************************************************
    public MMC(String fileURL){
        modelIN(fileURL, this);
        espaceNom = fileURL.replace(".txt", "");
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

    public MMC(MMC model){
        this.espaceNom = model.getEspaceNom();
        this.setA(model.getA());
        this.setB(model.getB());
        this.setPI(model.getPI());
        this.setEtats(model.getEtats());
        this.setSymboles(model.getSymboles());
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

//******************* Variables forward *************************************************
    protected double[][] getAlphas(String O) {
        return getAlphas(O, this);
    }
    protected double[][] getAlphas(String O, MMC model) {
        // Fonction de calcul et remplissage de la matrice des variables Forward
        int[] O_indice = getO_indices(O);
        double[][] alpha = new double[Etats.length][O_indice.length];
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
    protected double[][] getBetas(String O) {
        return getBetas(O, this);
    }
    protected double[][] getBetas(String O, MMC model) {
        // Fonction de calcul et remplissage de la matrice des variables Backward
        int[] O_indice = getO_indices(O);
        double[][] beta = new double[Etats.length][O_indice.length];
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
    public double evaluer(double[][] alpha, double[][] beta){
        double ev = 0.0;
        // Choix aléatoire de t barre
        int t_bar = (int) (Math.random() * alpha[0].length);
        //Calcul de Pr(O|Lambda)
        for (int i = 0; i < A.length; i++)
            ev += alpha[i][t_bar]*beta[i][t_bar];

        return ev;
    }
    public double evaluer(String O, MMC model){
        double ev = 0.0;
        double[][] alpha = getAlphas(O, model);
        double[][] beta = getBetas(O, model);
        return evaluer(alpha, beta);
    }
    public double evaluer(String O){ // Evaluation selon l'algo de Forward-Backward
        return evaluer(O, this);
    }

// ******************* Optimisation des paramètres du MMC **************************

    //---------------------- Création de la matrice des Xi -----------------------------------------
    public double[][][] getXiTab(String O, double[][] alpha, double[][] beta, MMC model){
        int[] O_i = getO_indices(O);
        double ev = evaluer(alpha, beta);
        double[][][] Xi = new double[A.length][A.length][O_i.length-1];
        for (int i = 0; i < Xi.length; i++){
            for (int j = 0; j < Xi[0].length; j++){
                for (int t = 0; t < Xi[0][0].length; t++){
                    Xi[i][j][t] = (alpha[i][t] * model.getA()[i][j] * model.getB()[j][O_i[t+1]] * beta[j][t+1]) / ev;
                }
            }
        }
        return Xi;
    }
    public double[][][] getXiTab(String O){
        double[][] alpha = getAlphas(O);
        double[][] beta = getBetas(O);
        return getXiTab(O, alpha, beta, this);
    }

    //----------------- Optimisation de Pi --------------------------------------------
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

    // ---------------- Optimisation de la matrice des transitions d'états ------------------------
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

    //-------------------- Optimisation de la matrice des observations --------------------------------------
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
//------------------- Baum-Welch mono séquence -----------------------------------------
    public void BaumWelchMono(String O, MMC model, int MaxIter, double epsilon){
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
            alpha = getAlphas(O, model);
            beta = getBetas(O, model);
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
            alpha_bar = getAlphas(O, model_bar);
            beta_bar = getBetas(O, model_bar);
            // Calcul de la valeur de l'évaluation du modèle initial et du modèle optimisé
            ev = evaluer(alpha, beta);
            ev_bar = evaluer(alpha_bar, beta_bar);
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
        System.out.println("Iterations: "+iter+"\nNouvelle matrice A");
        affichage2(model_bar.getA());
        System.out.println("Nouvelle matrice B");
        affichage2(model_bar.getB());
        System.out.println("Nouveau vecteur Pi");
        affichage1(model_bar.getPI());
        //------------ écriture sur fichier du nouveau model -----------------------------
        modelOUT(model.getEspaceNom()+"_optimized.txt", model_bar);
    }

//--------- Baum-Welch multi séquences --------------------------

//***************** Lecture et écriture de modèle **************************************
    //--------------- Lecture et parsing des données du modèle à partir d'un fichier ------------------------
    protected void modelIN(String fileURL, MMC M) {
        File model = new File(fileURL);
        FileReader fr;
        StringBuilder s = new StringBuilder();
        double[][] _A, _B;
        double[] _PI;
        String[] _Etats, _Symboles;
        int k;
        //------------------- Lecture du fichier source -----------------------------------
        try {
            fr = new FileReader(model);
            while((k = fr.read()) != -1){
                s.append((char) k);
            }
            fr.close();
        } catch (IOException e){
            e.printStackTrace();
        }
        //----------------- Allocation des dimensions aux différents tableaux -------------------------
        String[] sTab = s.toString().split("\n");
        _Etats = new String[Integer.parseInt(sTab[0].trim())];
        _Symboles = new String[Integer.parseInt(sTab[1].trim())];
        _A = new double[_Etats.length][_Etats.length];
        _B = new double[_Etats.length][_Symboles.length];
        _PI = new double[_Etats.length];
        //----------------- Remplissage des tableaux alloués -------------------------------------------
        _Etats = sTab[2].split(" ");
        _Symboles = sTab[3].split(" ");
        String[] tmp = sTab[4].split(" ");
        k = 0;
        for (int i = 0; i < _A.length; i++){
            for (int j = 0; j < _A[i].length; j++) {
                _A[i][j] = Double.parseDouble(tmp[k]);
                ++k;
            }
        }
        tmp = sTab[5].split(" ");
        k = 0;
        for (int i = 0; i < _B.length; i++){
            for (int j = 0; j < _B[i].length; j++) {
                _B[i][j] = Double.parseDouble(tmp[k]);
                ++k;
            }
        }
        tmp = sTab[6].split(" ");
        for (int i = 0; i < _PI.length; i++)
            _PI[i] = Double.parseDouble(tmp[i]);
        //----------------------- attribution des paramètres du modèle -----------------------------------
        M.setEtats(_Etats);
        M.setSymboles(_Symboles);
        M.setA(_A);
        M.setB(_B);
        M.setPI(_PI);
    }
    //-------------- écriture du modèle -----------------------------------------------
    protected void modelOUT(String outputFileURL, MMC model){
        FileWriter fw;
        try {
            fw = new FileWriter(outputFileURL);
            fw.write(""+model.getEtats().length+"\n"+model.getSymboles().length+"\n");
            litteralsWriter(model.getEtats(), fw);
            litteralsWriter(model.getSymboles(), fw);
            doubleMatrixWriter(model.getA(), fw);
            doubleMatrixWriter(model.getB(), fw);
            for (double d : model.getPI())
                fw.write(""+d+" ");
            fw.write("\b");

            fw.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private void litteralsWriter(String[] sTab, FileWriter fw) throws IOException {
        StringBuilder tmp = new StringBuilder();
        for (String s : sTab) tmp.append(s).append(" ");
        int idx = tmp.lastIndexOf(" ");
        tmp.deleteCharAt(idx);
        fw.write(""+tmp.toString());
    }
    private void doubleMatrixWriter(double[][] dTab, FileWriter fw) throws IOException {
        for (double[] ligne : dTab){
            for (double cellule : ligne)
                fw.write(""+cellule+" ");
            fw.write("\b");
        }
        fw.write("\n");
    }

//***************************** Affichages (test) ************************************

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
