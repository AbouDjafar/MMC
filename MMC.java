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
    public MMC(String fileURI){
        modelIN(fileURI, this);
        espaceNom = fileURI.replace(".txt", "");
    }

    public MMC(MMC model){
        this.espaceNom = model.getEspaceNom();
        this.setA(model.getA());
        this.setB(model.getB());
        this.setPI(model.getPI());
        this.setEtats(model.getEtats());
        this.setSymboles(model.getSymboles());
    }

    protected int[] get_Tab_indices(String S, int param){
    // Vecteur des indices de la suite d'observations. Au lieu d'avoir O = (o1 02...0n) ou Q = (q1 q2 ... qn) on aura O_indice = (3 0 2 ...) où
    // chaque chiffre correspond à l'indice du symbole ou de l'état correspondant dans le vecteur des Symboles S
    // param == 0 : Séquences et param == 0: Etats
        param = (param == 1) ? param : 0;
        String[] O_Tab = S.split(" ");
        int[] S_indice = new int[O_Tab.length];
        int k = 0, compteur;
        boolean trouve;
        String[][] SE = {Symboles, Etats};
        for (String o : O_Tab) {
            trouve = false;
            compteur = 0;
            while (!trouve && compteur < SE[param].length) {
                if (o.equals(SE[param][compteur].trim())) {
                    S_indice[k] = compteur;
                    ++k;
                    trouve = true;
                } else {
                    ++compteur;
                }
            }
        }
        return S_indice;
    }
//******************* Evaluation Naive ************************************************
    public double evaluationNaive(String O, String Q, MMC model){
        double produit_O = 1, produit_Q;
        double[] _PI = model.getPI();
        double[][] _A = model.getA(), _B = model.getB();
        int[] O_tab = get_Tab_indices(O, 0);
        int[] Q_tab = get_Tab_indices(Q, 1);
        produit_Q = _PI[Q_tab[0]];
        for (int i = 1; i < Q_tab.length; i++)
            produit_Q *= _A[Q_tab[i - 1]][Q_tab[i]];
        for (int i = 0; i < O_tab.length; i++)
            produit_O *= _B[Q_tab[i]][O_tab[i]];

        return (produit_Q * produit_O);
    }

//***************** Lecture et écriture de modèle **************************************
    //--------------- Lecture et parsing des données du modèle à partir d'un fichier ------------------------
    protected void modelIN(String fileURL, MMC M) {
        File model = new File(fileURL);
        BufferedReader br;
        StringBuilder s = new StringBuilder();
        double[][] _A, _B;
        double[] _PI;
        String[] _Etats, _Symboles;
        int k;
        //------------------- Lecture du fichier source -----------------------------------
        try {
            br = new BufferedReader(new FileReader(model));
            String str;
            while((str = br.readLine()) != null){
                s.append(str).append('\n');
            }
            br.close();
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
    protected void modelOUT(String outputFileURI, MMC model){
        BufferedWriter bw;
        String str;
        try {
            bw = new BufferedWriter(new FileWriter(outputFileURI));
            str = model.getEtats().length+"\n";
            bw.write(str);
            str = model.getSymboles().length+"\n";
            bw.write(str);
            litteralsWriter(model.getEtats(), bw);
            litteralsWriter(model.getSymboles(), bw);
            doubleDimMatrixWriter(model.getA(), bw);
            doubleDimMatrixWriter(model.getB(), bw);
            singleDimMatrixWriter(model.getPI(), bw);

            bw.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private void litteralsWriter(String[] sTab, BufferedWriter bw) throws IOException {
        StringBuilder tmp = new StringBuilder();
        for (String s : sTab) tmp.append(s).append(" ");
        tmp.deleteCharAt(tmp.length() - 1);
        tmp.append("\n");
        bw.write(""+tmp.toString());
    }
    private void doubleDimMatrixWriter(double[][] dTab, BufferedWriter bw) throws IOException {
        StringBuilder tmp = new StringBuilder();
        for (double[] ligne : dTab){
            for (double cellule : ligne)
                tmp.append(cellule).append(" ");
        }
        tmp.deleteCharAt(tmp.length() - 1);
        tmp.append("\n");
        bw.write(tmp.toString());
    }
    private void singleDimMatrixWriter(double[] sTab, BufferedWriter bw) throws IOException {
        StringBuilder tmp = new StringBuilder();
        for (double s : sTab) tmp.append(s).append(" ");
        tmp.deleteCharAt(tmp.length() - 1);
        tmp.append("\n");
        bw.write(""+tmp.toString());
    }

    @Override
    public String toString(){
        StringBuilder str = new StringBuilder();
        str.append("Etats: \n");
        str.append(affichage(Etats));
        str.append("Symboles: \n");
        str.append(affichage(Symboles));
        str.append("Matrice de transition d'états: \n");
        str.append(affichage2(A));
        str.append("Matrice d'observations: \n");
        str.append(affichage2(B));
        str.append("Matrice des états initiaux: \n");
        str.append(affichage1(PI));

        return str.toString();
    }

//***************************** Affichages (test) ************************************

    private String affichage3(double[][][] o){
        StringBuilder str = new StringBuilder();
        for (double[][] o1: o) {
            for (double[] o2 : o1){
                for (double o3 : o2)
                    str.append(o3).append(" ");
                str.append("\n");
            }
        }

        return str.toString();
    }

    private String affichage2(double[][] o){
        StringBuilder str = new StringBuilder();
        for (double[] oT: o) {
            for (Double oi : oT)
                str.append(oi).append(" ");
            str.append("\n");
        }

        return str.toString();
    }
    private String affichage(String[] o){
        StringBuilder str = new StringBuilder();
        for (String oi:o) {
            str.append(oi.trim()).append(" ");
        }
        str.append("\n");

        return str.toString();
    }
    private String affichage1(double[] o){
        StringBuilder str = new StringBuilder();
        for (double oi:o) {
            str.append(oi).append(" ");
        }
        str.append("\n");

        return str.toString();
    }

}
