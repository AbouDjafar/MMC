import java.util.Scanner;
import javax.swing.*;
import javax.swing.filechooser.FileView;

public class IHM {
    private static Scanner sc;

    public static void main(String[] args){
        MMC model = new MMC(System.getProperty("user.dir")+"/exo1_1.txt");
        /*BaumWelch BW = new BaumWelch(model);
        //String O = "Rire Dormir Rire Dormir Mentir Rire Transpirer Rire Demontrer Rire";
        //String O = "Visionner Dormir Danser Manger Jouer Dormir";
        //String Q = "Eric Dol Eric Dol Dol Eric Dol Eric Dol Eric";
        String[] O = {"Visionner Dormir Danser", "Manger Manger Dormir Manger", "Danser Jouer Dormir Danser", "Manger Danser Visionner Dormir Manger"};
        System.out.println("*** Affichage  ***");
        //double[][] alpha = new ForwardBackward().getAlphas(O, model), beta = new ForwardBackward().getBetas(O, model);
        //double[][][] Xi = BW.getXiTab(O, alpha, beta, model);
        //double[][] iGammas = BW.getIGammas(alpha, beta, model), iGammas2 = BW.getIGammas2(Xi);
        //System.out.println(BW.MonoSequence(O, model, 100, 0, false));
        System.out.println(BW.MultiSequence2(O, model, 100, 0, false));
        //System.out.println("Naive subjectif: "+model.evaluationNaive(O, Q, model));
        //Q = new Viterbi().viterbi(O, model);
        //System.out.println(Q);
        //System.out.println("Naive optimal: "+model.evaluationNaive(O, Q, model));
        //System.out.println("forward-backward: "+model.evaluer(O, model));*/

        int usr_ch;
        do {
            //TODO: Charger au préalale un modèle
            System.out.print("\n1. charger un nouveau modèle \n" +
                    "2. evaluer une séquence (Forward-Backward)\n" +
                    "3. optimiser le modele selon une sequence de symboles (Baum-Welch mono)\n" +
                    "4. optimiser le modele selon plusieurs sequences de symboles (Baum-Welch multi)\n" +
                    "5. ressortir la sequence d'etats (Viterbi) \n" +
                    "6. Afficher le modele \n" +
                    "0. Quitter \n" +
                    "-> Que voulez-vous faire [0-5]: ");
            sc = new Scanner(System.in);
            usr_ch = sc.nextInt();
            switch (usr_ch) {
                case 0:
                    System.out.println("Fin");
                break;
                case 1:
                    model = ajoutModel();
                break;
                case 2:
                    double ev = evaluation(model);
                    System.out.println("Pr[0|\\lambda] = "+ev);
                break;
                case 3:
                    optimisation_mono(model);
                break;
                case 4:
                    optimisation_multi(model);
                break;
                case 5:
                    System.out.println(viterbi(model));
                break;
                case 6:
                    System.out.println(model.toString());
                break;
                default:
                    System.err.println("un chiffre compris entre 0 et 6 svp");
                break;
            }
        }while (usr_ch != 0);
    }

    private static double evaluation(MMC model){
        String O = lectureSequence();
        return (new ForwardBackward().evaluer(O, model));
    }
    private static void optimisation_mono(MMC model){
        String O = lectureSequence();
        BaumWelch bw = new BaumWelch(model);
        System.out.println(bw.MonoSequence(O, model, 100, 0, false).toString());
    }
    private static void optimisation_multi(MMC model){
        System.out.println("Nombre de séquences à saisir: ");
        final int n = sc.nextInt();
        String[] O = new String[n];
        for (int i = 0; i < n; i++) {
            System.out.print((i + 1) + " - ");
            O[i] = lectureSequence();
        }
        BaumWelch bw = new BaumWelch(model);
        System.out.println(bw.MultiSequence2(O, model, 100, 0, false).toString());
    }
    private static String viterbi(MMC model){
        String O = lectureSequence();
        return (new Viterbi().viterbi(O, model));
    }
    private static MMC ajoutModel(){
        System.out.println("/!\\ Formatage des fichier d'entrée des MMC /!\\\n" +
                "    ligne 1: nombre total des etats (N)\n" +
                "    ligne 2: nombre total des symboles (T)\n" +
                "    ligne 3: lister les etats, separes par un espace (e1 e2 e3 e4 ... eN)\n" +
                "    ligne 4: lister les symboles, separes par un espace (s1 s2 s3 s4 ... sT)\n" +
                "    ligne 5: derouler la matrice de transition d'etats (A11 A12 .. A1N A21 A22 ... A2N A31 A32 ... A3N ... AN1 AN2 ... ANN)\n" +
                "    ligne 6: derouler la matrice d'observations (B11 B12 ...B1T B21 B22 ... B2T ... BN1 BN2 ... BNT)\n" +
                "    ligne 7: lister le vecteur des états initiaux (pi1 pi2 pi3 ... piN)" +
                "    \n(Le selecteur de fichier apparait souvent en arriere-plan)");

        JFileChooser selecteur = new JFileChooser();
        selecteur.setFileSelectionMode(JFileChooser.FILES_ONLY);
        selecteur.showOpenDialog(new JPanel());
        String fichier = selecteur.getSelectedFile().getAbsolutePath();
        return (new MMC(fichier));
    }

    private static String lectureSequence() {
        System.out.println("Entrez une séquence de longueur n (O_1 O_2 ... O_n):");
        sc = new Scanner(System.in);
        return sc.nextLine();
    }

}
