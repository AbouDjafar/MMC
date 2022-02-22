public class IHM {
    public static void main(String[] args){
        MMC model = new MMC(System.getProperty("user.dir")+"/exo1_1.txt");
        //String O = "Visionner Dormir Danser";
        String[] O = {"Visionner Dormir Danser", "Manger Manger Dormir Manger", "Danser Jouer Dormir Danser", "Manger Danser Visionner Dormir Manger"};
        /*double[][][] Xi = model.getXiTab(O);
        double[][] ABar = model.getOptimizedA(Xi, iGammas);
        double[][] alpha, beta, iGammas, iGammas2;
        alpha = model.getAlphas(O);
        beta = model.getBetas(O);
        iGammas = model.getIGammas(alpha, beta);
        iGammas2 = model.getIGammas2(Xi);*/
        System.out.println("*** Affichage  ***");
        //System.out.println(model.BaumWelchMonoSeq(O, model, 100, 0, false));
        System.out.println(model.BaumWelchMultiSeq(O, model, 100, 0, true));
    }
}
