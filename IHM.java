public class IHM {
    public static void main(String[] args){
        MMC model = new MMC(System.getProperty("user.dir")+"/groupeDol.txt");
        String O = "Rire Dormir Rire Dormir Mentir Rire Transpirer Rire Demontrer Rire";
        String Q = "Eric Dol Eric Dol Dol Eric Dol Eric Dol Eric";
        //String[] O = {"Visionner Dormir Danser", "Manger Manger Dormir Manger", "Danser Jouer Dormir Danser", "Manger Danser Visionner Dormir Manger"};
        /*double[][][] Xi = model.getXiTab(O);
        double[][] ABar = model.getOptimizedA(Xi, iGammas);
        double[][] alpha, beta, iGammas, iGammas2;
        alpha = model.getAlphas(O);
        beta = model.getBetas(O);
        iGammas = model.getIGammas(alpha, beta);
        iGammas2 = model.getIGammas2(Xi);*/
        System.out.println("*** Affichage  ***");
        //System.out.println(model.BaumWelchMonoSeq(O, model, 100, 0, false));
        //System.out.println(model.BaumWelchMultiSeq(O, model, 100, 0, true));
        System.out.println(model.evaluationNaive(O, Q, model));
        System.out.println(model.evaluer(O, model));
    }
}
