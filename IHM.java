public class IHM {
    public static void main(String[] args){
        MMC model = new MMC("C:\\Users\\Djafar\\Documents\\exo1_1.txt");
        String O = "Visionner Dormir Danser Manger Dormir";
        /*double[][][] Xi = model.getXiTab(O);
        double[][] ABar = model.getOptimizedA(Xi, iGammas);
        double[][] alpha, beta, iGammas, iGammas2;
        alpha = model.getAlphas(O);
        beta = model.getBetas(O);
        iGammas = model.getIGammas(alpha, beta);
        iGammas2 = model.getIGammas2(Xi);*/
        System.out.println("*** Affichage  ***");
        model.BaumWelchMono(O, 100, 0);
    }
}
