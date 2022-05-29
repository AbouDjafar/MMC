public class Viterbi {

    private double[] AMax(double[] tab){
        double max = 0;
        int argmax = -1;
        for (int i = 0; i < tab.length; i++){
            if (max < tab[i]){
                max = tab[i];
                argmax = i;
            }
        }

        return (new double[]{argmax, max});
    }

    public String viterbi(String O, MMC model){
        final int T = O.split(" ").length, N = model.getEtats().length;
        int[] O_Tab = model.get_Tab_indices(O, 0), chemin = new int[T];
        double[] _Pi = model.getPI(), tmp = new double[N], maxis; //tmp: juste un tableau tampon et maxis servira de support aux résultats de la fonction AMax pour avoir max er argmax
        double[][] _A = model.getA(), _B = model.getB(), delta; // delta matrice des proba de passer par un noeud du treillis
        int[][] psi; //matrice des indices des états (noeuds du treilli)
        delta = new double[N][T];
        psi = new int[N][T];
        for(int i = 0; i < N; i++) { // A l'initialisation
            delta[i][0] = _Pi[i] * _B[i][O_Tab[0]];
            psi[i][0] = i;
        }
        for (int t = 1; t < T; t++){ // on procède temps par temps
            for (int i = 0; i < N; i++){ // ensuite état par état
                for (int j = 0; j < N; j++) // chaque état est susceptible d'être la débouchée d'uns transition au temps t-1
                    tmp[j] = delta[j][t-1] * _A[i][j];
                maxis = AMax(tmp);
                delta[i][t] = maxis[1] * _B[i][O_Tab[t]];
                psi[i][t] = (int)maxis[0];
            }
        }
        //Construction du chemin optimal
        for (int i = 0; i < N - 1; i++)
            tmp[i] = delta[i][T-1];
        chemin[T-1] = (int)AMax(tmp)[0];
        for (int t = T - 2; t >= 0; t--)
            chemin[t] = psi[chemin[t+1]][t+1];
        StringBuilder Q = new StringBuilder();
        for (int i : chemin)
            Q.append(model.getEtats()[i].trim()).append(" ");

        return Q.toString();
    }

}
