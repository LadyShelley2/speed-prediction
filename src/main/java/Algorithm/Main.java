package Algorithm;

import Utils.FileUtil;
import org.jblas.DoubleMatrix;

public class Main {
    public static void main(String[] args){
        DataPicker dp = new DataPicker();
        dp.initFlows();

        double[][] M = dp.flows;
        double[][] T = dp.genToep();

        dp.output(M, new FileUtil(dp.BASE_URL+"M.txt"));

        Smatrix smatrix = new Smatrix(M);
        double[][] S = smatrix.calcS();

        dp.output(S, new FileUtil(dp.BASE_URL+"S.txt"));

        CSstALS cSstALS = new CSstALS(new DoubleMatrix(M),new DoubleMatrix(S),new DoubleMatrix(T));
        double[][] data = cSstALS.estimate();

        dp.output(data, new FileUtil(dp.BASE_URL+"estimate.txt"));
        dp.output(dp.test_base, new FileUtil(dp.BASE_URL+"test_base.txt"));
        dp.output(data[dp.HISTORY_SIZE-1], new FileUtil(dp.BASE_URL+"last_row.txt"));

        double mape = cSstALS.getMAPEValid(new DoubleMatrix(dp.test_base),new DoubleMatrix(data[dp.HISTORY_SIZE-1]));
        double rmse = cSstALS.getRMSE(new DoubleMatrix(dp.test_base),new DoubleMatrix(data[dp.HISTORY_SIZE-1]));

        System.out.println("mape is :"+mape);
        System.out.println("rmse is :"+rmse);
    }
}
