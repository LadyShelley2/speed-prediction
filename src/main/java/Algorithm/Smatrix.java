package Algorithm;

import org.jblas.*;

import java.util.Comparator;
import java.util.PriorityQueue;

public class Smatrix {

    private int historySize, topK =10;
    private double[][] flows;

    Smatrix(double[][] _flows){
        flows=_flows;
        historySize =_flows.length;
    }

    /**
     * 获得路段的线性回归矩阵，每一行代表一个路段与其它路段间的相关系数
     * @return
     */
    double[][] calcS(){
        int n=flows[0].length;
        double S[][] = new double[n][n];
        double[][]pear = getPearson();
        for(int i=0;i<n;++i){
            int[] top = selectTopk(pear, i);
            double x[][] = new double[historySize][topK];
            double y[] = new double[historySize];
            for(int j=0;j<top.length;++j)
                for(int k=0;k<historySize;++k)
                    x[k][j] = flows[k][top[j]];
            for(int k=0;k<historySize;++k)
                y[k] = flows[k][i];

            DoubleMatrix xm = new DoubleMatrix(x);
            DoubleMatrix ym= new DoubleMatrix(y);
            double[] coef=Solve.pinv(xm.transpose().mmul(xm)).mmul(xm.transpose()).mmul(ym).toArray();
            for(int j=0;j<top.length;++j)
                S[i][top[j]] = coef[j];
            S[i][i] = -1;
        }
        return S;
    }
    //根据相似函数矩阵close，选出与i最近的topK个
    int[] selectTopk(double[][] close,int i){
        int[] top = new int[topK];
        int n= flows[0].length;
        PriorityQueue<Integer> pq = new PriorityQueue<Integer>(new Comparator<Integer>(){ //建立小顶堆，存储i的K近邻
            @Override
            public int compare(Integer arg0, Integer arg1) {
                if(close[i][arg0]<close[i][arg1]) return -1;
                else if(close[i][arg0]>close[i][arg1]) return 1;
                return 0;
            }
        });
        for(int j=0;j<n;++j)
            if(j!=i&&pq.size()<topK)
                pq.add(j);
            else if(j!=i&&close[i][j]>pq.peek()){
                pq.remove();
                pq.add(j);
            }
        int k=0;
        for(Integer t:pq)
            top[k++] = t;
        return top;
    }
    //计算任意两列的     相关系数/距离
    double[][] getPearson(){
        int n = flows[0].length,m=flows.length;
        double pear[][] = new double[n][n];
        for(int i=0;i<n;++i){
            for(int j=0;j<i;++j){
                double x=0,y=0,xy=0,x2=0,y2=0;
                for(int k=0;k<m;++k){
                    x+=flows[k][i];
                    y+=flows[k][j];
                    xy+=flows[k][i]*flows[k][j];
                    x2+=flows[k][i]*flows[k][i];
                    y2+=flows[k][j]*flows[k][j];
                }
                pear[i][j] = pear[j][i]=(m*xy-x*y)/Math.sqrt((m*x2-x*x)*(m*y2-y*y));///dist[i][j];
            }
            pear[i][i]=1;
        }
        return pear;
    }
}
