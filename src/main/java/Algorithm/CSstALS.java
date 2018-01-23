package Algorithm;

import org.jblas.DoubleMatrix;
import org.jblas.MatrixFunctions;
import org.jblas.Solve;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.function.Predicate;
import java.util.stream.Collectors;

/**
 * 交替最小二乘求解时空约束的压缩感知
 * @author cyl
 *
 */
public class CSstALS {
    static int r=2,iteration=10;
    static double lamda=100.0,yita = 1,eps = 1e-10;
    private DoubleMatrix Ma = null,S=null,T=null;
    public CSstALS(DoubleMatrix M,DoubleMatrix S,DoubleMatrix T){
        Ma = M;
        this.S = S;
        this.T = T;
    }
    public double[][] estimate(){
        int m = Ma.rows,n = Ma.columns;
        double[][] B = new double[m][n];
        for(int i=0;i<m;++i)
            for(int j=0;j<n;++j)
                if(Ma.get(i, j)>=eps){
                    B[i][j]=1;
                }
        DoubleMatrix Mb = new DoubleMatrix(B);
        double v0 = Double.MAX_VALUE;
        DoubleMatrix Lo = null,Ro = null;

        DoubleMatrix L = DoubleMatrix.rand(m, r).mul(10);
        for(int iter = 0;iter<iteration;++iter){
            DoubleMatrix R = myLeastSquare(L,Ma,Mb,S,T);
            L = myLeastSquare(R,Ma.transpose(),Mb.transpose(),T,S);
            DoubleMatrix LR = L.mmul(R.transpose());
            double v = sumSquare(Mb.mul(LR).sub(Ma).toArray2())+lamda*(sumSquare(L.toArray2())+sumSquare(R.toArray2()))+
                    yita*sumSquare(T.mmul(LR).toArray2())+yita*sumSquare(LR.mmul(S.transpose()).toArray2());
            if(v<v0){
                Lo = L;Ro = R;
                v0 = v;
            }
        }
        return Lo.mmul(Ro.transpose()).toArray2();
    }
    //
    private static DoubleMatrix myLeastSquare(DoubleMatrix L,DoubleMatrix M,DoubleMatrix B,DoubleMatrix S,DoubleMatrix T){
        int m =  M.rows,n = M.columns;
        DoubleMatrix LT = L.transpose().mmul(T.transpose()).mmul(T).mmul(L).mul(yita).add(DoubleMatrix.eye(r).muli(lamda));
        DoubleMatrix Sts = S.transpose().mmul(S),Ltl = L.transpose().mmul(L).mul(yita);
        DoubleMatrix Mtl = M.transpose().mmul(L);
        double[][] equat = new double[n*r][n*r]; //n*r个方程，n*r个未知数
        double y[] = new double[n*r];
        for(int i=0;i<n;++i){
            DoubleMatrix LBL = L.transpose().mmul(DoubleMatrix.diag(B.getColumn(i))).mmul(L);
            for(int j=0;j<r;++j){
                int number = i*r+j;
                y[number] = Mtl.get(i, j);
                for(int k=0;k<r;++k){
                    equat[number][i*r+k] = LBL.get(k, j)+LT.get(k, j);
                }
                for(int k1=0;k1<n;++k1)
                    for(int k2 = 0;k2<r;++k2)
                        equat[number][k1*r+k2]+=Sts.get(i, k1)*Ltl.get(k2, j);
            }
        }
        double colR[] = Solve.pinv(new DoubleMatrix(equat)).mmul(new DoubleMatrix(y)).toArray();
        double[][] R = new double[n][r];
        for(int i=0;i<n;++i)
            for(int j=0;j<r;++j)
                R[i][j] = colR[i*r+j];
        return new DoubleMatrix(R);
    }

    private static double sumSquare(double[][] a){
        double s = 0;
        for(int i=0;i<a.length;++i)
            for(int j=0;j<a[0].length;++j)
                s+=a[i][j]*a[i][j];
        return s;
    }
    private static void printArray(double a[][]){
        for(int i=0;i<a.length;++i){
            for(int j=0;j<a[0].length;++j)
                System.out.print(a[i][j]+" ");
            System.out.println();
        }
    }
    public static double getMAPE(DoubleMatrix base, DoubleMatrix estimate){
        int m = base.rows, n = base.columns;
        return MatrixFunctions.abs(base.sub(estimate)).divi(base.add(DoubleMatrix.ones(m,n).mmul(eps))).sum()/base.length;
    }

    /**
     * 去除实际值中没有值得元素
     * @param base
     * @param estimate
     * @return
     */
    public static double getMAPEValid(DoubleMatrix base, DoubleMatrix estimate){

        int m = base.rows, n = base.columns;
        double[][] baseArr = base.toArray2();
        double[][] estimateArr = estimate.toArray2();

        double acc = 0.0;
        int counter =0;
        for(int i=0;i<m;i++){
            for(int j =0;j<n;j++){
                if(baseArr[i][j]<0.0) continue;
                if(baseArr[i][j]==0.0){
                    counter++;
                    continue;
                }
                acc+=Math.abs((estimateArr[i][j]-baseArr[i][j])/baseArr[i][j]);
                counter++;
            }
        }

        return acc/counter;
    }

    public double getRMSE(DoubleMatrix base, DoubleMatrix estimate){
//        System.out.println(base.sub(estimate));
        System.out.println(MatrixFunctions.pow(base.sub(estimate),2).sum());
        return MatrixFunctions.sqrt(MatrixFunctions.pow(base.sub(estimate),2).sum()/base.length);
    }
    public static void main(String[] args) {

        DoubleMatrix base = new DoubleMatrix(new double[][]{{1,2,3},{2,4,4}});

        DoubleMatrix estimate = new DoubleMatrix(new double[][]{{1,3,3},{2,3,4}});

        System.out.println(estimate.toArray2()[0]);
        System.out.println(CSstALS.getMAPE(base,estimate));
        System.out.println(CSstALS.getMAPEValid(base,estimate));

        double[][] a = {{1,2,3},{3,4,5}};
        double[] b = Arrays.copyOf(a[0],a[0].length);
        a[0][0]=-1;

        System.out.println("______________");

    }
}
