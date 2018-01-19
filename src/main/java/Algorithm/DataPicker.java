package Algorithm;

import Utils.FileUtil;
import Utils.TimeFormat;
import org.jblas.DoubleMatrix;

import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;
import java.util.Map;

public class DataPicker {
    public static String BASE_URL = "D:\\data\\speed-prediction\\";

    /*路段数，比结点数少1*/
    private static Map<String, Integer> segments = new HashMap<String, Integer>() {
        {
            put("G15", 740);
            put("S35", 117);
            put("G25", 304);
            put("G70", 212);
            put("G72", 224);
            put("G76", 233);
            put("G1501", 67);

            //另外两条高速数据异常
        }
    };

    private Map<String, Integer> segmentIndexAcc = new HashMap<String, Integer>();/*each road segment index in a row*/
    private int totalCount;
    private DateFormat df;
    private double[][] flows;
    private int flowsCurrIndex = 0;

    private static int HISTORY_SIZE = 4 * 24 * 7;
    private static int TIMESLOT = 5 * 60;
    private static int BEGIN_DATE = 20151101;
    private static int END_DATE = 20151104;
    private static int TEST_BEGIN_DATE = 20151103;
    private static String BEGIN_TIME_STRING = " 00:00:00";
    private static String END_TIME_STRING = " 24:00:00";


    DataPicker() {
        /*init segementIndexAcc*/
        int counter = 0;
        for (Map.Entry<String, Integer> entry : segments.entrySet()) {
            segmentIndexAcc.put(entry.getKey(), counter);
            counter += entry.getValue();
        }
        totalCount = segments.entrySet().stream().map(e -> e.getValue()).reduce(0, (acc, a) -> acc + a);
        flows = new double[HISTORY_SIZE][totalCount];
        df = new SimpleDateFormat("yyyy/mm/dd HH:MM:ss");
    }

    public void initFlows() {
        for (int i = BEGIN_DATE; i < END_DATE; i++) {
            long beginTime = TimeFormat.parse(i + BEGIN_TIME_STRING).getTime() / 1000;
            long endTime = TimeFormat.parse(i + END_TIME_STRING).getTime() / 1000;
            FileUtil fin = new FileUtil(BASE_URL + i + ".txt");
            String line = fin.readLine();

            for (long t = beginTime; t < endTime; t += TIMESLOT) {
                Arrays.fill(flows[flowsCurrIndex], -1);
                while (line != null && line.contains(df.format(new Date(t * 1000)))) {
                    insertItem(line, flowsCurrIndex);
                    flowsCurrIndex = (flowsCurrIndex + 1) % HISTORY_SIZE;
                }
            }
        }
        return;
    }

    private void insertItem(String line, int _flowsCurrIndex){
        String[] strs = line.split("\t");

        /*sample: 2016/11/01 00:15:00	G15	25	1	45*/
        int index = findIndex(strs[1],strs[2]);//
        double speed = Double.parseDouble(strs[3]);

        flows[_flowsCurrIndex][index]=speed;
        return;
    }

    private int findIndex(String roadId, String segmentId){
        return segmentIndexAcc.get(roadId).intValue()+Integer.parseInt(segmentId)-1;
    }

    private double[][] genToep(){
        double[][] Toep = new double[HISTORY_SIZE-1][HISTORY_SIZE];

        for(int i=0;i<HISTORY_SIZE-1;++i){
            Toep[i][i]=1;
            Toep[i][i+1]=-1;
        }

        return Toep;
    }

    public void output(double[][] flow, FileUtil fout) {
        int rowCount = flow.length;
        int columnCount = flow[0].length;

        String line = "";
        for (int i = 0; i < rowCount; i++) {
            line = "";
            for (int j = 0; j < columnCount; j++) {
                line += String.valueOf(flow[i][j]);
                line += "\t";
            }
            fout.writeLine(line);
        }
        fout.close();
    }

    public static void main(String[] args){
        DataPicker dp = new DataPicker();
        dp.initFlows();

        double[][] M = dp.flows;
        double[][] T = dp.genToep();

        dp.output(M, new FileUtil(BASE_URL+"M.txt"));

        Smatrix smatrix = new Smatrix(M);
        double[][] S = smatrix.calcS();

        dp.output(S, new FileUtil(BASE_URL+"S.txt"));

        CSstALS cSstALS = new CSstALS(new DoubleMatrix(M),new DoubleMatrix(S),new DoubleMatrix(T));
        double[][] data = cSstALS.estimate();

        dp.output(data, new FileUtil(BASE_URL+"estimate.txt"));
    }
}