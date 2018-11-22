package consoperators;

import java.text.DecimalFormat;

/***
 * N阶行列式求解
 * @author imlilu
 *
 */
//public class Test {


public class JacobianMatrixDeterminant {

    /*test
    public static void main(String[] args) {
//		double[][] test = {{2,1,-1},{4,-1,1},{201,102,-99}}; 			结果为-18
//		double[][] test = {{1,1,-1,3},{-1,-1,2,1},{2,5,2,4},{1,2,3,2}};  结果为33
//		double[][] test = {{1,0,-1,2},{-2,1,3,1},{0,1,0,-1},{1,3,4,-2}}; 结果为31
        //结果为12
        double[][] test = {{1,-1,2,-3,1},{-3,3,-7,9,-5},{2,0,4,-2,1},{3,-5,7,-14,6},{4,-4,10,-10,2}};
        double result;
        try {
            result = mathDeterminantCalculation(test);
            System.out.println(result);
        } catch (Exception e) {
            System.out.println("不是正确的行列式！！");
        }
    }
    */
    
    public static double Det(double[][] Matrix, int N)//计算n阶行列式（N=n-1）
    {
        int T0;
        int T1;
        int T2;
        double Num;
        int Cha;
        double[][] B;
        if (N > 0) {
            Cha = 0;
            B = new double[N][N];
            Num = 0;
            if (N == 1) {
                return Matrix[0][0] * Matrix[1][1] - Matrix[0][1] * Matrix[1][0];
            }
            for (T0 = 0; T0 <= N; T0++)
            {
                for (T1 = 1; T1 <= N; T1++)
                {
                    for (T2 = 0; T2 <= N - 1; T2++)
                    {
                        if (T2 == T0) {
                            Cha = 1;

                        }
                        B[T1 - 1][T2] = Matrix[T1][T2 + Cha];
                    }
                    Cha = 0;
                }
                Num = Num + Matrix[0][T0] * Det(B, N - 1) * Math.pow((-1), T0);
            }
            return Num;
        } else if (N == 0) {
            return Matrix[0][0];
        }

        return 0;

    }






        /***
         * 求行列式的算法
         * @param value 需要算的行列式
         * @return 计算的结果
         */
        public static double mathDeterminantCalculation(double[][] value) throws Exception{
            if (value.length == 1) {
                //当行列式为1阶的时候就直接返回本身
                return value[0][0];
            }else if (value.length == 2) {
                //如果行列式为二阶的时候直接进行计算
                return value[0][0]*value[1][1]-value[0][1]*value[1][0];
            }
            //当行列式的阶数大于2时
            double result = 1;
            for (int i = 0; i < value.length; i++) {
                //检查数组对角线位置的数值是否是0，如果是零则对该数组进行调换，查找到一行不为0的进行调换
                if (value[i][i] == 0) {
                    value = changeDeterminantNoZero(value,i,i);
                    result*=-1;
                }
                for (int j = 0; j <i; j++) {
                    //让开始处理的行的首位为0处理为三角形式
                    //如果要处理的列为0则和自己调换一下位置，这样就省去了计算
                    if (value[i][j] == 0) {
                        continue;
                    }
                    //如果要是要处理的行是0则和上面的一行进行调换
                    if (value[j][j]==0) {
                        double[] temp = value[i];
                        value[i] = value[i-1];
                        value[i-1] = temp;
                        result*=-1;
                        continue;
                    }
                    double  ratio = -(value[i][j]/value[j][j]);
                    value[i] = addValue(value[i],value[j],ratio);
                }
            }
            DecimalFormat df = new DecimalFormat(".##");
            return Double.parseDouble(df.format(mathValue(value,result)));
        }

        /**
         * 计算行列式的结果
         * @param value
         * @return
         */
        public static double mathValue(double[][] value,double result) throws Exception{
            for (int i = 0; i < value.length; i++) {
                //如果对角线上有一个值为0则全部为0，直接返回结果
                if (value[i][i]==0) {
                    return 0;
                }
                result *= value[i][i];
            }
            return result;
        }

        /***
         * 将i行之前的每一行乘以一个系数，使得从i行的第i列之前的数字置换为0
         * @param currentRow 当前要处理的行
         * @param frontRow i行之前的遍历的行
         * @param ratio 要乘以的系数
         * @return 将i行i列之前数字置换为0后的新的行
         */
        public static double[] addValue(double[] currentRow,double[] frontRow, double ratio)throws Exception{
            for (int i = 0; i < currentRow.length; i++) {
                currentRow[i] += frontRow[i]*ratio;
            }
            return currentRow;
        }

        /**
         * 指定列的位置是否为0，查找第一个不为0的位置的行进行位置调换，如果没有则返回原来的值
         * @param determinant 需要处理的行列式
         * @param line 要调换的行
         * @param row 要判断的列
         */
        public static double[][] changeDeterminantNoZero(double[][] determinant,int line,int row)throws Exception{
            for (int j = line; j < determinant.length; j++) {
                //进行行调换
                if (determinant[j][row] != 0) {
                    double[] temp = determinant[line];
                    determinant[line] = determinant[j];
                    determinant[j] = temp;
                    return determinant;
                }
            }
            return determinant;
        }


}
