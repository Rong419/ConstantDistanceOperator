package consoperators;

import java.text.DecimalFormat;


public class JacobianMatrixDeterminant {

    /*
    compute the determinant of n*n matrix
    N = n - 1
    Matrix is the input Jacobian matrix
     */
    public static double Determinant(double[][] Matrix, int N)
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
                Num = Num + Matrix[0][T0] * Determinant(B, N - 1) * Math.pow((-1), T0);
            }
            return Num;
        } else if (N == 0) {
            return Matrix[0][0];
        }

        return 0;

    }
        /*test the calculation

    public static void main(String[] args) {
//		double[][] test = {{2,1,-1},{4,-1,1},{201,102,-99}}; 			"return -18"
//		double[][] test = {{1,1,-1,3},{-1,-1,2,1},{2,5,2,4},{1,2,3,2}};  "return 33"
//		double[][] test = {{1,0,-1,2},{-2,1,3,1},{0,1,0,-1},{1,3,4,-2}}; "return 31"
        //"return 12"
        double[][] test = {{1,-1,2,-3,1},{-3,3,-7,9,-5},{2,0,4,-2,1},{3,-5,7,-14,6},{4,-4,10,-10,2}};
        double result;
        try {
            result = mathDeterminantCalculation(test);
            System.out.println(result);
        } catch (Exception e) {
            System.out.println("Wrong answer！！");
        }
    }
    */

        /*
         * Algorithm to compute the determinant
         * @param value: the matrix
         * @return: the output
         */

        public static double mathDeterminantCalculation(double[][] value) throws Exception{
            if (value.length == 1) {
                return value[0][0];
            }else if (value.length == 2) {
                return value[0][0]*value[1][1]-value[0][1]*value[1][0];
            }
            //for matrix size is larger than 2
            double result = 1;
            for (int i = 0; i < value.length; i++) {
                if (value[i][i] == 0) {
                    value = changeDeterminantNoZero(value,i,i);
                    result*=-1;
                }
                for (int j = 0; j <i; j++) {
                    if (value[i][j] == 0) {
                        continue;
                    }
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
         * Compute the result of determinant
         * @param value
         * @return
         */

        public static double mathValue(double[][] value,double result) throws Exception{
            for (int i = 0; i < value.length; i++) {
                if (value[i][i]==0) {
                    return 0;
                }
                result *= value[i][i];
            }
            return result;
        }

        /***
         *
         * @param currentRow
         * @param frontRow
         * @param ratio
         * @return
         */
        public static double[] addValue(double[] currentRow,double[] frontRow, double ratio)throws Exception{
            for (int i = 0; i < currentRow.length; i++) {
                currentRow[i] += frontRow[i]*ratio;
            }
            return currentRow;
        }

        /**
         * @param determinant
         * @param line
         * @param row
         */
        public static double[][] changeDeterminantNoZero(double[][] determinant,int line,int row)throws Exception{
            for (int j = line; j < determinant.length; j++) {
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
