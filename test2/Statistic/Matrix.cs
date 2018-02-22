using System;

namespace Statistic
{
    #region "Exception in the Library"
    class MatrixLibraryExceptions : ApplicationException
    { public MatrixLibraryExceptions(string message) : base(message) { } }

    // The Exceptions in this Class
    class MatrixNullException : ApplicationException
    {
        public MatrixNullException() :
            base("To do this operation, matrix can not be null")
        { }
    }
    class MatrixDimensionException : ApplicationException
    {
        public MatrixDimensionException() :
            base("Dimension of the two matrices not suitable for this operation !")
        { }
    }
    class MatrixNotSquare : ApplicationException
    {
        public MatrixNotSquare() :
            base("To do this operation, matrix must be a square matrix !")
        { }
    }
    class MatrixDeterminentZero : ApplicationException
    {
        public MatrixDeterminentZero() :
            base("Determinent of matrix equals zero, inverse can't be found !")
        { }
    }
    class VectorDimensionException : ApplicationException
    {
        public VectorDimensionException() :
            base("Dimension of matrix must be [3 , 1] to do this operation !")
        { }
    }
    class MatrixSingularException : ApplicationException
    {
        public MatrixSingularException() :
            base("Matrix is singular this operation cannot continue !")
        { }
    }
    #endregion
    public class Matrix
    {

        static public double Determinant(double[,] a)
        {
            var matrix1 = new MatrixLib.Matrix(a);
            double det1 = matrix1.Determinant;
            return det1;
        }

        static public double Sum(double[,] InputArray)
        {
            double sum = 0;
            for (int i = 0; i < InputArray.GetLength(0); i++)
            {
                for (int j = 0; j < InputArray.GetLength(1); j++)
                {
                    sum += InputArray[i, j];
                }
            }
            return sum;
        }

        static public double[,] SubMatrix(double[,] A, int row, int col)
        {
            double[,] rez = new double[A.GetLength(0) - 1, A.GetLength(1) - 1];
            for (int i1 = 0, i2 = 0; i1 < A.GetLength(0); i1++, i2++)
            {
                if (i1 == row)
                {
                    i2--;
                    continue;
                }
                for (int i3 = 0, i4 = 0; i3 < A.GetLength(1); i3++)
                {
                    if (i3 != col)
                    {
                        rez[i2, i4++] = A[i1, i3];
                    }
                }
            }
            return rez;
        }
        static public double[,] SwapRowsret(double[,] A, int row1, int row2)
        {
            double.TryParse("1+2", out double rez);
            double[,] Arez = new double[A.GetLength(0), A.GetLength(1)];
            for (int row = 0; row < A.GetLength(0); row++)
                for (int col = 0; col < A.GetLength(1); col++)
                    Arez[row, col] = A[row, col];
            for (int col = 0; col < A.GetLength(1); col++)
            {
                double swap = Arez[row1, col];
                Arez[row1, col] = Arez[row2, col];
                Arez[row2, col] = swap;
            }
            return Arez;
        }
        static public double[,] MultiplicMatrix(double[,] A, double[,] B)
        {
            double[,] rez = new double[A.GetLength(0), B.GetLength(1)];
            for (int i1 = 0; i1 < A.GetLength(0); i1++)
            {
                for (int i2 = 0; i2 < B.GetLength(1); i2++)
                {
                    double hlpdoubl = 0;
                    for (int i3 = 0; i3 < A.GetLength(1); i3++)
                        hlpdoubl += A[i1, i3] * B[i3, i2];
                    rez[i1, i2] = hlpdoubl;
                }
            }
            return rez;
        }
        static public double[,] Abs(double[,] A)
        {
            double[,] rez = new double[A.GetLength(0), A.GetLength(1)];
            for (int i1 = 0; i1 < A.GetLength(0); i1++)
                for (int i2 = 0; i2 < A.GetLength(1); i2++)
                    rez[i1, i2] = Math.Abs(A[i1, i2]);
            return rez;
        }
        static public double[,] TranspMatrix(double[,] A)
        {
            double[,] rez = new double[A.GetLength(1), A.GetLength(0)];
            for (int i1 = 0; i1 < A.GetLength(1); i1++)
                for (int i2 = 0; i2 < A.GetLength(0); i2++)
                    rez[i1, i2] = A[i2, i1];
            return rez;
        }
        static public double[,] InverseMatrix(double[,] A)
        {
            double[,] M = new double[A.GetLength(0), A.GetLength(1)];
            for (int col = 0; col < A.GetLength(1); col++)
                for (int row = 0; row < A.GetLength(0); row++)
                    M[row, col] = Math.Pow(-1, row + col) * Determinant(SubMatrix(A, row, col));
            M = TranspMatrix(M);
            double det = Determinant(A);
            M = MultiplicNumber(M, 1 / det);
            return M;
        }

        static public double[,] MultiplicNumber(double[,] A, double Num)
        {
            double[,] rez = new double[A.GetLength(0), A.GetLength(1)];
            for (int i1 = 0; i1 < A.GetLength(0); i1++)
            {
                for (int i2 = 0; i2 < A.GetLength(1); i2++)
                {
                    rez[i1, i2] = A[i1, i2] * Num;
                }
            }
            return rez;
        }

        static public double[,] Substraction(double[,] A, double[,] B)
        {

            double[,] rez = new double[A.GetLength(0), A.GetLength(1)];
            /*if (B.GetLength(0) == 1)
                for (int i1 = 0; i1 < A.GetLength(0); i1++)
                {
                    for (int i2 = 0; i2 < A.GetLength(1); i2++)
                    {
                        rez[i1, i2] = A[i1, i2] - B[0, i2];
                    }
                }*/
            if (A.GetLength(0) != B.GetLength(0) || A.GetLength(1) != B.GetLength(1))
                return new double[1, 1] { { 0 } };
            for (int i1 = 0; i1 < A.GetLength(0); i1++)
            {
                for (int i2 = 0; i2 < A.GetLength(1); i2++)
                {
                    rez[i1, i2] = A[i1, i2] - B[i1, i2];
                }
            }
            return rez;
        }

        static public double[,] Addition(double[,] A, double[,] B)
        {
            if (A.GetLength(0) != B.GetLength(0) || A.GetLength(1) != B.GetLength(1))
                return new double[1, 1] { { 0 } };
            double[,] rez = new double[A.GetLength(0), A.GetLength(1)];
            for (int i1 = 0; i1 < A.GetLength(0); i1++)
            {
                for (int i2 = 0; i2 < A.GetLength(1); i2++)
                {
                    rez[i1, i2] = A[i1, i2] + B[i1, i2];
                }
            }
            return rez;
        }

        public static void OwnVectors(double[,] startcoefficients, out double[] eightvalue ,out double[,] solution)
        {
            solution = new double[startcoefficients.GetLength(0), startcoefficients.GetLength(1)];
            eightvalue = new double[startcoefficients.GetLength(0)];
            for (int index = 0; index < startcoefficients.GetLength(0); index++)
                solution[index, index] = 1;
            double[,] coefficients  = startcoefficients.Clone() as double[,];
            int result = 1;
            int i, j, k;
            int maxI = 0, maxJ = 0;
            double max, fi, precision = 0.0000001;
            double[][] matricaPoworota;
            int numberOfEquation = coefficients.GetLength(0);
            matricaPoworota = new double[numberOfEquation][];
            for (i = 0; i < numberOfEquation; i++)
            {
                matricaPoworota[i] = new double[numberOfEquation];
            }
            double[][] temp;
            temp = new double[numberOfEquation][];
            for (i = 0; i < numberOfEquation; i++)
            {
                temp[i] = new double[numberOfEquation];
            }
            double fault = 0.0;
            for (i = 0; i < numberOfEquation; i++)
            {
                for (j = i + 1; j < numberOfEquation; j++)
                {
                    fault = fault + coefficients[i,j] * coefficients[i,j];
                }
            }
            fault = Math.Sqrt(2 * fault);
            while (fault > precision)
            {
                max = 0.0;
                for (i = 0; i < numberOfEquation; i++)
                {
                    for (j = i + 1; j < numberOfEquation; j++)
                    {
                        if (coefficients[i,j] > 0 && coefficients[i,j] > max)
                        {
                            max = coefficients[i,j];
                            maxI = i;
                            maxJ = j;
                        }
                        else if (coefficients[i,j] < 0 && -coefficients[i,j] > max)
                        {
                            max = -coefficients[i,j];
                            maxI = i;
                            maxJ = j;
                        }
                    }
                }
                for (i = 0; i < numberOfEquation; i++)
                {
                    for (j = 0; j < numberOfEquation; j++)
                    {
                        matricaPoworota[i][j] = 0;
                    }
                    matricaPoworota[i][i] = 1;
                }
                if (coefficients[maxI, maxI] == coefficients[maxJ,maxJ])
                {
                    matricaPoworota[maxI][maxI] = matricaPoworota[maxJ][maxJ] =
                    matricaPoworota[maxJ][maxI] = Math.Sqrt(2.0) / 2.0;
                    matricaPoworota[maxI][maxJ] = -Math.Sqrt(2.0) / 2.0;
                }
                else
                {
                    fi = 0.5 * Math.Atan((2.0 * coefficients[maxI, maxJ]) /
                    (coefficients[maxI, maxI] - coefficients[maxJ, maxJ]));
                    matricaPoworota[maxI][maxI] = matricaPoworota[maxJ][maxJ] = Math.Cos(fi);
                    matricaPoworota[maxI][maxJ] = -Math.Sin(fi);
                    matricaPoworota[maxJ][maxI] = Math.Sin(fi);
                }
                for (i = 0; i < numberOfEquation; i++)
                {
                    for (j = 0; j < numberOfEquation; j++)
                    {
                        temp[i][j] = 0.0;
                    }
                }
                for (i = 0; i < numberOfEquation; i++)
                {
                    for (j = 0; j < numberOfEquation; j++)
                    {
                        for (k = 0; k < numberOfEquation; k++)
                        {
                            temp[i][j] = temp[i][j] + matricaPoworota[k][i] * coefficients[k, j];
                        }
                    }
                }
                for (i = 0; i < numberOfEquation; i++)
                {
                    for (j = 0; j < numberOfEquation; j++)
                    {
                        coefficients[i,j] = 0.0;
                    }
                }
                for (i = 0; i < numberOfEquation; i++)
                {
                    for (j = 0; j < numberOfEquation; j++)
                    {
                        for (k = 0; k < numberOfEquation; k++)
                        {
                            coefficients[i, j] = coefficients[i,j] +
                            temp[i][k] * matricaPoworota[k][j];
                        }
                    }
                }
                fault = 0.0;
                for (i = 0; i < numberOfEquation; i++)
                {
                    for (j = i + 1; j < numberOfEquation; j++)
                    {
                        fault = fault + coefficients[i, j] * coefficients[i, j];
                    }
                }
                fault = Math.Sqrt(2 * fault);
                for (i = 0; i < numberOfEquation; i++)
                {
                    for (j = 0; j < numberOfEquation; j++)
                    {
                        temp[i][j] = 0.0;
                    }
                }
                for (i = 0; i < numberOfEquation; i++)
                {
                    for (j = 0; j < numberOfEquation; j++)
                    {
                        for (k = 0; k < numberOfEquation; k++)
                        {
                            temp[i][j] = temp[i][j] + solution[i, k] * matricaPoworota[k][j];
                        }
                    }
                }
                for (i = 0; i < numberOfEquation; i++)
                {
                    for (j = 0; j < numberOfEquation; j++)
                    {
                        solution[i,j] = temp[i][j];
                    }
                }
                result++;
            }
            for (int index = 0; index < eightvalue.Length; index++)
                eightvalue[index] = coefficients[index, index];
        }


        public static void SortEighten(ref double[] eightvalues, ref double[,] eightvectors)
        {
            int n = eightvalues.GetLength(0);
            double[][] eightvectors1 = new double[n][];
            for (int i = 0; i < n; i++)
            {
                eightvectors1[i] = new double[n];
                for (int j = 0; j < n; j++)
                    eightvectors1[i][j] = eightvectors[j, i];
            }
            Array.Sort(eightvalues, eightvectors1);
            Array.Reverse(eightvalues);
            Array.Reverse(eightvectors1);
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                    eightvectors[i, j] = eightvectors1[j][i];
            }
            /*
            List<PointF> eightvalues  = new List<PointF>();
            Array.Sort(eightvalues, eightvectors);*/

        }
    }
}
