package org.example.sandbox;

import java.text.DecimalFormat;

public class Lab2 {
    private static final DecimalFormat df = new DecimalFormat("0.00");

    public static void main(String[] args) {
        System.out.println("gaussMethod");
        printResult(gaussMethod(
                new double[][]{
                        {7, 2, 3, 0, 32},
                        {0, 3, 2, 6, 47},
                        {2, 5, 1, 0, 23},
                        {0, 1, 4, 2, 29}
                }));

        System.out.println("sqrtMethod");
        printResult(sqrtMethod(
                new double[][]{
                        {1, 2, 0},
                        {2, 2, 3},
                        {0, 3, 2}
                },
                new double[]{8, 21, 17}));

        System.out.println("simpleIterationMethod");
        printResult(simpleIterationMethod(
                new double[][]{
                        {3, 0, 0, 1},
                        {0, 6, 2, 0},
                        {0, 2, 3, 0},
                        {1, 0, 0, 4}
                },
                new double[]{7, 18, 13, 17}));

        System.out.println("zedellMethod");
        printResult(zedellMethod(
                new double[][]{
                        {6, 0, 2, 3},
                        {0, 4, 2, 1},
                        {2, 2, 5, 0},
                        {1, 1, 0, 3}
                },
                new double[]{24, 18, 21, 15}));
    }

    private static double[] gaussMethod(double[][] A_matrix) {
        for (int i = 0; i < A_matrix.length; i++) {
            for (int j = i + 1; j < A_matrix.length + 1; j++) {
                A_matrix[i][j] = A_matrix[i][j] / A_matrix[i][i];
                for (int k = i + 1; k < A_matrix.length; k++) {
                    A_matrix[k][j] = A_matrix[k][j] - A_matrix[k][i] * A_matrix[i][j];
                }
            }
            A_matrix[i][i] = 1;
            for (int k = i + 1; k < A_matrix.length; k++) {
                A_matrix[k][i] = 0;
            }
            System.out.println("Step = " + (i + 1));
            printM(A_matrix);
            System.out.println();
        }
        double[] result = new double[A_matrix.length];

        for (int i = A_matrix.length - 1; i >= 0; i--) {
            result[i] = A_matrix[i][A_matrix.length];
            for (int j = i + 1; j < A_matrix.length; j++) {
                result[i] -= A_matrix[i][j] * result[j];
            }
        }
        return result;
    }

    private static double[] sqrtMethod(double[][] a, double[] b) {
        double[][] d = new double[a.length][a[0].length];
        double[][] s = new double[a.length][a[0].length];
        for (int i = 0; i < a.length; i++) {
            double sum = 0;
            for (int p = 0; p < i; p++) {
                sum += s[p][i] * s[p][i] * d[p][p];
            }
            d[i][i] = Math.signum(a[i][i] - sum);
            s[i][i] = Math.sqrt(Math.abs(a[i][i] - sum));
            for (int j = i + 1; j < a[0].length; j++) {
                double sum2 = 0;
                for (int p = 0; p < i; p++) {
                    sum += s[p][i] * d[p][p] * s[p][j];
                }
                s[i][j] = (a[i][j] - sum2) / (d[i][i] * s[i][i]);
            }
        }

        System.out.println("S:");
        printM(s);
        System.out.println("D:");
        printM(d);

        double[][] std = multiplyM(getTranspose(s), d);
        double[][] stdPlus = new double[b.length][b.length + 1];
        for (int i = 0; i < b.length; i++) {
            System.arraycopy(std[i], 0, stdPlus[i], 0, b.length);
            stdPlus[i][b.length] = b[i];
        }

        for (int i = 0; i < stdPlus.length; i++) {
            double mod = stdPlus[i][i];
            for (int j = 0; j < stdPlus.length + 1; j++) {
                stdPlus[i][j] /= mod;
            }
        }
        double[] y = gaussMethod(stdPlus);
        System.out.println("Y:");
        for (double v : y) {
            System.out.print(df.format(v) + "\t");
        }
        System.out.println();
        System.out.println("result:");

        double[][] sPlus = new double[b.length][b.length + 1];
        for (int i = 0; i < b.length; i++) {
            System.arraycopy(s[i], 0, sPlus[i], 0, b.length);
            sPlus[i][b.length] = y[i];
        }
        return gaussMethod(sPlus);
    }

    private static double[] simpleIterationMethod(double[][] a_matrix, double[] b) {
        double[] result = new double[b.length];
        int n = 1;
        while (n <= 10) {
            double[] result1 = new double[result.length];
            System.arraycopy(result, 0, result1, 0, result1.length);
            for (int i = 0; i < b.length; i++) {
                result1[i] = b[i];
                for (int j = 0; j < b.length; j++)
                    if (j != i)
                        result1[i] -= a_matrix[i][j] * result[j];
                result1[i] /= a_matrix[i][i];
            }
            result = result1;
            System.out.println("Step " + n + ":");
            for (int i = 0; i < result.length; i++)
                System.out.println("X" + (i + 1) + " = " + df.format(result[i]) + ";");
            n++;
        }
        return result;
    }

    private static double[] zedellMethod(double[][] a_matrix, double[] b) {
        double[] result = new double[b.length];
        int n = 1;
        while (n <= 10) {
            for (int i = 0; i < b.length; i++) {
                result[i] = b[i];
                for (int j = 0; j < b.length; j++)
                    if (j != i)
                        result[i] -= a_matrix[i][j] * result[j];
                result[i] /= a_matrix[i][i];
            }
            System.out.println("Step " + n + ":");
            for (int i = 0; i < result.length; i++)
                System.out.println("X" + (i + 1) + " = " + df.format(result[i]) + ";");
            n++;
        }
        return result;
    }

    private static double[][] multiplyM(double[][] a, double[][] b) {
        double[][] result = new double[a.length][b[0].length];
        for (int row = 0; row < result.length; row++) {
            for (int col = 0; col < result[row].length; col++) {
                result[row][col] = multiplyMatricesCell(a, b, row, col);
            }
        }
        return result;
    }

    private static double multiplyMatricesCell(double[][] firstMatrix, double[][] secondMatrix, int row, int col) {
        double cell = 0;
        for (int i = 0; i < secondMatrix.length; i++) {
            cell += firstMatrix[row][i] * secondMatrix[i][col];
        }
        return cell;
    }

    public static double[][] getTranspose(double[][] a) {
        double[][] transpose = new double[a.length][a[0].length];
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < a[i].length; j++) {
                transpose[i][j] = a[j][i];
            }
        }
        return transpose;
    }

    private static void printResult(double[] doubles) {
        for (double aDouble : doubles) {
            System.out.println(df.format(aDouble));
        }
        System.out.println("---------------------");
    }

    private static void printM(double[][] a_matrix) {
        for (double[] line : a_matrix) {
            for (double val : line) {
                System.out.print(df.format(val) + "\t");
            }
            System.out.println();
        }
    }
}
