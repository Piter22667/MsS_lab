package test3;

import org.apache.commons.math3.linear.*;
import java.io.*;
import java.util.*;

public class ParameterIdentification {
    private static final int SYSTEM_SIZE = 6;
    private static final double DT = 0.2;
    private static final double T = 50.0;
    private static final double EPSILON = 0.0001;

    // Клас для параметрів системи
    public static class SystemParameters {
        public double c1, c2, c3, c4, m1, m2, m3;

        public SystemParameters(double c1, double c2, double c3, double c4,
                                double m1, double m2, double m3) {
            this.c1 = c1;
            this.c2 = c2;
            this.c3 = c3;
            this.c4 = c4;
            this.m1 = m1;
            this.m2 = m2;
            this.m3 = m3;
        }

        public double[] getBeta() {
            return new double[]{c1, m1, m2};
        }

        public void updateBeta(double[] deltaBeta) {
            c1 += deltaBeta[0];
            m1 += deltaBeta[1];
            m2 += deltaBeta[2];
        }
    }
        // Створення матриці  на основі рівнянь 5.1-5.3
        private static RealMatrix createMatrixA(SystemParameters params) {
            double[][] A = new double[SYSTEM_SIZE][SYSTEM_SIZE];

            // Перший рядок
            A[0][1] = 1;
            //  (з рівняння 5.1)
            A[1][0] = -(params.c2 + params.c1) / params.m1;
            A[1][2] = params.c2 / params.m1;
            // Третій рядок
            A[2][3] = 1;
            //  (з рівняння 5.2)
            A[3][0] = params.c2 / params.m2;
            A[3][2] = -(params.c2 + params.c3) / params.m2;
            A[3][4] = params.c3 / params.m2;
            // П'ятий рядок
            A[4][5] = 1;
            //  (з рівняння 5.3)
            A[5][2] = params.c3 / params.m3;
            A[5][4] = -(params.c4 + params.c3) / params.m3;
            return MatrixUtils.createRealMatrix(A);
        }

        // Метод Рунге-Кутти 4-го порядку
        private static double[] rungeKutta4(RealMatrix A, double[] y, double dt) {
            int n = y.length;
            double[] k1 = new double[n];
            double[] k2 = new double[n];
            double[] k3 = new double[n];
            double[] k4 = new double[n];
            double[] yTemp = new double[n];

            // k1 = dt * A * y
            k1 = multiplyMatrixVector(A, y);
            for (int i = 0; i < n; i++) k1[i] *= dt;
            // k2 = dt * A * (y + k1/2)
            for (int i = 0; i < n; i++) yTemp[i] = y[i] + k1[i]/2;
            k2 = multiplyMatrixVector(A, yTemp);
            for (int i = 0; i < n; i++) k2[i] *= dt;

            // k3 = dt * A * (y + k2/2)
            for (int i = 0; i < n; i++) yTemp[i] = y[i] + k2[i]/2;
            k3 = multiplyMatrixVector(A, yTemp);
            for (int i = 0; i < n; i++) k3[i] *= dt;

            // k4 = dt * A * (y + k3)
            for (int i = 0; i < n; i++) yTemp[i] = y[i] + k3[i];
            k4 = multiplyMatrixVector(A, yTemp);
            for (int i = 0; i < n; i++) k4[i] *= dt;

            // y(t+dt) = y(t) + (k1 + 2k2 + 2k3 + k4)/6
            double[] result = new double[n];
            for (int i = 0; i < n; i++) {
                result[i] = y[i] + (k1[i] + 2*k2[i] + 2*k3[i] + k4[i])/6;
            }

            return result;
        }

        private static double[] multiplyMatrixVector(RealMatrix A, double[] vector) {
            double[] result = new double[vector.length];
            for (int i = 0; i < A.getRowDimension(); i++) {
                result[i] = 0;
                for (int j = 0; j < A.getColumnDimension(); j++) {
                    result[i] += A.getEntry(i, j) * vector[j];
                }
            }
            return result;
        }


        private static RealMatrix computeDerivativeC1(SystemParameters params) {
            double[][] dA = new double[SYSTEM_SIZE][SYSTEM_SIZE];

            // Похідна по c1 впливає тільки на A[1][0]
            dA[1][0] = -1.0 / params.m1;

            return MatrixUtils.createRealMatrix(dA);
        }

        private static RealMatrix computeDerivativeM1(SystemParameters params) {
            double[][] dA = new double[SYSTEM_SIZE][SYSTEM_SIZE];

            // Похідна по m1 впливає на елементи в рівнянні 5.1
            dA[1][0] = (params.c2 + params.c1) / (params.m1 * params.m1);
            dA[1][2] = -params.c2 / (params.m1 * params.m1);
            return MatrixUtils.createRealMatrix(dA);
        }

        private static RealMatrix computeDerivativeM2(SystemParameters params) {
            double[][] dA = new double[SYSTEM_SIZE][SYSTEM_SIZE];
            // Похідна по m2 впливає на всі елементи в рівнянні 5.2
            dA[3][0] = -params.c2 / (params.m2 * params.m2);
            dA[3][2] = (params.c2 + params.c3) / (params.m2 * params.m2);
            dA[3][4] = -params.c3 / (params.m2 * params.m2);

            return MatrixUtils.createRealMatrix(dA);
        }

        private static double[][] computeSensitivityMatrix(RealMatrix A, double[] y, SystemParameters params) {
            double[][] U = new double[SYSTEM_SIZE][3];
            // Обчислюємо похідні від A*y по кожному параметру
            RealMatrix dC1 = computeDerivativeC1(params);
            RealMatrix dM1 = computeDerivativeM1(params);
            RealMatrix dM2 = computeDerivativeM2(params);

            // Заповнюємо матрицю чутливості
            for (int i = 0; i < SYSTEM_SIZE; i++) {
                U[i][0] = multiplyMatrixVector(dC1, y)[i];
                U[i][1] = multiplyMatrixVector(dM1, y)[i];
                U[i][2] = multiplyMatrixVector(dM2, y)[i];
            }
            return U;
        }


    private static double[][] computeURK(RealMatrix A, double[][] U_prev, double[][] A_deriv_beta, double dt) {
        int n = U_prev.length;
        int m = U_prev[0].length;
        double[][] result = new double[n][m];

        // k1 = dt * (A*U + A_deriv_beta)
        double[][] k1 = new double[n][m];
        double[][] AU = multiplyMatrices(A, U_prev);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                k1[i][j] = dt * (AU[i][j] + A_deriv_beta[i][j]);
            }
        }

        // k2 = dt * (A*(U + k1/2) + A_deriv_beta)
        double[][] k2 = new double[n][m];
        double[][] U_temp = addMatrices(U_prev, multiplyMatrixByScalar(k1, 0.5));
        AU = multiplyMatrices(A, U_temp);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                k2[i][j] = dt * (AU[i][j] + A_deriv_beta[i][j]);
            }
        }

        // k3 = dt * (A*(U + k2/2) + A_deriv_beta)
        double[][] k3 = new double[n][m];
        U_temp = addMatrices(U_prev, multiplyMatrixByScalar(k2, 0.5));
        AU = multiplyMatrices(A, U_temp);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                k3[i][j] = dt * (AU[i][j] + A_deriv_beta[i][j]);
            }
        }

        // k4 = dt * (A*(U + k3) + A_deriv_beta)
        double[][] k4 = new double[n][m];
        U_temp = addMatrices(U_prev, k3);
        AU = multiplyMatrices(A, U_temp);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                k4[i][j] = dt * (AU[i][j] + A_deriv_beta[i][j]);
            }
        }

        // Фінальний результат: U + (k1 + 2k2 + 2k3 + k4)/6
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                result[i][j] = U_prev[i][j] +
                        (k1[i][j] + 2*k2[i][j] + 2*k3[i][j] + k4[i][j])/6;
            }
        }

        return result;
    }




    private static double[][] multiplyMatrices(RealMatrix A, double[][] B) {
        int n = A.getRowDimension();
        int m = B[0].length;
        double[][] result = new double[n][m];

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                result[i][j] = 0;
                for (int k = 0; k < A.getColumnDimension(); k++) {
                    result[i][j] += A.getEntry(i, k) * B[k][j];
                }
            }
        }
        return result;
    }

    private static double[][] addMatrices(double[][] A, double[][] B) {
        int n = A.length;
        int m = A[0].length;
        double[][] result = new double[n][m];

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                result[i][j] = A[i][j] + B[i][j];
            }
        }
        return result;
    }

    private static double[][] multiplyMatrixByScalar(double[][] A, double scalar) {
        int n = A.length;
        int m = A[0].length;
        double[][] result = new double[n][m];

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                result[i][j] = A[i][j] * scalar;
            }
        }
        return result;
    }


    public static Result computeParameters(double[][] measurements, SystemParameters initialParams) {
        List<Double> qualityHistory = new ArrayList<>();
        int maxIterations = 100;
        double[][] calculated = new double[SYSTEM_SIZE][measurements[0].length];

        // Ініціалізація параметрів
        SystemParameters currentParams = new SystemParameters(
                initialParams.c1, initialParams.c2, initialParams.c3, initialParams.c4,
                initialParams.m1, initialParams.m2, initialParams.m3
        );

        // Початкові умови
        double[] y0 = new double[SYSTEM_SIZE];
        for (int i = 0; i < SYSTEM_SIZE; i++) {
            y0[i] = measurements[i][0];
            calculated[i][0] = y0[i];
        }

        // Основний цикл ітерацій
        for (int iter = 1; iter <= maxIterations; iter++) {
            // Ініціалізація для поточної ітерації
            double[][] integralInv = new double[3][3];
            double[][] integralUY = new double[3][1];
            double qualityIndicator = 0.0;
            double[] y = Arrays.copyOf(y0, SYSTEM_SIZE);
            double[][] U = new double[SYSTEM_SIZE][3];
            RealMatrix A = createMatrixA(currentParams);

            // Проходимо по всіх точках вимірювань
            for (int i = 1; i < measurements[0].length; i++) {
                // Розв'язання системи методом Рунге-Кутти
                y = rungeKutta4(A, y, DT);

                for (int j = 0; j < SYSTEM_SIZE; j++) {
                    calculated[j][i] = y[j];
                }

                // Обчислення матриці чутливості
                double[][] A_deriv_beta = computeSensitivityMatrix(A, y, currentParams);
                U = computeURK(A, U, A_deriv_beta, DT);

                // Обчислення різниці між вимірюваннями і розрахунками
                double[] diff = new double[SYSTEM_SIZE];
                for (int j = 0; j < SYSTEM_SIZE; j++) {
                    diff[j] = measurements[j][i] - y[j];
                }

                // Обчислення інтегралів
                for (int k = 0; k < 3; k++) {
                    for (int l = 0; l < 3; l++) {
                        double sum = 0;
                        for (int j = 0; j < SYSTEM_SIZE; j++) {
                            sum += U[j][k] * U[j][l];
                        }
                        integralInv[k][l] += sum * DT;
                    }

                    double sum = 0;
                    for (int j = 0; j < SYSTEM_SIZE; j++) {
                        sum += U[j][k] * diff[j];
                    }
                    integralUY[k][0] += sum * DT;
                }

                // Обчислення показника якості
                for (int j = 0; j < SYSTEM_SIZE; j++) {
                    qualityIndicator += diff[j] * diff[j] * DT;
                }
            }

            qualityHistory.add(qualityIndicator);

            // Виведення проміжних результатів
            System.out.println("Iteration " + iter);
            System.out.println("Quality indicator: " + qualityIndicator);
            System.out.println("Current parameters: c1=" + currentParams.c1 +
                    ", m1=" + currentParams.m1 + ", m2=" + currentParams.m2);

            SystemResponsePlotter.createPlot(measurements, calculated,
                    "system_response_iter_" + iter + ".png");

            // Перевірка умови збіжності
            if (qualityIndicator < EPSILON) {
                ConvergencePlotter.createPlot(qualityHistory, "convergence_plot.png");
                return new Result(currentParams, qualityIndicator, iter);
            }

            try {
                // Розв'язання системи для знаходження beta
                RealMatrix integralInvMatrix = MatrixUtils.createRealMatrix(integralInv);
                RealMatrix integralUYMatrix = MatrixUtils.createRealMatrix(integralUY);

                DecompositionSolver solver = new LUDecomposition(integralInvMatrix).getSolver();
                RealMatrix deltaBeta = solver.solve(integralUYMatrix);

                // Оновлення параметрів
                double[] deltaParams = new double[]{
                        deltaBeta.getEntry(0, 0),
                        deltaBeta.getEntry(1, 0),
                        deltaBeta.getEntry(2, 0)
                };
                currentParams.updateBeta(deltaParams);

            } catch (RuntimeException e) {
                System.out.println("Error solving the system: " + e.getMessage());
                ConvergencePlotter.createPlot(qualityHistory, "convergence_plot.png");
                SystemResponsePlotter.createPlot(measurements, calculated, "system_response_final.png");
                return new Result(currentParams, qualityIndicator, iter);
            }
        }

        ConvergencePlotter.createPlot(qualityHistory, "convergence_plot.png");
        SystemResponsePlotter.createPlot(measurements, calculated, "system_response_final.png");
        return new Result(currentParams, qualityHistory.get(qualityHistory.size() - 1), maxIterations);
    }

        public static void main(String[] args) {
            try {
                double[][] measurements = readDataFromFile("y1.txt");
//                validateData(measurements);
                SystemParameters initialParams = new SystemParameters(
                        0.1, 0.3, 0.2, 0.12, 11.0, 23.0, 18.0);

                Result result = computeParameters(measurements, initialParams);

                System.out.println("\nFinal Results:");
                System.out.println("Iterations: " + result.iterations);
                System.out.println("Quality indicator: " + result.quality);
                System.out.println("Parameters:");
                System.out.println("c1 = " + result.parameters.c1);
                System.out.println("m1 = " + result.parameters.m1);
                System.out.println("m2 = " + result.parameters.m2);

            } catch (Exception e) {
                e.printStackTrace();
            }
        }




    private static double[][] readDataFromFile(String filename) throws IOException {
        try (BufferedReader reader = new BufferedReader(new FileReader(filename))) {
            List<String[]> lines = new ArrayList<>();
            String line;

            while ((line = reader.readLine()) != null) {
                String[] values = line.trim().split("\\s+");
                lines.add(values);
            }

            if (lines.size() != SYSTEM_SIZE) {
                throw new IllegalArgumentException(
                        "Invalid input data: expected " + SYSTEM_SIZE + " rows, got " + lines.size());
            }
            int timePoints = lines.get(0).length;
            double[][] data = new double[SYSTEM_SIZE][timePoints];
            for (int i = 0; i < SYSTEM_SIZE; i++) {
                String[] values = lines.get(i);
                if (values.length != timePoints) {
                    throw new IllegalArgumentException(
                            "Inconsistent data length in row " + (i + 1));
                }
                for (int j = 0; j < timePoints; j++) {
                    try {
                        data[i][j] = Double.parseDouble(values[j]);
                    } catch (NumberFormatException e) {
                        throw new IllegalArgumentException(
                                "Invalid number format at row " + (i + 1) + ", column " + (j + 1));
                    }
                }
            }
            return data;
        }
    }
}