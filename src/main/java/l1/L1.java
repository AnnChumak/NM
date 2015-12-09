package l1;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.Arrays;

import static java.lang.Math.abs;

public class L1 {
	public static final double EPS = 1e-5;
	private static File log = new File("log1.txt");

	private static double[][] a = {{6.48, 1.04, 1.03, 1.21},
			{1.04, 3.5, 1.3, 0.16},
			{1.03, 1.3, 5.66, 2.1},
			{1.21, 0.16, 2.1, 5.88}};


	private static double[][] a_init = {{6.48, 1.04, 1.03, 1.21},
			{1.04, 3.5, 1.3, 0.16},
			{1.03, 1.3, 5.66, 2.1},
			{1.21, 0.16, 2.1, 5.88}};

	public strictfp static void main(String[] args) throws FileNotFoundException {
		// Поток вывода для лога
		PrintWriter out = new PrintWriter(log.getAbsoluteFile());

		double[][] aNew;    // Переменная для хранения преобразованной матрицы
		double[][] r;       // Матрица поворота
		double[][] h;       // Матрица для нахождения собственных векторов
		{
			h = new double[a.length][a.length];
			for (int i = 0; i < h[0].length; i++) {
				h[i][i] = 1.0;
			}
		}
		double acc;         // Точность
		int n = 0;          // Счетчик количества итераций
		do {
			// Находим наибольший по модулю недиагональный элемент матрицы и его позицию
			int posI = 0;
			int posJ = 1;
			double biggest = abs(a[0][1]);
			for (int i = 0; i < a.length; i++) {
				for (int j = i; j < a[0].length; j++) {
					if (i != j && abs(a[i][j]) > biggest) {
						biggest = abs(a[i][j]);
						posI = i;
						posJ = j;
					}
				}
			}

			// Формируем матрицу поворота
			double[][] m = new double[a.length][a[0].length];
			for (int k = 0; k < m.length; k++) {
				for (int l = 0; l < m.length; l++) {
					if (k == l) {
						m[k][l] = 1;
					} else {
						m[k][l] = 0;
					}
				}
			}
			double fi = 1.0 / 2 * Math.atan((2 * a[posI][posJ]) / (a[posI][posI] - a[posJ][posJ]));
			m[posI][posI] = Math.cos(fi);
			m[posJ][posJ] = m[posI][posI];
			m[posJ][posI] = Math.sin(fi);
			m[posI][posJ] = -m[posJ][posI];
			r = m;

			// Транспонируем матрицу
			double[][] rTr = new double[r[0].length][r.length];
			for (int i = 0; i < r.length; i++) {
				for (int j = 0; j < r[0].length; j++) {
					rTr[j][i] = r[i][j];
				}
			}

			// R(tr) * A * R
			double[][] aN = multiply(rTr, a);
			aNew = multiply(aN, r);

			// Умножаем результат, чтобы в дальнейшем найти собственные вектора
			h = multiply(h, r);

			// Пишем в лог
			out.printf("Iteration: %d%n", n++);
			out.printf("Matrix: %n" + toString(a));
			out.printf("Biggest number: %5.5f; at A[%d, %d]%n", biggest, posI, posJ);
			out.printf("New matrix: %n" + toString(aNew));
			out.print("Diagonal elements: ");
			for (int i = 0; i < a[0].length; i++) {
				out.printf("%5.5f ", a[i][i]);
			}
			out.println();
			out.printf("Quotients: %5.5f, %5.5f, %5.5f, %5.5f ---- %5.5f %n",
					m[posI][posI], m[posJ][posJ], m[posJ][posI], m[posI][posJ],
					m[posI][posI] * m[posI][posI] + m[posJ][posI] * m[posJ][posI]);
			double d = 0.0;
			double w = 0.0;
			for (int i = 0; i < a.length; i++) {
				for (int j = 0; j < a.length; j++) {
					if (i == j) {
						d += Math.pow(a[i][j], 2);
					} else {
						w += Math.pow(a[i][j], 2);
					}
				}
			}
			out.printf("Diagonal = %5.5f, NonDiagonal = %5.5f, Sum = %5.5f", d, w, d + w);
			out.printf("%n%n");

			// Перезаписываем матрицу для следующей итерации
			a = aNew;

			// Считаем погрешность
			acc = 0.0;
			for (int i1 = 0; i1 < a.length; i1++) {
				for (int j = i1; j < a[0].length; j++) {
					if (i1 != j) acc += Math.sqrt(a[i1][j] * a[i1][j]);
				}
			}
		} while (acc > EPS);
		out.close();

		double[][] h_res = new double[a.length][a.length];
		for (int i = 0; i < h.length; i++) {
			for (int j = 0; j < h[0].length; j++) {
				h_res[j][i] = h[i][j];
			}
		}
		h = h_res;

		// Выводим результаты пользователю
		System.out.printf("Eigen values: \t Eigen vectors:%n");
		for (int i = 0; i < a.length; i++) {
			System.out.printf("%5.5f \t\t [", a[i][i]);
			for (int j = 0; j < a.length; j++) {
				System.out.printf("%5.5f ", h[j][i]);
			}
			System.out.printf("]%n");
		}
		System.out.printf("%n");

		System.out.println(Arrays.toString(error(a, a_init, h, 0)));
		System.out.println(Arrays.toString(error(a, a_init, h, 1)));
		System.out.println(Arrays.toString(error(a, a_init, h, 2)));
		System.out.println(Arrays.toString(error(a, a_init, h, 3)));

		System.out.println();
		Step.degree(a_init);
	}

	static double[][] multiply(double[][] a, double[][] b) {
		double s;
		int n = a.length;
		double[][] result = new double[n][n];
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				s = 0;
				for (int k = 0; k < n; k++) {
					s += a[i][k] * b[k][j];
				}
				result[i][j] = s;
			}
		}
		return result;
	}

	static String toString(double[][] a) {
		StringBuilder sb = new StringBuilder();
		for (double[] anA : a) {
			sb.append("[ ");
			for (int j = 0; j < a[0].length; j++) {
				sb.append(String.format("%5.5f", anA[j]));
				if (j != a[0].length - 1) {
					sb.append(", ");
				}
			}
			sb.append(" ]\n");
		}
		return sb.toString();
	}

	static double[] error(double[][] a, double[][] a_init, double[][] h, int k) {
		double[] result = new double[a.length];
		double[] v = new double[a.length];
		for (int i = 0; i < h.length; i++) {
			v[i] = h[k][i];
		}
		double[] av = mult(a_init, v);
		double[] lv = new double[a.length];
		for (int i = 0; i < a.length; i++) {
			lv[i] = v[i] * a[k][k];
		}
		for (int i = 0; i < a.length; i++) {
			result[i] = av[i] - lv[i];
		}

		return result;
	}

	public static double[] mult(double[][] a, double[] b) {
		double s;
		int n = a.length;
		double[] result = new double[n];
		for (int i = 0; i < n; i++) {
			s = 0;
			for (int j = 0; j < n; j++) {
				s = s + a[i][j] * b[j];
			}
			result[i] = s;
		}
		return result;
	}
}