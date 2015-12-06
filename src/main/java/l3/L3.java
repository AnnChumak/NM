package l3;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;

public class L3 {
	// Наша функция
	static double y(double x) {
		return a * x * x + b * x + c + 1 / (d * x + e);
	}

	// Первая производная
	static double dy(double x) {
		return 2 * a * x + b - d / Math.pow(d * x + e, 2);
	}

	// Вторая производная
	static double ddy(double x) {
		return 2 * a + 2 * d * d / Math.pow(d * x + e, 3);
	}

	// A y
	static double f(double x) {
		return ddy(x) + p(x) * dy(x) + q(x) * y(x);
	}

	public static void tridiagonal(int n, double[] d1, double[] d2, double[] d3, double[] fa, double[] x) {
		double m;
		for (int i = 1; i < n; i++) {
			m = d1[i] / d2[i - 1];
			d2[i] -= m * d3[i - 1];
			fa[i] -= m * fa[i - 1];
		}
		x[n - 1] = fa[n - 1] / d2[n - 1];

		for (int i = n - 2; i >= 0; i--)
			x[i] = (fa[i] - d3[i] * x[i + 1]) / d2[i];
	}

	private static double a = -1.147,
			b = 2.642,
			c = -2.894,
			d = -1.599,
			e = 2.093;

	// Коэфициент при первой производной
	static double p(double x) {
		return 1 / x;
	}

	// Коэфициент при второй
	static double q(double x) {
		return 2.0;
	}

	public static void main(String[] args) throws FileNotFoundException {
		int n = 100;
		double x1 = 0.7, x2 = 1;
		double h = (x2 - x1) / n;
		double c1 = 2,
				c2 = 3,
				a1 = y(x2),
				a2 = c1 * dy(x1) + c2 * y(x1);

		double[] pa = new double[n + 1];
		double[] qa = new double[n + 1];
		double[] fa = new double[n + 1];
		double[] ya = new double[n + 1];
		double A[][] = new double[n + 1][n + 1];

		for (int i = 0; i <= n; i++)
			for (int j = 0; j <= n; j++)
				A[i][j] = 0;

		for (int i = 0; i <= n; i++) {
			ya[i] = y(x1 + i * h);
			pa[i] = p(x1 + i * h);
			qa[i] = q(x1 + i * h);
			fa[i] = f(x1 + i * h);

			if (i == 0) {
				A[i][i] = qa[i] - 2 / (h * h) + 2 * h * c2 * (1 / (h * h) - pa[0] / (2 * h)) / c1;
				A[i][i + 1] = 2 / (h * h);
				fa[i] += 2 * h * a2 * (1 / (h * h) - pa[i] / (2 * h)) / c1;
			} else if (i == n) {
				A[i][i - 1] = 0;
				A[i][i] = 1;
				fa[i] = a1;
			} else {
				A[i][i - 1] = 1 / (h * h) - pa[i] / (2 * h);
				A[i][i] = -2 / (h * h) + qa[i];
				A[i][i + 1] = 1 / (h * h) + pa[i] / (2 * h);
			}
		}

		double[] d1 = new double[n + 1];
		double[] d2 = new double[n + 1];
		double[] d3 = new double[n + 1];
		double[] x = new double[n + 1];

		for (int i = 0; i <= n; i++) {
			if (i == 0)
				d1[i] = 0;
			else
				d1[i] = A[i][i - 1];

			if (i == n)
				d3[i] = 0;
			else
				d3[i] = A[i][i + 1];

			d2[i] = A[i][i];
		}

		tridiagonal(n + 1, d1, d2, d3, fa, x);

		// Выводим отчет
		PrintWriter out = new PrintWriter(new File("out.txt").getAbsoluteFile());
		double max = Math.abs(ya[0] - x[0]);
		out.printf("x\t\t\t\t  y(x)\t\t\t  Approx\t\t e(x)%n");

		for (int i = 0; i <= n; i++) {

			if (max < Math.abs(ya[i] - x[i]))
				max = Math.abs(ya[i] - x[i]);

			out.printf("%.10f \t %.10f \t %.10f \t %.10f%n",
					x1 + i * h,
					ya[i],
					x[i],
					Math.abs(ya[i] - x[i]));
		}

		out.printf("%n");
		out.printf("%ne(h) = %.10f%n", max);
		out.close();
	}
}
