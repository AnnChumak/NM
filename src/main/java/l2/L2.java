package l2;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;

public class L2 {
	static double dy(double x, double y) {
		return (x * x) + (y / x);
	}

	static double solve(double x) {
		return ((x * x * x) / 2) - (x / 2);
	}

	static double[] RungeKutta(double h, double x0, double y0) {
		int n = (int) (1.0 / h);
		double xi, F1, F2, F3, F4;
		double[] y = new double[n + 2];
		y[0] = y0;
		for (int i = 0; i < n + 1; i++) {
			xi = x0 + i * h;
			F1 = dy(xi, y[i]);
			F2 = dy(xi + h / 2, y[i] + h * F1 / 2);
			F3 = dy(xi + h / 2, y[i] + h * F2 / 2);
			F4 = dy(xi + h, y[i] + h * F3);
			y[i + 1] = y[i] + h * (F1 + 2 * F2 + 2 * F3 + F4) / 6;
		}
		return y;
	}

	static double[] Adams(double h, double x0, double y0) {
		int n = (int) (1.0 / h);
		double xi, F1, F2, F3, F4;
		double[] y = new double[n + 2];
//		y.resize(n + 2); ??
		y[0] = y0;
		for (int i = 0; i < 4; i++) {
			xi = x0 + i * h;
			F1 = dy(xi, y[i]);
			F2 = dy(xi + h / 2, y[i] + h * F1 / 2);
			F3 = dy(xi + h / 2, y[i] + h * F2 / 2);
			F4 = dy(xi + h, y[i] + h * F3);
			y[i + 1] = y[i] + h * (F1 + 2 * F2 + 2 * F3 + F4) / 6;
		}
		F1 = dy(x0, y[0]);
		F2 = dy(x0 + h, y[1]);
		F3 = dy(x0 + 2 * h, y[2]);
		F4 = dy(x0 + 3 * h, y[3]);
		for (int i = 4; i < n + 2; i++) {
			y[i] = y[i - 1] + h * (55 * F4 - 59 * F3 + 37 * F2 - 9 * F1) / 24;
			F1 = F2;
			F2 = F3;
			F3 = F4;
			F4 = dy(x0 + i * h, y[i]);
		}
		return y;
	}

	static double max(double[] v, double a, double h) {
		double m = Math.abs(solve(a) - v[0]);
		for (int i = 0; i < v.length; i++)
			if (Math.abs(solve(a + i * h) - v[i]) > m)
				m = Math.abs(solve(a + i * h) - v[i]);

		return m;
	}

	public static void main(String[] args) throws FileNotFoundException {
		File file = new File("output.txt");
		PrintWriter out = new PrintWriter(file.getAbsoluteFile());

		double h = 0.01;
		double x0 = 1.0;

		double[] yr = RungeKutta(h, x0, solve(x0));
		double[] ya = Adams(h, x0, solve(x0));

		int n = (int) (1 / h);

		out.printf("x[i]\t y[i]\t     rk[i]\t     ad[i]\t     |y[i] - rk[i]|\t |y[i] - ad[i]|y[i]%n");
		for (int i = 1; i < n + 2; i++)
			out.printf("%.2f\t %.6f\t %.6f\t %.6f\t %.6f\t     %.6f\t %n",
					x0 + i * h,
					solve(x0 + i * h),
					yr[i],
					ya[i],
					Math.abs(solve(x0 + i * h) - yr[i]),
					Math.abs(solve(x0 + i * h) - ya[i]));

		out.printf("h\t\t     e(ad)\t\t e(rk)%n");
		for (int i = 1; i < 24; i++) {
			n = 10 * i;
			h = 1.0 / n;
			yr = RungeKutta(h, x0, solve(x0));
			ya = Adams(h, x0, solve(x0));
			out.printf("%.6f \t %.6f \t %.6f %n", h, max(yr, x0, h), max(ya, x0, h));
		}
		out.close();
	}
}

