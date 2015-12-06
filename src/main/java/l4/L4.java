package l4;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;

import static java.lang.Math.*;

public class L4 {

	static double s = 1.0 / 4;
	static double L = 1;

	static double u0(double x) {
		return x * cos(PI * x * 0.5);
	}

	static double real(double t, double x) {
		return u0(x) * cos(PI * t);
	}

	static double utt(double t, double x) {
		return (-PI * PI * u0(x) * cos(PI * t));
	}

	static double uxx(double t, double x) {
		return -2 * PI * sin(PI * x) * cos(PI * t) - x * PI * PI * cos(PI * x) * cos(PI * t);
	}

	static double ut0(double x) {
		return 0;
	}

	static double u1(double t) {
		return real(t, 0);
	}

	static double u2(double t) {
		return real(t, L);
	}

	static double f(double t, double x) {
		return utt(t, x) - uxx(t, x);
	}

	static void SLAR(double[] a1, double[] a2, double[] a3, double[] mf, double[] y, int n) {
		n++;
		double[] b1;
		double[] b2;
		b1 = new double[n];
		b2 = new double[n];
		b1[1] = -a2[0] / a3[0];
		b2[1] = mf[0] / a3[0];

		for (int i = 2; i < n; i++) {
			b1[i] = -a2[i - 1] / (a1[i - 1] * b1[i - 1] + a3[i - 1]);
			b2[i] = (mf[i - 1] - a1[i - 1] * b2[i - 1]) / (a1[i - 1] * b1[i - 1] + a3[i - 1]);
		}

		y[n - 1] = (mf[n - 1] - a1[n - 1] * b2[n - 1]) / (a1[n - 1] * b1[n - 1] + a3[n - 1]);
		for (int i = n - 2; i > -1; i--)
			y[i] = b1[i + 1] * y[i + 1] + b2[i + 1];
	}

	static void pde(int k, PrintWriter file, PrintWriter file1, PrintWriter file2, PrintWriter file3,
	                PrintWriter file4, int c) {
		int x1 = 0, x2 = (int) L, n = 10;
		double h = (double) (x2 - x1) / n;
		double[] mf = new double[n + 1];
		double[] uj1 = new double[n + 2];
		double[] uj2 = new double[n + 2];
		double[] uj3 = new double[n + 2];
		double[] a1 = new double[n + 1];
		double[] a2 = new double[n + 1];
		double[] a3 = new double[n + 1];
		double z = (double) 2 / k;
		//double y_n = new double[n + 1];
		//double uy_n = new double[n + 1];
		//double y = new double[n + 1];
		//double uy = new double[n + 1];

		k++;

		for (int i = 0; i < n + 1; i++)
			uj1[i] = u0(x1 + i * h);

		uj2[0] = u1(z);
		uj2[n] = u2(z);

		for (int i = 1; i < n; i++)
			uj2[i] = uj1[i] + z * ut0(x1 + i * h) + pow(z, 2) * ((uj1[i + 1] - 2 * uj1[i] + uj1[i - 1]) / pow(h, 2) + f(0, x1 + i * h)) / 2;

		a1[0] = a2[0] = a1[n] = a2[n] = 0;
		a3[0] = a3[n] = 1;
		double le = 0, ly = 0;

		if (c == 1) {
			file.printf("t= 0%n");
			for (int i = 0; i < n + 1; i++) {
				file.printf("u[0,%.4f]= %.4f \tu_real= %.4f", i * h, uj1[i], real(0, i * h));
				file.printf("\tu-y= %.4f%n", abs(real(0, i * h) - uj1[i]));
				if (uj1[i] > ly) ly = uj1[i];
				if (abs(real(0, i * h) - uj1[i]) > le) le = (abs(real(0, i * h) - uj1[i]));
			}
			file.println();
			file3.println(ly);
			file3.println(le);

			ly = le = 0;
			file.printf("t= %.4f%n", z);
			for (int i = 0; i < n + 1; i++) {
				file.printf("u[%.4f,%.4f]= %.4f \tu_real= %.4f", z, i * h, uj2[i], real(z, i * h));
				file.printf("\tu-y= %.4f%n", abs(real(z, i * h) - uj2[i]));
				if (uj2[i] > ly) ly = uj2[i];
				if (abs(real(z, i * h) - uj2[i]) > le) le = (abs(real(z, i * h) - uj2[i]));
			}
			file.println();
			file3.println(ly);
			file3.println(le);

		}

		for (int j = 2; j < k; j++) {
			ly = le = 0;
			mf[0] = u1(z * j);
			mf[n] = u2(z * j);

			for (int i = 1; i < n; i++) {
				a1[i] = a2[i] = s * pow(h, -2);
				a3[i] = -2 * s * pow(h, -2) - pow(z, -2);
				mf[i] = (-2 * uj2[i] + uj1[i]) * pow(z, -2) - ((1 - 2 * s) * (uj2[i + 1] - 2 * uj2[i] + uj2[i - 1]) + s * (uj1[i + 1] - 2 * uj1[i] + uj1[i - 1])) * pow(h, -2) - s * f(z * j, x1 + i * h) - (1 - 2 * s) * f(z * (j - 1), x1 + i * h) - s * f(z * (j - 2), x1 + i * h);
			}

			SLAR(a1, a2, a3, mf, uj3, n);

			if (c == 1) {
				file.printf("t= %.4f%n", z * j);
				for (int i = 0; i < n + 1; i++) {
					file.printf("u[%.4f,%.4f]= %.4f \tu_real= %.4f", z * j, i * h, uj3[i], real(z * j, i * h));
					file.printf("\tu-y= %.4f%n", abs(real(z * j, i * h) - uj3[i]));
					if (j == k - 1) {
						file1.println(abs(real(z * j, i * h) - uj3[i]));
						file2.println(uj3[i]);
					}
					if (uj3[i] > ly) ly = uj3[i];
					if (abs(real(z * j, i * h) - uj3[i]) > le) le = (abs(real(z * j, i * h) - uj3[i]));
				}
				file.println();
				file3.println(ly);
				file3.println(le);
			}

			for (int b = 0; b < n + 2; b++) {
				uj1[b] = uj2[b];
				uj2[b] = uj3[b];
				uj3[b] = 0;
			}
		}
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	public static void main(String[] args) throws FileNotFoundException {
		int k = 40;
		PrintWriter f = new PrintWriter(new File("out.txt").getAbsoluteFile());
		PrintWriter f1 = new PrintWriter(new File("(u-y)^(k).txt").getAbsoluteFile());
		PrintWriter f2 = new PrintWriter(new File("y^(k).txt").getAbsoluteFile());
		PrintWriter f3 = new PrintWriter(new File("(u-y)(t).txt").getAbsoluteFile());
		PrintWriter f4 = new PrintWriter(new File("y(t).txt").getAbsoluteFile());


		pde(k, f, f1, f2, f3, f4, 1);
		f.close();
		f1.close();
		f2.close();
		f3.close();
		f4.close();
	}
}
