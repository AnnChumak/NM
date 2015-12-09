package l1;

public class Step {
	public static void degree(double[][] a) {

		double[][] y = new double[a.length][1];
		double[][] newY = new double[a.length][1];
		double[][] b = new double[a.length][a.length];
		double[][] I = new double[a.length][a.length];
		double lamda = 0, newlamda = 0, eps = 0.0001;

		for (int i = 0; i < a.length; i++) {
			y[i][0] = 1;//Y0
		}

		newY = mult(a, y); //Y1
		lamda = newY[1][0] / y[1][0]; //ë1

		for (int i = 0; i < a.length; i++) {
			y[i][0] = newY[i][0]; //Y1
		}


		int z = 0;
		do {
			lamda = newlamda;//ë1
			newY = mult(a, y);//Y2
			newlamda = newY[1][0] / y[1][0];//ë2
			newY = norm(newY);
			for (int i = 0; i < a.length; i++) {
				y[i][0] = newY[i][0]; //Y1
			}

		} while (Math.abs(newlamda - lamda) > eps);

		System.out.println("Max eigenvalue = " + newlamda);
		double h = newlamda;
		for (int i = 0; i < newY.length; i++) {
			System.out.println(newY[i][0]);
		}

		for (int i = 0; i < a.length; i++) {
			for (int j = 0; j < a.length; j++) {
				if (i == j) {
					b[i][j] = a[i][j] - newlamda;
				} else {
					b[i][j] = a[i][j];
				}
			}
		}

		for (int i = 0; i < b.length; i++) {
			y[i][0] = 1;//Y0
		}

		newY = mult(b, y); //Y1
		lamda = newY[1][0] / y[1][0]; //ë1

		for (int i = 0; i < b.length; i++) {
			y[i][0] = newY[i][0]; //Y1
		}


		z = 0;
		do {
			lamda = newlamda;//ë1
			newY = mult(b, y);//Y2
			newlamda = newY[1][0] / y[1][0];//ë2
			newY = norm(newY);
			for (int i = 0; i < b.length; i++) {
				y[i][0] = newY[i][0]; //Y1
			}

		} while (Math.abs(newlamda - lamda) > eps);

		System.out.println("Min eigenvalue = " + (h + newlamda));
		for (int i = 0; i < newY.length; i++) {
			System.out.println(newY[i][0]);
		}


	}

	private static double[][] norm(double[][] newY) {
		return newY;
	}

	public static double[][] mult(double[][] a, double[][] b) {
		double s;
		int n = a[0].length;
		int m = b.length;
		double[][] result = new double[n][m];
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				s = 0.0;
				for (int k = 0; k < n; k++) {
					s += a[i][k] * b[k][j];
				}
				result[i][j] = s;
			}
		}
		return result;
	}

}


