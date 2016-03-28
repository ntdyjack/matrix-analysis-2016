package matrixanalysis;

import java.io.FileNotFoundException;

public class CholDecomp {
	
	private final Matrix L;
	private final Matrix V;
	
	public CholDecomp(Matrix V) {
		this.V = V;
		int N = V.N();
		this.L = new Matrix(N,N);
	}
	
	public Matrix L() { return this.L; }
	
	public void factor() {
		double sum = 0;
		int i, j, k;
		int N = this.V.N();
		
		for (i=0; i < N; ++i) { // What row are we looking at
			for (j=i; j < N; ++j) { // What column are we looking at
				sum = 0;
				if (i==j) {	
					for (k=0; k<i; ++k )  // Sum up previous entries
						sum += L.data[i][k] * L.data[i][k];
					L.data[i][i] = Math.sqrt(this.V.data[i][i] - sum);
				}
				else {
					for (k=0; k<i; ++k)
						sum +=  L.data[i][k] * L.data[j][k];
					L.data[j][i] = (this.V.data[i][j] - sum) / L.data[i][i];
				}
			}
		}
	}
	
	public static void main(String[] args) throws FileNotFoundException {
		Matrix A = Matrix.random(5,5);
		A = A.times(A.T());
		//A.show();
		//System.out.println();
		
		CholDecomp CD = new CholDecomp(A);
		CD.factor();
		//Matrix L = CD.L();
		//L.show();
		//System.out.println();
		//L.times(L.T()).show();
		//String "matrixanalysis/iris.csv";
		//CholDecomp CD = new CholDecomp(A);

	}
}
