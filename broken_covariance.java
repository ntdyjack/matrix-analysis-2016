package matrixanalysis;
import java.io.FileNotFoundException;
import java.lang.Double;
import java.lang.Math;

public class broken_covariance {
	
	private final Matrix A;
	private final int N;
	private final int M;
	public Matrix V;
	
	public broken_covariance(Matrix A) {
		this.A = A;
		this.N = A.N();
		this.M = A.M();
		this.V = new Matrix(M,M);
	}
	
	private int[] find(Matrix vector) {
		int N = this.N;
		int[] miss = new int[N];
		for(int i=0; i<N; ++i) {
			if (Double.isNaN(vector.data[0][i]))
				miss[i] = 1;
		}
		return miss;
	}
	
	private int[] add_missing(int[] miss1, int[] miss2) {
		int N = this.N;
		int[] miss = new int[N];
		for(int i=0; i<N; ++i) {
			miss[i] = miss1[i] + miss2[i];
		}
		return miss;
	}
	
	private double multiply_vec(Matrix vec1, Matrix vec2, int[] miss) {
		double entry=0;
		int N = this.N;
		for(int i=0; i<N; ++i) {
			if(miss[i] > 0)
				continue;
			else
				entry += vec1.data[0][i] * vec2.data[0][i];
		}
		return entry;
	}
	
	public void fix() {
		Matrix A = this.A;
		Matrix V = this.V;
		int N = this.N; int M = this.M;
//		Matrix vec1 = new Matrix(1,N); Matrix vec2 = new Matrix(1,N);
//		int[] miss1 = new int[N]; int[] miss2 = new int[N]; int[] miss = new int[N];
		
		for(int i=0; i<M; ++i) {
			for(int j=0; j<=i; ++j) {
				Matrix vec1 = A.getRow(i); Matrix vec2 = A.getRow(j);
				int[] miss1 = find(vec1); int[] miss2 = find(vec2);
				int[] miss = add_missing(miss1, miss2);
				
//				
//				System.out.println(miss.length);
//				for(int k=0; k<N; ++k) { 
//					System.out.println(miss[k]);
//					System.out.println(k);
//				}
//				
				
				V.data[i][j] = multiply_vec(vec1, vec2, miss);
				V.data[j][i] = V.data[i][j];
			}
		}
		this.V = V;
	}
	
	public static void main(String args[]) {
		CSVread R = new CSVread();
		try {
			Matrix [] Iris = R.go();
			int N = Iris[0].N(); int M = Iris[0].M();
			Matrix use = Iris[1];
			
			for(int i=0; i<N; ++i) {
				for(int j=0; j<M; ++j) {
					if(Math.random() < .1) {
						use.data[j][i] = Double.NaN;
					}
				}
			}
			
			broken_covariance A = new broken_covariance(use);
			A.fix();
			
			System.out.println("\n\n");
			
			A.V.show();
			System.out.println("\n\n");

			JacRot B = new JacRot(A.V);
			
			B.factor(2000, 0.00001);
			B.A.show();
			
			System.out.println("\n\n");
			//Iris[0].show();
			//System.out.println();
			//Iris[1].show();
			//System.out.println();
			//Iris[2].show();
		} catch(FileNotFoundException e) {
			System.out.println("Whoops");
		} finally {
			System.out.println("Wat");
		}
		

	}
	
}
