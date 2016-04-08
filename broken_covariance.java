package matrixanalysis;
import java.io.FileNotFoundException;
import java.lang.Double;

public class broken_covariance {
	// Takes in a N x M matrix and produces an M x M matrix
	//
	//
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
		double mean1 = 0; double mean2 = 0; int count = 0;
		for(int i=0; i<N; ++i) {
			if(miss[i] > 0)
				continue;
			else {
				count += 1;
				mean1 += vec1.data[0][i];
				mean2 += vec2.data[0][i];
				entry += vec1.data[0][i] * vec2.data[0][i];
			}
		}
		entry /= count;
		entry -= (mean1 / count) * (mean2 / count);
		return entry;
	}
	
	public void fix() {
		Matrix A = this.A;
		Matrix V = this.V;
		int M = this.M;
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
		try {
			CSVread read = new CSVread();
			//pass the file path, are there colnames? are there rownames?
			Matrix use = read.go("/home/nikwoj/workspace/matrix-analysis-2016/matrixanalysis/DSGaggregate.csv",true,true);
			use = use.T();
			int N = use.M();
			
//			for(int i=0; i<N; ++i) {
//				for(int j=0; j<M; ++j) {
//					if(Math.random() < .5) {
//						use.data[j][i] = Double.NaN;
//					}
//				}
//			}
			
//			double mean = 0; double count = 0;
//			for(int i=0; i<M; ++i) {
//				mean = 0; count = 0;
//				for(int j=0; j<N; ++j) {
//					if(!Double.isNaN(use.data[i][j]))
//						mean += use.data[i][j];
//						count += 1;
//				}
//				mean *= 1 / count;
//				for(int j=0; j<N; ++j) {
//					use.data[i][j] -= mean;
//				}
//			}
			

			broken_covariance A = new broken_covariance(use);
			A.fix();
			
			System.out.println("\n\n");
			A.V.show();
			Matrix V = A.V;
			System.out.println("\n\n");

			JacRot B = new JacRot(A.V);
			
			B.factor(2000, 0.00001);
			Matrix dd = B.A;
			
			dd.show();
			System.out.println("\n\n");
			
			
			for(int i=0; i<N; ++i) {
				if(dd.data[i][i] < 100) { 
					dd.data[i][i] = 100.0; 
				}
			}
			
			Matrix AA = B.D.T().times(dd.times(B.D));
			Matrix C = AA.plus(V.scale(-1));
			
			double sum = 0;
			for(int i=0; i<N; ++i) {
				for(int j=0; j<N; ++j) 
					sum += Math.abs(C.data[i][j]);
			}
			System.out.println(sum);
			C.show();
			
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
