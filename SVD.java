package matrixanalysis;
import java.lang.Math;

public class SVD {
	
	private final Matrix A;
	private  Matrix U;
	private  Matrix B;
	private  Matrix V;
	
	private final int N;
	private final int M;
	
	public SVD(Matrix A) {
		this.A = A;
		this.N = A.N();
		this.M = A.M();
		this.U = Matrix.identity(M);
		this.B = this.A;
		this.V = Matrix.identity(N);
	};
	
	public Matrix [] factor(){
		Matrix [] UBV = new Matrix [3];
		
    	SVD L = new SVD(A);
    	L.bidiagonalize();
		
		
		UBV[0] = U; UBV[1] = B; UBV[2] = V;
		return UBV;
	}
	
	public void bidiagonalize(){ // bidiaonalize A = UBV^t where U and V are the operations to zero the columns and rows (respectively)
		Matrix vk;
		Matrix wk;
		Matrix Hk;
		Matrix Kk;
		double sign;
		double alphak;
		
		for (int k=0; k<M; ++k){

			//zero out entries below diagonal of k,kth entry
			 vk = new Matrix(1,M);
			for(int i=k;i<M;++i){
				vk.data[0][i] = B.data[i][k];
			}
			 sign = (int) Math.signum(A.data[k][k]);
			 alphak = -sign*vk.l2norm();
			vk.data[0][k] += alphak;
			 wk = vk.scale(1/vk.l2norm());
			 Hk = Matrix.identity(M).minus(wk.outer(wk.T()).scale(2));
			//compute updated B and U
			B = Hk.times(B);
			U = U.times(Hk);
			
			//zero out entries right of diagonal of k,k+1th entry
			vk = new Matrix(1,N);
			for(int i=k+1;i<N;++i){
				vk.data[0][i] = B.data[k][i];
			}
			 sign = (int) Math.signum(A.data[k][k+1]);
			 alphak = -sign*vk.l2norm();
			vk.data[0][k+1] += alphak;
			wk = vk.scale(1/vk.l2norm());
			Kk = Matrix.identity(N).minus(wk.outer(wk.T()).scale(2));
			//compute updated B and V
			B = B.times(Kk);
			V = V.times(Kk);
		}
	}

	
    public static void main(String[] args) {
    	
      double[][] a = { {4, 3, 0, 2}, {2, 1, 2, 1}, {4, 4, 0, 3} };
      Matrix A = new Matrix(a);
    	
    	
    	SVD B = new SVD(A);
    	Matrix [] UBV = B.factor();
//    	B.U.show();
//    	System.out.println("\n");
//    	B.B.show();
//    	System.out.println("\n");
//    	B.V.show();
    	//B.U.times(B.B).times(B.V.T()).show();
    	
    }
	
}