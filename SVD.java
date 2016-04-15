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
	}
	
	public double [] rot(double f, double g){
		double c; double s; double r;
		if(f == 0){
			c = 0; s = 1; r = g;
		}
		else if(Math.abs(f)>Math.abs(g)){
			double t = Math.sqrt(1 + Math.pow(g/f,2));
			c = 1/t;
			s = t*c;
			r = f*t;
		}
		else{
			double t = Math.sqrt(1 + Math.pow(f/g,2));
			s = 1/t;
			c = t*s;
			r = g*t;
		}
		double [] vect = new double [3];
		vect [0] = c; vect [1] = s; vect [2] = r;
		return vect;
	}
	
	public Matrix [] factor(){
		Matrix [] UBV = new Matrix [3];
		double c;
		double s;
		double r;
		double rotation [];
		Matrix Q;
    	SVD L = new SVD(A);
    	L.bidiagonalize();
    	L.B.show();
    	// Demmel & Kahan zero-shift QR downward sweep
    	for( int k = 0; k < M-1; ++k){
    		rotation = rot(L.B.data[k][k],L.B.data[k][k+1]);
    		c = rotation[0]; s = rotation [1]; r = rotation[2];
    		
    		// construct matrix Q and multiply on the right by Q’
    		// this annihilates both B(k-1,k+1) and B(k,k+1)
    		// but makes B(k+1,k) non-zero
    		Q = Matrix.identity(N);
    		
    		
    		Q.data[k][k] = c;
    		Q.data[k+1][k+1] = c;
    		Q.data[k][k+1] = -s;
    		Q.data[k+1][k] = s;
    		
    		//L.B.show();
    		//System.out.println("\n");
    		//		(k:k+1,k:k+1)=[c s;-s c];
    		//L.B = L.B.times(Q.T());
    		
    		
    		// construct matrix Q and multiply on the left by Q
    		// This annihilates B(k+1,k) but makes B(k,k+1) and
    		// B(k,k+2) non-zero
    		
    		//L.B.show();
    		//System.out.println("\n");
    		
    		rotation = rot(L.B.data[k][k],L.B.data[k][k+1]);
    		c = rotation[0]; s = rotation [1]; r = rotation[2];
    		
    		Q = Matrix.identity(M);
    		
    		Q.data[k][k] = c;
    		Q.data[k+1][k+1] = c;
    		Q.data[k][k+1] = -s;
    		Q.data[k+1][k] = s;
    		
    		//System.out.println("\n");
    		//Q.show();
    		
    		L.B = Q.times(L.B);
    	}
		
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
		
		for (int k=0; k<M-1; ++k){

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
    	
      //double[][] a = { {4, 3, 0, 2}, {2, 1, 2, 1}, {4, 4, 0, 3} };
      //Matrix A = new Matrix(a);
    	Matrix A = Matrix.random(4, 4);
    	
    	
    	SVD B = new SVD(A);
    	B.factor();
    	//B.bidiagonalize();
    	//B.B.show();
    	//Matrix [] UBV = B.factor();
    	//UBV[1].show();
//    	System.out.println("\n");
//    	B.B.show();
//    	System.out.println("\n");
//    	B.V.show();
    	//B.U.times(B.B).times(B.V.T()).show();
    	
    }
	
}