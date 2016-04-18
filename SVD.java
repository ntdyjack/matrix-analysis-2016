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
		double c; double s; double r; double t; double t1;
		if(Math.abs(f) < 0.00001){
			c = 0; s = 1; r = g;
		}
		else if(Math.abs(f)>Math.abs(g)){
			t = g/f;
			t1 = Math.sqrt(1 + Math.pow(t,2));
			c = 1/t1;
			s = t*c;
			r = f*t1;
		}
		else{
			t = f/g;
			t1 = Math.sqrt(1 + Math.pow(t,2));
			s = 1/t1;
			c = t*s;
			r = g*t1;
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
    	
    	
    	for(int j = 0; j<500; ++j){
    	for( int k = 0; k < M-1; ++k){
    		
    		rotation = rot(L.B.data[k][k],L.B.data[k][k+1]);
    		c = rotation[0]; s = rotation [1]; r = rotation[2];
    		
    		// construct matrix Q and multiply on the right by Q’
    		// this annihilates both B(k-1,k+1) and B(k,k+1)
    		// but makes B(k+1,k) non-zero
    		Q = Matrix.identity(N);
    		
    		
    		Q.data[k][k] = c;
    		Q.data[k+1][k+1] = c;
    		Q.data[k][k+1] = s;
    		Q.data[k+1][k] = -s;
    		
    		L.B = L.B.times(Q.T());
    		L.V = Q.T().times(L.V);
    		
    		// construct matrix Q and multiply on the left by Q
    		// This annihilates B(k+1,k) but makes B(k,k+1) and
    		// B(k,k+2) non-zero
    		rotation = rot(L.B.data[k][k],L.B.data[k+1][k]);
    		c = rotation[0]; s = rotation [1]; r = rotation[2];
    		
    		Q = Matrix.identity(N);
    		
    		Q.data[k][k] = c;
    		Q.data[k+1][k+1] = c;
    		Q.data[k][k+1] = s;
    		Q.data[k+1][k] = -s;
    		
    		L.U = L.U.times(Q);
    		L.B = Q.times(L.B);
    	}
    	}
		UBV[0] = L.U; UBV[1] = L.B; UBV[2] = L.V;
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
			
			if(vk.iszero()==true){break;} // check if it's all zero and skip if yes
			
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
    	Matrix A = Matrix.random(3, 3);
    	A.show();
    	System.out.println("\n");
    	
    	SVD B = new SVD(A);
    	Matrix [] vals = B.factor();
    	
    	vals[0].show();
    	System.out.println("\n");
    	vals[1].show();
    	System.out.println("\n");
    	vals[2].show();
    	System.out.println("\n");
    	//Matrix S = vals[0].times(vals[1]);
    	//S.times(vals[2]).show();
    	//(vals[0].times(vals[1])).times(vals[2]).show();
    	
    }
	
}