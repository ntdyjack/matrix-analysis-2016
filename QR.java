package matrixanalysis;

final public class QR {
    
	private final Matrix V;
	private final Matrix Q;
	private final Matrix R;
	private final Matrix Eigen;
	
    public QR(Matrix A) {
    	this.V = A;
    	int N = A.N(); int M = A.M();
    	this.Q = new Matrix(M,N);
    	this.R = new Matrix(N,N);
    	this.Eigen = new Matrix(M,N);
    }
    
    
    public Matrix eigen(int nIter){
    	int i = 0;
		Matrix Eigen = this.Eigen;
		Eigen = this.R().times(this.Q());
    	while(i++ < nIter){
    		QR L = new QR(Eigen);
    		L.factor();
    		Eigen = L.R().times(L.Q());
    		//System.out.println(Eigen.max());
    	}
    	System.out.println("Maximum Iterations Reached \n");
    	return Eigen;
    }

	public void factor(){ // returns gramschmidt / QR decomposition
		//Matrix V = this;
		int M = V.M();
		int N = V.N();
		Matrix Q = this.Q;
		Matrix R = this.R;
		Matrix q = new Matrix(M, 1);  // A column of V to make orthonormal
		Matrix z = new Matrix(M, 1);  //
		
		double scaleFactor;
		
		int i; int k;
		for(i=0, k=0; i<N; i++, k++){ // loop over number of matrix
			q = V.getCol(i);  // grab the ith column of V
	
			for(int j=0; j<k; ++j) {
				z = Q.getCol(j);
				scaleFactor = V.getCol(i).inner(z); 
				z = z.scale(scaleFactor / z.inner(z) );
    			q = q.minus(z);
			}
			
			scaleFactor = Math.sqrt(q.inner(q));
			if (scaleFactor < 0.001){
				for(int j=0; j<k; ++j) {
					R.data[j][i] = Q.getCol(j).inner(V.getCol(i)); // set the off diagonal entries of the ith column of R
				}
				 continue; 
			}
			
			else Q.setCol(i, q.scale( 1 / scaleFactor));
			
			R.data[i][k] = Q.getCol(i).inner(V.getCol(i)); // set the diagonal entry of R
			
			for(int j=0; j<k; ++j) {
				R.data[j][i] = Q.getCol(j).inner(V.getCol(i)); // set the off diagonal entries of the ith column of R
			}
			
		}
	}
    
	public Matrix Q() { return this.Q; }
	
	public Matrix R() { return this.R; }
	
	
	public static void main(String[] args) {
		Matrix A = Matrix.random(3, 3);
		A.show();
		System.out.println("\n");
		QR B = new QR(A);
		B.factor();
//		Matrix Q = B.Q();
//		Matrix R = B.R();
//		
//		Q.show();
//		System.out.println("\n\n");
//		R.show();
//		System.out.println("\n\n");
//		Matrix z = Q.times(R);
//		z.show();
//		
//		System.out.println("\n\n");
		Matrix UT = B.eigen(1000); // this integer is the number of iterations
		UT.show();
		
	}
}
