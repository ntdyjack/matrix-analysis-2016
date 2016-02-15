package matrixanalysis;

final public class QR {
	
	public Matrix[] factor(Matrix V){ // returns gramschmidt / QR decomposition
		int[] shape = V.shape();
		int M = shape[0]; int N = shape[1];
		Matrix Q = new Matrix(M, N); // Declare Orthonormal Matrix Q (ie the gramschmidt)
		Matrix R = new Matrix(N, N); // Declare UT matrix of operations R
		Matrix q = new Matrix(M, 1);  // A column of V to make orthonormal
		Matrix z = new Matrix(M, 1);  //
		
		double scaleFactor;
		double entry;
		
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
					entry = Q.getCol(j).inner(V.getCol(i));
					R.setVal(j, i, entry);  // set the off diagonal entries of the ith column of R
				}
				 continue; 
			}
			
			else Q.setCol(i, q.scale( 1 / scaleFactor));
			
			entry = Q.getCol(i).inner(V.getCol(i));
			R.setVal(i, k, entry); // set the diagonal entry of R
			
			for(int j=0; j<k; ++j) {
				entry = Q.getCol(j).inner(V.getCol(i));
				R.setVal(j, i, entry); // set the off diagonal entries of the ith column of R
			}
		}
		
		Matrix[] a = new Matrix[2];
		a[0] = Q; a[1] = R;
		return (a);
	}
	
	public static void main(String[] args) {
		Matrix A = Matrix.random(4, 3);
		Matrix[] QQRR = new Matrix[2];
		QR test = new QR();
		
		A.show();
		
		QQRR = test.factor(A);
//		QQRR[0].show(); QQRR[1].show();
		System.out.println("\n\n");
		
		Matrix C = QQRR[0].times(QQRR[1]);
		C.show();
	}
}
