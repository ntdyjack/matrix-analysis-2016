package matrixanalysis;



/******************************************************************************
 *  Compilation:  javac Matrix.java
 *  Execution:    java Matrix
 *
 *  A bare-bones immutable data type for M-by-N matrices.
 *  
 *  Added functions (Nathan Dyjack + Nik Wojtalewicz)
 *  
 *  scalar multiplication of matrices .scalar()
 *  get specified column(s) and row(s) from a matrix .getCol()
 *  set values of specified rows/columns in a matrix .setCol()
 *  inner product of two column vectors .inner()
 *  A=QR factorization (the matrix Q is the gram-schmidt orthogonal matrix) .factor()
 *  
 *  
 *
 ******************************************************************************/

final public class Matrix {
    private final int M;             // number of rows
    private final int N;             // number of columns
    private final double[][] data;   // M-by-N array

    // create M-by-N matrix of 0's
    public Matrix(int M, int N) {
        this.M = M;
        this.N = N;
        data = new double[M][N];
    }

    // create matrix based on 2d array
    public Matrix(double[][] data) {
        M = data.length;
        N = data[0].length;
        this.data = new double[M][N];
        for (int i = 0; i < M; i++)
            for (int j = 0; j < N; j++)
                    this.data[i][j] = data[i][j];
    }

    // copy constructor
    private Matrix(Matrix A) { this(A.data); }

    // create and return a random M-by-N matrix with values between 0 and 1
    public static Matrix random(int M, int N) {
        Matrix A = new Matrix(M, N);
        for (int i = 0; i < M; i++)
            for (int j = 0; j < N; j++)
                A.data[i][j] = Math.random();
        return A;
    }

    // create and return the N-by-N identity matrix
    public static Matrix identity(int N) {
        Matrix I = new Matrix(N, N);
        for (int i = 0; i < N; i++)
            I.data[i][i] = 1;
        return I;
    }

    // swap rows i and j
    private void swap(int i, int j) {
        double[] temp = data[i];
        data[i] = data[j];
        data[j] = temp;
    }

    // create and return the transpose of the invoking matrix
    public Matrix T() {
        Matrix A = new Matrix(N, M);
        for (int i = 0; i < M; i++)
            for (int j = 0; j < N; j++)
                A.data[j][i] = this.data[i][j];
        return A;
    }

    // return C = A + B
    public Matrix plus(Matrix B) {
        Matrix A = this;
        if (B.M != A.M || B.N != A.N) throw new RuntimeException("Illegal matrix dimensions.");
        Matrix C = new Matrix(M, N);
        for (int i = 0; i < M; i++)
            for (int j = 0; j < N; j++)
                C.data[i][j] = A.data[i][j] + B.data[i][j];
        return C;
    }


    // return C = A - B
    public Matrix minus(Matrix B) {
        Matrix A = this;
        if (B.M != A.M || B.N != A.N) throw new RuntimeException("Illegal matrix dimensions.");
        Matrix C = new Matrix(M, N);
        for (int i = 0; i < M; i++)
            for (int j = 0; j < N; j++)
                C.data[i][j] = A.data[i][j] - B.data[i][j];
        return C;
    }

    // does A = B exactly?
    public boolean eq(Matrix B) {
        Matrix A = this;
        if (B.M != A.M || B.N != A.N) throw new RuntimeException("Illegal matrix dimensions.");
        for (int i = 0; i < M; i++)
            for (int j = 0; j < N; j++)
                if (A.data[i][j] != B.data[i][j]) return false;
        return true;
    }

    // return C = A * B
    public Matrix times(Matrix B) {
        Matrix A = this;
        if (A.N != B.M) throw new RuntimeException("Illegal matrix dimensions.");
        Matrix C = new Matrix(A.M, B.N);
        for (int i = 0; i < C.M; i++)
            for (int j = 0; j < C.N; j++)
                for (int k = 0; k < A.N; k++)
                    C.data[i][j] += (A.data[i][k] * B.data[k][j]);
        return C;
    }

    // return x = A^-1 b, assuming A is square and has full rank
    public Matrix solve(Matrix rhs) {
        if (M != N || rhs.M != N || rhs.N != 1)
            throw new RuntimeException("Illegal matrix dimensions.");

        // create copies of the data
        Matrix A = new Matrix(this);
        Matrix b = new Matrix(rhs);

        // Gaussian elimination with partial pivoting
        for (int i = 0; i < N; i++) {

            // find pivot row and swap
            int max = i;
            for (int j = i + 1; j < N; j++)
                if (Math.abs(A.data[j][i]) > Math.abs(A.data[max][i]))
                    max = j;
            A.swap(i, max);
            b.swap(i, max);

            // singular
            if (A.data[i][i] == 0.0) throw new RuntimeException("Matrix is singular.");

            // pivot within b
            for (int j = i + 1; j < N; j++)
                b.data[j][0] -= b.data[i][0] * A.data[j][i] / A.data[i][i];

            // pivot within A
            for (int j = i + 1; j < N; j++) {
                double m = A.data[j][i] / A.data[i][i];
                for (int k = i+1; k < N; k++) {
                    A.data[j][k] -= A.data[i][k] * m;
                }
                A.data[j][i] = 0.0;
            }
        }

        // back substitution
        Matrix x = new Matrix(N, 1);
        for (int j = N - 1; j >= 0; j--) {
            double t = 0.0;
            for (int k = j + 1; k < N; k++)
                t += A.data[j][k] * x.data[k][0];
            x.data[j][0] = (b.data[j][0] - t) / A.data[j][j];
        }
        return x;
   
    }

    // print matrix to standard output
    public void show() {
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) 
                System.out.printf("%9.4f ", data[i][j]);
            System.out.println();
        }
    }
    
    public void setCol(int colNum, Matrix Col) {
    	Matrix A = this;
    	if (colNum > A.N || colNum < 0) throw new RuntimeException("Column number "
    			+ "must be less than number of columns and greater than zero");
    	if (Col.N > 1) throw new RuntimeException("Column to add must be column vector");
    	
    	for (int i = 0; i < A.M; ++i) 
    		A.data[i][colNum] = Col.data[i][0];
    }
    
    public Matrix getCols( int[] columns ) {
    	Matrix col = new Matrix(M, columns.length);
    	for (int i=0; i < columns.length; ++i){
    		for (int j = 0; j < M; j++){
    			col.data[j][i] = data[j][columns[i]];
    		}
    	}
    	return col;	
    }

    public Matrix getCol( int column_num ) {
    	if(column_num > N || column_num < 0 ) throw new RuntimeException("Column Number "
    			+ "must be less than number of columns in A and bigger than 0");
		Matrix col = new Matrix(M, 1);

		for (int j = 0; j < M; j++){
			col.data[j][0] = data[j][column_num];
		}		
	return col;	
    }    
    
    public Matrix getRows( int[] rows )  {
		Matrix row = new Matrix(rows.length, N);
		for (int i=0; i < rows.length; ++i){
			for (int j = 0; j < N; j++){
				row.data[i][j] = data[rows[i]][j];
			}
		}
		return row;
    }    
    
    public Matrix getRow( int row_num )  {
		Matrix row = new Matrix(1, N);
		
		for (int j = 0; j < N; j++){
			row.data[0][j] = data[row_num][j];
		}
		return row;
    }    
    
    public Matrix scale(double a) {
    	Matrix V = this; Matrix B = new Matrix(V.M, V.N);
    	for(int i=0; i<V.M; ++i) {
    		for(int j=0; j<V.N; ++j) {
    			B.data[i][j] = a * V.data[i][j];
    		}
    	}
    	return B;
    }
    
    public double inner(Matrix W) {
    	Matrix V = this;
    	if (V.N != 1 && W.N != 1) throw new RuntimeException("Illegal Matrix Dimensions");
    	if (V.M != W.M) throw new RuntimeException("Dimensions must be same");
    	
    	double sum = 0;
    	for (int i=0; i < V.M; ++i) 
    		sum += V.data[i][0] * W.data[i][0];
    	return sum;
    }
    
    // Take base matrix and apply Gram Schimdt to its columns
    public Matrix GramSchmidt(){
    	Matrix V = this;
    	Matrix R = new Matrix(V.M,V.N); // Declare Orthonormal Matrix R
    	Matrix q = new Matrix(V.M, 1);  // A column of V to make orthonormal
    	Matrix z = new Matrix(V.M, 1);  // A column of R
    	double scaleFactor;
    	
    	int i; int k;
    	for(i=0, k=0; i<V.N; i++, k++){ // loop over number of matrix
    		q.setCol(0, V.getCol(i));

    		for(int j=0; j<k; ++j) {
    			scaleFactor = q.inner(R.getCol(j));
    			z.setCol(0, R.getCol(j));
    			q = q.minus(z.scale(scaleFactor)); // Ugh
    			
    		}
    		
    		scaleFactor = Math.sqrt(q.inner(q));
    		if (scaleFactor < 0.01) 
    			 continue;
    		else R.setCol(i, q.scale( 1 / scaleFactor));
    	}
    	return R;
    }

    
	public Matrix[] factor(){ // returns gramschmidt / QR decomposition
		Matrix V = this;
		Matrix Q = new Matrix(V.M,V.N); // Declare Orthonormal Matrix Q (ie the gramschmidt)
		Matrix R = new Matrix(V.M,V.N); // Declare UT matrix of operations R

		Matrix q = new Matrix(V.M, 1);  // A column of V to make orthonormal
		Matrix z = new Matrix(V.M, 1);  //
		
		double scaleFactor;
		
		int i; int k;
		for(i=0, k=0; i<V.N; i++, k++){ // loop over number of matrix
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
		
		Matrix[] a = new Matrix[2];
		a[0] = Q; a[1] = R;
		return (a);
	}
    
    
    // test client
    public static void main(String[] args) {

<<<<<<< HEAD
        //double[][] a = { {1, 2, 3, 4}, {4, 5, 6, 4}, {7, 8, 9, 4} };
        //Matrix A = new Matrix(a);
       Matrix A = Matrix.random(3,3);
=======
        double[][] a = { {1, 2, 3}, {4, 5, 6}, {7, 8, 9} };
        //Matrix A = new Matrix(a);
        Matrix A = Matrix.random(3,3);
>>>>>>> origin/master
        System.out.println("A =");
        A.show();
        System.out.println();
        
        
        Matrix[] L = A.factor(); // call QR factorization
        Matrix Q = L[0]; Matrix R = L[1];
        System.out.println("Q =");
		Q.show();
		System.out.println("R =");
		R.show();
        Matrix B = Q.times(R);
		System.out.println("A = QR");
        B.show();
    }
}
