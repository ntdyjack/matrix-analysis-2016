package matrixanalysis;
import java.lang.Math;


/******************************************************************************
 *  Compilation:  javac Matrix.java
 *  Execution:    java Matrix
 *
 *  A bare-bones immutable data type for M-by-N matrices.
 *  
 *  Added functions (Nathan Dyjack + Nik Wojtalewicz)
 *  
 *  scalar multiplication of matrices .scalar()
 *  get specified column(s) and row(s) from a matrix .getCol/Row()
 *  set values of specified rows/columns in a matrix .setCol/Row()...
 *  inner product of two column vectors .inner()
 *  A=QR factorization (the matrix Q is the gram-schmidt orthogonal matrix) .factor()
 *  A=LU factorization (lower/upper triangular matrices)
 *  
 *  To do : 
 *  Finish setting up LU class
 *   - Make backward solver
 *   - Make forward solver
 *   - Write det function
 *  
 *  Rename MaxLowDiag and maxLowDiag functions
 *
 ******************************************************************************/

import java.lang.Math;

final public class Matrix {
    private final int M;             // number of rows
    private final int N;             // number of columns
    public final double[][] data;   // M-by-N array

    
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
    
    
    public int M() { return this.M; }
    public int N() { return this.N; }
    
    public void setVal(int i, int j, double val) { this.data[i][j] = val; }
    
    public double val(int i, int j) { return this.data[i][j]; }
    
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

    //return the l2 norm of a row(or column) vector
    public double l2norm(){
    	double norm = 0;
    	if(M!=1 && N!=1) throw new RuntimeException("not a row or column vector");
    	
    	if(N>M){ //we have a row vector
    		for(int i =0; i<N; i++){norm += Math.pow(this.data[0][i],2);}
    	}
    	else{ //we have a column vector
    		for(int i =0; i<M; ++i){norm += Math.pow(this.data[i][0],2);}
    	}
    	return Math.sqrt(norm);
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

    public double max(){
    	Matrix A = this;
    	double max = 0;
        for (int i = 0; i < M; i++)
            for (int j = 0; j < N; j++)
            	if (Math.abs(A.data[i][j]) > max)
            		max = A.data[i][j];
        return max;
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
    
    
    public void setRow(int rowNum, Matrix Row) {
    	Matrix A = this;
    	if (rowNum > A.N || rowNum < 0) throw new RuntimeException("Column number "
    			+ "must be less than number of columns and greater than zero");
    	if (Row.M > 1) throw new RuntimeException("Column to add must be column vector");
    	
    	for (int i = 0; i < A.N; ++i) 
    		A.data[rowNum][i] = Row.data[0][i];
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
    
    public Matrix outer(Matrix W){
    	Matrix V = this;
    	if (V.M != W.N) throw new RuntimeException("Dimensions must be same");
    	Matrix OP = new Matrix(V.N,W.M);
    	for(int i = 0; i<V.N;++i){
    		for(int j = 0; j<W.M;++j){
    			OP.data[i][j] = V.data[0][i]*W.data[j][0];
    		}
    	}
    	
    	return OP;
    }
    
    public Matrix[] LU(){
    	Matrix V = this;
    	Matrix L = new Matrix(V.M,V.N);
    	Matrix U = new Matrix(V.M,V.N);
    	
    	U.setRow(0, V.getRow(0)); // set the first row of U
    	
    	L.setCol(0, V.getCol(0).scale(1 / U.val(0, 0)));
    	
//    	for(int i = 0; i<V.M; i++){ // set first column of L
//    		L.data[i][0] = V.data[i][0]/U.data[0][0];
//    	}
    	

    	int i; int j; int k; double sum = 0;
    	for(i=1;i<V.M;i++){
    		
    		if(i>1){ //set the off diagonals of L
    			for(j=1;j<i;j++){
    			sum=0;
    			k=0;
    			while(k<j){
    				sum += U.data[k][j]*L.data[i][k];
    				k++;
    			}
    			L.data[i][j]= (1/U.data[j][j])*(V.data[i][j] - sum);}
    		}
    		
    		if(i>0 && i<V.M-1){ // set the off diagonals of U
    			for(j=i+1;j<V.M;j++){
	    			sum = 0;
	    			k=0;
	    			while(k<i){
	    				sum += U.data[k][j]*L.data[i][k];
	    				k++;
	    			}
		    		U.data[i][j] = V.data[i][j] - sum;}
    		}

	    		sum=0; //j=i; // set the diagonal entry of U and L
	    		for( k=0; k<i; k++){ 
	    			sum += U.data[k][i]*L.data[i][k]; 
	    		}
	    		U.data[i][i]= V.data[i][i] - sum; //uij = aij - sum(k to i-1) ukj*lik
	    		L.data[i][i]= 1; //(1/U.data[j][j])*(V.data[i][j] - sum); //lij = 1/uij * (aij - sum(k to j-1) ukj*lik
	    		
    		}
		Matrix[] a = new Matrix[2];
		a[0] = L; a[1] = U;
		return (a);
    }
    
    public double MaxLowDiag() {
    	int N = this.N;
		Matrix V = this;
		double maxElem = 0;
		for (int i=0; i < N; ++i) {
			for (int j=i+1; j < N; ++j) {
				if (Math.abs(V.data[j][i]) > maxElem) {
					maxElem = Math.abs(V.data[j][i]);
				}
			}
		}
		return maxElem;
    }
    
	public int[] maxLowDiag() {
		int N = this.N;
		Matrix V = this;
		int[] max = new int[2];
		double maxElem = 0;
		for (int i=0; i < N; ++i) {
			for (int j=i+1; j < N; ++j) {
				if (Math.abs(V.data[j][i]) > maxElem) {
					maxElem = Math.abs(V.data[j][i]);
					max[0] = j; max[1] = i;
				}
			}
		}
		return max;
	}
	
    // test client
    public static void main(String[] args) {

//        double[][] a = { {1, 2, 3, 4}, {4, 5, 6, 4}, {7, 8, 9, 4} };
//        Matrix A = new Matrix(a);

//    	
//        System.out.println("A =");
//        A.show();
//        System.out.println();
//        
//        Matrix[] Q = A.LU(); // call QR factorization
//        Matrix L = Q[0]; Matrix U= Q[1];
//        System.out.println("L =");
//        L.show();
//        System.out.println();
//        System.out.println("U =");
//        U.show();
//        System.out.println();
//        Matrix B = L.times(U);
//        System.out.println("LU =");
//        B.show();
        
    	//Matrix A = Matrix.random(1,4);
        //double[][] a = { {1, 2, 3}};
        //Matrix A = new Matrix(a);
    	//A.show();
        //A.show();
        //A.outer(A.T()).show();
    	//System.out.println(A.l2norm());
    	//A.T().show();
    	//System.out.println(A.T().l2norm());
        
//        System.out.format("%d%n", B.shape()[0]);
    }
}
