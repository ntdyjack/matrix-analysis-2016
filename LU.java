package matrixanalysis;

public class LU {
	private final Matrix V;
	private final Matrix L;
	private final Matrix U;
	
    public LU(Matrix A) {
    	this.V = A;
    	int N = V.N(); int M = V.M();
    	this.L = new Matrix(M,N);
    	this.U = new Matrix(M,N);
    }
    
    public void factor(){
    	Matrix V = this.V;
    	Matrix L = this.L;
    	Matrix U = this.U;
    	int N = V.N(); int M = V.M();
    	
    	U.setRow(0, V.getRow(0)); // set the first row of U
    	
    	L.setCol(0, V.getCol(0).scale(1 / U.val(0, 0)));

    	int i; int j; int k; double sum = 0;
    	for(i=1; i<M; i++){
    		
    		if(i>1){ //set the off diagonals of L
    			for(j=1;j<i;j++){
    			sum=0;
    			k=0;
    			while( k<j ){
    				sum += U.data[k][j]*L.data[i][k]; // This is just an inner product
    				k++;
    			}
    			L.data[i][j]= (1/U.data[j][j])*(V.data[i][j] - sum);}
    		}
    		
    		if(i>0 && i<M-1){ // set the off diagonals of U
    			for(j=i+1;j<M;j++){
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
    }
    
	public Matrix L() { return this.L; }
	
	public Matrix U() { return this.U; }
	
	public Matrix backsub(Matrix vector){
	    Matrix U = this.U;
	    int N = U.N(); int M = U.M();
    	Matrix x = new Matrix(N, 1);
	    for (int j = N - 1; j >= 0; j--) {
	        double t = 0.0;
	        for (int k = j + 1; k < N; k++){
	            t += U.data[j][k] * x.data[k][0];
	            }
	        x.data[j][0] = (vector.data[j][0] - t) / U.data[j][j];
	    }
	    return x;
	}
	
	public Matrix forwardsub(Matrix vector){
	    Matrix L = this.L;
	    int N = L.N(); int M = L.M();
    	Matrix x = new Matrix(N, 1);
	    for (int j = 0; j < M; j++) {
	        double t = 0.0;
	        for (int k = 0; k <= j; k++){
	            t += L.data[j][k] * x.data[k][0];
	            //System.out.println(k);
	            }
	        x.data[j][0] = (vector.data[j][0] - t);
	        //x.show();
	        //System.out.println("\n");
	    }
	    return x;
	}
    
	public Matrix LUsolve(Matrix vector){
    	Matrix L = this.L;
    	LU current = this;
    	
	    int N = L.N(); //int M = L.M();
    	Matrix x = new Matrix(N, 1);
    	
    	Matrix y = current.forwardsub(vector);
    	x = current.backsub(y);
    	
    	return(x);
    	
	}
    
	public static void main(String[] args) {
		//Matrix A = Matrix.random(3, 3);
		//Matrix b = Matrix.random(3, 1);
		//double[][] a = {{5,9,8}, {0,2, 8},  {0,0,1} };
      double[][] a = {{5,0,0}, {2, 8, 0},  {3, 9,8} };
      Matrix A = new Matrix(a);;
      double[][] b = { {5}, {7}, {8}};
      Matrix B = new Matrix(b);
		
		A.show();
		System.out.println("\n");
		B.show();
		System.out.println("\n");
		
		LU Z = new LU(A);
		Z.factor();
		
		System.out.println("\n");

		A.times(Z.LUsolve(B)).show();
	}
    
}


