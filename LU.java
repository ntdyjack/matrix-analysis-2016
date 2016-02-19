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
    	
//    	for(int i = 0; i<V.M; i++){ // set first column of L
//    		L.data[i][0] = V.data[i][0]/U.data[0][0];
//    	}
    	

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
    }
    
    public Matrix solve(Matrix vector) {// back substitution
	    Matrix V = this.V;
	    Matrix L = this.L;
	    Matrix U = this.U;
	    int N = V.N(); int M = V.M();
    	Matrix x = new Matrix(N, 1);
	    for (int j = N - 1; j >= 0; j--) {
	        double t = 0.0;
	        for (int k = j + 1; k < N; k++)
	            t += A.data[j][k] * x.data[k][0];
	        x.data[j][0] = (b.data[j][0] - t) / A.data[j][j];
	    }
	    return x;
    }
}
