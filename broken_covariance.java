package matrixanalysis;
import java.lang.Math;
import java.lang.Double;

public class broken_covariance {
	
	private final Matrix A;
	private final int N;
	private final int M;
	public Matrix F;
	
	public broken_covariance(Matrix A) {
		this.A = A;
		this.N = A.N();
		this.M = A.M();
	}
	
	private int[] find(Matrix vector) {
		int M = this.M;
		int[] missing = new int[M];
		for(int i=0; i<M; ++i) {
			if (Double.isNaN(vector.data[i][0]))
				missing[i] = 1;
		}
		return missing;
	}
	
	private int[] add_missing(int[] missing1, int[] missing2) {
		int M = this.M;
		int[] missing = new int[M];
		for(int i=0; i<M; ++i) {
			missing[i] = missing1[i] + missing2[i];
		}
		return missing;
	}
	
	public void fix() {
		
	}
	
}
