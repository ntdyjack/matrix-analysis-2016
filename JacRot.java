package matrixanalysis;
import java.lang.Math;

public class JacRot {
	
	public Matrix A;
	public Matrix D;
	public final int N;
	
	public JacRot(Matrix A) {
		this.A = A;
		this.N = A.N();
		this.D = Matrix.identity(N);
	}
	
	
	public void factor(int maxiter, double tolerance) {
		Matrix A = this.A;
		Matrix D = this.D;
		int N = this.N; int i=0;
		int[] max = A.maxLowDiag();
		double theta;
		double aii; double aij; double ajj;
		
		aij = A.data[max[0]][max[1]];
		while (Math.abs(aij) > tolerance && i++ < maxiter) {
			Matrix R = Matrix.identity(N);
			aii = A.data[max[0]][max[0]]; aij = A.data[max[0]][max[1]]; ajj = A.data[max[1]][max[1]];
			if ( aii == ajj ) {
				theta = Math.PI / 2;
			}
			else {
				theta = (1/2.0) * Math.atan(2.0*aij / (aii-ajj));
			}
			R.data[max[0]][max[0]] = Math.cos(theta);
			R.data[max[1]][max[1]] = Math.cos(theta);
			R.data[max[0]][max[1]] = Math.sin(theta) * (-1);
			R.data[max[1]][max[0]] = Math.sin(theta);
			
			A = R.T().times(A.times(R));
			D = R.T().times(D);
			max = A.maxLowDiag();
		}
		this.A = A;
		this.D = D;
	}
	
	public static void main(String args[]) {
		Matrix A = Matrix.random(5, 5);
		Matrix B = A.times(A.T());
		
		JacRot JR = new JacRot(B);
		B.show();
		System.out.println("test");
		JR.factor(1000, 0.00001);
		
		JR.A.show();
		Matrix R = JR.D;
		R.show();
		System.out.println();
		R.T().times(R).show();
	}
}
