package matrixanalysis;


import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;


public class CSVread  {
	private final Matrix Set;
	private final Matrix Ver;
	private final Matrix Vir;
	
	public CSVread() {
    	this.Set = new Matrix(4,50);
    	this.Ver = new Matrix(4,50);
    	this.Vir = new Matrix(4,50);
	}
	
	public Matrix [] go() throws FileNotFoundException {

		Matrix [] Iris = new Matrix[3];
		
        Scanner scanner = new Scanner(new File("/home/nikwoj/workspace/matrix-analysis-2016/matrixanalysis/iris.csv"));
        scanner.useDelimiter(",|\n");
        int counter = 0;
        int colnum = 0;
        String name = "a";
        
        
        while(scanner.hasNext()){
        	counter ++;
        	
            String x = scanner.next();
        	if ( (counter % 5)==0){
        		 name = x;
        		 colnum++;
        		 continue;
        	}
        	
        	if (counter >5 & ( name.equals("Iris-setosa") | name.equals("name")) ){
        		Set.data[(counter % 5)-1][colnum-2] = Double.parseDouble(x);
        	} else if (name.equals("Iris-versicolor")){
        		Ver.data[(counter % 5)-1][colnum-52] = Double.parseDouble(x);
        	} else if (name.equals("Iris-virginica")) {
        		Vir.data[(counter % 5)-1][colnum-102] = Double.parseDouble(x);
        	};
        }
    	scanner.close();
    	Iris[0] = Set; Iris[1] = Ver; Iris[2] = Vir;
        return Iris;   
	}
        

		
	}