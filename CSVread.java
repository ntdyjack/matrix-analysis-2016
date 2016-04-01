package matrixanalysis;

//reads CSV files into a Matrix object

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;

public class CSVread  {
	//private final Matrix Mat;
	//private final Matrix Ver;
	//private final Matrix Vir;
	
	public CSVread() {
    	//this.Mat = new Matrix(1,1);
    	//this.Ver = new Matrix(4,50);
    	//this.Vir = new Matrix(4,50);
	}
	
	private Double returnval(String val){
		double ans;
		
		if(val.equals("NA") | val.equals("na")|val.equals("Na")){ans = Double.NaN;}
		else{ans = Double.valueOf(val);}
		
		
		return ans;
	}
	
	
	public Matrix go(String Path, Boolean colnames, Boolean rownames) throws FileNotFoundException {

		
        Scanner scanner = new Scanner(new File(Path));
        scanner.useDelimiter(",");
        
        int ncol=0; int nrow = 0;
        
        //searches by comma separations until the first carriage return
        // this gives us the number of columns
        boolean m = true;
        while(m==true){
        	ncol++;
        	String [] x = scanner.next().split("\n");
        	if(x.length==2){m = false;}
        }   
        
        
        //reset the scanner
        // read the number of rows
         scanner = new Scanner(new File(Path));
        scanner.useDelimiter("\n");
        while(scanner.hasNext()){
        	nrow++;
        	scanner.next();
        }
        
        if(colnames == true){nrow--;}
        if(rownames == true){ncol--;}
        
        
        //declare the matrix
        Matrix  Mat = new Matrix(nrow,ncol);
        int currentrow = 0; int currentcol = 0;
        scanner = new Scanner(new File(Path));
        //System.out.println(ncol);
        if(colnames == true){
        	scanner.useDelimiter("\n");
        	scanner.next();
        }
        scanner.useDelimiter(",|\n");
        
        if(rownames == true){
        	while(currentrow<nrow){
        		if(currentcol==0){scanner.next();}
        		String x = scanner.next();
        		Mat.data[currentrow][currentcol] = returnval(x);
        		currentcol++;
        		if( currentcol == ncol ){currentrow++; currentcol=0;}
        	}
        	
        } else {
        	while(currentrow<nrow){ 
        		String x = scanner.next();
        		Mat.data[currentrow][currentcol] = returnval(x);
        		currentcol++;
        		if( currentcol == ncol ){currentrow++;}
        	}
        	
        }
        
        
    	scanner.close();
        return Mat;   
	}
        
	public static void main(String[] args) throws FileNotFoundException {
		CSVread read = new CSVread();
		//pass the file path, are there colnames? are there rownames?
		Matrix x = read.go("matrixanalysis/DSGaggregate.csv",true,true);
		System.out.println(x.data[14150][2]);
	}
	
	}
		