import java.math.*;
import java.util.Arrays;
import java.util.Formatter;

public class Lab1 {
	public static double[] gaussSelectionByColumn(int matrixOrder,double [][] matrixA,double [] columnF){
		//������ ���
		double [] columnX = new double[matrixOrder];
        int maxIndex;
        double maxInColumn;
        double temp;
		for (int i=0; i<matrixOrder-1; i++) {
        	maxIndex=i;
        	maxInColumn=Math.abs(matrixA[i][i]);
        	for (int j=i; j<matrixOrder; j++){
        		if(maxInColumn<Math.abs(matrixA[j][i]))
        		{
        			maxInColumn=Math.abs(matrixA[j][i]);
        			maxIndex=j;
        		}
        		if(maxInColumn==0)
        			throw new ArithmeticException("������� ������� - ������� ������ ���� ����������");
        	}
        	for (int k=i; k<matrixOrder; k++){
        		temp=matrixA[i][k];
        		matrixA[i][k]=matrixA[maxIndex][k];
        		matrixA[maxIndex][k]=temp;
        	}
        	 temp=columnF[i];
        	 columnF[i]=columnF[maxIndex];
        	 columnF[maxIndex]=temp;
        	 for (int j=i+1; j<matrixOrder; j++)
        		 matrixA[i][j]= matrixA[i][j]/ matrixA[i][i];
        	 columnF[i]= columnF[i]/matrixA[i][i];
        	 matrixA[i][i]=1.0;
        	 for (int j=i+1; j<matrixOrder; j++)
        	 {
        		 for (int k=i+1; k<matrixOrder; k++)
        			 matrixA[j][k]= matrixA[j][k]+matrixA[i][k]*matrixA[j][i]*(-1.0);
        		 columnF[j]= columnF[j]+ columnF[i]*matrixA[j][i]*(-1.0);
    			 matrixA[j][i]=0;
        	 }
        }    
			columnF[matrixOrder-1]=columnF[matrixOrder-1]/matrixA[matrixOrder-1][matrixOrder-1];
            matrixA[matrixOrder-1][matrixOrder-1]=1;     
             
             //�������� ���
		 columnX[matrixOrder-1]=columnF[matrixOrder-1];
         for (int i=matrixOrder-2; i>-1; i--) {
         	columnX[i]=columnF[i];
         	for (int k=i; k<matrixOrder-1; k++)
         	columnX[i]=columnX[i]+matrixA[i][k+1]*columnX[k+1]*(-1.0);
         }
             return columnX;
	}
	
	
	public static void main(String[] args) {
		int matrixOrder = 10;
		double [][] matrixA = new double[matrixOrder] [matrixOrder];
		double [] columnX = new double[matrixOrder];
		double [] exact�olumnX = new double[matrixOrder];
		double [] columnF = new double[matrixOrder];
		
		 //������ ������� ����
		 for (int i=0; i<matrixOrder; i++)
			 exact�olumnX[i]=i+1;
		
		//��������� ������� A ���������� ������� �� ��������� �� -100 �� 100
        for (int i=0; i<matrixOrder; i++)
             for (int j=0; j<matrixOrder; j++)
            	 matrixA[i][j] = 100-(Math.random()*200);

        
        //������ ������� f
        for(int i = 0; i < matrixOrder; i++)
            for(int j = 0; j < matrixOrder; j++)
                columnF[i] += exact�olumnX[j] * matrixA[i][j];
        
        //�����
        try {
        	columnX=gaussSelectionByColumn(matrixOrder,Arrays.stream(matrixA).map(double[]::clone).toArray(double[][]::new),columnF);
        }catch(ArithmeticException e){
         	System.out.println(e.getMessage());
         	System.exit(0);
         }
        
        //�����������
        double normXX=Math.abs(exact�olumnX[0]-columnX[0]);
        double normX = Math.abs(exact�olumnX[0]);
        for (int i=1; i<matrixOrder; i++) {
        	if(Math.abs(exact�olumnX[i]-columnX[i])>normXX)
        		normXX=Math.abs(exact�olumnX[i]-columnX[i]);
        	if(Math.abs(exact�olumnX[i])>normX)
        		normX=Math.abs(exact�olumnX[i]);
        }
        double accuracy=normXX/normX*100;
        
        //���������� ��������
        double [][] inverseMatrixA = new double[matrixOrder][matrixOrder];
        double [][] unitMatrix = new double[matrixOrder][matrixOrder];
        for (int i=0; i<matrixOrder; i++)
        	for (int j=0; j<matrixOrder; j++) {
        		if(i==j)
        			unitMatrix[i][j]=1;
        		else
        			unitMatrix[i][j]=0;
        	}  
        double [] tempColumn = new double[matrixOrder];
        for (int i=0; i<matrixOrder; i++) {
        	tempColumn=unitMatrix[i];
        	 try {
        		 tempColumn=gaussSelectionByColumn(matrixOrder,Arrays.stream(matrixA).map(double[]::clone).toArray(double[][]::new),tempColumn);
             }catch(ArithmeticException e){
              	System.out.println(e.getMessage());
              	System.exit(0);
              }
        	 for (int j=0; j<matrixOrder; j++)
        		 inverseMatrixA[j][i]=tempColumn[j];
        }
        //������������ ������� � ��������
        double[][] mulMatr=new double[matrixOrder][matrixOrder];
        for (int i=0; i<matrixOrder; i++)
            for (int j=0; j<matrixOrder; j++)
              for (int k=0; k<matrixOrder; k++)
            	  mulMatr[i][j] += matrixA[i][k] * inverseMatrixA[k][j];
        
        //������ �����������
        System.out.println("������� � � ������� f");
        for (int i=0; i<matrixOrder; i++) {
        	Formatter f = new Formatter();
            for (int j=0; j<matrixOrder; j++) {
				f.format("%8.4f ", matrixA[i][j]);
            }
            f.format("| %8.4f%n", columnF[i]);
            System.out.print(f);                                            
       }
        System.out.println();   
        System.out.println("������ ������� �������                   ������ ������������� �������");
    	Formatter f = new Formatter();
        for (int j=0; j<matrixOrder; j++)
			f.format("%19.1f %40.16f%n", exact�olumnX[j], columnX[j]);
		 System.out.print(f);      
         System.out.print("������������� ����������� =");
         f = new Formatter();
         f.format("%19.16f%%%n", accuracy);
         System.out.print(f);  
         System.out.println("�������� ������� �");
         for (int i=0; i<matrixOrder; i++) {
         	 f = new Formatter();
             for (int j=0; j<matrixOrder; j++) {
 				f.format("%8.4f ", inverseMatrixA[i][j]);
             }
             System.out.println(f);                                            
        }
         System.out.println("������������ �*�^(-1)");
         for (int i=0; i<matrixOrder; i++) {
         	 f = new Formatter();
             for (int j=0; j<matrixOrder; j++) {
 				f.format("%12.4e ", mulMatr[i][j]);
             }
             System.out.println(f);   
        }
         
         
         //������ �������� ������� ������� �������� �������������� ���������
         System.out.println("������� �������     ������������� ����������� ");
         f = new Formatter();
         final int k=6+1;
         for (int i=k; i<=101+k; i+=10) {
        	matrixOrder = i;
     		matrixA = new double[matrixOrder] [matrixOrder];
     		columnX = new double[matrixOrder];
     		exact�olumnX = new double[matrixOrder];
     		columnF = new double[matrixOrder];
     		 for (int i1=0; i1<matrixOrder; i1++)
     			 exact�olumnX[i1]=i1+1;
             for (int i1=0; i1<matrixOrder; i1++)
                  for (int j=0; j<matrixOrder; j++)
                 	 matrixA[i1][j] = 100-(Math.random()*200);
             for(int i1 = 0; i1 < matrixOrder; i1++)
                 for(int j = 0; j < matrixOrder; j++)
                     columnF[i1] += exact�olumnX[j] * matrixA[i1][j];
             try {
             	columnX=gaussSelectionByColumn(matrixOrder,Arrays.stream(matrixA).map(double[]::clone).toArray(double[][]::new),columnF);
             }catch(ArithmeticException e){
              	System.out.println(e.getMessage());
              	System.exit(0);
              }
             
             //�����������
             normXX=Math.abs(exact�olumnX[0]-columnX[0]);
             normX = Math.abs(exact�olumnX[0]);
             for (int i1=1; i1<matrixOrder; i1++) {
             	if(Math.abs(exact�olumnX[i1]-columnX[i1])>normXX)
             		normXX=Math.abs(exact�olumnX[i1]-columnX[i1]);
             	if(Math.abs(exact�olumnX[i1])>normX)
             		normX=Math.abs(exact�olumnX[i1]);
             }
             accuracy=normXX/normX*100;
             f.format("%5d%36.14f%%%n",i, accuracy);
         }
         System.out.print(f); 
	}
}