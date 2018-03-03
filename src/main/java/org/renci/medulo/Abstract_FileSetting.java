package org.renci.medulo;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.stream.Collectors;

import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.stat.StatUtils;
import org.renci.medulo.tools.CHATBufferedFileReader;
import org.renci.medulo.tools.CHATBufferedFileWriter;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public abstract class Abstract_FileSetting {
	
	private String pcaLoadings = "GeneLoadings.VehCells.PCs.txt";
	private String controlLouvainIds = "VehicleCells.Louvain.IDs.txt";
	private String controlPCScores = "VehicleCells.PC.Scores.txt";
	private String testPCScores = "VismodegibCells.PC.Scores.txt";
	private String controlScaledExression = "VehicleCells.Scaled.UMIs.txt";
	private String controlTSNEScores = "VehicleCells.tSNE.Scores.txt";
	private String testScaledExpression = "VismodegibCells.Scaled.UMIs.txt";
	
	static final Logger logger = LoggerFactory.getLogger(Abstract_FileSetting.class);

	
	public BigFrame calcPC(BigFrame exp, BigFrame pcLoadings) {
		BigFrame out = new BigFrame("newPCs", "", pcLoadings.getXAxisName(), exp.getXAxisName(),"\t");
		out.setxAxisColNames(pcLoadings.xAxisColNames);
		out.setxAxisColNames2Index(pcLoadings.getxAxisColNames2Index());
		out.setyAxisRowNames(exp.getxAxisColNames());
		out.setyAxisRowNames2Index(exp.getxAxisColNames2Index());
		Double [][] outData = new Double [out.getyAxisRowNames().length][out.getxAxisColNames().length];
		String outCurCell;
		int inCurCellExpCol;
		int inCurPCLoadCol;
		String outCurPC;
		int curInGeneExpY;
		int curInGeneLoadY;
		double curScore= 0d;
		for(int y =0; y< outData.length; y++ ) {
			outCurCell=out.getyAxisRowNames()[y];
			inCurCellExpCol = exp.getxAxisColNames2Index().get(outCurCell);
			for(int x =0; x< outData[y].length; x++ ) {
				outCurPC = out.getxAxisColNames()[x];
				inCurPCLoadCol = pcLoadings.getxAxisColNames2Index().get(outCurPC);
				curScore = 0d;
				for(String g : pcLoadings.getyAxisRowNames() ) {
					curInGeneExpY=exp.getyAxisRowNames2Index().get(g);
					curInGeneLoadY = pcLoadings.getyAxisRowNames2Index().get(g);
					try {
						curScore += exp.getData()[curInGeneExpY][inCurCellExpCol] * pcLoadings.getData()[curInGeneLoadY][inCurPCLoadCol];
					} catch (Exception e) {
						e.printStackTrace();
					}
				}
				outData[y][x]=curScore;
			}
		}
		out.setData(outData);
		return out;
	}

	public class kNearestNeighbor {
		/**
		 * Adapted from kNearestNeighbor
		 * Created by Leon Aburime on 7/15/2014.
		 * Built using IntelliJ
		 */
		RealMatrix m;
		RealMatrix sm;

		//Value we are trying to predict whether default will happen
	    private int k;
	    private int numberPC;
	    private BigFrame tSNE;
	    private BigFrame cluster;
	    private BigFrame controlPC;
	    private boolean standardize;
	    private boolean allowAmbigousCluster;
	    
	    StringBuffer out = new StringBuffer();
	    StringBuffer out2 = new StringBuffer();
	    
	    public kNearestNeighbor(BigFrame controlPC, int numberPC,int k, BigFrame tSNE, BigFrame cluster,boolean standardize, boolean allowAmbigousCluster) {
	    	this.controlPC=controlPC;
	        double [][] input= new double[controlPC.getData().length][controlPC.getData()[0].length];
	    	for(int row = 0; row< input.length; row++) {
	    		for(int col=0;col<input[0].length ; col++) {
	    			input[row][col]=controlPC.getData()[row][col].doubleValue();
	    		}
	    	}
	    	this.m = MatrixUtils.createRealMatrix(input);
	    	if(standardize) {
	    		sm = m.copy();
	    		standardizeMatrix();
	    	}
	    	this.k = k;
	    	this.numberPC=numberPC;
	    	this.tSNE=tSNE;
	    	this.cluster=cluster;
	    	this.standardize = standardize;
	    	this.allowAmbigousCluster = allowAmbigousCluster;
	    }

	    public  String  calculate_KNN(Double [] predict, String cellID) {
	    	out.delete(0, out.capacity());
	    	out.append("\""+cellID +"\"").append("\t");
	    	double [] temp = new double [predict.length];
	    	for(int i = 0;i<predict.length;i++)temp[i]=predict[i];
	    	if(standardize)temp = 	 standardizeArray( m, temp );//Standardize prediction values
	        List<Integer> answers;
	        if(standardize) {
	        	answers= euclideanDistance( temp,sm,cellID);
	    	}else {
	    		answers= euclideanDistance( temp,sm,cellID);
	    	}	
	        Double[] t = calcAverageTsne(answers);
	        for(int i =0;i<t.length;i++)out.append(t[i].toString()).append("\t");
	        out.append(calcClass(answers)).append("\t");
	        return out.toString();
	    }


	    private String calcClass(List<Integer> answers) {
	    	Map<Long,Integer> map = new LinkedHashMap<Long, Integer>();
	    	for(int c : answers) {
				String curCellPCControl = controlPC.getyAxisRowNames()[c];
				int rowCellCluster= cluster.getyAxisRowNames2Index().get(curCellPCControl);
				Long key = Math.round(cluster.getData()[rowCellCluster][0]);
				if(!map.containsKey(key)) {
					map.put(key, 1);
				}else {
					map.put(key, map.get(key)+1);
				}
	    	}
			int maxCt = 0;
			for(Long key:map.keySet()) if(maxCt< map.get(key))maxCt= map.get(key);
			out2.delete(0, out2.capacity());
			for(Long key:map.keySet()) {
				if(maxCt== map.get(key)) {
					if(allowAmbigousCluster) {
						if(out2.length()>0) {
							out2.append(",").append(key.toString());
						}else {
							out2.append(key.toString());
						}
					}else {
						out2.append(key.toString());
						break;
					}
				}
			}
			return out2.toString();
		}

		private Double[] calcAverageTsne(List<Integer> answers) {
			Double [] out = new Double [tSNE.getxAxisColNames().length];
			for(int col = 0;col<tSNE.getxAxisColNames().length;col++) {
				out[col]=0d;
				for(int c : answers) {
					String curCellPCControl = controlPC.getyAxisRowNames()[c];
					int rowCellTsne= tSNE.getyAxisRowNames2Index().get(curCellPCControl);
					out[col]+=tSNE.getData()[rowCellTsne][col];
				}
			}
			for(int col=0;col< out.length;col++)out[col] = out[col]/(double)answers.size();
			return out;
		}

		//Perform Euclidean distance formula to find out the distance
	    //between our prediction value and each row in the matrix
	    public List <Integer> euclideanDistance( double[] y, RealMatrix sm2, String cellID){

	        Map <Double,Integer>map = new HashMap<Double,Integer>();
	        //Lets turn out 'y' value or label into vector for easier math operations
	        RealVector Y = MatrixUtils.createRealVector( y);

	        String currentRefCellName;		
	        if(standardize) {
		        for (int i=0; i<sm.getRowDimension(); i++){
		        	currentRefCellName=controlPC.getyAxisRowNames()[1];
		        	if(!currentRefCellName.equalsIgnoreCase(cellID)) {
			            RealVector vec = sm.getRowVector(i);
			            RealVector sub1 = vec.subtract( Y );	
			            double [] sub  = Arrays.copyOfRange(sub1.toArray(),0,numberPC);
			            //Take square root of sum of square values that were subtracted a line above
			            Double distance = Math.sqrt(StatUtils.sumSq(sub));
			            //Use the distance to each data point(or row) as key with the 'default' option as value
			            map.put( distance  , i );
		        	}
		        }
	        }else {
	        	for (int i=0; i<m.getRowDimension(); i++){
		        	currentRefCellName=controlPC.getyAxisRowNames()[1];
		        	if(!currentRefCellName.equalsIgnoreCase(cellID)) {
			            RealVector vec = m.getRowVector(i);
			            RealVector sub1 = vec.subtract( Y );	
			            double [] sub  = Arrays.copyOfRange(sub1.toArray(),0,numberPC);
			            //Take square root of sum of square values that were subtracted a line above
			            Double distance = Math.sqrt(StatUtils.sumSq(sub));
			            //Use the distance to each data point(or row) as key with the 'default' option as value
			            map.put( distance  , i );
		        	}
		        }
	        }
	        //Now lets sort the map's keys into a set
	        SortedSet<Double> keys = new TreeSet<Double>(map.keySet());
	        List<Integer> neighbors = new ArrayList<Integer>();

	        //For each key add the values in that order into the list
	        for (Double key : keys) {
	            neighbors.add( map.get(key));
	            if(neighbors.size()>=k) {
	            	break;
	            }
	        }
	        return neighbors;
	    }

	    //Standardize our matrix i.e. ..
	    //(each_value - min_val_in_column)/(max_val_in_column - min_val_in_column
	    public void standardizeMatrix(){
	        for (int i=0; i<m.getColumnDimension(); i++){
	            //Get each column as Vector
	            RealVector vec = m.getColumnVector(i);

	            //Lets get the max and the min of the column
	            double max = vec.getMaxValue();
	            double min = vec.getMinValue();

	            vec.mapSubtractToSelf(min);
	            vec.mapDivideToSelf(max - min);

	            //Now lets reset the value back into the matrix
	            sm.setColumnVector(i, vec);
	        }
	    }

	    //Standardize our vector by the matrix i.e. ..
	    //(each_array_value - min_val_in_column)/(max_val_in_column - min_val_in_column
	    public double[] standardizeArray(RealMatrix m, double[] arr) {
	        //Create new array to store the our standardized values in
	        double[] new_arr = new double[arr.length];

	        //Iterate each of the columns of the matrix to find max and min.
	        //Lets use those for calculations on each individual entry of the vector
	        for (int i = 0; i < m.getColumnDimension(); i++) {
	            //Get each column as Vector
	            RealVector vec = m.getColumnVector(i);

	            //Lets get the max and the min of the column
	            double max = vec.getMaxValue();
	            double min = vec.getMinValue();

	            //Store each value in the vector
	            new_arr[i] = (arr[i] - min) / (max - min);
	        }

	        return new_arr;
	    }


	}
	
	public class BigFrame{
		private String frameName;
		private String fileName;
		private String XAxisName;
		private String YAxisName;
		private String delimiter;
		private String [] xAxisColNames;
		private Map<String,Integer> xAxisColNames2Index = new LinkedHashMap<String,Integer>();
		private String [] yAxisRowNames; 
		private Map<String,Integer> yAxisRowNames2Index = new LinkedHashMap<String,Integer>();
		private Double[][] data;

		public BigFrame(String name, String fileName, String xAxisName, String yAxisName, String delimiter) {
			super();
			this.frameName= name;
			this.fileName = fileName;
			XAxisName = xAxisName;
			YAxisName = yAxisName;
			this.delimiter = delimiter;
		}
		public void load(boolean testDimensions) {
			if(testDimensions) {
				if(!testDimentions()) {
					logger.error("The file: " + fileName + " has an irregular shape");
					System.exit(-1);
				}
			}
			logger.info("Loading file: " + fileName + " into " + frameName);
			File inFile = new File(fileName);
			CHATBufferedFileReader in = new CHATBufferedFileReader(inFile);
			xAxisColNames = in.nextLine().split(delimiter);
			for(int i=0;i<xAxisColNames.length; i++)xAxisColNames[i]=removeQuotes(xAxisColNames[i]);
			for(int i=0;i<xAxisColNames.length; i++)xAxisColNames2Index.put(xAxisColNames[i],i);
			if(xAxisColNames.length>=2) {
				System.out.println("Length header for columns \t" + xAxisColNames.length + " first two col header names \t" + removeQuotes(xAxisColNames[0]) + "\t" + removeQuotes(xAxisColNames[1]));
			}else{
				System.out.println("Length header for columns \t" + xAxisColNames.length + " first two col header names \t" + removeQuotes(xAxisColNames[0] ));
			}
			String [] line = in.nextLine().split(delimiter);
			System.out.println("Beginning first line data \t"+(line.length-1) + "\t" + removeQuotes(line[0]) + "\t" + line[1]);
			String nextline = null;
			int lineCt=0;
			while((nextline=in.nextLine())!=null){
				lineCt++;
			}
			in.close();
			data = new Double [lineCt +1][xAxisColNames.length];
			logger.info("Matrix size: " +xAxisColNames.length + " X " + (lineCt+1));
			in = new CHATBufferedFileReader(inFile);
			nextline = in.nextLine();
			lineCt=-1;
			while((nextline=in.nextLine())!=null){
				lineCt++;
//				String currentLineNumber = "Current line number :" + lineCt;
				String [] s = nextline.split(delimiter);
//				String currentLineLength = "Current line length :" + s.length;
				yAxisRowNames2Index.put(removeQuotes(s[0]), lineCt);
//				String currentGene = s[0];
//				String inputValue;
				for(int i =1; i<s.length;i++) {
//					inputValue = "At i: " + (i-1) + " row: " + lineCt + " value: " + s[i]; 
					try {
						data[lineCt][i - 1] = Double.valueOf(removeQuotes(s[i]));
					} catch (Exception e) {
	/*					System.out.println(matrixSize);
						System.out.println(currentLineNumber);
						System.out.println(currentLineLength);
						System.out.println(currentGene);
						System.out.println(inputValue);
						System.out.println("");
	*/					e.printStackTrace();
						System.exit(-1);
					}
				}
				yAxisRowNames = new String[yAxisRowNames2Index.size()];
				for(String r :yAxisRowNames2Index.keySet())yAxisRowNames[yAxisRowNames2Index.get(r)]=r;
				if(lineCt%1000==0)logger.info("Lines read:\t" + lineCt);
			}
			in.close();
			logger.info("Loaded file: " + fileName + " into " + frameName	);
		}
		public boolean testDimentions() {
			logger.info("Test dimensions for file " + fileName);
			File inFile = new File(fileName);
			CHATBufferedFileReader in = new CHATBufferedFileReader(inFile);
			String [] line1 = in.nextLine().split("\t");
			if(line1.length>2) {
				System.out.println("First 2 entries header\t"+  line1.length + "\t" + removeQuotes(line1[0]) + "\t" + removeQuotes(line1[1]));
			}else{
				System.out.println("First entry header\t"+  line1.length + "\t" + removeQuotes(line1[0]) );
			}
			String [] line = in.nextLine().split("\t");
			System.out.println("Length first line and first 2 elements\t "+ (line.length-1) + "\t" + removeQuotes(line[0]) + "\t" + line[1]);
			String nextline = null;
			int pcaIndex=-1;
			String lastMarker=null;
			int lastlength = 0;
			boolean out = true;
			in.close();
			in = new CHATBufferedFileReader(inFile);
			in.nextLine();
			while((nextline=in.nextLine())!=null){
				pcaIndex++;
				String [] s = nextline.split("\t");
				if(line1.length!=s.length-1) {
					if(s.length > line1.length + 1) {
						System.out.println(pcaIndex + "\t" + s[0] +"\t" + s.length + "\t" + s[line1.length] + "\t" +lastMarker + "\t" +lastlength);
					}else {
						System.out.println(pcaIndex + "\t" + s[0] +"\t" + s.length + "\t" + s[s.length-1] + "\t" +lastMarker + "\t" +lastlength);
					}
					out=false;
				}else{
					lastMarker = s[0];
					lastlength=s.length;
				}
				if(pcaIndex%1000==0)logger.info("Lines read:\t" + pcaIndex);
			}
			logger.info("Lines read:\t" + pcaIndex);
			in.close();
			logger.info("Finished loading into " + frameName);
			return out;
		}
		public void dump(String outFileName) {
			logger.info("Dumping " +  frameName + " into " + outFileName);
			File outFile  = new File(outFileName);
			CHATBufferedFileWriter o = new CHATBufferedFileWriter();
			o.open(outFile.getAbsolutePath());
			o.writeString(Arrays.stream(xAxisColNames).map(a -> addQuotes(a)).collect(Collectors.joining("\t")));
			for(int l=0;l<yAxisRowNames.length;l++) {
				o.writeString(addQuotes(yAxisRowNames[l]) + "\t" + Arrays.stream(data[l]).map(a ->a.toString()).collect(Collectors.joining("\t")));
			}
			o.close();
			logger.info("Finishe dumping " +  frameName + " into " + outFileName);
		}
		public String removeQuotes(String in ) {
			return in.replace("\"", "");
		}
		public String addQuotes(String in ) {
			return "\"" + in + "\"";
		}
		public String getFileName() {
			return fileName;
		}
		public void setFileName(String fileName) {
			this.fileName = fileName;
		}
		public String getFrameName() {
			return frameName;
		}
		public void setFrameName(String frameName) {
			this.frameName = frameName;
		}
		public String getXAxisName() {
			return XAxisName;
		}
		public void setXAxisName(String xAxisName) {
			XAxisName = xAxisName;
		}
		public String getYAxisName() {
			return YAxisName;
		}
		public void setYAxisName(String yAxisName) {
			YAxisName = yAxisName;
		}
		public String getDelimiter() {
			return delimiter;
		}
		public void setDelimiter(String delimiter) {
			this.delimiter = delimiter;
		}
		public String[] getxAxisColNames() {
			return xAxisColNames;
		}
		public void setxAxisColNames(String[] xAxisColNames) {
			this.xAxisColNames = xAxisColNames;
		}
		public Map<String, Integer> getxAxisColNames2Index() {
			return xAxisColNames2Index;
		}
		public void setxAxisColNames2Index(Map<String, Integer> xAxisColNames2Index) {
			this.xAxisColNames2Index = xAxisColNames2Index;
		}
		public String[] getyAxisRowNames() {
			return yAxisRowNames;
		}
		public void setyAxisRowNames(String[] yAxisRowNames) {
			this.yAxisRowNames = yAxisRowNames;
		}
		public Map<String, Integer> getyAxisRowNames2Index() {
			return yAxisRowNames2Index;
		}
		public void setyAxisRowNames2Index(Map<String, Integer> yAxisRowNames2Index) {
			this.yAxisRowNames2Index = yAxisRowNames2Index;
		}
		public Double[][] getData() {
			return data;
		}
		public void setData(Double[][] data) {
			this.data = data;
		}
	}
	
	public String getPcaLoadings() {
		return pcaLoadings;
	}

	public void setPcaLoadings(String pcaLoadings) {
		this.pcaLoadings = pcaLoadings;
	}

	public String getControlLouvainIds() {
		return controlLouvainIds;
	}

	public void setControlLouvainIds(String controlLouvainIds) {
		this.controlLouvainIds = controlLouvainIds;
	}

	public String getControlPCScores() {
		return controlPCScores;
	}

	public void setControlPCScores(String controlPCScores) {
		this.controlPCScores = controlPCScores;
	}

	public String getTestPCScores() {
		return testPCScores;
	}

	public void setTestPCScores(String testPCScores) {
		this.testPCScores = testPCScores;
	}

	public String getControlScaledExression() {
		return controlScaledExression;
	}

	public void setControlScaledExression(String controlScaledExression) {
		this.controlScaledExression = controlScaledExression;
	}

	public String getControlTSNEScores() {
		return controlTSNEScores;
	}

	public void setControlTSNEScores(String controlTSNEScores) {
		this.controlTSNEScores = controlTSNEScores;
	}

	public String getTestScaledExpression() {
		return testScaledExpression;
	}

	public void setTestScaledExpression(String testScaledExpression) {
		this.testScaledExpression = testScaledExpression;
	}

	public static Logger getLogger() {
		return logger;
	}

}
