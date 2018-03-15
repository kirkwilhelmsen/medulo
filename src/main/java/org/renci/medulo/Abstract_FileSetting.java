package org.renci.medulo;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
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
		public void load(boolean testDimensions,Set<String> allGenes) {
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
			Set<String> geneSet = new LinkedHashSet<String>();
			if(allGenes.contains(removeQuotes(line[0])))geneSet.add(removeQuotes(line[0]));
			System.out.println("Beginning first line data \t"+(line.length-1) + "\t" + removeQuotes(line[0]) + "\t" + line[1]);
			String nextline = null;
			boolean infoLinewritten =false;
			while((nextline = in.nextLine())!=null)	{
				String [] l = nextline.split(delimiter);
				if(allGenes.contains(removeQuotes(l[0])))geneSet.add(removeQuotes(l[0]));
				if(!infoLinewritten && geneSet.size()>0 && geneSet.size()%20==0) {
					logger.debug(geneSet.size() +" expected genes found");
					infoLinewritten = true;
				}
				if(infoLinewritten && geneSet.size()>0 && (geneSet.size()-1)%20==0) {
					infoLinewritten = false;
				}
			}
			in.close();
			List<String> geneList = new ArrayList<String>();
			geneList.addAll(geneSet);
			Collections.sort(geneList);
			yAxisRowNames = geneList.toArray(new String [geneList.size()]);
			for(int i = 0; i<yAxisRowNames.length;i++)yAxisRowNames2Index.put(geneList.get(i),i);
			data = new Double [yAxisRowNames.length][xAxisColNames.length];
			
			logger.info("Matrix size:  " + geneSet.size()+ " X " +xAxisColNames.length);
			in = new CHATBufferedFileReader(inFile);
			nextline = in.nextLine();
			int lineCt=0;
			int linesInterrogated=0;
			while((nextline=in.nextLine())!=null){
				linesInterrogated++;
				String [] s = nextline.split(delimiter);
				if(yAxisRowNames2Index.containsKey(removeQuotes(s[0]))){
					lineCt++;
					for(int i =1; i<s.length;i++) {
	//					inputValue = "At i: " + (i-1) + " row: " + lineCt + " value: " + s[i]; 
						try {
							data[yAxisRowNames2Index.get(removeQuotes(s[0]))][i - 1] = Double.valueOf(removeQuotes(s[i]));
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
				}
				if(linesInterrogated%1000==0)logger.info("Lines read:\t" + linesInterrogated +"\tLines used:\t" + lineCt);
			}
			in.close();
			logger.info(lineCt +" gene expression data found of the " + allGenes.size() + " expected");
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

	public class Evolution{
		public double scoreThreshold;
		public double probMutation;
		public double fractChange;
		public boolean pareto =true;
		public boolean sortAngel = false;
		public double paretoAngle =5;
		public double angle;
		public double wt0;
		public double wt1;
		public int minProgeney;
		//public int genTarget;
		public int genomes2Save;
		Map<String,Set<String>> ggIn ;
		Map<String,Set<String>> g2gr = new LinkedHashMap<String,Set<String>>();
		double[][] expGeneByCell;
		public BigFrame exp;
		public EvolutionModel startModel = new EvolutionModel();
		public List<EvolutionModel> newModels = new ArrayList<EvolutionModel>();
		public List<EvolutionModel> models = new ArrayList<EvolutionModel>();
		Random r = new Random(System.currentTimeMillis());
 		public String [] grOrder;
 		public String [] gOrder;
 		public String [] cOrder;
// 		public CompairEvolutionModel cmp;
// 		public CompairEvolutionModelAngle cmpa;
 		public StringBuffer line = new StringBuffer();
 		public boolean verbose1=false;
 		public Evolution(Map<String, Set<String>> gg, BigFrame exp/*, CompairEvolutionModel cmp, CompairEvolutionModelAngle cmpa*/) {
			super();
			this.ggIn = gg;
			this.exp = exp;
//			this.cmp= cmp;
			for(String gr:gg.keySet()) {
				for(String g: gg.get(gr)) {
					if(!g2gr.containsKey(g))g2gr.put(g, new LinkedHashSet<String>());
					g2gr.get(g).add(gr);
				}
			}
			gOrder = exp.yAxisRowNames2Index.keySet().toArray(new String[exp.yAxisRowNames2Index.size()]);
			grOrder = gg.keySet().toArray(new String[gg.keySet().size()]);
			cOrder =exp.getxAxisColNames2Index().keySet().toArray(new String[exp.getxAxisColNames2Index().keySet().size()]);
			expGeneByCell = new double [gOrder.length][cOrder.length];
			for(int g = 0;g<gOrder.length;g++) {
				for(int c=0;c<cOrder.length;c++) {
					expGeneByCell[g][c] = exp.data[exp.yAxisRowNames2Index.get(gOrder[g])][exp.getxAxisColNames2Index().get(cOrder[c])];
				}
			}
			exp=null;
			System.gc();
		}
		public void makeInitalModel() {
			startModel.betas=new double[gOrder.length][grOrder.length]; 
			for(int gr = 0;gr<grOrder.length;gr++) {
				for(int g= 0;g<gOrder.length;g++) {
					startModel.betas[g][gr] =g2gr.get(gOrder[g]).contains(grOrder[gr])?1d:-1d;
				}
			}
			normalizeBetas(startModel);
			calcModel(startModel);
			if(verbose1)System.out.println( "metrix scores for initial model:\t" +  startModel.metrix[0].toString() + "\t" + startModel.metrix[1].toString());
			models.add(startModel);
		}
		private void normalizeBetas(EvolutionModel curModel) {
			double sumBetas;
			for(int gr = 0;gr<grOrder.length;gr++) {
				sumBetas=0;
				for(int g= 0;g<gOrder.length;g++) {
					sumBetas += Math.abs(startModel.betas[g][gr]);
				}
				sumBetas = sumBetas/((double)gOrder.length);
				for(int g= 0;g<gOrder.length;g++) {
					startModel.betas[g][gr] = startModel.betas[g][gr]/sumBetas;
				}
			}
		}
		private void calcModel(EvolutionModel curModel) {
			curModel.scores = new double[cOrder.length][grOrder.length]; 
			curModel.identity = new String [cOrder.length];
			for(int gr = 0;gr<grOrder.length;gr++) {
				for(int c =0;c<cOrder.length;c++ ) {
					double curScore = 0;
					for(int g =0; g<gOrder.length;g++) {
						try {
							curScore += curModel.betas[g][gr] * expGeneByCell[g][c];
						} catch (Exception e) {
//							logger.warn("There is some unexpected missing data (likely expression) for cell:\t" + c + " gene:\t" + g +"\t" + gOrder[g]);
						}
					}
					curModel.scores[c][gr]=curScore;
				}
			}
			//normalize by group
			double cMax;
			Double [] metrix = {0d,0d};
			for(int m=0;m<grOrder.length;m++) {
				cMax = Double.NEGATIVE_INFINITY;
				for( int c = 0;c<curModel.scores.length;c++) {
					cMax=Double.max(cMax, curModel.scores[c][m]);
				}
				for( int c = 0;c<curModel.scores.length;c++) {
					curModel.scores[c][m]= curModel.scores[c][m]/cMax;
				}
			}
			int intGr = -9;
			double cellTotal;
			for( int c = 0;c<curModel.scores.length;c++) {
				cMax=Double.NEGATIVE_INFINITY;
				cellTotal=0;
				for(int gr = 0;gr<grOrder.length;gr++) {
					cellTotal+= curModel.scores[c][gr];
					if(curModel.scores[c][gr]>cMax) {
						intGr = gr;
						cMax = curModel.scores[c][gr];
					}
				}
				if(cMax>scoreThreshold) {
					curModel.identity[c]=grOrder[intGr];
					cellTotal-= cMax;
					metrix[0] +=cMax;
					metrix[1] -= cellTotal;
				}else {
					curModel.identity[c]="null";
					metrix[0] += 1-cMax;
					metrix[1] -=cellTotal;
				}
			}
			curModel.metrix= metrix;
			curModel.euclidAndAngle = new Double[2];
			curModel.euclidAndAngle[0]=Math.sqrt(Math.pow(metrix[0],2)+Math.pow(metrix[1], 2));
			curModel.euclidAndAngle[1]=Math.atan(metrix[0]/metrix[1])/Math.PI*180d;
		}
		public void makeModels(int x) {
			logger.info("Making new " + x + " for each existing saved model");
			int ct =0;
			if(x<minProgeney)x = minProgeney;
			for(EvolutionModel cur :models ) {
				for(int rep=0;rep<x;rep++) {
					ct++;
					EvolutionModel newM = new EvolutionModel();
					newM.betas = cur.betas.clone();
					for(int gr=0;gr<grOrder.length;gr++) {
						for(int c=0;c<exp.getyAxisRowNames().length;c++) {
							if(r.nextFloat()<probMutation) {
									newM.betas[c][gr] += r.nextGaussian()*fractChange;
							}
						}
					}
					normalizeBetas(newM);
					calcModel(newM);
					newModels.add(newM);
					if(verbose1)System.out.println( "metrix scores for new model:\t" + ct + "\t" + newM.metrix[0].toString() + "\t" + newM.metrix[1].toString());
					if(ct%100==0)logger.info(ct + " new models made");
				}
			}
			logger.info("Finish making new models");
		}
		public class EvolutionModel implements Comparable<EvolutionModel>{
			double[][] betas;
			double scores[][];
			String [] identity;
			Double [] metrix;
			Double [] euclidAndAngle;
			@Override
			public int compareTo(EvolutionModel o) {
				if(!sortAngel) {
					Double temp1 =(this.metrix[0]*wt0 + this.metrix[1]*wt1);
					Double temp2 =(o.metrix[0]*wt0 + o.metrix[1]*wt1);
					return -temp1.compareTo(temp2);
				}else {
					int x = this.euclidAndAngle[1].compareTo(o.euclidAndAngle[1]);
					if(x==0) {	
						return this.euclidAndAngle[0].compareTo(o.euclidAndAngle[0]);
					}else {
						return x;
					}
				}
			}
		}
		public void calcNewModels() {
			for(EvolutionModel cur:newModels) {
				calcModel(cur);
				models.add(cur);
			}
			newModels.clear();
		}
		public void evolve(int gens,int progeny) {
			for(int reps = 0;reps<gens;reps++) {
				logger.info("Starting generation: " + reps);
				int p2 = (int)((double) progeny * (double)genomes2Save/(double) models.size());
				if(p2<10)p2=10;
				makeModels(p2);
				selectNextGen();
			}
		}
		public void selectNextGen() {
			if (!pareto) {
				models.addAll(newModels);
				newModels.clear();
				sortAngel=false;
				Collections.sort(models);
				for (int i = 0; i < genomes2Save; i++) {
					newModels.add(models.get(i));
				}
				models.clear();
				for (EvolutionModel ev : newModels) {
					/*				System.out.println(ev.metrix[0]);
									System.out.println(ev.metrix[1]);
									System.out.println(startModel.metrix[0]);
									System.out.println(startModel.metrix[1]);
									System.out.println(cmp.compare(ev, startModel));*/
					if (ev.compareTo(startModel) <= 0 || r.nextFloat() < 0.1)
						models.add(ev);
				}
				Collections.sort(models);
				newModels.clear();
				logger.info(models.size() + " models saved");
				logger.info("metrix scores for initial model\t" + startModel.metrix[0].toString() + "\t"
						+ startModel.metrix[1].toString());
				logger.info("metrix scores for best model:\t" + models.get(0).metrix[0].toString() + "\t"
						+ models.get(0).metrix[1].toString());
			}else {
				sortAngel=true;
				models.addAll(newModels);
				newModels.clear();
				Collections.sort(models);
				int ct = (int)((double)genomes2Save/((double)90/(double)paretoAngle));
				if(ct<10)ct=10;
				for(double a=0;a<90;a+=paretoAngle) {
					List<Double> euclidDist = new ArrayList<Double>();
					for(int i=0; i< models.size() ; i++) {
						if( models.get(i).euclidAndAngle[1].doubleValue()>=a && models.get(i).euclidAndAngle[1].doubleValue()<=a + paretoAngle) {
							euclidDist.add(models.get(i).euclidAndAngle[0]);
						}
					}
					if(euclidDist.size()>0) {
						Collections.sort(euclidDist);
						double eu;
						if(euclidDist.size()<=ct) {
							eu = Double.MIN_VALUE;
						}else {
							eu = euclidDist.get(euclidDist.size()-ct);
						}
						logger.info("For angle " + a + " to " + (a +paretoAngle) + " has " + euclidDist.size() + " of these up to " + ct + " will be saved");
						for(int i=0; i< models.size() ; i++) {
							if( models.get(i).euclidAndAngle[0]>= eu &&models.get(i).euclidAndAngle[1].doubleValue()>a && models.get(i).euclidAndAngle[1].doubleValue()<=(a + paretoAngle)) {
								newModels.add(models.get(i));
							}
						}
					}
				}
				sortAngel=false;
				models.clear();
				models.addAll(newModels);
				newModels.clear();
				Collections.sort(models);
				newModels.clear();
				logger.info(models.size() + " models saved");
				logger.info("metrix scores for initial model based on metrix\t" + startModel.metrix[0].toString() + "\t"
						+ startModel.metrix[1].toString());
				logger.info("metrix scores for best model:\t" + models.get(0).metrix[0].toString() + "\t"
						+ models.get(0).metrix[1].toString());
			}
		}
		public void write(File file) {
			CHATBufferedFileWriter out = new CHATBufferedFileWriter();
			out.open(file.getAbsolutePath());
			if(pareto) {
				sortAngel=true;
			}else {
				sortAngel=false;
			}
			Collections.sort(models);
			for(int m=0;m<models.size();m++) {
				EvolutionModel cur = models.get(m);
				writeModel(m, out,cur);
			}
			out.close();
		}
		private void writeModel(int m, CHATBufferedFileWriter out, EvolutionModel cur) {
			line.delete(0, line.capacity());
			line.append("Start\t").append(m);
			out.writeString(line.toString(),true);
			line.delete(0, line.capacity());
			line.append(cur.metrix[0]).append("\t").append(cur.metrix[1]);
			out.writeString(line.toString(),true);
			line.delete(0, line.capacity());
			line.append(cur.euclidAndAngle[0]).append("\t").append(cur.euclidAndAngle[1]);
			out.writeString(line.toString(),true);
			for(int c=0;c<cOrder.length;c++) {
				line.delete(0, line.capacity());
				line.append(cOrder[c]).append("\t");
				for(int gr=0;gr<grOrder.length;gr++) {
					line.append(cur.scores[c][gr]).append("\t");
				}
				line.append(cur.identity[c]);
				
			}
			
			out.writeString("Betas");
			
			for(int g=0;g<gOrder.length;g++) {
				line.delete(0, line.capacity());
				line.append(gOrder[g]).append("\t");
				for(int gr=0;gr<grOrder.length;gr++) {
					line.append(cur.betas[g][gr]).append("\t");
				}
				out.writeString(line.toString(),true);
			}
			line.delete(0, line.capacity());
			line.append("End\t").append(m);
			out.writeString(line.toString(),true);

		}
		public void setScoreThreshold(double scoreThreshold) {
			this.scoreThreshold = scoreThreshold;
		}
		public void setFractioChange(double d) {
			this.fractChange=d;
		}
		public void setGenomes2Save(int genomes2Save) {
			this.genomes2Save = genomes2Save;
		}
		public void setProbMutation(double probMutation) {
			this.probMutation = probMutation;
		}
		public void setPareto(boolean pareto) {
			this.pareto = pareto;
		}
		public void setAngle(double angle) {
			this.angle = angle;
		}
		public void setWt0(double wt0) {
			this.wt0 = wt0;
		}
		public void setWt1(double wt1) {
			this.wt1 = wt1;
		}
		public void setMinProgeney(int minProgeney) {
			this.minProgeney = minProgeney;
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
