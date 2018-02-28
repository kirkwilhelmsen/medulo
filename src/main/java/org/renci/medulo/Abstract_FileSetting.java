package org.renci.medulo;

import java.io.File;
import java.util.LinkedHashMap;
import java.util.Map;

import org.renci.medulo.tools.CHATBufferedFileReader;
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
	
	private String [] pcNames = {};
	private Map<String,double[]> pcLoadings = new LinkedHashMap<String,double[]>();
	private String [] controlCellId;
	private Map<String,Integer> controlExpGenes = new LinkedHashMap<String,Integer>();
	private float [][] controlExp;//[cell][gene]
	private Map<String,float[]> controlPCs= new LinkedHashMap<String, float []>();
	private String[] controlPCId;
	private String [] testCellId;
	Map<String,Integer> testExpGenes = new LinkedHashMap<String,Integer>();
	float [][] testExp;//[cell][gene]
	private Map<String,float[]> testPCs= new LinkedHashMap<String, float []>();
	private String[] testPCId;
	
	static final Logger logger = LoggerFactory.getLogger(Abstract_FileSetting.class);

	
	public void loadPCALoadings() {
		logger.info("Loading PCA loadings in " + pcaLoadings);
		File inPCALoading = new File(pcaLoadings);
		CHATBufferedFileReader in = new CHATBufferedFileReader(inPCALoading);
		String line = in.nextLine();
//		System.out.println(line);
		pcNames = line.split("\t");
		for(int i=0;i<pcNames.length; i++)pcNames[i]=removeQuotes(pcNames[i]);
		while((line=in.nextLine())!=null){
//			System.out.println(line);
			String [] f = line.split("\t");
			double [] l = new double [f.length-1];
			for(int i=1;i<(f.length); i++)l[i-1] = Double.valueOf(f[i]);
			pcLoadings.put(removeQuotes(f[0]), l);
		}
		logger.info("Finished loading PCA loadings");
	}
	public void loadBigMatrix(String fileName, boolean isControl, boolean isExpresion, String delimiter) {
		logger.info("Loading control scaled expression in : " +fileName);
		File inScaledExression = new File(fileName);
		CHATBufferedFileReader in = new CHATBufferedFileReader(inScaledExression);
		if(isControl) {
			if(isExpresion) {
				controlCellId = in.nextLine().split(delimiter);
				for(int i=0;i<controlCellId.length; i++)controlCellId[i]=removeQuotes(controlCellId[i]);
				System.out.println("Length header for columns \t" + controlCellId.length + " first two col header names \t" + controlCellId[0] + "\t" + controlCellId[1]);
			}else {
				controlPCId = in.nextLine().split(delimiter);
				for(int i=0;i<controlPCId.length; i++)controlPCId[i]=removeQuotes(controlPCId[i]);
				System.out.println("Length header for columns \t" + controlPCId.length + " first two col header names \t" + controlPCId[0] + "\t" + controlPCId[1]);
			}
		}else {
			if(isExpresion) {
				testCellId = in.nextLine().split(delimiter);
				for(int i=0;i<testCellId.length; i++)testCellId[i]=removeQuotes(testCellId[i]);
				System.out.println("Length header for columns \t" + testCellId.length + " first two col header names \t" + testCellId[0] + "\t" + testCellId[1]);
			}else {
				testPCId = in.nextLine().split(delimiter);
				for(int i=0;i<testPCId.length; i++)testPCId[i]=removeQuotes(testPCId[i]);
				System.out.println("Length header for columns \t" + testPCId.length + " first two col header names \t" + testPCId[0] + "\t" + testPCId[1]);
			}
		}
		String [] line = in.nextLine().split(delimiter);
		System.out.println("Beginning first line data \t"+(line.length-1) + "\t" + removeQuotes(line[0]) + "\t" + line[1]);
		String nextline = null;
		int lineCt=0;
		while((nextline=in.nextLine())!=null){
			lineCt++;
		}
		in.close();
		if(isControl) {
			if(isExpresion) {
				controlExp = new float [controlCellId.length][lineCt +1];
				logger.info("Matrix size: " +controlCellId.length + " X " + lineCt+1);
			}else {
				logger.info("Matrix size: " +controlPCId.length + " X " + lineCt+1);
			}
		}else {
			if(isExpresion) {
				testExp = new float [testCellId.length][lineCt +1];
				logger.info("Matrix size: " +testCellId.length + " X " + lineCt+1);
			}else {
				logger.info("Matrix size: " +testPCId.length + " X " + lineCt+1);
			}
		}
		in = new CHATBufferedFileReader(inScaledExression);
		nextline = in.nextLine();
		lineCt=-1;
		while((nextline=in.nextLine())!=null){
			lineCt++;
//			String currentLineNumber = "Current line number :" + lineCt;
			String [] s = nextline.split(delimiter);
//			String currentLineLength = "Current line length :" + s.length;
			if(isControl) {
				if(isExpresion) {
					controlExpGenes.put(removeQuotes(s[0]), lineCt);
				}else {
					float[] p = new float[controlPCId.length];
					for(int i=1;i<s.length;i++)p[i-1]=Float.valueOf(s[i]);
					controlPCs.put(removeQuotes(s[0]), p);
				}
			}else {
				if(isExpresion) {
					testExpGenes.put(removeQuotes(s[0]), lineCt);
				}else {
					float[] p = new float[testPCId.length];
					for(int i=1;i<s.length;i++)p[i-1]=Float.valueOf(s[i]);
					testPCs.put(removeQuotes(s[0]), p);
				}
			}
//			String currentGene = s[0];
//			String inputValue;
			for(int i =1; i<s.length;i++) {
//				inputValue = "At i: " + (i-1) + " row: " + lineCt + " value: " + s[i]; 
				try {
					if(isControl) {
						if(isExpresion) controlExp[i - 1][lineCt] = Float.valueOf(s[i]);
					}else {
						if(isExpresion) testExp[i - 1][lineCt] = Float.valueOf(s[i]);
					}
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
			if(lineCt%1000==0)logger.info("Lines read:\t" + lineCt);
		}
		in.close();
		logger.info("Finish loading control scaled expression");
	}
	
	public boolean testDimentionsExpresion(String stringFileName, String delimiter) {
		logger.info("Test dimensions scaled expression in " + stringFileName);
		File inScaledExression = new File(stringFileName);
		CHATBufferedFileReader in = new CHATBufferedFileReader(inScaledExression);
		String [] cellId = in.nextLine().split(delimiter);
		System.out.println("Width of table and first header col names \t"+ cellId.length + "\t" + removeQuotes(cellId[0]) + "\t" + removeQuotes(cellId[1]));
		String [] line = in.nextLine().split(delimiter);
		System.out.println("Length first line and first 2 elements\t "+(line.length-1) + "\t" + removeQuotes(line[0]) + "\t" + line[1]);
		String nextline = null;
		int geneIndex=0;
		String lastMarker=null;
		int lastlength = 0;
		boolean out = true;
		
		while((nextline=in.nextLine())!=null){
			geneIndex++;
			String [] s = nextline.split(delimiter);
			if(cellId.length!=(s.length-1)) {
				if(s.length > controlCellId.length + 1) {
					System.out.println(geneIndex + "\t" + s[0] +"\t" + s.length + "\t" + s[controlCellId.length] + "\t" +lastMarker + "\t" +lastlength);
				}else {
					System.out.println(geneIndex + "\t" + s[0] +"\t" + s.length + "\t" + s[s.length-1] + "\t" +lastMarker + "\t" +lastlength);
				}
				out=false;
			}else{
				lastMarker = s[0];
				lastlength=s.length;
			}
			if(geneIndex%1000==0)logger.info("Lines read:\t" + geneIndex);
		}
		in.close();
		logger.info("Finish test dimensions scaled expression in " + stringFileName);
		return out;
	}
	public boolean testDimensionPCScores(String stringFileName, String delimiter) {
		logger.info("Test dimensions PC scores file in " + stringFileName);
		File inPCScores = new File(stringFileName);
		CHATBufferedFileReader in = new CHATBufferedFileReader(inPCScores);
		String [] pcId = in.nextLine().split(delimiter);
		System.out.println("Width of table and first header col names \t"+ pcId.length + "\t" + removeQuotes(pcId[0]) + "\t" + removeQuotes(pcId[1]));
		String [] line = in.nextLine().split(delimiter);
		System.out.println("Length first line and first 2 elements\t "+(line.length-1) + "\t" + removeQuotes(line[0]) + "\t" + line[1]);
		String nextline = null;
		int cellIndex=0;
		String lastCell=null;
		int lastlength = 0;
		boolean out = true;
		
		while((nextline=in.nextLine())!=null){
			cellIndex++;
			String [] s = nextline.split(delimiter);
			if(pcId.length!=(s.length-1)) {
				if(s.length > controlCellId.length + 1) {
					System.out.println(cellIndex + "\t" + s[0] +"\t" + s.length + "\t" + s[pcId.length] + "\t" +lastCell + "\t" +lastlength);
				}else {
					System.out.println(cellIndex + "\t" + s[0] +"\t" + s.length + "\t" + s[s.length-1] + "\t" +lastCell + "\t" +lastlength);
				}
				out=false;
			}else{
				lastCell = s[0];
				lastlength=s.length;
			}
			if(cellIndex%1000==0)logger.info("Lines read:\t" + cellIndex);
		}
		in.close();
		logger.info("Finish test dimensions PC scores in " + stringFileName);
		return out;
	}
	public boolean testDimentionsPCALoads() {
		logger.info("Test dimensions data in PCA loading file");
		File inPCALoading = new File(pcaLoadings);
		CHATBufferedFileReader in = new CHATBufferedFileReader(inPCALoading);
		String [] pcaID = in.nextLine().split("\t");
		System.out.println("First 2 entries header\t"+  pcaID.length + "\t" + removeQuotes(pcaID[0]) + "\t" + removeQuotes(pcaID[1]));
		String [] line = in.nextLine().split("\t");
		System.out.println("Length first line and first 2 elements\t "+ (line.length-1) + "\t" + removeQuotes(line[0]) + "\t" + line[1]);
		String nextline = null;
		int pcaIndex=-1;
		String lastMarker=null;
		int lastlength = 0;
		boolean out = true;
		in.close();
		in = new CHATBufferedFileReader(inPCALoading);
		in.nextLine();
		while((nextline=in.nextLine())!=null){
			pcaIndex++;
			String [] s = nextline.split("\t");
			if(pcaID.length!=s.length-1) {
				if(s.length > pcaID.length + 1) {
					System.out.println(pcaIndex + "\t" + s[0] +"\t" + s.length + "\t" + s[pcaID.length] + "\t" +lastMarker + "\t" +lastlength);
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
		logger.info("Finish test dimensions data in PCA loading file ");
		return out;
	}
	public String removeQuotes(String in ) {
		return in.replace("\"", "");
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
	public String[] getPcNames() {
		return pcNames;
	}
	public void setPcNames(String[] pcNames) {
		this.pcNames = pcNames;
	}
	public Map<String, double[]> getPcLoadings() {
		return pcLoadings;
	}
	public void setPcLoadings(Map<String, double[]> pcLoadings) {
		this.pcLoadings = pcLoadings;
	}
	public String[] getControlCellId() {
		return controlCellId;
	}
	public void setControlCellId(String[] controlCellId) {
		this.controlCellId = controlCellId;
	}
	public Map<String, Integer> getControlExpGenes() {
		return controlExpGenes;
	}
	public void setControlExpGenes(Map<String, Integer> controlExpGenes) {
		this.controlExpGenes = controlExpGenes;
	}
	public float[][] getControlExp() {
		return controlExp;
	}
	public void setControlExp(float[][] controlExp) {
		this.controlExp = controlExp;
	}
	public String[] getTestCellId() {
		return testCellId;
	}
	public void setTestCellId(String[] testCellId) {
		this.testCellId = testCellId;
	}
	public Map<String, Integer> getTestExpGenes() {
		return testExpGenes;
	}
	public void setTestExpGenes(Map<String, Integer> testExpGenes) {
		this.testExpGenes = testExpGenes;
	}
	public float[][] getTestExp() {
		return testExp;
	}
	public void setTestExp(float[][] testExp) {
		this.testExp = testExp;
	}
	public String getTestPCScores() {
		return testPCScores;
	}
	public void setTestPCScores(String testPCScores) {
		this.testPCScores = testPCScores;
	}
	public Map<String, float[]> getControlPCs() {
		return controlPCs;
	}
	public void setControlPCs(Map<String, float[]> controlPCs) {
		this.controlPCs = controlPCs;
	}
	public String[] getControlPCId() {
		return controlPCId;
	}
	public void setControlPCId(String[] controlPCId) {
		this.controlPCId = controlPCId;
	}
	public Map<String, float[]> getTestPCs() {
		return testPCs;
	}
	public void setTestPCs(Map<String, float[]> testPCs) {
		this.testPCs = testPCs;
	}
	public String[] getTestPCId() {
		return testPCId;
	}
	public void setTestPCId(String[] testPCId) {
		this.testPCId = testPCId;
	}
	public static Logger getLogger() {
		return logger;
	}
	
	

}
