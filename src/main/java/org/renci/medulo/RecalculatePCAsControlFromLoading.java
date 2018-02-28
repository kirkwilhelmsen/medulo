package org.renci.medulo;

import java.util.Map;

public class RecalculatePCAsControlFromLoading extends Abstract_FileSetting implements Runnable {
	int cellCountToCheck = 10;
	
	@Override
	public void run() {
		boolean proceed= true;
/*		if(!testDimentionsPCALoads()) {
			logger.error("Problem with PCA loads file");
			proceed = false;
		}
		if(!testDimentionsExpresion(getTestScaledExpression()," ")) {
			logger.error("Problem with test scaled expression file");
			proceed = false;
		}
		if(!testDimentionsExpresion(getControlScaledExression()," ")) {
			logger.error("Problem with control scaled expression file");
			proceed = false;
		}
		if(!testDimensionPCScores(getControlPCScores(),"\t")) {
			logger.error("Problem with control PC scores file");
			proceed = false;
		};*/
		if(!proceed) {
			System.exit(-1);
		}
//		loadPCALoadings();
//		loadBigMatrix(getControlScaledExression(),true,true,"\t");
		loadBigMatrix(getControlPCScores(),true,false,"\t");
		Map<String,float[]> sc = getControlPCs();
		String[] id = getControlPCId();
//		loadBigMatrix(getTestScaledExpression(), false, false," ");
	}

	public static void main(String[] args) {
		RecalculatePCAsControlFromLoading rcpca = new RecalculatePCAsControlFromLoading();
		rcpca.run();
	}

}
