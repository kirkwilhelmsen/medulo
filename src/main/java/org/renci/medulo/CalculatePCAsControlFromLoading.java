package org.renci.medulo;

public class CalculatePCAsControlFromLoading extends Abstract_FileSetting implements Runnable {
	int cellCountToCheck = 10;
	
	@Override
	public void run() {
	/*	BigFrame controlExp = new BigFrame("controlExp",getControlScaledExression(), "control_cell", "gene", "\t");
		controlExp.load(false);
		BigFrame pcLoadings = new BigFrame("pcLoading",getPcaLoadings(), "PC", "gene", "\t");
		pcLoadings.load(false);
		BigFrame reCalControlPC = calcPC(controlExp, pcLoadings);
		reCalControlPC.setFrameName("reCalControlPC");
		reCalControlPC.dump("recalc-" + getControlPCScores())*/;		
		BigFrame testExp = new BigFrame("testExp",getTestScaledExpression(), "test_cell", "gene", "\t");
		testExp.load(true);
		BigFrame pcLoadings = new BigFrame("pcLoading",getPcaLoadings(), "PC", "gene", "\t");
		pcLoadings.load(false);
		BigFrame testPC = calcPC(testExp, pcLoadings);
		testPC.setFrameName("reCalControlPC");
		testPC.dump( getTestPCScores());
	}

	public static void main(String[] args) {
		CalculatePCAsControlFromLoading rcpca = new CalculatePCAsControlFromLoading();
		rcpca.setPcaLoadings("GeneLoadings.VehCells.PCs.txt");
		rcpca.setControlLouvainIds("VehicleCells.Louvain.IDs.txt");
		rcpca.setControlPCScores("VehicleCells.PC.Scores.txt");
		rcpca.setTestPCScores("VismodegibCells.PC.Scores.txt");
		rcpca.setControlScaledExression("VehicleCells.Scaled.UMIs.txt");
		rcpca.setControlTSNEScores("VehicleCells.tSNE.Scores.txt");
		rcpca.setTestScaledExpression("VismodegibCells.Scaled.UMIs.txt");
		rcpca.run();
	}

}
