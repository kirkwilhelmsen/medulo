package org.renci.medulo;

import org.renci.medulo.tools.CHATBufferedFileWriter;

public class CalculateKNNClusteAndTsne extends Abstract_FileSetting implements Runnable {
	int cellCountToCheck = 10;
	
	@Override
	public void run() {
		BigFrame cluster = new BigFrame("cluster",getControlLouvainIds(), "control_cell", "cluster", "\t");
		cluster.load(true);
		BigFrame tsne = new BigFrame("tsne",getControlTSNEScores(), "control_cell", "tsne", "\t");
		tsne.load(true);
		BigFrame testPC = new BigFrame("testPC",getTestPCScores(), "PC", "test_cell", "\t");
		testPC.load(false);
		BigFrame controlPC = new BigFrame("controlPC",getControlPCScores(), "PC", "control_cell", "\t");
		controlPC.load(false);
		kNearestNeighbor knn = new kNearestNeighbor(controlPC, 14, 3, tsne, cluster, true,true);
		CHATBufferedFileWriter out = new CHATBufferedFileWriter();
		out.open("controlCells.tSNE.cluster.Std_notamb_3knn.txt");
		logger.info("Opening knn output file");
		out.writeString("tsne1\ttsne2\tx");
		for(int c=0;c<testPC.getData().length;c++) {
			out.writeString(knn.calculate_KNN(controlPC.getData()[c],controlPC.getyAxisRowNames()[c]));
		}
		out.close();
		logger.info("Closing knn output file");
		
	}

	public static void main(String[] args) {
		CalculateKNNClusteAndTsne rcpca = new CalculateKNNClusteAndTsne();
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
