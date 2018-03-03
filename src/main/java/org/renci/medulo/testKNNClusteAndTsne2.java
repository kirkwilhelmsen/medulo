package org.renci.medulo;

import org.renci.medulo.tools.CHATBufferedFileWriter;

/*
I recommend running unstandardized
*/
public class testKNNClusteAndTsne2 extends Abstract_FileSetting implements Runnable {
	int cellCountToCheck = 10;
	
	@Override
	public void run() {
		BigFrame cluster = new BigFrame("cluster",getControlLouvainIds(), "control_cell", "cluster", "\t");
		cluster.load(false);
		BigFrame tsne = new BigFrame("tsne",getControlTSNEScores(), "control_cell", "tsne", "\t");
		tsne.load(false);
		BigFrame controlPC = new BigFrame("controlPC",getControlPCScores(), "PC", "control_cell", "\t");
		controlPC.load(false);
		kNearestNeighbor knn = new kNearestNeighbor(controlPC, 14, 10, tsne, cluster,false,true);
		CHATBufferedFileWriter out = new CHATBufferedFileWriter();
		out.open("testKNN_unstand_ambig.txt");
		out.writeString("tsne1\ttsne2\tx");
		for(int c=0;c<10;c++) {
			String currentCell = controlPC.getyAxisRowNames()[c];
			int tSNERow = tsne.getyAxisRowNames2Index().get(currentCell);
			int clusterRow = cluster.getyAxisRowNames2Index().get(currentCell);
			out.writeString(knn.calculate_KNN(controlPC.getData()[c],controlPC.getyAxisRowNames()[c]) + "\t" + tsne.getData()[tSNERow][0].toString()+ "\t" + tsne.getData()[tSNERow][1].toString() + "\t" + (Math.round(cluster.getData()[clusterRow][0])));
		}
		out.close();
		
	}

	public static void main(String[] args) {
		testKNNClusteAndTsne2 rcpca = new testKNNClusteAndTsne2();
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
