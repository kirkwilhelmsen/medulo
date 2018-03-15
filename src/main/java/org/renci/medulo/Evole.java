package org.renci.medulo;

import java.io.File;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.Map;
import java.util.Set;

import org.renci.medulo.tools.CHATBufferedFileReader;

public class Evole extends Abstract_FileSetting implements Runnable {
	File groupStart;
	
	@Override
	public void run() {
		CHATBufferedFileReader in = new CHATBufferedFileReader(groupStart);
		String l;
		Set<String> allGenes = new LinkedHashSet<String>();
		Map<String, Set<String>> gg= new LinkedHashMap<String, Set<String>>();
		while((l=in.nextLine())!=null) {
			String [] f = l.split("\t");
			if(!gg.containsKey(f[1]))gg.put(f[1], new LinkedHashSet<String>());
			gg.get(f[1]).add(f[0]);
			allGenes.add(f[0]);
		}
		in.close();
		/*int ct=0;
		for(String gr:gg.keySet())ct += gg.get(gr).size();
		if(allGenes.size()!=ct) {
			logger.error("Genes are on more than one list");
			System.exit(-1);
		}*/
		BigFrame controlExp = new BigFrame("controlExp",getControlScaledExression(), "control_cell", "gene", "\t");
		controlExp.load(false,allGenes);
		Evolution ev = new Evolution(gg, controlExp/*,new CompairEvolutionModel(10, 1), new CompairEvolutionModelAngle()*/);
		ev.setGenomes2Save(100);
		ev.setFractioChange(.5d);
		ev.setProbMutation(1);
		ev.setScoreThreshold(.3);
		ev.setMinProgeney(100);
		ev.setWt0(1);
		ev.setWt1(.1);
		ev.setPareto(true);
		ev.setAngle(5);
		ev.makeInitalModel();
		ev.makeModels(1000);
		ev.selectNextGen();
		ev.evolve(10,10);
		ev.write(new File("curResults,txt"));
		
	}

	public static void main(String[] args) {
		Evole rcpca = new Evole();
		rcpca.setControlScaledExression("VehicleCells.Scaled.UMIs.txt");
		rcpca.setTestScaledExpression("VismodegibCells.Scaled.UMIs.txt");
		rcpca.groupStart = new File("cell.tsv");
		rcpca.run();
	}

}
