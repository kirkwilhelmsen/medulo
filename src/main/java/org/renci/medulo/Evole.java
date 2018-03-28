package org.renci.medulo;

public class Evole extends Abstract_FileSetting implements Runnable {
	String []args;
	public void run() {
		Evolution ev = new Evolution();
		ev.setScaledExpresion("VehicleCells.Scaled.UMIs.txt");
		ev.setGroupStart("cell.tsv");
		ev.setGenomes2Save(Integer.valueOf(args[0]));
		ev.setFractioChange(Double.valueOf(args[1]));
		ev.setProbMutation(Double.valueOf(args[2]));
		ev.setMinProgeney(Integer.valueOf(args[3]));
		ev.setReps(Integer.valueOf(args[4]));
		ev.setWt0(Double.valueOf(args[5]));
		ev.setWt1(Double.valueOf(args[6]));
		ev.setPareto(Boolean.valueOf(args[7]));
		ev.setAngle(Double.valueOf(args[8]));
		ev.setOutFile(args[9]);
		ev.run();
	}
	
	public static void main(String[] args) {
		Evole rcpca = new Evole();
		rcpca.args=args;
		rcpca.run();
	}

}
