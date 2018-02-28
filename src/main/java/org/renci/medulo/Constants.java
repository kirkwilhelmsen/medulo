package org.renci.medulo;


public interface Constants {
	
	

    public static final Integer BUFFER_SIZE_JUMBO = Double.valueOf(Math.pow(2, 20)).intValue();
    public static final Integer BUFFER_SIZE_LARGE = Double.valueOf(Math.pow(2, 16)).intValue();
    public static final Integer BUFFER_SIZE_MEDIUM = Double.valueOf(Math.pow(2, 12)).intValue();
    public static final Integer BUFFER_SIZE_SMALL = Double.valueOf(Math.pow(2, 8)).intValue();
    public static final Integer BUFFER_SIZE_TINY = Double.valueOf(Math.pow(2, 4)).intValue();

    public static final String TAB = "\t";

    public static final String NEW_LINE = "\n";

    public static final String COMMA = ",";
    
    public static final String DASH = "-";

    public static final char DOUBLE_QUOTE = '"';
    
}

