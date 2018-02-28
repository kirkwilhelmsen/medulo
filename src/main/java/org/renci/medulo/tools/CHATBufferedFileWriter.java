package org.renci.medulo.tools;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class CHATBufferedFileWriter {

	private static final Logger logger = LoggerFactory.getLogger(CHATBufferedFileWriter.class);
			
    public BufferedWriter bufferedWriter = null;

    private String filename = null;

    public void open(String pathAndFileName) {
        filename = pathAndFileName;
        this.open(false);
    }
    
    public void openWithAppend(String pathAndFileName) {
        filename = pathAndFileName;
        this.open(true);
    }

    private void open(boolean append) {
        try {
            File file = new File(filename);
            FileWriter fileWriter = new FileWriter(file, append);
            bufferedWriter = new BufferedWriter(fileWriter);
        } catch (IOException e) {
        	logger.error("{}, {} file could not be opened", e.getMessage(), filename, e);
            System.exit(-1);
        }
    }

    public void flush() {
        try {
            bufferedWriter.flush();
        } catch (IOException e) {
            logger.error(e.getMessage(), e);
            System.exit(-1);
        }
    }

    public void close() {
        try {
            bufferedWriter.close();
        } catch (IOException e) {
            logger.error(e.getMessage(), e);
        }
    }

    public void writeString(String lineToWrite) {
        try {
            bufferedWriter.append(lineToWrite);
            bufferedWriter.newLine();
        } catch (IOException e) {
            logger.error(e.getMessage(), e);
            System.exit(-1);
        } catch (Exception e) {
            logger.error(e.getMessage(), e);
            System.exit(-1);
        }
    }

	public void writeString(String writeString, boolean echo) {
		writeString(writeString);
		
		if(echo)
			logger.info("writeString echo {}", writeString);
		
	}

	public String getFilename() {
		return filename;
	}

	public void setFilename(String filename) {
		this.filename = filename;
	}

}