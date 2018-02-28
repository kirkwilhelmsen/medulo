package org.renci.medulo.tools;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

import org.renci.medulo.Constants;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class CHATBufferedFileReader {

	private static final Logger logger = LoggerFactory.getLogger(CHATBufferedFileReader.class);
	
    public BufferedReader reader = null;

    public FileReader ff = null;
    

    public CHATBufferedFileReader(File f) {
        try {
            ff = new FileReader(f);
            reader = new BufferedReader(ff, Constants.BUFFER_SIZE_MEDIUM);
        } catch (FileNotFoundException e) {
            logger.error("{}, {} file could not be opened", e.getMessage(), f.getAbsolutePath(), e);
            System.exit(0);
        }
    }

    public void close() {
        try {
           	reader.close();
            reader = null;
           	ff.close();
           	ff = null;
        } catch (IOException e) {
        	logger.error(e.getMessage(), e);
            System.exit(0);
        }
    }

    public String nextLine() {
        try {
            return reader.readLine();
        } catch (IOException e) {
        	logger.error(e.getMessage(), e);
            System.exit(0);
        }
        return null;
    }

}