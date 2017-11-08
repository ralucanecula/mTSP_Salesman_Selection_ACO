package aco;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

/**
 * This is a data reader for reading the input data set that will be partitioned into k clusters
 */
public class DataReader {
    /** file name */
    private String fileName = "mtsp51.txt";
    
    /** number of attributes */
    private int nrAttributes = 2;
    
    /** separator */
    private String separator = "\t";

    public DataReader(String fileName) {
        this.fileName = fileName;
    }

    /**
     * Read the data from the given file
     * @return the data from the file
     */
    public Data read() {
    	String[] str = fileName.split("\\.");
    	
		//we are dealing with a input file in TSP format 
		if (str[1].equals("tsp")) {
			boolean foundCoordSection = false;
			separator = " ";
			File file = new File("input/tsplib/" + fileName);
			BufferedReader in = null;
	        ArrayList<Record> records = new ArrayList<Record>();
	        int lineNumber = 0;
	        Record depotCity = null;
	        Data d;
		        
	        if (file.exists()) {
	            try {
	                in = new BufferedReader(new FileReader(file));
	                String line = in.readLine();
	                while (line != null) {
	                	if (line.startsWith("EOF")) {
					    	break;
					    }
				
					    if (!foundCoordSection) {
							if (line.startsWith("NAME")) {
							} else if (line.startsWith("COMMENT")) {
							} else if (line.startsWith("TYPE") && !line.contains("TSP")) {
							    System.err.println("Not a TSP Tsp.instance in TSPLIB format !!");
							    System.exit(1);
							} else if (line.startsWith("DIMENSION")) {		
							} else if (line.startsWith("DISPLAY_DATA_TYPE")) {
							} else if (line.startsWith("EDGE_WEIGHT_TYPE")) {
							}
					    } else {
							if (lineNumber == 0) {
		                		depotCity = readDepotCity(line, separator);
		                	} 
							else {
		                		readRecord(records, line, separator);
		                	}
							lineNumber++;						
					    }
				
					    if (line.startsWith("NODE_COORD_SECTION")) {
					    	foundCoordSection = true;
					    }
				
					    line = in.readLine();
	                }
	                in.close();
	            } catch (FileNotFoundException ignored) {
	            } catch (IOException e) {
	                System.out.println("Error occurred while reading file: " + file + " " + e.getMessage());
	            }
	        } else {
	            return null;
	        }
	        d = new Data(records);
	        d.setDepotCity(depotCity);
	        return d;
        }
		            
		//we are dealing with a mTSP input file
		else  if (str[1].equals("txt")) {
			separator = "\t";
			File file = new File("input/" + fileName);
	        BufferedReader in = null;
	        ArrayList<Record> records = new ArrayList<Record>();
	        int lineNumber = 0;
	        Record depotCity = null;
	        Data d;
	        
	        if (file.exists()) {
	            try {
	                in = new BufferedReader(new FileReader(file));
	                String line = in.readLine();
	                while (line != null) {
	                	if (lineNumber == 0) {
	                		depotCity = readDepotCity(line, separator);
	                	}
	                	else {
	                		readRecord(records, line, separator);
	                	}
	                    line = in.readLine();
	                    lineNumber++;
	                }
	                in.close();
	            } catch (FileNotFoundException ignored) {
	            } catch (IOException e) {
	                System.out.println("Error occurred while reading file: " + file + " " + e.getMessage());
	            }
	        } else {
	            return null;
	        }
	
	        d = new Data(records);
	        d.setDepotCity(depotCity);
	        return d;
		}
        return null;
    }

    /**
     * Read the first line from the file representing the depot (home) city
     * @param record first record from the file
     * @param line string to extract the record from
     */
    private Record readDepotCity(String line, String inputFormat) {
        String[] strRecord = line.split(separator);
        ArrayList<String> strRecordFiltered = new ArrayList<String>();
        
        //some entries from strRecord vector may be empty space, so remove them from initial record vector
        //since they not contain any useful information
        for (int i = 0; i < strRecord.length; i++) {
        	if (strRecord[i].equals("")) {
        		continue;
        	}
        	else {
        		strRecordFiltered.add(strRecord[i]);
        	}
        }
        
        strRecord = strRecordFiltered.toArray(new String[0]);
        
        double[] recordValues = new double[nrAttributes];
        Record depotRecord = null;
        int id = 0;
        
        try {
            id = Integer.parseInt(strRecord[0].trim());
        } catch (NumberFormatException e) {
            System.out.println("NumberFormatException " + e.getMessage() + " strRecord[0]=" + strRecord[0] + " line=" + line);
        }
        
        for (int i = 0; i < nrAttributes; i++) {
            try {
                recordValues[i] = Double.parseDouble(strRecord[i + 1].trim());
            } catch (NumberFormatException e) {
                System.out.println("NumberFormatException " + e.getMessage() + " i= " + i
                		+ " strRecord[i]=" + strRecord[i] + " line=" + line);
            }
        }
        depotRecord = new Record(line, recordValues, id);
        return depotRecord;
    }
    
    /**
     * Read a single record and store it in records
     * @param records records list
     * @param line string to extract the record from
     */
    private void readRecord(ArrayList<Record> records, String line, String inputFormat) {
        String[] strRecord = line.split(separator);
        ArrayList<String> strRecordFiltered = new ArrayList<String>();
        
        //some entries from strRecord vector may be empty space, so remove them from initial record vector
        //since they not contain any useful information
        for (int i = 0; i < strRecord.length; i++) {
        	if (strRecord[i].equals("")) {
        		continue;
        	}
        	else {
        		strRecordFiltered.add(strRecord[i]);
        	}
        }
        
        strRecord = strRecordFiltered.toArray(new String[0]);
        
        double[] record = new double[nrAttributes];
        int id = 0;
        
        try {
            id = Integer.parseInt(strRecord[0].trim());
        } catch (NumberFormatException e) {
            System.out.println("NumberFormatException " + e.getMessage() + " strRecord[0]=" + strRecord[0] + " line=" + line);
        }
        
        for (int i = 0; i < nrAttributes; i++) {
            try {
                record[i] = Double.parseDouble(strRecord[i + 1].trim());
            } catch (NumberFormatException e) {
                System.out.println("NumberFormatException " + e.getMessage() + " i= " + i
                		+ " strRecord[i]=" + strRecord[i] + " line=" + line);
            }
        }
        Record r = new Record(line, record, id);
        records.add(r);
    }

	public String getSeparator() {
		return separator;
	}

	public void setSeparator(String separator) {
		this.separator = separator;
	}
    
    
    
}
