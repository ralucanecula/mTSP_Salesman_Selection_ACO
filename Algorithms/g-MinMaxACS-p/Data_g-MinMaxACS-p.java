package aco;

import java.util.ArrayList;

/**
 * Represent the read data.
 */
public class Data {
    /** Data records */
    private ArrayList<Record> records;
    
    /** the depot city (home city) from the mTSP problem such that all the salespersons have to start 
     * and end their tour at depot (home city). */
    public Record depotCity;

    public Data(ArrayList<Record> records) {
        this.records = records;
    }

    public ArrayList<Record> getRecords() {
        return records;
    }

	public Record getDepotCity() {
		return depotCity;
	}

	public void setDepotCity(Record depotCity) {
		this.depotCity = depotCity;
	}
    
    
   

}
