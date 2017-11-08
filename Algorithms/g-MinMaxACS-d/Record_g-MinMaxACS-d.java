package aco;

/**
 * A single data records read from the file.
 */
public class Record {
    /** Actual record as in the file */
    private String actual;
    
    /** Double array as data */
    private double[] data;
    
    /** id of the record */
    private int id;
    
    public Record(String actual, double[] data) {
        this.actual = actual;
        this.data = data;
    }
    
    public Record(String actual, double[] data, int id) {
        this(actual, data);
        this.id = id;
    }

    /**
     * Get the actual string read from the file
     * @return the actual record as in the file
     */
    public String getActual() {
        return actual;
    }

    /**
     * Get the data as a double array
     * @return data as a double array
     */
    public double[] getData() {
        return data;
    }

	public int getId() {
		return id;
	}

	public void setId(int id) {
		this.id = id;
	}
    
    
   
    
}
