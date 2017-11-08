package aco;

import java.util.ArrayList;

import aco.Tsp.Problem;

/*class that hold multiple instances of TSP problem and fill them with data, according the formed clusters
  given by the kMeans clustering algorithm 
*/
public class MTsp {
	
	//instance of mTSP problem
	public static Problem instance;
	
	//depot (start city)
	public static Tsp.Point depotCity;
	
	//contains the id of cities assigned to each cluster 
	public static ArrayList<Integer>[] clusterCities;
	
	//number of cities from the mTSP problem (e.g. for mtsp51.txt input file, there will be 51 cities)
	//the number of cities should exclude the depot city, so there will be 50 cities
	public static int n;
	
	//number of salesmen from the mTSP problem (it should be equal with the number of obtained clusters)
	public static int m;
	
	public MTsp() {}
	
	public MTsp(int nrCities, int nrSalesmen) {
		n = nrCities;
		m = nrSalesmen;
	}
	
	public MTsp(ArrayList<Record> records, Record depotCity_, int nrCities, int nrSalesmen) {
		this(nrCities, nrSalesmen);
		
		instance = new Tsp.Problem();
		instance.n = nrCities;
		instance.name = "mTSP with " + nrCities + " cities and " + nrSalesmen  + " salesmen"; 
		instance.nodes = new Tsp.Point[records.size()];
		
		//set coordinates for depot city
    	depotCity =  new Tsp.Point();
        double[] coords = depotCity_.getData();
    	depotCity.x = coords[0];
    	depotCity.y = coords[1];
    	double coordsCity[];
		
		for (int i = 0;  i < records.size(); i++) {
			instance.nodes[i] = new Tsp.Point();
			Record r = records.get(i);
        	coordsCity = r.getData();	
    		instance.nodes[i].x = coordsCity[0];
    		instance.nodes[i].y = coordsCity[1];	
		}
			
	}
    
    public void printDepotCity() {
    	 System.out.print("Depot (home) city coordinates" + ": (" + depotCity.x + ", " + depotCity.y + ")");	
    	 System.out.println();
    }
    
    public void printRecords(ArrayList<Record> records) {
    	System.out.println("Records are:");
    	for (int i = 0;  i < records.size(); i++) {
    		System.out.println(records.get(i).getActual());
    	}
    }
	
}
