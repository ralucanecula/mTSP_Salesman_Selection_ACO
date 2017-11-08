package aco;

import java.util.ArrayList;

import aco.Ants.Ant;

/**
 * ACO algorithms for the TSP
 * 
 * This code is based on the ACOTSP project of Thomas Stuetzle.
 * It was initially ported from C to Java by Adrian Wilke.
 * 
 * Project website: http://adibaba.github.io/ACOTSPJava/
 * Source code: https://github.com/adibaba/ACOTSPJava/
 */

/* this Eclipse project represents the ACO algorithm adapted from the mTSP_ACO global Eclipse project,
 * in which the quantity of the laid pheromone depends on the total cost of the tours and is the same 
 * on all the visited edges. Also there is no clustering information used (k and l lower bounds or 
 * probabilities given by the fuzzy c-Means clustering algorithm)  
 * In this version the laid pheromone depends on the lowest cost of the longest tour -> min(max(subtours))
 * for each of the n_ants (usually equal to 10) candidate solutions, we keep the length of the longest tour
 * and best so far ant (global best) is chosen among the obtained solutions in an iteration such that the chosen
 * one will have the lowest longest cost; this solution will be used for the global pheromone update;
 * the quantity of the pheromone will be according to the total cost and will be computed in the same 
 * way on all the edges  
 * 
 * 
 * for solving the mTSP problem and it's an improved one that actually add the depot city in the tour that each salesman travel,
 * and pheromone is deposited on the edges incident (that enters or leave) with the depot city 
 * this approach should be a better one, since in nature this is how things actually happen: the ant 
 * leaves/deposit pheromone anywhere it travels in the search for food
 * 
 * the only difference with respect to this second approach (see FuzzyClustering_mTSP_ACO project in Eclipse)
 * is that in the transition formula for an ant it doesn't uses the probabilities given by fuzzy c-Means
 */

/* contains the main entry point for the implementation of ACO for solving the mTSP problem
 */
public class TSP_ACO {
	
	/*
     * ################################################
     * ########## ACO algorithms for the TSP ##########
     * ################################################
     * 
     * Version: 1.0
     * File: main.c
     * Author: Thomas Stuetzle
     * Purpose: main routines and control for the ACO algorithms
     * Check: README and gpl.txt
     * Copyright (C) 2002 Thomas Stuetzle
     */

    /***************************************************************************
     * Program's name: acotsp
     * 
     * Ant Colony Optimization algorithms (AS, ACS, EAS, RAS, MMAS, BWAS) for the
     * symmetric TSP
     * 
     * Copyright (C) 2004 Thomas Stuetzle
     * 
     * This program is free software; you can redistribute it and/or modify
     * it under the terms of the GNU General Public License as published by
     * the Free Software Foundation; either version 2 of the License, or
     * (at your option) any later version.
     * 
     * This program is distributed in the hope that it will be useful,
     * but WITHOUT ANY WARRANTY; without even the implied warranty of
     * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
     * GNU General Public License for more details.
     * 
     * You should have received a copy of the GNU General Public License
     * along with this program; if not, write to the Free Software
     * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
     * 
     * email: stuetzle no@spam informatik.tu-darmstadt.de
     * mail address: Universitaet Darmstadt
     * Fachbereich Informatik
     * Hochschulstr. 10
     * D-64283 Darmstadt
     * Germany
     ***************************************************************************/
	
	/** File name to be used for input data set */
	//private static String fileName = "rat99.tsp";
	private static String fileName = "mtsp100.txt";
	
	/** the number of salesmen from the mTSP instance)*/
    private static int m = 10;
    
    /** name of the TSP input instance from TSPLIB **/
    public static String instanceName = new String();
    
    /** indicates whether local search is used in the ACO algorithm **/
    public static boolean ls_flag = true; 


	//checks whether termination condition is met
    static boolean termination_condition() {   
    	//return (((InOut.n_tours >= InOut.max_tours) && (Timer.elapsed_time() >= InOut.max_time)) || (Ants.best_so_far_ant.tour_length <= InOut.optimal));
    	//return ((InOut.n_tours >= InOut.max_tours) && (Timer.elapsed_time() >= InOut.max_time));
    	return (InOut.iteration >= InOut.max_iterations);
    }
    
    /*static int chooseSalesman(double[] probabilities) {
    	int i;
    	double rnd, partialSum = 0.;
    	
    	rnd = Utilities.random01();
    	i = 0;
    	partialSum = probabilities[i];
    	
    	 This loop always stops because last value from probabilities vector is a huge value 
    	while (partialSum <= rnd) {
			i++;
			partialSum += probabilities[i];
	    }
    	return i;
    }*/
    
    //check if there is still an ant with left cities to visit
    static boolean isDone() {
    	boolean done = true;
    	
    	for (int k = 0; k < Ants.n_ants; k++) {
    		if (Ants.ants[k].toVisit > 0) {
    			return false;
    		}
    	}
    	
    	return done;
    }
    
    //choose in a probabilistic way (according to the length of its subtour) the next salesman 
    //that will visit a city; when all the subtours are empty choose in a random manner the salesman
    static int chooseSalesmanProbabilistic(Ant a, double[] fitnessValues, double totalFitnes) {
    	//test if all the ants'subtours are empty
    	boolean isEmpty = true;
    	int nextSalesman;
    	double rnd, partialSum = 0.;
    	double[] probabilities = new double[MTsp.m + 1];
		int j;
	
    	for (int i = 0; i < MTsp.m; i++) {
    		if (a.tour_lengths[i] > 0) {
    			isEmpty = false;
    			break;
    		}
    	}
    	//when all the subtours are empty and no city is visited choose randomly the salesman
    	if (isEmpty) {
    		nextSalesman = (int)(Math.random() * MTsp.m);
    	}
    	else {
    		for (int i = 0; i < MTsp.m; i++) { 
    			probabilities[i] = fitnessValues[i]/totalFitnes;
    		}
    		probabilities[MTsp.m] = Double.POSITIVE_INFINITY;	
        	
        	rnd = Utilities.random01();
        	j = 0;
        	partialSum = probabilities[j];
        	
        	/* This loop always stops because last value from probabilities vector is a huge value */
        	while (partialSum <= rnd) {
    			j++;
    			partialSum += probabilities[j];
    	    }
        	nextSalesman = j;
    	}
    	return nextSalesman;
    }

    //manage the solution construction phase (construct a set of complete and closed tours, 
    //each tour for one of the m salesmen)
    static void construct_solutions() {
		int k; /* counter variable */
		int step; /* counter of the number of construction steps */
		int salesman;
		double[] fitnessValues = new double[MTsp.m];
		double totalFitness = 0.0;
		
		/* Mark all cities as unvisited */
		for (k = 0; k < Ants.n_ants; k++) {
		    Ants.ant_empty_memory(Ants.ants[k]);
		}
	
		step = 0;
		/* Place the ants on same initial city, which is the depot city */
		for (k = 0; k < Ants.n_ants; k++) {
			//Ants.place_ant(Ants.ants[k], step);
			for (int i = 0; i < MTsp.m; i++) {
				//place each ant on the depot city, which is the start city of each tour
				// -1 is a special marker for the deport city, so that it's not be confused with the rest of the cities
				// all the rest of the cities are represented by integer values > 0
				Ants.ants[k].tours[i].add(-1);   
			}
		}
	
		while (!isDone()) {
		    for (k = 0; k < Ants.n_ants; k++) {
		    	totalFitness = 0;
		    	if (Ants.ants[k].toVisit > 0) {
		    		//choose for each ant in a probabilistic way the salesman that will visit the next city
		    		//the salesman is selected in a probabilistic way according to some type of roulette wheel selection: 
		    		//in which a high probability corresponds to a subtour of lowest length
		    		
		    		//compute lengths for all partial subtours constructed so far 
			    	for (int i1 = 0; i1 < MTsp.m; i1++) {
			    		Ants.ants[k].tour_lengths[i1] = Tsp.compute_tour_length(Ants.ants[k].tours[i1]);
			    		fitnessValues[i1] = 1.0 / (Ants.ants[k].tour_lengths[i1] + 0.1); 
			    		totalFitness += fitnessValues[i1];
			    	}
		    		salesman = chooseSalesmanProbabilistic(Ants.ants[k], fitnessValues, totalFitness); 
			    	
					Ants.neighbour_choose_and_move_to_next(Ants.ants[k], salesman);
					if (Ants.acs_flag)
					    Ants.local_acs_pheromone_update(Ants.ants[k], salesman);
		    	}
		    	
		    }
		}
	
		double longestTourLength;
		int idLongestTour = 0;
		for (k = 0; k < Ants.n_ants; k++) {
			longestTourLength = Double.MIN_VALUE;
			for (int i = 0; i < MTsp.m; i++) {
				step = Ants.ants[k].tours[i].size();
				Ants.ants[k].tours[i].add(step, -1);
				
				Ants.ants[k].tour_lengths[i] = Tsp.compute_tour_length_(Ants.ants[k].tours[i]);
				Ants.ants[k].total_tour_length += Ants.ants[k].tour_lengths[i];
				if (longestTourLength < Ants.ants[k].tour_lengths[i]) {
					longestTourLength = Ants.ants[k].tour_lengths[i];
					idLongestTour = i;
				}
				
			    if (Ants.acs_flag)
			    	Ants.local_acs_pheromone_update(Ants.ants[k], i);
			}
			Ants.ants[k].longest_tour_length = longestTourLength;
			Ants.ants[k].indexLongestTour = idLongestTour;
			//Ants.ants[k].total_tour_length = Tsp.compute_tour_lengths(Ants.ants[k].tours);
			Ants.ants[k].costObjectives[0] = Ants.ants[k].total_tour_length;
			Ants.ants[k].costObjectives[1] = Ants.computeToursAmplitude(Ants.ants[k]);
		}
		InOut.n_tours += (Ants.n_ants * MTsp.m); //each ant constructs a complete and closed tour
    }

    //initialize variables appropriately when starting a trial
    static void init_try() {

		Timer.start_timers();
		InOut.time_used = Timer.elapsed_time();
		InOut.time_passed = InOut.time_used;
	
		/* Initialize variables concerning statistics etc. */
		InOut.n_tours = 1 * MTsp.m;
		InOut.iteration = 1;
		Ants.best_so_far_ant.total_tour_length = Integer.MAX_VALUE;
		InOut.found_best = 0;
		
		/*
		 * Initialize the Pheromone trails, only if ACS is used, Ants.pheromones
		 * have to be initialized differently
		 */
		if (!(Ants.acs_flag)) {
		    Ants.trail_0 = 1. / ((Ants.rho) * Ants.nn_tour());
		    /*
		     * in the original papers on Ant System it is not exactly defined what the
		     * initial value of the Ants.pheromones is. Here we set it to some
		     * small constant, analogously as done in MAX-MIN Ant System.
		     */
		    Ants.init_pheromone_trails(Ants.trail_0);
		}
		if (Ants.acs_flag) {
		    Ants.trail_0 = 1. / ((double) MTsp.n * (double) Ants.nn_tour());
		    Ants.init_pheromone_trails(Ants.trail_0);
		}
	
		/* Calculate combined information Ants.pheromone times heuristic information */
		Ants.compute_total_information();
		
    }
    
    //apply to each solution constructed by a ant, a local search phase to further improve the solution
    //the local search procedure is composed of 3 operators applied sequentially in this order: relocation, exchange and 2-opt
    //as a result all ants of the colony will have locally optimal tours
    static void apply_local_search() {
    	for (int k = 0; k < Ants.n_ants; k++) {
    		Ants.ants[k] = local_search(Ants.ants[k]);
    	}
    }
    
    static Ant local_search(Ant a) {  	
		//apply relocation local search operator 
		//a = relocation(a);  
		
		//apply exchange local search operator 
		a = exchangeIterated(a);    
		
		//apply 2-opt heuristic on each of the m routes; 2-opt tries to improve single route
		//by replacing its two non-adjacent edges by two other edges
		for (int l = 0; l < m; l++) {
			LocalSearch.two_opt_first(a.tours[l]);
		}
		
		//compute new distances and update longest tour
		double sum = 0.0;
		double longestTourLength = Double.MIN_VALUE;
		int idLongestTour = 0;
		for (int l = 0; l < m; l++) {
			a.tour_lengths[l] = Tsp.compute_tour_length(a.tours[l]);
			sum += a.tour_lengths[l];
			if (longestTourLength < a.tour_lengths[l]) {
				longestTourLength = a.tour_lengths[l];
				idLongestTour = l;
			}
		}
		a.total_tour_length = sum;
		a.longest_tour_length = longestTourLength;
		a.indexLongestTour = idLongestTour;
    	
		return a;
    }
    
    static Ant exchange(Ant a) {
    	int k, city1, city2, idLongestTour;
    	double longestTourLength;
    	int x1, x2;
    	
    	boolean found = false;
    	for (int index = 0; index < m; index++) {   
    		int lastPos =  a.tours[index].size() - 1;
			x1 = a.tours[index].get(0);
			x2 = a.tours[index].get(lastPos);
			if (x1 != -1 || x2 != -1) {
				found = true;
				break;
			}
    	}
    	if (found) {
    		System.out.println("Found not standard subtour..");
    	}
    	
    	Ant improvedAnt = new Ant();
    	improvedAnt.tours = (ArrayList<Integer>[])new ArrayList[MTsp.m];
    	improvedAnt.tour_lengths = new double[MTsp.m];
	    for (int j = 0; j < MTsp.m; j++) {
	    	improvedAnt.tours[j] = new ArrayList<Integer>();
	    	improvedAnt.tour_lengths[j] = 0;
	    }
	    improvedAnt.visited = new boolean[MTsp.n];
	    improvedAnt.toVisit = MTsp.n;
	    improvedAnt.longest_tour_length = Double.MAX_VALUE;
	    
	    improvedAnt.costObjectives = new double[2];
	    for (int indexObj = 0; indexObj < 2; indexObj++) {
	    	improvedAnt.costObjectives[indexObj] = 0;
    	}
	    
	    Ant temp = new Ant();
	    temp.tours = (ArrayList<Integer>[])new ArrayList[MTsp.m];
	    temp.tour_lengths = new double[MTsp.m];
	    for (int j = 0; j < MTsp.m; j++) {
	    	temp.tours[j] = new ArrayList<Integer>();
	    	temp.tour_lengths[j] = 0;
	    }
	    temp.visited = new boolean[MTsp.n];
	    temp.toVisit = MTsp.n;
	    temp.longest_tour_length = Double.MAX_VALUE;
		
	    temp.costObjectives = new double[2];
	    for (int indexObj = 0; indexObj < 2; indexObj++) {
	    	temp.costObjectives[indexObj] = 0;
    	}
    	
    	//Ant iterBestAnt = Ants.ants[indexAnt];
    	Ants.copy_from_to(a, improvedAnt);
    	Ants.copy_from_to(a, temp);
    	
    	k = a.indexLongestTour;
    	//iterate over cities from the longest tour except the first and last city from the tour which is the depot
    	for (int i = 1; i < a.tours[k].size() - 1; i++) {   
    		city1 = a.tours[k].get(i);
    		//iterate over the rest of other tours
    		for (int l = 0; l < m; l++) {
    			if (l != k) {
    				for (int j = 1; j < a.tours[l].size() - 1; j++) {
	    				city2 = a.tours[l].get(j);
	    				//swap city1 and city2 in the 2 considered subtours
	    				temp.tours[k].set(i, city2);
	    				temp.tours[l].set(j, city1);
	    				
	    				//compute the new subtour lengths and determine longest subtour
	    				temp.tour_lengths[k] = Tsp.compute_tour_length(temp.tours[k]);
	    				temp.tour_lengths[l] = Tsp.compute_tour_length(temp.tours[l]);
	    				
	    				longestTourLength = Double.MIN_VALUE;
	    				idLongestTour = 0;
	    				for (int index = 0; index < MTsp.m; index++) {
	    					if (longestTourLength < temp.tour_lengths[index]) {
	    						longestTourLength = temp.tour_lengths[index];
	    						idLongestTour = index;
	    					}
	    				}
	    				//verify if we obtain an improvement after this exchange operation
	    				if (longestTourLength < improvedAnt.longest_tour_length) {
	    					//System.out.println("Improvement found..");
	    					temp.longest_tour_length = longestTourLength;
	    					temp.indexLongestTour = idLongestTour;
	    					Ants.copy_from_to(temp, improvedAnt);
	    				}
	    				
	    				//restore temp to its initial value
	    				Ants.copy_from_to(a, temp);
    				}
    			}
    			
    		}
    	}
    	
    	return improvedAnt;
    }
    
    static Ant exchangeIterated(Ant a) {
    	int k, city1, city2, idLongestTour;
    	double longestTourLength;
    	boolean foundImprovement = true;
    	int startValue, endValue;
        ArrayList<Integer> result;
     	
     	for (int index = 0; index < m; index++) {   
     		int lastPos =  a.tours[index].size() - 1;
     		startValue = a.tours[index].get(0);
     		endValue = a.tours[index].get(lastPos);
 			if (startValue != -1 || endValue != -1) {
 				result = new ArrayList<Integer>(lastPos + 1);
 				result = rotateTour(a.tours[index]);	
 				a.tours[index] = new ArrayList<Integer>(result);
 			}
     	}
    	
    	Ant improvedAnt = new Ant();
    	improvedAnt.tours = (ArrayList<Integer>[])new ArrayList[MTsp.m];
    	improvedAnt.tour_lengths = new double[MTsp.m];
	    for (int j = 0; j < MTsp.m; j++) {
	    	improvedAnt.tours[j] = new ArrayList<Integer>();
	    	improvedAnt.tour_lengths[j] = 0;
	    }
	    improvedAnt.visited = new boolean[MTsp.n];
	    improvedAnt.toVisit = MTsp.n;
	    improvedAnt.longest_tour_length = Double.MAX_VALUE;
	    
	    improvedAnt.costObjectives = new double[2];
	    for (int indexObj = 0; indexObj < 2; indexObj++) {
	    	improvedAnt.costObjectives[indexObj] = 0;
    	}
	    
	    Ant temp = new Ant();
	    temp.tours = (ArrayList<Integer>[])new ArrayList[MTsp.m];
	    temp.tour_lengths = new double[MTsp.m];
	    for (int j = 0; j < MTsp.m; j++) {
	    	temp.tours[j] = new ArrayList<Integer>();
	    	temp.tour_lengths[j] = 0;
	    }
	    temp.visited = new boolean[MTsp.n];
	    temp.toVisit = MTsp.n;
	    temp.longest_tour_length = Double.MAX_VALUE;
		
	    temp.costObjectives = new double[2];
	    for (int indexObj = 0; indexObj < 2; indexObj++) {
	    	temp.costObjectives[indexObj] = 0;
    	}
    	
    	Ants.copy_from_to(a, improvedAnt);
    	Ants.copy_from_to(a, temp);
    	
    	while (foundImprovement) {
    		foundImprovement = false;
    		
    		Ants.copy_from_to(improvedAnt, a);
    		Ants.copy_from_to(a, temp);
    		
    		k = a.indexLongestTour;
	    	//iterate over cities from the longest tour except the first and last city from the tour which is the depot
	    	for (int i = 1; i < a.tours[k].size() - 1; i++) {   
	    		city1 = a.tours[k].get(i);
	    		//iterate over the rest of other tours
	    		for (int l = 0; l < m; l++) {
	    			if (l != k) {
	    				for (int j = 1; j < a.tours[l].size() - 1; j++) {
		    				city2 = a.tours[l].get(j);
		    				//swap city1 and city2 in the 2 considered subtours
		    				temp.tours[k].set(i, city2);
		    				temp.tours[l].set(j, city1);
		    				
		    				//compute the new subtour lengths and determine longest subtour
		    				temp.tour_lengths[k] = Tsp.compute_tour_length(temp.tours[k]);
		    				temp.tour_lengths[l] = Tsp.compute_tour_length(temp.tours[l]);
		    				
		    				longestTourLength = Double.MIN_VALUE;
		    				idLongestTour = 0;
		    				for (int index = 0; index < MTsp.m; index++) {
		    					if (longestTourLength < temp.tour_lengths[index]) {
		    						longestTourLength = temp.tour_lengths[index];
		    						idLongestTour = index;
		    					}
		    				}
		    				//verify if we obtain an improvement after this exchange operation
		    				if (longestTourLength < improvedAnt.longest_tour_length) {
		    					//System.out.println("Improvement found..");
		    					temp.longest_tour_length = longestTourLength;
		    					temp.indexLongestTour = idLongestTour;
		    					Ants.copy_from_to(temp, improvedAnt);
		    					foundImprovement = true;
		    				}
		    				
		    				//restore temp to its initial value
		    				Ants.copy_from_to(a, temp);
	    				}
	    			}
	    			
	    		}
	    	}
    	}
    	
    	
    	return improvedAnt;
    }
    
    static Ant relocation(Ant a) {
    	int k, city1, idLongestTour;
    	double longestTourLength;
        int x1, x2;
    	
    	boolean found = false;
    	for (int index = 0; index < m; index++) {   
    		int lastPos =  a.tours[index].size() - 1;
			x1 = a.tours[index].get(0);
			x2 = a.tours[index].get(lastPos);
			if (x1 != -1 || x2 != -1) {
				found = true;
				break;
			}
    	}
    	if (found) {
    		System.out.println("Found not standard subtour..");
    	}
    	
    	Ant improvedAnt = new Ant();
    	improvedAnt.tours = (ArrayList<Integer>[])new ArrayList[MTsp.m];
    	improvedAnt.tour_lengths = new double[MTsp.m];
	    for (int j = 0; j < MTsp.m; j++) {
	    	improvedAnt.tours[j] = new ArrayList<Integer>();
	    	improvedAnt.tour_lengths[j] = 0;
	    }
	    improvedAnt.visited = new boolean[MTsp.n];
	    improvedAnt.toVisit = MTsp.n;
	    improvedAnt.longest_tour_length = Double.MAX_VALUE;
	    
	    improvedAnt.costObjectives = new double[2];
	    for (int indexObj = 0; indexObj < 2; indexObj++) {
	    	improvedAnt.costObjectives[indexObj] = 0;
    	}
	    
	    Ant temp = new Ant();
	    temp.tours = (ArrayList<Integer>[])new ArrayList[MTsp.m];
	    temp.tour_lengths = new double[MTsp.m];
	    for (int j = 0; j < MTsp.m; j++) {
	    	temp.tours[j] = new ArrayList<Integer>();
	    	temp.tour_lengths[j] = 0;
	    }
	    temp.visited = new boolean[MTsp.n];
	    temp.toVisit = MTsp.n;
	    temp.longest_tour_length = Double.MAX_VALUE;
		
	    temp.costObjectives = new double[2];
	    for (int indexObj = 0; indexObj < 2; indexObj++) {
	    	temp.costObjectives[indexObj] = 0;
    	}
    	
    	Ants.copy_from_to(a, improvedAnt);
    	Ants.copy_from_to(a, temp);
    	
    	k = a.indexLongestTour;
    	//iterate over cities from the longest tour except the first and last city from the tour which is the depot
    	for (int i = 1; i < a.tours[k].size() - 1; i++) {   
    		city1 = a.tours[k].get(i);
    		temp.tours[k].remove(i);
    		
    		//iterate over the rest of other tours
    		for (int l = 0; l < m; l++) {
    			if (l != k) {
    				for (int j = 1; j < a.tours[l].size(); j++) {
	    				temp.tours[l].add(j, city1);
	    				
	    				//compute the new subtour lengths and determine longest subtour
	    				temp.tour_lengths[k] = Tsp.compute_tour_length(temp.tours[k]);
	    				temp.tour_lengths[l] = Tsp.compute_tour_length(temp.tours[l]);
	    				
	    				longestTourLength = Double.MIN_VALUE;
	    				idLongestTour = 0;
	    				for (int index = 0; index < MTsp.m; index++) {
	    					if (longestTourLength < temp.tour_lengths[index]) {
	    						longestTourLength = temp.tour_lengths[index];
	    						idLongestTour = index;
	    					}
	    				}
	    				//verify if we obtain an improvement after this exchange operation
	    				if (longestTourLength < improvedAnt.longest_tour_length) {
	    					//System.out.println("Improvement found..");
	    					temp.longest_tour_length = longestTourLength;
	    					temp.indexLongestTour = idLongestTour;
	    					Ants.copy_from_to(temp, improvedAnt);
	    				}
	    				
	    				//restore temp 
	    				//Ants.copy_from_to(iterBestAnt, temp);
	    				temp.tours[l].remove(j);
    				}
    			}
    			
			}
    		Ants.copy_from_to(a, temp);
    	}

    	return improvedAnt;
    }
    
    static ArrayList<Integer> rotateTour(ArrayList<Integer> t) {
    	ArrayList<Integer> result = new ArrayList<Integer> (t.size());
    	int pos = 0;
    	
    	for (int i = 0; i < t.size(); i++) {
    		if (t.get(i) == -1) {
    			pos = i;
    			break;
    		}
    	}
    	for (int i = pos; i < t.size(); i++) {
    		result.add(t.get(i));
    	}
    	for (int i = 1; i <= pos; i++) {
    		result.add(t.get(i));
    	}
    	
    	return result;
    }
    
    static Ant relocationIterated(Ant a) {
    	int k, city1, idLongestTour;
    	double longestTourLength;
    	boolean foundImprovement = true;
        int startValue, endValue;
        ArrayList<Integer> result;
    	
    	for (int index = 0; index < m; index++) {   
    		int lastPos =  a.tours[index].size() - 1;
    		startValue = a.tours[index].get(0);
    		endValue = a.tours[index].get(lastPos);
			if (startValue != -1 || endValue != -1) {
				result = new ArrayList<Integer>(lastPos + 1);
				result = rotateTour(a.tours[index]);	
				a.tours[index] = new ArrayList<Integer>(result);
			}
    	}
    	
    	Ant improvedAnt = new Ant();
    	improvedAnt.tours = (ArrayList<Integer>[])new ArrayList[MTsp.m];
    	improvedAnt.tour_lengths = new double[MTsp.m];
	    for (int j = 0; j < MTsp.m; j++) {
	    	improvedAnt.tours[j] = new ArrayList<Integer>();
	    	improvedAnt.tour_lengths[j] = 0;
	    }
	    improvedAnt.visited = new boolean[MTsp.n];
	    improvedAnt.toVisit = MTsp.n;
	    improvedAnt.longest_tour_length = Double.MAX_VALUE;
	    
	    improvedAnt.costObjectives = new double[2];
	    for (int indexObj = 0; indexObj < 2; indexObj++) {
	    	improvedAnt.costObjectives[indexObj] = 0;
    	}
	    
	    Ant temp = new Ant();
	    temp.tours = (ArrayList<Integer>[])new ArrayList[MTsp.m];
	    temp.tour_lengths = new double[MTsp.m];
	    for (int j = 0; j < MTsp.m; j++) {
	    	temp.tours[j] = new ArrayList<Integer>();
	    	temp.tour_lengths[j] = 0;
	    }
	    temp.visited = new boolean[MTsp.n];
	    temp.toVisit = MTsp.n;
	    temp.longest_tour_length = Double.MAX_VALUE;
		
	    temp.costObjectives = new double[2];
	    for (int indexObj = 0; indexObj < 2; indexObj++) {
	    	temp.costObjectives[indexObj] = 0;
    	}
    	
    	Ants.copy_from_to(a, improvedAnt);
    	Ants.copy_from_to(a, temp);
    	
    	while (foundImprovement) {
    		foundImprovement = false;
    		
    		Ants.copy_from_to(improvedAnt, a);
    		Ants.copy_from_to(a, temp);
    		
    		k = a.indexLongestTour;
	    	//iterate over cities from the longest tour except the first and last city from the tour which is the depot
	    	for (int i = 1; i < a.tours[k].size() - 1; i++) {   
	    		city1 = a.tours[k].get(i);
	    		temp.tours[k].remove(i);
	    		
	    		//iterate over the rest of other tours
	    		for (int l = 0; l < m; l++) {
	    			if (l != k) {
	    				for (int j = 1; j < a.tours[l].size(); j++) {
		    				temp.tours[l].add(j, city1);
		    				
		    				//compute the new subtour lengths and determine longest subtour
		    				temp.tour_lengths[k] = Tsp.compute_tour_length(temp.tours[k]);
		    				temp.tour_lengths[l] = Tsp.compute_tour_length(temp.tours[l]);
		    				
		    				longestTourLength = Double.MIN_VALUE;
		    				idLongestTour = 0;
		    				for (int index = 0; index < MTsp.m; index++) {
		    					if (longestTourLength < temp.tour_lengths[index]) {
		    						longestTourLength = temp.tour_lengths[index];
		    						idLongestTour = index;
		    					}
		    				}
		    				//verify if we obtain an improvement after this exchange operation
		    				if (longestTourLength < improvedAnt.longest_tour_length) {
		    					//System.out.println("Improvement found..");
		    					temp.longest_tour_length = longestTourLength;
		    					temp.indexLongestTour = idLongestTour;
		    					Ants.copy_from_to(temp, improvedAnt);
		    					foundImprovement = true;
		    				}
		    				
		    				//restore temp 
		    				//Ants.copy_from_to(iterBestAnt, temp);
		    				temp.tours[l].remove(j);
	    				}
	    			}
	    			
				}
	    		Ants.copy_from_to(a, temp);
	    	}
    	}
    	

    	return improvedAnt;
    }
    
    //manage some statistical information about the trial, especially if a new best solution
    //(best-so-far) is found and adjust some parameters if a new best solution is found
    static void update_statistics(boolean saveIterCosts) {
		int iteration_best_ant;
	
		iteration_best_ant = Ants.find_best(); /* iteration_best_ant is a global variable */
		
		//Ants.ants[iteration_best_ant] = local_search(Ants.ants[iteration_best_ant]);
		
		/*Ant a = Ants.ants[iteration_best_ant];
		//apply exchange local search operator for the iteration best ant
	    a = exchangeIterated(a);    
		
		//apply 2-opt heuristic on each of the m routes; 2-opt tries to improve single route
		//by replacing its two non-adjacent edges by two other edges
		for (int l = 0; l < m; l++) {
			LocalSearch.two_opt_first(a.tours[l]);
		}
		
		//compute new distances and update longest tour
		double sum = 0.0;
		double longestTourLength = Double.MIN_VALUE;
		int idLongestTour = 0;
		for (int l = 0; l < m; l++) {
			a.tour_lengths[l] = Tsp.compute_tour_length(a.tours[l]);
			sum += a.tour_lengths[l];
			if (longestTourLength < a.tour_lengths[l]) {
				longestTourLength = a.tour_lengths[l];
				idLongestTour = l;
			}
		}
		a.total_tour_length = sum;
		a.longest_tour_length = longestTourLength;
		a.indexLongestTour = idLongestTour;
		
        if (a.longest_tour_length < Ants.best_so_far_ant.longest_tour_length) {
			
		    InOut.time_used = Timer.elapsed_time();  //best solution found after time_used 
		    Ants.copy_from_to(a, Ants.best_so_far_ant);
	
		    InOut.found_best = InOut.iteration;	 
		    System.out.println("Iter: " + InOut.iteration + " Best ant -> longest tour=" + Ants.best_so_far_ant.longest_tour_length);
		}*/
	
		if (Ants.ants[iteration_best_ant].longest_tour_length < Ants.best_so_far_ant.longest_tour_length) {
	
		    InOut.time_used = Timer.elapsed_time();  //best solution found after time_used 
		    Ants.copy_from_to(Ants.ants[iteration_best_ant], Ants.best_so_far_ant);
	
		    InOut.found_best = InOut.iteration;	
		    
		    System.out.println("Iter: " + InOut.iteration + " Best ant -> longest tour=" + Ants.best_so_far_ant.longest_tour_length);
		}
		
		/*if (!saveIterCosts) {
			//compute non-dominated set of solutions (iteration non-dominated front)
			ParetoFront.iterationPareto.clear();
			Ant copyAnt;
			for (int i = 0; i < Ants.n_ants; i++) {
				copyAnt = Ants.copyAnt(Ants.ants[i]);
				ParetoFront.paretoUpdateWithSolution(ParetoFront.iterationPareto, copyAnt);
			}
			
			//update BestSoFarPareto external set
			ParetoFront.paretoUpdate(ParetoFront.bestSoFarPareto, ParetoFront.iterationPareto);
		}*/
		
    }

    //occasionally compute some statistics
   //at every 5 iterations save the value of the longest cost of the solution/tour given by the best so far ant
    static void search_control_and_statistics(ArrayList<Double> iterLongestCost, ArrayList<Integer> iterNumber, boolean saveDetailedOutput, int trial)
    {
    	double longestCost;
    	
    	if (saveDetailedOutput) {
    		if ((InOut.iteration % 5) == 0) {
			    //System.out.println("TSP(" + tspIndex + "): best tour length so far " + Ants.best_so_far_ant[tspIndex].tour_length + ", iteration: " + InOut.iteration);
	    		longestCost = Ants.best_so_far_ant.longest_tour_length;
	    		
	    		iterLongestCost.add(longestCost);
	    		if (trial == 0) {
	    			iterNumber.add(InOut.iteration);
	    		}
    		}
    	}
    }

    //manage global Ants.pheromone deposit for Ant System
    static void as_update() {
		int k;
	
		for (k = 0; k < Ants.n_ants; k++)
		    Ants.global_update_pheromone(Ants.ants[k]);
    }

    //manage global Ants.pheromone deposit for Ant Colony System
    static void acs_global_update() {
    	Ants.global_acs_pheromone_update(Ants.best_so_far_ant);
    }

    //manage global Ants.pheromone trail update for the ACO algorithms
	static void pheromone_trail_update()
	{
		//Simulate the Ants.pheromone evaporation of all Ants.pheromones; this is not necessary for ACS
		if (Ants.as_flag) {
			/* evaporate all Ants.pheromone trails */
			Ants.evaporation();
		}
	
		/* Next, apply the Ants.pheromone deposit for the various ACO algorithms */
		if (Ants.as_flag)
		    as_update();
		else if (Ants.acs_flag)
		    acs_global_update();
	
	
		/*
		 * Compute combined information Ants.pheromone times heuristic info after
		 * the Ants.pheromone update for all ACO algorithms except ACS; in the ACS case
		 * this is already done in the Ants.pheromone update procedures of ACS
		 */
		if (Ants.as_flag) {
			Ants.compute_total_information();
	    }
	}

	public static void main(String[] args) {
		boolean saveIterCosts = false;
		for (int trial = 0; trial < 10; trial++) {
			long startTime = System.currentTimeMillis();
			
	    	//read the data from the input file
			DataReader reader = new DataReader(fileName);
	        //read the data from the file
	        Data d = reader.read();
	        
	        String[] str = fileName.split("\\.");
			//we are dealing with a input file in TSP format 
			if (str[1].equals("tsp")) {
				instanceName = str[0];
			}
	        
	    	int n = d.getRecords().size();  //nr of cities except the depot city
	    	ArrayList<Record> records = d.getRecords();
	    	Record depotCity = d.getDepotCity();
	    	
	        MTsp mtsp = new MTsp(records, depotCity, n, m);
	        mtsp.printDepotCity();
	        //mtsp.printRecords(records);
	        
            InOut.init_program(args, trial);
			
			int[][][] result = new int[2][][];
			result = Tsp.compute_nn_lists();
		    MTsp.instance.nn_list = result[0];
		    MTsp.instance.nn_list_all = result[1];
			
			Ants.pheromone = new double[MTsp.n + 1][MTsp.n + 1];
			Ants.total = new double[MTsp.n + 1][MTsp.n + 1];
	
			init_try(); 
	
			ArrayList<Double> iterLongestCost = null;
			ArrayList<Integer> iterNumber = null;
			if (saveIterCosts) {
				//for saving detailed cost values at each 5 iterations
				iterLongestCost = new ArrayList<Double>();
				iterNumber = new ArrayList<Integer>();
			}
			
			InOut.iteration = 0;
		    while (!termination_condition()) {
				construct_solutions();
				if (ls_flag) {
					apply_local_search();
				}
				update_statistics(saveIterCosts);
				pheromone_trail_update();
				search_control_and_statistics(iterLongestCost, iterNumber, saveIterCosts, trial);
				InOut.iteration++;
		    }
		    
            /*System.out.println("Before local search->longest tour length:" + Ants.best_so_far_ant.longest_tour_length);
		    
		    //apply relocation local search operator for the iteration best ant
		    Ants.best_so_far_ant = relocationIterated(Ants.best_so_far_ant);  
			
			//apply exchange local search operator for the iteration best ant
		    Ants.best_so_far_ant = exchangeIterated(Ants.best_so_far_ant);    
			
			//apply 2-opt heuristic on each of the m routes; 2-opt tries to improve single route
			//by replacing its two non-adjacent edges by two other edges
			for (int l = 0; l < m; l++) {
				LocalSearch.two_opt_first(Ants.best_so_far_ant.tours[l]);
			}
			
			//compute new distances and update longest tour
			double sum = 0.0;
			double longestTourLength = Double.MIN_VALUE;
			int idLongestTour = 0;
			for (int l = 0; l < m; l++) {
				Ants.best_so_far_ant.tour_lengths[l] = Tsp.compute_tour_length(Ants.best_so_far_ant.tours[l]);
				sum += Ants.best_so_far_ant.tour_lengths[l];
				if (longestTourLength < Ants.best_so_far_ant.tour_lengths[l]) {
					longestTourLength = Ants.best_so_far_ant.tour_lengths[l];
					idLongestTour = l;
				}
			}
			Ants.best_so_far_ant.total_tour_length = sum;
			Ants.best_so_far_ant.longest_tour_length = longestTourLength;
			Ants.best_so_far_ant.indexLongestTour = idLongestTour;
		    
		    System.out.println("After local search->longest tour length:" + Ants.best_so_far_ant.longest_tour_length);*/
		    
		    //Utilities.setIterLongestCost(iterLongestCost);
		    //Utilities.setIterNumber(iterNumber);
		    
		    InOut.exit_try(trial);
		    if (!saveIterCosts) {
		    	////Utilities.writeParetoSet(ParetoFront.bestSoFarPareto, trial);
			    //System.out.println("Reached " + InOut.iteration + " iterations");
			    //Utilities.writeExcel(MTsp.n, MTsp.m, totalLength);
			    ////Utilities.writeParetoSolutions(ParetoFront.bestSoFarPareto);
			    ////ParetoFront.bestSoFarPareto.clear();
		    }
		    ////Utilities.writeResultsExcel(trial, saveIterCosts);	
		    
	       long endTime = System.currentTimeMillis();
		   double difference = (endTime - startTime)/1000.0;
		   System.out.println("\nRun #" + (trial + 1)  + " Elapsed seconds: " + difference);   
		}
		
    }
}
