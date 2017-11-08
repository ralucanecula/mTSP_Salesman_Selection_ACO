package aco;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.Reader;

import aco.Tsp.Problem;

/**
 * ACO algorithms for the TSP
 * 
 * This code is based on the ACOTSP project of Thomas Stuetzle.
 * It was initially ported from C to Java by Adrian Wilke.
 * 
 * Project website: http://adibaba.github.io/ACOTSPJava/
 * Source code: https://github.com/adibaba/ACOTSPJava/
 */
public class InOut {

	  /*
     * ################################################
     * ########## ACO algorithms for the TSP ##########
     * ################################################
     * 
     * Version: 1.0
     * File: InOut.c
     * Author: Thomas Stuetzle
     * Purpose: mainly input / output / statistic routines
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
	
    enum Distance_type {EUC_2D, CEIL_2D, GEO, ATT};

    static Distance_type distance_type;

    static int best_in_try;
    static int best_found_at;
    static double time_best_found;
    static double time_total_run;

    static int n_try; /* try counter */
    static int n_tours; /* counter of number constructed tours */
    static int iteration; /* iteration counter */
    static int max_tries; /* maximum number of independent tries */
    static int max_tours; /* maximum number of tour constructions in one try */
    public static int max_iterations; /* maximum number of iterations */

    static double max_time; /* maximal allowed run time of a try */
    static double time_used; /* time used until some given event */
    static double time_passed; /* time passed until some moment */
    static int optimal; /* optimal solution or bound to find */

    static int found_best; /* iteration in which best solution is found for each instance of the TSP problem */

    static String inputFile;
    static int opt;

    static void set_default_as_parameters() {
		/* number of ants (-1 means MTsp.instance.n size) and number of nearest neighbors in tour construction */
		Ants.n_ants = -1;
		Ants.nn_ants = 20;
		
		Ants.alpha = 1.0;
		Ants.beta = 2.0;
		Ants.rho = 0.5;
		Ants.q_0 = 0.0;
    }

    static void set_default_acs_parameters() {
		/* number of ants (-1 means MTsp.instance.n size) and number of nearest neighbors in tour construction */
		Ants.n_ants = 15;
		Ants.nn_ants = 15;
 
		Ants.alpha = 1.0;
		Ants.beta = 2.0;
		Ants.rho = 0.1;
		Ants.q_0 = 0.9;
    }
    
    //set default parameter settings
    static void set_default_parameters() {
		/* number of ants and number of nearest neighbors in tour construction */
		Ants.n_ants = 25;
		Ants.nn_ants = 20; 
	
		Ants.alpha = 1.0;
		Ants.beta = 2.0;
		Ants.rho = 0.5;
		Ants.q_0 = 0.0;
		max_tries = 10;
		max_tours = 10 * 20;
		Utilities.seed = (int) System.currentTimeMillis();
		//10 seconds allowed for running time; it is used in the termination condition of ACO
		max_time = 15.0 * MTsp.m; 
		//maximum number of allowed iterations until the ACO algorithm stops
		max_iterations = 2200;
		optimal = 1;
		Ants.u_gb = Integer.MAX_VALUE;
		Ants.acs_flag = false;
		Ants.as_flag = false;
		distance_type = Distance_type.EUC_2D;
    }
    
    public static float average(int[] array) {
    	int sum = 0;

    	
    	for (int  i = 0; i < array.length; i++)
    		sum += array[i];
    	
    	return (float)sum / (float)array.length;
    }
    	 
	public static float variance(int[] array){
		double var = 0;

		double average = average(array);
		for (int  i = 0; i < array.length; i++)
		   var += (array[i] - average)*(array[i] - average);
		
		return (float)var/(float)array.length;
    }

    // save some statistical information on a trial once it finishes
    static void exit_try(int ntry) {
    	//String[] toursText = new String[MTsp.m];
    	//int[] nrCities = new int[MTsp.m];
    	//double[] subtoursCost = new double[MTsp.m];
    	double totalCost, longestSubtour, amplitude;
    	
    	System.out.println("\nRun #" + (ntry + 1) + ": " + ParetoFront.bestSoFarPareto.size() + " solutions found in the best so far Pareto set");
    	
    	System.out.println("\nBest so far ant: final tours traveled by each salesman are:");
		for (int index = 0; index < MTsp.m; index++) {
			//System.out.print("\ntour(" + index + "): ");
			//print depot city
			System.out.print("1 ");
			StringBuilder sb = new StringBuilder();
			sb.append("1 ");
			
			int tourLength =  Ants.best_so_far_ant.tours[index].size();
			int count = 0;
			for (int i = 1; i < tourLength - 1; i++) {
				int city = Ants.best_so_far_ant.tours[index].get(i);
				city = city + 2;  //so as to correspond to the city indexes from the mtsp input file
				System.out.print(city + " ");
				sb.append(city + " ");
			}
			count = tourLength - 2; //excluding the 2 depots cities from the start and end of tour
			System.out.print("1  (#" + count + ") Total cost: " + Ants.best_so_far_ant.tour_lengths[index] + "\n");
			sb.append("1  (#" + count + ")\n");
			//subtoursCost[index] = Ants.best_so_far_ant.tour_lengths[index];
			//toursText[index] = sb.toString();
			//nrCities[index] = count;
		}
		//System.out.println("Variance of tour costs is: " + variance(costs));

		double length = 0;
		System.out.println("\nBest solution is " + Ants.best_so_far_ant.total_tour_length + ", longest tour is " + Ants.best_so_far_ant.longest_tour_length + ", found at iteration " + found_best);
	    length = Ants.best_so_far_ant.total_tour_length;
	    totalCost = length;
	    longestSubtour = Ants.best_so_far_ant.longest_tour_length;
	    amplitude = Ants.best_so_far_ant.costObjectives[1];
	    
		System.out.println("\nTotal sum of tour lenghts is: " + length);
		
		/*Utilities.setTours(toursText);
		Utilities.setNrCities(nrCities);
		Utilities.setSubtoursCost(subtoursCost);*/
		
		Utilities.setTotalCost(totalCost);
		Utilities.setLongestSubtour(longestSubtour);
		Utilities.setAmplitude(amplitude);
		
    }

    //initialize the program
    static void init_program(String[] args, int runNumber) {
	
		set_default_parameters();
		Parse.parse_commandline(args, runNumber);
		
		//compute distance matrix between cities and allocate ants
    	if (Ants.n_ants < 0)
    		Ants.n_ants = MTsp.n;
 
    	MTsp.instance.distance = Tsp.compute_distances();
    	Ants.allocate_ants();
		
    }

}
