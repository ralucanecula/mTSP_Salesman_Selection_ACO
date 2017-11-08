package aco;

import java.util.ArrayList;

/**
 * ACO algorithms for the TSP
 * 
 * This code is based on the ACOTSP project of Thomas Stuetzle.
 * It was initially ported from C to Java by Adrian Wilke.
 * 
 * Project website: http://adibaba.github.io/ACOTSPJava/
 * Source code: https://github.com/adibaba/ACOTSPJava/
 */


public class Ants {
	
	 /*
     * ################################################
     * ########## ACO algorithms for the TSP ##########
     * ################################################
     * 
     * Version: 1.0
     * File: ants.c
     * Author: Thomas Stuetzle
     * Purpose: implementation of procedures for ants' behaviour
     * Check: README.txt and legal.txt
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

    static class Ant {
    	//for each of the m salesmen an ant will construct a tour, so that a candidate solution constructed by
    	//an ant will be represented by a list of tours, one for each salesman
    	ArrayList<Integer>[] tours;
		boolean[] visited;
		double[] tour_lengths;
		double total_tour_length;
		double longest_tour_length;
		//cities left to be visited by an ant (initially toVisit = n, which is the number of cities from the mTSP instance)
		int toVisit;
		//stores the cost of each solution according to the considered objectives (2 in this case)
		double costObjectives[];
    }

    public static final int MAX_ANTS = 1024;
    public static final int MAX_NEIGHBOURS = 512;

    static Ant ants[];
    static Ant best_so_far_ant;

    static double pheromone[][];
    static double total[][];  //keeps heuristic information times pheromone for each arc

    static double prob_of_selection[];

    static int n_ants; /* number of ants */
 
    static int nn_ants; /* length of nearest neighbor lists for the ants' solution construction */

    static double rho; /* parameter for evaporation */
    static double alpha; /* importance of trail */
    static double beta; /* importance of heuristic evaluate */
    static double q_0; /* probability of best choice in tour construction */

    static boolean as_flag; /* ant system */
    static boolean acs_flag; /* ant colony system (ACS) */

    static int u_gb; /* every u_gb iterations update with best-so-far ant */

    static double trail_0; /* initial pheromone level in ACS */

    static double HEURISTIC(int m, int n) {
    	return (1.0 / (double) MTsp.instance.distance[m][n]);
    }

    //allocate the memory for the ant colony, the best-so-far ant
    static void allocate_ants() {
		int i, j;
	
		ants = new Ant[n_ants];
	
		for (i = 0; i < n_ants; i++) {
		    ants[i] = new Ant();
		    ants[i].tours = (ArrayList<Integer>[])new ArrayList[MTsp.m];
		    ants[i].tour_lengths = new double[MTsp.m];
		    for (j = 0; j < MTsp.m; j++) {
		    	ants[i].tours[j] = new ArrayList<Integer>();
		    	ants[i].tour_lengths[j] = 0;
		    }
		    ants[i].visited = new boolean[MTsp.n];
		    ants[i].toVisit = MTsp.n;
		    ants[i].costObjectives = new double[2];
		    for (int indexObj = 0; indexObj < 2; indexObj++) {
		    	ants[i].costObjectives[indexObj] = 0;
	    	}
		}
		
		best_so_far_ant = new Ant();
		best_so_far_ant.tours = (ArrayList<Integer>[])new ArrayList[MTsp.m];
		best_so_far_ant.tour_lengths = new double[MTsp.m];
	    for (j = 0; j < MTsp.m; j++) {
	    	best_so_far_ant.tours[j] = new ArrayList<Integer>();
	    	best_so_far_ant.tour_lengths[j] = 0;
	    }
		best_so_far_ant.visited = new boolean[MTsp.n];
		best_so_far_ant.toVisit = MTsp.n;
		best_so_far_ant.longest_tour_length = Double.MAX_VALUE;
		
		best_so_far_ant.costObjectives = new double[2];
	    for (int indexObj = 0; indexObj < 2; indexObj++) {
	    	best_so_far_ant.costObjectives[indexObj] = 0;
    	}
	
		prob_of_selection = new double[nn_ants + 1];
		for (i = 0; i < nn_ants + 1; i++) {
		    prob_of_selection[i] = Double.POSITIVE_INFINITY;
		}
    }

    // find the best ant of the current iteration (the one with the lowest cost of its longest tour)
    static int find_best()
    {
    	double min;
		int k, k_min;
	
		min = ants[0].longest_tour_length;
		k_min = 0;
		for (k = 1; k < n_ants; k++) {
		    if (ants[k].longest_tour_length < min) {
				min = ants[k].longest_tour_length;
				k_min = k;
		    }
		}
		return k_min;
    }

    //initialize pheromone trails
    //matricea cu urmele de feromoni trebuie sa se faca relativ la toate cele n orase
    static void init_pheromone_trails(double initial_trail)
    {
		int i, j;
	
		/* Initialize pheromone trails */
		for (i = 0; i < (MTsp.n + 1); i++) {
		    for (j = 0; j <= i; j++) {
				pheromone[i][j] = initial_trail;
				pheromone[j][i] = initial_trail;
				total[i][j] = initial_trail;
				total[j][i] = initial_trail;
		    }
		}
    }

    //implements the pheromone trail evaporation
    static void evaporation()
    {
		int i, j;
		
		for (i = 0; i < MTsp.n + 1; i++) {
		    for (j = 0; j <= i; j++) {
				pheromone[i][j] = (1 - rho) * pheromone[i][j];
				pheromone[j][i] = pheromone[i][j];
		    }
		}
    }

    //reinforces edges used in ant k's solution
    static void global_update_pheromone(Ant a)
    {
		int i, j, h, k, size;
		double d_tau;
	
		d_tau = 1.0 / (double) a.total_tour_length;
		for (i = 0; i < MTsp.m; i++) {
			size = a.tours[i].size();
			for (k = 0; k < size - 1; k++) {
			    j = a.tours[i].get(k); 
			    h = a.tours[i].get(k + 1);
			    
			    j++;
	            h++;
	            
			    pheromone[j][h] += d_tau;
			    pheromone[h][j] = pheromone[j][h];
			}
		}
		
    }

    //calculates heuristic info times pheromone for each arc
    static void compute_total_information()
    {
		int i, j;
	
		for (i = 0; i < (MTsp.n + 1); i++) {
		    for (j = 0; j < i; j++) {
				total[i][j] = Math.pow(pheromone[i][j], alpha) * Math.pow(HEURISTIC(i, j), beta);
				total[j][i] = total[i][j];
		    }
		}
    }

    //empty the ants's memory regarding visited cities
    static void ant_empty_memory(Ant a)
    {
        int i, j;

		a.total_tour_length = 0;
		a.longest_tour_length = Integer.MAX_VALUE;
		
		for (int indexObj = 0; indexObj < 2; indexObj++) {
	    	a.costObjectives[indexObj] = 0;
	    }
		
		//clear all the elements (cities) from the m tours of an ant
		for (i = 0; i < MTsp.m; i++) {
		    a.tours[i].clear();
		    a.tour_lengths[i] = 0;
		}
		
		for (j = 0; j < MTsp.n; j++) {
		    a.visited[j] = false;
		}
		a.toVisit = MTsp.n;
		
    }
    
    //create a copy of the ant a and return the created copy at the output
    static Ant copyAnt(Ant a) {
    	//first create an empty ant
    	Ant copy = new Ant();
    	copy.tours = (ArrayList<Integer>[])new ArrayList[MTsp.m];
    	copy.tour_lengths = new double[MTsp.m];
	    for (int j = 0; j < MTsp.m; j++) {
	    	copy.tours[j] = new ArrayList<Integer>();
	    	copy.tour_lengths[j] = 0;
	    }
	    copy.visited = new boolean[MTsp.n];
	    copy.toVisit = MTsp.n;
	    
	    copy.costObjectives = new double[2];
	    //copy.weights = new double[TSP_ACO.k];
    	
	    //then copy the information from the ant a
	    copy.total_tour_length = a.total_tour_length;
		for (int indexObj = 0; indexObj < 2; indexObj++) {
			copy.costObjectives[indexObj] = a.costObjectives[indexObj];
    	}
		for (int i = 0; i < MTsp.m; i++) {
			copy.tour_lengths[i] = a.tour_lengths[i];
			int size = a.tours[i].size();
			copy.tours[i] = new ArrayList<Integer>(size);
			for (int j = 0; j < size; j++) {
				int elem = a.tours[i].get(j);
				copy.tours[i].add(elem);
			}		
		}
    	
    	return copy;
    }

    //choose for an ant as the next city the one with maximal value of heuristic information times pheromone
    static void choose_best_next(Ant a, int indexSalesman)
    {
		int city, current_city, next_city;
		double value_best;
	
		next_city = MTsp.n;
		int lastPos = a.tours[indexSalesman].size() - 1;
		current_city = a.tours[indexSalesman].get(lastPos);
    	current_city++;
    	
		value_best = -1.; /* values in total matrix are always >= 0.0 */
		for (city = 0; city < MTsp.n; city++) {
		    if (a.visited[city])
			; /* city already visited, do nothing */
		    else {
				if (total[current_city][city + 1] > value_best) {
				    next_city = city;
				    value_best = total[current_city][city + 1];
				}
		    }
		}
		a.tours[indexSalesman].add(next_city);
		a.visited[next_city] = true;
		a.toVisit--;
    }

    //chooses for an ant as the next city the one with maximal value of heuristic information times pheromone
    static void neighbour_choose_best_next(Ant a, int indexSalesman) {
		int i, current_city, next_city, help_city;
		double value_best, help;
	
		next_city = MTsp.n;   //next_city = Integer.MAX_VALUE;
		int lastPos = a.tours[indexSalesman].size() - 1;
		current_city = a.tours[indexSalesman].get(lastPos);
        current_city++;
 
		value_best = -1.;  //values in total matrix are always >= 0.0 
		for (i = 0; i < nn_ants; i++) {
		    help_city = MTsp.instance.nn_list[current_city][i];
		    if (a.visited[help_city - 1])
		    	; // city already visited, do nothing 
		    else {
				help = total[current_city][help_city];
				if (help > value_best) {
				    value_best = help;
				    next_city = help_city - 1;
				}
		    }
		}
		
		if (next_city == MTsp.n)
		    // all cities in nearest neighbor list were already visited 
		    choose_best_next(a, indexSalesman);
		else {
			a.tours[indexSalesman].add(next_city);
		    a.visited[next_city] = true;
		    a.toVisit--;
		}
    }
    
    static void choose_closest_next(Ant a, int indexSalesman)
    {
		int current_city, next_city, city;
	    double  min_distance;
		
		next_city = MTsp.n;
		int lastPos = a.tours[indexSalesman].size() - 1;
		current_city = a.tours[indexSalesman].get(lastPos);
		current_city++;
		
		min_distance = Integer.MAX_VALUE;  //Search shortest edge 
		for (city = 0; city < MTsp.n; city++) {
		    if (a.visited[city])
			;  //city already visited 
		    else {
				if (MTsp.instance.distance[current_city][city + 1] < min_distance) {
				    next_city = city;
				    min_distance = MTsp.instance.distance[current_city][city + 1];
				}
		    }
		}
		a.tours[indexSalesman].add(next_city);
		a.visited[next_city] = true; 
		a.toVisit--;
    }

    //Choose for an ant probabilistically a next city among all unvisited cities in the current city's candidate list
    static void neighbour_choose_and_move_to_next(Ant a, int indexSalesman) {
		int i, help, city;
		int current_city;
		double rnd, partial_sum = 0., sum_prob = 0.0;
		double prob_ptr[];
	
		if ((q_0 > 0.0) && (Utilities.random01() < q_0)) {
		    /*
		     * with a probability q_0 make the best possible choice
		     * according to pheromone trails and heuristic information, this corresponds to exploitation
		     */
		    /*
		     * we first check whether q_0 > 0.0, to avoid the very common case
		     * of q_0 = 0.0 to have to compute a random number, which is
		     * expensive computationally
		     */
		    neighbour_choose_best_next(a, indexSalesman);
		    return;
		}
	
		prob_ptr = prob_of_selection;
	
	    /* current_city city of ant k */
		int lastPos = a.tours[indexSalesman].size() - 1;
		current_city = a.tours[indexSalesman].get(lastPos);
    	current_city++;
    
		for (i = 0; i < nn_ants; i++) {
			city = MTsp.instance.nn_list[current_city][i];
		    if (a.visited[city - 1])
		    	prob_ptr[i] = 0.0; /* city already visited */
		    else {
				prob_ptr[i] = total[current_city][city];
				sum_prob += prob_ptr[i];
		    }
		}
	
		if (sum_prob <= 0.0) {
		    /* All cities from the candidate set are tabu (are already visited) */
		    choose_best_next(a, indexSalesman);
		} else {
		    /*
		     * at least one neighbor is eligible, choose one according to the
		     * selection probabilities
		     */
		    rnd = Utilities.random01();
		    rnd *= sum_prob;
		    i = 0;
		    partial_sum = prob_ptr[i];
		    /* This loop always stops because prob_ptr[nn_ants] == HUGE_VAL */
		    while (partial_sum <= rnd) {
				i++;
				partial_sum += prob_ptr[i];
		    }
		    /*
		     * This may very rarely happen because of rounding if rnd is close to 1.
		     */
		    if (i == nn_ants) {
				neighbour_choose_best_next(a, indexSalesman);
				return;
		    }
		    help = MTsp.instance.nn_list[current_city][i];
		    a.tours[indexSalesman].add(help - 1);
		    a.visited[help - 1] = true;
		    a.toVisit--;
		}
    }

    //reinforces the edges used in ant's solution as in ACS
    static void global_acs_pheromone_update(Ant a) {
		int i, j, h, k, size;
		double d_tau;
	
		d_tau = 1.0 / (double) a.longest_tour_length;
	
		for (i = 0; i < MTsp.m; i++) {
			size = a.tours[i].size();
			for (k = 0; k < size - 1; k++)  {
			    j = a.tours[i].get(k);
			    h = a.tours[i].get(k + 1);	
			    
			    j++;
			    h++;
		
			    pheromone[j][h] = (1. - rho) * pheromone[j][h] + rho * d_tau;
			    pheromone[h][j] = pheromone[j][h];
		
			    total[h][j] = Math.pow(pheromone[h][j], alpha) * Math.pow(HEURISTIC(h, j), beta);
			    total[j][h] = total[h][j];
			}
		}
		
    }

    //removes some pheromone on edge just passed by the ant
    static void local_acs_pheromone_update(Ant a, int indexSalesman) {
		int h, j;
		
		int lastPos = a.tours[indexSalesman].size() - 1;
		j = a.tours[indexSalesman].get(lastPos);
		h = a.tours[indexSalesman].get(lastPos - 1);	
		
		j++;
		h++;

		/* still additional parameter has to be introduced */
		pheromone[h][j] = (1. - 0.1) * pheromone[h][j] + 0.1 * trail_0;
		pheromone[j][h] = pheromone[h][j];
		total[h][j] = Math.pow(pheromone[h][j], alpha) * Math.pow(HEURISTIC(h, j), beta);
		total[j][h] = total[h][j];
    }

    //copy solution from ant a1 into ant a2
    static void copy_from_to(Ant a1, Ant a2) {
		int i, j;
	
		Ants.ant_empty_memory(a2);
		
		a2.total_tour_length = a1.total_tour_length;
		a2.longest_tour_length = a1.longest_tour_length;
		for (int indexObj = 0; indexObj < 2; indexObj++) {
		   a2.costObjectives[indexObj] = a1.costObjectives[indexObj];
    	}
		for (i = 0; i < MTsp.m; i++) {
			a2.tour_lengths[i] = a1.tour_lengths[i];
			int size = a1.tours[i].size();
			a2.tours[i] = new ArrayList<Integer>(size);
			for (j = 0; j < size; j++) {
				int elem = a1.tours[i].get(j);
				a2.tours[i].add(elem);
			}
			
		}
    }
    
    static double computeToursAmplitude(Ant a) {
    	double min, max;
		int i;

		min = a.tour_lengths[0];
		max = a.tour_lengths[0];
		for (i = 1; i < a.tours.length; i++) {
		    if (a.tour_lengths[i] < min) {
				min = a.tour_lengths[i];
		    }
		    if (a.tour_lengths[i] > max) {
				max = a.tour_lengths[i];
		    }	    
		}
		
		return (max - min);
    }

    //generate some nearest neighbor tour and compute tour length
    static double nn_tour() {
    	int k, step, salesman;
    	double sum = 0;
    	
    	ant_empty_memory(ants[0]);
    	step = 0;
    	
    	for (int i = 0; i < MTsp.m; i++) {
    		//place the ant on the depot city, which is the start city of each tour
			// -1 is a special marker for the deport city, so that it's not be confused with the rest of the cities
			// all the rest of the cities are represented by integer values > 0
    		ants[0].tours[i].add(-1);  	
    	}
		
		while (ants[0].toVisit > 0) {   //there are still left cities to be visited	
    		//choose in a random manner with equal probability which salesman to consider next, that will visit a city
	    	salesman = (int)(Math.random() * MTsp.m);
	    	choose_closest_next(ants[0], salesman);	 
		}
		
		double longestTourLength = Double.MIN_VALUE;
		for (int i = 0; i < MTsp.m; i++) {
			step = ants[0].tours[i].size();
			ants[0].tours[i].add(step, -1);
			
			ants[0].tour_lengths[i] = Tsp.compute_tour_length_(ants[0].tours[i]);
			ants[0].total_tour_length += ants[0].tour_lengths[i];
			if (longestTourLength < ants[0].tour_lengths[i]) {
				longestTourLength = ants[0].tour_lengths[i];
			}
		}
		ants[0].longest_tour_length = longestTourLength;
		//sum = ants[0].total_tour_length;
		ant_empty_memory(ants[0]);
		
		//return sum;
		System.out.println("Initial (nearest neighbour tour) longest tour length:" + longestTourLength);
		return longestTourLength;
    }


    //generate some nearest neighbor tour and compute tour length
    /*static int nn_tour_(Partition[] clusters) 
    {
		int phase, sum = 0;
		
		//for each of the m crisp partitions obtained with fuzzy cMeans clustering algorithm, that define a 
		//TSP problem, a nearest neighbor tour should be constructed and at the end, the total sum of these m
		//tour lengths should be computed and returned as the length of the nn_tour
	    for (int i = 0; i < clusters.length; i++) {	    	
	    	ant_empty_memory(ants[0]);
	    	phase = 0;  //counter of the construction steps 
			//place_ant(i, ants[0], phase);
	    	ants[0].tours[i].add(phase, -1);
		
			while (phase < MTsp.clusterCities[i].size()) {
			    phase++;
			    choose_closest_next_(i, ants[0], phase);
			}
			phase = MTsp.clusterCities[i].size() + 1;
			ants[0].tours[i].add(phase, -1);
			//we assume the tour contains at the first and last position the depot city, so that the total length of the tour 
			//should take into account the depot city
			InOut.n_tours += 1;
			ants[0].total_tour_length = Tsp.compute_tour_length(ants[0].tours[i]);
		
			sum += ants[0].total_tour_length;
			ant_empty_memory(ants[0]);
	    }
		
		return sum;
    }*/


}
