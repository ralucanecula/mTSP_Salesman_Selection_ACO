package aco;

import java.util.ArrayList;
import java.util.HashMap;

/**
 * ACO algorithms for the TSP
 * 
 * This code is based on the ACOTSP project of Thomas Stuetzle.
 * It was initially ported from C to Java by Adrian Wilke.
 * 
 * Project website: http://adibaba.github.io/ACOTSPJava/
 * Source code: https://github.com/adibaba/ACOTSPJava/
 */
public class LocalSearch {
	/*
     * 
     * ################################################
     * ########## ACO algorithms for the TSP ##########
     * ################################################
     * 
     * Version: 1.0
     * File: ls.c
     * Author: Thomas Stuetzle
     * Purpose: implementation of local search routines
     * Check: README and gpl.txt
     * Copyright (C) 1999 Thomas Stuetzle
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
	
	//maximal depth of nearest neighbour lists used in the local search
	static int nn_ls = 20; 
  
	//flag indicating whether don't look bits are used
	static boolean dlb_flag = true; 
	
	/*
     * FUNCTION: generate a random permutation of the integers 0 .. n-1
     * INPUT: length of the array
     * OUTPUT: pointer to the random permutation
     * COMMENTS: only needed by the local search procedures
     */
    static int[] generate_random_permutation(ArrayList<Integer> t) {
		int i, help, node, tot_assigned = 0;
		double rnd;
		int[] r, permArray;
	
		int n = t.size() - 1;
		r = new int[n];
		permArray = new int[n];
	
		for (i = 0; i < n; i++)
		    r[i] = i;
	
		for (i = 0; i < n; i++) {
		    /* find (randomly) an index for a free unit */
		    rnd = Utilities.random01();
		    node = (int) (rnd * (n - tot_assigned));
		    assert (i + node < n);
		    help = r[i];
		    r[i] = r[i + node];
		    r[i + node] = help;
		    tot_assigned++;
		}
		for (i = 0; i < n; i++) {
			permArray[i] = t.get(r[i]);
		}
		return permArray;
    }

    /*
     * FUNCTION: 2-opt a tour
     * INPUT: pointer to the tour that undergoes local optimization
     * OUTPUT: none
     * (SIDE)EFFECTS: tour is 2-opt
     * COMMENTS: the neighbouAnts.rhood is scanned in random order (this need
     * not be the best possible choice). Concerning the speed-ups used
     * here consult, for example, Chapter 8 of
     * Holger H. Hoos and Thomas Stuetzle,
     * Stochastic Local Search---Foundations and Applications,
     * Morgan Kaufmann Publishers, 2004.
     * or some of the papers online available from David S. Johnson.
     */
    static void two_opt_first(ArrayList<Integer> tour) {
		boolean gotoExchange = false;
	
		int c1, c2; /* cities considered for an exchange */
		int s_c1, s_c2; /* successor cities of c1 and c2 */
		int p_c1, p_c2; /* predecessor cities of c1 and c2 */
		int pos_c1, pos_c2; /* positions of cities c1, c2 */
		int i, j, h, l;
		int help;
		boolean improvement_flag;
		int countImprovement;
		int h1 = 0, h2 = 0, h3 = 0, h4 = 0;
		double radius; /* radius of nn-search */
		double gain = 0.0;
		int[] random_vector;
		HashMap<Integer, Integer> pos; /* positions of cities in tour */
		HashMap<Integer, Boolean> dlb;  /* contains don't look bits */
		
		pos = new HashMap<Integer, Integer>(tour.size() - 1);
		dlb = new HashMap<Integer, Boolean>(tour.size() - 1);
		for (i = 0; i < tour.size() - 1; i++) {
		    pos.put(tour.get(i), i);
		    dlb.put(tour.get(i), false);
		}
	
		improvement_flag = true;
		countImprovement = 0;
		random_vector = generate_random_permutation(tour);
	
		while (improvement_flag) {
			if (countImprovement > 50) {
				//System.out.println("Exiting from loop..");
				break;
			}
		    improvement_flag = false;
	
		    for (l = 0; l < tour.size() - 1; l++) {
				c1 = random_vector[l];
				if (dlb_flag && dlb.get(c1))
				    continue;
				pos_c1 = pos.get(c1);
				s_c1 = tour.get(pos_c1 + 1);
				radius = MTsp.instance.distance[c1 + 1][s_c1 + 1];
		
				/* First search for c1's nearest neighbours, use successor of c1 */
				for (h = 0; h < nn_ls; h++) {
				    c2 = MTsp.instance.nn_list_all[c1 + 1][h]; /* exchange partner, determine its position */
				    //make sure that c2 is in the same subtour with c1; if it isn't choose another element
				    if (pos.get(c2 - 1) == null) {
				    	continue;
				    }
				    if (radius > MTsp.instance.distance[c1 + 1][c2]) {
						s_c2 = tour.get(pos.get(c2 - 1) + 1);
						gain = -radius + MTsp.instance.distance[c1 + 1][c2] + MTsp.instance.distance[s_c1 + 1][s_c2 + 1]
							- MTsp.instance.distance[c2][s_c2 + 1];
						if (gain < 0) {
						    h1 = c1;
						    h2 = s_c1;
						    h3 = c2 - 1;
						    h4 = s_c2;
						    gotoExchange = true;
						    break;
						}
				    } else
				    	break;
				}
		
				if (gotoExchange) {
				    /* Search one for next c1's h-nearest neighbours, use predecessor c1 */
				    if (pos_c1 > 0)
				    	p_c1 = tour.get(pos_c1 - 1);
				    else
				    	p_c1 = tour.get(tour.size() - 2);
				    radius = MTsp.instance.distance[p_c1 + 1][c1 + 1];
				    for (h = 0; h < nn_ls; h++) {
						c2 = MTsp.instance.nn_list_all[c1 + 1][h]; /* exchange partner, determine its position */
						//make sure that c2 is in the same subtour with c1; if it isn't choose another element
					    if (pos.get(c2 - 1) == null) {
					    	continue;
					    }
						if (radius > MTsp.instance.distance[c1 + 1][c2]) {
						    pos_c2 = pos.get(c2 - 1);
						    if (pos_c2 > 0)
						    	p_c2 = tour.get(pos_c2 - 1);
						    else
						    	p_c2 = tour.get(tour.size() - 2);
						    if (p_c2 == c1)
						    	continue;
						    if (p_c1 == (c2 - 1))
						    	continue;
						    gain = -radius + MTsp.instance.distance[c1 + 1][c2] + MTsp.instance.distance[p_c1 + 1][p_c2 + 1]
							    - MTsp.instance.distance[p_c2 + 1][c2];
						    if (gain < 0) {
								h1 = p_c1;
								h2 = c1;
								h3 = p_c2;
								h4 = c2 - 1;
								gotoExchange = true;
								break;
						    }
						} else
						    break;
				    }
				}
		
				if (!gotoExchange) {
				    /* No exchange */
				    dlb.put(c1, true);
				    continue;
				}
		
				if (gotoExchange) {
				    gotoExchange = false;
				    improvement_flag = true;
				    countImprovement++;
				    //System.out.println("countImprovement=" + countImprovement);
				    dlb.put(h1, false);
				    dlb.put(h2, false);
				    dlb.put(h3, false);
				    dlb.put(h4, false);
				    /* Now perform move */
				    if (pos.get(h3) < pos.get(h1)) {
						help = h1;
						h1 = h3;
						h3 = help;
						help = h2;
						h2 = h4;
						h4 = help;
				    }
				    if (pos.get(h3) - pos.get(h2) < (tour.size() - 1) / 2 + 1) {
						/* reverse inner part from pos[h2] to pos[h3] */
						i = pos.get(h2);
						j = pos.get(h3);
						//System.out.println("i=" + i + " j=" + j);
						if (i >= j) {
							//System.out.println("Got here..");
							//countImprovement--;
						}
						while (i < j) {
						    c1 = tour.get(i);
						    c2 = tour.get(j);
						    tour.set(i, c2);
						    tour.set(j, c1);
						    pos.put(c1, j);
						    pos.put(c2, i);
						    i++;
						    j--;
						}
						//System.out.println("Improvement applied->case 1..");
				    } else {
						/* reverse outer part from pos[h4] to pos[h1] */
						i = pos.get(h1);
						j = pos.get(h4);
						if (j > i)
						    help = tour.size() - (j - i);
						else
						    help = (i - j) + 1;
						help = help / 2;
						for (h = 0; h < help; h++) {
						    c1 = tour.get(i);
						    c2 = tour.get(j);
						    tour.set(i, c2);
						    tour.set(j, c1);
						    pos.put(c1, j);
						    pos.put(c2, i);
						    i--;
						    j++;
						    if (i < 0)
						    	i = tour.size() - 2;
						    if (j >= tour.size() - 1)
						    	j = 0;
						}
						tour.set(tour.size() - 1, tour.get(0));
						//System.out.println("Improvement applied->case 2..");
				    }
				} else {
				    dlb.put(c1, true);
				}
	
		    }
		}

    }
}
