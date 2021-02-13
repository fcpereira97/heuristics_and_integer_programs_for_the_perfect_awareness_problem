/*
File: random_plus_greedy.cpp
Description: GRASP construction phase for PAP with Random plus Greedy strategie
Author: Felipe de C. Pereira

Based on the article "Effective Heuristics for the Perfect Awareness Problem"
Copyright (C) 2020 Felipe de C. Pereira, Pedro J. de Rezebde, Cid C. de Souza
*/

#include <bits/stdc++.h>
#include "my_struct.h"
#include "propagation.h"
#include "priority_queue.h"

using namespace std;

// GRASP construction phase with Random plus Greedy strategie
void random_plus_greedy_construction(int n_vertices, vector<Vertex*> &vertices, vector<Vertex*> &seed_set,
	int &sol_value, int &n_rounds, int &n_aware, int n_seeds_per_insertion, double random_steps_rate, 
	int &draws_altruistic, int &draws_philanthropic, int &steps, int &tot_spreaders, int &tot_edges_explored,
	double &rlc_cl_ratio_mean, mt19937 &rng)
{
	int round, random_steps; // Control variables
	queue<Vertex*> new_seeds; // New selected seeds
	vector<Vertex*> cl(n_vertices); // Candidate list (CL)
	double alpha; // Greediness parameter
	vector<int> null_array;

	// Initialize variables
	n_aware = 0;
	round = 0;
	alpha = 1.0;
	random_steps = ceil(random_steps_rate * n_vertices);

	// Erase the propagation process
	erase_propagation(n_vertices, vertices, 3);

	// Turn to seeds all vertices that has degree equals to 0
	initialize_natural_seeds(n_vertices, vertices, seed_set, n_aware);

	// Build CL as a priority queue due to the benefits of the CL elements
	cl.assign(vertices.begin(), vertices.end());
	make_heap(cl.begin(), cl.end(), comp_pipeline);

	// Update vertices positions in priority queue
	for(int i = 0; i < (int)cl.size(); i++)
		cl[i]-> cl_pos = i;

	// While the seed set is still a infeasible solution do
	while(n_aware < n_vertices)
	{
		// For the first random_steps steps, we have alpha = 1.0, that is, total randomness
		// For the remaining steps, we have alpha = 0.0, that is, total greediness
		if(steps > random_steps - 1)
			alpha = 0.0;

		// For as many seeds per selection as wanted
		for(int i = 0; i < n_seeds_per_insertion; i++)
		{
			// Selected vertice variable
			Vertex *v_chosen;

			if(alpha > 0)
			{
				// Draw a random vertex from cl
				v_chosen = cl[rng() % cl.size()];
			}
			else
			{
				// Get the best vertex from cl
				v_chosen = cl[0];

				// Update statistics about ties
				if( (cl.size() >= 2 && v_chosen-> n_unaware_neighs == cl[1]-> n_unaware_neighs) ||
					(cl.size() >= 3 && v_chosen-> n_unaware_neighs == cl[2]-> n_unaware_neighs))
				{
					draws_altruistic++;
					if( (cl.size() >= 2 && v_chosen-> n_unaware_neighs == cl[1]-> n_unaware_neighs
						&& v_chosen-> n_spreader_neighs_gain == cl[1]-> n_spreader_neighs_gain) ||
						(cl.size() >= 3 && v_chosen-> n_unaware_neighs == cl[2]-> n_unaware_neighs
						&& v_chosen-> n_spreader_neighs_gain == cl[2]-> n_spreader_neighs_gain) )
						draws_philanthropic++;
				}
			}	

			// If the selected vertex is not spreader
			if(v_chosen-> state != 2)
			{
				// If the selected vertex is unaware
				if(v_chosen-> state == 0)
				{
					n_aware++; // Update the total number of aware vertices
					
					// Update the number of unaware neighs of the selected vertex neighbors
					decrease_n_unaware_neighs(v_chosen, cl, null_array, 3);
				}

				// Turn selected vertex into a seed
				// Update the cl and the benefits
				v_chosen-> state = 2;
				v_chosen-> isSeed = true;
				if(v_chosen -> lack_threshold == 1)
					decrease_n_spreader_neighs_gain(v_chosen, cl, 3);
				v_chosen-> lack_threshold = 0;
				seed_set.push_back(v_chosen);
				new_seeds.push(v_chosen);
				delete_node_from_heap(cl, v_chosen -> cl_pos, comp_pipeline);
			}
		}

		// Continue the propagation with the new seeds 
		propagate_from_a_state_construction_phase(n_vertices, new_seeds, cl, null_array,
			n_aware, round, tot_edges_explored, tot_spreaders, 3);
		steps++;		
	}		

	//Propagate from initial state
	propagate_from_initial_state(n_vertices, vertices, seed_set, n_aware, round);

	// Update the solution value
	sol_value = seed_set.size(); 
	n_rounds = round;
}