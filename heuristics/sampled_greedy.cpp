/*
File: sampled_greedy.cpp
Description: GRASP construction phase for PAP with Sampled Greedy strategie
Author: Felipe de C. Pereira

Based on the article "Effective Heuristics for the Perfect Awareness Problem"
Copyright (C) 2020 Felipe de C. Pereira, Pedro J. de Rezebde, Cid C. de Souza
*/

#include <bits/stdc++.h>
#include "my_struct.h"
#include "propagation.h"

using namespace std;

// GRASP construction phase with Sampled Greedy strategie
void sampled_greedy_construction(int n_vertices, vector<Vertex*> &vertices, vector<Vertex*> &seed_set,
	int &sol_value,	int &n_rounds, int &n_aware, int n_seeds_per_insertion, double sample_size_rate,
	int &draws_altruistic, int &draws_philanthropic, int &steps, int &tot_spreaders, int &tot_edges_explored,
	double &rlc_cl_ratio_mean, mt19937 &rng)
{
	int round, sample_size; // Control variables
	queue<Vertex*> new_seeds; // New selected seeds
	vector<Vertex*> cl(n_vertices); // Candidate list (CL)
	sample_size = ceil(sample_size_rate * n_vertices); // Size of random sample
	vector<int> null_array;

	// Initialize variables
	n_aware = 0;
	round = 0;

	// Erase the propagation process
	erase_propagation(n_vertices, vertices, 4);

	// Turn to seeds all vertices that has degree equals to 0
	initialize_natural_seeds(n_vertices, vertices, seed_set, n_aware);

	// Build CL
	cl.assign(vertices.begin(), vertices.end());
	for(int i = 0; i < int(cl.size()); i++)
		cl[i]-> cl_pos = i;

	// While the seed set is still a infeasible solutions do
	while(n_aware < n_vertices)
	{	
		// For as many seeds per selection as wanted
		for(int i = 0; i < n_seeds_per_insertion; i++)
		{
			// Selected vertice variable
			Vertex *v_chosen;
			int v_chosen_index;

			// Tie counters
			int partial_draws_philanthropic = 0;
			int partial_draws_altruistic = 0;
			
			// Sample the first element from cl
			v_chosen_index = rng() % cl.size();
			v_chosen = cl[v_chosen_index];
			swap(cl[v_chosen_index], cl[cl.size()-1]);
			swap(cl[v_chosen_index]-> cl_pos, cl[cl.size()-1]-> cl_pos);

			// Sample sample_size elements from cl and select the one with greater benefit
			for(int i = 1; i < min(sample_size, int(cl.size())); i++)
			{
				// Sample another element from cl
				int j = rng() % (cl.size() - i);
				swap(cl[j], cl[cl.size()-1-i]);
				swap(cl[j]->cl_pos, cl[cl.size()-1-i]-> cl_pos);

				// Apply tiebreaker criteria
				if(cl[j]-> n_unaware_neighs > v_chosen-> n_unaware_neighs)
				{
					partial_draws_altruistic = 0;
					partial_draws_philanthropic = 0;
					v_chosen = cl[j];
				}
				else if(cl[j]-> n_unaware_neighs == v_chosen-> n_unaware_neighs 
					&& cl[j]-> n_spreader_neighs_gain > v_chosen-> n_spreader_neighs_gain)
				{
					partial_draws_altruistic = 1;
					partial_draws_philanthropic = 0;
					v_chosen = cl[j];
				}
				else if(cl[j]-> n_unaware_neighs == v_chosen-> n_unaware_neighs 
					&& cl[j]-> n_spreader_neighs_gain == v_chosen-> n_spreader_neighs_gain 
					&& cl[j]-> lack_threshold < v_chosen-> lack_threshold)
				{
					partial_draws_altruistic = 1;
					partial_draws_philanthropic = 1;
					v_chosen = cl[j];
				}
			}

			// Update statistics about ties
			draws_altruistic += partial_draws_altruistic;
			draws_philanthropic += partial_draws_philanthropic;
			
			// Save benefits to use in local search phase
			v_chosen-> n_unaware_neighs_old = v_chosen-> n_unaware_neighs;
			v_chosen-> n_spreader_neighs_gain_old = v_chosen-> n_spreader_neighs_gain;
			v_chosen-> lack_threshold_old = v_chosen-> lack_threshold;

			// If the vertex selected is not spreader
			if(v_chosen-> state != 2)
			{
				// If the vertex selected is unaware
				if(v_chosen-> state == 0)
				{
					n_aware++; // Update total number of aware vertices
					
					// Update the number of unaware neighs of the selected vertex neighbors
					decrease_n_unaware_neighs(v_chosen, cl, null_array, 4);
				}

				// Turn selected vertex into a seed
				// Update the cl and the benefits
				v_chosen-> state = 2;
				v_chosen-> isSeed = true;
				if(v_chosen -> lack_threshold == 1)
					decrease_n_spreader_neighs_gain(v_chosen, cl, 4);
				v_chosen-> lack_threshold = 0;
				seed_set.push_back(v_chosen);
				new_seeds.push(v_chosen);
				delete_from_cl_SG(v_chosen, cl);
			}
		}

		// Continue the propagation with the new seeds 

		propagate_from_a_state_construction_phase(n_vertices, new_seeds, cl, null_array,
			n_aware, round, tot_edges_explored, tot_spreaders, 4);
		steps++;
	}		

	//Propagate from initial state
	propagate_from_initial_state(n_vertices, vertices, seed_set, n_aware, round);

	// Update the solution value
	sol_value = seed_set.size(); 
	n_rounds = round;
}