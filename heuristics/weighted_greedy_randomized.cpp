/*
File: weighted_greedy_randomized.cpp
Description: GRASP construction phase for PAP with Weighted Greedy Randomized strategie
Author: Felipe de C. Pereira

Based on the article "Effective Heuristics for the Perfect Awareness Problem"
Copyright (C) 2020 Felipe de C. Pereira, Pedro J. de Rezebde, Cid C. de Souza
*/

#include <bits/stdc++.h>
#include "my_struct.h"
#include "propagation.h"

using namespace std;

// GRASP construction phase with Weighted Greedy Randomized strategie
void weighted_greedy_randomized_construction(int n_vertices, vector<Vertex*> &vertices,
	vector<Vertex*> &seed_set, int &sol_value, 	int &n_rounds, int &n_aware,
	int n_seeds_per_insertion, 	double alpha, int &steps, int &tot_spreaders,
	int &tot_edges_explored, double &rlc_cl_ratio_mean, mt19937 &rng)
{
	int round, rcl_size, min_benefit, max_benefit, max_degree; // Control variables
	queue<Vertex*> new_seeds; // New selected seeds
	vector<Vertex*> cl(n_vertices); // Candidate list (CL)

	// Auxiliar array used to keep CL sorted
	// cl_aux[x] = y indicates that the last element in cl with benefit x is at position y
	vector<int> cl_aux;

	// Initialize variables
	n_aware = 0;
	round = 0;
	max_degree = 0;

	// Erase the propagation process
	erase_propagation(n_vertices, vertices, 2);

	// Turn to seeds all vertices that has degree equals to 0
	initialize_natural_seeds(n_vertices, vertices, seed_set, n_aware);

	// Build CL and sort its elements by number unware neighbors
	cl.assign(vertices.begin(), vertices.end());
	sort(cl.begin(), cl.end(), compare_two_vertices_by_unaware_neighs);

	// Find maximum degree in CL elements
	for(int i = 0; i < int(cl.size()); i++)
		if(vertices[i]-> degree > max_degree)
			max_degree = vertices[i]-> degree;

	// Configure cl_aux
	cl_aux.resize(max_degree+1);
	fill(cl_aux.begin(), cl_aux.end(), -1);
	for(int i = 0; i < int(cl.size()) - 1; i++)
	{
		cl[i]-> cl_pos = i;
		if(cl[i]-> n_unaware_neighs != cl[i+1]-> n_unaware_neighs)
			cl_aux[cl[i]-> n_unaware_neighs] = i;
	}
	cl[cl.size()-1]-> cl_pos = cl.size()-1;
	cl_aux[cl[cl.size()-1]-> n_unaware_neighs] = cl.size()-1;

	// Fix cl_aux for benefit values that do not correspond with any vertex
	for(int i = cl_aux.size() - 2; i >= 0; i--)
		if(cl_aux[i] == -1)
			cl_aux[i] = cl_aux[i+1];

	// While the seed set is still a infeasible solution do
	while(n_aware < n_vertices)
	{	
		// Using alpha parameter, calculates the maximum and the minimum benefit values
		max_benefit = cl[0]-> n_unaware_neighs;
		min_benefit = max_benefit - int(alpha * (max_benefit - cl[cl.size() - 1] -> n_unaware_neighs));

		// Determine the range of the RCL
		rcl_size = cl_aux[min_benefit] + 1;
		rlc_cl_ratio_mean += double(rcl_size) / double(cl.size());

		// For as many seeds per selection as wanted
		for(int i = 0; i < n_seeds_per_insertion; i++)
		{	
			// For each vertex in RCL we calculate the weight for draw
			vector<int> weights(rcl_size);
			for(int j = 0; j < rcl_size; j++)
				weights[j] = vertices[j]-> n_spreader_neighs_gain;

			// Normalize the weights
			discrete_distribution<int> distribution(weights.begin(), weights.end());

			// Draw a vertex from rcl using the weight of a vertex v as the probability of v be selected
			Vertex *v_chosen = cl[distribution(rng)];

			// If the selected vertex is not spreader
			if(v_chosen-> state != 2)
			{
				// If the selected vertex is unaware
				if(v_chosen-> state == 0)
				{
					n_aware++; // Update the total number of aware vertices
					
					// Update the number of unaware neighs of the selected vertex neighbors
					decrease_n_unaware_neighs(v_chosen, cl, cl_aux, 2);
				}

				// Turn selected vertex into a seed
				// Update the cl and the benefits
				v_chosen-> state = 2;
				v_chosen-> isSeed = true;
				if(v_chosen -> lack_threshold == 1)
					decrease_n_spreader_neighs_gain(v_chosen, cl, 2);
				v_chosen-> lack_threshold = 0;
				seed_set.push_back(v_chosen);
				new_seeds.push(v_chosen);
				delete_from_cl_GR_WGR(v_chosen, cl, cl_aux);
			}
		}

		// Continue the propagation with the new seeds 
		propagate_from_a_state_construction_phase(n_vertices, new_seeds, cl, cl_aux,
			n_aware, round, tot_edges_explored, tot_spreaders, 2);
		steps++;
	}		

	//Propagate from initial state
	propagate_from_initial_state(n_vertices, vertices, seed_set, n_aware, round);

	// Update the solution value
	sol_value = seed_set.size(); 
	n_rounds = round;
}