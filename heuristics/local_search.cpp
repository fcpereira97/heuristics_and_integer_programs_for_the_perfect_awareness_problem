/*
File: local_search.cpp
Description: GRASP local search phase for PAP
Author: Felipe de C. Pereira

Based on the article "Effective Heuristics for the Perfect Awareness Problem"
Copyright (C) 2020 Felipe de C. Pereira, Pedro J. de Rezebde, Cid C. de Souza
*/

#include <bits/stdc++.h>
#include "my_struct.h"
#include "propagation.h"

using namespace std;

// Delete from the solution every seed that has the number of neighboring seeds
// more or equal than its threshold
void clean_solution(int n_vertices, vector<Vertex*> &vertices, vector<Vertex*> &seed_set,
	int &sol_value,  int &n_rounds, int &n_aware)
{
	Vertex *seed;
	// For each seed
	for (int i = 0; i < (int)seed_set.size(); i++)
	{
		seed = seed_set[i];

		// If the seed has a number of neighboring seeds more or equal than its threshold
		if(seed-> degree >= 1 && seed-> n_seed_neighs >= seed-> threshold)
		{
			seed-> isSeed = false;
			// Decrease the number of neighboring seeds of its neighbors
			for (list<Vertex*>::iterator neighbor = seed-> neighbors.begin(); neighbor != seed-> neighbors.end(); ++neighbor)
				(*neighbor)-> n_seed_neighs--;

			// Remove the seed from the seed set
	    	seed_set.erase(seed_set.begin() + i);
	    	i--;
	    }
	}

	sol_value = seed_set.size();
}

// Identify and delete exceeding seeds by simulating the propagation
// with a subset of the seed set that grows progressively
void local_search_progressive(int n_vertices, vector<Vertex*> &vertices, vector<Vertex*> &seed_set,
	int &sol_value, int &n_rounds, int &n_aware, string construction_phase_flag,
	double construction_phase_parameter)
{
	int n_aware_aux, round, half_index, last_index; // Control variables

	// ### Begin of the first step ###
	half_index = floor(0.5 * (int)seed_set.size()); // Take the element at the half of seed set

	// Propagate with second half
	erase_propagation(n_vertices, vertices, -1);
	queue<Vertex*> seed_set_aux;
	for (int j = half_index + 1; j < (int)seed_set.size(); j++)
	{
		seed_set_aux.push(seed_set[j]);
		seed_set[j]->state = 2;
	}
	n_aware_aux = seed_set_aux.size();
	round = 0;
	propagate_from_a_state_general(n_vertices, seed_set_aux, n_aware_aux, round);

	// If the removal of a seed v (from the first half)
	// does not affected the feasibility of the solution,
	// then remove v from the seed set
	for (int j = 0; j <= half_index; j++)
	{
		// If v become a spreader due to the propagation
		// or if the subset of seeds is feasible
		if(seed_set[j]->state == 2  || n_aware_aux == n_vertices) 
		{
			seed_set.erase(seed_set.begin() + j);
			j--;
			half_index--;
		}
	}
	// ### End of the first step ###

	// ### Begin of the second step ###
	last_index = half_index; // Last index of first part
	half_index = floor((double)last_index/2); /// Half index of first part

	while(last_index > 0)
	{
		// Propagate with additional seeds from half_index to last_index
		erase_propagation(n_vertices, vertices, -1);
		queue<Vertex*> seed_set_aux_1;
		for (int j = half_index + 1; j < (int)seed_set.size(); j++)
		{
			seed_set_aux_1.push(seed_set[j]);
			seed_set[j]->state = 2;
		}
		n_aware_aux = seed_set_aux_1.size();
		round = 0;
		propagate_from_a_state_general(n_vertices, seed_set_aux_1, n_aware_aux, round);

		// If the removal of a seed v does not affected the feasibility of the solution,
		// then remove v from the seed set
		for (int j = 0; j <= half_index; j++)
		{
			// If v become a spreader due to the propagation
			// or if the subset of seeds is feasible
			if(seed_set[j]->state == 2  || n_aware_aux == n_vertices)
			{
				seed_set.erase(seed_set.begin() + j);
				j--;
				last_index--;
				half_index--;
			}
		}

		// Propagate with additional seeds from first_index to half_index
		erase_propagation(n_vertices, vertices, -1);
		queue<Vertex*> seed_set_aux_2;
		for (int j = 0; j <= half_index; j++)
		{
			seed_set_aux_2.push(seed_set[j]);
			seed_set[j]->state = 2;
		}
		for (int j = last_index + 1; j < (int)seed_set.size(); j++)
		{
			seed_set_aux_2.push(seed_set[j]);
			seed_set[j]->state = 2;
		}
		n_aware_aux = seed_set_aux_2.size();
		round = 0;
		propagate_from_a_state_general(n_vertices, seed_set_aux_2, n_aware_aux, round);

		// If the removal of a seed v does not affected the feasibility of the solution,
		// then remove v from the seed set
		for (int j = half_index + 1; j <= last_index; j++)
		{
			// If v become a spreader due to the propagation
			// or if the subset of seeds is feasible
			if(seed_set[j]->state == 2 || n_aware_aux == n_vertices)
			{
				seed_set.erase(seed_set.begin() + j);
				j--;
				last_index--;
			}
		}

		// Update control variables
		last_index = half_index;
		half_index = floor((double)last_index/2);
	}
	// ### End of the second step ###

	//Propagate from initial state
	propagate_from_initial_state(n_vertices, vertices, seed_set, n_aware_aux, round);

	// Update the solution value
	sol_value = seed_set.size();
	n_rounds = round;
	n_aware = n_aware_aux;
}

// Evaluates, seed by seed, whether it is possible to remove the seed v
// If so, v is deleted from the solution
void local_search_one_by_one(int n_vertices, vector<Vertex*> &vertices, vector<Vertex*> &seed_set,
	int &sol_value, int &n_rounds, int &n_aware)
{
	// Control variables
	Vertex *v1, *v2;
	int n_aware_aux, round;

	for (int i = 0; i < (int)seed_set.size(); i++)
	{	
		v1 = seed_set[i]; // Select a candidate to be removed

		// Erase propagation
		erase_propagation(n_vertices, vertices, -1);

		// Propagate from the initial state without the selected seed
		queue<Vertex*> seed_set_aux;

		for (int j = 0; j < (int)seed_set.size(); j++)
		{
			v2 = seed_set[j];
			if(v1-> id != v2->id)
			{
				seed_set_aux.push(v2);
				v2->state = 2;
			}
		}
		n_aware_aux = seed_set_aux.size();
		round = 0;
		propagate_local_search_one_by_one(n_vertices, seed_set_aux, n_aware_aux, round, v1);

		// If the removal of a seed v does not affected the feasibility of the solution,
		// then remove v from the seed set
		if(n_aware_aux == n_vertices || v1->state == 2)
		{
			seed_set.erase(seed_set.begin() + i);
			i--;
		}
	}

	//Propagate from initial state
	propagate_from_initial_state(n_vertices, vertices, seed_set, n_aware_aux, round);

	// Update the solution value
	sol_value = seed_set.size();
	n_rounds = round;
	n_aware = n_aware_aux;
}

// Evaluates, whether it is possible to remove an element from a subset S' of the seed seed S
// If so, these elements are deleted from S
void local_search_by_sets(int n_vertices, vector<Vertex*> &vertices, vector<Vertex*> &seed_set,
	int &sol_value, int &n_rounds, int &n_aware, double division_parameter)
{
	// Control variables
	int n_aware_aux, round, first_subset_seed, last_subset_seed;
	int subset_size = ceil(division_parameter * seed_set.size());
	first_subset_seed = 0;
	last_subset_seed = subset_size - 1;

	// If the calculated subset size is not equal 0
	if(subset_size != 0)
	{
		// While the first element of the selected subset S' does not exceed the array
		while(first_subset_seed < (int)seed_set.size())
		{
			// Propagate with S - S'
			erase_propagation(n_vertices, vertices, -1);
			queue<Vertex*> seed_set_aux;
			for (int j = 0; j < first_subset_seed; j++)
			{
				seed_set_aux.push(seed_set[j]);
				seed_set[j]->state = 2;
			}
			for (int j = last_subset_seed + 1; j < (int)seed_set.size(); j++)
			{
				seed_set_aux.push(seed_set[j]);
				seed_set[j]->state = 2;
			}
			n_aware_aux = seed_set_aux.size();
			round = 0;
			propagate_from_a_state_general(n_vertices, seed_set_aux, n_aware_aux, round);

			// If the removal of a seed v does not affected the feasibility of the solution,
			// then remove v from the seed set
			for (int j = first_subset_seed; j <= last_subset_seed && j < (int)seed_set.size(); j++)
			{
				if(seed_set[j]->state == 2 || n_aware_aux == n_vertices)
				{
					seed_set.erase(seed_set.begin() + j);
					j--;
					last_subset_seed--;
				}
			}

			// Get next subset
			first_subset_seed = last_subset_seed + 1;
			last_subset_seed = first_subset_seed + subset_size - 1;
		}
	}

	// Propagate from initial state
	propagate_from_initial_state(n_vertices, vertices, seed_set, n_aware_aux, round);

	// Update the solution value
	sol_value = seed_set.size();
	n_rounds = round;
	n_aware = n_aware_aux;
}
