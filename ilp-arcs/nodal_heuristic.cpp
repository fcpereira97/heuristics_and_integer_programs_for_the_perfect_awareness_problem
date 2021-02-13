#include <bits/stdc++.h>
#include "my_struct.h"

using namespace std;

void nodal_erase_propagation(int n_vertices, vector<Vertex*> *vertices)
{
	for (int i = 0; i < n_vertices; ++i)
	{
		(*vertices)[i]-> n_spreader_neighs = 0;
		(*vertices)[i]-> state = 0;
		(*vertices)[i]-> isSeed = false;
		(*vertices)[i]-> round_became_spreader = -1;
	}
}

// Starting from a snapshot of the propagation state, this method simulates continue the propagation with new spreaders
void nodal_propagate_from_a_state(int n_vertices, queue<Vertex*> * next_spreaders, int *n_aware, int *round)
{
	// Spreaders that will possibly nodal_propagate_from_a_state at this round
	int n_next_spreaders;

	// While exists a vertex that is still unaware or while at least a new vertex become spreader at the round, do
	while(!(*next_spreaders).empty())
	{
		// The vertices that become spreaders in last round are considered the current ones to nodal_propagate_from_a_state
		n_next_spreaders = (*next_spreaders).size();
		(*round)++; // Count a new round
		
		// For each spreader that nodal_propagate_from_a_states at this round, do
		while(n_next_spreaders > 0)
		{
			// Get the next vertex
			Vertex * spreader = (*next_spreaders).front();
			(*next_spreaders).pop();

			// For each neigh of the spreader do
			for (list<Vertex*>::iterator neigh = spreader-> neighbors.begin(); neigh != spreader-> neighbors.end(); ++neigh)
			{
				// Number of spreaders neighs of the neigh is increased
				(*neigh)-> n_spreader_neighs++;

    			// If the neigh is unaware, it becomes aware
    			if((*neigh)-> state == 0)
    			{
    				(*neigh)-> state = 1;
    				(*n_aware)++; // Count a new aware vertex
    			}

	    		// If the neigh is aware and its threshold was reached, it becomes a spreader
	    		if((*neigh)-> state == 1 && (*neigh)-> n_spreader_neighs >= (*neigh)-> threshold)
	    		{
	    			(*neigh)-> state = 2;
	    			(*neigh)-> round_became_spreader = (*round);
	    			(*next_spreaders).push((*neigh)); // Insert the neigh in the list of vertex to nodal_propagate_from_a_state at the next round
	    		}
			}
			n_next_spreaders--;
		}
	}
}

// Starting from round 0 of the propagation state, this method simulates all the propagation
void nodal_propagate_from_initial_state(int n_vertices, vector<Vertex*> *vertices, vector<Vertex*> *seed_set, int *n_aware, int *round)
{
	// Erase the propagation
	nodal_erase_propagation(n_vertices, vertices);

	// Configure a auxiliar seed set
	queue<Vertex*> seed_set_aux;
	for (int i = 0; i < (int)(*seed_set).size(); i++)
	{
		(*seed_set)[i]-> state = 2;
		(*seed_set)[i]-> isSeed = true;
		(*seed_set)[i]-> round_became_spreader = 0;
		seed_set_aux.push((*seed_set)[i]);
	}

	// Use the auxiliar seed set to nodal_propagate_from_a_state from the initial state of the graph
	(*n_aware) = (*seed_set).size();
	(*round) = 0;
	nodal_propagate_from_a_state(n_vertices, &seed_set_aux, n_aware, round);
}

// This local search deletes from the solution all seed that if removed won't turn the solution infeasible
void nodal_local_search(int n_vertices, vector<Vertex*> *vertices, vector<Vertex*> * seed_set, int *sol_value)
{
	Vertex *v1, *v2;
	int n_aware, round, iterations_limit;
	iterations_limit = (*seed_set).size();

	for (int i = 0; i < iterations_limit; i++)
	{	
		v1 = (*seed_set)[i]; // Select a candidate to be removed

		// Erase propagation
		nodal_erase_propagation(n_vertices, vertices);

		// nodal_propagate_from_a_state from the initial state of the graph without the selected candidate within the seeders
		queue<Vertex*> seed_set_aux;

		for (int j = 0; j < (int)(*seed_set).size(); j++)
		{
			v2 = (*seed_set)[j];
			if(v1-> id != v2->id)
			{
				seed_set_aux.push(v2);
				v2->state = 2;
			}
		}
		n_aware = seed_set_aux.size();
		round = 0;
		nodal_propagate_from_a_state(n_vertices, &seed_set_aux, &n_aware, &round);

		// If the removal does not affected the feasibility of the solution, remove the seed from the seed set
		if(n_aware == n_vertices)
		{
			(*seed_set).erase((*seed_set).begin() + i);
			i--;
			iterations_limit--;
		}
	}

	//Propagate from initial state
	nodal_propagate_from_initial_state(n_vertices, vertices, seed_set, &n_aware, &round);

	// Update the solution value
	*sol_value = (*seed_set).size();
}

void nodal_heuristic(int n_vertices, vector<Vertex*> *vertices, int *sol_value)
{
	vector<Vertex*> seed_set;
	vector<Vertex*> vertices_aux(n_vertices);
	vertices_aux.assign((*vertices).begin(), (*vertices).end());

	sort(vertices_aux.begin(), vertices_aux.end(), compare_two_vertices_by_var_value);
	
	int next_seed, round, n_aware;
	next_seed = round = n_aware = 0;

	queue<Vertex*> new_seeds; // Next vertices to nodal_propagate_from_a_state

	nodal_erase_propagation(n_vertices, vertices);	// Erase the propagation scheme

	// While the seed set is still a infeasible solutions, do
	while(n_aware < n_vertices)
	{	
		Vertex *v_chosen = vertices_aux[next_seed];
		next_seed++;

		// If the vertex selected is not spreader
		if(v_chosen-> state != 2)
		{
			// If the vertex selected is unaware
			if(v_chosen-> state == 0)
				n_aware++; // Update total number of aware vertices

			// Turn selected vertex into a seed
			v_chosen-> state = 2;
			v_chosen-> isSeed = true;
			seed_set.push_back(v_chosen);
			new_seeds.push(v_chosen);

			// Continue the propagation with the new seeds 
			nodal_propagate_from_a_state(n_vertices, &new_seeds, &n_aware, &round);
		}
	}		

	//Propagate from initial state
	nodal_propagate_from_initial_state(n_vertices, vertices, &seed_set, &n_aware, &round);

	// Update the solution value
	*sol_value = (seed_set).size();
	//cout << "SOL_VALUE = " << *sol_value << endl;

	nodal_local_search(n_vertices, vertices, &seed_set, sol_value);
}

