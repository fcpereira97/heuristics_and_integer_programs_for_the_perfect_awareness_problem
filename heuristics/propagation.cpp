/*
File: propagation.cpp
Description: Propagation process of PAP used by both construction and local search phases of GRASP
Author: Felipe de C. Pereira

Based on the article "Effective Heuristics for the Perfect Awareness Problem"
Copyright (C) 2020 Felipe de C. Pereira, Pedro J. de Rezebde, Cid C. de Souza
*/

#include <bits/stdc++.h>
#include "my_struct.h"
#include "priority_queue.h"

using namespace std;

// Erase vertices status related to the propagation process
void erase_propagation(int n_vertices, vector<Vertex*> &vertices, int construction_phase_flag)
{
	// For each vertex v
	for (int i = 0; i < n_vertices; ++i)
	{
		// Erase status data of v
		vertices[i]-> n_spreader_neighs = 0;
		vertices[i]-> n_seed_neighs = 0;
		vertices[i]-> n_unaware_neighs = vertices[i]->degree;
		vertices[i]-> n_spreader_neighs_gain = 0;
		vertices[i]-> lack_threshold = vertices[i]->threshold;
		vertices[i]-> state = 0;
		vertices[i]-> isSeed = false;
		vertices[i]-> id_area = 0;
	}

	// If the method was called in construction phase with strategie WGR, RG or SG
	if(construction_phase_flag >= 2)
	{
		// For each vertex v with threshold 1 initialize the  
		// gain of neighboring spreaders of each neighbors of v
		for (int i = 0; i < n_vertices; ++i)
			if(vertices[i]-> threshold == 1)
				for (list<Vertex*>::iterator neigh = vertices[i]-> neighbors.begin(); neigh != vertices[i]-> neighbors.end(); ++neigh)
					(*neigh)-> n_spreader_neighs_gain++;
	}
}

// Move the vertex to the correct position at CL keeping the CL sorted
// The vertex is moved to the range of vertices in CL with a lower benefit
// Used in construction phase stategies GR and WGR
void move_in_cl_GR_WGR(Vertex* vertex, vector<Vertex*> &cl, vector<int> &cl_aux)
{
	// Swap the vertex in CL with the righest element
	// of the CL with same benefit value of the vertex
	cl[vertex-> cl_pos] = cl[cl_aux[vertex-> n_unaware_neighs]];
	cl[vertex-> cl_pos]-> cl_pos = vertex->cl_pos;
	cl[cl_aux[vertex-> n_unaware_neighs]] = vertex;
	vertex-> cl_pos = cl_aux[vertex-> n_unaware_neighs];

	// Update the index of the righest element with same benefit
	cl_aux[vertex-> n_unaware_neighs]--;
}

// Remove a vertex from CL, keeping cl sorted
// Used in construction phase stategies GR and WGR
void delete_from_cl_GR_WGR(Vertex *vertex, vector<Vertex*> &cl, vector<int> &cl_aux)
{
	// Decrease the benefit of the vertex and
	// move it in CL sequentially to the last position
	while(vertex-> n_unaware_neighs >= 0)
	{
		move_in_cl_GR_WGR(vertex, cl, cl_aux);
		vertex-> n_unaware_neighs--;
	}
	cl.erase(cl.end()-1); // Delete the vertex from CL
}

// Remove a vertex from CL
// Used in construction phase stategies SG
void delete_from_cl_SG(Vertex* vertex, vector<Vertex*> &cl)
{
	// Swap the vertex in CL with the one at the last position
	cl[vertex-> cl_pos] = cl[cl.size() - 1];
	cl[vertex-> cl_pos]-> cl_pos = vertex->cl_pos;
	cl[cl.size() - 1] = vertex;
	vertex-> cl_pos = cl.size() - 1;

	cl.pop_back(); // Delete the vertex from CL
}

// For a new aware vertex v, decrease the number of unaware neighbors
// of the neighbors of v
void decrease_n_unaware_neighs(Vertex *new_aware_vertex, vector<Vertex*> &cl,
	vector<int> &cl_aux, int construction_phase_flag)
{
	// For each neighbor u of the new aware vertex v
	for (list<Vertex*>::iterator neighbor = new_aware_vertex-> neighbors.begin(); neighbor != new_aware_vertex-> neighbors.end(); ++neighbor)
	{
		// If u is not a spreader, then move it to the correct position at CL

		if((*neighbor)-> state != 2 && construction_phase_flag < 3)	// For GR or WGR			
			move_in_cl_GR_WGR((*neighbor), cl, cl_aux);

		(*neighbor)-> n_unaware_neighs--; // Update the number of unaware neighbors

		if((*neighbor)-> state != 2 && construction_phase_flag == 3) // For RG
			heapify(cl, (*neighbor)-> cl_pos, comp_pipeline);
	}
}

// For a new spreader vertex v, decrease the gain
// of neighboring spreaders of each neighbors of v
void decrease_n_spreader_neighs_gain(Vertex* new_spreader_vertex, vector<Vertex*> &cl,
	int construction_phase_flag)
{
	// For each neighbor u of new spreader vertex v
	for (list<Vertex*>::iterator neighbor = new_spreader_vertex-> neighbors.begin(); neighbor != new_spreader_vertex-> neighbors.end(); ++neighbor)
	{
		(*neighbor)-> n_spreader_neighs_gain--; // Decrease the gain of neighboring spreaders of u
		
		// In case of RG strategie, fix the priority queue
		if((*neighbor)-> state != 2 && construction_phase_flag == 3)
			heapify(cl, (*neighbor)-> cl_pos, comp_pipeline);
	}
}

// For a new almost spreader vertex v, increase the gain
// of neighboring spreaders of each neighbors of v
void increase_n_spreader_neighs_gain(Vertex* new_almost_spreader_vertex, vector<Vertex*> &cl,
	int construction_phase_flag)
{
	// For each neighbor u of new spreader almost vertex v
	for (list<Vertex*>::iterator neighbor = new_almost_spreader_vertex-> neighbors.begin(); neighbor != new_almost_spreader_vertex-> neighbors.end(); ++neighbor)
	{
		(*neighbor)-> n_spreader_neighs_gain++; // Increase the gain of neighboring spreaders of u
		
		// In case of RG strategie, fix the priority queue
		if((*neighbor)-> state != 2 && construction_phase_flag == 3)
			heapify(cl, (*neighbor)-> cl_pos, comp_pipeline);
	}
}

// Continue the propagation process from the last propagation state using new spreaders 
// Used in construction phase
void propagate_from_a_state_construction_phase(int n_vertices, queue<Vertex*> &next_spreaders,
	vector<Vertex*> &cl, vector<int> &cl_aux, int &n_aware, int &round, int &tot_edges_explored,
	int &tot_spreaders, int construction_phase_flag)
{
	// Number of spreaders that will propagate at the current round
	int n_next_spreaders;

	// While exists a vertex that is still unaware and
	// there is at least one spreader that will propagate at the current round
	while(n_aware < n_vertices && !next_spreaders.empty())
	{
		n_next_spreaders = next_spreaders.size();
		round++; // Count a new round
		
		// For each spreader that will propagate at this round
		while(n_next_spreaders > 0)
		{
			// Get the next spreader v
			Vertex *spreader = next_spreaders.front();
			next_spreaders.pop();

			// Update statistics
			tot_edges_explored += spreader -> degree + 1;
			tot_spreaders++;
			
			// For each neighnbor u of v
			for (list<Vertex*>::iterator neighbor = spreader-> neighbors.begin(); neighbor != spreader-> neighbors.end(); ++neighbor)
			{
				// Update the status of u
				(*neighbor)-> n_spreader_neighs++;
				(*neighbor)-> lack_threshold--;			

				// In case of RG, fix the priority queue
				if((*neighbor)-> state != 2 && construction_phase_flag == 3) 
					heapify(cl, (*neighbor)-> cl_pos, comp_pipeline);

    			// If v is actually a seed (that was selected) then the
    			// number of neighboring seeds of the u is increased
    			if(spreader-> isSeed)
    				(*neighbor)-> n_seed_neighs++;

    			// If the u is unaware then it becomes aware
    			if((*neighbor)-> state == 0)
    			{
    				(*neighbor)-> state = 1;
    				n_aware++;
    				decrease_n_unaware_neighs((*neighbor), cl, cl_aux, construction_phase_flag);
    			}

    			// If u is almost spreader, then updates the status of the neighbors of u
    			if(construction_phase_flag >= 2 && (*neighbor)-> lack_threshold == 1)
    				increase_n_spreader_neighs_gain(*neighbor, cl, construction_phase_flag);
    			
	    		// If u is aware and its threshold was reached, it becomes a spreader
	    		if((*neighbor)-> state == 1 && (*neighbor)-> n_spreader_neighs >= (*neighbor)-> threshold)
	    		{
	    			// Turn u into a spreader
	    			(*neighbor)-> state = 2;
	    			(*neighbor)-> lack_threshold = 0;
	    			if(construction_phase_flag >= 2)
	    				decrease_n_spreader_neighs_gain(*neighbor, cl, construction_phase_flag);

	    			// Insert u in the list of vertices that will propagate at the next round
	    			next_spreaders.push((*neighbor)); 

	    			// Remove u from CL
	    			if(construction_phase_flag < 3)
	    				delete_from_cl_GR_WGR(*neighbor, cl, cl_aux); 
	    			else if(construction_phase_flag == 3)
	    				delete_node_from_heap(cl, (*neighbor)-> cl_pos, comp_pipeline);
	    			else if(construction_phase_flag == 4)
	    				delete_from_cl_SG(*neighbor, cl);
	    		}
			}
			n_next_spreaders--;
		}
	}
}

// Continue the propagation process from the last propagation state using new spreaders 
// For general purposes
void propagate_from_a_state_general(int n_vertices, queue<Vertex*> &next_spreaders, int &n_aware,
	int &round)
{
	// Number of spreaders that will propagate at the current round
	int n_next_spreaders;

	// While exists a vertex that is still unaware and
	// there is at least one spreader that will propagate at the current round
	while(n_aware < n_vertices && !next_spreaders.empty())
	{
		n_next_spreaders = next_spreaders.size();
		round++; // Count a new round
		
		// For each spreader that will propagate at this round
		while(n_next_spreaders > 0)
		{
			// Get the next spreader v
			Vertex * spreader = next_spreaders.front();
			next_spreaders.pop();

			// For each neighnbor u of v
			for (list<Vertex*>::iterator neigh = spreader-> neighbors.begin(); neigh != spreader-> neighbors.end(); ++neigh)
			{
				// Update the status of u
				(*neigh)-> n_spreader_neighs++;

    			// If the spreader is a seed, the number of seed neighs of the neigh is increased
    			if(spreader-> isSeed)
    				(*neigh)-> n_seed_neighs++;

    			// If the u is unaware then it becomes aware
    			if((*neigh)-> state == 0)
    			{
    				(*neigh)-> state = 1;
    				n_aware++;
    			}

	    		// If u is aware and its threshold was reached, it becomes a spreader
	    		if((*neigh)-> state == 1 && (*neigh)-> n_spreader_neighs >= (*neigh)-> threshold)
	    		{	
	    			// Turn u into a spreader
	    			// Insert u in the list of vertices that will propagate at the next round
	    			(*neigh)-> state = 2;
	    			next_spreaders.push((*neigh));
	    		}
			}
			n_next_spreaders--;
		}
	}
}

// Simulate a complete propagation process, starting from round 0
void propagate_from_initial_state(int n_vertices, vector<Vertex*> &vertices, vector<Vertex*> &seed_set,
	int &n_aware, int &round)
{
	// Erase the propagation
	erase_propagation(n_vertices, vertices, -1);

	// Configure an auxiliar seed set and initial state of the seeds
	queue<Vertex*> seed_set_aux;
	for (int i = 0; i < (int)seed_set.size(); i++)
	{
		seed_set[i]-> state = 2;
		seed_set[i]-> isSeed = true;
		seed_set[i]-> lack_threshold = 0;
		seed_set_aux.push(seed_set[i]);
	}

	// Use the auxiliar seed set to start the propagation from round 0
	n_aware = seed_set.size();
	round = 0;
	propagate_from_a_state_general(n_vertices, seed_set_aux, n_aware, round);
}

// Turn into seeds all vertices that has degree 0
void initialize_natural_seeds(int n_vertices, vector<Vertex*> &vertices, vector<Vertex*> &seed_set,
	int &n_aware)
{
	// For each vertex v
	for(int i = 0; i < n_vertices; i++)
	{
		if(vertices[i]-> degree == 0) // If v has degree 0
		{
			// Configure the initial state of v as seed and insert v into the seed set
			vertices[i]-> state = 2;
			vertices[i]-> isSeed = true;
			vertices[i]-> lack_threshold = 0;
			seed_set.push_back(vertices[i]);
			n_aware++;
		}
	}
}

// Continue the propagation process from the last propagation state using the
// original seed set except for a seed that was removed
// Used in local search
void propagate_local_search_one_by_one(int n_vertices, queue<Vertex*> &next_spreaders, int &n_aware,
	int &round, Vertex *removed_seed)
{
	// Number of spreaders that will propagate at the current round
	int n_next_spreaders;

	// While exists a vertex that is still unaware and
	// there is at least one spreader that will propagate at the current round and
	// the removed seed did not become a spreader
	while(n_aware < n_vertices && !next_spreaders.empty() && removed_seed-> state < 2)
	{
		n_next_spreaders = next_spreaders.size();
		round++; // Count a new round
		
		// For each spreader that will propagate at this round
		while(n_next_spreaders > 0)
		{
			// Get the next spreader v
			Vertex * spreader = next_spreaders.front();
			next_spreaders.pop();

			// For each neighnbor u of v
			for (list<Vertex*>::iterator neigh = spreader-> neighbors.begin(); neigh != spreader-> neighbors.end(); ++neigh)
			{
				// Update the status of u
				(*neigh)-> n_spreader_neighs++;

    			// If the spreader is a seed, the number of seed neighs of the neigh is increased
    			if(spreader-> isSeed)
    				(*neigh)-> n_seed_neighs++;

    			// If the u is unaware then it becomes aware
    			if((*neigh)-> state == 0)
    			{
    				(*neigh)-> state = 1;
    				n_aware++;
    			}

	    		// If u is aware and its threshold was reached, it becomes a spreader
	    		if((*neigh)-> state == 1 && (*neigh)-> n_spreader_neighs >= (*neigh)-> threshold)
	    		{
	    			// Turn u into a spreader
	    			// Insert u in the list of vertices that will propagate at the next round
	    			(*neigh)-> state = 2;
	    			next_spreaders.push((*neigh));
	    		}
			}
			n_next_spreaders--;
		}
	}
}