/*
File: propagation.h
Description: Header file of propagation.cpp
Author: Felipe de C. Pereira

Based on the article "Effective Heuristics for the Perfect Awareness Problem"
Copyright (C) 2020 Felipe de C. Pereira, Pedro J. de Rezebde, Cid C. de Souza
*/

#ifndef PROPAGATION_H
#define PROPAGATION_H

// Erase vertices status related to the propagation process
void erase_propagation(int n_vertices, vector<Vertex*> &vertices, int construction_phase_flag);

// Move the vertex to the correct position at CL keeping the CL sorted
// The vertex is moved to the range of vertices in CL with a lower benefit
// Used in construction phase stategies GR and WGR
void move_in_cl_GR_WGR(Vertex* vertex, vector<Vertex*> &cl, vector<int> &cl_aux);

// Remove a vertex from CL, keeping cl sorted
// Used in construction phase stategies GR and WGR
void delete_from_cl_GR_WGR(Vertex* vertex, vector<Vertex*> &cl, vector<int> &cl_aux);

// Remove a vertex from CL
// Used in construction phase stategies SG
void delete_from_cl_SG(Vertex* vertex, vector<Vertex*> &cl);

// For a new aware vertex v, decrease the number of unaware neighbors
// of the neighbors of v
void decrease_n_unaware_neighs(Vertex* new_aware_vertex, vector<Vertex*> &cl, vector<int> &cl_aux,
 int construction_phase_flag);

// For a new spreader vertex v, decrease the gain
// of neighboring spreaders of each neighbors of v
void decrease_n_spreader_neighs_gain(Vertex* new_spreader_vertex, vector<Vertex*> &cl,
	int construction_phase_flag);

// For a new almost spreader vertex v, increase the gain
// of neighboring spreaders of each neighbors of v
void increase_n_spreader_neighs_gain(Vertex* new_almost_spreader_vertex, vector<Vertex*> &cl,
	int construction_phase_flag);

// Continue the propagation process from the last propagation state using new spreaders 
// Used in construction phase
void propagate_from_a_state_construction_phase(int n_vertices, queue<Vertex*> &next_spreaders,
	vector<Vertex*> &cl, vector<int> &cl_aux, int &n_aware, int &round, int &tot_edges_explored,
	int &tot_spreaders, int construction_phase_flag);

// Continue the propagation process from the last propagation state using new spreaders 
// For general purposes
void propagate_from_a_state_general(int n_vertices, queue<Vertex*> &next_spreaders, int &n_aware,
	int &round);

// Simulate a complete propagation process, starting from round 0
void propagate_from_initial_state(int n_vertices, vector<Vertex*> &vertices,
	vector<Vertex*> &seed_set, int &n_aware, int &round);

// Turn into seeds all vertices that has degree 0
void initialize_natural_seeds(int n_vertices, vector<Vertex*> &vertices,
	vector<Vertex*> &seed_set,int &n_aware);

// Continue the propagation process from the last propagation state using the
// original seed set except for a seed that was removed
// Used in local search
void propagate_local_search_one_by_one(int n_vertices, queue<Vertex*> &next_spreaders,
	int &n_aware, int &round, Vertex *removed_seed);
	
#endif