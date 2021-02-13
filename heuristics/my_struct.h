#ifndef MYSTRUCT_H
#define MYSTRUCT_H

#include <bits/stdc++.h>

using namespace std;

/* Representation of a vertex
* id - number assinged to the vertex during the graph loading
* threshold - threshold assinged to the vertex
* degree - number of edges that contain the vertex
* n_spreader_neighs - number of neighbors of the vertex that are spreaders
* n_seed_neighs - number of neighbors of the vertex that are seeders
* n_aware_neighs - number of neighbors of the vertex that are aware
* n_unaware_neighs - number of neighbors of the vertex that are unaware
* state - 0 is unaware, 1 is aware and 2 is spreader
* isSeed - true if it is a seed, else false
* neighbors - list of references to neighboring vertices
*/
struct Vertex
{
	int id;
	int threshold;
	int degree;
	int n_spreader_neighs;
	int n_seed_neighs;
	int state;
	int group_id;
	int community_id;
	bool isSeed;
	int n_spreader_neighs_gain;
	int n_unaware_neighs;
	int lack_threshold;
	int n_spreader_neighs_gain_old;
	int n_unaware_neighs_old;
	int lack_threshold_old;
	int cl_pos;
	int ls_benefit;
	int id_area;

	list<int> ids_before_preprocessing;
	list<Vertex*> neighbors;
	list<Vertex*> neighbors_backup;
};

typedef struct Vertex Vertex;

// Compare vertices by number of unaware neighbors
bool compare_two_vertices_by_unaware_neighs(Vertex* a, Vertex* b);

bool compare_two_vertices_by_lack_threshold(Vertex* a, Vertex* b);

bool compare_two_vertices_by_spreader_neighs_gain(Vertex* a, Vertex* b);

bool compare_two_vertices_by_state(Vertex* a, Vertex* b);

bool compare_two_vertices_by_ls_gr_criteria(Vertex* a, Vertex* b);

// Return the vertex with best benefit considering the pipeline
bool comp_pipeline(Vertex *v1, Vertex *v2);

bool comp_pipeline_by_ls_sg_criteria(Vertex *v1, Vertex *v2);

#endif