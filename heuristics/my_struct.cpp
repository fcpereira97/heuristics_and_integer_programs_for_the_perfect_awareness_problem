#include <bits/stdc++.h>
#include "my_struct.h"

using namespace std;

// Compare vertices by number of unaware neighbors
bool compare_two_vertices_by_unaware_neighs(Vertex* a, Vertex* b) 
{ 
    return a-> n_unaware_neighs > b-> n_unaware_neighs;
}

// Compare vertices by number of unaware neighbors
bool compare_two_vertices_by_lack_threshold(Vertex* a, Vertex* b) 
{ 
    return a-> lack_threshold > b-> lack_threshold;
}

// Compare vertices by number of unaware neighbors
bool compare_two_vertices_by_spreader_neighs_gain(Vertex* a, Vertex* b) 
{ 
    return a-> n_spreader_neighs_gain > b-> n_spreader_neighs_gain;
}

// Compare vertices by number of state neighbors
bool compare_two_vertices_by_state(Vertex* a, Vertex* b) 
{ 
    return a-> state < b-> state;
}

bool compare_two_vertices_by_ls_gr_criteria(Vertex* a, Vertex* b)
{
	return a-> ls_benefit < b-> ls_benefit;
}

// Return the vertex with best benefit considering the pipeline
bool comp_pipeline(Vertex *v1, Vertex *v2)
{
	return( (v1-> n_unaware_neighs < v2-> n_unaware_neighs) ||
		(v1-> n_unaware_neighs == v2-> n_unaware_neighs && v1-> n_spreader_neighs_gain < v2-> n_spreader_neighs_gain) ||
		(v1-> n_unaware_neighs == v2-> n_unaware_neighs && v1-> n_spreader_neighs_gain == v2-> n_spreader_neighs_gain && v1-> lack_threshold < v2-> lack_threshold) ||
		(v1-> n_unaware_neighs == v2-> n_unaware_neighs && v1-> n_spreader_neighs_gain == v2-> n_spreader_neighs_gain && v1-> lack_threshold == v2-> lack_threshold && v1-> id < v2-> id));
}

bool comp_pipeline_by_ls_sg_criteria(Vertex *v1, Vertex *v2)
{
	return( (v1-> n_unaware_neighs_old < v2-> n_unaware_neighs_old) ||
		(v1-> n_unaware_neighs_old == v2-> n_unaware_neighs_old && v1-> n_spreader_neighs_gain_old < v2-> n_spreader_neighs_gain_old) ||
		(v1-> n_unaware_neighs_old == v2-> n_unaware_neighs_old && v1-> n_spreader_neighs_gain_old == v2-> n_spreader_neighs_gain_old && v1-> lack_threshold_old < v2-> lack_threshold_old) ||
		(v1-> n_unaware_neighs_old == v2-> n_unaware_neighs_old && v1-> n_spreader_neighs_gain_old == v2-> n_spreader_neighs_gain_old && v1-> lack_threshold_old == v2-> lack_threshold_old && v1-> id < v2-> id));
}