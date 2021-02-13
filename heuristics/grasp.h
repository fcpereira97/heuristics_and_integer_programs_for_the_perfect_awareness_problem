#ifndef GRASP_H
#define GRASP_H

// Main method of grasp algorithm
void grasp(int original_n_vertices, int original_n_edges, int n_components,
	vector<vector<Vertex*>> *components, string construction_phase_flag,
	double construction_phase_parameter, int n_seeds_per_insertion, string local_search_phase_flag,
	double local_search_rate, double time_begin, int min_time, int time_limit, int iteration_limit,
	int target,	string instance_name, FILE* output_file, FILE* sol_file, mt19937 &rng,
	vector<vector<Vertex*>> *best_seed_set_comp, int heap_of_seed_sets_size);

#endif