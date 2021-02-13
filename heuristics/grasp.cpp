#include <bits/stdc++.h>
#include <sys/resource.h>
#include "my_struct.h"
#include "propagation.h"
#include "greedy_randomized.h"
#include "weighted_greedy_randomized.h"
#include "random_plus_greedy.h"
#include "sampled_greedy.h"
#include "local_search.h"

using namespace std;

// Macro for getting runtime execution
// this macro produces a time equal to the one produced by clock(),
// but does not suffer from the wraparound problem of clock()
extern int getrusage();
#define CPUTIME(ruse) (getrusage(RUSAGE_SELF,&ruse),ruse.ru_utime.tv_sec + ruse.ru_stime.tv_sec + 1e-6 * (ruse.ru_utime.tv_usec + ruse.ru_stime.tv_usec))
struct rusage grb_ruse_grasp; // Global variable used on time counter

// Print each vertex from the graph by listing its id, degree and ides of its neighs
void print_graph(int n_vertices, vector<Vertex*> *vertices){
	for(int i = 0; i < n_vertices; i++)
	{
		cout << "Vertex " << (*vertices)[i]-> id << endl;
		cout << "Degree: " << (*vertices)[i]-> degree << endl;
		cout << "Threshold: " << (*vertices)[i]-> threshold << endl;
		cout << "Neighs: ";
		for (list<Vertex*>::iterator it = (*vertices)[i]-> neighbors.begin(); it != (*vertices)[i]-> neighbors.end(); ++it)
    		cout << (*it)-> id << " ";
    	cout << endl << endl;
	}
}

void sort_seeds_gr_criteria(vector<Vertex*> * seed_set)
{
	Vertex *v;

	for(int i = 0; i < (int)(*seed_set).size(); i++)
	{
		v = (*seed_set)[i];
		v-> ls_benefit = 0;

		for (list<Vertex*>::iterator neigh = v-> neighbors.begin(); neigh != v-> neighbors.end(); ++neigh)
			if((*neigh)-> state == 1 && (*neigh)-> n_spreader_neighs < 2)
				v-> ls_benefit++;
	}
	sort((*seed_set).begin(), (*seed_set).end(), compare_two_vertices_by_ls_gr_criteria);
}

bool heap_compare(vector<Vertex*> &seed_set_1, vector<Vertex*> &seed_set_2)
{
	return (int)(seed_set_2).size() < (int)(seed_set_1).size(); 
}

// Heapify the subtree rooted at a given node
void heapify(vector<vector<Vertex*>> *heap, function<bool (vector<Vertex*>&, vector<Vertex*>&)> comp_func)
{
	int root = 0;
	bool stop = false;
	while(!stop)
	{
		int left, right, best;

		// Get children indices
		left = 2*root + 1;
		right = 2*root + 2;

		best = root;

		// Get the best child of root
		if(left < (int)(*heap).size() && comp_func((*heap)[root], (*heap)[left]))
			best = left;
		if(right < (int)(*heap).size() && comp_func((*heap)[best], (*heap)[right]))
			best = right;

		if(best != root)
		{
			swap((*heap)[root], (*heap)[best]); // Swap best child and root
			root = best;
		}
		else
		{
			stop = true;
		}
	}
}

// Main method of grasp algorithm
void grasp(int original_n_vertices, int original_n_edges, int n_components,
	vector<vector<Vertex*>> *components, string construction_phase_flag,
	double construction_phase_parameter, int n_seeds_per_insertion, string local_search_phase_flag,
	double local_search_rate, double time_begin, int min_time, int time_limit, int iteration_limit,
	int target,	string instance_name, FILE* output_file, FILE* sol_file, mt19937 &rng,
	vector<vector<Vertex*>> *best_seed_set_comp, int heap_of_seed_sets_size)
{
	//double mean_it_best_solution = 0;
	//int it_counter_at_1000 = 0;

	// Solution variables
	int best_sol_value;
	best_sol_value = INT_MAX;
	vector<vector<Vertex*>> seed_set_comp(n_components);
	vector<vector<vector<Vertex*>>> heap_of_seed_sets_comp(n_components);
	vector<int> best_sol_value_comp(n_components, INT_MAX);

	int sol_const, sol_ls, sol_value, comp_sol_value, n_rounds, comp_n_rounds, n_aware, comp_n_aware, tot_spreaders, tot_edges_explored, steps, it_best_sol, sol_const_heap;
	sol_value = comp_sol_value = n_rounds = comp_n_rounds = n_aware = comp_n_aware = tot_spreaders = tot_edges_explored = steps = it_best_sol = sol_const_heap = 0;

	double total_time, mean_time_const, mean_time_ls, time_const, time_ls, time_const_start, time_const_end, time_ls_start, time_ls_end;
	total_time = mean_time_const = mean_time_ls = 0;

	double mean_sol, mean_sol_const, rlc_cl_ratio_mean, mean_rlc_cl_ratio_mean, mean_tot_spreaders, mean_steps, mean_tot_edges_explored;
	mean_sol = mean_sol_const = rlc_cl_ratio_mean = mean_rlc_cl_ratio_mean = mean_tot_spreaders = mean_steps = mean_tot_edges_explored = 0;

	double best_gain_ls, mean_gain_ls, mean_draws_altruistic, mean_draws_philanthropic;
	best_gain_ls = mean_gain_ls = mean_draws_altruistic = mean_draws_philanthropic = 0;

	int draws_altruistic, draws_philanthropic;
	draws_altruistic = draws_philanthropic = 0;

	int n_iteration = 1;

	// Loop of GRASP iterations
	while(best_sol_value > target && ((CPUTIME(grb_ruse_grasp) < min_time) || (CPUTIME(grb_ruse_grasp) < time_limit && n_iteration <= iteration_limit)))
	{
		sol_value = 0;
		n_rounds = 0;
		n_aware = 0;
		sol_const = 0;
		sol_const_heap = 0;
		sol_ls = 0;
		time_const = 0;
		time_ls = 0;
		draws_altruistic = 0;
		draws_philanthropic = 0;
		tot_spreaders = 0;
		tot_edges_explored = 0;
		rlc_cl_ratio_mean = 0;
		steps = 0;

		int i = 0;
		while(i < n_components)
		{
			//cout << "component " << i + 1 << " of " << n_components << endl;
			seed_set_comp[i].clear();
			time_const_start = CPUTIME(grb_ruse_grasp);

			// Construction phase
			if(construction_phase_flag == "GR")
			{
				greedy_randomized_construction((int)(*components)[i].size(), (*components)[i],
					seed_set_comp[i], comp_sol_value, comp_n_rounds, comp_n_aware,
					n_seeds_per_insertion, construction_phase_parameter, steps, tot_spreaders,
					tot_edges_explored, rlc_cl_ratio_mean, rng);
			}
			else if (construction_phase_flag == "WGR")
			{
				weighted_greedy_randomized_construction((int)(*components)[i].size(), (*components)[i],
					seed_set_comp[i], comp_sol_value, comp_n_rounds, comp_n_aware,
					n_seeds_per_insertion, construction_phase_parameter, steps, tot_spreaders,
					tot_edges_explored, rlc_cl_ratio_mean, rng);
			}
			else if (construction_phase_flag == "RG")
			{
				random_plus_greedy_construction((int)(*components)[i].size(), (*components)[i],
					seed_set_comp[i], comp_sol_value, comp_n_rounds, comp_n_aware, n_seeds_per_insertion,
					construction_phase_parameter, draws_altruistic, draws_philanthropic, steps, tot_spreaders,
					tot_edges_explored, rlc_cl_ratio_mean, rng);
			}
			else if(construction_phase_flag == "SG")
			{
				sampled_greedy_construction((int)(*components)[i].size(), (*components)[i],
					seed_set_comp[i], comp_sol_value, comp_n_rounds, comp_n_aware, n_seeds_per_insertion,
					construction_phase_parameter, draws_altruistic, draws_philanthropic, steps, tot_spreaders,
					tot_edges_explored, rlc_cl_ratio_mean, rng);
			}
			
			if(CPUTIME(grb_ruse_grasp) >= time_limit)
				break;

			clean_solution((int)(*components)[i].size(), (*components)[i], seed_set_comp[i],
				comp_sol_value, comp_n_rounds, comp_n_aware);

			sol_const += comp_sol_value;
			time_const_end = CPUTIME(grb_ruse_grasp);
			time_const += time_const_end - time_const_start;

			time_ls_start = CPUTIME(grb_ruse_grasp);

			if(construction_phase_flag == "GR" || construction_phase_flag == "WGR")
				sort_seeds_gr_criteria(&seed_set_comp[i]);
			else if(construction_phase_flag == "SG")
				sort((seed_set_comp[i]).begin(), (seed_set_comp[i]).end(), comp_pipeline_by_ls_sg_criteria);

			if(heap_of_seed_sets_size == 0)
			{
				local_search_progressive((int)(*components)[i].size(), (*components)[i],
					seed_set_comp[i], comp_sol_value, comp_n_rounds, comp_n_aware,
					construction_phase_flag, construction_phase_parameter);

				if(local_search_phase_flag == "OBO")
					local_search_one_by_one((int)(*components)[i].size(), (*components)[i],
						seed_set_comp[i], comp_sol_value, comp_n_rounds, comp_n_aware);
				else if(local_search_phase_flag == "SBS")
					local_search_by_sets((int)(*components)[i].size(), (*components)[i],
						seed_set_comp[i], comp_sol_value, comp_n_rounds, comp_n_aware,
						local_search_rate);
			}
			else if(n_iteration <= heap_of_seed_sets_size)
			{
				heap_of_seed_sets_comp[i].push_back(seed_set_comp[i]);

				if(n_iteration == heap_of_seed_sets_size)
					make_heap(heap_of_seed_sets_comp[i].begin(), heap_of_seed_sets_comp[i].end(), heap_compare);
			}
			else
			{
				sol_const_heap += (int)heap_of_seed_sets_comp[i].front().size();
				swap(seed_set_comp[i], heap_of_seed_sets_comp[i].front());
				heapify(&heap_of_seed_sets_comp[i], heap_compare);
				local_search_progressive((int)(*components)[i].size(), (*components)[i],
					seed_set_comp[i], comp_sol_value, comp_n_rounds, comp_n_aware,
					construction_phase_flag, construction_phase_parameter);

				if(local_search_phase_flag == "OBO")
					local_search_one_by_one((int)(*components)[i].size(), (*components)[i],
						seed_set_comp[i], comp_sol_value, comp_n_rounds, comp_n_aware);
				else if(local_search_phase_flag == "SBS")
					local_search_by_sets((int)(*components)[i].size(), (*components)[i],
						seed_set_comp[i], comp_sol_value, comp_n_rounds, comp_n_aware,
						local_search_rate);
			}
			

			sol_ls += comp_sol_value;
			time_ls_end = CPUTIME(grb_ruse_grasp);
			time_ls += time_ls_end - time_ls_start;

			if(CPUTIME(grb_ruse_grasp) >= time_limit)
				break;

			sol_value += comp_sol_value;
			n_aware += comp_n_aware;

			i++;
		}

		rlc_cl_ratio_mean /= steps;

		if(heap_of_seed_sets_size == 0)
		{
			mean_sol += sol_value;
			mean_sol_const += sol_const;
			mean_gain_ls += sol_const - sol_ls;
		}
		else if(n_iteration > heap_of_seed_sets_size)
		{
			mean_sol += sol_value;
			mean_sol_const += sol_const_heap;
			mean_gain_ls += sol_const_heap - sol_ls;
		}
		mean_time_const += time_const;
		mean_time_ls += time_ls;
		mean_draws_altruistic += draws_altruistic;
		mean_draws_philanthropic += draws_philanthropic;
		mean_tot_spreaders += tot_spreaders;
		mean_steps += steps;
		mean_tot_edges_explored += tot_edges_explored;
		mean_rlc_cl_ratio_mean += rlc_cl_ratio_mean;

		if(sol_ls - sol_const_heap > best_gain_ls)
			best_gain_ls = sol_ls - sol_const_heap;

		//Update best solutions per component
		bool new_incumbent = false;
		for(int i = 0; i < n_components; i++)
		{
			if((int)seed_set_comp[i].size() < best_sol_value_comp[i])
			{
				new_incumbent = true;
				best_sol_value_comp[i] = seed_set_comp[i].size();
				(*best_seed_set_comp)[i].clear();
				(*best_seed_set_comp)[i].assign(seed_set_comp[i].begin(), seed_set_comp[i].end());
			}
		}

		if(new_incumbent)
		{
			it_best_sol = n_iteration;
			best_sol_value = 0;
			for(int i = 0; i < n_components; i++)
				best_sol_value += best_sol_value_comp[i];

			printf("[Iteration %d] Solution value = '%d' (%d) [%d] {%d} -- New incumbent solution found!\n", n_iteration, sol_const, sol_const_heap, sol_value, best_sol_value);
		}
		else
		{
			printf("[Iteration %d] Solution value = '%d' (%d) [%d] {%d}\n", n_iteration, sol_const, sol_const_heap, sol_value, best_sol_value);
		}

		/*
		if(n_iteration % 1000 == 0)
		{
			cout << "It best = " << it_best_sol << endl;
			mean_it_best_solution += (it_best_sol- it_counter_at_1000*1000);
			best_sol_value = INT_MAX;
			it_counter_at_1000++;
		}
		*/
		n_iteration++;
	}

	if(n_iteration > iteration_limit)
		n_iteration = iteration_limit;
	
	total_time = CPUTIME(grb_ruse_grasp) - time_begin;
	if(total_time > time_limit)
		total_time = time_limit;

	mean_sol /= (n_iteration - heap_of_seed_sets_size);
	mean_sol_const /= (n_iteration - heap_of_seed_sets_size);
	mean_time_const /= n_iteration;
	mean_time_ls /= (n_iteration - heap_of_seed_sets_size);
	mean_gain_ls /= (n_iteration - heap_of_seed_sets_size);
	mean_draws_philanthropic /= n_iteration;
	mean_draws_altruistic /= n_iteration;
	mean_tot_spreaders /= n_iteration;
	mean_steps /= n_iteration;
	mean_rlc_cl_ratio_mean /= n_iteration;
	mean_tot_edges_explored /= n_iteration;
	//mean_it_best_solution /= 30;

	// Print info on console
	//printf("\nNumber of iterations = %d\n", n_iteration);
	//printf("Mean of solution values = %lf\n\n", mean_sol);
	printf("\nBest solution value: %d\n", best_sol_value);
	printf("Best solution found at iteration: %d\n", it_best_sol);
	printf("Total elapsed time: %.2lf\n\n", total_time);

	fprintf(output_file, "%s, %s, %d, %d, %.10lf, %.10lf, %.10lf, %d, %d, %.10lf, %.10lf, %.10lf, %d, %.10lf, %.10lf, %.10lf, %.10lf, %.10lf, %.10lf\n",
		instance_name.c_str(), construction_phase_flag.c_str(), original_n_vertices, original_n_edges, mean_sol_const,
		mean_sol, mean_gain_ls, best_sol_value, it_best_sol, mean_time_const, mean_time_ls, total_time, n_iteration,
		mean_steps, mean_rlc_cl_ratio_mean, mean_tot_spreaders, mean_tot_edges_explored, mean_draws_altruistic,
		mean_draws_philanthropic);

	// Print on solution file
	fprintf(sol_file, "%d\n", best_sol_value);
	for(int i = 0; i < n_components; i++)
	{
		for(int j = 0; j < (int)(*best_seed_set_comp)[i].size(); j++)
		{
			fprintf(sol_file, "%d\n", (*best_seed_set_comp)[i][j]-> id);
		}
	}

	//Print on ttt-plot file
	if(target > 0)
	{
		FILE *ttt_plot_file;
		transform (construction_phase_flag.begin(), construction_phase_flag.end(), construction_phase_flag.begin(), ::tolower);
		string ttt_plot_path = "ttt-plot_" + construction_phase_flag + "/" + instance_name + ".ttt";
		ttt_plot_file = fopen(ttt_plot_path.c_str(), "a");
		fprintf(ttt_plot_file, "%.8lf\n", total_time);
		fclose(ttt_plot_file);
	}
}