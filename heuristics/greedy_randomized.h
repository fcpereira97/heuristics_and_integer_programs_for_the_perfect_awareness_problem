/*
File: greedy_randomized.h
Description: Header file of greedy_randomized.cpp
Author: Felipe de C. Pereira

Based on the article "Effective Heuristics for the Perfect Awareness Problem"
Copyright (C) 2020 Felipe de C. Pereira, Pedro J. de Rezebde, Cid C. de Souza
*/

#ifndef GREEDY_RANDOMIZED_H
#define GREEDY_RANDOMIZED_H

// GRASP construction phase with Greedy Randomized strategie
void greedy_randomized_construction(int n_vertices, vector<Vertex*> &vertices, vector<Vertex*> &seed_set,
	int &sol_value, int &n_rounds, int &n_aware, int n_seeds_per_insertion, double alpha, int &steps,
	int &tot_spreaders, int &tot_edges_explored, double &rlc_cl_ratio_mean, mt19937 &rng);

#endif