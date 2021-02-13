/*
File: random_plus_greedy.h
Description: Header file of random_plus_greedy.cpp
Author: Felipe de C. Pereira

Based on the article "Effective Heuristics for the Perfect Awareness Problem"
Copyright (C) 2020 Felipe de C. Pereira, Pedro J. de Rezebde, Cid C. de Souza
*/

#ifndef RANDOM_PLUS_GREEDY_H
#define RANDOM_PLUS_GREEDY_H

// GRASP construction phase with Random plus Greedy strategie
void random_plus_greedy_construction(int n_vertices, vector<Vertex*> &vertices, vector<Vertex*> &seed_set,
	int &sol_value, int &n_rounds, int &n_aware, int n_seeds_per_insertion, double random_steps_rate, 
	int &draws_altruistic, int &draws_philanthropic, int &steps, int &tot_spreaders, int &tot_edges_explored,
	double &rlc_cl_ratio_mean, mt19937 &rng);

#endif