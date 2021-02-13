/*
File: local_search.h
Description: Header file of local_search.cpp
Author: Felipe de C. Pereira

Based on the article "Effective Heuristics for the Perfect Awareness Problem"
Copyright (C) 2020 Felipe de C. Pereira, Pedro J. de Rezebde, Cid C. de Souza
*/

#ifndef LOCAL_SEARCH_H
#define LOCAL_SEARCH_H

// Delete from the solution every seed that has the number of neighboring seeds
// more or equal than its threshold
void clean_solution(int n_vertices, vector<Vertex*> &vertices,
	vector<Vertex*> &seed_set,	int &sol_value,  int &n_rounds, int &n_aware);

// Identify and delete exceeding seeds by simulating the propagation
// with a subset of the seed set that grows progressively
void local_search_progressive(int n_vertices, vector<Vertex*> &vertices,
	vector<Vertex*> &seed_set, int &sol_value, int &n_rounds, int &n_aware,
	string construction_phase_flag, double construction_phase_parameter);

// Evaluates, seed by seed, whether it is possible to remove the seed v
// If so, v is deleted from the solution
void local_search_one_by_one(int n_vertices, vector<Vertex*> &vertices,
	vector<Vertex*> &seed_set, int &sol_value, int &n_rounds, int &n_aware);

// Evaluates, whether it is possible to remove an element from a subset S' of the seed seed S
// If so, these elements are deleted from S
void local_search_by_sets(int n_vertices, vector<Vertex*> &vertices,
	vector<Vertex*> &seed_set, int &sol_value, int &n_rounds, int &n_aware,
	double division_parameter);
#endif