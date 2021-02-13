/*
File: preprocessing_by_collapse.h
Description: Header file of preprocessing_by_collapse.cpp
Author: Felipe de C. Pereira

Based on the article "Effective Heuristics for the Perfect Awareness Problem"
Copyright (C) 2020 Felipe de C. Pereira, Pedro J. de Rezebde, Cid C. de Souza
*/

#ifndef PREPROCESSING_BY_COLLAPSE_H
#define PREPROCESSING_BY_COLLAPSE_H

using namespace std;

// We denote by a group, the set of vertices of a maximal connected subgraph in which every vertex has threshold equals to 1
// Preprocessing routine that collapses each group into one new vertex
void preprocessing_by_collapse(int &n_vertices, int &n_edges, vector<Vertex*> &vertices, double threshold_rate, string instance_name,
	string sol_folder_path);
#endif