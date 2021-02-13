/*
File: preprocessing_by_components.h
Description: Header file of preprocessing_by_components.cpp
Author: Felipe de C. Pereira

Based on the article "Effective Heuristics for the Perfect Awareness Problem"
Copyright (C) 2020 Felipe de C. Pereira, Pedro J. de Rezebde, Cid C. de Souza
*/

#ifndef PREPROCESSING_BY_COMPONENTS_H
#define PREPROCESSING_BY_COMPONENTS_H

using namespace std;

// Identify the connected components
void preprocessing_by_components(int n_vertices, int &n_components, vector<Vertex*> &vertices, vector<vector<Vertex*>> &components);

#endif