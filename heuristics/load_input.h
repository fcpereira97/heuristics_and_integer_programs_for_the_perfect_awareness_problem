/*
File: load_input.h
Description: Header file of load_input.cpp
Author: Felipe de C. Pereira

Based on the article "Effective Heuristics for the Perfect Awareness Problem"
Copyright (C) 2020 Felipe de C. Pereira, Pedro J. de Rezebde, Cid C. de Souza
*/

#ifndef LOAD_INPUT_H
#define LOAD_INPUT_H

// Load instance from input file
void load_input(FILE* input_file, int &n_vertices, int &n_edges, vector<Vertex*> &vertices, double threshold_rate);

#endif