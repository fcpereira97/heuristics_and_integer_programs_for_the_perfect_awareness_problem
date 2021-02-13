/*
File: load_input.cpp
Description: Load a given instance of the PAP from an input file
Author: Felipe de C. Pereira

Based on the article "Effective Heuristics for the Perfect Awareness Problem"
Copyright (C) 2020 Felipe de C. Pereira, Pedro J. de Rezebde, Cid C. de Souza
*/

#include <bits/stdc++.h>
#include "my_struct.h"

using namespace std;

// Compare two edges by the lowest ids of its vertices
bool compare_two_edges_by_endpoints_ids(pair<int, int> edge_x, pair<int, int> edge_y)
{
	if(edge_x.first == edge_y.first)
		return edge_x.second < edge_y.second;
	else
		return edge_x.first < edge_y.first;
}

// Load the number of vertices and edges of the graph
void load_graph_size(FILE* input_file, int &n_vertices, int &n_edges)
{	
	char trash[50];
	if(!fscanf(input_file, "%s %s %d %d", trash, trash, &n_vertices, &n_edges))
		printf("Error reading the input file.\n");
}

// Load edges
void load_edges(FILE* input_file, int &n_vertices, int &n_edges, vector<Vertex*> &vertices)
{
	int next_index = 0; // Identifier to be assigned to the next discovered vertex
	int n_non_duplicated_edges = 0; // Number of non-duplicated edges
	int v1_index, v2_index, v1_id, v2_id; // Auxiliar variables
	vector<pair<int, int>> edges; // Array of edges

	// Auxiliar hash map where the key is the original id from input file
	//and the value is the index of vertice in array of vertices
	unordered_map<int, int> id_index_map;

	// Load each edge from input. 
	for(int i = 0; i < n_edges; i++) 
	{
		if(!fscanf(input_file, "%d %d", &v1_id, &v2_id))
			printf("Error reading the input file.\n");

		// Discard edge if its a loop
		// Notice that if a vertex is the only one vertes in its connected component, then this vertex is discarded
		if(v1_id != v2_id)
		{
			// The first endpoint (source) is the one with less original id
			edges.push_back(make_pair(min(v1_id, v2_id), max(v1_id, v2_id)));
		}
	}

	// Sort	edges by its endpoints ids
	sort(edges.begin(), edges.end(), compare_two_edges_by_endpoints_ids);

	// For each edge loaded
	for(int i = 0; i < (int)edges.size(); i++)
	{
		// If it is the first edge or it is not a non-duplicated edge
		// Notice that it discards the edge if it is duplicated
		if(i == 0 || (edges[i].first != edges[i-1].first) || (edges[i].second != edges[i-1].second))
		{
			// Get the endpoints' ids
			v1_id = edges[i].first;
			v2_id = edges[i].second;

			// If a new vertex was discovered, assign the next_index to it
			if(id_index_map.find(v1_id) == id_index_map.end())
			{
				id_index_map[v1_id] = next_index;
				next_index++;
				vertices[id_index_map[v1_id]]-> id = v1_id;
			}
			if(id_index_map.find(v2_id) == id_index_map.end())
			{
				id_index_map[v2_id] = next_index;
				next_index++;
				vertices[id_index_map[v2_id]]-> id = v2_id;
			}

			// Auxiliar variables receive the ids assigned to the vertices
			v1_index = id_index_map[v1_id];
			v2_index = id_index_map[v2_id];

			// Update the list of neighbors of both of endpoints of the edge
			if(v1_index != v2_index)
			{
				vertices[v1_index]-> neighbors.push_back(vertices[v2_index]);
				vertices[v2_index]-> neighbors.push_back(vertices[v1_index]);
			}

			// Increase the number of non-duplicated edges
			n_non_duplicated_edges++;
		}
	}
	n_edges = n_non_duplicated_edges; // Fix the real number of edges
	n_vertices = next_index; // Fix the real number of vertices
	vertices.resize(n_vertices); // Fix the vertices array due to the real number of vertices
}

// Load instance from input file
void load_input(FILE* input_file, int &n_vertices, int &n_edges, vector<Vertex*> &vertices, double threshold_rate)
{
	// Load graph size
	load_graph_size(input_file, n_vertices, n_edges);

	// Create vertices
	vertices.resize(n_vertices);
	for(int i = 0; i < n_vertices; i++)
	{
		vertices[i] = new Vertex;
		vertices[i]-> degree = 0;
		vertices[i]-> n_spreader_neighs = 0;
		vertices[i]-> state = 0;
	}

	// Load edges list
	load_edges(input_file, n_vertices, n_edges, vertices);

	// Set degrees and thresholds
	for(int i = 0; i < n_vertices; i++)
	{
		vertices[i]-> degree = vertices[i]-> neighbors.size();
		vertices[i]-> threshold = ceil((double)vertices[i]-> degree * threshold_rate);
	}
}