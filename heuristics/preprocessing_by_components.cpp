/*
File: preprocessing_by_components.cpp
Description: Preprocessing routine that identifies connected components
Author: Felipe de C. Pereira

Based on the article "Effective Heuristics for the Perfect Awareness Problem"
Copyright (C) 2020 Felipe de C. Pereira, Pedro J. de Rezebde, Cid C. de Souza
*/

#include <bits/stdc++.h>
#include <sys/resource.h>
#include "my_struct.h"

using namespace std;

// Compare the size of components
bool sort_by_component_size (vector<Vertex*> &v1, vector<Vertex*> &v2)
{
	return ((int)v1.size() > (int)v2.size());
}

// Label each vertex with a null group
void initialize_groups(int n_vertices, vector<Vertex*> &vertices)
{
	for(int i = 0; i < n_vertices; i++)
	{
		vertices[i]-> group_id = -1;
	}
}

// Perform a bfs for finding all vertices of the component
void bfs(Vertex *root, int group_id, vector<Vertex*> &component)
{
	queue<Vertex*> unexplored_vertices; // bfs queue
	unexplored_vertices.push(root); // Add root
	Vertex *vertex;

	root-> group_id = group_id; // Give the new label to the root
	component.push_back(root); // Add root to component

	while(!unexplored_vertices.empty()) // While queue is not empty
	{
		// Get next unexplored vertex v 
		vertex = unexplored_vertices.front();
		unexplored_vertices.pop();
		
		// For all neighbors of v
		for (list<Vertex*>::iterator neighbor = vertex-> neighbors.begin(); neighbor != vertex-> neighbors.end(); ++neighbor)
		{
	    	if((*neighbor)->group_id == -1) // If neighbor has not been labeled yet
	    	{
	    		(*neighbor)->group_id = group_id; // Label neighbor
	    		component.push_back(*neighbor); // Insert neighbor into component
	    		unexplored_vertices.push(*neighbor); // Insert neighbor into queue 
	    	}
		}
	}
}

// Identify the connected components
// Each group corresponds to a connected componet
// For each vertex not labeled, we create a new group and
// perform a bfs to find the entire group
void preprocessing_by_components(int n_vertices, int &n_components, vector<Vertex*> &vertices, vector<vector<Vertex*>> &components)
{
	int n_groups = 0;
	initialize_groups(n_vertices, vertices); //Label all vertex with a null component value

	for (int i = 0; i < n_vertices; ++i) 
	{
		if(vertices[i]-> group_id == -1) // For each vertex not labeled
		{
			vector<Vertex*> component; // Create new component
			n_groups++; // Increase the number of components detected

			bfs(vertices[i], n_groups-1, component); // Find all vertices of the new component by bfs

			components.push_back(component); // Add the new component to the array of components
		}
	}

	n_components = (int)components.size(); // Update number of components
	sort(components.begin(), components.end(), sort_by_component_size); // Sort array of components
}