/*
File: preprocessing_by_collapse.cpp
Description: Preprocessing routine that collapses neighboring vertices that have threshold equal to 1
Author: Felipe de C. Pereira

Based on the article "Effective Heuristics for the Perfect Awareness Problem"
Copyright (C) 2020 Felipe de C. Pereira, Pedro J. de Rezebde, Cid C. de Souza
*/

#include <bits/stdc++.h>
#include "my_struct.h"

using namespace std;

// Label each vertex with a null group id (-1)
void initialize_vertices(int n_vertices, vector<Vertex*> &vertices)
{
	for(int i = 0; i < n_vertices; i++)
		vertices[i]-> group_id = -1;
}

// Perform a Depth-first search in order to label the vertices
void dfs(Vertex* vertex, int group_id)
{
	vertex-> group_id = group_id;

	for (list<Vertex*>::iterator neighbor = vertex-> neighbors.begin(); neighbor != vertex-> neighbors.end(); ++neighbor)
		// If a neighbor with threshold 1 has not yet been labeled then continue the dfs over it
    	if((*neighbor)-> threshold == 1 && (*neighbor)-> group_id == -1) 
    		dfs(*neighbor, group_id);
}

// Identify groups of connected vertices with threshold equals to 1
void identify_groups(int n_vertices, vector<Vertex*> &vertices, int &n_groups)
{
	// For each vertex v
	for(int i = 0; i < n_vertices; i++)
		if(vertices[i]-> threshold == 1) // If v has threshold 1
			if(vertices[i]-> group_id == -1) // If v has not been labeled yet
			{
				n_groups++; // Count new group
				dfs(vertices[i], n_groups -1); // Run dfs to find other vertices in this new group
			}
}


// For each group G, include in G any vertex that has neighbors only in G
void expand_groups(int n_vertices, vector<Vertex*> &vertices)
{
	bool expansion_occurred; // Flag that indicates if an expansion has occured

	do
	{
		expansion_occurred = false; 

		// For each vertex v
		for(int i = 0; i < n_vertices; i++) 
		{
			// If v has not been labeled yet and v has at least one neighbor
			if(vertices[i]-> group_id == -1 && !vertices[i]->neighbors.empty())
			{
				// Get the group if of the first neighbor of v
				int group_id = vertices[i]->neighbors.front()->group_id; 

				// If the group if is not null
				if(group_id != -1)
				{
					// Verify if all neighbors of v belong to the same group
					bool expanded_this_group = true;
					for (list<Vertex*>::iterator neighbor = vertices[i]-> neighbors.begin(); neighbor != vertices[i]-> neighbors.end(); ++neighbor)
		   			{
				    	if((*neighbor)-> group_id != group_id)
				    	{
				    		expanded_this_group = false;
				    		break;
				    	}
		   			}

		   			// If the group can be expanded
		   			if(expanded_this_group)
		   			{
		   				vertices[i]-> group_id = group_id; // Insert v to the group
		   				expansion_occurred = true; // Flag indicates that at least one expansion has ocurred
		   			}
				}
			}
		}
	} while(expansion_occurred);// Repeat until no expansion has ocurred, that is, no vertice was inserted into a group
}

// Fix adjacency lists of the vertices
// Insert into the adjacency lists the new vertices that represents the groups
// Delete from the adjacency lists the neighbors that were collapsed into a group
void fix_neighbors(int n_vertices, vector<Vertex*> &vertices, vector<Vertex*> &colapsed_vertices)
{
	for(int i = 0; i < n_vertices; i++) // For each vertex v
	{
		if(vertices[i]-> group_id == -1) // If v does not belong to a group
		{
			list<Vertex*> neighbors_aux; // Auxiliar adjacency list of v

			// For each neighbor u of v
			for (list<Vertex*>::iterator neighbor = vertices[i]-> neighbors.begin(); neighbor != vertices[i]-> neighbors.end(); ++neighbor)
			{
				int neighbor_group_id = (*neighbor)-> group_id; // Get the id of group that the u belongs to

				if(neighbor_group_id != -1) // If u belong to a group
				{
					// Insert the vertex that represents the group into the auxiliar adjacency list of v
					neighbors_aux.push_back(colapsed_vertices[neighbor_group_id]);
					// Insert v into the adjacency list of the vertice that represents the group
					colapsed_vertices[neighbor_group_id]-> neighbors.push_back(vertices[i]);
				}
				else // If u does not belong o a group
				{
					neighbors_aux.push_back(*neighbor); // Insert u into the auxiliar adjacency list of v 
				}
			}
			// Rewrite the adjacency list of v
			vertices[i]->neighbors.assign(neighbors_aux.begin(), neighbors_aux.end());
		}
	}
}

// Fix array of vertices
// Add new vertices that represent ew vertices that represents the groups
// Delete the vertices that were collapsed into a group
// For each vertex produced by collapse, build a list of original ids of the collapsed vertices
void fix_vertices(int &n_vertices, vector<Vertex*> &vertices, int n_groups, vector<Vertex*> &colapsed_vertices)
{
	list<Vertex*> vertices_aux; // Auxiliar list of vertices

	// For each original vertex v
	for(int i = 0; i < n_vertices; i++)
	{
		int v_group_id = vertices[i]-> group_id; // Get the id of the group that v belongs to
		if(v_group_id == -1) // If v does not belong to a group
		{
			// Insert v's original id into the list of v's ids before the preprocessing
			vertices[i]-> ids_before_preprocessing.push_back(vertices[i]-> id); 
			vertices_aux.push_back(vertices[i]); // Add v in the auxiliar array
		}
		else
		{
			// Insert v's original id into the list of the u's ids before the preprocessing
			// Where u is the vertex that represents a group
			colapsed_vertices[v_group_id]-> ids_before_preprocessing.push_back(vertices[i]-> id);
		}
	}

	// For each group, add the vertex that represents it into the auxiliar array
	for(int i = 0; i < n_groups; i++)
	{
		colapsed_vertices[i]-> degree = colapsed_vertices[i]-> neighbors.size();
		vertices_aux.push_back(colapsed_vertices[i]);
	}

	// Rewrite array of vertices
	vertices.assign(vertices_aux.begin(), vertices_aux.end());
	n_vertices = vertices.size(); // Fix the number of vertices

	for(int i = 0; i < n_vertices; i++)
		vertices[i]-> id = i; // Set new ids for the vertices
}

// Write a file that maps the old ids (before preprocessing) into new ids (after preprocessinf)
// A vertex that was produced by a collapse has its id mapped into a list of old ids
void print_mapping_file(int n_vertices, vector<Vertex*> &vertices, string instance_name, string sol_folder_path)
{
	// Create map file
	FILE *mapping_file;
	string mapping_file_path = sol_folder_path + instance_name + ".map";
	mapping_file = fopen(mapping_file_path.c_str(), "w");

	// For each vertex v
	for(int i = 0; i < n_vertices; i++)
	{
		fprintf(mapping_file, "%d: ", vertices[i]-> id); // Print the new id of v
		
		// For each old id before the preprocessing
		for (list<int>::iterator old_id = vertices[i]-> ids_before_preprocessing.begin(); old_id != vertices[i]-> ids_before_preprocessing.end(); ++old_id)
			fprintf(mapping_file, "%d ", *old_id); // Print the old id

		fprintf(mapping_file, "\n");
	}
	fclose(mapping_file);
}

// We denote by a group, the set of vertices of a maximal connected subgraph in which every vertex has threshold equals to 1
// Preprocessing routine that collapses each group into one new vertex
void preprocessing_by_collapse(int &n_vertices, int &n_edges, vector<Vertex*> &vertices, double threshold_rate, string instance_name,
	string sol_folder_path)
{
	int n_groups = 0; // Number of groups
	initialize_vertices(n_vertices, vertices); // Label each vertex with a null group identifier (-1) 
	identify_groups(n_vertices, vertices, n_groups); // Identify groups
	expand_groups(n_vertices, vertices); // For each group G, include in G any vertex that has neighbors only in G
	
	// Create a array of vertices representing each identified group
	// Each group will be collapsed into a single vertex
	vector<Vertex*> colapsed_vertices(n_groups); 
	for (int i = 0; i < n_groups; ++i)
	{
		colapsed_vertices[i] = new Vertex;
		colapsed_vertices[i]-> threshold = 1;
	}

	// Fix graph due to the collapse
	fix_neighbors(n_vertices, vertices, colapsed_vertices);
	fix_vertices(n_vertices, vertices, n_groups, colapsed_vertices);

	// Fix the number of edges
	n_edges = 0;
	for(int i = 0; i < n_vertices; i++)
		n_edges += vertices[i]-> degree;
	n_edges = n_edges/2;

	// Print map file of ids
	print_mapping_file(n_vertices, vertices, instance_name, sol_folder_path);
}