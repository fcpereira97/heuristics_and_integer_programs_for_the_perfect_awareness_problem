/*
File: priority_queue.cpp
Description: Methods for fixing a priority queue implemented as heap
Author: Felipe de C. Pereira

Based on the article "Effective Heuristics for the Perfect Awareness Problem"
Copyright (C) 2020 Felipe de C. Pereira, Pedro J. de Rezebde, Cid C. de Souza
*/

#include <bits/stdc++.h>
#include "my_struct.h"

using namespace std;

// Raise an element in heap
int heapify_up(vector<Vertex*> &heap, int node_index, function<bool (Vertex*, Vertex*)> comp_function)
{
	int parent_index = (node_index - 1)/2; // Get parent index

	// While the priority of the node is greater then its parent
	while(node_index > 0 && comp_function(heap[parent_index], heap[node_index]))
	{
		// Swap node and its parent
		swap(heap[node_index]-> cl_pos, heap[parent_index]-> cl_pos);
		swap(heap[node_index], heap[parent_index]);
		node_index = parent_index;
		parent_index = (node_index - 1)/2;
	}

	return node_index; // Return the final index of the element in the array
}

// Heapify the subtree rooted at a given node
// Down an element in heap
void heapify_down(vector<Vertex*> &heap, int root_index, function<bool (Vertex*, Vertex*)> comp_func)
{
	bool stop = false;
	while(!stop)
	{
		int left, right, best;

		// Get children indexes
		left = 2*root_index + 1;
		right = 2*root_index + 2;

		// Get the best child of root
		best = root_index;
		if(left < (int)heap.size() && comp_func(heap[root_index], heap[left]))
			best = left;
		if(right < (int)heap.size() && comp_func(heap[best], heap[right]))
			best = right;

		// If the best child priority is greater than the priority of the root do
		if(best != root_index)
		{
			// Swap best child and root
			swap(heap[root_index]-> cl_pos, heap[best]-> cl_pos);
			swap(heap[root_index], heap[best]); 
			root_index = best;
		}
		else // Else the heap is already correct
		{
			stop = true;
		}
	}
}

// Fix the heap when a node priority has been modified
void heapify(vector<Vertex*> &heap, int node, function<bool (Vertex*, Vertex*)> comp_func)
{
	// 
	node = heapify_up(heap, node, comp_func);
	heapify_down(heap, node, comp_func);
}

// Delete a node from heap
void delete_node_from_heap(vector<Vertex*> &heap, int node, function<bool (Vertex*, Vertex*)> comp_func)
{
	// Swap node to be deleted with last node of heap
	swap(heap[node]-> cl_pos, heap.back()-> cl_pos);
	swap(heap[node], heap.back());

	// Remove last node
	heap.pop_back();

	// Heapify from the original position of the deleted node
	heapify(heap, node, comp_func);
}