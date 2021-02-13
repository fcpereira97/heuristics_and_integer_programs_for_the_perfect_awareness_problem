/*
File: priority_queue.h
Description: Header of priority_queue.cpp
Author: Felipe de C. Pereira

Based on the article "Effective Heuristics for the Perfect Awareness Problem"
Copyright (C) 2020 Felipe de C. Pereira, Pedro J. de Rezebde, Cid C. de Souza
*/

#ifndef PRIORITY_QUEUE_H
#define PRIORITY_QUEUE_H

// Fix the heap when a node priority has been modified
void heapify(vector<Vertex*> &heap, int node_index, function<bool(Vertex*, Vertex*)> comp_func);

// Delete a node from heap
void delete_node_from_heap(vector<Vertex*> &heap, int node_index, function<bool (Vertex*, Vertex*)> comp_func);

#endif