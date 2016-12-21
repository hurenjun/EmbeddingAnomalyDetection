/*
 * Graph class firstly load a graph from input file,
 * then give a community structure initialization according to type:
 *		0: random
 *		1: heuristic (graph partition)
 *
 * cmties[i] : the ith community
 *		cmties[i][0]: size of ith community
 *		cmties[i][1]--cmties[i][size]: nodes of ith community
 */

#pragma once
#include "Edge.h"
#include <string>
#include <stddef.h>

using namespace std;

class Graph
{
public:
	Graph(int, string &, int);
	~Graph();
	void printCmties();
	int getNodeCount() const;
	int getEdgeCount() const;
	int getCmtyCount() const;
	Edge* getEdges() const;
	int ** getCmties() const;

private:
	int nodeCount, edgeCount, cmtyCount;
	Edge *edges;
	int **cmties;	

	void loadGraph(string &);
	void randonPartition();
	int findCmty(int, int[]);
	void heuristicPartition();
	void metisPartititon();
	void randomEdgeset();
	int getInterGroupEdgeCount(int *);
};


