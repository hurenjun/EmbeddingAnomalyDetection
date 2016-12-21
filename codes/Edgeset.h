/*
 * Edgeset: store a array of either edges or non-edges
 * and implement the stress computing method for corresponding 
 * edgeset or non-edgeset
 */

#pragma once
#include "Edge.h"
#include "Graph.h"
#include "Vector.h"
#include <fstream>
#include <stddef.h>

using namespace std;

class Edgeset{
public:
	Edgeset(const Graph & );
	Edgeset(const Edgeset *);
	Edgeset(const Edgeset &, int **, int, int, int);
	~Edgeset();
	void printSet();
	Edge getEdge(int) const;
	int getSize() const;
	double getBalance() const;
	double getStress(Vector *) const;
	void printStress(Vector *, ofstream &);
private:
	int type = 0;
	bool sampling = false;
	double balance = 0;
	Edge *set = NULL;
	// edgeCount: edge count of the original graph
	// Not Size Of Set
	int size, nodeCount, edgeCount, cmtyCount;

	void randomNonedgeSet(const Edgeset &, int **);
	void getFullNonedgeSet(const Edgeset &);
	bool containEdge(int, int) const;
	int getSamplingSize();
	void copyEdge(Edge *);
};
