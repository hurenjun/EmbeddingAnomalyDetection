// Class Embedding implement the main procedure of gradient descent 
// Things it do: 
//		1). give each node a vector
//		2). optimize the objective function, whose input is all variables in vectors
// Gradient Descent: line search, backtracking that holdd Armijo rule 

#pragma once
#include "Edgeset.h"
#include "Vector.h"
class Embedding
{
public:
	Embedding(int, double, double, double, string &, int, int);
	Embedding(int, double, double, double, string &, int);
	~Embedding();	
	void gradientDescent();
	int getNodeCount() const;
	int getEdgeCount() const;
	int getCmtyCount() const;
	Edgeset* getEdges() const;
	Vector* getVectors() const;
	void PrintVectors() const;
private:
	int nodeCount, edgeCount, cmntCount;
	int para_k;
	Edgeset *edges = NULL, *nonedges = NULL;
	Edge *doubleEdge, *doubleNonedge;
	Vector *vectors = NULL;
	int **cmnts = NULL;
	int * node_comm = NULL;
	int iteration = 0;
	double *tempDrct = NULL;
	double balance = 0;
	double alpha = 0.01, beta = 0.8, t = 1, eps = 0.001;

	void initialVectors();
	double moveStepSize(int, double);
	void calDescentDirection();
	void iterateVectors(double);
	double objectiveFunction();
	void copyNextToCurrent();
	void objectiveFunctionDetail();
};

