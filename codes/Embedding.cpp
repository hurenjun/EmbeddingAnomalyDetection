#include "Embedding.h"
#include "AnomalyDetection.h"
#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include <time.h>
#include <algorithm>

using namespace std;

Embedding::Embedding(int sign, double alpha, double beta, double eps, string &path, int dimension, int i_k) {
	this->alpha = alpha;
	this->beta = beta;
	this->eps = eps;
	this->para_k = i_k;
	/*
	* parameter is the method of graph partition
	* 	0: random method
	*  1: method recommended in the article
	*/
	Graph g(sign % 3, path, dimension);

	this->nodeCount = g.getNodeCount();
	this->edgeCount = g.getEdgeCount();
	this->cmntCount = g.getCmtyCount();
	cmnts = g.getCmties();
	//int ** comm_temp = g.getCmties();
	//cmnts = new int*[this->cmntCount];
	//for (int i = 0; i < this->cmntCount; i++) {
	//	cmnts[i] = new int[comm_temp[i][0] + 1];
	//	cmnts[i][0] = comm_temp[i][0];
	//	for (int j = 1; j <= comm_temp[i][0]; j++) {
	//		cmnts[i][j] = comm_temp[i][j];
	//	}
	//}
	this->initialVectors();	

	this->edges = new Edgeset(g);
	this->nonedges = new Edgeset(*this->edges, this->cmnts, this->nodeCount, this->edgeCount, this->cmntCount);
	this->doubleEdge = new Edge[2 * edgeCount];
	this->doubleNonedge = new Edge[2 * edgeCount];
	for (int i = 0; i < edgeCount; i++) {
		doubleEdge[2 * i].setPoint(edges->getEdge(i).getStart(), edges->getEdge(i).getEnd());
		doubleEdge[2 * i + 1].setPoint(edges->getEdge(i).getEnd(), edges->getEdge(i).getStart());
		doubleNonedge[2 * i].setPoint(nonedges->getEdge(i).getStart(), nonedges->getEdge(i).getEnd());
		doubleNonedge[2 * i + 1].setPoint(nonedges->getEdge(i).getEnd(), nonedges->getEdge(i).getStart());
	}
	sort(doubleEdge, doubleEdge + 2 * edgeCount);
	sort(doubleNonedge, doubleNonedge + 2 * edgeCount);
	this->balance = this->nonedges->getBalance();	

	tempDrct = new double[cmntCount];
}

Embedding::Embedding(int sign, double alpha, double beta, double eps, string &path, int dimension) {
	this->alpha = alpha;
	this->beta = beta;
	this->eps = eps;
	/*
	* parameter is the method of graph partition
	* 	0: random method
	*  1: method recommended in the article
	*/
	Graph g(sign % 3, path, dimension);

	this->nodeCount = g.getNodeCount();
	this->edgeCount = g.getEdgeCount();
	this->cmntCount = g.getCmtyCount();
	cmnts = g.getCmties();
	this->initialVectors();

	this->edges = new Edgeset(g);
	this->nonedges = new Edgeset(*this->edges, this->cmnts, this->nodeCount, this->edgeCount, this->cmntCount);
	this->doubleEdge = new Edge[2 * edgeCount];
	this->doubleNonedge = new Edge[2 * edgeCount];
	for (int i = 0; i < edgeCount; i++) {
		doubleEdge[2 * i].setPoint(edges->getEdge(i).getStart(), edges->getEdge(i).getEnd());
		doubleEdge[2 * i + 1].setPoint(edges->getEdge(i).getEnd(), edges->getEdge(i).getStart());
		doubleNonedge[2 * i].setPoint(nonedges->getEdge(i).getStart(), nonedges->getEdge(i).getEnd());
		doubleNonedge[2 * i + 1].setPoint(nonedges->getEdge(i).getEnd(), nonedges->getEdge(i).getStart());
	}
	sort(doubleEdge, doubleEdge + 2 * edgeCount);
	sort(doubleNonedge, doubleNonedge + 2 * edgeCount);
	this->balance = this->nonedges->getBalance();

	tempDrct = new double[cmntCount];
}

Embedding::~Embedding() {
	delete edges;			edges = NULL;
	delete nonedges;		nonedges = NULL;
	delete[] vectors;		vectors = NULL;
	delete[] tempDrct;		tempDrct = NULL;
	delete[] doubleEdge;	doubleEdge = NULL;
	delete[] doubleNonedge;	doubleNonedge = NULL;
	if (node_comm != NULL) { delete[] node_comm; node_comm = NULL; }
	//if (cmnts != NULL) {
	//	for (int i = 0; i < cmntCount; i++) { delete[] cmnts[i]; cmnts[i] = NULL; }
	//	delete[] cmnts; cmnts = NULL;
	//}
}

void Embedding::initialVectors() {
	this->vectors = new Vector[this->nodeCount];
	node_comm = new int[this->nodeCount];

	for (int i = 0; i < this->cmntCount; i++) {
		for (int j = 1; j <= this->cmnts[i][0]; j++) {
			node_comm[this->cmnts[i][j]] = i;
		}
	}

	//int p = min((int)ceil(log2(nodeCount * 1.0)), cmntCount);
	int p = min((int)ceil(edgeCount * 2.0 / nodeCount), cmntCount);		// avg degree
	//int p = para_k;
	int q = min(p / 4, cmntCount - p);
	cout << p << " " << q << endl;
	for (int i = 0; i < this->nodeCount; i++) {
		this->vectors[i].initialVector(node_comm[i], p, q);
		//this->vectors[i].initialVector(nodeCmtyID[i], cmntCount / 20, cmntCount / 20 / 3);
	}
}

void Embedding::gradientDescent() {
	double currentOF = 0, lastOF = 0;
	iteration = 0;

	lastOF = this->objectiveFunction();

	while (true) {
		iteration++;
		if (nonedges != NULL) {
			delete nonedges;		nonedges = NULL;
		}
		this->nonedges = new Edgeset(*this->edges, this->cmnts, this->nodeCount, this->edgeCount, this->cmntCount);
		for (int i = 0; i < edgeCount; i++) {
			doubleNonedge[2 * i].setPoint(nonedges->getEdge(i).getStart(), nonedges->getEdge(i).getEnd());
			doubleNonedge[2 * i + 1].setPoint(nonedges->getEdge(i).getEnd(), nonedges->getEdge(i).getStart());
		}
		sort(doubleNonedge, doubleNonedge + 2 * edgeCount);
		this->balance = this->nonedges->getBalance();

		this->calDescentDirection();
		currentOF = this->moveStepSize(1, lastOF);

		this->copyNextToCurrent();
		//cout << "Iteration " << iteration << ": " << currentOF << endl;
		if (fabs(lastOF - currentOF) / currentOF < eps || iteration == 50) break;
		else lastOF = currentOF;
	}

	cout << "Iteration Round Quantity: " << iteration << endl;
	//objectiveFunctionDetail();
}


double Embedding::moveStepSize(int method, double lastOF) {
	//clock_t start = clock();
	switch (method) {
	case 0: {
				// Constant Move, can not be used in practice for 0.01 will 
				// break the descent condition in later iterations.
				this->iterateVectors(0.01);
				return this->objectiveFunction();
	}
	case 1: {
				//backtracking line search holding Armijo rule 
				double t = 1;   //initial step size
				double pg = 0;
				for (int i = 0; i < this->nodeCount; i++) {
					pg += this->vectors[i].getVectorSquareSum(2);
				}
				this->iterateVectors(t);
				while (true) {
					double curOF = this->objectiveFunction();
					if (curOF < lastOF - alpha * t * pg) {
						return curOF;
					}
					else {
						t *= beta;
						this->iterateVectors(t);
						if (t < 1e-12) {
							this->iterateVectors(0);
							return lastOF;
						}
					}
				}
				break;
	}
	default: {
				 return 0;
	}
	}
}

void Embedding::calDescentDirection() {
	int edgeIndex = 0, nonedgeIndex = 0;
	int t;

	for (int i = 0; i < nodeCount; i++) {
		for (int j = 0; j < cmntCount; j++)
			tempDrct[j] = 0;
		while (edgeIndex < 2 * edgeCount && doubleEdge[edgeIndex].getStart() == i) {
			// for current edge:
			// 		direction = 2 * |Xs-Xt| * d(|Xs - Xt|)
			// 		|Xs - Xt| = sqrt((Xs1 - Xt1)^2 + ... + (Xsd - Xtd)^2)
			// 		(vectors with d dimension)
			t = doubleEdge[edgeIndex].getEnd();
			edgeIndex++;
			this->vectors[i].directionPlus(vectors[t].getVectorCurP(), 2, tempDrct);
		}
		while (nonedgeIndex < 2 * edgeCount && doubleNonedge[nonedgeIndex].getStart() == i) {
			 // so we will sampling a set of m non - edge relationships
			 //
			 // for current non - edge
			 // 		direction = 2 * (| Xs - Xt | -1) * d(| Xs - Xt | )
			 // | Xs - Xt | = sqrt((Xs1 - Xt1) ^ 2 + ... + (Xsd - Xtd) ^ 2)
			 //		(vectors with d dimension)
			 
			t = doubleNonedge[nonedgeIndex].getEnd();
			nonedgeIndex++;
			double lvd = vectors[i].getLengthOfMinus(vectors[t].getVectorNextP());
			if (lvd < 1e-13) {
				continue;
			}
			this->vectors[i].directionPlus(vectors[t].getVectorCurP(), 2 * (lvd - 1) * this->balance / lvd, tempDrct);
		}
		vectors[i].directionFinal(tempDrct, cmntCount);
	}
}

void Embedding::iterateVectors(double stepSize) {
	for (int i = 0; i < this->nodeCount; i++) {
		this->vectors[i].iterateVector(stepSize);
	}
}

double Embedding::objectiveFunction() {
	double of = 0;
	of += this->edges->getStress(this->vectors);
	of += this->nonedges->getStress(this->vectors);
	return of;	
}

void Embedding::copyNextToCurrent() {
	for (int i = 0; i < this->nodeCount; i++) {
		vectors[i].copyNextToCur();
	}
}

void Embedding::objectiveFunctionDetail() {
	string writeFile = "Objective Function Detail.txt";
	ofstream oFile;
	oFile.open(writeFile);
	oFile << "--------------Objective Function Deatils Start------------" << endl;
	this->edges->printStress(vectors, oFile);
	this->nonedges->printStress(vectors, oFile);

	oFile << "--------------Objective Function Deatils End------------" << endl;
	oFile.close();
}

void Embedding::PrintVectors() const {
	string writeFile = "embedding.txt";
	ofstream oFile;
	oFile.open(writeFile);
	for (int i = 0; i < nodeCount; i++) {
		oFile << i << ": ";
		vectors[i].printVector(oFile);
	}
	oFile.close();
}

int Embedding::getNodeCount() const {
	return nodeCount;
}

int Embedding::getEdgeCount() const {
	return edgeCount;
}

int Embedding::getCmtyCount() const {
	return cmntCount;
}

Edgeset* Embedding::getEdges() const {
	return edges;
}

Vector* Embedding::getVectors() const {
	return vectors;
}