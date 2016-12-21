#include "Edgeset.h"
#include "Vector.h"
#include <algorithm>
#include <iostream>

//constructor of edge set
Edgeset::Edgeset(const Graph & g) {
	this->size = g.getEdgeCount();
	this->nodeCount = g.getNodeCount();
	this->cmtyCount = g.getCmtyCount();
	this->edgeCount = g.getEdgeCount();
	this->set = new Edge[this->size];
	this->copyEdge(g.getEdges());
	sort(set, set + size);
}

//constructor of double edge set
Edgeset::Edgeset(const Edgeset *single) {
	/*
	* This constructor 
	* regard the network as a directed one, so double edges and sort it
	* it will be quite convenient when extracting edges that connecting with some node
	*
	*  we generate Edge e here by 'setPoint(s, t)' method which do not have a intra-sort process
	*/
	this->size = single->size * 2;
	this->nodeCount = single->nodeCount;
	this->cmtyCount = single->cmtyCount;
	this->edgeCount = single->edgeCount;
	this->set = new Edge[this->size];

	for (int i = 0; i < single->size; i++) {
		this->set[i * 2].setPoint(single->set[i].getStart(), single->set[i].getEnd());
		this->set[i * 2 + 1].setPoint(single->set[i].getEnd(), single->set[i].getStart());
	}
	sort(set, set + size);
}

//constructor of non-edge set
Edgeset::Edgeset(const Edgeset & edges, int ** comm, int nodeCount, int edgeCount, int cmtyCount) {
	this->type = 1;
	this->nodeCount = nodeCount;
	this->edgeCount = edgeCount;
	this->cmtyCount = cmtyCount;
	//this->sampling = (this->nodeCount <= 1000) ? false : true;
	this->sampling = true;
	if (!this->sampling) {
		this->size = this->nodeCount * (this->nodeCount - 1) / 2 - edges.size;
	}
	else {
		this->size = this->getSamplingSize();
	}

	this->set = new Edge[this->size];

	if (this->sampling) {
		this->randomNonedgeSet(edges, comm);
	}
	else {
		this->getFullNonedgeSet(edges);
	}

	this->balance = edges.size * 1.0 / this->size;
}

Edgeset::~Edgeset() {
	if (set != NULL) {
		delete[] set;
		set = NULL;
	}	
}

double Edgeset::getStress(Vector *vectors) const {
	double result = 0, len;
	int s, t;
	int *x = NULL; double *y = NULL;

	switch (this->type) {
	case 0:
		for (int i = 0; i < this->size; i++) {
			s = this->set[i].getStart();
			t = this->set[i].getEnd();
			len = vectors[s].getLengthOfMinus(vectors[t].getVectorNextP());
			result += len * len;
		}
		break;
	case 1:
		for (int i = 0; i < this->size; i++) {
			s = this->set[i].getStart();
			t = this->set[i].getEnd();
			len = vectors[s].getLengthOfMinus(vectors[t].getVectorNextP()) - 1;
			result += this->balance * len * len;
		}
		break;
	}

	return result;
}

void Edgeset::printStress(Vector *vectors, ofstream &oFile) {
	int s, t;
	double len;
	switch (this->type) {
	case 0:
		for (int i = 0; i < this->size; i++) {
			s = this->set[i].getStart();
			t = this->set[i].getEnd();
			//lvd: length of vectors Difference			 
			len = vectors[s].getLengthOfMinus(vectors[t].getVectorNextP());
			oFile << s << " " << t << " " << len * len << endl;
			if (i % 100 == 9) {
				oFile.flush();
			}
		}
		break;
	case 1:
		for (int i = 0; i < this->size; i++) {
			s = this->set[i].getStart();
			t = this->set[i].getEnd();
			//lvd: length of vectors Difference			 
			double len = vectors[s].getLengthOfMinus(vectors[t].getVectorNextP()) - 1;
			oFile << s << " " << t << " " << this->balance * len * len << endl;
			if (i % 100 == 9) {
				oFile.flush();
			}
		}
		break;
	}
}

Edge Edgeset::getEdge(int k) const{
	return this->set[k];
}

int Edgeset::getSize() const{
	return this->size;
}

double Edgeset::getBalance() const{
	return this->balance;
}

void Edgeset::copyEdge(Edge *src) {
	for (int i = 0; i < size; i++) {
		this->set[i].setPointWithSort(src[i].getStart(), src[i].getEnd());
	}
}

void Edgeset::randomNonedgeSet(const Edgeset & edges, int ** comm) {
	for (int i = 0; i < this->size;) {
		//int k = (rand() << 16 | rand()) % cmtyCount;
		//int s = (rand() << 16 | rand()) % comm[k][0]; s = comm[k][s + 1];
		//int t = (rand() << 16 | rand()) % comm[k][0]; t = comm[k][t + 1];
		int s = (rand() << 16 | rand()) % nodeCount;
		int t = (rand() << 16 | rand()) % nodeCount;
		if (s == t) {
			continue;
		}
		if (!edges.containEdge(s, t)) {
			this->set[i].setPointWithSort(s, t);
			i++;
		}
	}
}

void Edgeset::getFullNonedgeSet(const Edgeset & edges) {
	int nonedgeIndex = 0;
	for (int s = 0; s < this->nodeCount; s++) {
		for (int t = s + 1; t < this->nodeCount; t++) {
			if (!edges.containEdge(s, t)) {
				this->set[nonedgeIndex++].setPoint(s, t);
			}
		}
	}
}

bool Edgeset::containEdge(int s, int t) const{
	int left = 0;
	int right = this->size - 1;
	int mid;
	Edge e(s, t);
	while (left <= right) {
		mid = (left + right) >> 1;
		if (this->set[mid] == e) {
			return true;
		}
		else if (this->set[mid] < e) {
			left = mid + 1;
		}
		else {
			right = mid - 1;
		}
	}
	return false;
}

int Edgeset::getSamplingSize() {
	return edgeCount;
}

void Edgeset::printSet() {
	for (int i = 0; i < size; i++) {
		std::cout << set[i].getStart() << " " << set[i].getEnd() << endl;
	}
}