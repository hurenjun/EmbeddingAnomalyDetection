// Class Vector implement a vector of each node. 
// Vector only maintain p (p << d) active dimensions, and other q alternate dimensions.
// 
// Main data structure in Vector:
//		vectorDrct: direction vector, a.k.a ¡®partial derivative¡¯ of objective function under current vectors
//		vectorCur: current vector				
//		vectorNext: next current, normalized vector of (vectorCur - stepSize * vectorDrct)
// More aboud direction vector:
//		we will implement two sections to preserve direction: 1) a set of (p + q), just same with 
//		(p +q) elements in vectorCur; 2) a heap with top (p + q) min elements 
// Vector has four common-used operations:
//		1. vector minus: (V_i - V_j) where V_i & V_j are vectors
//		2. update direction vector by (V_i - V_j)
//		3. update vectorNext by vectorCur - stepSize * vectorDrct
//		4. vector normalize (if necessary)
//		

#pragma once
#include <fstream>

#define MAX_P_SIZE			(24)
#define MAX_Q_SIZE			(6)
#define MAX_HEAP_SIZE		(30)

using namespace std;

struct VectorItem {
	int index;
	double weight;
	VectorItem(int index = -1, double weight = 0) {
		this->index = index;
		this->weight = weight;
	}

	void setValue(int index, double weight) {
		this->index = index;
		this->weight = weight;
	}

	void exchange(VectorItem & other) {
		int tempInt;
		tempInt = index; index = other.index; other.index = tempInt;
		double tempDbl;
		tempDbl = weight; weight = other.weight; other.weight = tempDbl;
	}

	bool operator < (const VectorItem &vi) const{
		if (weight > vi.weight)  {
			return true;
		}
		else{
			return false;
		}
	}
};

typedef struct Node {
	VectorItem vi;
	struct Node *next;
	
	Node(int index = -1, double weight = 0) {
		this->vi.setValue(index, weight);
		this->next = NULL;
	}
} LinkNode, *LinkList;

class Vector {
public:
	Vector();
	~Vector();
	void initialVector(int, int, int);
	void copyNextToCur();
	void directionPlus(VectorItem *, double, double *);
	void directionFinal(double *, const int &);
	void iterateVector(double);
	VectorItem* getVectorCurP();
	VectorItem* getVectorNextP();
	double getVectorSquareSum(int);
	double getLengthOfMinus(VectorItem *);
	void printVector(ofstream &);
	int IndexOfMaxW() const;
	int GetDomiDim(int *, double *) const;
private:	
	static const double sqrt2d2;
	// sort in index
	VectorItem vectorCurP[MAX_P_SIZE], vectorCurQ[MAX_Q_SIZE];
	VectorItem vectorNextP[MAX_P_SIZE], vectorNextQ[MAX_Q_SIZE];
	VectorItem vectorDrctP[MAX_P_SIZE], vectorDrctQ[MAX_Q_SIZE];
	LinkList vectorDrctOther = NULL;	
	LinkNode otherSet[MAX_HEAP_SIZE];
	double maxWeight;
	int p, q;
	void regularize();
	void adjustMaxHeapUpward(int, int, VectorItem *);
	void adjustMaxHeapDownward(int, int, VectorItem *);
	void adjustMinHeapUpward(int, int, VectorItem *);
	void adjustMinHeapDownward(int, int, VectorItem *);
	void swapOtherIndex(int, int);
	int findSmallerSon(int, int, VectorItem*);
	int findLargerSon(int, int, VectorItem *);
	void insertMinHeap(VectorItem*, int &, int, double);
	void insertMaxHeap(VectorItem*, int &, int, double);
};

