#include "Vector.h"
#include <math.h>
#include <iostream>
#include <climits>
#include <algorithm>

const double Vector::sqrt2d2 = sqrt(2) * 0.5;

Vector::Vector() {
}

Vector::~Vector() {
}

void Vector::initialVector(int index, int p, int q) {
	this->p = p;
	this->q = q;
	vectorCurP[0].setValue(index, sqrt2d2);
	vectorCurP[1].index = -1;
	vectorNextP[0].setValue(index, sqrt2d2);
	vectorNextP[1].index = -1;
	vectorCurQ[0].index = -1;
	vectorNextQ[0].index = -1;
}

// vectorDrct += (x - y) * coe
// x = vectorCur
// y is also vectorCur of another node
void Vector::directionPlus(VectorItem * y, double coe, double *tempDrct) {
	VectorItem *x = vectorCurP;
	for (int i = 0; i < p; i++) {
		if (x[i].index == -1)	break;
		tempDrct[x[i].index] += x[i].weight * coe;
	}
	for (int i = 0; i < p; i++) {
		if (y[i].index == -1)	break;
		tempDrct[y[i].index] -= y[i].weight * coe;
	}
}

void Vector::directionFinal(double *tempDrct, const int &cmtySize) {
	int i;
	for (i = 0; i < p + 1; i++) {
		if (vectorCurP[i].index == -1)	break;
		vectorDrctP[i].setValue(vectorCurP[i].index, tempDrct[vectorCurP[i].index]);
		tempDrct[vectorCurP[i].index] = INT_MAX;		
	}
	vectorDrctP[i].index = -1;
	for (i = 0; i < q + 1; i++) {
		if (vectorCurQ[i].index == -1)	break;
		vectorDrctQ[i].setValue(vectorCurQ[i].index, tempDrct[vectorCurQ[i].index]);
		tempDrct[vectorCurQ[i].index] = INT_MAX;		
	}
	vectorDrctQ[i].index = -1;
	VectorItem heap[MAX_HEAP_SIZE];
	int heapSize = 0;
	for (i = 0; i < cmtySize; i++) {
		if (tempDrct[i] < 0)	insertMaxHeap(heap, heapSize, i, tempDrct[i]);
	}
	if (heapSize == 0)	vectorDrctOther = NULL;
	else {
		otherSet[0].vi.setValue(heap[0].index, heap[0].weight);
		otherSet[0].next = NULL;
		vectorDrctOther = &otherSet[0];
	}
	for (i = 1; i < heapSize; i++) {
		otherSet[i].vi.setValue(heap[i].index, heap[i].weight);
		otherSet[i].next = NULL;
		otherSet[i - 1].next = &otherSet[i];
	}	
}

//vectorNext = vectorCur - VectorDrct * stepSize
void Vector::iterateVector(double stepSize) {
	// get the top (p + q) max elememts
	// use a min heap
	VectorItem heap[MAX_HEAP_SIZE];
	int heapSize = 0;
	int index;
	double weight;

	for (int i = 0; i < p; i++) {
		if (vectorCurP[i].index == -1) {
			break;
		}
		index = vectorCurP[i].index;
		weight = vectorCurP[i].weight - vectorDrctP[i].weight * stepSize;
		insertMinHeap(heap, heapSize, index, weight);
	}
	for (int i = 0; i < q; i++) {
		if (vectorCurQ[i].index == -1) {
			break;
		}
		index = vectorCurQ[i].index;
		weight = vectorCurQ[i].weight - vectorDrctQ[i].weight * stepSize;
		insertMinHeap(heap, heapSize, index, weight);
	}

	LinkList l = vectorDrctOther;
	while (l != NULL) {
		insertMinHeap(heap, heapSize, l->vi.index, -l->vi.weight * stepSize);
		l = l->next;
	}

	//for (int i = 0; i < otherSize; i++) {
	//	insertMinHeap(heap, heapSize, vectorDrctOther[i].index, -vectorDrctOther[i].weight * stepSize);
	//}

	// sort by weight
	sort(heap, heap + heapSize);

	for (int i = 0; i < min(p, heapSize); i++) {
		for (int j = min(p, heapSize) - 1; j > i; j--) {
			if (heap[j - 1].index > heap[j].index) {
				heap[j].exchange(heap[j - 1]);
			}
		}
	}

	for (int i = p; i < heapSize; i++) {
		for (int j = heapSize - 1; j > i; j--) {
			if (heap[j - 1].index > heap[j].index) {
				heap[j].exchange(heap[j - 1]);
			}
		}
	}
	int i;
	for (i = 0; i < min(p, heapSize); i++) {
		vectorNextP[i].setValue(heap[i].index, heap[i].weight);
	}
	vectorNextP[i].index = -1;

	for (i = p; i < heapSize; i++) {
		vectorNextQ[i - p].setValue(heap[i].index, heap[i].weight);
	}
	vectorNextQ[i - p].index = -1;

	this->regularize();
}

void Vector::copyNextToCur() {
	int i;
	for (i = 0; i < p; i++) {
		if (vectorNextP[i].index == -1) {
			break;
		}
		vectorCurP[i].setValue(vectorNextP[i].index, vectorNextP[i].weight);
	}
	vectorCurP[i].index = -1;
	for (i = 0; i < q; i++) {
		if (vectorNextQ[i].index == -1) {
			break;
		}
		vectorCurQ[i].setValue(vectorNextQ[i].index, vectorNextQ[i].weight);
	}
	vectorCurQ[i].index = -1;
}

void Vector::printVector(ofstream &oFile) {
	double avg = 0;
	int k = 0;
	for (int i = 0; i < p; i++) {
		if (vectorCurP[i].index == -1) {
			break;
		}
		//oFile << '(' << vectorCurP[i].index << ", " << vectorCurP[i].weight << ") ";
		avg += vectorCurP[i].weight;
		k++;
	}
	avg /= k;
	oFile << "avg=" << avg << "\t";
	for (int i = 0; i < p; i++) {
		if (vectorCurP[i].index == -1) break;
		if (vectorCurP[i].weight > avg) oFile << '(' << vectorCurP[i].index << ", " << vectorCurP[i].weight << ") ";
	}
	oFile << endl;
	oFile.flush();
}

VectorItem* Vector::getVectorCurP() {
	return vectorCurP;
}

VectorItem * Vector::getVectorNextP() {
	return vectorNextP;
}

double Vector::getVectorSquareSum(int type) {
	double sum = 0;
	if (type == 0) {
		for (int i = 0; i < p; i++) {
			if (vectorNextP[i].index == -1) {
				break;
			}
			sum += vectorCurP[i].weight * vectorCurP[i].weight;
		}
	}
	else if (type == 1){
		for (int i = 0; i < p; i++) {
			if (vectorNextP[i].index == -1) {
				break;
			}
			sum += vectorNextP[i].weight * vectorNextP[i].weight;
		}
	}
	else {
		for (int i = 0; i < p; i++) {
			if (vectorDrctP[i].index == -1) {
				break;
			}
			sum += vectorDrctP[i].weight * vectorDrctP[i].weight;
		}
		for (int i = 0; i < q; i++) {
			if (vectorDrctQ[i].index == -1) {
				break;
			}
			sum += vectorDrctQ[i].weight * vectorDrctQ[i].weight;
		}
		LinkList ll = vectorDrctOther;
		while (ll != NULL) {
			sum += ll->vi.weight * ll->vi.weight;
			ll = ll->next;
		}
	}
	return sum;
}

// return the length of vector(vextorNextP - y);
double Vector::getLengthOfMinus(VectorItem *y) {
	double sum = 0;
	VectorItem *x = vectorNextP;
	int i, j;
	for (i = 0, j = 0; x[i].index != -1 && y[j].index != -1;) {
		if (x[i].index == y[j].index) {
			sum += (x[i].weight - y[j].weight) * (x[i].weight - y[j].weight);
			++i;
			++j;
		}
		else if (x[i].index < y[j].index) {
			sum += (x[i].weight)* (x[i].weight);
			++i;
		}
		else {
			sum += (y[j].weight)* (y[j].weight);
			++j;
		}
	}
	for (; x[i].index != -1; ++i) {
		sum += (x[i].weight)* (x[i].weight);
	}
	for (; y[j].index != -1; ++j) {
		sum += (y[j].weight)* (y[j].weight);
	}

	return sqrt(sum);
}

// heap vectorDrctOther and DrctOtherIndex are corressponding ids
// adjust heap upward from pos 
// para. pos: index of array, which has changed
void Vector::adjustMaxHeapUpward(int pos, int size, VectorItem *heap) {
	int father = (pos - 1) / 2;
	while (pos != 0 && heap[pos].weight > heap[father].weight) {
		// need to adjust
		// exchange pos and father
		heap[pos].exchange(heap[father]);
		//swapOtherIndex(heap[pos].index, heap[father].index);
		pos = father;
		father = (pos - 1) / 2;
	}
}

// heap vectorDrctOther and DrctOtherIndex are corressponding ids
// adjust heap download from pos
// para. pos: index of array, which has changed
void Vector::adjustMaxHeapDownward(int pos, int size, VectorItem *heap) {
	int son = findLargerSon(pos, size, heap);
	while (son != -1) {
		heap[pos].exchange(heap[son]);
		//swapOtherIndex(heap[pos].index, heap[son].index);
		pos = son;
		son = findLargerSon(pos, size, heap);
	}
}

void Vector::adjustMinHeapUpward(int pos, int size, VectorItem *heap) {
	int father = (pos - 1) / 2;
	while (pos != 0 && heap[pos].weight < heap[father].weight) {
		heap[pos].exchange(heap[father]);
		pos = father;
		father = (pos - 1) / 2;
	}
}

void Vector::adjustMinHeapDownward(int pos, int size, VectorItem *heap) {
	int son = findSmallerSon(pos, size, heap);
	while (son != -1) {
		heap[pos].exchange(heap[son]);
		pos = son;
		son = findSmallerSon(pos, size, heap);
	}
}

// exchange the value of follosing two keys
// MyMap otherIndex
void Vector::swapOtherIndex(int id1, int id2) {
	//int temp = otherIndex[id1];
	//otherIndex[id1] = otherIndex[id2];
	//otherIndex[id2] = temp;
}

// find a larger son of current node
// return -1 of not exist
// para. pos: current node 
int Vector::findLargerSon(int pos, int size, VectorItem *heap) {
	int leftSon = 2 * pos + 1;
	int rightSon = 2 * pos + 2;
	double leftSonWeight = leftSon <size ? heap[leftSon].weight : -INT_MAX;
	double rightSonWeight = rightSon < size ? heap[rightSon].weight : -INT_MAX;
	if (leftSonWeight > heap[pos].weight && leftSonWeight >= rightSonWeight) {
		return leftSon;
	}
	if (rightSonWeight > heap[pos].weight && rightSonWeight > leftSonWeight) {
		return rightSon;
	}
	return -1;
}

int Vector::findSmallerSon(int pos, int size, VectorItem* heap) {
	int leftSon = 2 * pos + 1;
	int rightSon = 2 * pos + 2;
	double leftSonWeight = leftSon < size ? heap[leftSon].weight : INT_MAX;
	double rightSonWeight = rightSon < size ? heap[rightSon].weight : INT_MAX;
	if (leftSonWeight < heap[pos].weight && leftSonWeight <= rightSonWeight) {
		return leftSon;
	}
	if (rightSonWeight < heap[pos].weight && rightSonWeight < leftSonWeight) {
		return rightSon;
	}
	return -1;
}

// insert VectotItem(index, weight) into maxheap heap
// size is the current size of maxheap
// may do not insert
void Vector::insertMinHeap(VectorItem* heap, int & size, int index, double weight) {
	if (weight <= 0) {		// remove element that is less than 0
		return;
	}
	if (size < p + q) {		// heap is not full, add in current VectorItem directly
		heap[size].index = index;
		heap[size].weight = weight;
		size++;
		adjustMinHeapUpward(size - 1, size, heap);
	}
	else if (weight > heap[0].weight){		//heap is full, and current VectorItem > min VectorItem in heap
		heap[0].index = index;
		heap[0].weight = weight;
		adjustMinHeapDownward(0, size, heap);
	}
}

void Vector::insertMaxHeap(VectorItem *heap, int &size, int index, double weight) {
	if (size < p + q) {
		heap[size].index = index;
		heap[size].weight = weight;
		size++;
		adjustMaxHeapUpward(size - 1, size, heap);
	}
	else if (weight < heap[0].weight) {
		heap[0].index = index;
		heap[0].weight = weight;
		adjustMaxHeapDownward(0, size, heap);
	}
}

// Regularize criteria: vector length no longer than sqrt(2) / 2
// The vectot that will be regularized is NextP, not CurP
void Vector::regularize() {
	double length = sqrt(getVectorSquareSum(1));
	if (length > sqrt2d2) {
		double coe = sqrt2d2 / length;

		for (int i = 0; i < p; i++) {
			if (vectorNextP[i].index == -1) {
				break;
			}
			vectorNextP[i].weight *= coe;
		}
		for (int i = 0; i < q; i++) {
			if (vectorNextQ[i].index == -1) {
				break;
			}
			vectorNextQ[i].weight *= coe;
		}
	}
}

int Vector::IndexOfMaxW() const {
	int index = vectorCurP[0].index;
	double maxW = vectorCurP[0].weight;
	
	for (int i = 1; i < p; i++) {
		if (vectorCurP[i].index == -1) break;
		if (vectorCurP[i].weight > maxW) {
			maxW = vectorCurP[i].weight;
			index = vectorCurP[i].index;
		}
	}
	return index;
}

int Vector::GetDomiDim(int * dim, double * w) const {
	double avg = 0;
	int k = 0;
	for (int i = 0; i < p; i++) {
		if (vectorCurP[i].index == -1) {
			break;
		}
		avg += vectorCurP[i].weight;
		k++;
	}
	avg /= k;
	int n_dim = 0;
	for (int i = 0; i < p; i++) {
		if (vectorCurP[i].index == -1) break;
		if (vectorCurP[i].weight > avg) {
			dim[n_dim] = vectorCurP[i].index;
			w[n_dim] = vectorCurP[i].weight;
			n_dim++;
		}
	}
	return n_dim;
}