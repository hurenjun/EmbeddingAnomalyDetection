
/*
 * Constructor of Edge will automaticly give both start and end -1.
 * Create an edge with start < end should use setPointWithSort, 
 * while create an edge with no constrains should use a setPoint after constructor.
 */
#pragma once

struct Edge
{
	int start, end;
	void intraSort() {
		if (this->start > this->end) {
			int temp = this->start;
			this->start = this->end;
			this->end = temp;
		}
	}

	Edge() {
		this->start = -1;
		this->end = -1;
	}
	Edge(const int s, const int t) {
		this->start = s;
		this->end = t;
		this->intraSort();
	}

	int getStart() const {
		return this->start;
	}
	int getEnd() const {
		return this->end;
	}

	void setPoint(const int start, const int end) {
		this->start = start;
		this->end = end;
	}
	void setPointWithSort(const int start, const int end) {
		this->start = start;
		this->end = end;
		this->intraSort();
	}

	bool operator < (const Edge & e) const {
		if ((start < e.start) || (start == e.start && end < e.end)) {
			return true;
		}
		else{
			return false;
		}
	}
	bool operator == (const Edge & e) const {
		if (start == e.start && end == e.end) {
			return true;
		}
		return false;
	}
}; 