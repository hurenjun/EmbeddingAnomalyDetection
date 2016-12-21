#include "Graph.h"
#include "metis.h"
#include <iostream>
#include <fstream>
#include <algorithm>

Graph::Graph(int method, string &path, int cmtyCount) {
	this->cmtyCount = cmtyCount;
	loadGraph(path);

	switch (method) {
	case 0:
		randonPartition();
		break;
	case 1:
		heuristicPartition();
		break;
	case 2:
		metisPartititon();
		break;
	default:
		break;
	}
}

Graph::~Graph() {
	delete[] edges;
	edges = NULL;
	for (int i = 0; i < cmtyCount; i++) {
		delete[] cmties[i];
		cmties[i] = NULL;
	}
	delete[] cmties;
	cmties = NULL;
}

void Graph::loadGraph(string & path) {
	ifstream inFile;
	inFile.open(path);
	if (!inFile) {
		cout << "File Processing Error£¡" << endl;
		exit(1);
	}
	inFile >> nodeCount;
	inFile >> edgeCount;

	edges = new Edge[edgeCount];
	cmties = new int*[cmtyCount];

	int s, t;
	for (int i = 0; i < edgeCount; i++) {
		/*
		 * We recommend the node ID start from 0
		 * but in some input files, they start from 1, so we minus 1 when reading
		 *	minus: jure datasets, example dataset
		 *	no minus: synthetic dataset
		 */
		inFile >> s;
		inFile >> t;
		edges[i].setPointWithSort(s, t);
		//edges[i].setPointWithSort(s - 1, t - 1);
	}
	inFile.close();
}

void Graph::randonPartition() {
	int* nodeCmty = new int[nodeCount];
	int* cmtySize = new int[cmtyCount];
	memset(nodeCmty, 0, nodeCount * sizeof(int));
	memset(cmtySize, 0, cmtyCount * sizeof(int));	

	for (int i = 0; i < nodeCount; i++) {
		int currentnodeCmty = rand() % cmtyCount;
		nodeCmty[i] = currentnodeCmty;
		cmtySize[currentnodeCmty]++;
	}

	for (int i = 0; i < cmtyCount; i++) {
		cmties[i] = new int[cmtySize[i] + 1];
		//cout << "size: " << sizeof(cmties[i]) << "cmty size: " << cmtySize[i] << endl;
		cmties[i][0] = cmtySize[i];
		for (int j = 0, k = 1; j < nodeCount; j++) {
			if (nodeCmty[j] == i) {
				cmties[i][k++] = j;
			}
		}
	}

	delete[] nodeCmty;
	delete[] cmtySize;
}

int Graph::findCmty(int k, int nodeCmty[]) {
	if (k == nodeCmty[k]) {
		return k;
	}
	nodeCmty[k] = findCmty(nodeCmty[k], nodeCmty);
	return nodeCmty[k];
}

void Graph::heuristicPartition() {
	//System.out.println("Start Partition");	
	// nodeUsed: for each node, 0 not used; 1 used
	// nodeCmty: node's community ID
	int interGroupEdge = edgeCount + 1;
	int *bestNodeCmty = new int[nodeCount];
	int *nodeUsed = new int[nodeCount];
	int *nodeCmty = new int[nodeCount];
	int *cmtySize = new int[nodeCount];

	for (int ite = 0; ite < 5; ite++) {		
		// nodeRest: the number of nodes that haven't be merged into a community
		// cmtyRest: the number of community rest
		int nodeRest = nodeCount;
		int cmtyRest = cmtyCount;
		
		memset(nodeUsed, 0, nodeCount * sizeof(int));

		int round = 0;
		while (cmtyRest > 1) {
			/*
			* Initialize each node i's community ID by i.
			*
			* An edge's presence will merge two communities(A, B) into one,
			* through changing IDs of nodes in community B to IDs of nodes
			* in community A.
			*/
			for (int i = 0; i < nodeCount; i++) {
				if (nodeUsed[i] == 0) {
					nodeCmty[i] = i;
					cmtySize[i] = 1;
				}
			}

			this->randomEdgeset();

			int edgeIndex;
			for (edgeIndex = 0; edgeIndex < edgeCount; edgeIndex++) {
				int start = edges[edgeIndex].getStart();
				int end = edges[edgeIndex].getEnd();
				if (nodeUsed[start] == 1 || nodeUsed[end] == 1) {
					continue;
				}

				int startCmty = findCmty(start, nodeCmty);
				int endCmty = findCmty(end, nodeCmty);
				if (startCmty == endCmty) {
					continue;
				}
				nodeCmty[endCmty] = startCmty;
				cmtySize[startCmty] += cmtySize[endCmty];
				cmtySize[endCmty] = 0;

				if (cmtySize[startCmty] >= nodeRest / cmtyRest) {
					nodeRest -= cmtySize[startCmty];
					cmtyRest--;
					for (int j = 0; j < nodeCount; j++) {
						if (findCmty(j, nodeCmty) == startCmty) {
							nodeUsed[j] = 1;
						}
					}
					break;
				}
			}

			/*
			* remove the largest component when no community >= threshold
			*/
			if (edgeIndex == edgeCount) {
				int maxcmtySize = 0;
				int maxcmtyIndex = -1;
				for (int i = 0; i < nodeCount; i++) {
					int cmty = findCmty(i, nodeCmty);
					if (nodeUsed[i] == 0 && cmtySize[cmty] > maxcmtySize) {
						maxcmtySize = cmtySize[cmty];
						maxcmtyIndex = cmty;
					}
				}
				nodeRest -= maxcmtySize;
				cmtyRest--;
				for (int i = 0; i < nodeCount; i++) {
					if (findCmty(i, nodeCmty) == maxcmtyIndex) {
						nodeUsed[i] = 1;
					}
				}
			}
		}

		/*
		* merge the rest nodes into last community
		*/
		int lastcmtyID = -1;
		for (int i = 0; i < nodeCount; i++) {
			if (nodeUsed[i] == 0) {
				lastcmtyID = i;
				break;
			}
		}
		for (int i = 0; i < this->nodeCount; i++) {
			if (nodeUsed[i] == 0) {
				nodeCmty[i] = lastcmtyID;
			}
		}

		int inter = getInterGroupEdgeCount(nodeCmty);
		if (inter < interGroupEdge) {
			interGroupEdge = inter;
			for (int i = 0; i < nodeCount; i++)
				bestNodeCmty[i] = nodeCmty[i];
		}
	}

	for (int i = 0; i < nodeCount; i++)
		nodeCmty[i] = bestNodeCmty[i];
	/*
	* extract communities
	*/
	int* cmtyIDs = new int[cmtyCount];
	int tempIndex = 0;
	for (int i = 0; i < this->nodeCount; i++) {
		int j;
		for (j = 0; j < tempIndex; j++) {
			if (nodeCmty[i] == cmtyIDs[j]) {
				break;
			}
		}
		if (j == tempIndex) {
			cmtyIDs[tempIndex++] = nodeCmty[i];
		}
	}
	for (int i = 0; i < this->cmtyCount; i++) {
		int cmtySize = 0;
		for (int j = 0; j < this->nodeCount; j++) {
			if (nodeCmty[j] == cmtyIDs[i]) {
				cmtySize++;
			}
		}
		cmties[i] = new int[cmtySize + 1];
		cmties[i][0] = cmtySize;
		for (int j = 0, k = 1; j < this->nodeCount; j++) {
			if (nodeCmty[j] == cmtyIDs[i]) {
				cmties[i][k++] = j;
			}
		}
	}

	delete[] nodeUsed;
	delete[] nodeCmty;
	delete[] cmtySize;
	delete[] cmtyIDs;
	delete[] bestNodeCmty;
}

int Graph::getInterGroupEdgeCount(int *nodeCmty) {
	int result = 0;
	for (int i = 0; i < edgeCount; i++) {
		int s = edges[i].getStart();
		int t = edges[i].getEnd();
		if (findCmty(s, nodeCmty) != findCmty(t, nodeCmty)) result++;
	}
	return result;
}

void Graph::metisPartititon() {
	Edge *doubleEdge = new Edge[edgeCount * 2];
	for (int i = 0; i < edgeCount; i++) {
		doubleEdge[i * 2].setPoint(edges[i].getStart(), edges[i].getEnd());
		doubleEdge[i * 2 + 1].setPoint(edges[i].getEnd(), edges[i].getStart());
	}
	sort(doubleEdge, doubleEdge + 2 * edgeCount);

	idx_t nvtxs = nodeCount, ncon = 1, nparts = cmtyCount, objval;
	idx_t *xadj = new idx_t[nodeCount + 1];
	idx_t *adjncy = new idx_t[edgeCount * 2];
	idx_t *part = new idx_t[nodeCount];

	int edgeIndex = 0;
	for (int i = 0; i < nodeCount; i++) {
		xadj[i] = edgeIndex;
		while (edgeIndex < 2 * edgeCount && doubleEdge[edgeIndex].getStart() == i) {
			adjncy[edgeIndex] = doubleEdge[edgeIndex].getEnd();
			edgeIndex++;
		}
	}
	xadj[nodeCount] = edgeIndex;
	int metisResult = METIS_PartGraphKway(&nvtxs, &ncon, xadj, adjncy,
		NULL, NULL, NULL, &nparts, NULL,
		NULL, NULL, &objval, part);

	if (metisResult != METIS_OK) {
		cout << "METIS ERROR!";
		delete[] xadj;
		delete[] adjncy;
		delete[] part;
		delete[] doubleEdge;
		exit(1);
	}

	int* cmtySize = new int[cmtyCount];
	memset(cmtySize, 0, sizeof(int)* cmtyCount);

	for (int i = 0; i < nodeCount; i++)
		cmtySize[part[i]]++;
	for (int i = 0; i < cmtyCount; i++) {
		cmties[i] = new int[cmtySize[i] + 1];
		cmties[i][0] = cmtySize[i];
		cmtySize[i] = 1;
	}
	for (int i = 0; i < nodeCount; i++) {
		//cmties[i] = new int[cmtySize[i] + 1];
		//cout << "size: " << sizeof(cmties[i]) << "cmty size: " << cmtySize[i] << endl;
		//cmties[i][0] = cmtySize[i];
		//for (int j = 0, k = 1; j < nodeCount; j++) {
		//	if (part[j] == i) {
		//		cmties[i][k++] = j;
		//	}
		//}
		cmties[part[i]][cmtySize[part[i]]] = i;
		cmtySize[part[i]]++;
	}

	delete[] xadj;
	delete[] adjncy;
	delete[] part;
	delete[] doubleEdge;
	delete[] cmtySize;
}

void Graph::randomEdgeset() {
	for (int i = 0; i < this->edgeCount; i++) {
		int j = rand() % edgeCount;
		int tempStart = edges[i].getStart();
		int tempEnd = edges[i].getEnd();
		edges[i].setPoint(this->edges[j].getStart(), this->edges[j].getEnd());
		edges[j].setPoint(tempStart, tempEnd);
	}
}

void Graph::printCmties() {
	for (int i = 0; i < this->cmtyCount; i++) {
		cout << "Community " <<  i << ": ";
		for (int j = 1; j <= cmties[i][0]; j++) {
			cout<< this->cmties[i][j] << " ";
			if (j % 100 == 0) {
				cout << endl;
			}
		}
		cout << endl;
	}
}

int Graph::getNodeCount() const{
	return nodeCount;
}

int Graph::getEdgeCount() const{
	return edgeCount;
}

int Graph::getCmtyCount() const{
	return cmtyCount;
}

Edge* Graph::getEdges() const{
	return edges;
}

int** Graph::getCmties() const{
	return cmties;	
}