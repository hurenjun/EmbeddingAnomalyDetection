// Class AnomalyDetection load the embedding vectors of a graph
// and then use a specific metric for defined anomaly detection
//
// Matric Z-value > 3:
//		compute the z-value of all node-stresses, and labelled nodes whose 
//		node-stress z-value > 0 as anomalies

#include "Embedding.h"
#include "Edgeset.h"
#include "Vector.h"
#include <set>

struct Sorted_Dou_Int_Pair {
	double d;
	int i;
	void SetValue(const int i_i, const double i_d) {
		i = i_i;
		d = i_d;
	}
	bool operator < (const Sorted_Dou_Int_Pair & i_other) const {	// > actually
		if (d > i_other.d)
			return true;
		return false;
	}
};

struct IID {
	int s, t; double d;
	void SetValue(const int ss, const int tt, const double dd) {
		s = ss; t = tt; d = dd;
	}
	bool operator < (const IID & i_other) const {	// > actually
		if (d > i_other.d)
			return true;
		return false;
	}
};

class AnomalyDetection
{
public:
	AnomalyDetection(Embedding *, const double);
	~AnomalyDetection();
	void GetPreRecCurve(const string);
	double GetF1(const string, const double);
	void TopAnomalies(const string, const double);
	void TopK(const int);
	void RemoveEdge(const string, const string);
	void RemoveNode(const string, const string, const double);
	void GraphRewrite(const string, const string, const double);
private:
	int nodeCount, edgeCount, cmntCount;
	double theta;
	// create own edgeset for anomaly detection
	// and only copy a pointer of vectors
	Edgeset *edges = NULL;
	Vector *vectors = NULL;
	double *nodeStress = NULL;
	double * avg_node_stress = NULL;
	int *deg, *link_comm, *larg_comm;
	int ** node_comm = NULL;
	double **node_comm_w = NULL;
	int * ncc = NULL;
	double avgStress, stdDev, min_total, min_avg, max_total, max_avg;
	Sorted_Dou_Int_Pair * anomaly_deg = NULL;
	int stress_dis[105];
	void detectAnomaly();
	double getNodeStress(int);
	int getEdgeIndex(int, int) const;
	void PrintDetailNodeStress(const int, ofstream &);
	void GetNodeComm();
	double LargNeiComm(const int) const;
	void PrintLargNeiComm(const int, ofstream &) const;
	double LenOfVD(const int, const int) const;
};
