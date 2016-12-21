#include "AnomalyDetection.h"
#include "Edge.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <algorithm>
#include <set>
#include <map>

using namespace std;

AnomalyDetection::AnomalyDetection(Embedding *embd, const double i_theta) {
	this->nodeCount = embd->getNodeCount();
	this->edgeCount = embd->getEdgeCount() * 2;
	this->cmntCount = embd->getCmtyCount();
	this->theta = i_theta;
	this->vectors = embd->getVectors();
	this->edges = new Edgeset(embd->getEdges());
	this->nodeStress = new double[this->nodeCount];
	this->avg_node_stress = new double[this->nodeCount];
	this->deg = new int[this->nodeCount];
	for (int k = 0; k < nodeCount; k++) {
		int left = this->getEdgeIndex(k, 0);
		int right = this->getEdgeIndex(k, this->nodeCount - 1);
		deg[k] = right - left;
	}
	this->link_comm = new int[this->nodeCount];
	this->larg_comm = new int[this->nodeCount];
	this->anomaly_deg = new Sorted_Dou_Int_Pair[this->nodeCount];
	this->detectAnomaly();
	
	//this->printNodeStress();
}

AnomalyDetection::~AnomalyDetection() {
	delete edges;			edges = NULL;
	delete[] nodeStress;	nodeStress = NULL;
	delete[] avg_node_stress;	avg_node_stress = NULL;
	if (deg != NULL) { delete[] deg; deg = NULL; }
	if (link_comm != NULL)	 { delete[] link_comm; link_comm = NULL; }
	if (anomaly_deg != NULL) { delete[] anomaly_deg; anomaly_deg = NULL; }
	if (larg_comm != NULL) { delete[] larg_comm; larg_comm = NULL; }
	if (ncc != NULL) { delete[] ncc; ncc = NULL; }
	if (node_comm != NULL) {
		for (int i = 0; i < nodeCount; i++) { 
			delete[] node_comm[i]; 
			delete[] node_comm_w[i];
		}
		delete[] node_comm;
		delete[] node_comm_w;
	}
}

void AnomalyDetection::detectAnomaly() {
	//double totalStress = 0;
	//memset(stress_dis, 0, sizeof(int)* 100);
	//for (int i = 0; i < this->nodeCount; i++) {
	//	this->getNodeStress(i);
	//	totalStress += this->nodeStress[i];
	//}
	//avgStress = totalStress / this->nodeCount;
	
	//ofstream stre_dis_ofs("stress_dis.dat");
	//double cum_prob = 0;
	//for (int i = 0; i < 101; i++) {
	//	cum_prob += stress_dis[i] * 1.0 / edgeCount;
	//	stre_dis_ofs << "[" << i * 1.0 / 100 << ", " << (i + 1) * 1.0 / 100 << "]: " << stress_dis[i] / 2
	//		<< "\tprob=" << stress_dis[i] * 1.0 / edgeCount << "\tcum_prob=" << cum_prob << endl;
	//}
	//stre_dis_ofs.close();

	//stdDev = 0;
	//for (int i = 0; i < this->nodeCount; i++) {
	//	stdDev += (this->nodeStress[i] - avgStress) * (this->nodeStress[i] - avgStress);
	//}
	//stdDev = sqrt(stdDev / this->nodeCount);

	GetNodeComm();
	for (int i = 0; i < this->nodeCount; i++) {
		//if (nodeStress[i] > 1.5)
		//	//anomaly_deg[i].SetValue(i, abs(nodeStress[i] - avgStress) / stdDev);
		//	anomaly_deg[i].SetValue(i, avg_node_stress[i]);
		//else anomaly_deg[i].SetValue(i, 0);
		//anomaly_deg[i].SetValue(i, larg_comm[i] * 1.0 / deg[i]);
		anomaly_deg[i].SetValue(i, LargNeiComm(i));
	}
	sort(anomaly_deg, anomaly_deg + nodeCount);

	//double totalStress = 0;
	//for (int i = 0; i < nodeCount; i++) {
	//	totalStress += anomaly_deg[i].d;	
	//}
	//avgStress = totalStress / this->nodeCount;
	//stdDev = 0;
	//for (int i = 0; i < this->nodeCount; i++) {
	//	stdDev += (anomaly_deg[i].d - avgStress) * (anomaly_deg[i].d - avgStress);
	//}
	//stdDev = sqrt(stdDev / this->nodeCount);
	//cout << "AVG. " << avgStress << endl;
	//cout << "STD DEV " << stdDev << endl;
}

double AnomalyDetection::getNodeStress(int k) {
	double ns = 0;
	int left = this->getEdgeIndex(k, 0);
	int right = this->getEdgeIndex(k, this->nodeCount - 1);
	int s, t;
	double len;
	map<int, int> nb_comm;
	for (int i = left; i < right; i++) {
		s = this->edges->getEdge(i).getStart();
		t = this->edges->getEdge(i).getEnd();
		len = vectors[s].getLengthOfMinus(vectors[t].getVectorNextP());
		//if (len * len > 0.1) 
		stress_dis[(int)floor(len * len * 100)]++;
		ns += len * len;
		int comm = vectors[t].IndexOfMaxW();
		if (nb_comm.find(comm) == nb_comm.end()) nb_comm.insert(pair<int, int>(comm, 1));
		else nb_comm.find(comm)->second++;
	}
	this->nodeStress[k] = ns;
	this->avg_node_stress[k] = ns / (right - left);
	this->link_comm[k] = nb_comm.size();
	this->deg[k] = right - left;
	int max_comm = 0;
	for (map<int, int>::iterator iter = nb_comm.begin(); iter != nb_comm.end(); ++iter) {
		if (iter->second > max_comm) max_comm = iter->second;
	}
	this->larg_comm[k] = max_comm;
	//return ns;
	return ns / (right - left);
}

void AnomalyDetection::PrintDetailNodeStress(const int k, ofstream & ofs) {
	ofs << k << ": ";
	double ns = 0;
	int left = this->getEdgeIndex(k, 0);
	int right = this->getEdgeIndex(k, this->nodeCount - 1);
	int s, t;
	double len;
	map<int, int> nb_comm;
	for (int i = left; i < right; i++) {
		s = this->edges->getEdge(i).getStart();
		t = this->edges->getEdge(i).getEnd();
		len = vectors[s].getLengthOfMinus(vectors[t].getVectorNextP());
		//ofs << vectors[t].IndexOfMaxW() << " ";
		int index = vectors[t].IndexOfMaxW();
		map<int, int>::iterator iter = nb_comm.find(index);
		if (iter == nb_comm.end()) nb_comm.insert(pair<int, int>(index, 1));
		else iter->second++;
		//if (len * len > 0.1) {
			//ns += len * len;
			//ofs << "<" << t << "," << len * len << "> ";
		//}
	}

	//ofs << "avg=" << ns / (right - left) << " total=" << ns << endl;
	ofs << "total: " << nb_comm.size() << "\t";
	for (map<int, int>::iterator iter = nb_comm.begin(); iter != nb_comm.end(); ++iter) {
		ofs << "<" << iter->first << ", " << iter->second << "> ";
	}
	ofs << endl;
	min_total = min(ns, min_total);
	min_avg = min(ns / (right - left), min_avg);
	max_total = max(ns, max_total);
	max_avg = max(ns / (right - left), max_avg);
}

int AnomalyDetection::getEdgeIndex(int s, int t) const {
	Edge e(-1, -1);
	e.setPoint(s, t);
	if (e == this->edges->getEdge(this->edgeCount - 1)) {
		return this->edgeCount;
	}

	int left = 0;
	int right = this->edgeCount - 1;
	while (left < right) {
		int mid = (left + right) / 2;
		if (e == this->edges->getEdge(mid)) {
			return mid;
		}
		if (e < this->edges->getEdge(mid)) {
			right = mid;
		}
		else {
			left = mid + 1;
		}
	}
	return left;
}

void AnomalyDetection::GetPreRecCurve(const string filename) {
	fstream anomaly(filename + "-anomaly.txt");
	if (!anomaly) {
		cout << "can not open anomaly file: " + filename + "-anomaly.txt" << endl;
		for (int i = 0; i <= 10; i++) {
			cout << i << " " << anomaly_deg[i * 1000].d << endl;
		}
		return;
	}
	int anomaly_cnt = 0;
	set<int> anml;
	string line;
	while (anomaly) {
		getline(anomaly, line);
		if (line.length() == 0)
			break;
		anomaly_cnt++;
		anml.insert(stoi(line));
	}
	anomaly.close();

	int found = 0;
	//double rec = 0.05;
	int rec = anml.size() / 20;
	double f1 = 0;
	for (int i = 0; i < nodeCount; i++) {
		if (anml.find(anomaly_deg[i].i) != anml.end()) {
			found++;
			//if (found >= anomaly_cnt * rec) {
			if (found >=  rec) {
				double recall = rec * 1.0 / anomaly_cnt;
				double prec = found * 1.0 / (i + 1);
				cout << recall << " " << prec << " " 
					<< (anomaly_deg[i].d - avgStress) / stdDev << endl;
				rec += anml.size() / 20;
				f1 = max(f1, 2 * recall * prec / (recall + prec));
			}
		}
	}
	cout << "F1 Score: " << f1 << endl;
}

void AnomalyDetection::TopAnomalies(const string filename, const double theta) {
	fstream anomaly(filename + "-anomaly.txt");
	if (!anomaly) {
		cout << "can not open anomaly file: anomaly.txt" << endl;
	}
	int anomaly_cnt = 0;
	set<int> anml;
	string line;
	while (anomaly) {
		getline(anomaly, line);
		if (line.length() == 0)
			break;
		anomaly_cnt++;
		anml.insert(stoi(line));
	}
	anomaly.close();

	ofstream stress_ofs(filename + "-stress.txt");
	//stress_ofs << "hub nodes" << endl;
	//min_total = 10000;	min_avg = 10000;
	//max_total = 0; max_avg = 0;
	//for (set<int>::iterator iter = anml.begin(); iter != anml.end(); ++iter) {
	//	if (*iter < 100) {
	//		//PrintDetailNodeStress(*iter, stress_ofs);
	//		PrintLargNeiComm(*iter, stress_ofs);
	//	}
	//}
	////stress_ofs << "min_total=" << min_total << " min_avg=" << min_avg << endl;
	////stress_ofs << "max_total=" << max_total << " max_avg=" << max_avg << endl;
	//stress_ofs << "outlies" << endl;
	//min_total = 10000;	min_avg = 10000;
	//max_total = 0; max_avg = 0;
	//for (set<int>::iterator iter = anml.begin(); iter != anml.end(); ++iter) {
	//	if (*iter > 100) {
	//		//PrintDetailNodeStress(*iter, stress_ofs);
	//		PrintLargNeiComm(*iter, stress_ofs);
	//	}
	//}
	//stress_ofs << "min_total=" << min_total << " min_avg=" << min_avg << endl;
	//stress_ofs << "max_total=" << max_total << " max_avg=" << max_avg << endl;
	int k = 0;
	stress_ofs << "top 5000" << endl;
	min_total = 10000; 	min_avg = 10000;
	max_total = 0; max_avg = 0;
	for (int i = 0; i < nodeCount; i++) {
		if (anml.find(anomaly_deg[i].i) != anml.end()) {
			stress_ofs << "yes ";
		}
		k++;
			//PrintDetailNodeStress(anomaly_deg[i].i, stress_ofs);
		stress_ofs << anomaly_deg[i].d << " " << deg[anomaly_deg[i].i] << "; ";
		PrintLargNeiComm(anomaly_deg[i].i, stress_ofs);
		if (anomaly_deg[i].d < theta) break;
		
	}
	//stress_ofs << "min_total=" << min_total << " min_avg=" << min_avg << endl;
	//stress_ofs << "max_total=" << max_total << " max_avg=" << max_avg << endl;
	stress_ofs.close();
}

void AnomalyDetection::TopK(const int k) {
	ofstream top_ofs("topk.txt");
	for (int i = 0; i < k; i++) {
		top_ofs << anomaly_deg[i].i << " " << anomaly_deg[i].d << " "; // << endl;
		PrintLargNeiComm(anomaly_deg[i].i, top_ofs);
	}
	top_ofs.close();
}

void AnomalyDetection::GetNodeComm() {
	int * domi_dim = new int[cmntCount];
	double * domi_w = new double[cmntCount];
	ncc = new int[nodeCount];
	node_comm = new int*[nodeCount];
	node_comm_w = new double*[nodeCount];
	for (int i = 0; i < nodeCount; i++) {
		ncc[i] = vectors[i].GetDomiDim(domi_dim, domi_w);
		node_comm[i] = new int[ncc[i]];
		node_comm_w[i] = new double[ncc[i]];
		for (int j = 0; j < ncc[i]; j++) { 
			node_comm[i][j] = domi_dim[j]; 
			node_comm_w[i][j] = domi_w[j];
		}
	}
	
	delete[] domi_dim;
	delete[] domi_w;
}

double AnomalyDetection::LargNeiComm(const int s) const {
	int left = this->getEdgeIndex(s, 0);
	int right = this->getEdgeIndex(s, this->nodeCount - 1);
	//map<int, double> comm_pow;
	double * comm_pow = new double[cmntCount];
	//for (int i = 0; i < cmntCount; i++) { comm_pow[i] = 0; }
	memset(comm_pow, 0, cmntCount * sizeof(double));
	int t, tar_comm; double stren; 
	//map<int, double>::iterator iter;
	for (int i = left; i < right; i++) {
		t = edges->getEdge(i).getEnd();
		//stren = 1 - vectors[s].getLengthOfMinus(vectors[t].getVectorNextP());
		stren = 1 - LenOfVD(s, t);
		double total_w = 0;
		for (int j = 0; j < ncc[t]; j++) {
			total_w += node_comm_w[t][j];
		}
		for (int j = 0; j < ncc[t]; j++) {
			tar_comm = node_comm[t][j];
			comm_pow[tar_comm] += stren * node_comm_w[t][j] / total_w;
			//comm_pow[tar_comm] += stren * node_comm_w[t][j];
			//iter = comm_pow.find(tar_comm);
			//if (iter == comm_pow.end()) comm_pow.insert(pair<int, double>(tar_comm, stren * node_comm_w[t][j] / total_w));
			//else iter->second += stren * node_comm_w[t][j] / total_w;
		}
	}
	//double avg_pow = 0;
	//for (iter = comm_pow.begin(); iter != comm_pow.end(); ++iter) { avg_pow += iter->second; }
	//avg_pow /= comm_pow.size();
	//int valid_comm = 0;
	//for (iter = comm_pow.begin(); iter != comm_pow.end(); ++iter) { 
	//	if (iter->second > avg_pow) valid_comm++; 
	//}
	//return valid_comm;
	sort(comm_pow, comm_pow + cmntCount);
	//double base = comm_pow[cmntCount - 1];
	//double avg_pow = 0;
	//int valid = 0;
	//for (int i = 0; i < cmntCount; i++) {
	//	avg_pow += comm_pow[i];
	//	valid += (comm_pow[i] > 0) ? 1 : 0;
	//}
	//avg_pow /= valid;
	//double score = 0;
	//for (int i = cmntCount - 1; i >= 0; i--) {
	//	//if (comm_pow[i] >= base / 7)
	//	if (comm_pow[i] >= avg_pow)
	//		//score += log(1 + comm_pow[i]) / log(1 + base);
	//		score += (comm_pow[i]) / base;
	//}
	double max_comm = comm_pow[cmntCount - 1];
	double base = max_comm * theta;
	double score = 0;
	for (int i = cmntCount - 1; i >= 0; i--) {
		if (comm_pow[i] >= base) score += comm_pow[i] / max_comm;
	}
	delete[] comm_pow;
	return score;
}

void AnomalyDetection::PrintLargNeiComm(const int s, ofstream & ofs) const {
	ofs << s << ": ";
	int left = this->getEdgeIndex(s, 0);
	int right = this->getEdgeIndex(s, this->nodeCount - 1);
	map<int, double> comm_pow;
	int t, tar_comm; double stren; map<int, double>::iterator iter;
	for (int i = left; i < right; i++) {
		t = edges->getEdge(i).getEnd();
		//stren = 1 - vectors[s].getLengthOfMinus(vectors[t].getVectorNextP());
		stren = 1 - LenOfVD(s, t);
		double total_w = 0;
		for (int j = 0; j < ncc[t]; j++) {
			total_w += node_comm_w[t][j];
		}
		for (int j = 0; j < ncc[t]; j++) {
			tar_comm = node_comm[t][j];
			iter = comm_pow.find(tar_comm);
			if (iter == comm_pow.end()) comm_pow.insert(pair<int, double>(tar_comm, stren * node_comm_w[t][j] / total_w));
			else iter->second += stren * node_comm_w[t][j] / total_w;
			//if (iter == comm_pow.end()) comm_pow.insert(pair<int, double>(tar_comm, stren * node_comm_w[t][j]));
			//else iter->second += stren * node_comm_w[t][j];
		}
	}
	double max_pow = 0, avg_pow = 0;;
	for (iter = comm_pow.begin(); iter != comm_pow.end(); ++iter) {
		if (iter->second > max_pow) max_pow = iter->second;
		avg_pow += iter->second;
	}
	avg_pow /= comm_pow.size();
	double base = max_pow * theta;
	ofs << max_pow << " " << base << "; ";
	for (iter = comm_pow.begin(); iter != comm_pow.end(); ++iter) { 
		//if (iter->second >= max_pow / 7)
		//if (iter->second >= avg_pow)
		if (iter->second >= base)
			ofs << "<" << iter->first << ", " << iter->second << ">\t";
	}
	ofs << endl;
}

double AnomalyDetection::LenOfVD(const int s, const int t) const {
	double sum = 0;
	int i, j;
	for (i = 0, j = 0; i < ncc[s] && j < ncc[t];) {
		if (node_comm[s][i] == node_comm[t][j]) {
			sum += (node_comm_w[s][i] - node_comm_w[t][j]) * (node_comm_w[s][i] - node_comm_w[t][j]);
			++i;
			++j;
		}
		else if (node_comm[s][i] < node_comm[t][j]) {
			sum += (node_comm_w[s][i])* (node_comm_w[s][i]);
			++i;
		}
		else {
			sum += (node_comm_w[t][j])* (node_comm_w[t][j]);
			++j;
		}
	}
	for (; i < ncc[s]; ++i) {
		sum += (node_comm_w[s][i])* (node_comm_w[s][i]);
	}
	for (; j < ncc[t]; ++j) {
		sum += (node_comm_w[t][j])* (node_comm_w[t][j]);
	}

	return sqrt(sum);
}

void AnomalyDetection::RemoveEdge(const string source, const string target) {
	ifstream in(source);
	int n, m, s, t, mm;
	in >> n >> m;
	IID * edges = new IID[m];
	for (int i = 0; i < m; i++) {
		in >> s >> t;
		edges[i].SetValue(s, t, LenOfVD(s, t));
	}
	in.close();
	sort(edges, edges + m);
	double avg = 0;
	for (int i = 0; i < m; i++) { avg += edges[i].d; }
	avg /= m;
	double std_dev = 0;
	for (int i = 0; i < m; i++) { std_dev += (edges[i].d - avg) * (edges[i].d - avg); }
	std_dev = sqrt(std_dev / m);
	cout << "avg & std_dev: " << avg << " " << std_dev << endl;
	int start = 0;
	//for (int i = 0; i < 100; i++) {
	//	cout << edges[i].s << " " << edges[i].t << " " << edges[i].d << "\t";
	//}
	cout << endl;
	while ((edges[start].d - avg) / std_dev > 1.5) { start++; }
	cout << "start: " << start << endl;
	ofstream out(target);
	in >> n >> m;
	out << n << " " << m - start << endl;
	for (start; start < m; start++) {
		out << edges[start].s << " " << edges[start].t << endl;
	}
	out.close();
	delete[] edges; edges = NULL;
}

void AnomalyDetection::RemoveNode(const string source, const string target, const double thre) {
	//int * mapping = new int[nodeCount];
	//fill(mapping, mapping + nodeCount, 0);
	//for (int i = 0; i < nodeCount; i++) {
	//	if (anomaly_deg[i].d >= thre) mapping[anomaly_deg[i].i] = -1;
	//}
	//int valid = 0;
	//for (int i = 0; i < nodeCount; i++) {
	//	if (mapping[i] != -1) {
	//		mapping[i] = valid;
	//		valid++;
	//	}
	//}

	//ifstream in(source);
	//ofstream out(target);
	//int n, m, s, t;
	//in >> n >> m;
	//out << valid << endl;
	//cout << "valid node: " << valid << endl;
	//for (int i = 0; i < m; i++) {
	//	in >> s >> t;
	//	if (mapping[s] >= 0 && mapping[t] >= 0) out << mapping[s] << " " << mapping[t] << endl;
	//}
	//in.close();
	//out.close();
	//delete[] mapping; mapping = NULL;

	ofstream out(target);
	int anomaly = 0, valid = 0;
	for (int i = 0; i < nodeCount; i++) {
		if (anomaly_deg[i].d >= thre) anomaly++;
		else valid++;
	}
	cout << "valid node: " << valid << endl;
	out << anomaly << endl;
	for (int i = 0; i < nodeCount; i++) {
		if (anomaly_deg[i].d >= thre) out << anomaly_deg[i].i << endl;
	}
	out.close();
}

void AnomalyDetection::GraphRewrite(const string source, const string target, const double thre) {
	ifstream source_ifs(source);
	int n, m, s, t;
	source_ifs >> n >> m;
	int * mapping = new int[n];
	fill(mapping, mapping + n, 0);
	for (int i = 0; i < n; i++) {
		if (anomaly_deg[i].d >= thre) mapping[anomaly_deg[i].i] = -1;
	}
	int valid = 0;
	for (int i = 0; i < n; i++) {
		if (mapping[i] == 0) {
			mapping[i] = valid;
			valid++;
		}
	}
	int validEdge = 0;
	for (int i = 0; i < m; i++) {
		source_ifs >> s >> t;
		if (mapping[s] != -1 && mapping[t] != -1)
			validEdge++;
	}
	source_ifs.close();
	source_ifs.open(source);
	source_ifs >> n >> m;
	ofstream target_ifs(target);
	target_ifs << valid << " " << validEdge << endl;
	for (int i = 0; i < m; i++) {
		source_ifs >> s >> t;
		if (mapping[s] != -1 && mapping[t] != -1)
			target_ifs << mapping[s] << " " << mapping[t] << endl;
	}
	source_ifs.close();
	target_ifs.close();

	delete[] mapping; mapping = NULL;
}

double AnomalyDetection::GetF1(const string filename, const double thre) {
	fstream anomaly(filename + "-anomaly.txt");
	if (!anomaly) {
		cout << "can not open anomaly file: anomaly.txt" << endl;
		return -1;
	}
	int anomaly_cnt = 0;
	set<int> anml;
	string line;
	while (anomaly) {
		getline(anomaly, line);
		if (line.length() == 0)
			break;
		anomaly_cnt++;
		anml.insert(stoi(line));
	}
	anomaly.close();

	int found = 0, detect = 0;
	double f1 = 0; double best_sep = -1;
	for (int i = 0; i < nodeCount; i++) {
		if (anomaly_deg[i].d >= thre) {
			detect++;
			if (anml.find(anomaly_deg[i].i) != anml.end()) found++;
		}
	}
	double prec = found * 1.0 / detect;
	double recall = found * 1.0 / anomaly_cnt;
	//if (prec + recall > 0 && 2 * prec * recall / (prec + recall) > f1) {
	//	f1 = 2 * prec * recall / (prec + recall);
	//	best_sep = anomaly_deg[i].d;
	//}
	//cout << "best sep line: " << best_sep << endl;
	f1 = 2 * prec * recall / (prec + recall);
	return f1;
	//cout << "detect count: " << detect << endl;
	//return 2 * prec * recall / (prec + recall);
}