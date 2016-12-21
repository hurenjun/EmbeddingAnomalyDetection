#include "Graph.h"
#include "Edgeset.h"
#include "Vector.h"
#include "Embedding.h"
#include "AnomalyDetection.h"
#include <time.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <string>

using namespace std;

void printAnomaly(string & target, vector<int> & anomalies);
int main(int argc, char *argv[]) {
	if (argc != 7) {
		cout << "Usage: filename d thre eps pa rg" << endl;
		cout << "\tfilename: filename of network, (no '.txt' suffix);" << endl;
		cout << "\td: number of dimensions, n/500 by default;" << endl;
		cout << "\tthre: parameter thre, used in AScore for detecting anomalies;" << endl;
		cout << "\teps: parameter eps, stop condition of gradient descent, 0.001 by default;" << endl;
		cout << "\tpa: binary number, 1 print anomaly, 0 not;" << endl;
		cout << "\trg: binary number, 1 rewrite graph by deleting anomalies, 0 not;" << endl;
		exit(1);
	}

	cout << "filename: " << argv[1] << endl;
	string filename(argv[1]);
	int method = 2;
	cout << "d = " << argv[2] << endl;
	string dStr(argv[2]);
	int d = atoi(argv[2]);
	cout << "thre = " << argv[3] << endl;
	string thre_str(argv[3]);
	double thre = atof(argv[3]);
	cout << "eps = " << argv[4] << endl;
	string epsStr(argv[4]);
	double eps = atof(argv[4]);
	cout << "pa = " << argv[5] << endl;
	string paStr(argv[5]);
	double pa = atof(argv[5]);
	cout << "rg = " << argv[6] << endl;
	string rgStr(argv[6]);
	double rg = atof(argv[6]);
	
	srand((unsigned)time(NULL));
	time_t start = clock();
	try {
		Embedding embd(method, 0.04, 0.1, eps, filename + ".txt", d);
		embd.gradientDescent();
		AnomalyDetection ad(&embd, 0.10);
		//ad.GetPreRecCurve(filename);
		cout << "F1: " << ad.GetF1(filename, thre) << endl;
		cout << "Total time: " << (clock() - start) / 1000 << endl;

		if (pa == 1) {	// rewrite graph by removing edges with large stress
			string target = filename + "-" + dStr + '-' + thre_str + "-pa.txt";
			ad.RemoveNode(filename + ".txt", target, thre);
		}
		if (rg == 1) { 
			string target = filename + "-" + dStr + '-' + thre_str + "-rg.txt";
			ad.GraphRewrite(filename + ".txt", target, thre);
		}
		cout << endl;
	}
	catch (exception e) {
		cerr << e.what() << endl;
		exit(1);
	}
}
