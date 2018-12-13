# EmbeddingAnomalyDetection

# Codes for paper "An Embedding Approach to Anomaly Detection.", ICDE, 2016

## Usage: command line parameters
filename: filename of network (no '.txt' suffix), e.g., input 'network' for 'network.txt';

d: number of dimensions, n/500 by default;

thre: parameter thre, used in AScore for detecting anomalies;

eps: parameter eps, stop condition of gradient descent, 0.001 by default;

pa: binary number, 1 print anomaly, 0 not;

rg: binary number, 1 rewrite graph by deleting anomalies, 0 not;

## network file
first line: n m (#nodes & #edges)

following m lines: s t (end points of an edge) (The indices of nodes start with 0. Only one edge of each node pair needs to be included in the edge list.)

## ground-truth of anomalies
If networks have ground-truth of anomalies, the filename of the ground-truth should by [filename]-anomaly.txt

E.g., the network filename is 'network.txt', the ground-truth should be 'network-anomaly.txt'

## pa file (print anomaly)
first line: #anomalies

following k lines: node id of an anomaly

## rg file (rewrite graph)
The format is the same to network file.

The node indices are reordered, i.e., indices of anomalies are used by other nodes.

E.g., original network has 3 edges: <0,1> <0,2> <1,2>

if 0 is detected as an anomaly, the rg network should only have 1 edge: <0,1>, where the remaining nodes are reordered.

## external library
We use the METIS library for graph partitioning.

The deployment of METIS for MS Visual Studio in x64 platform under Release mode is as follows: 

1. Open project Property Page;

2. Configuration Properties -> VC++ Directories, add the directory containing 'metis.h' & 'metis.lib' into "Include Directories" and "Library Directories";

3. Configuration Properties -> Linker -> Input, add metis.lib into "Additional Dependencies" 

if MSVS reports the unresolved external symbol error, please refer to the following page for help.
https://stackoverflow.com/questions/30412951/unresolved-external-symbol-imp-fprintf-and-imp-iob-func-sdl2/36504365#36504365

For Linux OS users, please follow the guides in the homepage of METIS.
http://glaros.dtc.umn.edu/gkhome/metis/metis/overview

## baselines 
1. ABC: Adaptive Betweenness Centrality 

Reference: Yuichi Yoshada. Almost Linear-Time Algorithms for Adaptive Betweenness Centrality using Hupergraph Sketches. In KDD, 2014.

2. OddBall

Reference: Leman Akoglu, Mary McGlohon, Christos Faloutsos. oddball: Spotting Anomalies in Weighted Graphs. In PAKDD, 2010.

3. MDS

Reference: V. de Silva and J. B. Tenenbaum. Global versus local methods in nonlinear dimensionality reduction. In NIPS, 2002.
