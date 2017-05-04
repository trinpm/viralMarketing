/*
 * graph.h
 *
 *  Created on: May 20, 2016
 *      Author: trinpm
 */

#ifndef GRAPH_H_
#define GRAPH_H_

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <set>
#include <map>
#include <queue>
#include <iterator>
#include <algorithm>
#include <math.h>
#include "IGraph.h"
//#include "MatlabEngine.h"
#include <cstring>
#include "GraphLabel.h"
#include "GraphLabelSlow.h"
#include "PythonWrapper.h"

using namespace std;

typedef int d_type;

#define	FOR(i, n) for(int i = 0; i < (int)(n); ++i)

struct Pair {
	//used for logistic regression impl.
	int u_id;
	int v_id;
	int label; //estimated label
	int real_label; // real label for pair<u_id,v_id>, used for computing error rate
};

struct Feature {
	int rnode_size, ynode_size;
	int rnode_outdeg;
	int total_RR_edges;
	float local_indeg = 0.0, local_outdeg = 0.0, local_bc = 0.0, local_cc = 0.0,
			local_eig = 0.0, local_pr = 0.0, local_clc = 0.0;
	float glob_indeg = 0.0, glob_outdeg = 0.0, glob_bc = 0.0, glob_cc = 0.0,
			glob_eig = 0.0, glob_pr = 0.0, glob_clc = 0.0;
	int total_links2_gnodes = 0, total_links2_rnodes = 0, total_links2_ynodes =
			0;
	int m_prime = 0;
	int k = 0;
	int m = 0;
	int u_2_ynodes = 0;
	float avgDeg2Hops = 0.0;
};

class Graph {

private:

	int n; //#nodes
	int m; //#edges

	bool directedFlag;

	vector<int> nodes; //nodeid starting from 0 to n - 1
	vector<int> edges;
	vector<int> colors; //1: green, 2: yellow, 3: red

	vector<vector<int>> adjList;

	vector<int> head;
	vector<int> tail;

	vector<vector<float>> fMatrix; //feature matrix
	vector<vector<float>> fMatrix_pairwise; //feature matrix for pairwise method

	//global node centrality
	vector<float> res_glob_indeg, res_glob_outdeg, res_glob_bc, res_glob_cc,
			res_glob_eig, res_glob_pr, res_glob_clc;

public:

	Graph() {
	}
	;

	static string A0, A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11, A12, A13,
			A14, A15, A16, A17, A18, A19, A20, A21, A22, A23, A24, A25, A26,
			A27, A28;

	static string A1_v, A2_v, A3_v, A4_v, A5_v, A6_v, A7_v, A8_v, A9_v, A10_v,
			A11_v, A12_v, A13_v, A14_v, A15_v, A16_v, A17_v, A18_v, A19_v,
			A20_v, A21_v, A22_v, A23_v;

	int getNumberOfEdges() {
		return m;
	}
	;

	int getNumberOfNodes() {
		return n;
	}
	;

	vector<int> getEdges() {
		return edges;
	}
	;

	void clearPairwiseMatrix() {
		fMatrix_pairwise.clear();
	}
	;
	void clearfMatrix() {
		fMatrix.clear();
	}
	;

	//init Hung's code
	map<d_type, int> getBenefitMapping(string gfile, vector<d_type> rNodes,
			int k, bool slowFlag);

	int getBenefitOfNode(int uid);

	/**
	 * Load graph in this format:
	 * "n m" - n is #nodes, m is #edges
	 * "u1 v1"
	 * "u2 v2"
	 */
	void loadGraphFromFile(string filename, bool dFlag);

	void readVectorFromFile(string fname, vector<int> &v, bool flag);

	void makeAdjList(vector<int> head, vector<int> tail);

	int addEdge(int u, int v);

	void addNode();

	void doSamplingByBFS(float ratio, vector<int> &YNodes, vector<int> &RNodes,
			vector<int> &edgeList);

	void doSamplingByDFS(float ratio, vector<int> &YNodes, vector<int> &RNodes,
			vector<int> &edgeList);

	void doSamplingByRW(float ratio, vector<int> &YNodes, vector<int> &RNodes,
			vector<int> &edgeList);

	void BFS(int u, int &sampled, const int budget, vector<int> &edgeList);

	void DFS(int u, int &sampled, const int budget, bool &breakFlag,
			vector<int> &edgeList);

	void RW(int u, const vector<vector<float>> w, int &sampled,
			const int budget, vector<int> &edgeList);

	void generateFeatureMatrix(string gfile, vector<int> YNodes,
			vector<int> RNodes, vector<int> sampledEdgeList, int dFlag,
			bool pairwiseFlag, int k_path, int c, float alpha, bool powLawFlag,
			bool slowFlag, bool BCFlag, bool deepWalkFlag, bool CCFlag);

	void getFeatureVectorForRedNodes(vector<float> &rFeature, int u,
			const vector<int> &RNodes, const vector<int> &YNodes,
			int rnode_outdeg, const vector<float> &res_indeg,
			const vector<float> &res_outdeg, const vector<float> &res_bc,
			const vector<float> &res_cc, const vector<float> &res_eig,
			const vector<float> &res_pr, const vector<float> &res_clc, int m,
			const int total_rnode_deg, const int total_RR_edges,
			const int total_ynode_deg, IGraph &g_prime, bool BCFlag,
			bool dFlag, bool CCFlag);

	void getFeatureVector(vector<vector<float>> &rFeatures,
			vector<float> &yFeature, int u, const vector<int> &RNodes,
			const vector<int> &YNodes, int rnode_outdeg,
			const vector<float> &res_indeg, const vector<float> &res_outdeg,
			const vector<float> &res_bc, const vector<float> &res_cc,
			const vector<float> &res_eig, const vector<float> &res_pr,
			const vector<float> &res_clc, int m, const int total_rnode_deg,
			const int total_RR_edges, const int total_ynode_deg,
			IGraph &g_prime, bool BCFlag, bool dFlag, bool deepWalkFlag, bool CCFlag);

	void getCommonFeatures(vector<int> &sampledEdges, vector<int> &RNodes,
			vector<int> &YNodes, vector<float> &res_indeg,
			vector<float> &res_outdeg, int &rnode_outdeg, int &total_RR_edges,
			int &total_rnode_deg, int &total_ynode_deg, int &m, bool dFlag);

	void predictHairNode(vector<int> edges, vector<int> sampledEdges,
			vector<int> RNodes, vector<int> YNodes, string lib_path,
			string mdl_path, bool dFlag, bool BCFlag);

	void getGlobalFeature(vector<int> edges, const int dFlag, bool BCFlag, bool CCFlag);

	void getGlobalDegree(vector<int> edges, const int dFlag);

	void getLocalFeature(IGraph &igraph,
			vector<float> &res_indeg, vector<float> &res_outdeg,
			vector<float> &res_bc, vector<float> &res_cc,
			vector<float> &res_eig, vector<float> &res_pr,
			vector<float> &res_clc, const int dFlag, bool BCFlag, bool CCFlag);

//	void getLocalFeature(IGraph igraph, vector<float> &res_indeg,
//			vector<float> &res_outdeg, vector<float> &res_bc,
//			vector<float> &res_cc, vector<float> &res_eig,
//			vector<float> &res_pr, vector<float> &res_clc, const int dFlag);

	vector<vector<float>> getFeatureMatrix();

	vector<vector<float>> getPairWiseFeatureMatrix();

	int getTotalDegreeOfNodes(vector<int> nodes, vector<float> degV);

	int getTotalRedRedEdges(const vector<int> RNodes,
			const vector<int> sampledEdges);

	int getTotalColoredEdges(const int u, const vector<int> edges,
			const vector<int> RNodes, const vector<int> YNodes,
			const bool dFlag);

	int getTotalColoredEdges(const int u, const vector<int> edges,
			const vector<int> nodes, const bool dFlag);

	vector<int> convertAdjList2EdgeList(vector<vector<int>> aList);

	void generateWeights(vector<vector<float>> &w);

	void getInDegree(vector<int> &inDeg, vector<vector<int>> aList);

	void getOutDegree(vector<int> &outDeg, vector<vector<int>> aList);

	float getRN01();

	vector<int> getNeighbors(int u);

	set<int> getNeighbors(int u, IGraph &g, bool dFlag);

	void resetColors();

	vector<int> getNodes();

	bool find(vector<int> v, int e);

	bool contains(queue<int> Q, int e);

	/**
	 * @brief remove duplicated edges in case of undirected graph
	 * duplicated edge:
	 * ...
	 * 3 7 (e1)
	 * ...
	 * 7 3 (e2)
	 * ...
	 */
	void removeDuplicatedEdges(vector<int> &edges);

	bool checkEdgeExist(vector<int> edges, int u, int v, bool dFlag);

	void convertSet2Vector();

	vector<int> concat2Vectors(vector<int> A, vector<int> B);

	void removeElemInVector(int elem, vector<int> &v);

	void shuffleVector(vector<int> &v);

	void printEdgeList(vector<int> eList);

	void printQueue(queue<int> Q);

	template<class T>
	void printVector(vector<T> v, bool flag);

	void printAdjList();

	template<class T>
	void printMatrix(vector<vector<T>> matrix);

	void printSet(set<int> s);

	void printMap(map<int, double> map);

//	template<typename K, typename V>
//	void printMapping(const std::map<K, V>& mapping);

	/**
	 * @brief: write feature matrix to text file
	 */
	void write2File(vector<vector<float>> matrix, string fname);

	void write2File(string fname, vector<int> v);

	void write2File(string fname, vector<float> v);

	void write2File(string fname, vector<int> v, vector<int> rnodes,
			vector<int> ynodes, bool headerFlag);

	void doAdaptiveProbing(vector<int> edges, vector<int> sampledEdges,
			vector<int> RNodes, vector<int> YNodes, int probe_budget,
			map<string, double> coef_map, int ival, bool dFlag,
			bool logisticFlag, int k_path, int c, bool randomForestRegressFlag,
			bool randomForestClassFlag, bool rf_feature, string rfModelPath,
			bool BCFlag, bool CCFlag);

	void doAdaptiveProbing2(vector<int> edges, vector<int> sampledEdges,
			vector<int> RNodes, vector<int> YNodes, int budget,
			map<string, double> coef_map, int ival, bool dFlag, bool BCFlag, bool CCFlag, bool inverse);

	void doGreedyAdaptiveProbing(vector<int> edges, vector<int> sampledEdges,
			vector<int> RNodes, vector<int> YNodes, int budget, int ival,
			bool dFlag);

	void doGreedyNonAdaptiveProbing(vector<int> edges, vector<int> sampledEdges,
			vector<int> RNodes, vector<int> YNodes, int budget, int ival,
			bool dFlag);

	void doRandomAdaptiveProbing(vector<int> edges, vector<int> sampledEdges,
			vector<int> RNodes, vector<int> YNodes, int budget, int ival,
			bool dFlag);

	void doNonAdaptiveProbing(vector<int> edges, vector<int> sampledEdges,
			vector<int> RNodes, vector<int> YNodes, int budget,
			map<string, double> coef_map, int ival, bool dFlag,
			bool logisticFlag, int k_path, int c, bool randomForestRegressFlag,
			bool randomForestClassifyFlag, bool rf_feature, string rfModelPath,
			bool BCFlag, bool CCFlag);

	void doNonAdaptiveProbing2(vector<int> edges, vector<int> sampledEdges,
			vector<int> RNodes, vector<int> YNodes, int budget, int ival,
			bool dFlag, string lib_path, string mdl_path, bool BCFlag,
			bool deepWalkFlag, bool CCFlag);

	void doNonAdaptiveProbing3(vector<int> edges, vector<int> sampledEdges,
			vector<int> RNodes, vector<int> YNodes, int budget,
			map<string, double> coef_map, int ival, bool dFlag, bool BCFlag, bool CCFlag);

	void doNonAdaptiveProbingDeepWalk(vector<int> edges,
			vector<int> sampledEdges, vector<int> RNodes, vector<int> YNodes,
			int budget, int ival, map<int, double> node2gain, bool dFlag,
			bool BCFlag);

	void doRandomNonAdaptiveProbing(vector<int> edges, vector<int> sampledEdges,
			vector<int> RNodes, vector<int> YNodes, int budget, int ival,
			bool dFlag);

	set<int> getNeighbors_kpath(int s, int node_count, int k_path,
			vector<int> RNodes, IGraph &g, vector<int> &distance, bool dFlag);

	void getNodeInfo(int u, std::map<string, double> &info, int rnode_size,
			int ynode_size, int rnode_outdeg, int total_rnode_deg,
			int total_ynode_deg, int total_RR_edges, int m,
			vector<float> &res_indeg, vector<float> &res_outdeg,
			vector<float> &res_bc, vector<float> &res_cc, vector<float> &res_eig,
			vector<float> &res_pr, vector<float> &res_clc,
			vector<float> &res_glob_indeg, vector<float> &res_glob_outdeg,
			vector<float> &res_glob_bc, vector<float> &res_glob_cc,
			vector<float> &res_glob_eig, vector<float> &res_glob_pr,
			vector<float> &res_glob_clc, vector<int> &edges, vector<int> &RNodes,
			vector<int> &YNodes, vector<int> &sampledEdges, bool dFlag);

	void getNodeInfoForLogRegression(std::map<string, double> &info,
			vector<double> f);

	vector<double> getNodeInfoForRFClassification(
			std::map<string, double> &info, vector<double> f, bool rf_feature);

	void probeNode(int u, vector<int> &gp_nodes, vector<int> &sampledEdges,
			vector<int> &RNodes, vector<int> &YNodes, IGraph &g, bool dFlag);

	int getTotalLinks2GreenNodes(int u, vector<int> RNodes, vector<int> YNodes);

	int getTotalLinks2Nodes(int u, vector<int> nodes, bool dflag);

//	void probeNode2(int u, vector<int> &gp_nodes, vector<int> &sampledEdges,
//				vector<int> &RNodes, vector<int> &YNodes, vector<int> &GNodes, IGraph &g, bool dFlag);

	double estimateOutcome(map<string, double> coefs_map,
			map<string, double> info, bool logisticFlag);

	float totalDegree2Hops(set<int> nodes, vector<int> YNodes,
			vector<int> RNodes, vector<float> degree);

	float getCosineSimilarity(vector<float> a, vector<float> b);

	float getEuclideanSimilarity(vector<float> a, vector<float> b);

	float getSumSimiDeg(vector<vector<float>> &rFeatures,
			vector<float> &yFeature, string simiType);

};

#endif /* GRAPH_H_ */
