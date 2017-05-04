/*
 * graph.cpp
 *
 *  Created on: May 20, 2016
 *      Author: trinpm
 */

#include "graph.h"
#include <math.h>

string Graph::A0 = "INTERCEPT";

//A1~A4
string Graph::A1 = "rnode_size"; 		//A1
string Graph::A2 = "ynode_size"; 	 	//A2
string Graph::A3 = "rnode_outdeg"; 		//A3
string Graph::A4 = "total_RR_edges"; 	//A4

//A5~6
string Graph::A5 = "m"; 	//A5 = estimated #edges in G (use local information)
string Graph::A6 = "k"; 				//A6 = estimated degree of u in G

//A7~9
string Graph::A7 = "total_links2_gnodes"; //A7 - #links to green nodes, not used as a feature of node, degree dist. assumed
string Graph::A8 = "m_prime"; //A8 = estimated #edges in G (with degree distribution)
string Graph::A9 = "total_links2_ynodes"; //A9 - not used as a feature of node, degree dist. assumed

//A10~16
string Graph::A10 = "u_2_ynodes"; //A10 - Estimated #links from u to yellow nodes
string Graph::A11 = "outdeg";  			//A11
string Graph::A12 = "bc"; 				//A12
string Graph::A13 = "cc"; 				//A13
string Graph::A14 = "eig"; 				//A14
string Graph::A15 = "pr"; 				//A15
string Graph::A16 = "clc"; 				//A16

//A17~23
string Graph::A17 = "glob_indeg";		//A17
string Graph::A18 = "glob_outdeg";		//A18
string Graph::A19 = "glob_bc";			//A19
string Graph::A20 = "glob_cc";			//A20
string Graph::A21 = "glob_eig";			//A21
string Graph::A22 = "glob_pr"; 			//A22
string Graph::A23 = "glob_clc"; 		//A23

string Graph::A24 = "deg2hops";
string Graph::A25 = "cosin_simi";
string Graph::A26 = "euclid_simi";

//TODO to be removed
//string Graph::A24 = "xf1"; 		//A23

// FOR LOGISTIC REGRESSION
//A1~A4
string Graph::A1_v = "rnode_size_v"; 		//A1
string Graph::A2_v = "ynode_size_v"; 	 	//A2
string Graph::A3_v = "rnode_outdeg_v"; 		//A3
string Graph::A4_v = "total_RR_edges_v"; 	//A4

//A5~6
string Graph::A5_v = "m_v"; //A5 = estimated #edges in G (use local information)
string Graph::A6_v = "k_v"; 				//A6

//A7~9
string Graph::A7_v = "total_links2_gnodes_v"; //A7 - not used as feature of node
string Graph::A8_v = "m_prime_v"; //A8 = estimated #edges in G (with degree distribution)
string Graph::A9_v = "total_links2_ynodes_v"; //A9 - not used as feature of node

//A10~16
string Graph::A10_v = "u_2_ynodes"; //A10 - Estimated #links from u to yellow nodes
string Graph::A11_v = "outdeg_v";  			//A11
string Graph::A12_v = "bc_v"; 				//A12
string Graph::A13_v = "cc_v"; 				//A13
string Graph::A14_v = "eig_v"; 				//A14
string Graph::A15_v = "pr_v"; 				//A15
string Graph::A16_v = "clc_v"; 				//A16

//A17~23
string Graph::A17_v = "glob_indeg_v";		//A17
string Graph::A18_v = "glob_outdeg_v";		//A18
string Graph::A19_v = "glob_bc_v";			//A19
string Graph::A20_v = "glob_cc_v";			//A20
string Graph::A21_v = "glob_eig_v";			//A21
string Graph::A22_v = "glob_pr_v"; 			//A22
string Graph::A23_v = "glob_clc_v"; 		//A23

//string Graph::A24 = "label"; 		//A24 - label based on u only
//string Graph::A25 = "label_1_path";	//A25 - label based on neighbors in 1-path
//string Graph::A26 = "label_2_path";	//A26 - label based on neighbors in 2-path
//string Graph::A27 = "label_3_path";	//A27 - label based on neighbors in 3-path
//string Graph::A28 = "label_4_path";	//A28 - label based on neighbors in 4-path

map<d_type, int> Graph::getBenefitMapping(string gfile, vector<d_type> rNodes,
		int k, bool slowFlag) {

	//convert string to char*
	string str = gfile;
	char *cstr = new char[str.length() + 1];
	strcpy(cstr, str.c_str());

	map<d_type, int> benefit_mapping;
	if (slowFlag == false) {
		GraphLabel graphLabel(cstr, rNodes, k);
		cout << "estimating ben. in k-hop with fast version..." << endl;
		benefit_mapping = graphLabel.benefitMapping();
	} else {
		GraphLabelSlow graphLabelSlow(cstr, rNodes, k);
		cout << "estimating ben. in k-hop with slow version..." << endl;
		benefit_mapping = graphLabelSlow.benefitMapping();
	}

//	for (map<d_type, int>::iterator it = benefit_mapping.begin();
//			it != benefit_mapping.end(); ++it)
//		cout << it->first << " => " << it->second << '\n';

	return benefit_mapping;
}

int Graph::getBenefitOfNode(d_type uid) {
//	return benefit_mapping[uid];
}

void Graph::loadGraphFromFile(string fname, bool directed) {

	directedFlag = directed;

	std::ifstream file(fname);

	file >> n >> m;

	cout << "Reading... " << fname << "\t" << n << " nodes, " << m << " edges, "
			<< endl;

	//vector<long> edge_hash; //used for detecting multiple edges in undirected graph format, i.e., edge <1,2> and edge<2,1> should be considered as 1 edge in undirected graph

//vector<int> head(m), tail(m);
	FOR (i, m)
	{
		int u, v;
		file >> u >> v;
		head.push_back(u);
		tail.push_back(v);
		edges.push_back(u);
		edges.push_back(v);
//		file >> head[i] >> tail[i];
	}
	file.close();

// init. adj list
	FOR (i, n)
	{
		vector<int> vec;
		adjList.push_back(vec);

		colors.push_back(1);
		nodes.push_back(i);
	}

	makeAdjList(head, tail);

}

void Graph::readVectorFromFile(string fname, vector<int> &v, bool flag) {

	std::ifstream file(fname);

	int n;

	if (flag) {
		// used to read edge list
		int m;
		file >> m >> n;
		n = n * 2; //since vector of edges is 1D only
	} else {
		// used to read rfile, yfile
		file >> n;
	}

	FOR (i, n)
	{
		int u;
		file >> u;
		v.push_back(u);
	}

	file.close();

}

// COMPLEXITY: f(m) = O(m): 0 < f(m) < c*g(m)
void Graph::makeAdjList(vector<int> head, vector<int> tail) {

	cout << "makeAdjList, directedFlag: " << directedFlag << endl;

	for (unsigned int i = 0; i < head.size(); i++) {
		int u = head[i];
		int v = tail[i];

//		cout << "u-v: " << u << "-" << v << endl;

		if (directedFlag
				|| (directedFlag == false && find(adjList[u], v) == false))
			adjList[u].push_back(v);

		if (directedFlag == false && find(adjList[v], u) == false)
			adjList[v].push_back(u);
	}

}

bool Graph::find(vector<int> v, int e) {
	if (std::find(v.begin(), v.end(), e) != v.end())
		return true; // v contains e
	else
		return false; // v doesn't contain e
}

void Graph::doSamplingByBFS(float ratio, vector<int> &YNodes,
		vector<int> &RNodes, vector<int> &edgeList) {

	int budget = n * ratio;
//	cout << "#Nodes = " << n << endl;
	cout << "#Budgets = " << budget << endl;
	int sampled = 0;

	for (int i = 0; i < budget;) {
		int idx, seed;
		do {

			idx = rand() % n;
			seed = nodes[idx];
//			cout << "idx = " << idx << endl;
//			cout << "seed = " << seed << endl;

		} while (find(RNodes, seed) == true);

		// do BFS for seedId
		BFS(seed, sampled, budget, edgeList);

		FOR (j, colors.size())
		{
			if (colors[j] == 2 && find(YNodes, j) == false) {
				YNodes.push_back(j);
			} else if (colors[j] == 3 && find(RNodes, j) == false) {
				RNodes.push_back(j);
			}
		}

		i = sampled;
	}

	//	cout << "YNodes: " << endl;
//	printVector(YNodes, false);
//	cout << "RNodes: " << endl;
//	printVector(RNodes, false);

//	sampledNodes = concat2Vectors(YNodes, RNodes);
//	cout << "sampledNodes: " << endl;
//	printVector(sampledNodes, false);

//	cout << "output edgelist: " << endl;
	if (directedFlag == false)
		removeDuplicatedEdges(edgeList);
//	printEdgeList(edgeList);

//	cout << "#Nodes = " << n << endl;
//	cout << "#Budgets = " << budget << endl;
//	cout << "#Sampled = " << sampled << endl;
//	cout << "#Sampled edges = " << (edgeList.size() / 2) << endl;

	resetColors();
}

// COMPLEXITY: O(n + m)
// @brief: for BFS sampling: #budget ==  #redNodes + #yellowNodes;
//		  					 #sampledNodes ==  (#redNodes + #yellowNodes)
//							 #sampledNodes equals (or just abit larger than) #budget
void Graph::BFS(int s, int &sampled, const int budget, vector<int> &edgeList) {

//	cout << "seed = " << s << endl;

	queue<int> Q;
	Q.push(s);

	while (!Q.empty() && sampled < budget) {

		int u = Q.front();
		Q.pop();

		if (colors[u] == 1)
			sampled++;
		colors[u] = 3;

//		cout << "*u = " << u << endl;
		FOR (i, adjList[u].size())
		{
			int v = adjList[u][i];

			if (colors[v] == 1) {

//				cout << "v = " << v << endl;
//				cout << "u-v : " << u << "-" << v << endl;
				edgeList.push_back(u);
				edgeList.push_back(v);

				Q.push(v);
				colors[v] = 2;

				sampled++;
			} else {
				// even if node v is traversed, the edge u-v still needs to be stored
//				cout << "v = " << v << endl;
//				cout << "u-v : " << u << "-" << v << endl;
				edgeList.push_back(u);
				edgeList.push_back(v);
			}

		}
//		cout << "sampled = " << sampled << endl;
	}
}

void Graph::doSamplingByDFS(float ratio, vector<int> &YNodes,
		vector<int> &RNodes, vector<int> &edgeList) {
//vector<int> Graph::doSamplingByDFS(float ratio) {

	int budget = n * ratio;
//	cout << "#Nodes = " << n << endl;
	cout << "#Budgets = " << budget << endl;
	int sampled = 0;
	bool flag = false;

	for (int i = 0; i < budget;) {
		int idx, seed;
		do {

			idx = rand() % n;
			seed = nodes[idx];
//			cout << "idx = " << idx << endl;
			cout << "seed = " << seed << endl;
//			cout << endl;

		} while (find(RNodes, seed) == true);

		// do DFS for seedId
		DFS(seed, sampled, budget, flag, edgeList);

		FOR (j, colors.size())
		{
			if (colors[j] == 2 && find(YNodes, j) == false) {
				YNodes.push_back(j);
			} else if (colors[j] == 3 && find(RNodes, j) == false) {
				RNodes.push_back(j);
			}
		}

		i = sampled;
	}

//	cout << "#Nodes = " << n << endl;
//	cout << "#Budgets = " << budget << endl;
//	cout << "#Sampled = " << sampled << endl;
//
//	cout << "YNodes: " << endl;
//	printVector(YNodes, false);
//	cout << "RNodes: " << endl;
//	printVector(RNodes, false);

//	cout << "output edgelist: " << endl;
	if (directedFlag == false)
		removeDuplicatedEdges(edgeList);
//	printEdgeList(edgeList);

	resetColors();
}

// O(n + m)
// @brief: for DFS sampling: #budget can equal #redNodes; or
//                           #sampledNodes == (#red nodes + #yellow nodes)
//							 #sampledNodes can be much more larger than #budget
void Graph::DFS(int u, int &sampled, const int budget, bool &flag,
		vector<int> &edgeList) {

	if (sampled > budget) {
		return;
	}

	if (colors[u] == 1)
		sampled++;

	colors[u] = 2;

	FOR(i, adjList[u].size())
	{
		int v = adjList[u][i];

		if (colors[v] == 1) {
			colors[v] = 2;
			sampled++;
			edgeList.push_back(u);
			edgeList.push_back(v);

			DFS(v, sampled, budget, flag, edgeList);
		} else {
			edgeList.push_back(u);
			edgeList.push_back(v);
		}

	}
	colors[u] = 3;

	if (sampled > budget)
		return;

}

void Graph::doSamplingByRW(float ratio, vector<int> &YNodes,
		vector<int> &RNodes, vector<int> &edgeList) {

	vector<vector<float>> w;
	vector<vector<int>> aList;

//init. nb, w, r
	FOR(i, n)
	{
		vector<int> t;
		aList.push_back(t);

		vector<float> t2;
		w.push_back(t2);
	}

//generate weight 0~1 for all edges in graph
	generateWeights(w);

//	vector<int> outDeg(n);
	//get out-degree of node in G - it is also degree of node (for undirected node)
//	getOutDegree(outDeg, adjList);

	int budget = n * ratio;
//	cout << "#Nodes = " << n << endl;
	cout << "#Budgets = " << budget << endl;
	int sampled = 0;

	vector<int> outDegSample(n);

	for (int i = 0; i < budget;) {

		int idx, seed;
		do {

			idx = rand() % n;
			seed = nodes[idx];
//			cout << "idx = " << idx << endl;
//			cout << "seed = " << seed << endl;

		} while (find(RNodes, seed) == true);

		RW(seed, w, sampled, budget, edgeList);

		FOR(j, colors.size())
		{
			//cout << "j=" << j << endl;
			if (colors[j] == 3 && find(RNodes, j) == false) {
				RNodes.push_back(j);
			} else if (colors[j] == 2 && find(YNodes, j) == false) {
				YNodes.push_back(j);
			}
		}

		i = sampled;
	}

//	cout << "#Nodes = " << n << endl;
//	cout << "#Budgets = " << budget << endl;
//	cout << "#Sampled = " << sampled << endl;
//
//	cout << "YNodes: " << endl;
//	printVector(YNodes, false);
//	cout << "RNodes: " << endl;
//	printVector(RNodes, false);

//	printAdjList(aList);
//	cout << "output edgelist: " << endl;
	if (directedFlag == false) {
		removeDuplicatedEdges(edgeList);
	}
//	printEdgeList(edgeList);

	resetColors();

}
// edgeList
// 1 2 0.5
// 1 3 0.6
// u(1): 2 3
// w(1): 0.5 0.6

// @brief: for DFS sampling: #budget == #redNodes;
//                           #sampledNodes == (#red nodes + #yellow nodes)
//							 #sampledNodes can be much more larger than #budget
void Graph::RW(int s, const vector<vector<float>> w, int &sampled,
		const int budget, vector<int> &edgeList) {

//	cout << "start from seed = " << s << endl;
	vector<int> RNode_candidates;
	queue<int> Q;
	Q.push(s);

	// start randomly collecting red nodes
	int estimated_budget = 0;
	while (!Q.empty() && sampled <= budget) {

		int u = Q.front();
		Q.pop();

		if (colors[u] != 3) {

			if (find(RNode_candidates, u) == false) {
				RNode_candidates.push_back(u);
				sampled++;

				FOR(i, adjList[u].size())
				{
					int v = adjList[u][i];
					float rn = getRN01();

					if (rn < w[u][i]) {
						if (contains(Q, v) == false && colors[v] != 3) {
							estimated_budget = adjList[v].size();
							Q.push(v);
						}
					}
				}
			}
		}
	}

//	printVector(RNode_candidates, true);

	//start probing rnodes
	FOR(i, RNode_candidates.size())
	{
		int u = RNode_candidates[i];

		// get all neighbors of u
		FOR(j, adjList[u].size())
		{
			int v = adjList[u][j];
			if (colors[v] == 1)
				colors[v] = 2; // turn v to yellow node
			edgeList.push_back(u);
			edgeList.push_back(v);
		}

		colors[u] = 3; // turn u to red nodes
	}

}

void Graph::getGlobalFeature(vector<int> edges, const int dFlag, bool BCFlag,
		bool CCFlag) {

	IGraph igraph;

	igraph.loadGraphFromEdgeList(edges, dFlag);

	//compute global node centrality
	res_glob_indeg = igraph.computeNodeCentrality(property::IN_DEGREE, dFlag);
//	cout << "1. IN_DEGREE (GLOB.)" << endl;
//	printVector(res_glob_indeg, true);
	res_glob_outdeg = igraph.computeNodeCentrality(property::OUT_DEGREE, dFlag);
//	cout << "2. OUT_DEGREE (GLOB.)" << endl;
//	printVector(res_glob_outdeg, true);
	if (BCFlag) {
		res_glob_bc = igraph.computeNodeCentrality(property::BETWEENNESS_CENT,
				dFlag);
//		cout << "3. BETWEENNESS_CENT (GLOB.) DISABLED" << endl;
//			printVector(res_glob_bc, true);
	}

	if (CCFlag) {
		res_glob_cc = igraph.computeNodeCentrality(property::CLOSENESS_CENT,
				dFlag);
	}

//	cout << "4. CLOSENESS_CENT (GLOB.)" << endl;
//	printVector(res_glob_cc, true);
	res_glob_eig = igraph.computeNodeCentrality(property::EIGEN_CENT, dFlag);
//	cout << "5. EIGEN_CENT (GLOB.)" << endl;
//	printVector(res_glob_eig, true);
	res_glob_pr = igraph.computeNodeCentrality(property::PAGERANK_CENT, dFlag);
//	cout << "6. PAGERANK_CENT (GLOB.)" << endl;
//	printVector(res_glob_pr, true);

	//compute global cluster coeff.
	res_glob_clc = igraph.computeClusCoef(dFlag);
//	cout << "7. CLUS_COEFF (GLOB.)" << endl;
//	printVector(res_glob_clc, true);
}

void Graph::getGlobalDegree(vector<int> edges, const int dFlag) {

	IGraph igraph;

	igraph.loadGraphFromEdgeList(edges, dFlag);

	//compute global node centrality
	res_glob_indeg = igraph.computeNodeCentrality(property::IN_DEGREE, dFlag);
//	cout << "1. IN_DEGREE (GLOB.) (f17):" << endl;
//	printVector(res_glob_indeg, true);
	res_glob_outdeg = igraph.computeNodeCentrality(property::OUT_DEGREE, dFlag);
	//	cout << "2. OUT_DEGREE (GLOB.) (f18):" << endl;
	//	printVector(res_glob_outdeg, true);
}

//void Graph::getLocalFeature(IGraph &igraph, vector<float> &res_indeg,
//		vector<float> &res_outdeg, vector<float> &res_bc, vector<float> &res_cc,
//		vector<float> &res_eig, vector<float> &res_pr, vector<float> &res_clc,
//		const int dFlag, bool BCFlag) {
//
//	//IGraph igraph;
//	//igraph.loadGraphFromEdgeList(sampledEdges, dFlag);
//
//	//compute local node centrality
//	res_indeg = igraph.computeNodeCentrality(property::IN_DEGREE, dFlag);
////	cout << "res_indeg" << endl;
//	res_outdeg = igraph.computeNodeCentrality(property::OUT_DEGREE, dFlag);
////	cout << "res_outdeg" << endl;
//	if (BCFlag) {
//		//enable or disable calculating BC for performance issue
//		res_bc = igraph.computeNodeCentrality(property::BETWEENNESS_CENT,
//				dFlag);
//	} else {
//		vector<float> res_bc_0(1); //init value 0 for any element in vector
//		res_bc = res_bc_0;
//	}
//
////	cout << "res_bc" << endl;
//	res_cc = igraph.computeNodeCentrality(property::CLOSENESS_CENT, dFlag);
////	cout << "res_cc" << endl;
//	res_eig = igraph.computeNodeCentrality(property::EIGEN_CENT, dFlag);
////	cout << "res_eig" << endl;
//	res_pr = igraph.computeNodeCentrality(property::PAGERANK_CENT, dFlag);
////	cout << "res_pr" << endl;
//
//	//compute local cluster coeff.
//	res_clc = igraph.computeClusCoef(dFlag);
//
//}

void Graph::getLocalFeature(IGraph &igraph, vector<float> &res_indeg,
		vector<float> &res_outdeg, vector<float> &res_bc, vector<float> &res_cc,
		vector<float> &res_eig, vector<float> &res_pr, vector<float> &res_clc,
		const int dFlag, bool BCFlag, bool CCFlag) {

//compute local node centrality
	res_indeg = igraph.computeNodeCentrality(property::IN_DEGREE, dFlag);
	res_outdeg = igraph.computeNodeCentrality(property::OUT_DEGREE, dFlag);

	cout << "** 1-1" << endl;
	if (BCFlag)
		res_bc = igraph.computeNodeCentrality(property::BETWEENNESS_CENT,
				dFlag);
	else {
		vector<float> res_bc_0(1); //init value 0 for any element in vector
		res_bc = res_bc_0;
	}

//	cout << "CCFlag: " << CCFlag << endl;
	if (CCFlag)
		res_cc = igraph.computeNodeCentrality(property::CLOSENESS_CENT, dFlag);
	else {
		vector<float> res_cc_0(1);
		res_cc = res_cc_0;
	}

	cout << "** 1-2" << endl;
	res_eig = igraph.computeNodeCentrality(property::EIGEN_CENT, dFlag);
	cout << "** 1-3" << endl;
	res_pr = igraph.computeNodeCentrality(property::PAGERANK_CENT, dFlag);
	cout << "** 1-4" << endl;

	//compute local cluster coeff.
	res_clc = igraph.computeClusCoef(dFlag);
	cout << "** 1-5" << endl;
}

void Graph::generateFeatureMatrix(string gfile, vector<int> YNodes,
		vector<int> RNodes, vector<int> sampledEdges, int dFlag, bool pwFlag,
		int k_path, int c, float alpha, bool powLawFlag, bool slowFlag,
		bool BCFlag, bool deepWalkFlag, bool CCFlag) {

	IGraph g;		// load original graph
	g.loadGraphFromEdgeList(edges, dFlag);

	IGraph g_prime; // sample graph
	g_prime.loadGraphFromEdgeList(sampledEdges, dFlag);

//	printEdgeList(sampledEdges);

	vector<float> res_indeg, res_outdeg, res_bc, res_cc, res_eig, res_pr,
			res_clc;

	//1. get node centrality of all nodes in sampled graph
	getLocalFeature(g_prime, res_indeg, res_outdeg, res_bc, res_cc, res_eig,
			res_pr, res_clc, dFlag, BCFlag, CCFlag);

	int rnode_outdeg = 0;
	int total_RR_edges = 0;
	int total_rnode_deg = 0, total_ynode_deg = 0;
	int m = 0;
	getCommonFeatures(sampledEdges, RNodes, YNodes, res_indeg, res_outdeg,
			rnode_outdeg, total_RR_edges, total_rnode_deg, total_ynode_deg, m,
			dFlag);

//  LABEL BY HUNG'S PROGRAM
//	vector<map<d_type, int>> vmapping;
//	if (powLawFlag) {
//		vector<int> k_hops;
//		k_hops.push_back(1);
//		k_hops.push_back(2);
//		k_hops.push_back(3);
//
//		for (int kid = 0; kid < k_hops.size(); kid++) {
//			int k = k_hops[kid];
//			cout << "k = " << k << endl;
//			map<d_type, int> ben_map = getBenefitMapping(gfile, RNodes, k, slowFlag);
//			vmapping.push_back(ben_map);
//		}
//	}

	//feature vector for red nodes
	vector<vector<float>> redFeatures;
	//feature vector for red nodes
	FOR(i, RNodes.size())
	{
		vector<float> rFeature;
		Feature fstr;

		int u = RNodes[i];
		getFeatureVectorForRedNodes(rFeature, u, RNodes, YNodes, rnode_outdeg,
				res_indeg, res_outdeg, res_bc, res_cc, res_eig, res_pr, res_clc,
				m, total_rnode_deg, total_RR_edges, total_ynode_deg, g_prime,
				BCFlag, dFlag, CCFlag);

		//name of red node corresponds to this feature
		rFeature.push_back(u);

		redFeatures.push_back(rFeature);
	}

//	cout << "print red feature matrix: " << endl;
//	printMatrix(redFeatures);

	//feature vector for yellow nodes
	FOR(i, YNodes.size())
	{
		vector<float> yFeature;
		int u = YNodes[i];
//		cout << "node u: " << u << endl;

		getFeatureVector(redFeatures, yFeature, u, RNodes, YNodes, rnode_outdeg,
				res_indeg, res_outdeg, res_bc, res_cc, res_eig, res_pr, res_clc,
				m, total_rnode_deg, total_RR_edges, total_ynode_deg, g_prime,
				BCFlag, dFlag, deepWalkFlag, CCFlag);

		// LABEL BY HUNG'S PROGRAM
		//		if (powLawFlag) {
		//			FOR(m, vmapping.size())
		//			{
		//				map<d_type, int> ben_map = vmapping[m];
		//				int u_ben = ben_map[u];
		//				feature.push_back(u_ben);
		//			}
		//		}

		//add to feature matrix
		fMatrix.push_back(yFeature);
	}

//	cout << "print feature matrix: " << endl;
//	printMatrix(fMatrix);

	if (pwFlag) {
		//in case feature matrix is used for pairwise ranking
		for (int i = 0; i < fMatrix.size() - 1; i++) {
			vector<float> f1 = fMatrix[i];
			for (int j = i + 1; j < fMatrix.size(); j++) {
				vector<float> f2 = fMatrix[j];
				vector<float> f;
				f.insert(f.end(), f1.begin(), f1.end());
				f.insert(f.end(), f2.begin(), f2.end());

				fMatrix_pairwise.push_back(f);
			}
		}

		cout << "size of feature matrix (pairwise):" << fMatrix_pairwise.size()
				<< endl;
//		cout << "print out feature matrix (pairwise):" << endl;
//		printMatrix (fMatrix_pairwise);
	}

//	cout << endl;
//	cout << "print out feature matrix:" << endl;
//	printMatrix (fMatrix);

}

void Graph::getCommonFeatures(vector<int> &sampledEdges, vector<int> &RNodes,
		vector<int> &YNodes, vector<float> &res_indeg,
		vector<float> &res_outdeg, int &rnode_outdeg, int &total_RR_edges,
		int &total_rnode_deg, int &total_ynode_deg, int &m, bool dFlag) {
	int rnode_indeg = 0;
	int ynode_indeg = 0;
	int ynode_outdeg = 0;
	rnode_indeg = getTotalDegreeOfNodes(RNodes, res_indeg);
	rnode_outdeg = getTotalDegreeOfNodes(RNodes, res_outdeg);
	ynode_indeg = getTotalDegreeOfNodes(YNodes, res_indeg);
	ynode_outdeg = getTotalDegreeOfNodes(YNodes, res_outdeg);
	total_RR_edges = getTotalRedRedEdges(RNodes, sampledEdges);

	if (dFlag) {
		total_rnode_deg = rnode_indeg + rnode_outdeg;
		total_ynode_deg = ynode_indeg + ynode_outdeg;
	} else {
		total_rnode_deg = rnode_outdeg;
		total_ynode_deg = ynode_outdeg;
	}
	m = pow(total_rnode_deg, 2) / (2 * total_RR_edges); // TODO: calculate m without degree distribution info
}

float Graph::getSumSimiDeg(vector<vector<float>> &rFeatures,
		vector<float> &yFeature, string simiType) {
//	cout << "---" << endl;
	float sumSimiDeg = 0.0;
	//traverse all features of red nodes:
	for (int j = 0; j < rFeatures.size(); j++) {
		vector<float> rFeature = rFeatures[j];

		vector<float> rVector, yVector;

//		for (int idx = 0; idx < rFeatures[0].size()-1; idx++) { // get only local node centrality to measure similarity
		for (int idx = 10; idx < 16; idx++) { // get only local node centrality to measure similarity
			rVector.push_back(rFeature[idx]);
			yVector.push_back(yFeature[idx]);
		}

		//	printVector(rVector, true);
		//	printVector(yVector, true);
		float simi = 0.0;
		if (simiType == "cos") {
			simi = getCosineSimilarity(rVector, yVector);
//			cout << "Cosine similarity: " << simi << endl;
		} else if (simiType == "euclid") {
			simi = getEuclideanSimilarity(rVector, yVector);
//			cout << "Euclidean similarity: " << simi << endl;
		}
//		sumSimiDeg += simi;
		int rNode = rFeature[rFeature.size() - 1]; //name of current Red Node
		sumSimiDeg += simi * res_glob_indeg[rNode];
	}
	return sumSimiDeg;
}

float Graph::getCosineSimilarity(vector<float> a, vector<float> b) {

	float cos_similarity;

	float nominator = 0.0;
	float a_sum = 0.0;
	float b_sum = 0.0;
	for (int i = 0; i < a.size(); i++) {
		nominator += a[i] * b[i];

		a_sum += a[i] * a[i];
		b_sum += b[i] * b[i];
	}

	cos_similarity = nominator / (sqrt(a_sum) * sqrt(b_sum));

	return cos_similarity;
}

float Graph::getEuclideanSimilarity(vector<float> a, vector<float> b) {

	float diff_sum = 0.0;

	for (int i = 0; i < a.size(); i++) {
		float diff = a[i] - b[i];
		diff_sum += diff * diff;
	}

	return sqrt(diff_sum);
}

void Graph::getFeatureVector(vector<vector<float>> &rFeatures,
		vector<float> &yFeature, int u, const vector<int> &RNodes,
		const vector<int> &YNodes, int rnode_outdeg,
		const vector<float> &res_indeg, const vector<float> &res_outdeg,
		const vector<float> &res_bc, const vector<float> &res_cc,
		const vector<float> &res_eig, const vector<float> &res_pr,
		const vector<float> &res_clc, int m, const int total_rnode_deg,
		const int total_RR_edges, const int total_ynode_deg, IGraph &g_prime,
		bool BCFlag, bool dFlag, bool deepWalkFlag, bool CCFlag) {

	Feature f;

	f.rnode_size = RNodes.size();
	f.ynode_size = YNodes.size();
	f.rnode_outdeg = rnode_outdeg;
	f.total_RR_edges = total_RR_edges;

	f.local_indeg = res_indeg[u];
	f.local_outdeg = res_outdeg[u];
	if (BCFlag)
		f.local_bc = res_bc[u];
	else
		f.local_bc = 0.0;

	if (CCFlag)
		f.local_cc = res_cc[u];
	else
		f.local_cc = 0.0;

	f.local_eig = res_eig[u];
	f.local_pr = res_pr[u];
	f.local_clc = res_clc[u];

	//global node centrality
	f.glob_indeg = res_glob_indeg[u];
	f.glob_outdeg = res_glob_outdeg[u];
	if (BCFlag)
		f.glob_bc = res_glob_bc[u]; //DISABLE
	else
		f.glob_bc = 0.0;

	if (CCFlag)
		f.glob_cc = res_glob_cc[u]; //DISABLE
	else
		f.glob_cc = 0.0;

	f.glob_eig = res_glob_eig[u];
	f.glob_pr = res_glob_pr[u];
	f.glob_clc = res_glob_clc[u];

	f.total_links2_gnodes = getTotalLinks2GreenNodes(u, RNodes, YNodes); //option 1: #links from u to green nodes (degree dist. assumed)
	f.total_links2_rnodes = getTotalLinks2Nodes(u, RNodes, dFlag); //option 2: #links from u to red nodes
	f.total_links2_ynodes = getTotalLinks2Nodes(u, YNodes, dFlag); //option 3: #links from u to yellow nodes (degree dist. assumed)

	f.m_prime = getNumberOfEdges(); //TODO: calculate m with degree distribution available

	//estimate m
	if (total_RR_edges != 0)
		f.k = (f.total_links2_rnodes * total_rnode_deg) / total_RR_edges;

	if (total_rnode_deg != 0)
		f.u_2_ynodes = (total_ynode_deg * f.total_links2_rnodes)
				/ total_rnode_deg; // estimate # links from u to yellow nodes

	set<int> nb_2hops;
	nb_2hops = getNeighbors(u, g_prime, dFlag);
	//cout << "2 hops nb of u: " << endl;
	//printSet(nb_2hops);
	float totalDeg2Hops = totalDegree2Hops(nb_2hops, YNodes, RNodes, res_indeg); // total degree of neighbors of u in 2 hops
	f.avgDeg2Hops = (float) totalDeg2Hops / (float) nb_2hops.size();

	//add to feature vector (23 features)

	if (deepWalkFlag) {
		// todo: just to test deep walk
		yFeature.push_back(u);
		yFeature.push_back(f.total_links2_gnodes);
	} else {

		//A1~A4
		yFeature.push_back(f.rnode_size);	//A1
		yFeature.push_back(f.ynode_size);  	//A2
		yFeature.push_back(f.rnode_outdeg); 	//A3
		yFeature.push_back(f.total_RR_edges);	//A4

		//A5~6
		yFeature.push_back(m); 				//A5
		yFeature.push_back(f.k); 				//A6

		//A7~9
		yFeature.push_back(f.total_links2_gnodes);  //A7
		yFeature.push_back(f.m_prime);  			  //A8
		yFeature.push_back(f.total_links2_ynodes);  //A9

		//A10~16
		yFeature.push_back(f.u_2_ynodes); 	//A10 - u_2_ynodes
		yFeature.push_back(f.local_outdeg);  	//A11
		yFeature.push_back(f.local_bc); 		//A12
		yFeature.push_back(f.local_cc); 		//A13
		yFeature.push_back(f.local_eig); 		//A14
		yFeature.push_back(f.local_pr); 		//A15
		yFeature.push_back(f.local_clc); 		//A16

		//A17~23
		yFeature.push_back(f.glob_indeg);  //A17
		yFeature.push_back(f.glob_outdeg); //A18
		yFeature.push_back(f.glob_bc); 	//A19
		yFeature.push_back(f.glob_cc); 	//A20
		yFeature.push_back(f.glob_eig); 	//A21
		yFeature.push_back(f.glob_pr); 	//A22
		yFeature.push_back(f.glob_clc); 	//A23

		yFeature.push_back(f.avgDeg2Hops); //A24

		float sumCosSimiDeg = getSumSimiDeg(rFeatures, yFeature, "cos");
		yFeature.push_back(sumCosSimiDeg); //A25

		float sumEuclidSimiDeg = getSumSimiDeg(rFeatures, yFeature, "euclid");
		yFeature.push_back(sumEuclidSimiDeg); //A26
	}

}

void Graph::getFeatureVectorForRedNodes(vector<float> &rFeature, int u,
		const vector<int> &RNodes, const vector<int> &YNodes, int rnode_outdeg,
		const vector<float> &res_indeg, const vector<float> &res_outdeg,
		const vector<float> &res_bc, const vector<float> &res_cc,
		const vector<float> &res_eig, const vector<float> &res_pr,
		const vector<float> &res_clc, const int m, const int total_rnode_deg,
		const int total_RR_edges, const int total_ynode_deg, IGraph &g_prime,
		bool BCFlag, bool dFlag, bool CCFlag) {

	Feature f;

	f.rnode_size = RNodes.size();
	f.ynode_size = YNodes.size();
	f.rnode_outdeg = rnode_outdeg;
	f.total_RR_edges = total_RR_edges;

	f.local_indeg = res_indeg[u];
	f.local_outdeg = res_outdeg[u];
	if (BCFlag)
		f.local_bc = res_bc[u];
	else
		f.local_bc = 0.0;
	f.local_cc = res_cc[u];
	f.local_eig = res_eig[u];
	f.local_pr = res_pr[u];
	f.local_clc = res_clc[u];

	//global node centrality
	f.glob_indeg = res_glob_indeg[u];
	f.glob_outdeg = res_glob_outdeg[u];
	if (BCFlag)
		f.glob_bc = res_glob_bc[u]; //DISABLE
	else
		f.glob_bc = 0.0;

	if (CCFlag)
		f.glob_cc = res_glob_cc[u]; //DISABLE
	else
		f.glob_cc = 0.0;

	f.glob_eig = res_glob_eig[u];
	f.glob_pr = res_glob_pr[u];
	f.glob_clc = res_glob_clc[u];

	f.total_links2_gnodes = getTotalLinks2GreenNodes(u, RNodes, YNodes); //option 1: #links from u to green nodes (degree dist. assumed)
	f.total_links2_rnodes = getTotalLinks2Nodes(u, RNodes, dFlag); //option 2: #links from u to red nodes
	f.total_links2_ynodes = getTotalLinks2Nodes(u, YNodes, dFlag); //option 3: #links from u to yellow nodes (degree dist. assumed)

	f.m_prime = getNumberOfEdges(); //TODO: calculate m with degree distribution available

	//estimate m

	if (total_RR_edges != 0)
		f.k = (f.total_links2_rnodes * total_rnode_deg) / total_RR_edges;

	if (total_rnode_deg != 0)
		f.u_2_ynodes = (total_ynode_deg * f.total_links2_rnodes)
				/ total_rnode_deg; // estimate # links from u to yellow nodes

	set<int> nb_2hops;
	nb_2hops = getNeighbors(u, g_prime, dFlag);
	//cout << "2 hops nb of u: " << endl;
	//printSet(nb_2hops);
	float totalDeg2Hops = totalDegree2Hops(nb_2hops, YNodes, RNodes, res_indeg); // total degree of neighbors of u in 2 hops
	f.avgDeg2Hops = (float) totalDeg2Hops / (float) nb_2hops.size();

	//add to feature vector (23 features)

	//A1~A4
	rFeature.push_back(f.rnode_size); 		//A1
	rFeature.push_back(f.ynode_size);  		//A2
	rFeature.push_back(f.rnode_outdeg); 	//A3
	rFeature.push_back(f.total_RR_edges);	//A4

	//A5~6
	rFeature.push_back(m); 					//A5
	rFeature.push_back(f.k); 				//A6

	//A7~9
	rFeature.push_back(f.total_links2_gnodes);   //A7
	rFeature.push_back(f.m_prime);  			 //A8
	rFeature.push_back(f.total_links2_ynodes);   //A9

	//A10~16
	rFeature.push_back(f.u_2_ynodes); 		//A10 - u_2_ynodes
	rFeature.push_back(f.local_outdeg);  	//A11
	rFeature.push_back(f.local_bc); 		//A12
	rFeature.push_back(f.local_cc); 		//A13
	rFeature.push_back(f.local_eig); 		//A14
	rFeature.push_back(f.local_pr); 		//A15
	rFeature.push_back(f.local_clc); 		//A16

	//A17~23
	rFeature.push_back(f.glob_indeg);	//A17
	rFeature.push_back(f.glob_outdeg);	//A18
	rFeature.push_back(f.glob_bc); 		//A19
	rFeature.push_back(f.glob_cc);		//A20
	rFeature.push_back(f.glob_eig); 	//A21
	rFeature.push_back(f.glob_pr);		//A22
	rFeature.push_back(f.glob_clc); 	//A23

	rFeature.push_back(f.avgDeg2Hops);  //A24
}

double Graph::estimateOutcome(map<string, double> coefs_map,
		map<string, double> info, bool logisticFlag) {

	if (coefs_map.size() == info.size() + 1) {
		double result = 0;
		result = coefs_map[A0];

		std::map<string, double>::iterator it;
		for (it = coefs_map.begin(); it != coefs_map.end(); ++it) {
			string key = it->first;

			if (key != A0) {
				result = result + coefs_map[key] * info[key];

//				cout << "key: " << key << ", coef: " << coefs_map[key]
//						<< ", info: " << info[key] << endl;
//				cout << "result: " << result << endl;
			}
		}

		if (logisticFlag) {
			// result for logistic regression
			double tmp;
			tmp = 1 / (1 + exp(-result));
//			cout << "logistic(result) = " << tmp << endl;

			if (tmp >= 0.5)
				return 1.0;
			else
				return 0.0;
		}

		return result; // result for linear regression
	}
	cout << "size of map not matched!" << endl;
	return 0;
}

int Graph::getTotalLinks2GreenNodes(int u, vector<int> RNodes,
		vector<int> YNodes) {
	// get all neighbors of u
	int gnode_total = 0;
	FOR(j, adjList[u].size())
	{
		int v = adjList[u][j];
		if (find(RNodes, v) == true)
			continue;
		if (find(YNodes, v) == true)
			continue;
		gnode_total++;
	}
//	cout << "u: " << u << ", gnode_total: " << gnode_total << endl;
	return gnode_total;
}

//COMPLEXITY: sizeOf(neighborsOfu) * sizeOf(nodes)
int Graph::getTotalLinks2Nodes(int u, vector<int> nodes, bool dFlag) {
	// get all neighbors of u
	int r_links = 0;
	FOR(j, adjList[u].size())
	{
		int v = adjList[u][j];
		if (find(nodes, v) == true)
			r_links++;
	}
	return r_links;
}

vector<int> Graph::getNeighbors(int u) {

	vector<int> neighbors_of_u;
	FOR(j, adjList[u].size())
	{
		int v = adjList[u][j];
		neighbors_of_u.push_back(v);
	}

	return neighbors_of_u;
}

void Graph::doGreedyAdaptiveProbing(vector<int> edges, vector<int> sampledEdges,
		vector<int> RNodes, vector<int> YNodes, int budget, int ival,
		bool dFlag) {

	cout << "doGreedyAdaptiveProbing" << endl;

	// get list of nodes of gprime
	vector<int> gp_nodes;		// gprime's nodes
	gp_nodes = concat2Vectors(RNodes, YNodes);
	cout << "gp_nodes's size: " << gp_nodes.size() << endl;

	cout << "#nodes in GPrime: " << RNodes.size() + YNodes.size() << endl;
	cout << "budget: " << budget << endl;

	cout << "GREEDY PROBING" << endl;
	cout << "RNode's size: " << RNodes.size() << endl;
	cout << "YNode's size: " << YNodes.size() << endl;

	IGraph g;		// load original graph
	g.loadGraphFromEdgeList(edges, dFlag);

	for (int k = 1; k <= budget; k++) {

		vector<double> du_out;
		map<int, double> u_map;

		FOR(i, YNodes.size())
		{
			int u = YNodes[i];

			int total_green_edges = 0;

			//check if neighbor of u is in green nodes:

			total_green_edges = getTotalLinks2GreenNodes(u, RNodes, YNodes);
//			total_green_edges = getTotalColoredEdges(u, edges, RNodes, YNodes, dFlag); //option 1: #links from u to green nodes (degree dist. assumed)

			int outcome = total_green_edges;
			du_out.push_back(outcome);
			u_map[u] = outcome;
		}

		// sort YNodes based on node's outcome

		cout << "du_out: " << endl;
		printVector(du_out, true);
		// sort du_out
		std::sort(du_out.begin(), du_out.end(), std::greater<double>());
		double du_out_max = du_out[0];
		cout << "du_out_max: " << du_out_max << endl;

		for (int i = 0; i < YNodes.size(); i++) {

			int u = YNodes[i];
			double du_out = u_map[u];

			// PROBE node u - which has highest outcome
			if (find(RNodes, u) == false && du_out == du_out_max) {
				cout << "start probing u: " << u << endl;
				probeNode(u, gp_nodes, sampledEdges, RNodes, YNodes, g, dFlag);

				break;
			}
		}

		if (k == 1) {
			cout << "k = " << k << ", #Nodes in augmented Gprime: "
					<< gp_nodes.size() << endl;
		} else if (k % ival == 0) {

			cout << "k = " << k << ", #Nodes in augmented Gprime: "
					<< gp_nodes.size() << endl;

		}

		shuffleVector(YNodes);

	}
}

void Graph::doGreedyNonAdaptiveProbing(vector<int> edges,
		vector<int> sampledEdges, vector<int> RNodes, vector<int> YNodes,
		int budget, int ival, bool dFlag) {

	cout << "doGreedyNonAdaptiveProbing" << endl;

	// get list of nodes of gprime
	vector<int> gp_nodes;		// gprime's nodes
	gp_nodes = concat2Vectors(RNodes, YNodes);
	cout << "gp_nodes's size: " << gp_nodes.size() << endl;

	cout << "#nodes in GPrime: " << RNodes.size() + YNodes.size() << endl;
	cout << "budget: " << budget << endl;

	cout << "GREEDY NON-ADAPTIVE PROBING" << endl;
//	cout << "GNode's size: " << GNodes.size() << endl;
	cout << "RNode's size: " << RNodes.size() << endl;
	cout << "YNode's size: " << YNodes.size() << endl;

	IGraph g;		// load original graph
	g.loadGraphFromEdgeList(edges, dFlag);

	vector<double> du_out;
	map<int, double> u_map;

	FOR(i, YNodes.size())
	{
		int u = YNodes[i];
		int total_green_edges = 0;

		//cout << "finding green links... for " << u << endl;

		total_green_edges = getTotalLinks2GreenNodes(u, RNodes, YNodes);
		//	total_green_edges = getTotalColoredEdges(u, edges, RNodes, YNodes, dFlag); //option 1: #links from u to green nodes (degree dist. assumed)

		int outcome = total_green_edges;
		du_out.push_back(outcome);
		u_map[u] = outcome;
	}

	// sort du_out
	std::sort(du_out.begin(), du_out.end(), std::greater<double>());
	vector<double> du_out_sorted(du_out.begin(), du_out.begin() + budget);

	int k = 0;
	vector<int> probed_nodes;
	vector<int> YNodes_org = YNodes;
	if (YNodes_org.size() > budget) {

		while (k < budget) {

			for (int i = 0; i < YNodes_org.size(); i++) {

				double du_out = du_out_sorted[k];

				int u = YNodes_org[i];
				double outcome_u = u_map[u];

				if (du_out == outcome_u && find(RNodes, u) == false) {
					// probe node u
					cout << "start probing u: " << u << endl;

					probed_nodes.push_back(u);

					probeNode(u, gp_nodes, sampledEdges, RNodes, YNodes, g,
							dFlag);

					k = k + 1;

					if (k == 1) {
						cout << "k = " << k << ", #Nodes in augmented Gprime: "
								<< gp_nodes.size() << endl;
					} else if (k % ival == 0) {
						cout << "k = " << k << ", #Nodes in augmented Gprime: "
								<< gp_nodes.size() << endl;
					}
				}
			}
		}
	} else
		cerr << "YNodes.size() < budget" << endl;
}

void Graph::doRandomAdaptiveProbing(vector<int> edges, vector<int> sampledEdges,
		vector<int> RNodes, vector<int> YNodes, int budget, int ival,
		bool dFlag) {

	cout << "doRandomAdaptiveProbing" << endl;

	// get list of nodes of gprime
	vector<int> gp_nodes;		// gprime's nodes
	gp_nodes = concat2Vectors(RNodes, YNodes);

	cout << "#nodes in GPrime: " << RNodes.size() + YNodes.size() << endl;
	cout << "budget: " << budget << endl;

	IGraph g;		// load original graph
	g.loadGraphFromEdgeList(edges, dFlag);

	for (int k = 1; k <= budget; k++) {

		vector<double> du_out;
		map<int, double> u_map;

		int idx = rand() % YNodes.size();
//		cout << "random idx: " << idx << endl;
		int u = YNodes[idx]; //randomly select a node in Ynodes to probe

		// PROBE node u
		if (find(RNodes, u) == false) {
//			cout << "start probing u: " << u << endl;
			probeNode(u, gp_nodes, sampledEdges, RNodes, YNodes, g, dFlag);
		}

		if (k == 1) {
			cout << "k = " << k << ", #Nodes in augmented Gprime: "
					<< gp_nodes.size() << endl;
		} else if (k % ival == 0) {

			cout << "k = " << k << ", #Nodes in augmented Gprime: "
					<< gp_nodes.size() << endl;
		}

		shuffleVector(YNodes);
	}
}

void Graph::doRandomNonAdaptiveProbing(vector<int> edges,
		vector<int> sampledEdges, vector<int> RNodes, vector<int> YNodes,
		int budget, int ival, bool dFlag) {

	cout << "doRandomNonAdaptiveProbing" << endl;

	// get list of nodes of gprime
	vector<int> gp_nodes;		// gprime's nodes
	gp_nodes = concat2Vectors(RNodes, YNodes);

	cout << "RANDOM PROBING" << endl;
	cout << "#nodes in GPrime: " << RNodes.size() + YNodes.size() << endl;
	cout << "budget: " << budget << endl;

	IGraph g;		// load original graph
	g.loadGraphFromEdgeList(edges, dFlag);

	//randomly get k yellow nodes
	vector<int> nodes;
	for (int k = 1; k <= budget;) {

		int idx = rand() % YNodes.size();
		int u = YNodes[idx]; //randomly select a node in Ynodes to probe

		if (find(nodes, u) == false) {
			nodes.push_back(u);
			k++;
		}
	}

	for (int k = 1; k <= budget; k++) {
		int u = nodes[k - 1];

		// PROBE node u
		if (find(RNodes, u) == false) {
			cout << "start probing u: " << u << endl;
			probeNode(u, gp_nodes, sampledEdges, RNodes, YNodes, g, dFlag);
		}

		if (k == 1) {
			cout << "k = " << k << ", #Nodes in augmented Gprime: "
					<< gp_nodes.size() << endl;
		} else if (k % ival == 0) {

			cout << "k = " << k << ", #Nodes in augmented Gprime: "
					<< gp_nodes.size() << endl;
		}
	}

}

/**
 *
 */
void Graph::doAdaptiveProbing(vector<int> edges, vector<int> sampledEdges,
		vector<int> RNodes, vector<int> YNodes, int budget,
		map<string, double> coef_map, int ival, bool dFlag, bool logisticFlag,
		int k_path, int c, bool randomForestRegressFlag,
		bool randomForestClassifyFlag, bool rf_feature, string rfModelPath,
		bool BCFlag, bool CCFlag) {
	/*
	 // get list of nodes of gprime
	 vector<int> gp_nodes;		// gprime's nodes
	 gp_nodes = concat2Vectors(RNodes, YNodes);

	 cout << "#nodes in GPrime: " << RNodes.size() + YNodes.size() << endl;
	 cout << "budget: " << budget << endl;

	 IGraph g;		// load original graph
	 g.loadGraphFromEdgeList(edges, dFlag);
	 int node_count = g.get_n(); // get #nodes in G

	 //TODO: disable random forest
	 MatlabEngine matlabEngine;
	 if (randomForestRegressFlag || randomForestClassifyFlag) {
	 //init. matlab engine
	 matlabEngine.initEngine();
	 matlabEngine.loadModel(rfModelPath);
	 }

	 int rnode_size = RNodes.size();
	 int ynode_size = YNodes.size();

	 int total_RR_edges = getTotalRedRedEdges(RNodes, sampledEdges);

	 //get global festure of node u
	 getGlobalDegree(edges, dFlag);

	 vector<float> error_rate_v;
	 vector<float> f1_score_v;

	 // start probing
	 for (int k = 1; k <= budget; k++) {

	 //		cout << "** 0" << endl;
	 // cout << "k = " << k << endl;
	 // load sampledEdges to igraph object g

	 // get local feature of g
	 vector<float> res_indeg, res_outdeg, res_bc, res_cc, res_eig, res_pr,
	 res_clc;
	 IGraph gprime;
	 gprime.loadGraphFromEdgeList(sampledEdges, dFlag);
	 getLocalFeature(gprime, res_indeg, res_outdeg, res_bc, res_cc, res_eig,
	 res_pr, res_clc, dFlag, BCFlag, CCFlag);

	 //		cout << "BC score: " << endl;
	 ////		printVector(res_bc, true);
	 //		write2File("/home/student/workspace/GraphSampling/data/all/BFS/set4/check/tella5_bc_score", res_bc);
	 //
	 //		cout << "CLC score: " << endl;
	 //		write2File("/home/student/workspace/GraphSampling/data/all/BFS/set4/check/tella5_clc_score", res_clc);
	 ////		printVector(res_clc, true);
	 //
	 //		cout << "PR score: " << endl;
	 //		write2File("/home/student/workspace/GraphSampling/data/all/BFS/set4/check/tella5_pr_score", res_pr);
	 ////		printVector(res_pr, true);
	 //
	 //		cout << "CC score: " << endl;
	 //		write2File("/home/student/workspace/GraphSampling/data/all/BFS/set4/check/tella5_cc_score", res_cc);
	 ////		printVector(res_cc, true);

	 int rnode_indeg = 0;
	 int rnode_outdeg = 0;

	 int ynode_indeg = 0;
	 int ynode_outdeg = 0;

	 rnode_indeg = getTotalDegreeOfNodes(RNodes, res_indeg);
	 rnode_outdeg = getTotalDegreeOfNodes(RNodes, res_outdeg);

	 ynode_indeg = getTotalDegreeOfNodes(YNodes, res_indeg);
	 ynode_outdeg = getTotalDegreeOfNodes(YNodes, res_outdeg);

	 int total_rnode_deg;
	 int total_ynode_deg;
	 if (dFlag) {
	 total_rnode_deg = rnode_indeg + rnode_outdeg;
	 total_ynode_deg = ynode_indeg + ynode_outdeg;
	 } else {
	 total_rnode_deg = rnode_outdeg;
	 total_ynode_deg = ynode_outdeg;
	 }
	 int m;
	 m = pow(total_rnode_deg, 2) / (2 * total_RR_edges); // TODO: estimate m without degree distribution info

	 //	cout << "12. m':" << m << endl << endl;

	 vector<int> du;
	 vector<double> du_out;
	 map<int, double> u_map;
	 map<int, int> id_name_map;

	 vector<vector<double>> matrix;
	 vector<vector<double>> datapoint_rf;
	 //		cout << "YNodes.size:" << YNodes.size() << endl;
	 //		printVector(YNodes, true);
	 FOR(i, YNodes.size())
	 {
	 int u = YNodes[i];

	 //cout << "checking node u = " << u << endl;

	 std::map<string, double> info_u;

	 //			cout << "** 1" << endl;
	 //get info of node u
	 getNodeInfo(u, info_u, rnode_size, ynode_size, rnode_outdeg,
	 total_rnode_deg, total_ynode_deg, total_RR_edges, m,
	 res_indeg, res_outdeg, res_bc, res_cc, res_eig, res_pr,
	 res_clc, res_glob_indeg, res_glob_outdeg, res_glob_bc,
	 res_glob_cc, res_glob_eig, res_glob_pr, res_glob_clc, edges,
	 RNodes, YNodes, sampledEdges, dFlag);
	 //			cout << "** 2" << endl;

	 if (logisticFlag || randomForestClassifyFlag) {

	 //using logistic regression
	 //info[A24] = res_glob_outdeg[u] - res_outdeg[u]; //TODO get alpha from command line

	 id_name_map[i] = u; //mapping from nodes' id to internal index

	 vector<double> feature; // get features of node u

	 feature.push_back(info_u[A1]);
	 feature.push_back(info_u[A2]);
	 feature.push_back(info_u[A3]);
	 feature.push_back(info_u[A4]);

	 feature.push_back(info_u[A5]);	//m
	 feature.push_back(info_u[A6]);

	 feature.push_back(info_u[A7]);	//total_green_edges;
	 feature.push_back(info_u[A8]);	//m_prime;
	 feature.push_back(info_u[A9]); 	//total_links2_ynodes;
	 feature.push_back(info_u[A10]);	//u_2_ynodes;

	 feature.push_back(info_u[A11]);
	 feature.push_back(info_u[A12]);
	 feature.push_back(info_u[A13]);
	 feature.push_back(info_u[A14]);
	 feature.push_back(info_u[A15]);
	 feature.push_back(info_u[A16]);

	 feature.push_back(info_u[A17]);	//glob_indeg;
	 feature.push_back(info_u[A18]);
	 feature.push_back(info_u[A19]);
	 feature.push_back(info_u[A20]);
	 feature.push_back(info_u[A21]);
	 feature.push_back(info_u[A22]);
	 feature.push_back(info_u[A23]);

	 //print vector f:
	 //cout << "feature vector of node " << u << ":" << endl;
	 //printVector(feature, true);

	 matrix.push_back(feature);

	 } else if (randomForestRegressFlag) {

	 // using random forest
	 vector<double> datapoint;

	 if (rf_feature == 0) {
	 // RF with local features
	 datapoint.push_back(info_u[A1]);
	 datapoint.push_back(info_u[A2]);
	 datapoint.push_back(info_u[A3]);
	 datapoint.push_back(info_u[A4]);
	 datapoint.push_back(info_u[A5]);
	 datapoint.push_back(info_u[A6]);
	 datapoint.push_back(info_u[A11]);
	 datapoint.push_back(info_u[A12]);
	 datapoint.push_back(info_u[A13]);
	 datapoint.push_back(info_u[A14]);
	 datapoint.push_back(info_u[A15]);
	 datapoint.push_back(info_u[A16]);
	 } else if (rf_feature == 1) {
	 // RF with global features
	 datapoint.push_back(info_u[A1]);
	 datapoint.push_back(info_u[A2]);
	 datapoint.push_back(info_u[A3]);
	 datapoint.push_back(info_u[A4]);
	 datapoint.push_back(info_u[A10]);
	 datapoint.push_back(info_u[A11]);
	 datapoint.push_back(info_u[A12]);
	 datapoint.push_back(info_u[A13]);
	 datapoint.push_back(info_u[A14]);
	 datapoint.push_back(info_u[A15]);
	 datapoint.push_back(info_u[A16]);
	 datapoint.push_back(info_u[A17]);
	 }
	 datapoint_rf.push_back(datapoint);

	 //				printVector(datapoint, true);
	 //				double outcome = 0.0;
	 //				outcome = matlabEngine.predict(datapoint);
	 ////				cout << "outcome of " << u << ": " << outcome << endl;
	 //				du_out.push_back(outcome);
	 //				u_map[u] = outcome;

	 // estimate outcome of node u
	 } else {

	 //using linear regression
	 double total_outcome = 0.0; // outcome from kpath neighbors

	 //				//print info_map
	 //				cout << "y node: " << u << endl;
	 //				map<string, double>::iterator it;
	 //				for (it = info_u.begin(); it != info_u.end(); ++it) {
	 //					string key = it->first;
	 //					double value = it->second;
	 //
	 //					cout << "key: " << key << ", value: " << value << endl;
	 //				}

	 // estimate outcome of node u
	 //				cout << "bc score of u: " << info_u["bc"] << endl;
	 double outcome_u = estimateOutcome(coef_map, info_u,
	 logisticFlag);
	 //				cout << "outcome of u: " << outcome_u << endl;

	 total_outcome = outcome_u;

	 //cout << "outcome of " << u << ": " << outcome << endl;
	 du_out.push_back(total_outcome);
	 u_map[u] = total_outcome;

	 }
	 }

	 //start ranking and probing
	 if (logisticFlag || randomForestClassifyFlag) {

	 int max = 0;
	 for (int i = 1; i < matrix.size(); i++) {

	 //concatenate info of matrix[max], and matrix[i]
	 vector<double> f1 = matrix[max];
	 vector<double> f2 = matrix[i];

	 vector<double> f;
	 f.insert(f.end(), f1.begin(), f1.end());
	 f.insert(f.end(), f2.begin(), f2.end());

	 // get local feature of node u
	 std::map<string, double> info;
	 double outcome = 0.0;
	 if (logisticFlag) {

	 getNodeInfoForLogRegression(info, f);
	 outcome = estimateOutcome(coef_map, info, logisticFlag); //binary outcome

	 } else if (randomForestClassifyFlag) {

	 //TODO: disable random forest
	 vector<double> datapoint;
	 datapoint = getNodeInfoForRFClassification(info, f,
	 rf_feature);
	 outcome = matlabEngine.predict(datapoint); //binary outcome
	 }
	 if (outcome == 0) { // gain(matrix[max]) < gain(matrix[i])
	 max = i;
	 }
	 }

	 int idx = max;

	 //from id_name_map, get back node id of u in graph
	 int u; // node has highest rank
	 u = id_name_map[idx];
	 //cout << "internal index of u: " << idx << endl;
	 //cout << "node u that has highest rank: " << u << endl;

	 if (find(RNodes, u) == false) {

	 //				cout << "start probing u: " << u << endl;
	 probeNode(u, gp_nodes, sampledEdges, RNodes, YNodes, g, dFlag);

	 //cout << "** 8" << endl;
	 }
	 } else if (randomForestRegressFlag) {

	 //TODO: disable random forest
	 // random forest fore regression
	 int idx_u = matlabEngine.predict(datapoint_rf);
	 int u = YNodes[idx_u - 1]; //in matlab index of array starting from 1

	 // PROBE node u - which has highest predicted outcome
	 if (find(RNodes, u) == false) {

	 //cout << "start probing u: " << u << endl;
	 probeNode(u, gp_nodes, sampledEdges, RNodes, YNodes, g, dFlag);
	 }

	 } else {
	 // linear regression

	 // sort YNodes based on node's outcome

	 // sort du_out
	 std::sort(du_out.begin(), du_out.end(), std::less<double>());
	 double du_out_max = du_out[du_out.size() - 1];
	 //		cout << "du_out_max: " << du_out_max << endl;

	 map<int, double>::iterator it;
	 for (it = u_map.begin(); it != u_map.end(); it++) {

	 int u = it->first;
	 double du_out = it->second;

	 // PROBE node u - which has highest outcome
	 if (find(RNodes, u) == false && du_out == du_out_max) {

	 //					cout << "***********************************" << endl;
	 //					cout << "start probing u: " << u << endl;

	 probeNode(u, gp_nodes, sampledEdges, RNodes, YNodes, g,
	 dFlag);
	 //					cout << "** 3" << endl;

	 break;

	 }
	 }
	 }

	 if (k == 1) {
	 cout << "k = " << k << ", #Nodes in augmented Gprime: "
	 << gp_nodes.size() << endl;

	 } else if (k % ival == 0) {
	 cout << "k = " << k << ", #Nodes in augmented Gprime: "
	 << gp_nodes.size() << endl;

	 }

	 }

	 float error_rate_sum = 0.0;
	 float f1_score_sum = 0.0;
	 FOR(i,error_rate_v.size())
	 {
	 error_rate_sum += error_rate_v[i];
	 f1_score_sum += f1_score_v[i];
	 }

	 cout << "Error rate (mean) (#errors/#data) of this model: "
	 << error_rate_sum * 100 / error_rate_v.size() << "%" << endl;
	 cout << "F1 score (mean) of this model: "
	 << f1_score_sum / f1_score_v.size() << endl;

	 //	cout << "YNodes: " << endl;
	 //	printVector(YNodes, true);
	 //	cout << "RNodes: " << endl;
	 //	printVector(RNodes, true);
	 //	cout << "Nodes in Gprime: " << endl;
	 //	printVector(gp_nodes, true);

	 if (randomForestRegressFlag || randomForestClassifyFlag) {
	 //close matlab engine
	 //TODO: disable random forest
	 matlabEngine.closeEngine();
	 }*/
}

void Graph::doAdaptiveProbing2(vector<int> edges, vector<int> sampledEdges,
		vector<int> RNodes, vector<int> YNodes, int budget,
		map<string, double> coef_map, int ival, bool dFlag, bool BCFlag,
		bool CCFlag, bool inverse) {

	cout << "doAdaptiveProbing2" << endl;

	// get list of nodes of gprime
	vector<int> gp_nodes;		// gprime's nodes
	gp_nodes = concat2Vectors(RNodes, YNodes);

	cout << "#nodes in GPrime: " << RNodes.size() + YNodes.size() << endl;
	cout << "budget: " << budget << endl;

	IGraph g;		// load original graph
	g.loadGraphFromEdgeList(edges, dFlag);
	int node_count = g.get_n();

	int rnode_size = RNodes.size();
	int ynode_size = YNodes.size();

	int total_RR_edges = getTotalRedRedEdges(RNodes, sampledEdges);

	getGlobalDegree(edges, dFlag);

	cout << "coef_map" << endl;
	map<string, double>::iterator it;
	for (it = coef_map.begin(); it != coef_map.end(); ++it) {
		string key = it->first;
		double value = it->second;

		cout << "feature: " << key << ", coef: " << value << endl;
	}

	// start probing
	for (int k = 1; k <= budget; k++) {

		cout << "Ynodes.size: " << YNodes.size() << endl;
		cout << "Rnodes.size: " << RNodes.size() << endl;

		if (0) {
			cout << "sampledEdges.size: " << sampledEdges.size() / 2 << endl;
			set<int> nodes;
			for (int j = 0; j < sampledEdges.size(); j++) {
				nodes.insert(sampledEdges[j]);
			}
			cout << "set.size: " << nodes.size() << endl;
			cout << "nodes in set:" << endl;
			printSet(nodes);
		}

		cout << "k = " << k << endl;
		cout << "** 0" << endl;

		IGraph gprime;
		gprime.loadGraphFromEdgeList(sampledEdges, dFlag);

		// get local feature of g_prime
		vector<float> res_indeg, res_outdeg, res_bc, res_cc, res_eig, res_pr,
				res_clc;

		cout << "** 1" << endl;
		getLocalFeature(gprime, res_indeg, res_outdeg, res_bc, res_cc, res_eig,
				res_pr, res_clc, dFlag, BCFlag, CCFlag);

		int rnode_indeg = 0;
		int rnode_outdeg = 0;

		int ynode_indeg = 0;
		int ynode_outdeg = 0;

		rnode_indeg = getTotalDegreeOfNodes(RNodes, res_indeg);
		rnode_outdeg = getTotalDegreeOfNodes(RNodes, res_outdeg);
		cout << "** 2" << endl;
		ynode_indeg = getTotalDegreeOfNodes(YNodes, res_indeg);
		cout << "** 3" << endl;
		ynode_outdeg = getTotalDegreeOfNodes(YNodes, res_outdeg);
		cout << "** 4" << endl;

		int total_rnode_deg;
		int total_ynode_deg;
		if (dFlag) {
			total_rnode_deg = rnode_indeg + rnode_outdeg;
			total_ynode_deg = ynode_indeg + ynode_outdeg;
		} else {
			total_rnode_deg = rnode_outdeg;
			total_ynode_deg = ynode_outdeg;
		}
		int m;
		m = pow(total_rnode_deg, 2) / (2 * total_RR_edges); // TODO: estimate m without degree distribution info

		//	cout << "12. m':" << m << endl << endl;

		vector<int> du;
		vector<double> du_out;
		map<int, double> u_map;

		FOR(i, YNodes.size())
		{
			int u = YNodes[i];

			//cout << "checking node u = " << u << endl;

			std::map<string, double> info_u;

			//get info of node u
			getNodeInfo(u, info_u, rnode_size, ynode_size, rnode_outdeg,
					total_rnode_deg, total_ynode_deg, total_RR_edges, m,
					res_indeg, res_outdeg, res_bc, res_cc, res_eig, res_pr,
					res_clc, res_glob_indeg, res_glob_outdeg, res_glob_bc,
					res_glob_cc, res_glob_eig, res_glob_pr, res_glob_clc, edges,
					RNodes, YNodes, sampledEdges, dFlag);

//				//print info_map
//				cout << "y node: " << u << endl;
//				map<string, double>::iterator it;
//				for (it = info_u.begin(); it != info_u.end(); ++it) {
//					string key = it->first;
//					double value = it->second;
//
//					cout << "key: " << key << ", value: " << value << endl;
//				}

			// estimate outcome of node u
			double outcome_u = estimateOutcome(coef_map, info_u, false);

			if (0) {
				cout << "local degree of " << u << ": " << res_indeg[u] << endl;
				cout << "pagerank of " << u << ": " << res_pr[u] << endl;
				cout << "eig. value of " << u << ": " << res_eig[u] << endl;
				cout << "closeness of " << u << ": " << res_cc[u] << endl;
				cout << "clus. coefs of " << u << ": " << res_clc[u] << endl;
				cout << "glb degree of " << u << ": " << res_glob_indeg[u]
						<< endl;
				cout << "--outcome of " << u << ": " << outcome_u << endl;
			}

			du_out.push_back(outcome_u);
			u_map[u] = outcome_u;
		}

		if (1) {
			cout << "du_out: " << endl;
			printVector(du_out, false);
		}

		// sort YNodes based on node's outcome

		// sort du_out
		std::sort(du_out.begin(), du_out.end(), std::greater<double>());

		double du_out_max;
		if (inverse == false)
			du_out_max = du_out[0];
		else
			du_out_max = du_out[du_out.size() - 1]; //select node with lowest ranking

		cout << "du_out_max: " << du_out_max << endl;

		cout << "** 5" << endl;

		for (int i = 0; i < YNodes.size(); i++) {

			int u = YNodes[i];
			double du_out = u_map[u];

			// PROBE node u - which has highest outcome
			if (find(RNodes, u) == false && du_out == du_out_max) {

//					cout << "***********************************" << endl;
				cout << "start probing u: " << u << endl;

				probeNode(u, gp_nodes, sampledEdges, RNodes, YNodes, g, dFlag);

				break;
			}
		}

		cout << "** 5" << endl;

		if (k == 1) {
			cout << "k = " << k << ", #Nodes in augmented Gprime: "
					<< gp_nodes.size() << endl;

		} else if (k % ival == 0) {
			cout << "k = " << k << ", #Nodes in augmented Gprime: "
					<< gp_nodes.size() << endl;

		}
		// shuffle updated YNodes
		shuffleVector(YNodes);
	}

}

void Graph::shuffleVector(vector<int> &v) {
	cout << "shuffle vector" << endl;
	std::random_device rd;
	std::mt19937 g(rd());
	std::shuffle(v.begin(), v.end(), g);
}

/**
 *
 */
void Graph::predictHairNode(vector<int> edges, vector<int> sampledEdges,
		vector<int> RNodes, vector<int> YNodes, string lib_path,
		string mdl_path, bool dFlag, bool BCFlag) {
	/*
	 // get list of nodes of gprime
	 vector<int> gp_nodes;		// gprime's nodes
	 gp_nodes = concat2Vectors(RNodes, YNodes);

	 cout << "#nodes in GPrime: " << RNodes.size() + YNodes.size() << endl;

	 IGraph g;		//load original graph
	 g.loadGraphFromEdgeList(edges, dFlag);

	 IGraph g_prime; // sample graph
	 g_prime.loadGraphFromEdgeList(sampledEdges, dFlag);

	 //init. Python engine
	 //	lib_path = "/home/student/PycharmProjects/IsotonicRegression";
	 //	mdl_path = "/home/student/PycharmProjects/IsotonicRegression/model/ca_RF.mdl";
	 //	mdl_path = "/home/student/PycharmProjects/IsotonicRegression/model/as70_RF.mdl";
	 PythonWrapper pythonWrapper;
	 pythonWrapper.init();
	 pythonWrapper.setPath(lib_path);
	 pythonWrapper.importModule("Pickle");
	 int loaded = 0;
	 loaded = pythonWrapper.loadModel("loadModel", mdl_path);
	 if (loaded != 1) {
	 cerr << "cannot load trained model" << endl;
	 } else {
	 cout << "trained model loaded!" << endl;
	 }

	 //get global feature of node u
	 //	getGlobalDegree(edges, dFlag);
	 getGlobalFeature(edges, dFlag, false);

	 // get local feature of g
	 vector<float> res_indeg, res_outdeg, res_bc, res_cc, res_eig, res_pr,
	 res_clc;
	 getLocalFeature(sampledEdges, g_prime, res_indeg, res_outdeg, res_bc,
	 res_cc, res_eig, res_pr, res_clc, dFlag, BCFlag);

	 int rnode_outdeg = 0;
	 int total_RR_edges = 0;
	 int total_rnode_deg = 0, total_ynode_deg = 0;
	 int m = 0;
	 getCommonFeatures(sampledEdges, RNodes, YNodes, res_indeg, res_outdeg,
	 rnode_outdeg, total_RR_edges, total_rnode_deg, total_ynode_deg, m,
	 dFlag);

	 //feature vector for red nodes
	 vector<vector<float>> redFeatures;
	 FOR(i, RNodes.size())
	 {
	 vector<float> rFeature;
	 Feature fstr;

	 int u = RNodes[i];
	 getFeatureVectorForRedNodes(rFeature, u, RNodes, YNodes, rnode_outdeg,
	 res_indeg, res_outdeg, res_bc, res_cc, res_eig, res_pr, res_clc,
	 m, total_rnode_deg, total_RR_edges, total_ynode_deg, g_prime,
	 BCFlag, dFlag);

	 rFeature.push_back(u); //name of red node corresponds to this feature
	 redFeatures.push_back(rFeature);
	 }

	 vector<float> error_rate_v;
	 vector<float> f1_score_v;

	 vector<int> groundTruth;
	 vector<int> pred_res;

	 cout << "Size of YNodes: " << YNodes.size() << endl;
	 FOR(i, YNodes.size())
	 {
	 vector<float> yFeature;
	 Feature fstr;

	 int u = YNodes[i];
	 //		cout << "node u: " << u << endl;

	 getFeatureVector(redFeatures, yFeature, u, RNodes, YNodes, rnode_outdeg,
	 res_indeg, res_outdeg, res_bc, res_cc, res_eig, res_pr, res_clc,
	 m, total_rnode_deg, total_RR_edges, total_ynode_deg, g_prime,
	 BCFlag, dFlag);

	 //select features used for prediction
	 vector<double> X;
	 FOR(j, yFeature.size())
	 {
	 if ((j == 10) || (j >= 12 && j <= 15) || (j >= 23)) {
	 X.push_back(yFeature[j]);
	 }
	 }

	 //get ground-truth from global graph:
	 if (res_glob_outdeg[u] == 1)
	 groundTruth.push_back(1); //u is hairNode
	 else
	 groundTruth.push_back(0);

	 //		cout << "0" << endl;
	 //		cout << "X.size: " << X.size() << endl;
	 float pred = pythonWrapper.predict("predict", X);
	 cout << "---" << endl;
	 cout << "node u: " << u << endl;
	 cout << "gt: " << groundTruth[i] << endl;
	 cout << "pred: " << pred << endl;
	 pred_res.push_back(pred);
	 }

	 //get error rate
	 int correct_cnt = 0;
	 if (groundTruth.size() == pred_res.size()) {
	 FOR(i, groundTruth.size())
	 {
	 if (groundTruth[i] == pred_res[i]) {
	 correct_cnt++;
	 }
	 }
	 } else {
	 cerr << "Different size btw groundtruth and predicted vector!" << endl;
	 }

	 cout << "Error rate (%): " << correct_cnt * 100 / groundTruth.size()
	 << endl;

	 //	float error_rate_sum = 0.0;
	 //	float f1_score_sum = 0.0;
	 //	FOR(i,error_rate_v.size())
	 //	{
	 //		error_rate_sum += error_rate_v[i];
	 //		f1_score_sum += f1_score_v[i];
	 //	}
	 //
	 //	cout << "Error rate (mean) (#errors/#data) of this model: "
	 //			<< error_rate_sum * 100 / error_rate_v.size() << "%" << endl;
	 //	cout << "F1 score (mean) of this model: "
	 //			<< f1_score_sum / f1_score_v.size() << endl;

	 pythonWrapper.clean();
	 */
}
//vector<vector<float>> Graph::getRedFeatures(vector<int> &sampledEdges,
//		vector<int> &RNodes, vector<int> &YNodes, IGraph &g_prime,
//		vector<float> &res_indeg, vector<float> &res_outdeg,
//		vector<float> &res_bc, vector<float> &res_cc, vector<float> &res_eig,
//		vector<float> &res_pr, vector<float> &res_clc, const int dFlag,
//		bool BCFlag) {
//
//	int rnode_indeg = 0;
//	int rnode_outdeg = 0;
//	int ynode_indeg = 0;
//	int ynode_outdeg = 0;
//
//	rnode_indeg = getTotalDegreeOfNodes(RNodes, res_indeg);
//	rnode_outdeg = getTotalDegreeOfNodes(RNodes, res_outdeg);
//
//	ynode_indeg = getTotalDegreeOfNodes(YNodes, res_indeg);
//	ynode_outdeg = getTotalDegreeOfNodes(YNodes, res_outdeg);
//
//	int rnode_size = RNodes.size();
//	int ynode_size = YNodes.size();
//
//	int total_RR_edges = getTotalRedRedEdges(RNodes, sampledEdges);
//
//	int total_rnode_deg, total_ynode_deg;
//	if (dFlag) {
//		total_rnode_deg = rnode_indeg + rnode_outdeg;
//		total_ynode_deg = ynode_indeg + ynode_outdeg;
//	} else {
//		total_rnode_deg = rnode_outdeg;
//		total_ynode_deg = ynode_outdeg;
//	}
//
//	int m;
//	m = pow(total_rnode_deg, 2) / (2 * total_RR_edges); // TODO: calculate m without degree distribution info
//
//	vector<vector<float>> redFeatures;
//	//feature vector for red nodes
//	FOR(i, RNodes.size())
//	{
//		vector<float> feature;
//		Feature fstr;
//
//		int u = RNodes[i];
//		//		cout << "node u: " << u << endl;
//
//		getFeatureVector(u, fstr, RNodes, YNodes, rnode_outdeg, res_indeg,
//				res_outdeg, res_bc, res_cc, res_eig, res_pr, res_clc,
//				total_rnode_deg, total_RR_edges, total_ynode_deg, g_prime,
//				BCFlag, dFlag);
//
//
//
//		//add to feature vector (23 features)
//
//		//A1~A4
//		feature.push_back(fstr.rnode_size); 	//A1
//		feature.push_back(fstr.ynode_size);  	//A2
//		feature.push_back(fstr.rnode_outdeg); 	//A3
//		feature.push_back(fstr.total_RR_edges);	//A4
//
//		//A5~6
//		feature.push_back(fstr.m); 				//A5
//		feature.push_back(fstr.k); 				//A6
//
//		//A7~9
//		feature.push_back(fstr.total_links2_gnodes);  //A7
//		feature.push_back(fstr.m_prime);  			  //A8
//		feature.push_back(fstr.total_links2_ynodes);  //A9
//
//		//A10~16
//		feature.push_back(fstr.u_2_ynodes); 	//A10 - u_2_ynodes
//		feature.push_back(fstr.local_outdeg);  	//A11
//		feature.push_back(fstr.local_bc); 		//A12
//		feature.push_back(fstr.local_cc); 		//A13
//		feature.push_back(fstr.local_eig); 		//A14
//		feature.push_back(fstr.local_pr); 		//A15
//		feature.push_back(fstr.local_clc); 		//A16
//
//		//A17~23
//		feature.push_back(fstr.glob_indeg);  //A17
//		feature.push_back(fstr.glob_outdeg); //A18
//		feature.push_back(fstr.glob_bc); 	//A19
//		feature.push_back(fstr.glob_cc); 	//A20
//		feature.push_back(fstr.glob_eig); 	//A21
//		feature.push_back(fstr.glob_pr); 	//A22
//		feature.push_back(fstr.glob_clc); 	//A23
//
//		//A24
//		feature.push_back(fstr.avgDeg2Hops);
//
//		//name of red node corresponds to this feature
//		feature.push_back(u);
//
//		redFeatures.push_back(feature);
//	}
//	return redFeatures;
//}

void Graph::getNodeInfoForLogRegression(std::map<string, double> &info,
		vector<double> f) {
// get local feature of node u

	info[A1] = f[0];
	info[A2] = f[1];
	info[A3] = f[2];
	info[A4] = f[3];

	info[A5] = f[4];
	info[A6] = f[5];
	info[A7] = f[6];	//0, total_green_edges;
	info[A8] = f[7];	//m_prime;
	info[A9] = f[8]; 	//0, total_links2_ynodes;

	info[A10] = f[9];	//u_2_ynodes;
	info[A11] = f[10];
	info[A12] = f[11];
	info[A13] = f[12];
	info[A14] = f[13];
	info[A15] = f[14];
	info[A16] = f[15];

	info[A17] = f[16];	//0, glob_indeg;
	info[A18] = f[17];
	info[A19] = f[18];
	info[A20] = f[19];
	info[A21] = f[20];
	info[A22] = f[21];
	info[A23] = f[22];

	info[A1_v] = f[23];
	info[A2_v] = f[24];
	info[A3_v] = f[25];
	info[A4_v] = f[26];

	info[A5_v] = f[27];
	info[A6_v] = f[28];
	info[A7_v] = f[29];	//0, total_green_edges;
	info[A8_v] = f[30];	//m_prime;
	info[A9_v] = f[31];	//0, total_links2_ynodes;

	info[A10_v] = f[32]; //u_2_ynodes
	info[A11_v] = f[33];
	info[A12_v] = f[34];
	info[A13_v] = f[35];
	info[A14_v] = f[36];
	info[A15_v] = f[37];
	info[A16_v] = f[38];

	info[A17_v] = f[39]; //0, glob_indeg;
	info[A18_v] = f[40];
	info[A19_v] = f[41];
	info[A20_v] = f[42];
	info[A21_v] = f[43];
	info[A22_v] = f[44];
	info[A23_v] = f[45];
}

vector<double> Graph::getNodeInfoForRFClassification(
		std::map<string, double> &info, vector<double> f, bool rf_feature) {

	vector<double> datapoint;
	if (rf_feature == 0) {
		// get local feature of node u
		info[A1] = f[0];
		info[A2] = f[1];
		info[A3] = f[2];
		info[A4] = f[3];
		info[A5] = f[4];
		info[A6] = f[5];
		info[A11] = f[10];
		info[A12] = f[11];
		info[A13] = f[12];
		info[A14] = f[13];
		info[A15] = f[14];
		info[A16] = f[15];

		info[A1_v] = f[23];
		info[A2_v] = f[24];
		info[A3_v] = f[25];
		info[A4_v] = f[26];
		info[A5_v] = f[27];
		info[A6_v] = f[28];
		info[A11_v] = f[33];
		info[A12_v] = f[34];
		info[A13_v] = f[35];
		info[A14_v] = f[36];
		info[A15_v] = f[37];
		info[A16_v] = f[38];

		datapoint.push_back(info[A1]);
		datapoint.push_back(info[A2]);
		datapoint.push_back(info[A3]);
		datapoint.push_back(info[A4]);
		datapoint.push_back(info[A5]);
		datapoint.push_back(info[A6]);
		datapoint.push_back(info[A11]);
		datapoint.push_back(info[A12]);
		datapoint.push_back(info[A13]);
		datapoint.push_back(info[A14]);
		datapoint.push_back(info[A15]);
		datapoint.push_back(info[A16]);

		datapoint.push_back(info[A1_v]);
		datapoint.push_back(info[A2_v]);
		datapoint.push_back(info[A3_v]);
		datapoint.push_back(info[A4_v]);
		datapoint.push_back(info[A5_v]);
		datapoint.push_back(info[A6_v]);
		datapoint.push_back(info[A11_v]);
		datapoint.push_back(info[A12_v]);
		datapoint.push_back(info[A13_v]);
		datapoint.push_back(info[A14_v]);
		datapoint.push_back(info[A15_v]);
		datapoint.push_back(info[A16_v]);

	} else if (rf_feature == 1) {
		info[A1] = f[0];
		info[A2] = f[1];
		info[A3] = f[2];
		info[A4] = f[3];
		info[A10] = f[9];	//u_2_ynodes;
		info[A11] = f[10];
		info[A12] = f[11];
		info[A13] = f[12];
		info[A14] = f[13];
		info[A15] = f[14];
		info[A16] = f[15];
		info[A17] = f[16];	//0, glob_indeg;

		info[A1_v] = f[23];
		info[A2_v] = f[24];
		info[A3_v] = f[25];
		info[A4_v] = f[26];
		info[A10_v] = f[32]; //u_2_ynodes
		info[A11_v] = f[33];
		info[A12_v] = f[34];
		info[A13_v] = f[35];
		info[A14_v] = f[36];
		info[A15_v] = f[37];
		info[A16_v] = f[38];
		info[A17_v] = f[39];

		datapoint.push_back(info[A1]);
		datapoint.push_back(info[A2]);
		datapoint.push_back(info[A3]);
		datapoint.push_back(info[A4]);
		datapoint.push_back(info[A10]);
		datapoint.push_back(info[A11]);
		datapoint.push_back(info[A12]);
		datapoint.push_back(info[A13]);
		datapoint.push_back(info[A14]);
		datapoint.push_back(info[A15]);
		datapoint.push_back(info[A16]);
		datapoint.push_back(info[A17]);

		datapoint.push_back(info[A1_v]);
		datapoint.push_back(info[A2_v]);
		datapoint.push_back(info[A3_v]);
		datapoint.push_back(info[A4_v]);
		datapoint.push_back(info[A10_v]);
		datapoint.push_back(info[A11_v]);
		datapoint.push_back(info[A12_v]);
		datapoint.push_back(info[A13_v]);
		datapoint.push_back(info[A14_v]);
		datapoint.push_back(info[A15_v]);
		datapoint.push_back(info[A16_v]);
		datapoint.push_back(info[A17_v]);
	}

	return datapoint;
}

void Graph::doNonAdaptiveProbing(vector<int> edges, vector<int> sampledEdges,
		vector<int> RNodes, vector<int> YNodes, int budget,
		map<string, double> coef_map, int ival, bool dFlag, bool logisticFlag,
		int k_path, int c, bool randomForestRegressFlag,
		bool randomForestClassifyFlag, bool rf_feature, string rfModelPath,
		bool BCFlag, bool CCFlag) {

// get list of nodes of gprime
	vector<int> gp_nodes;		// gprime's nodes
	gp_nodes = concat2Vectors(RNodes, YNodes);

	cout << "#nodes in GPrime: " << RNodes.size() + YNodes.size() << endl;
	cout << "budget: " << budget << endl;

	IGraph g;	// load original graph
	g.loadGraphFromEdgeList(edges, dFlag);
	int node_count = g.get_n(); // get #nodes in G

	IGraph g_prime;	// load original graph
	g_prime.loadGraphFromEdgeList(edges, dFlag);

	int rnode_size = RNodes.size();
	int ynode_size = YNodes.size();

	int total_RR_edges = getTotalRedRedEdges(RNodes, sampledEdges);

//get global festure of node u
//	getGlobalFeature(edges, dFlag);
	getGlobalDegree(edges, dFlag);

// get local feature of g
	vector<float> res_indeg, res_outdeg, res_bc, res_cc, res_eig, res_pr,
			res_clc;
	getLocalFeature(g_prime, res_indeg, res_outdeg, res_bc, res_cc, res_eig,
			res_pr, res_clc, dFlag, BCFlag, CCFlag);

	int rnode_indeg = 0;
	int rnode_outdeg = 0;

	int ynode_indeg = 0;
	int ynode_outdeg = 0;

	rnode_indeg = getTotalDegreeOfNodes(RNodes, res_indeg);
	rnode_outdeg = getTotalDegreeOfNodes(RNodes, res_outdeg);

	ynode_indeg = getTotalDegreeOfNodes(YNodes, res_indeg);
	ynode_outdeg = getTotalDegreeOfNodes(YNodes, res_outdeg);

	int total_rnode_deg;
	int total_ynode_deg;
	if (dFlag) {
		total_rnode_deg = rnode_indeg + rnode_outdeg;
		total_ynode_deg = ynode_indeg + ynode_outdeg;
	} else {
		total_rnode_deg = rnode_outdeg;
		total_ynode_deg = ynode_outdeg;
	}
	int m;
	if (total_RR_edges != 0)
		m = pow(total_rnode_deg, 2) / (2 * total_RR_edges); // TODO: calculate m without degree distribution info

	vector<int> du;
	vector<double> du_out;
	map<int, double> u_map;
	map<int, int> id_name_map;

	vector<vector<double>> matrix;
	FOR(i, YNodes.size())
	{
		int u = YNodes[i];

//			cout << "checking node u = " << u << endl;

		std::map<string, double> info_u;
		//get info of node u
		getNodeInfo(u, info_u, rnode_size, ynode_size, rnode_outdeg,
				total_rnode_deg, total_ynode_deg, total_RR_edges, m, res_indeg,
				res_outdeg, res_bc, res_cc, res_eig, res_pr, res_clc,
				res_glob_indeg, res_glob_outdeg, res_glob_bc, res_glob_cc,
				res_glob_eig, res_glob_pr, res_glob_clc, edges, RNodes, YNodes,
				sampledEdges, dFlag);

		if (logisticFlag) {

			// using logistic regression
			//				info[A24] = res_glob_outdeg[u] - res_outdeg[u]; //TODO get alpha from command line

			id_name_map[i] = u; //mapping from nodes' id to internal index

			vector<double> feature; // get features of node u

			feature.push_back(info_u[A1]);
			feature.push_back(info_u[A2]);
			feature.push_back(info_u[A3]);
			feature.push_back(info_u[A4]);

			feature.push_back(info_u[A5]);	//m
			feature.push_back(info_u[A6]);

			feature.push_back(info_u[A7]);	//total_green_edges;
			feature.push_back(info_u[A8]);	//m_prime;
			feature.push_back(info_u[A9]); 	//total_links2_ynodes;
			feature.push_back(info_u[A10]);	//u_2_ynodes;

			feature.push_back(info_u[A11]);
			feature.push_back(info_u[A12]);
			feature.push_back(info_u[A13]);
			feature.push_back(info_u[A14]);
			feature.push_back(info_u[A15]);
			feature.push_back(info_u[A16]);

			feature.push_back(info_u[A17]);	//glob_indeg;
			feature.push_back(info_u[A18]);
			feature.push_back(info_u[A19]);
			feature.push_back(info_u[A20]);
			feature.push_back(info_u[A21]);
			feature.push_back(info_u[A22]);
			feature.push_back(info_u[A23]);

			//print vector f:
			//				cout << "feature vector of node " << u << ":" << endl;
			//				printVector(feature, true);

			matrix.push_back(feature);

		} else {

			//using linear regression
			double total_outcome = 0; // outcome from kpath neighbors

			// estimate outcome of node u
			double outcome_u = estimateOutcome(coef_map, info_u, logisticFlag);
			total_outcome = outcome_u;

//				 cout << "outcome of " << u << ": " << outcome << endl;
			du_out.push_back(total_outcome);
			u_map[u] = total_outcome;
		}
	}

	if (logisticFlag) {

		vector<Pair> pairs;

		// make pair-wise matrix
		// COMPLEXITY: sizeOf(matrix) Choose 2;
		for (int i = 0; i < matrix.size() - 1; i++) {
			vector<double> f1 = matrix[i];
			//				cout << "feature vector 1:" << endl;
			//				printVector(f1, true);
			for (int j = i + 1; j < matrix.size(); j++) {
				vector<double> f2 = matrix[j];
				//					cout << "feature vector 2:" << endl;
				//					printVector(f2, true);

				vector<double> f;
				f.insert(f.end(), f1.begin(), f1.end());
				f.insert(f.end(), f2.begin(), f2.end());

				// get local feature of node u
				std::map<string, double> info;
				info[A1] = f[0];
				info[A2] = f[1];
				info[A3] = f[2];
				info[A4] = f[3];

				info[A5] = f[4];
				info[A6] = f[5];
				info[A7] = f[6];	//0, total_green_edges;
				info[A8] = f[7];	//m_prime;
				info[A9] = f[8]; 	//0, total_links2_ynodes;

				info[A10] = f[9];	//u_2_ynodes;
				info[A11] = f[10];
				info[A12] = f[11];
				info[A13] = f[12];
				info[A14] = f[13];
				info[A15] = f[14];
				info[A16] = f[15];

				info[A17] = f[16];	//0, glob_indeg;
				info[A18] = f[17];
				info[A19] = f[18];
				info[A20] = f[19];
				info[A21] = f[20];
				info[A22] = f[21];
				info[A23] = f[22];

				info[A1_v] = f[23];
				info[A2_v] = f[24];
				info[A3_v] = f[25];
				info[A4_v] = f[26];

				info[A5_v] = f[27];
				info[A6_v] = f[28];
				info[A7_v] = f[29];	//0, total_green_edges;
				info[A8_v] = f[30];	//m_prime;
				info[A9_v] = f[31];	//0, total_links2_ynodes;

				info[A10_v] = f[32]; //u_2_ynodes
				info[A11_v] = f[33];
				info[A12_v] = f[34];
				info[A13_v] = f[35];
				info[A14_v] = f[36];
				info[A15_v] = f[37];
				info[A16_v] = f[38];

				info[A17_v] = f[39]; //0, glob_indeg;
				info[A18_v] = f[40];
				info[A19_v] = f[41];
				info[A20_v] = f[42];
				info[A21_v] = f[43];
				info[A22_v] = f[44];
				info[A23_v] = f[45];

				// estimate outcome with f(.)
				double outcome = estimateOutcome(coef_map, info, logisticFlag); //binary outcome
				f.push_back(outcome);

				//add to pairs
				Pair p;
				p.u_id = i; //u_id: id of node u internally
				p.v_id = j; //v_id: id of node v internally
				p.label = outcome;

				//get real label for u_id, v_id
				if (info[A7] > info[A7_v])
					p.real_label = 1;
				else
					p.real_label = 0;

				pairs.push_back(p);
			}
		}

//		cout << "#pairs: " << pairs.size() << endl;
//		for (int i = 0; i < pairs.size(); i++) {
//			int u_id = pairs[i].u_id;
//			int v_id = pairs[i].v_id;
//			cout << "pairs (" << u_id << "," << v_id << "), estimated label: "
//					<< pairs[i].label << endl;
//			cout << "pairs (" << u_id << "," << v_id << "), real label: "
//					<< pairs[i].real_label << endl;
//		}

		// init. rank_matrix
		vector<vector<int>> rank_matrix; // size of rank_matrix: YNodes.size() * YNodes.size()
		FOR(i, YNodes.size())
		{
			vector<int> v(YNodes.size());
			rank_matrix.push_back(v);
		}

		int error_count = 0;
		int TP = 0;
		int FP = 0;
		int FN = 0;
		int TN = 0;

		FOR(i, pairs.size())
		{
			int u_id = pairs[i].u_id;
			int v_id = pairs[i].v_id;

			//calculate error rate of training model
			if (pairs[i].label != pairs[i].real_label)
				error_count++;

			//calculate TP, FP, FN, TN
			if (pairs[i].label == 1 && pairs[i].real_label == 1)
				TP++;
			else if (pairs[i].label == 0 && pairs[i].real_label == 1)
				FN++;
			else if (pairs[i].label == 1 && pairs[i].real_label == 0)
				FP++;
			else if (pairs[i].label == 0 && pairs[i].real_label == 0)
				TN++;

			//prep. for ranking
			if (pairs[i].label == 0) {
				rank_matrix[u_id][v_id] = 0;
				rank_matrix[v_id][u_id] = 1;
			} else if (pairs[i].label == 1) {
				rank_matrix[u_id][v_id] = 1;
				rank_matrix[v_id][u_id] = 0;
			}
		}

//		cout << "YNode.size: " << YNodes.size() << endl;
//		cout << "#pairs: " << pairs.size() << endl;
//		cout << "error_count: " << error_count << endl;
		float error_rate = (float(error_count) / float(pairs.size()));
		cout << "Error rate (#errors/#data) of this model: " << error_rate
				<< "%" << endl;

		float precision = float(TP) / float(TP + FP);
		float recall = float(TP) / float(TP + FN);
		float f1_score = (2 * (precision * recall)) / (precision + recall);
		cout << "F1 score of this model: " << f1_score << endl;

		vector<int> rank_v;  // size of rank_v: YNodes.size()
		// print rank_matrix
		for (int a = 0; a < rank_matrix.size(); a++) {
			int rank = 0;
			for (int b = 0; b < rank_matrix[a].size(); b++) {
//				cout << rank_matrix[a][b] << " ";
				rank += rank_matrix[a][b];
			}
			// cout << endl;
			rank_v.push_back(rank);
		}

		// print vector rank:
//		cout << "print vector rank: " << endl;
//		printVector(rank_v, true);

		// copy rank_v to tmp
		vector<int> tmp(rank_v.begin(), rank_v.end());

		// sort tmp
		std::sort(tmp.begin(), tmp.end(), std::greater<double>()); //desc sort

		// get top tanks based on budget
		vector<int> top_ranks(tmp.begin(), tmp.begin() + budget);

//		cout << "print vector top_ranks: " << endl;
//		printVector(top_ranks, true);

		//print out id_name_map:
//		cout << "id_name_map: " << endl;
//		std::map<int, int>::iterator it;
//		for (it = id_name_map.begin(); it != id_name_map.end(); ++it) {
//			int idx = it->first;
//			int name = it->second;
//			cout << idx << " " << name << endl;
//		}

		FOR(k, top_ranks.size())
		{
			FOR(j, rank_v.size())
			{
				if (rank_v[j] == top_ranks[k]) {
					int idx = j;
					int u;
					u = id_name_map[idx]; //get name of node from its idx
					if (find(RNodes, u) == false) {

//						cout << "start probing u: " << u << endl;
						probeNode(u, gp_nodes, sampledEdges, RNodes, YNodes, g,
								dFlag);
					}
					if ((k + 1) == 1) {
						cout << "k = " << (k + 1)
								<< ", #Nodes in augmented Gprime: "
								<< gp_nodes.size() << endl;
					} else if ((k + 1) % ival == 0) {
						cout << "k = " << (k + 1)
								<< ", #Nodes in augmented Gprime: "
								<< gp_nodes.size() << endl;
					}
					break;
				}
			}
		}

	} else {
		// linear regression
		// sort YNodes based on node's outcome

		// print map
//		cout << "print map: " << endl;
//		map<int, double>::iterator it;
//		for (it = u_map.begin(); it != u_map.end(); it++) {
//			int u = it->first;
//			double du_out = it->second;
//			cout << "u: " << u << ", du_out: " << du_out << endl;
//		}

		// sort du_out
		std::sort(du_out.begin(), du_out.end(), std::greater<double>()); //desc sort
//		cout << "print du_out: " << endl;
//		printVector(du_out, true);

		vector<double> du_out_sorted(du_out.begin(), du_out.begin() + budget);
//		cout << "print du_out_sorted: " << endl;
//		printVector(du_out_sorted, true);

		// start probing
		int k = 0;
		while (k < budget) {

			bool flag = false;
			map<int, double>::iterator it;
			for (it = u_map.begin(); it != u_map.end(); it++) {

				int u = it->first;
				double du_out = it->second;

				// PROBE node u - which has highest outcome
				if (find(RNodes, u) == false && du_out == du_out_sorted[k]) {

//					cout << "start probing u: " << u << endl;

					flag = true;

					probeNode(u, gp_nodes, sampledEdges, RNodes, YNodes, g,
							dFlag);

					k = k + 1;

					if (k == 1) {
						cout << "k = " << k << ", #Nodes in augmented Gprime: "
								<< gp_nodes.size() << endl;
					} else if (k % ival == 0) {
						cout << "k = " << k << ", #Nodes in augmented Gprime: "
								<< gp_nodes.size() << endl;
					}

				}
			}

			if (flag == false)
				k = budget; // if size(YNodes) < budget: it's time to terminate
		}
	}

//	cout << "YNodes: " << endl;
//	printVector(YNodes, true);
//	cout << "RNodes: " << endl;
//	printVector(RNodes, true);
//	cout << "Nodes in Gprime: " << endl;
//	printVector(gp_nodes, true);

}

void Graph::doNonAdaptiveProbing2(vector<int> edges, vector<int> sampledEdges,
		vector<int> RNodes, vector<int> YNodes, int budget, int ival,
		bool dFlag, string lib_path, string mdl_path, bool BCFlag,
		bool deepWalkFlag, bool CCFlag) {

	cout << "doNonAdaptiveProbing2" << endl;

// get list of nodes of gprime
	vector<int> gp_nodes;		// gprime's nodes
	gp_nodes = concat2Vectors(RNodes, YNodes);

	cout << "#nodes in GPrime: " << RNodes.size() + YNodes.size() << endl;
	cout << "budget: " << budget << endl;

	IGraph g; //load original graph
	g.loadGraphFromEdgeList(edges, dFlag);

	IGraph g_prime;	// load sample graph
	g_prime.loadGraphFromEdgeList(sampledEdges, dFlag);

	// init. Python engine
	PythonWrapper pythonWrapper;
	pythonWrapper.init();
	pythonWrapper.setPath(lib_path);
	pythonWrapper.importModule("Pickle");
	int loaded = 0;
	loaded = pythonWrapper.loadModel("loadModel", mdl_path);
	if (loaded != 1) {
		cerr << "cannot load trained model" << endl;
	} else {
		cout << "trained model loaded!" << endl;
	}

	// get features
	getGlobalFeature(edges, dFlag, false, false);

	vector<float> res_indeg, res_outdeg, res_bc, res_cc, res_eig, res_pr,
			res_clc;
	getLocalFeature(g_prime, res_indeg, res_outdeg, res_bc, res_cc, res_eig,
			res_pr, res_clc, dFlag, BCFlag, CCFlag);

	int rnode_outdeg = 0;
	int total_RR_edges = 0;
	int total_rnode_deg = 0, total_ynode_deg = 0;
	int m = 0;
	getCommonFeatures(sampledEdges, RNodes, YNodes, res_indeg, res_outdeg,
			rnode_outdeg, total_RR_edges, total_rnode_deg, total_ynode_deg, m,
			dFlag);

	// get feature vector for red nodes
	vector<vector<float>> redFeatures;
	FOR(i, RNodes.size())
	{
		vector<float> rFeature;
		Feature fstr;

		int u = RNodes[i];
		getFeatureVectorForRedNodes(rFeature, u, RNodes, YNodes, rnode_outdeg,
				res_indeg, res_outdeg, res_bc, res_cc, res_eig, res_pr, res_clc,
				m, total_rnode_deg, total_RR_edges, total_ynode_deg, g_prime,
				BCFlag, dFlag, CCFlag);

		rFeature.push_back(u); //name of red node corresponds to this feature
		redFeatures.push_back(rFeature);
	}

	// get feature vector for yellow nodes + outcome
	vector<int> du;
	vector<double> du_out;
	map<int, double> u_map;

	FOR(i, YNodes.size())
	{
		vector<float> yFeature;
		int u = YNodes[i];

		cout << "node u = " << u << endl;

		getFeatureVector(redFeatures, yFeature, u, RNodes, YNodes, rnode_outdeg,
				res_indeg, res_outdeg, res_bc, res_cc, res_eig, res_pr, res_clc,
				m, total_rnode_deg, total_RR_edges, total_ynode_deg, g_prime,
				BCFlag, dFlag, deepWalkFlag, CCFlag);

		//select features used for prediction
		vector<double> X;
		if (yFeature.size() == 26) {
			FOR(j, yFeature.size())
			{
				if (deepWalkFlag) {
					if ((j >= 0 && j <= 5) || (j == 10)
							|| (j >= 12 && j <= 15)) {
						X.push_back(yFeature[j]);
					}
				} else {
					//Node centrality
					if ((j >= 0 && j <= 5) || (j == 10)
							|| (j >= 13 && j <= 15)) {
						X.push_back(yFeature[j]);
					}
				}
			}
		} else {
			cout << "size of feature is not correct!" << endl;
			return;
		}

		cout << "X.size: " << X.size() << endl;
		printVector(X, true);

		//using linear regression
		double total_outcome = 0.0; // outcome from kpath neighbors

		// estimate outcome of node u
		double outcome_u = pythonWrapper.predict("predict", X);
		total_outcome = outcome_u;

		cout << "outcome of " << u << ": " << outcome_u << endl;
		du_out.push_back(total_outcome);
		u_map[u] = total_outcome;
	}

	// sort du_out
	std::sort(du_out.begin(), du_out.end(), std::greater<double>()); //desc sort

	vector<double> du_out_sorted(du_out.begin(), du_out.begin() + budget);

	int k = 0;
	vector<int> probed_nodes;
	vector<int> YNodes_org = YNodes;
	if (YNodes_org.size() > budget) {

		while (k < budget) {

			for (int i = 0; i < YNodes_org.size(); i++) {
				//				if (k >= budget)
				//					break;

				double du_out = du_out_sorted[k];

				int u = YNodes_org[i];
				double outcome_u = u_map[u];

				if (du_out == outcome_u && find(RNodes, u) == false) {
					// probe node u
					cout << "start probing u: " << u << endl;

					probed_nodes.push_back(u);

					probeNode(u, gp_nodes, sampledEdges, RNodes, YNodes, g,
							dFlag);

					k = k + 1;

					if (k == 1) {
						cout << "k = " << k << ", #Nodes in augmented Gprime: "
								<< gp_nodes.size() << endl;
					} else if (k % ival == 0) {
						cout << "k = " << k << ", #Nodes in augmented Gprime: "
								<< gp_nodes.size() << endl;
					}
				}
			}
		}
	} else
		cerr << "YNodes.size() < budget" << endl;
}

/**
 *
 */
void Graph::doNonAdaptiveProbingDeepWalk(vector<int> edges,
		vector<int> sampledEdges, vector<int> RNodes, vector<int> YNodes,
		int budget, int ival, map<int, double> node2gain, bool dFlag,
		bool BCFlag) {

	cout << "doNonAdaptiveProbingDeepWalk" << endl;

// get list of nodes of gprime
	vector<int> gp_nodes;		// gprime's nodes
	gp_nodes = concat2Vectors(RNodes, YNodes);

	cout << "#nodes in GPrime: " << RNodes.size() + YNodes.size() << endl;
	cout << "budget: " << budget << endl;

	IGraph g;	// load original graph
	g.loadGraphFromEdgeList(edges, dFlag);
	int node_count = g.get_n(); // get #nodes in G

	int rnode_size = RNodes.size();
	int ynode_size = YNodes.size();

	vector<int> du;
	vector<double> du_out;
	map<int, double> u_map;

	FOR(i, YNodes.size())
	{
		int u = YNodes[i];
		cout << "u = " << u << endl;

// estimate outcome of node u
		double outcome_u = node2gain[u];

		cout << "outcome of " << u << ": " << outcome_u << endl;
		du_out.push_back(outcome_u);
		u_map[u] = outcome_u;
	}

//	printMap(node2gain);

// sort du_out
	std::sort(du_out.begin(), du_out.end(), std::greater<double>()); //desc sort

	vector<double> du_out_sorted(du_out.begin(), du_out.begin() + budget);

	cout << "du_out_sorted: " << endl;
	printVector(du_out_sorted, false);

// start probing
	int k = 0;
	vector<int> probed_nodes;
	vector<int> YNodes_org = YNodes;
	if (YNodes_org.size() > budget) {

		while (k < budget) {

			for (int i = 0; i < YNodes_org.size(); i++) {

				double du_out = du_out_sorted[k];

				int u = YNodes_org[i];
				double outcome_u = u_map[u];

				if (du_out == outcome_u && find(RNodes, u) == false) {
					// probe node u
					cout << "start probing u: " << u << endl;

					probed_nodes.push_back(u);

					probeNode(u, gp_nodes, sampledEdges, RNodes, YNodes, g,
							dFlag);

					k = k + 1;

					if (k == 1) {
						cout << "k = " << k << ", #Nodes in augmented Gprime: "
								<< gp_nodes.size() << endl;
					} else if (k % ival == 0) {
						cout << "k = " << k << ", #Nodes in augmented Gprime: "
								<< gp_nodes.size() << endl;
					}
				}
			}
		}
	} else
		cerr << "YNodes.size() < budget" << endl;

//		cout << "probed nodes in the order: " << endl;
//		printVector(probed_nodes, true);

//	cout << "YNodes: " << endl;
//	printVector(YNodes, true);do
//	cout << "RNodes: " << endl;
//	printVector(RNodes, true);
//	cout << "Nodes in Gprime: " << endl;
//	printVector(gp_nodes, true);

}

// probe in the order:
//void Graph::doNonAdaptiveProbingDeepWalk(vector<int> edges,
//		vector<int> sampledEdges, vector<int> RNodes, vector<int> YNodes,
//		int budget, int ival, map<int, double> node2gain, bool dFlag,
//		bool BCFlag) {
//
//	cout << "doNonAdaptiveProbingDeepWalk" << endl;
//
//// get list of nodes of gprime
//	vector<int> gp_nodes;		// gprime's nodes
//	gp_nodes = concat2Vectors(RNodes, YNodes);
//
//	cout << "#nodes in GPrime: " << RNodes.size() + YNodes.size() << endl;
//	cout << "budget: " << budget << endl;
//
//	IGraph g;	// load original graph
//	g.loadGraphFromEdgeList(edges, dFlag);
//	int node_count = g.get_n(); // get #nodes in G
//
//	int rnode_size = RNodes.size();
//	int ynode_size = YNodes.size();
//
//	vector<int> du;
//	vector<double> du_out;
//	map<int, double> u_map;
//
//	FOR(i, YNodes.size())
//	{
//		int u = YNodes[i];
//		cout << "u = " << u << endl;
//
//// estimate outcome of node u
//		double outcome_u = node2gain[u];
//
//		cout << "outcome of " << u << ": " << outcome_u << endl;
//		du_out.push_back(outcome_u);
//		u_map[u] = outcome_u;
//	}
//
////	printMap(node2gain);
//
//// sort du_out
//	std::sort(du_out.begin(), du_out.end(), std::greater<double>()); //desc sort
//
//	vector<double> du_out_sorted(du_out.begin(), du_out.begin() + budget);
//
//// start probing
//	int k = 0;
//	vector<int> probed_nodes;
//	while (k < budget) {
//
//		bool flag = false;
//		map<int, double>::iterator it;
//		for (it = u_map.begin(); it != u_map.end(); it++) {
//
//			int u = it->first;
//			double du_out = it->second;
//
//			// PROBE node u - which has highest outcome
//			if (find(RNodes, u) == false && du_out == du_out_sorted[k]) {
//
//				cout << "start probing u: " << u << endl;
//
//				flag = true;
//
//				probed_nodes.push_back(u);
//				probeNode(u, gp_nodes, sampledEdges, RNodes, YNodes, g, dFlag);
//
//				k = k + 1;
//
//				if (k == 1) {
//					cout << "k = " << k << ", #Nodes in augmented Gprime: "
//							<< gp_nodes.size() << endl;
//				} else if (k % ival == 0) {
//					cout << "k = " << k << ", #Nodes in augmented Gprime: "
//							<< gp_nodes.size() << endl;
//				}
//			}
//		}
//
//		if (flag == false)
//			k = budget; // if size(YNodes) < budget: it's time to terminate
//	}
//
////		cout << "probed nodes in the order: " << endl;
////		printVector(probed_nodes, true);
//
////	cout << "YNodes: " << endl;
////	printVector(YNodes, true);do
////	cout << "RNodes: " << endl;
////	printVector(RNodes, true);
////	cout << "Nodes in Gprime: " << endl;
////	printVector(gp_nodes, true);
//
//}

void Graph::getNodeInfo(int u, std::map<string, double> &info, int rnode_size,
		int ynode_size, int rnode_outdeg, int total_rnode_deg,
		int total_ynode_deg, int total_RR_edges, int m,
		vector<float> &res_indeg, vector<float> &res_outdeg,
		vector<float> &res_bc, vector<float> &res_cc, vector<float> &res_eig,
		vector<float> &res_pr, vector<float> &res_clc,
		vector<float> &res_glob_indeg, vector<float> &res_glob_outdeg,
		vector<float> &res_glob_bc, vector<float> &res_glob_cc,
		vector<float> &res_glob_eig, vector<float> &res_glob_pr,
		vector<float> &res_glob_clc, vector<int> &edges, vector<int> &RNodes,
		vector<int> &YNodes, vector<int> &sampledEdges, bool dFlag) {

//	cout << "** 4-1" << endl;

	float local_indeg = 0, local_outdeg = 0, local_bc = 0.0, local_cc = 0.0,
			local_eig = 0.0, local_pr = 0.0, local_clc = 0.0;
	float glob_indeg = 0.0, glob_outdeg = 0.0, glob_bc = 0.0, glob_cc = 0.0,
			glob_eig = 0.0, glob_pr = 0.0, glob_clc = 0.0;

	int total_green_edges = 0, total_links2_rnodes = 0, total_links2_ynodes = 0;

//global node centrality
	local_indeg = res_indeg[u];
	local_outdeg = res_outdeg[u];
	local_bc = res_bc[u];
	local_cc = res_cc[u];
	local_eig = res_eig[u];
	local_pr = res_pr[u];
	local_clc = res_clc[u];

//global node centrality
	glob_indeg = res_glob_indeg[u];
	glob_outdeg = res_glob_outdeg[u];
//	glob_bc = res_glob_bc[u];
//	glob_cc = res_glob_cc[u];
//	glob_eig = res_glob_eig[u];
//	glob_pr = res_glob_pr[u];
//	glob_clc = res_glob_clc[u];

//total_green_edges = getTotalColoredEdges(u, edges, RNodes, YNodes, dFlag); //option 1: #links from u to green nodes (degree dist. assumed
//	total_green_edges = getTotalLinks2GreenNodes(u, RNodes, YNodes);

//	cout << "** 4-2" << endl;
	total_links2_rnodes = getTotalLinks2Nodes(u, RNodes, dFlag);
//	cout << "** 4-3" << endl;
//	total_links2_rnodes = getTotalColoredEdges(u, sampledEdges, RNodes, dFlag); //option 2: #links from u to red nodes

//total_links2_ynodes = getTotalColoredEdges(u, edges, YNodes, dFlag); //option 3: #links from u to yellow nodes (degree dist. assumed)

	int m_prime; // TODO: calculate m with degree distribution info
//	m_prime = (res_glob_outdeg[u] * total_rnode_deg)
//			/ (total_links2_rnodes * 2);
	m_prime = getNumberOfEdges();

	int k;
	if (total_RR_edges != 0)
		k = (total_links2_rnodes * total_rnode_deg) / total_RR_edges;

	int u_2_ynodes = (total_ynode_deg * total_links2_rnodes) / total_rnode_deg; // estimate # links from u to yellow nodes

// get local feature of node u
// TODO: to be updated if needed

	info[A1] = rnode_size;
	info[A2] = ynode_size;
	info[A3] = rnode_outdeg;
	info[A4] = total_RR_edges;
	info[A5] = m;

	info[A6] = k;
	info[A7] = total_green_edges; //total_green_edges;
	info[A8] = m_prime; //m_prime;
	info[A9] = total_links2_ynodes; //total_links2_ynodes;

	info[A10] = u_2_ynodes; //u_2_ynodes;
	info[A11] = local_outdeg;
	info[A12] = local_bc;
	info[A13] = local_cc;
	info[A14] = local_eig;
	info[A15] = local_pr;
	info[A16] = local_clc;

	info[A17] = glob_indeg;
	info[A18] = glob_outdeg;
	info[A19] = glob_bc;
	info[A20] = glob_cc;
	info[A21] = glob_eig;
	info[A22] = glob_pr;
	info[A23] = glob_clc;

//TODO to be removed
//info[A24] = local_pr * local_clc;

}

set<int> Graph::getNeighbors(int u, IGraph &g, bool dFlag) {

	set<int> nb_2hops;

	vector<int> u_nb = g.getNeighbors(u, dFlag);
//cout << "neighbor of u: ";
	FOR(i, u_nb.size())
	{
		int v = u_nb[i];
		//cout << v << " ";
		nb_2hops.insert(v);

		//get nb of v:
		vector<int> v_nb = g.getNeighbors(v, dFlag);
		v_nb.erase(std::remove(v_nb.begin(), v_nb.end(), u), v_nb.end()); //nb of v can be u, so remove u

		FOR(j, v_nb.size())
		{
			int w = v_nb[j];
			nb_2hops.insert(w);
		}
	}

//	cout << endl;
	return nb_2hops;

}

float Graph::totalDegree2Hops(set<int> nodes, vector<int> YNodes,
		vector<int> RNodes, vector<float> degree) {
	set<int>::iterator it;
	float totalDeg = 0.0;
	for (it = nodes.begin(); it != nodes.end(); it++) {
		int u = *it;

		if (find(RNodes, u)) {
			totalDeg += degree[u];
		} else if (find(YNodes, u)) {
			totalDeg += degree[u] * 1.1;
		}
	}
	return totalDeg;

}

set<int> Graph::getNeighbors_kpath(int s, int node_count, int k_path,
		vector<int> RNodes, IGraph &g, vector<int> &distance, bool dFlag) {

	set<int> kpath_nbors_of_s;

//	cout << "kpath starts from " << s << endl;

//init vector
	vector<int> colors;
	FOR(i, node_count)
	{
		colors.push_back(1);
		distance.push_back(-1);
	}

//	kpath_nbors_of_s.insert(s); //add u to kpath_nbors_of_u, in case k_path == 0, size(kpath_nbors_of_u) == 1
	colors[s] = 2; //gray nodes in BFS
	distance[s] = 0;

	queue<int> q;
	q.push(s);

	while (q.empty() == false) {

		int s = q.front();
		q.pop();

		vector<int> s_nbors = g.getNeighbors(s, dFlag);
		FOR(j, s_nbors.size())
		{
			int v = s_nbors[j];
			if (colors[v] == 1 && find(RNodes, v) == false) {
				colors[v] = 2;
				distance[v] = distance[s] + 1;
				if (distance[v] <= k_path) {
					q.push(v);
					kpath_nbors_of_s.insert(v);
				}
			}
		}

		colors[s] = 3;
	}

//print out kpath_nbors_of_u:
//	cout << "kpath_nbors_of " << s << endl;
//	for (set<int>::iterator it = kpath_nbors_of_s.begin();
//			it != kpath_nbors_of_s.end(); ++it) {
//		cout << *it << " ";
//	}
//	cout << endl;

	return kpath_nbors_of_s;

}

void Graph::probeNode(int u, vector<int> &gp_nodes, vector<int> &sampledEdges,
		vector<int> &RNodes, vector<int> &YNodes, IGraph &g, bool dFlag) {

//	cout << "start probing u: " << u << endl;

// traverse all neighbors v of u
//	vector<int> u_nbors = g.getNeighbors(u, dFlag);
	vector<int> u_nbors = getNeighbors(u);

	cout << "size of neibors of u: " << u_nbors.size() << endl;
// printVector(u_nbors, false);

	int count = 0;
	FOR(j, u_nbors.size())
	{
		int v = u_nbors[j];
		// check if v exists in G' or not, if not, then add v to g
		if (find(gp_nodes, v) == false) {
			gp_nodes.push_back(v);
			count++;
// cout << "added " << v << " to Gprime" << endl;
		}

		// add edges (u-v) to g
		if (checkEdgeExist(sampledEdges, u, v, dFlag) == false) {
			sampledEdges.push_back(u);
			sampledEdges.push_back(v);
//	cout << " added edge: " << u << " " << v
//	<< " to Gprime " << endl;
		}

// if v not in RNodes and YNodes, then add v to yellow nodes (candidate nodes)
		if (find(RNodes, v) == false && find(YNodes, v) == false)
			YNodes.push_back(v);
	}
	cout << "added " << count << " nodes to " << u << endl;

// add node u to RNodes
	if (find(RNodes, u) == false)
		RNodes.push_back(u);

// remove node u from YNodes
	removeElemInVector(u, YNodes);

//	cout << "done probing " << u << "!" << endl;
//	cout << "YNodes: " << endl;
//	printVector(YNodes, true);
//	cout << "RNodes: " << endl;
//	printVector(RNodes, true);
//	cout << "Nodes in Gprime: " << endl;
//	printVector(gp_nodes, true);

}

void Graph::removeElemInVector(int elem, vector<int> &v) {

	int idx;
	FOR(i, v.size())
	{
		if (v[i] == elem) {
			idx = i;
			break;
		}
	}
	v.erase(v.begin() + idx);

}

vector<vector<float>> Graph::getFeatureMatrix() {
	return fMatrix;
}

vector<vector<float>> Graph::getPairWiseFeatureMatrix() {
	return fMatrix_pairwise;
}

bool Graph::checkEdgeExist(vector<int> edges, int u, int v, bool dFlag) {

//	cout << "check edge exist: " << u << " " << v << endl;
	for (int i = 0; i < edges.size(); i += 2) {
		if (dFlag && edges[i] == u && edges[i + 1] == v) {
			return true;
		} else if (dFlag == false
				&& ((edges[i] == u && edges[i + 1] == v)
						|| (edges[i] == v && edges[i + 1] == u))) {
			return true;
		}
	}
	return false;

}

int Graph::getTotalDegreeOfNodes(vector<int> nodes, vector<float> degV) {
	int total = 0;

	int idx;
	FOR(i, nodes.size())
	{
		idx = nodes[i];
		total = total + degV[idx];
	}

	return total;
}

int Graph::getTotalRedRedEdges(const vector<int> RNodes,
		const vector<int> sampledEdges) {
//	printVector(sampledEdges, false);
	int total = 0;

	for (unsigned int i = 0; i < sampledEdges.size(); i = i + 2) {
		int u = sampledEdges[i];
		int v = sampledEdges[i + 1];

		if (find(RNodes, u) == true && find(RNodes, v) == true) {
			total++; //red-red edges found
		}
	}
	return total;
}

int Graph::getTotalColoredEdges(const int u, const vector<int> edges,
		const vector<int> RNodes, const vector<int> YNodes, const bool dFlag) {

	int total = 0;

	for (unsigned int i = 0; i < edges.size(); i += 2) {

		if (dFlag == false) {
//UNDIRECTED GRAPH
			if ((edges[i] == u && find(RNodes, edges[i + 1]) == false
					&& find(YNodes, edges[i + 1]) == false)
					|| (edges[i + 1] == u && find(RNodes, edges[i]) == false
							&& find(YNodes, edges[i]) == false)) {
				//u links to green node only (not in red or yellow set)
//				cout << edges[i] << " " << edges[i + 1] << endl;
				total++; //edge found
			}
		} else {
//DIRECTED GRAPH
			if (edges[i] == u && find(RNodes, edges[i + 1]) == false
					&& find(YNodes, edges[i + 1]) == false) {
				//u points to green node only (not in red or yellow set)
//				cout << edges[i] << " " << edges[i + 1] << endl;
				total++; //edge starting from node u found
			}
		}

	}
//	cout << endl;
	return total;
}

int Graph::getTotalColoredEdges(const int u, const vector<int> edges,
		const vector<int> nodes, const bool dFlag) {

	int total = 0;

	for (unsigned int i = 0; i < edges.size(); i += 2) {

		if (dFlag == false) {
//UNDIRECTED GRAPH
//			if ((edges[i] == u && find(nodes, edges[i + 1]) == true)
//					|| (edges[i + 1] == u && find(nodes, edges[i]) == true)) {
//u links to red/yellow node only(not in yellow/red set)
			if (edges[i] == u && find(nodes, edges[i + 1]) == true) {
//				cout << edges[i] << " " << edges[i + 1] << endl;
				total++;			//edge found
				continue;
			}
			if (edges[i + 1] == u && find(nodes, edges[i]) == true) {
				total++; //edge found
				continue;
			}

		} else {
//DIRECTED GRAPH
			if (edges[i] == u && find(nodes, edges[i + 1]) == true) {
				//u points to red/yellow node only (not in yellow/red set)
//				cout << edges[i] << " " << edges[i + 1] << endl;
				total++; //edge starting from node u found
			}
		}

	}
//	cout << endl;
	return total;
}

bool Graph::contains(queue<int> Q, int e) {
	while (!Q.empty()) {
		int a = Q.front();
		Q.pop();
		if (a == e)
			return true;
	}
	return false;
}

vector<int> Graph::concat2Vectors(vector<int> A, vector<int> B) {

	vector<int> C;
	C.reserve(A.size() + B.size());
	C.insert(C.end(), A.begin(), A.end());
	C.insert(C.end(), B.begin(), B.end());

	return C;
}

void Graph::getInDegree(vector<int> &inDeg, vector<vector<int>> aList) {

	FOR(i, inDeg.size())
	{
		inDeg[i] = 0;
	}

	FOR(i, nodes.size())
	{
		FOR(j, aList.size())
		{
			FOR(k, aList[j].size())
			{
				if (aList[j][k] == i) {
					inDeg[i] = inDeg[i] + 1;
					break;
				}
			}
		}
	}
}

void Graph::generateWeights(vector<vector<float>> &w) {

//generate weight 0~1 for all edges in graph
	FOR(i, head.size())
	{
		float rw = getRN01();

		int u = head[i];
		int v = tail[i];
		w[u].push_back(rw);

//		cout << "u-v-w: " << u << " " << v << " " << rw << endl;

		if (directedFlag == false) {
			w[v].push_back(rw);
//			cout << "v-u-w: " << v << " " << u << " " << rw << endl;
		}
	}

}

void Graph::getOutDegree(vector<int> &outDeg, vector<vector<int>> aList) {
	FOR(i, aList.size())
	{
		outDeg[i] = aList[i].size();
	}
}

void Graph::printQueue(queue<int> Q) {
	while (!Q.empty()) {
		cout << Q.front() << " ";
		Q.pop();
	}
	cout << endl;
}

void Graph::removeDuplicatedEdges(vector<int> &edges) {

	vector<int> idx; // get ids of edges to be removed
	for (int i = 0; i < edges.size(); i = i + 2) {
		for (int j = i + 2; j < edges.size(); j = j + 2) {
			if (edges[i] == edges[j + 1] && edges[i + 1] == edges[j]) {
				idx.push_back(j);
				idx.push_back(j + 1);
			}
		}
	}

// Sort the vector using predicate and std::sort
	std::sort(idx.begin(), idx.end(), std::less<int>());

	for (int i = 0; i < idx.size(); i++) {
//		cout << idx[i] << " ";
		edges.erase(edges.begin() + idx[i] - i);
	}

}

//get random number in the range of 0~1
float Graph::getRN01() {
//cout << "Seed: " << time(NULL) << endl;
	float val = (float) rand() / (float) RAND_MAX;
	return val;

}

void Graph::resetColors() {
	FOR(i, colors.size())
	{
		colors[i] = 1;
	}
}

vector<int> Graph::getNodes() {
	return nodes;
}

vector<int> Graph::convertAdjList2EdgeList(vector<vector<int>> aList) {

	vector<int> eList;
	FOR(i, aList.size())
	{
		FOR(j, aList[i].size())
		{
			eList.push_back(i);
			eList.push_back(aList[i][j]);
		}
	}
	return eList;

}

void Graph::printEdgeList(vector<int> eList) {
	for (int i = 0; i < eList.size(); i = i + 2)
		cout << eList[i] << " " << eList[i + 1] << endl;
	cout << endl;
}

void Graph::printSet(set<int> s) {

	set<int>::iterator it;
	for (it = s.begin(); it != s.end(); it++) {
		cout << *it << " ";
	}
	cout << endl;
}

void Graph::printAdjList() {

	FOR (i, adjList.size())
	{
		cout << "List of adjacent nodes of head " << i << ": ";
		for (unsigned int j = 0; j < adjList[i].size(); j++) {
			cout << adjList[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl << endl;
}

template<class T>
void Graph::printMatrix(vector<vector<T>> matrix) {

	FOR (i, matrix.size())
	{
		cout << i << ": ";
		for (unsigned int j = 0; j < matrix[i].size(); j++) {
			cout << matrix[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl << endl;
}

template<class T>
void Graph::printVector(vector<T> v, bool flag) {

	cout << "vector's size = " << v.size() << endl;
	if (flag) {
		FOR (i, v.size())
		{
			cout << i << " ";
		}
		cout << endl;
	}

	FOR (i, v.size())
	{
		cout << v[i] << " ";
	}
	cout << endl;
	cout << endl;
}

void Graph::printMap(map<int, double> map) {
	std::map<int, double>::iterator it;
	for (it = map.begin(); it != map.end(); ++it) {
		cout << it->first << "-" << it->second << endl;
	}
}

//template<typename K, typename V>
//void Graph::printMapping(const std::map<K,V>& mapping) {
//	typename std::map<K,V>::const_iterator it;
//	for (it = mapping.begin(); it != mapping.end(); ++it) {
//		cout << it->first << "-" << it->second << endl;
//	}
//}

void Graph::write2File(vector<vector<float>> matrix, string fname) {

	cout << "write feature matrix to file..." << endl;

//	ofstream o(fname);
	ofstream ofs(fname, std::ofstream::app); // append mode

	FOR (i, matrix.size())
	{
		for (unsigned int j = 0; j < matrix[i].size(); j++) {
			ofs << matrix[i][j];
			if (j < matrix[i].size() - 1)
				ofs << ",";
//power_edgelist
		}
		ofs << endl;
	}
	cout << "wrote feature matrix to: " << fname << endl;

//close stream
	ofs.close();
}

void Graph::write2File(string fname, vector<float> v) {

	cout << "write sampled graph: " << fname << endl;

	ofstream o(fname);

	o << v.size();
	o << endl;

	for (int i = 0; i < v.size();) {
//		if (flag == true) {
//			o << v[i] << " " << v[i + 1];
//			o << endl;
//			i = i + 2;
//		} else {
		o << v[i] << endl;
		i = i + 1;
//		}
	}
	cout << "wrote vector to: " << fname << endl;

//close stream
	o.close();

}

void Graph::write2File(string fname, vector<int> v) {

	cout << "write sampled graph: " << fname << endl;

	ofstream o(fname);

	o << v.size();
	o << endl;

	for (int i = 0; i < v.size();) {
//		if (flag == true) {
//			o << v[i] << " " << v[i + 1];
//			o << endl;
//			i = i + 2;
//		} else {
		o << v[i] << " ";
		i = i + 1;
//		}
	}
	cout << "wrote vector to: " << fname << endl;

//close stream
	o.close();

}

void Graph::write2File(string fname, vector<int> v, vector<int> rnodes,
		vector<int> ynodes, bool headerFlag) {
//header is needed in default; the one used with deepwalk having no header

	cout << "write sampled graph: " << fname << endl;

	ofstream o(fname);

	int sampled_nodes = rnodes.size() + ynodes.size();

	if (headerFlag) {
		o << sampled_nodes << " ";
		o << v.size() / 2;
		o << endl;
	}

	for (int i = 0; i < v.size(); i = i + 2) {
//		if (flag == true) {
		o << v[i] << " " << v[i + 1];
		o << endl;
//		} else {
//			o << v[i] << " ";
//			i = i + 1;
//		}
	}
	cout << "wrote vector to: " << fname << endl;

//close stream
	o.close();

}

//void Graph::doGreedyNonAdaptiveProbing(vector<int> edges,
//		vector<int> sampledEdges, vector<int> RNodes, vector<int> YNodes,
//		int budget, int ival, bool dFlag) {
//
//	cout << "doGreedyNonAdaptiveProbing" << endl;
//
//	// get list of nodes of gprime
//	vector<int> gp_nodes;		// gprime's nodes
//	gp_nodes = concat2Vectors(RNodes, YNodes);
//	cout << "gp_nodes's size: " << gp_nodes.size() << endl;
//
//	cout << "#nodes in GPrime: " << RNodes.size() + YNodes.size() << endl;
//	cout << "budget: " << budget << endl;
//
//	cout << "GREEDY NON-ADAPTIVE PROBING" << endl;
////	cout << "GNode's size: " << GNodes.size() << endl;
//	cout << "RNode's size: " << RNodes.size() << endl;
//	cout << "YNode's size: " << YNodes.size() << endl;
//
//	IGraph g;		// load original graph
//	g.loadGraphFromEdgeList(edges, dFlag);
//
//	vector<double> du_out;
//	map<int, double> u_map;
//
//	FOR(i, YNodes.size())
//	{
//		int u = YNodes[i];
//		int total_green_edges = 0;
//
//		//cout << "finding green links... for " << u << endl;
//
//		total_green_edges = getTotalLinks2GreenNodes(u, RNodes, YNodes);
//		//	total_green_edges = getTotalColoredEdges(u, edges, RNodes, YNodes, dFlag); //option 1: #links from u to green nodes (degree dist. assumed)
//
//		int outcome = total_green_edges;
//		du_out.push_back(outcome);
//		u_map[u] = outcome;
//	}
//
//	// sort du_out
//	std::sort(du_out.begin(), du_out.end(), std::greater<double>());
//	vector<double> du_out_sorted(du_out.begin(), du_out.begin() + budget);
//
//	//start probing
//	int k = 0;
//	while (k < budget) {
//
////		cout << "GNode's size: " << GNodes.size() << endl;
//
//		map<int, double>::iterator it;
//		for (it = u_map.begin(); it != u_map.end(); it++) {
//
//			int u = it->first;
//			double du_out = it->second;
//
//			// PROBE node u
//			if (find(RNodes, u) == false && du_out == du_out_sorted[k]) {
//
////				cout << "start probing u: " << u << endl;
//
//				probeNode(u, gp_nodes, sampledEdges, RNodes, YNodes, g, dFlag);
//
//				k = k + 1;
//				if (k == 1) {
//					cout << "k = " << k << ", #Nodes in augmented Gprime: "
//							<< gp_nodes.size() << endl;
//				} else if (k % ival == 0) {
//					cout << "k = " << k << ", #Nodes in augmented Gprime: "
//							<< gp_nodes.size() << endl;
//				}
//			}
//		}
//	}
//}

//// doNonAdaptiveProbing for node centrality - (don't use python model)
//void Graph::doNonAdaptiveProbing3(vector<int> edges, vector<int> sampledEdges,
//		vector<int> RNodes, vector<int> YNodes, int budget,
//		map<string, double> coef_map, int ival, bool dFlag, bool BCFlag,
//		bool CCFlag) {
//
//	cout << "doNonAdaptiveProbing3" << endl;
//
//	// get list of nodes of gprime
//	vector<int> gp_nodes;		// gprime's nodes
//	gp_nodes = concat2Vectors(RNodes, YNodes);
//
//	cout << "#nodes in GPrime: " << RNodes.size() + YNodes.size() << endl;
//	cout << "budget: " << budget << endl;
//
//	IGraph g;	// load original graph
//	g.loadGraphFromEdgeList(edges, dFlag);
//	int node_count = g.get_n(); // get #nodes in G
//
//	IGraph g_prime;	// load sampled graph
//	g_prime.loadGraphFromEdgeList(sampledEdges, dFlag);
//
//	int rnode_size = RNodes.size();
//	int ynode_size = YNodes.size();
//
//	int total_RR_edges = getTotalRedRedEdges(RNodes, sampledEdges);
//
//	// get global feature of node u
//	getGlobalDegree(edges, dFlag);
//
//	// get local feature of g
//	vector<float> res_indeg, res_outdeg, res_bc, res_cc, res_eig, res_pr,
//			res_clc;
//	getLocalFeature(g_prime, res_indeg, res_outdeg, res_bc, res_cc, res_eig,
//			res_pr, res_clc, dFlag, BCFlag, CCFlag);
//
//	int rnode_indeg = 0;
//	int rnode_outdeg = 0;
//
//	int ynode_indeg = 0;
//	int ynode_outdeg = 0;
//
//	rnode_indeg = getTotalDegreeOfNodes(RNodes, res_indeg);
//	rnode_outdeg = getTotalDegreeOfNodes(RNodes, res_outdeg);
//
//	ynode_indeg = getTotalDegreeOfNodes(YNodes, res_indeg);
//	ynode_outdeg = getTotalDegreeOfNodes(YNodes, res_outdeg);
//
//	int total_rnode_deg;
//	int total_ynode_deg;
//	if (dFlag) {
//		total_rnode_deg = rnode_indeg + rnode_outdeg;
//		total_ynode_deg = ynode_indeg + ynode_outdeg;
//	} else {
//		total_rnode_deg = rnode_outdeg;
//		total_ynode_deg = ynode_outdeg;
//	}
//	int m;
//	if (total_RR_edges != 0)
//		m = pow(total_rnode_deg, 2) / (2 * total_RR_edges); // TODO: calculate m without degree distribution info
//
//	vector<int> du;
//	vector<double> du_out;
//	map<int, double> u_map;
//	map<int, int> id_name_map;
//
//	// print map
//	cout << "print coef map: " << endl;
//	map<string, double>::iterator it;
//	for (it = coef_map.begin(); it != coef_map.end(); it++) {
//		string feature = it->first;
//		double value = it->second;
//		cout << "feature: " << feature << ", value: " << value << endl;
//	}
//
//	cout << "res_indeg: " << endl;
//	printVector(res_indeg, true);
//
//	//start ranking
//	FOR(i, YNodes.size())
//	{
//		int u = YNodes[i];
//
//		cout << "checking node u = " << u << endl;
//
//		std::map<string, double> info_u;
//		//get info of node u
//		getNodeInfo(u, info_u, rnode_size, ynode_size, rnode_outdeg,
//				total_rnode_deg, total_ynode_deg, total_RR_edges, m, res_indeg,
//				res_outdeg, res_bc, res_cc, res_eig, res_pr, res_clc,
//				res_glob_indeg, res_glob_outdeg, res_glob_bc, res_glob_cc,
//				res_glob_eig, res_glob_pr, res_glob_clc, edges, RNodes, YNodes,
//				sampledEdges, dFlag);
//
//		//using linear regression
//		double total_outcome = 0; // outcome from kpath neighbors
//
//		// estimate outcome of node u
//		double outcome_u = estimateOutcome(coef_map, info_u, false);
//		total_outcome = outcome_u;
//
//		cout << "local degree of " << u << ": " << res_indeg[u] << endl;
//		cout << "pagerank of " << u << ": " << res_pr[u] << endl;
//		cout << "eig. value of " << u << ": " << res_eig[u] << endl;
//		cout << "closeness of " << u << ": " << res_cc[u] << endl;
//		cout << "clus. coefs of " << u << ": " << res_clc[u] << endl;
//		cout << "glb degree of " << u << ": " << res_glob_indeg[u] << endl;
//		cout << "--outcome of " << u << ": " << outcome_u << endl;
//
//		du_out.push_back(total_outcome);
//		u_map[u] = total_outcome;
//	}
//
//	// linear regression
//	// sort YNodes based on node's outcome
//
//	// print map
////		cout << "print map: " << endl;
////		map<int, double>::iterator it;
////		for (it = u_map.begin(); it != u_map.end(); it++) {
////			int u = it->first;
////			double du_out = it->second;
////			cout << "u: " << u << ", du_out: " << du_out << endl;
////		}
//
//	cout << "print du_out: " << endl;
//	printVector(du_out, true);
//
//	// sort du_out
//	std::sort(du_out.begin(), du_out.end(), std::greater<double>()); //desc sort
//
//	vector<double> du_out_sorted(du_out.begin(), du_out.begin() + budget);
//	cout << "print du_out_sorted: " << endl;
//	printVector(du_out_sorted, true);
//
//	// start probing
//	int k = 0;
//	vector<int> probed_nodes;
//	vector<int> YNodes_org = YNodes;
//	if (YNodes_org.size() > budget) {
//
//		while (k < budget) {
//
//			for (int i = 0; i < YNodes_org.size(); i++) {
//
//				double du_out = du_out_sorted[k];
//
//				int u = YNodes_org[i];
//				double outcome_u = u_map[u];
//
//				if (du_out == outcome_u && find(RNodes, u) == false) {
//					// probe node u
//					cout << "start probing u: " << u << endl;
//
//					probed_nodes.push_back(u);
//
//					probeNode(u, gp_nodes, sampledEdges, RNodes, YNodes, g,
//							dFlag);
//
//					k = k + 1;
//
//					if (k == 1) {
//						cout << "k = " << k << ", #Nodes in augmented Gprime: "
//								<< gp_nodes.size() << endl;
//					} else if (k % ival == 0) {
//						cout << "k = " << k << ", #Nodes in augmented Gprime: "
//								<< gp_nodes.size() << endl;
//					}
//				}
//			}
//		}
//	} else
//		cerr << "YNodes.size() < budget" << endl;
//
////	// start probing
////	int k = 0;
////	while (k < budget) {
////
////		bool flag = false;
////		map<int, double>::iterator it;
////		for (it = u_map.begin(); it != u_map.end(); it++) {
////
////			int u = it->first;
////			double du_out = it->second;
////
////			// PROBE node u - which has highest outcome
////			if (find(RNodes, u) == false && du_out == du_out_sorted[k]) {
////
////				cout << "start probing u: " << u << endl;
////
////				flag = true;
////
////				probeNode(u, gp_nodes, sampledEdges, RNodes, YNodes, g, dFlag);
////
////				k = k + 1;
////
////				if (k == 1) {
////					cout << "k = " << k << ", #Nodes in augmented Gprime: "
////							<< gp_nodes.size() << endl;
////				} else if (k % ival == 0) {
////					cout << "k = " << k << ", #Nodes in augmented Gprime: "
////							<< gp_nodes.size() << endl;
////				}
////
////			}
////		}
////
////		if (flag == false)
////			k = budget; // if size(YNodes) < budget: it's time to terminate
//
//}
