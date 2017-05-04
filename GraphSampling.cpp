#include "IGraph.h"
#include "graph.h"
#include "option.h"
#include "omp.h"
#include <time.h>
#include <sys/time.h>

template<class T>
inline vector<T> iVector_to_vector(igraph_vector_t *v) {
	int s = igraph_vector_size(v);
	vector<T> result(s);
	for (int i = 0; i < s; i++)
		result[i] = VECTOR(*v)[i];
	return result;
}

//void print_vector(igraph_vector_t *v, FILE *f) {
//	long int i;
//	for (i = 0; i < igraph_vector_size(v); i++) {
//		cout << VECTOR(*v)[i] << " ";
//		//fprintf(f, " %4.2f", );
//	}
//	cout << endl;
//	//fprintf(f, "\n");
//}

template<class T>
void printVector(vector<T> v) {

	cout << "vector's size = " << v.size() << endl;
	FOR (i, v.size())
	{
		cout << i << " ";
	}
	cout << endl;
	FOR (i, v.size())
	{
		cout << v[i] << " ";
	}
	cout << endl << endl;
}

vector<float> iVector_to_vector(igraph_vector_t *v) {
	int s = igraph_vector_size(v);
	vector<float> result(s);
	for (int i = 0; i < s; i++)
		result[i] = VECTOR(*v)[i];
	return result;
}

double P(int x, double c, double alpha) {
	double p = c * pow(x, -alpha);
	return p;
}

vector<float> getRatioByPowerLaw(int numberOfSamples, int V, int smin, int smax,
		float plalpha) {

//double alpha = 0.5;
//double alpha = 1;
//double alpha = 0.5;
	double alpha = plalpha;
	cout << "power law with alpha = " << alpha << endl;
	double sum = 0.0;
//int min_size = 40; int max_size = 1000;
	int min_size = smin;
	int max_size = smax;
	for (int x = min_size; x < max_size; x++) {
		sum = sum + pow(x, -alpha);
	}
	double c = 1 / sum;
	cout << "constant c = " << c << endl;

	vector<int> S;
	int count = 0;
	do {
		double Q_prev = 0.0;
		double Q_curr = 0.0;
		float y = (float) rand() / (float) RAND_MAX;
		for (int x = min_size; x < max_size; x++) {
			Q_curr = Q_prev + P(x, c, alpha);
			if (Q_prev < y && y <= Q_curr) {
				S.push_back(x);
				count++;
			}
			Q_prev = Q_curr;
		}
	} while (count < numberOfSamples);

//sort budget:
	std::sort(S.begin(), S.end(), std::less<int>());

	cout << "print budget: " << endl;
	int est_datapoints = 0;
	FOR(i, S.size())
	{
		cout << S[i] << " ";
		est_datapoints = est_datapoints + S[i];
	}
	cout << endl;

	cout << "estimated #datapoints: " << est_datapoints << endl;
	cout << "V: " << V << endl;
//convert to ratio
	vector<float> ratios;
	for (int i = 0; i < S.size(); i++) {
		float r = (float) S[i] / (float) V;
		ratios.push_back(r);
	}
//	printVector(ratios);

	return ratios;
}

void doSampling(string nfile, string sampling_method, int epoches, float iratio,
		bool dFlag, bool oFlag, bool pwFlag, int k_path, int c, float alpha,
		bool output_sample, bool vary_size, bool powLawFlag,
		int numberOfThreads, int min, int max, bool slowFlag, float plalpha,
		bool BCFlag, bool deepWalkFlag, bool CCFlag) {

	Graph orgGraph;
	orgGraph.loadGraphFromFile(nfile, dFlag);

	cout << "getting global node centrality... " << endl;
	orgGraph.getGlobalFeature(orgGraph.getEdges(), dFlag, BCFlag, CCFlag);

	cout << "init ratio: " << iratio << endl;
// do sampling
	int f_matrix_size = 0;
	int pairwise_matrix_size = 0;
	float ratio = iratio;

	vector<float> powLawRatio;
	if (powLawFlag) {
		powLawRatio = getRatioByPowerLaw(epoches, orgGraph.getNumberOfNodes(),
				min, max, plalpha);
		//std::sort(powLawRatio.begin(), powLawRatio.end(), std::less<float>()); //a->z
		//cout << "ratios: " << endl;
		//printVector(powLawRatio);
	}

	vector<Graph> graphPool;
	FOR(i, numberOfThreads)
	{
		Graph g = orgGraph;
		graphPool.push_back(g);
	}
	omp_set_num_threads(numberOfThreads);

	int total_steps = epoches;
	int step_completed = 0;

	cout << "total_steps: " << total_steps << endl;
	cout << "step_completed: " << step_completed << endl;

#pragma omp parallel for
	FOR(i, epoches)
	{
		int tid = omp_get_thread_num();
		cout << "thread id: " << tid << endl;

		Graph graph = graphPool[tid];

		cout << "epoch: " << i << endl;
		vector<int> YNodes, RNodes;
		vector<int> sampledEdgeList;

		if (powLawFlag) {

			ratio = powLawRatio[i];

			if (sampling_method == "BFS") {
				cout << "---------SAMPLING BY BFS---------" << endl;
				graph.doSamplingByBFS(ratio, YNodes, RNodes, sampledEdgeList);
			}
			if (sampling_method == "DFS") {
				cout << "---------SAMPLING BY DFS---------" << endl;
				graph.doSamplingByDFS(ratio, YNodes, RNodes, sampledEdgeList);
			}
			if (sampling_method == "RW") {
				cout << "---------SAMPLING BY RW---------" << endl;
				graph.doSamplingByRW(ratio, YNodes, RNodes, sampledEdgeList);
			}

		} else {

			if (vary_size) {
				int ratio_vary = rand() % 20;
				if (ratio_vary != 0)
					ratio = ratio_vary * iratio;
				cout << "ratio_vary: " << ratio_vary << endl;
				cout << "ratio: " << ratio << endl;
			}

			if (sampling_method == "BFS") {
				cout << "---------SAMPLING BY BFS---------" << endl;
				graph.doSamplingByBFS(ratio, YNodes, RNodes, sampledEdgeList);
			}
			if (sampling_method == "DFS") {
				cout << "---------SAMPLING BY DFS---------" << endl;
				graph.doSamplingByDFS(ratio, YNodes, RNodes, sampledEdgeList);
			}
			if (sampling_method == "RW") {
				cout << "---------SAMPLING BY RW---------" << endl;
				graph.doSamplingByRW(ratio, YNodes, RNodes, sampledEdgeList);
			}
		}

		if (output_sample == false) { // "-ws" flag
			// get features
			graph.generateFeatureMatrix(nfile, YNodes, RNodes, sampledEdgeList,
					dFlag, pwFlag, k_path, c, alpha, powLawFlag, slowFlag,
					BCFlag, deepWalkFlag, CCFlag);

			//write last sampled graph to files
			if (i == epoches - 1) {
				cout << "#Nodes in last sampled Gprime: "
						<< YNodes.size() + RNodes.size() << endl;
				cout << "#Edges in last sampled Gprime: "
						<< sampledEdgeList.size() / 2 << endl;
				cout << "#RNodes in last sampled Gprime: " << RNodes.size()
						<< endl;
				cout << "#YNodes in last sampled Gprime: " << YNodes.size()
						<< endl;

				string yfile = nfile + "_YNodes_"
						+ to_string(static_cast<long long>(i));
				string rfile = nfile + "_RNodes_"
						+ to_string(static_cast<long long>(i));
				string efile = nfile + "_Sampled_"
						+ to_string(static_cast<long long>(i));

				graph.write2File(yfile, YNodes);
				graph.write2File(rfile, RNodes);
				bool headerFlag = true;
				graph.write2File(efile, sampledEdgeList, RNodes, YNodes,
						headerFlag);
			}

			f_matrix_size += graph.getFeatureMatrix().size();
			pairwise_matrix_size += graph.getPairWiseFeatureMatrix().size();

#pragma omp critical
			{
				if (oFlag) { //"-o" flag write matrix?

					step_completed++;
					cout << "Progress: " << (step_completed * 100) / total_steps
							<< "%" << endl;
					//exclusively write matrix to files from each thread
					string ofile;
					if (pwFlag) {
						ofile = nfile + "_" + sampling_method
								+ "_pairwise_fmatrix";
						graph.write2File(graph.getPairWiseFeatureMatrix(),
								ofile); // append mode
						cout << "#rows in feature matrix (pairwise): "
								<< pairwise_matrix_size << endl;
						graph.clearfMatrix();
						graph.clearPairwiseMatrix();
					} else {
//						//write feature matrix to file

						string ofile;
						if (deepWalkFlag) {
							ofile = nfile + "_" + sampling_method + "_gain_"
									+ to_string(static_cast<long long>(i));
							graph.write2File(graph.getFeatureMatrix(), ofile); // append mode
							cout << "#rows in feature matrix: " << f_matrix_size
									<< endl;
							graph.clearfMatrix();

//							string yfile = nfile + "_YNodes_"
//									+ to_string(static_cast<long long>(i));
//							string rfile = nfile + "_RNodes_"
//									+ to_string(static_cast<long long>(i));
//							string efile = nfile + "_Sampled_"
//									+ to_string(static_cast<long long>(i));
//							graph.write2File(yfile, YNodes);
//							graph.write2File(rfile, RNodes);
//							bool headerFlag = false;
//							graph.write2File(efile, sampledEdgeList, RNodes,
//									YNodes, headerFlag);
						} else {
							ofile = nfile + "_" + sampling_method + "_fmatrix";

							graph.write2File(graph.getFeatureMatrix(), ofile); // append mode
							cout << "#rows in feature matrix: " << f_matrix_size
									<< endl;
							graph.clearfMatrix();
						}
					}
				}
			}
		} else if (output_sample) {
#pragma omp critical
			{
				//TODO: write samples to file
				cout << "#Nodes in sampled Gprime: "
						<< YNodes.size() + RNodes.size() << endl;
				cout << "#Edges in sampled Gprime: "
						<< sampledEdgeList.size() / 2 << endl;
				cout << "#RNodes in sampled Gprime: " << RNodes.size() << endl;
				cout << "#YNodes in sampled Gprime: " << YNodes.size() << endl;

				/*string yfile = nfile + "_YNodes_" + to_string(i);
				 string rfile = nfile + "_RNodes_" + to_string(i);
				 string efile = nfile + "_Sampled_" + to_string(i);*/
				string yfile = nfile + "_YNodes_"
						+ to_string(static_cast<long long>(i));
				string rfile = nfile + "_RNodes_"
						+ to_string(static_cast<long long>(i));
				string efile = nfile + "_Sampled_"
						+ to_string(static_cast<long long>(i));

				graph.write2File(yfile, YNodes);
				graph.write2File(rfile, RNodes);
				bool headerFlag = true;
				graph.write2File(efile, sampledEdgeList, RNodes, YNodes,
						headerFlag);
			}
		}
	}
}

map<int, double> loadGain(string fname) {
	map<int, double> mapping;

	std::ifstream in(fname);
	if (!in) {
		cout << "cannot load file " << fname << endl;
	}

	while (in.is_open()) {
		double u;
		double gain;
		while (in >> u >> gain) {
			int v = (int) u;
//			cout << v << "-" << gain << endl;
			mapping[v] = gain;
		}
		in.close();
	}

	return mapping;
}

void doProbing(string nfile, string rfile, string yfile, string efile,
		std::map<string, double> coefs_map, int probe_budget, int k, bool dFlag,
		bool logisticFlag, int k_path, int c, bool adapFlag,
		bool greedyProbeFlag, bool randomProbeFlag,
		bool randomForestRegressFlag, bool randomForestClassifyFlag,
		bool randomForestFeature, string randomForestModel, bool BCFlag,
		string lib_path, string mdl_path, string gainFromFile,
		bool deepWalkFlag, bool nodeCentralityFlag, bool CCFlag, bool inverse) {

	Graph g;
	g.loadGraphFromFile(nfile, dFlag); // build internal adjacency list

// TODO
//	string yfile = nfile + "_YNodes";
//	string rfile = nfile + "_RNodes";
//	string efile = nfile + "_Sampled";
	string sfile = nfile;

	vector<int> YNodes, RNodes;
	vector<int> sampledEdges, edges;

	g.readVectorFromFile(yfile, YNodes, false);
	g.readVectorFromFile(rfile, RNodes, false);
	g.readVectorFromFile(efile, sampledEdges, true);
	g.readVectorFromFile(sfile, edges, true);

	cout << "#Nodes in Gprime: " << YNodes.size() + RNodes.size() << endl;
	cout << "#Edges in Gprime: " << sampledEdges.size() / 2 << endl;
	cout << "#RNodes in Gprime: " << RNodes.size() << endl;
	cout << "#YNodes in Gprime: " << YNodes.size() << endl;

//	g.printVector(YNodes, false);
//	g.printVector(RNodes, false);
//	g.printEdgeList(sampledEdges);
//	g.printEdgeList(edges);

	cout << "read coefs_map:" << endl;
//	std::map<string, double>::iterator it;
//	for (it = coefs_map.begin(); it != coefs_map.end(); ++it) {
//		string key = it->first;
//		double coef = it->second;
//		cout << key << " " << coef << endl;
//	}

	bool predictHair = false;

	cout << "YNodes - before shuffling : " << endl;
	printVector(YNodes);
	for (int i = 0; i < 3; i++)
		g.shuffleVector(YNodes);
	cout << "YNodes - after shuffling: " << endl;
	printVector(YNodes);

	if (adapFlag == false) {
		if (greedyProbeFlag) {
			g.doGreedyNonAdaptiveProbing(edges, sampledEdges, RNodes, YNodes,
					probe_budget, k, dFlag);
		} else if (deepWalkFlag) {
			map<int, double> node2gain;
			node2gain = loadGain(gainFromFile);
			//		for (std::map<int, double>::const_iterator it = node2gain.begin();
			//				it != node2gain.end(); ++it) {
			//			cout << it->first << "-" << it->second << endl;
			//		}
			//		g.printMapping(node2gain);

			g.doNonAdaptiveProbingDeepWalk(edges, sampledEdges, RNodes, YNodes,
					probe_budget, k, node2gain, dFlag, BCFlag);
		} else if (randomProbeFlag) {
			g.doRandomNonAdaptiveProbing(edges, sampledEdges, RNodes, YNodes,
					probe_budget, k, dFlag);

		} else if (nodeCentralityFlag) {
//			g.doNonAdaptiveProbing3(edges, sampledEdges, RNodes, YNodes,
//					probe_budget, coefs_map, k, dFlag, BCFlag, CCFlag);
			g.doNonAdaptiveProbing2(edges, sampledEdges, RNodes, YNodes,
					probe_budget, k, dFlag, lib_path, mdl_path, BCFlag,
					false, CCFlag);

		} else {
			g.doNonAdaptiveProbing2(edges, sampledEdges, RNodes, YNodes,
					probe_budget, k, dFlag, lib_path, mdl_path, BCFlag,
					deepWalkFlag, CCFlag);
		}

	} else {
		if (predictHair)
			g.predictHairNode(edges, sampledEdges, RNodes, YNodes, lib_path,
					mdl_path, dFlag, BCFlag);
		else if (greedyProbeFlag)
			g.doGreedyAdaptiveProbing(edges, sampledEdges, RNodes, YNodes,
					probe_budget, k, dFlag);
		else if (randomProbeFlag)
			g.doRandomAdaptiveProbing(edges, sampledEdges, RNodes, YNodes,
					probe_budget, k, dFlag);
		else if (nodeCentralityFlag) {
			g.doAdaptiveProbing2(edges, sampledEdges, RNodes, YNodes,
					probe_budget, coefs_map, k, dFlag, BCFlag, CCFlag, inverse);
		} else if (adapFlag)
			g.doAdaptiveProbing(edges, sampledEdges, RNodes, YNodes,
					probe_budget, coefs_map, k, dFlag, logisticFlag, k_path, c,
					randomForestRegressFlag, randomForestClassifyFlag,
					randomForestFeature, randomForestModel, BCFlag, CCFlag);
		else
			g.doNonAdaptiveProbing(edges, sampledEdges, RNodes, YNodes,
					probe_budget, coefs_map, k, dFlag, logisticFlag, k_path, c,
					randomForestRegressFlag, randomForestClassifyFlag,
					randomForestFeature, randomForestModel, BCFlag, CCFlag);
	}
}

double get_wall_time() {
	struct timeval time;
	if (gettimeofday(&time, NULL)) {
		//  Handle error
		return 0;
	}
	return (double) time.tv_sec + (double) time.tv_usec * .000001;
}

int main(int argc, char ** argv) {

//	srand(5);

	OptionParser op(argc, argv);
	if (!op.validCheck()) {
		cout << "Parameters error!" << endl;
		return -1;
	}

	string action = op.getPara("-a");
	cout << action << endl;

//	if (action == "sampling") {

	string tmp;

	tmp = op.getPara("-f"); //graph to be sampled
	string nFile;
	if (tmp != "") {
		nFile = tmp;
	}

	tmp = op.getPara("-sm");
	string sampling_method;
	if (tmp != "") {
		sampling_method = tmp; // sampling method: BFS, DFS, RW
	}

	tmp = op.getPara("-varysize"); //vary sample size
	bool vary_size = false;
	if (tmp != "") {
		if (tmp == "TRUE")
			vary_size = true;
		else if (tmp == "FALSE")
			vary_size = false;
	}

	int epoches = 1;
	tmp = op.getPara("-epoches"); // #generated assampling data
	if (tmp != "")
		epoches = stoi(tmp);

	tmp = op.getPara("-ratio"); // ratio of sampled nodes/nodes
	float ratio;
	if (tmp != "")
		ratio = stof(tmp);

	tmp = op.getPara("-dFlag");
	bool dFlag;
	if (tmp != "") {
		if (tmp == "TRUE")
			dFlag = IGRAPH_DIRECTED;
		else if (tmp == "FALSE")
			dFlag = IGRAPH_UNDIRECTED;
	}

	tmp = op.getPara("-o"); // write feature matrix to file or not?
	bool oFlag;
	if (tmp != "") {
		if (tmp == "TRUE")
			oFlag = true;
		else if (tmp == "FALSE")
			oFlag = false;
	}

	tmp = op.getPara("-ws"); // write samples?
	bool sFlag;
	if (tmp != "") {
		if (tmp == "TRUE")
			sFlag = true;
		else if (tmp == "FALSE")
			sFlag = false;
	}

	tmp = op.getPara("-pw"); // run for logistic regression
	bool pairwiseFlag = false;
	if (tmp != "") {
		if (tmp == "TRUE")
			pairwiseFlag = true;
		else if (tmp == "FALSE")
			pairwiseFlag = false;
	}

	tmp = op.getPara("-powlaw"); // make budget obeys power law distribution
	bool powlawFlag;
	if (tmp != "") {
		if (tmp == "TRUE")
			powlawFlag = true;
		else if (tmp == "FALSE")
			powlawFlag = false;
	}

	tmp = op.getPara("-plalpha"); // make budget obeys power law distribution
	float plalpha;
	if (tmp != "")
		plalpha = stof(tmp);

	tmp = op.getPara("-slow");
	bool slowFlag;
	if (tmp != "") {
		if (tmp == "TRUE")
			slowFlag = true;
		else if (tmp == "FALSE")
			slowFlag = false;
	}

	tmp = op.getPara("-min"); //get min size of samples
	int min;
	if (tmp != "")
		min = stoi(tmp);

	tmp = op.getPara("-max"); //get max size of samples
	int max;
	if (tmp != "")
		max = stoi(tmp);

	tmp = op.getPara("-k_path");
	int k_path;
	if (tmp != "")
		k_path = stoi(tmp);

	tmp = op.getPara("-c");
	int c;
	if (tmp != "")
		c = stoi(tmp);

	tmp = op.getPara("-alpha");
	float alpha;
	if (tmp != "")
		alpha = stoi(tmp);

	tmp = op.getPara("-thread");
	int numberOfThreads;
	if (tmp != "") {
		numberOfThreads = stoi(tmp);
	}

	tmp = op.getPara("-BCFlag");
	bool BCFlag = true;
	if (tmp != "") {
		if (tmp == "TRUE")
			BCFlag = true;
		else if (tmp == "FALSE")
			BCFlag = false;

		cout << "BCFlag: " << BCFlag << endl;
	}

	tmp = op.getPara("-CCFlag");
	bool CCFlag = true;
	if (tmp != "") {
		if (tmp == "TRUE")
			CCFlag = true;
		else if (tmp == "FALSE")
			CCFlag = false;

		cout << "CCFlag: " << CCFlag << endl;
	}

	tmp = op.getPara("-DWFlag");
	bool DWFlag = false;
	if (tmp != "") {
		if (tmp == "TRUE")
			DWFlag = true;
		else if (tmp == "FALSE")
			DWFlag = false;
		cout << "DWFlag: " << DWFlag << endl;
	}

	tmp = op.getPara("-sRand"); //get max size of samples
	if (tmp != "") {
		srand(stoi(tmp));
	} else
		srand(5); //default seed rand

	cout << nFile << " " << sampling_method << " " << epoches << " " << ratio
			<< " " << dFlag << " " << oFlag << " " << pairwiseFlag
			<< " numberOfThreads: " << numberOfThreads << endl;

//	} else if (action == "probing") {

//	string tmp;

	tmp = op.getPara("-adap"); // write feature matrix to file or not?
	bool adapFlag = false;
	if (tmp != "") {
		if (tmp == "TRUE")
			adapFlag = true;
		else if (tmp == "FALSE")
			adapFlag = false;
	}

//	tmp = op.getPara("-f"); //path to original graph
//	string nFile;
//	if (tmp != "") {
//		nFile = tmp;
//		cout << nFile << endl;
//	}

	tmp = op.getPara("-rfile"); //path to red nodes of sample graph
	string rFile;
	if (tmp != "") {
		rFile = tmp;
		cout << rFile << endl;
	}

	tmp = op.getPara("-yfile"); //path to yellow nodes of sample graph
	string yFile;
	if (tmp != "") {
		yFile = tmp;
		cout << yFile << endl;
	}

	tmp = op.getPara("-efile"); //path to sample graph
	string eFile;
	if (tmp != "") {
		eFile = tmp;
		cout << eFile << endl;
	}

	tmp = op.getPara("-greedy"); // enable global greedy
	bool greedyProbeFlag = false;
	if (tmp != "") {
		if (tmp == "TRUE")
			greedyProbeFlag = true;
		else if (tmp == "FALSE")
			greedyProbeFlag = false;
	}

	tmp = op.getPara("-random"); // enable global greedy
	bool randomProbeFlag = false;
	if (tmp != "") {
		if (tmp == "TRUE")
			randomProbeFlag = true;
		else if (tmp == "FALSE")
			randomProbeFlag = false;
	}

	tmp = op.getPara("-randomForest"); // use random forest
	bool randomForestRegressFlag = false;
	bool randomForestClassifyFlag = false;
	if (tmp != "") {
		if (tmp == "regress")
			randomForestRegressFlag = true;
		else if (tmp == "classify")
			randomForestClassifyFlag = true;
	}

	tmp = op.getPara("-rfFeature"); // use local/global feature for random forest
	bool randomForestFeature = false;
	if (tmp != "") {
		if (tmp == "TRUE")
			randomForestFeature = true;
		else if (tmp == "FALSE")
			randomForestFeature = false;
	}

	tmp = op.getPara("-rfModel"); // path to random forest model
	string randomForestModel;
	if (tmp != "") {
		randomForestModel = tmp;
	}

	tmp = op.getPara("-budget"); //get probing budget
	int budget;
	if (tmp != "") {
		budget = stoi(tmp);
		cout << budget << endl;
	}

	tmp = op.getPara("-i"); //interval to print out result
	int ival;
	if (tmp != "") {
		ival = stoi(tmp);
		cout << ival << endl;
	}

//	tmp = op.getPara("-dFlag");
//	bool dFlag;
//	if (tmp != "") {
//		if (tmp == "TRUE")
//			dFlag = IGRAPH_DIRECTED;
//		else if (tmp == "FALSE")
//			dFlag = IGRAPH_UNDIRECTED;
//
//		cout << dFlag << endl;
//	}
//	tmp = op.getPara("-pw"); // write feature matrix to file or not?
//	bool pairwiseFlag = false;
//	if (tmp != "") {
//		if (tmp == "TRUE")
//			pairwiseFlag = true;
//		else if (tmp == "FALSE")
//			pairwiseFlag = false;
//		cout << pairwiseFlag << endl;
//	}
//
//	tmp = op.getPara("-k_path");
//	int k_path;
//	if (tmp != "")
//		k_path = stoi(tmp);
//
//	tmp = op.getPara("-c");
//	int c;
//	if (tmp != "")
//		c = stoi(tmp);

	tmp = op.getPara("-hairNodeFlag");
	bool hairFlag;
	if (tmp != "") {
		if (tmp == "TRUE")
			hairFlag = true;
		else if (tmp == "FALSE")
			hairFlag = false;

		cout << "hairFlag: " << hairFlag << endl;
	}

	tmp = op.getPara("-pythonLibPath");
	string pLibPath;
	if (tmp != "") {
		pLibPath = tmp;
		cout << "pLibPath: " << pLibPath << endl;
	}

	tmp = op.getPara("-pythonMdlPath");
	string pMdlPath;
	if (tmp != "") {
		pMdlPath = tmp;
		cout << "pythonMdlPath: " << pMdlPath << endl;
	}

	tmp = op.getPara("-gainFile");
	string gainFile;
	if (tmp != "") {
		gainFile = tmp;
		cout << "gainFile: " << gainFile << endl;
	}

	tmp = op.getPara("-nodeCentFlag");
	bool nodeCentFlag = false;
	if (tmp != "") {
		if (tmp == "TRUE")
			nodeCentFlag = true;
		else if (tmp == "FALSE")
			nodeCentFlag = false;
		cout << "nodeCentFlag: " << nodeCentFlag << endl;
	}

	tmp = op.getPara("-inverse");
	bool inverseProbing = false;
	if (tmp != "") {
		if (tmp == "TRUE")
			inverseProbing = true;
		else if (tmp == "FALSE")
			inverseProbing = false;
		cout << "inverseProbing: " << inverseProbing << endl;
	}

	std::map<string, double> coefs_map;
	/*
	if (nodeCentFlag) {

		string x0, x1, x2, x3, x4, x5, x6, x7;
		string x8, x9, x10, x11, x12, x13, x14;
		string x15, x16, x17, x18, x19, x20, x21, x22, x23;
		string x24, x25;

		//linear regression
		x0 = op.getPara("-x0");
		x1 = op.getPara("-x1");
		x2 = op.getPara("-x2");
		x3 = op.getPara("-x3");
		x4 = op.getPara("-x4");
		x5 = op.getPara("-x5");
		x6 = op.getPara("-x6");
		x7 = op.getPara("-x7");
		x8 = op.getPara("-x8");
		x9 = op.getPara("-x9");
		x10 = op.getPara("-x10");
		x11 = op.getPara("-x11");
		x12 = op.getPara("-x12");
		x13 = op.getPara("-x13");
		x14 = op.getPara("-x14");
		x15 = op.getPara("-x15");
		x16 = op.getPara("-x16");
		x17 = op.getPara("-x17");
		x18 = op.getPara("-x18");
		x19 = op.getPara("-x19");
		x20 = op.getPara("-x20");
		x21 = op.getPara("-x21");
		x22 = op.getPara("-x22");
		x23 = op.getPara("-x23");

		//TODO to be removed
		//x24 = op.getPara("-x24");

		coefs_map[Graph::A0] = stod(x0);		//intercept

		coefs_map[Graph::A1] = stod(x1);
		coefs_map[Graph::A2] = stod(x2);
		coefs_map[Graph::A3] = stod(x3);
		coefs_map[Graph::A4] = stod(x4);
		coefs_map[Graph::A5] = stod(x5);		//m
		coefs_map[Graph::A6] = stod(x6);
		coefs_map[Graph::A7] = stod(x7);		//not used as feature
		coefs_map[Graph::A8] = stod(x8);		//m_prime
		coefs_map[Graph::A9] = stod(x9);		//not used as feature

		coefs_map[Graph::A10] = stod(x10);		//u_2_ynodes (estimated)
		coefs_map[Graph::A11] = stod(x11);
		coefs_map[Graph::A12] = stod(x12);
		coefs_map[Graph::A13] = stod(x13);
		coefs_map[Graph::A14] = stod(x14);
		coefs_map[Graph::A15] = stod(x15);
		coefs_map[Graph::A16] = stod(x16);
		coefs_map[Graph::A17] = stod(x17);		//not used as feature
		coefs_map[Graph::A18] = stod(x18);
		coefs_map[Graph::A19] = stod(x19);
		coefs_map[Graph::A20] = stod(x20);
		coefs_map[Graph::A21] = stod(x21);
		coefs_map[Graph::A22] = stod(x22);
		coefs_map[Graph::A23] = stod(x23);

		//TODO to be removed
		//coefs_map[Graph::A24] = stod(x24);
		//cout << "x24: " << x24 << endl;

		if (pairwiseFlag) {
			string x24, x25, x26, x27, x28, x29;
			string x30, x31, x32, x33, x34, x35, x36, x37, x38;
			string x39, x40, x41, x42, x43, x44, x45, x46;

			//logistic regression:
			x24 = op.getPara("-x24");
			x25 = op.getPara("-x25");
			x26 = op.getPara("-x26");
			x27 = op.getPara("-x27");
			x28 = op.getPara("-x28");
			x29 = op.getPara("-x29");
			x30 = op.getPara("-x30");
			x31 = op.getPara("-x31");
			x32 = op.getPara("-x32");
			x33 = op.getPara("-x33");
			x34 = op.getPara("-x34");
			x35 = op.getPara("-x35");
			x36 = op.getPara("-x36");
			x37 = op.getPara("-x37");
			x38 = op.getPara("-x38");
			x39 = op.getPara("-x39");
			x40 = op.getPara("-x40");
			x41 = op.getPara("-x41");
			x42 = op.getPara("-x42");
			x43 = op.getPara("-x43");
			x44 = op.getPara("-x44");
			x45 = op.getPara("-x45");
			x46 = op.getPara("-x46");

			coefs_map[Graph::A1_v] = stod(x24);
			coefs_map[Graph::A2_v] = stod(x25);
			coefs_map[Graph::A3_v] = stod(x26);
			coefs_map[Graph::A4_v] = stod(x27);
			coefs_map[Graph::A5_v] = stod(x28);		//m
			coefs_map[Graph::A6_v] = stod(x29);
			coefs_map[Graph::A7_v] = stod(x30);		//not used as feature
			coefs_map[Graph::A8_v] = stod(x31);		//m_prime
			coefs_map[Graph::A9_v] = stod(x32);		//not used as feature
			coefs_map[Graph::A10_v] = stod(x33);	//u_2_ynodes (estimated)

			coefs_map[Graph::A11_v] = stod(x34);
			coefs_map[Graph::A12_v] = stod(x35);
			coefs_map[Graph::A13_v] = stod(x36);
			coefs_map[Graph::A14_v] = stod(x37);
			coefs_map[Graph::A15_v] = stod(x38);
			coefs_map[Graph::A16_v] = stod(x39);
			coefs_map[Graph::A17_v] = stod(x40);	//not used as feature
			coefs_map[Graph::A18_v] = stod(x41);
			coefs_map[Graph::A19_v] = stod(x42);
			coefs_map[Graph::A20_v] = stod(x43);
			coefs_map[Graph::A21_v] = stod(x44);
			coefs_map[Graph::A22_v] = stod(x45);
			coefs_map[Graph::A23_v] = stod(x46);
		}
	}
	*/

	double wall1 = get_wall_time();

	if (action == "sampling")
		doSampling(nFile, sampling_method, epoches, ratio, dFlag, oFlag,
				pairwiseFlag, k_path, c, alpha, sFlag, vary_size, powlawFlag,
				numberOfThreads, min, max, slowFlag, plalpha, BCFlag, DWFlag,
				CCFlag);
	else if ((action == "probing")) {
		doProbing(nFile, rFile, yFile, eFile, coefs_map, budget, ival, dFlag,
				pairwiseFlag, k_path, c, adapFlag, greedyProbeFlag,
				randomProbeFlag, randomForestRegressFlag,
				randomForestClassifyFlag, randomForestFeature,
				randomForestModel, BCFlag, pLibPath, pMdlPath, gainFile, DWFlag,
				nodeCentFlag, CCFlag, inverseProbing);
	}
	double wall2 = get_wall_time();

	cout << "elapsed time for probing (mins): " << (wall2 - wall1) / 60 << endl;

//	} else
//		cout << "invalid parameters" << endl;

}
