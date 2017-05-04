/*
 * IGraph.cpp
 *
 *  Created on: May 31, 2016
 *      Author: trinpm
 */

#include "IGraph.h"

void IGraph::loadGraphFromEdgeList(vector<int> edgeList, iBool dFlag) {

	vector<iReal> edges;
	FOR(i, edgeList.size())
	{
		edges.push_back(edgeList[i]);
//		if (i % 2 == 0)
//			cout << edgeList[i];
//		else
//			cout << " " << edgeList[i] << endl;

	}

	cout << "dFlag: " << dFlag << endl;

	iVector v;
	igraph_vector_view(&v, (iReal*) (&edges[0]), edges.size());
	igraph_create(&g, &v, 0, dFlag ? IGRAPH_DIRECTED : IGRAPH_UNDIRECTED);

	//get #nodes
	n = igraph_vcount(&g);
	cout << "igraph-#nodes: " << n << endl;

	//get #edges
	m = igraph_ecount(&g);
	cout << "igraph-#edges: " << m << endl;

	cout << "graph loaded by igraph" << endl;

}

void IGraph::iterateVertices() {

	igraph_vs_t vs;
	igraph_vit_t vit;

	igraph_integer_t size;

	igraph_vs_seq(&vs, 0, n);

	igraph_vit_create(&g, vs, &vit);

	igraph_vs_size(&g, &vs, &size);

//	printf("%li", (long int) size);
	cout << "vertex sequence's size: " << size << endl;

	while (!IGRAPH_VIT_END(vit)) {
		cout << (long int) IGRAPH_VIT_GET(vit) << " ";
		IGRAPH_VIT_NEXT(vit);

	}
	cout << endl;

	igraph_vit_destroy(&vit);
	igraph_vs_destroy(&vs);

}

vector<int> IGraph::getNeighbors(int u, iBool dFlag) {

	iVector ivec;
	igraph_vector_init(&ivec, n);

	if (dFlag)
		igraph_neighbors(&g, &ivec, u, IGRAPH_OUT);
	else
		igraph_neighbors(&g, &ivec, u, IGRAPH_ALL);

	igraph_vector_sort(&ivec);

	vector<int> v;

	v = iVector_to_vector2(&ivec);

//	printVector(v);

	igraph_vector_destroy(&ivec);

	return v;
}

vector<float> IGraph::computeNodeCentrality(int nodeCentrality,
iBool dFlag) {

	iReal cent;
	iVector ires;
	vector<float> res;

	igraph_vector_init(&ires, 0);

	switch (nodeCentrality) {

	case IN_DEGREE:
		igraph_centralization_degree(&g, /*res=*/&ires,
		/*mode=*/IGRAPH_IN, IGRAPH_NO_LOOPS, /*value*/&cent, /*theoretical_max=*/
		0,
		/*normalized=*/0);
		break;

	case OUT_DEGREE:
		igraph_centralization_degree(&g, /*res=*/&ires,
		/*mode=*/IGRAPH_OUT, IGRAPH_NO_LOOPS, /*value*/&cent, /*theoretical_max=*/
		0,
		/*normalized=*/0);
		break;

	case BETWEENNESS_CENT:

		igraph_centralization_betweenness(&g, /*res=*/&ires, dFlag, /*nobigint=*/
		0, /*value*/&cent,
		/*theoretical_max=*/0,
		/*normalized=*/0);

//		igraph_betweenness_estimate(
//		/*graph*/&g,
//		/*res=*/&ires,
//		/*vids=*/ igraph_vss_all(),
//		/*directed=*/ dFlag,
//		/*cutoff=*/ 3,
//		/*weights=*/ 0,
//		/*nobigint=*/ 1);

		break;

	case CLOSENESS_CENT:
		//normalized version has different value with one calculated in R
		igraph_closeness(&g, &ires, igraph_vss_all(), IGRAPH_OUT, /*weights*/0, /*normalization=*/
		0);
//		res = iVector_to_vector(&ires);
		break;

	case PAGERANK_CENT:
		igraph_pagerank(&g, IGRAPH_PAGERANK_ALGO_PRPACK, &ires, /*value*/0,
				igraph_vss_all(), dFlag, 0.85, /*weights*/0, /*options*/0);
		break;

	case EIGEN_CENT:
		igraph_arpack_options_t arpack_options;
		igraph_arpack_options_init(&arpack_options);
		igraph_centralization_eigenvector_centrality(&g, /*vector=*/&ires,
		/*value=*/0,
		/*directed=*/dFlag,
		/*scale=*/0, &arpack_options, /*value*/&cent,
		/*theoretical_max=*/0,
		/*normalization=*/0);
		break;
	}

	res = iVector_to_vector(&ires);

	igraph_vector_destroy(&ires);

	return res;
}

vector<float> IGraph::computeClusCoef(iBool dflag) {

	iVector ires;
	igraph_vector_init(&ires, 0);
	igraph_transitivity_local_undirected(&g, &ires, igraph_vss_all(),
			IGRAPH_TRANSITIVITY_ZERO);

	vector<float> res;
	res = iVector_to_vector(&ires);

	igraph_vector_destroy(&ires);

	return res;
}

template<class T>
void IGraph::printVector(vector<T> v) {

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
	cout << endl;
}

vector<float> IGraph::iVector_to_vector(igraph_vector_t *v) {
	int s = igraph_vector_size(v);
	vector<float> result(s);
	for (int i = 0; i < s; i++)
		result[i] = VECTOR(*v)[i];
	return result;
}

vector<int> IGraph::iVector_to_vector2(igraph_vector_t *v) {
	int s = igraph_vector_size(v);
	vector<int> result(s);
	for (int i = 0; i < s; i++)
		result[i] = VECTOR(*v)[i];
	return result;
}

