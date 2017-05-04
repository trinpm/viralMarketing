/*
 * IGraph.h
 *
 *  Created on: May 31, 2016
 *      Author: trinpm
 */

#ifndef IGRAPH_H_
#define IGRAPH_H_

#define	iGraph		igraph_t
#define	iVector		igraph_vector_t
#define	iInteger	igraph_integer_t
#define iBool		igraph_bool_t
#define	iReal		igraph_real_t
#define	iMatrix		igraph_matrix_t
#define iSelector 	igraph_vs_t

#define	FOR(i, n) for(int i = 0; i < (int)(n); ++i)

#include <igraph.h>
#include <string>
#include <iostream>
#include <vector>

using namespace std;

typedef enum {
	IN_DEGREE = 1,
	OUT_DEGREE = 2,
	BETWEENNESS_CENT = 3,
	CLOSENESS_CENT = 4,
	PAGERANK_CENT = 5,
	EIGEN_CENT = 6
} property;

class IGraph {

//	friend class Graph; // allow Graph access private/protected members of IGraph

	iGraph g;

	int n; //#nodes
	int m; //#edges

public:

	IGraph() {
		n = 0;
		m = 0;
	}
	;

	IGraph(IGraph &obj){
		//Here, obj is a reference to an object that is being used to initialize another object.
	}; // copy constructor - used as copying an object to pass it as an argument to a function

	int get_n() {
		return n;
	}

	void iterateVertices();

	vector<int> getNeighbors(int u, iBool dFlag);

	/**
	 * load graph from edgeList
	 * edgeList format: edges[] = {u1, v1, u2, v2,..., un, vn};
	 */
	void loadGraphFromEdgeList(vector<int> edgeList, iBool directed);

	/**
	 * @brief get centrality of nodes in graph
	 */
	vector<float> computeNodeCentrality(int nodeCentrality,
	iBool dflag);

	/**
	 * @briefg compute clustering coefficient
	 */
	vector<float> computeClusCoef(iBool dflag);

	/**
	 * @brief print out STL vector
	 */
	template<class T>
	void printVector(vector<T> v);

	/**
	 * @brief convert igraphVector to STL vector
	 */
	vector<float> iVector_to_vector(igraph_vector_t *v);

	vector<int> iVector_to_vector2(igraph_vector_t *v);

	~IGraph() {
		igraph_destroy(&g);
	}
	;
};

#endif /* IGRAPH_H_ */
