/*
 * GraphLabelSampling.h
 *
 *  Created on: Sep 10, 2016
 *      Author: trinpm
 */

#ifndef GRAPHLABELSLOW_H_
#define GRAPHLABELSLOW_H_

#include <iostream>
#include <fstream>
#include <queue>
#include <stack>
#include <vector>
#include <map>

typedef int d_type;

using namespace std;

class GraphLabelSlow{
	private:
		vector<vector<d_type> > edgelist;
		vector<d_type> nodeDegrees;
		vector<d_type> selectedNodes;
		d_type numNodes;
		d_type numEdges;
		d_type numSelected;
		d_type numObserved;
		vector<bool> observedNodes;
		vector<d_type> marginalDegrees;
		vector<vector<bool> > markvisited;
		vector<vector<bool> > markobserved;
		vector<d_type> accumulatedBenefit;
		vector<d_type> distance;
		vector<d_type> predecessor;
		vector<d_type> redNodes;
                vector<d_type> yellowNodes;
                int k;

	public:
        GraphLabelSlow(char * graphfile, vector<d_type> & probedNodes, int budget);
		void setRedNodes(vector<d_type> & red);
		void bfs(int budget);
		d_type getBestNode(int & budget);
		void probeNode(d_type u, int & k);
		d_type getTotalBenefit();
                map<d_type, int> benefitMapping();
};


#endif /* GRAPHLABELSLOW_H_ */
