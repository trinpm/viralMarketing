#ifndef _GRAPHLBL_H_
#define _GRAPHLBL_H_

#include <iostream>
#include <fstream>
#include <queue>
#include <vector>
#include <map>

typedef int d_type;

using namespace std;

class GraphLabel{
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
		vector<d_type> distance;
		vector<d_type> predecessor;
		vector<d_type> redNodes;
		vector<d_type> yellowNodes;
		int k;

	public:
		GraphLabel(char * graphfile, vector<d_type> & probedNodes, int budget);
		void setRedNodes(vector<d_type> & red);
		void bfs();
		d_type getBestNode(int & budget);
		d_type getBestNode_2(int & budget);
		void probeNode(d_type u);
		d_type getTotalBenefit();

		map<d_type, int> benefitMapping();
};

#endif
