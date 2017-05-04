/*
 * GraphLabelSlow.cpp
 *
 *  Created on: Sep 10, 2016
 *      Author: trinpm
 */

#include "GraphLabelSlow.h"

/*
 * real the graph and partially observed subgraph
 */
GraphLabelSlow::GraphLabelSlow(char * graphfile, vector<d_type> & red,
		int budget) {

	ifstream graphf(graphfile);
	d_type tmp1, tmp2;

	graphf >> numNodes >> numEdges;
	edgelist = vector<vector<d_type> >(numNodes + 1);
	nodeDegrees = vector<d_type>(numNodes + 1, 0);

	for (int i = 0; i < numEdges; ++i) {
		graphf >> tmp1 >> tmp2;
		edgelist[tmp1].push_back(tmp2);
		edgelist[tmp2].push_back(tmp1);
		nodeDegrees[tmp1]++;
		nodeDegrees[tmp2]++;
	}
	graphf.close();

	redNodes = red;

	k = budget;

	observedNodes = vector<bool>(numNodes + 1, false);

	/*for (d_type u : red) {
	 selectedNodes.push_back(u);
	 observedNodes[u] = true;
	 }*/

	for (int idx = 0; idx < red.size(); idx++) {
		int u = red[idx];
		selectedNodes.push_back(u);
		observedNodes[u] = true;
	}

	numSelected = selectedNodes.size();
	for (int i = 0; i < numSelected; ++i) {
		d_type tmpnode1 = selectedNodes[i];
		for (int j = 0; j < nodeDegrees[tmpnode1]; ++j) {
			d_type tmpnode2 = edgelist[tmpnode1][j];
			if (!observedNodes[tmpnode2]) {
				yellowNodes.push_back(tmpnode2);
				observedNodes[tmpnode2] = true;
			}
		}
	}

}

void GraphLabelSlow::setRedNodes(vector<d_type> & red) {
	selectedNodes.clear();
	observedNodes = vector<bool>(numNodes + 1, false);
	marginalDegrees = nodeDegrees;
	markvisited = vector<vector<bool> >(numNodes + 1,
			vector<bool>(numNodes + 1, false));
	markobserved = vector<vector<bool> >(numNodes + 1,
			vector<bool>(numNodes + 1, false));
	accumulatedBenefit = vector<d_type>(numNodes + 1, 0);
	distance = vector<d_type>(numNodes + 1, numNodes + 1);
	predecessor = vector<d_type>(numNodes + 1, numNodes + 1);

	/*for (d_type u : red) {
	 selectedNodes.push_back(u);
	 observedNodes[u] = true;
	 }*/

	for (int idx = 0; idx < red.size(); idx++) {
		int u = red[idx];
		selectedNodes.push_back(u);
		observedNodes[u] = true;
	}

	numSelected = selectedNodes.size();
	for (int i = 0; i < numSelected; ++i) {
		d_type tmpnode1 = selectedNodes[i];
		for (int j = 0; j < nodeDegrees[tmpnode1]; ++j) {
			d_type tmpnode2 = edgelist[tmpnode1][j];
			marginalDegrees[tmpnode2]--;
			if (!observedNodes[tmpnode2]) {
				observedNodes[tmpnode2] = true;
				for (int k = 0; k < nodeDegrees[tmpnode2]; ++k) {
					marginalDegrees[edgelist[tmpnode2][k]]--;
				}
			}
		}
	}

	numObserved = 0;

	//for (d_type v : selectedNodes) {
	for (int i = 0; i < selectedNodes.size(); i++) {
		int v = selectedNodes[i];
		markvisited[selectedNodes[0]][v] = true;
		markobserved[selectedNodes[0]][v] = true;
		/*for (d_type u : edgelist[v]) {
		 cout << "u1: " << u << endl;
		 markobserved[selectedNodes[0]][u] = true;
		 }*/
		for (int j = 0; j < edgelist[v].size(); j++) {
			int u = edgelist[v][j];
			markobserved[selectedNodes[0]][u] = true;
		}
	}
	predecessor[selectedNodes[0]] = -1;
	distance[selectedNodes[0]] = 0;

	for (int i = 1; i < selectedNodes.size(); ++i) {
		markvisited[selectedNodes[i]] = markvisited[selectedNodes[0]];
		markobserved[selectedNodes[i]] = markobserved[selectedNodes[0]];
		predecessor[selectedNodes[i]] = -1;
		distance[selectedNodes[i]] = 0;
	}

}

/*
 * Run BFS to set the distances and predecessors for nodes
 */
void GraphLabelSlow::bfs(int budget) {
	queue<d_type> nextvisit;

	for (int i = 0; i < numSelected; ++i) {
		nextvisit.push(selectedNodes[i]);
	}

	d_type tmpnode, tmpnode1;
	d_type benefit0 = 0;
	d_type tmpBenefit;

	while (!nextvisit.empty()) {
		tmpnode = nextvisit.front();
		nextvisit.pop();
		if (distance[tmpnode] > budget)
			continue;

		benefit0 = accumulatedBenefit[tmpnode];

		//for (d_type node : edgelist[tmpnode]) {
		for (int i = 0; i < edgelist[tmpnode].size(); i++) {
			int node = edgelist[tmpnode][i];

			tmpnode1 = tmpnode;
			while (distance[tmpnode1] > 0 && predecessor[tmpnode1] != node) {
				tmpnode1 = predecessor[tmpnode1];
			}

			if (!markvisited[tmpnode][node] && distance[tmpnode1] == 0) {
				tmpBenefit = benefit0;
				//for (d_type node1 : edgelist[node]) {
				for (int j = 0; j < edgelist[node].size(); j++) {
					int node1 = edgelist[node][j];
					if (!markobserved[tmpnode][node1]) {
						tmpBenefit += 1;
					}
				}

				if ((double) tmpBenefit / (double) (distance[tmpnode] + 1)
						> (double) accumulatedBenefit[node]
								/ (double) distance[node]) {
					accumulatedBenefit[node] = tmpBenefit;
					distance[node] = distance[tmpnode] + 1;
					predecessor[node] = tmpnode;
					nextvisit.push(node);
					markvisited[node] = markvisited[tmpnode];
					markvisited[node][node] = true;
					markobserved[node] = markobserved[tmpnode];
					//for (d_type node1 : edgelist[node]) {
					for (int idx = 0; idx < edgelist[node].size(); idx++) {
						int node1 = edgelist[node][idx];
						markobserved[node][node1] = true;
					}
				}
			}

		}
	}
}

/*
 * Return the node with highest ratio of marginal benefit to distance from observed nodes
 */
d_type GraphLabelSlow::getBestNode(int & budget) {
	d_type bestNode = -1;
	double bestValue = -1;
	for (int i = 0; i <= numNodes; ++i) {
		if (distance[i] != 0 && distance[i] <= budget
				&& (double) accumulatedBenefit[i] / (double) distance[i]
						> bestValue) {
			bestValue = (double) accumulatedBenefit[i] / (double) distance[i];
			bestNode = i;
		}
	}

	return bestNode;
}

/*
 * Probe node u and update the distances, predecessor and other related variables
 */
void GraphLabelSlow::probeNode(d_type u, int & k) {
	if (distance[u] == 0)
		return;

	vector<d_type> updates;
	d_type tmpnode = u;
	k -= distance[tmpnode];
	while (distance[tmpnode] > 0) {
		accumulatedBenefit[tmpnode] = 0;
		updates.push_back(tmpnode);
		numSelected++;
		selectedNodes.push_back(tmpnode);
		marginalDegrees[tmpnode] = 0;
		tmpnode = predecessor[tmpnode];
	}

	//for (d_type tmp : updates) {
	for (int idx = 0; idx < updates.size(); idx++) {
		int tmp = updates[idx];
		if (!observedNodes[tmp]) {
			observedNodes[tmp] = true;
			numObserved++;
		}
		for (int i = 0; i < nodeDegrees[tmp]; ++i) {
			d_type tmp1 = edgelist[tmp][i];
			if (!observedNodes[tmp1]) {
				observedNodes[edgelist[tmp][i]] = true;
				numObserved++;
			}
		}
	}

	markvisited = vector<vector<bool> >(numNodes + 1,
			vector<bool>(numNodes + 1, false));
	markobserved = vector<vector<bool> >(numNodes + 1,
			vector<bool>(numNodes + 1, false));
	accumulatedBenefit = vector<d_type>(numNodes + 1, 0);
	distance = vector<d_type>(numNodes + 1, numNodes + 1);
	predecessor = vector<d_type>(numNodes + 1, numNodes + 1);

	//for (d_type v : selectedNodes) {
	for (int i = 0; i < selectedNodes.size(); i++) {
		int v = selectedNodes[i];
		markvisited[selectedNodes[0]][v] = true;
		markobserved[selectedNodes[0]][v] = true;
		/*for (d_type u : edgelist[v]) {
		 markobserved[selectedNodes[0]][u] = true;
		 }*/
		for (int j = 0; j < edgelist[v].size(); j++) {
			int u = edgelist[v][j];
			markobserved[selectedNodes[0]][u] = true;
		}
	}
	predecessor[selectedNodes[0]] = -1;
	distance[selectedNodes[0]] = 0;

	for (int i = 1; i < selectedNodes.size(); ++i) {
		markvisited[selectedNodes[i]] = markvisited[selectedNodes[0]];
		markobserved[selectedNodes[i]] = markobserved[selectedNodes[0]];
		predecessor[selectedNodes[i]] = -1;
		distance[selectedNodes[i]] = 0;
	}
}

d_type GraphLabelSlow::getTotalBenefit() {
	return numObserved;
}

map<d_type, int> GraphLabelSlow::benefitMapping() {
	map<d_type, int> m;
	vector<d_type> tmp;
	d_type vmax;
	int t;
	setRedNodes(redNodes);
	vector<int> yBenefit(yellowNodes.size(), 0);
	for (int i = 0; i < yellowNodes.size(); ++i) {
		yBenefit[i] = marginalDegrees[yellowNodes[i]];
	}
	int marginalBen = 0;
	for (int i = 0; i < yellowNodes.size(); ++i) {
		tmp = redNodes;
		tmp.push_back(yellowNodes[i]);
		setRedNodes(tmp);
		t = k - 1;
		marginalBen = -1;
		while (t > 0 && marginalBen < getTotalBenefit()) {
			marginalBen = getTotalBenefit();
			bfs(t);
			vmax = getBestNode(t);
			probeNode(vmax, t);
		}
		m.insert(
				pair<d_type, int>(yellowNodes[i],
						yBenefit[i] + getTotalBenefit()));
	}
	return m;
}
