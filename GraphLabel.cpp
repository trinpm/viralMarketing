#include "GraphLabel.h"

/*
 * real the GraphLabel and partially observed subGraphLabel
 */
GraphLabel::GraphLabel(char * GraphLabelfile, vector<d_type> & red,
		int budget) {

	string str(GraphLabelfile);
//	cout << "GraphLabelfile: " << str << endl;

	ifstream GraphLabelf(GraphLabelfile);
	d_type tmp1, tmp2;

	GraphLabelf >> numNodes >> numEdges;
	edgelist = vector<vector<d_type> >(numNodes + 1);
	nodeDegrees = vector<d_type>(numNodes + 1, 0);

	for (int i = 0; i < numEdges; ++i) {
		GraphLabelf >> tmp1 >> tmp2;
		edgelist[tmp1].push_back(tmp2);
		edgelist[tmp2].push_back(tmp1);
		nodeDegrees[tmp1]++;
		nodeDegrees[tmp2]++;
	}
	GraphLabelf.close();

	redNodes = red;

	k = budget;

	observedNodes = vector<bool>(numNodes + 1, false);

	/*for (d_type u : red) {
	 selectedNodes.push_back(u);
	 observedNodes[u] = true;
	 }*/

	for (int i = 0; i < red.size(); i++) {
		int u = red[i];
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

/*
 * Set sed nodes
 */
void GraphLabel::setRedNodes(vector<d_type> & red) {
	selectedNodes.clear();
	observedNodes = vector<bool>(numNodes + 1, false);
	marginalDegrees = nodeDegrees;
	distance = vector<d_type>(numNodes + 1, numNodes + 1);
	predecessor = vector<d_type>(numNodes + 1, numNodes + 1);

	/*for (d_type u : red) {
	 selectedNodes.push_back(u);
	 observedNodes[u] = true;
	 }*/

	for (int i = 0; i < red.size(); i++) {
		int u = red[i];
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
}

/*
 * Run BFS to set the distances and predecessors for nodes
 */
void GraphLabel::bfs() {
	vector<bool> visited(numNodes + 1, false);
	queue<d_type> nextvisit;

	for (int i = 0; i < numSelected; ++i) {
		nextvisit.push(selectedNodes[i]);
		distance[selectedNodes[i]] = 0;
		visited[selectedNodes[i]] = true;
	}

	d_type tmpnode;
	while (!nextvisit.empty()) {
		tmpnode = nextvisit.front();
		nextvisit.pop();
		for (int i = 0; i < nodeDegrees[tmpnode]; ++i) {
			if (!visited[edgelist[tmpnode][i]]) {
				distance[edgelist[tmpnode][i]] = distance[tmpnode] + 1;
				predecessor[edgelist[tmpnode][i]] = tmpnode;
				nextvisit.push(edgelist[tmpnode][i]);
				visited[edgelist[tmpnode][i]] = true;
			}
		}
	}
}

/*
 * Return the node with highest ratio of marginal benefit to distance from observed nodes
 */
d_type GraphLabel::getBestNode(int & budget) {
	d_type bestNode = -1;
	double bestValue = -1;
	for (int i = 0; i <= numNodes; ++i) {
		if (distance[i] != 0 && distance[i] <= budget
				&& (double) marginalDegrees[i] / (double) distance[i]
						> bestValue) {
			bestValue = (double) marginalDegrees[i] / (double) distance[i];
			bestNode = i;
		}
	}

	budget -= distance[bestNode];

	return bestNode;
}

/*
 * Second version of the previous function
 */
d_type GraphLabel::getBestNode_2(int & budget) {
	d_type bestNode = -1;
	double bestValue = -1;

	int tmp_bestValue = -1;
	d_type tmp_node;
	for (int i = 0; i < numNodes + 1; ++i) {
		if (distance[i] > 0 && distance[i] <= budget) {
			vector<bool> tmp_observedNodes = observedNodes;
			tmp_bestValue = 0;
			tmp_node = i;
			while (distance[tmp_node] > 0) {
				for (d_type u : edgelist[tmp_node]) {
					if (!tmp_observedNodes[u]) {
						tmp_bestValue++;
						tmp_observedNodes[u] = true;
					}
				}
				tmp_node = predecessor[tmp_node];
			}
			if ((double) tmp_bestValue / (double) distance[i] > bestValue) {
				bestValue = (double) tmp_bestValue / (double) distance[i];
				bestNode = i;
			}
		}
	}

	budget -= distance[bestNode];

	return bestNode;
}

/*
 * Probe node u and update the distances, predecessor and other related variables
 */
void GraphLabel::probeNode(d_type u) {
	if (distance[u] == 0)
		return;

	queue<d_type> updateNodes;
	d_type tmpnode = u;
	while (distance[tmpnode] != 0) {
		updateNodes.push(tmpnode);
		numSelected++;
		selectedNodes.push_back(tmpnode);
		distance[tmpnode] = 0;
		if (!observedNodes[tmpnode]) {
			observedNodes[tmpnode] = true;
			numObserved++;
			for (int i = 0; i < nodeDegrees[tmpnode]; ++i) {
				if (distance[edgelist[tmpnode][i]] > 0) {
					marginalDegrees[edgelist[tmpnode][i]]--;
				}
			}
		}
		for (int i = 0; i < nodeDegrees[tmpnode]; ++i) {
			d_type tmpnode1 = edgelist[tmpnode][i];
			if (!observedNodes[tmpnode1]) {
				observedNodes[edgelist[tmpnode][i]] = true;
				numObserved++;
				for (int j = 0; j < nodeDegrees[tmpnode1]; ++j) {
					if (distance[edgelist[tmpnode1][j]] > 0) {
						marginalDegrees[edgelist[tmpnode1][j]]--;
					}
				}
			}
		}
		tmpnode = predecessor[tmpnode];
	}

	while (!updateNodes.empty()) {
		tmpnode = updateNodes.front();
		updateNodes.pop();
		for (int i = 0; i < nodeDegrees[tmpnode]; ++i) {
			if (distance[edgelist[tmpnode][i]] > distance[tmpnode] + 1) {
				distance[edgelist[tmpnode][i]] = distance[tmpnode] + 1;
				predecessor[edgelist[tmpnode][i]] = tmpnode;
				updateNodes.push(edgelist[tmpnode][i]);
			}
		}
	}
}

d_type GraphLabel::getTotalBenefit() {
	return numObserved;
}

map<d_type, int> GraphLabel::benefitMapping() {
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
		bfs();
		t = k - 1;
		marginalBen = -1;
		while (t > 0 && marginalBen < getTotalBenefit()) {
			marginalBen = getTotalBenefit();
			vmax = getBestNode(t);
			probeNode(vmax);
		}
		m.insert(
				pair<d_type, int>(yellowNodes[i],
						yBenefit[i] + getTotalBenefit()));
	}
	return m;
}

//map<d_type, int> GraphLabel::benefitMappingSlow() {

//}
