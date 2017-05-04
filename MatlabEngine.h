/*
 * MatlabEngine.h
 *
 *  Created on: Aug 7, 2016
 *      Author: trinpm
 */

#ifndef MATLABENGINE_H_
#define MATLABENGINE_H_

#include "engine.h"
#include <iostream>
#include <string>
#include <sstream>
#include <cstring>
#include <vector>

using namespace std;

class MatlabEngine {

	Engine *ep;
	mxArray *T = NULL;

public:
	int initEngine();

	void loadModel(string modelPath);

	double predict(vector<double> datapoint);

	double predict(vector<vector< double>> datapoint);

	vector<double> predictVector(vector<vector< double>> datapoint);

	double parseResult(string str);

	double parseResult_single(string str);

	vector<double> parseVectorResult(string result);

	void closeEngine();

};

#endif /* MATLABENGINE_H_ */
