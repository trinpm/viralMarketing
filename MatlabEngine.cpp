/*
 * MatlabEngine.cpp
 *
 *  Created on: Aug 7, 2016
 *      Author: trinpm
 */

#include "MatlabEngine.h"
#include <stdio.h>

#define BUFSIZE 2560

int MatlabEngine::initEngine() {
	//start engine
	if (!(ep = engOpen(""))) {
		cout << "cannot start MATLAB engine" << endl;
		return 0;
	}
	cout << "started MATLAB engine" << endl;
	return 1;
}

void MatlabEngine::loadModel(string modelPath) {
	string load_model = "load('" + modelPath + "');";
	engEvalString(ep, load_model.c_str());
	cout << "loaded random forest model: " << modelPath << endl;
}

double MatlabEngine::predict(vector<double> datapoint) {

	// convert vector of datapoint to array
	int ARR_SIZE = datapoint.size();
	double data[ARR_SIZE];
	std::copy(datapoint.begin(), datapoint.end(), data);

	// send data to Matlab
	mxArray *T = NULL;
	T = mxCreateDoubleMatrix(1, ARR_SIZE, mxREAL);
	memcpy((void *) mxGetPr(T), (void *) data, sizeof(data));

	// catch any output to buffer
	char buffer[BUFSIZE + 1];
	buffer[BUFSIZE] = '\0';
	engOutputBuffer(ep, buffer, BUFSIZE);

	// run matlab command
	engPutVariable(ep, "T", T);
	engEvalString(ep, "predict = forest_model.predict(T);");
	engEvalString(ep, "predict = str2num(predict{1,1});");
	engEvalString(ep, "sprintf('predict is: %d', predict)");
//	cout << buffer << endl;

	double predict = parseResult_single(buffer);
//	cout << "predict: " << predict << endl;

// free variables
	mxDestroyArray(T);

	return predict;
}

double MatlabEngine::parseResult_single(string result) {
	stringstream ss;
	ss << result;
	double predict;
	std::string temp;
	while (std::getline(ss, temp, ' ')) {
		if (stringstream(temp) >> predict) {
			break;
		}
	}

//	cout << "predict: " << predict << endl;
	return predict;
}

double MatlabEngine::predict(vector<vector<double>> datapoint) {

	// convert vector of datapoint to array
	int row = datapoint.size();
	int col = datapoint[0].size();

//	cout << "array of data: " << endl;
	double data[row][col];
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			data[i][j] = datapoint[i][j];
//			cout << data[i][j] << " ";
		}
//		cout << endl;
	}

	// send data to Matlab
	mxArray *T = NULL;
	T = mxCreateDoubleMatrix(col, row, mxREAL);
	memcpy((void *) mxGetPr(T), (void *) data, sizeof(data));

	// catch any output to buffer
	char buffer[BUFSIZE + 1];
	buffer[BUFSIZE] = '\0';
	engOutputBuffer(ep, buffer, BUFSIZE);

	// run matlab command
	engPutVariable(ep, "T", T);
	engEvalString(ep, "T = T';");
	engEvalString(ep, "result = forest_model.predict(T)");
	engEvalString(ep, "[M,I] = max(result)");
	int idx = parseResult(buffer);

// free variables
	mxDestroyArray(T);

	return idx;
}

double MatlabEngine::parseResult(string result) {
	stringstream ss;
	ss << result;
	double number;
	std::string temp;
	int index;
	while (std::getline(ss, temp, ' ')) {
		if (stringstream(temp) >> number) {
			index = number;
		}
	}

//	cout << "index: " << index << endl;
	return index;
}

vector<double> MatlabEngine::predictVector(vector<vector<double>> datapoint) {

	// convert vector of datapoint to array
	int row = datapoint.size();
	int col = datapoint[0].size();

//	cout << "array of data: " << endl;
	double data[row][col];
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			data[i][j] = datapoint[i][j];
//			cout << data[i][j] << " ";
		}
//		cout << endl;
	}

	// send data to Matlab
	mxArray *T = NULL;
	T = mxCreateDoubleMatrix(col, row, mxREAL);
	memcpy((void *) mxGetPr(T), (void *) data, sizeof(data));

	// catch any output to buffer
	char buffer[BUFSIZE + 1];
	buffer[BUFSIZE] = '\0';
	engOutputBuffer(ep, buffer, BUFSIZE);

	// run matlab command
	engPutVariable(ep, "T", T);
	engEvalString(ep, "T = T';");
	engEvalString(ep, "result = forest_model.predict(T)");

//	cout << "buffer: " << buffer << endl;
	vector<double> results;
	results = parseVectorResult(buffer);

	// free variables
	mxDestroyArray(T);

	return results;
}

vector<double> MatlabEngine::parseVectorResult(string result) {
	stringstream ss;
	ss << result;
	double number;
	std::string temp;

	vector<double> results;
	while (std::getline(ss, temp, ' ')) {
		if (stringstream(temp) >> number) { //return 0 if temp is not number
			results.push_back(number);
		}
	}

//	cout << "index: " << index << endl;
	return results;
}

void MatlabEngine::closeEngine() {
	cout << "closed MATLAB engine" << endl;
	//close engine
	engEvalString(ep, "close;");
	engClose(ep);
}
