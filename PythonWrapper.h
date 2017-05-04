/*
 * PythonWrapper.h
 *
 *  Created on: Mar 9, 2017
 *      Author: trinpm
 */

#include <iostream>
#include <vector>
#include <string>
#include <Python.h>

using namespace std;

#ifndef PYTHONWRAPPER_H_
#define PYTHONWRAPPER_H_

class PythonWrapper {
private:
	PyObject *pModule;

public:
	PythonWrapper();

	void init();

	void setPath(string modulePath);

	void importModule(string modulePath);

	int loadModel(string func, string arg);

	double predict(string func, vector<double> datapoint);

	void clean();

};

#endif /* PYTHONWRAPPER_H_ */
