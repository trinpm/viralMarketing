/*
 * PythonWrapper.cpp
 *
 *  Created on: Mar 9, 2017
 *      Author: trinpm
 */

#include "PythonWrapper.h"

PythonWrapper::PythonWrapper() {

}

void PythonWrapper::init() {
	Py_Initialize();
}

void PythonWrapper::setPath(string modulePath) {
	PyObject *sys = PyImport_ImportModule("sys");
	PyObject *path = PyObject_GetAttrString(sys, "path");
	PyList_Append(path, PyString_FromString(modulePath.c_str()));
}

void PythonWrapper::importModule(string module) {
	PyObject *temp;
	temp = PyString_FromString(module.c_str());
	pModule = PyImport_Import(temp);

	if (pModule != NULL) {
		cout << "loaded Python module: " << module << endl;
	} else {
		PyErr_Print();
	}

	Py_DECREF(temp);
}

int PythonWrapper::loadModel(string func, string arg) {
	PyObject *pFunc, *pValue, *pArgs, *temp;

	int success = 0;
	if (pModule != NULL) {
		pFunc = PyObject_GetAttrString(pModule, func.c_str());

		if (pFunc && PyCallable_Check(pFunc)) {

			//set arg to be passed to Python
			pArgs = PyTuple_New(1);
			pValue = PyString_FromString(arg.c_str());

			//insert arg to tuple
			PyTuple_SetItem(pArgs, 0, pValue);

			//call function with tuples containing args
			pValue = PyObject_CallObject(pFunc, pArgs);

			if (pValue != NULL) {
				success = 1;
				cout << "Result of call: " << PyFloat_AsDouble(pValue) << endl;
			} else {
				cerr << "cannot load model!!" << endl;
				PyErr_Print();
			}
		} else {
			cerr << "cannot load function" << endl;
			PyErr_Print();
		}
	}
	Py_DECREF(pFunc);
	Py_DECREF(pValue);
	Py_DECREF(pArgs);
	return success;
}

double PythonWrapper::predict(string func, vector<double> dp) {
	double pred = 0.0;
	PyObject *pVec, *pValue;

//	cout << "1" << endl;
	//transfer c++ vector to python tuple
	pVec = PyTuple_New(dp.size());
	for (int i = 0; i < dp.size(); i++) {
		pValue = PyFloat_FromDouble(dp[i]);

		if (pValue != NULL) {
			PyTuple_SetItem(pVec, i, pValue);
		} else {
			cerr << "cannot convert array value" << endl;
			return 0.0;
		}
	}
//	cout << "2" << endl;
	PyObject *pArgs;
	pArgs = PyTuple_New(1);
	// assign pvec to tuple
	PyTuple_SetItem(pArgs, 0, pVec);

//	cout << "3" << endl;
	PyObject *pFunc;
	pFunc = PyObject_GetAttrString(pModule, func.c_str());
	if (pFunc && PyCallable_Check(pFunc)) {
		// call predict function
		pValue = PyObject_CallObject(pFunc, pArgs);
		if (pValue != NULL) {
			pred = PyFloat_AsDouble(pValue);
			cout << "Result of call: " << PyFloat_AsDouble(pValue) << endl;
		}
		else {
			cerr << "no result returned from Python" << endl;
			PyErr_Print();
		}
	} else {
		cerr << "cannot load function" << endl;
		PyErr_Print();
	}

	Py_DECREF(pVec);
	Py_DECREF(pValue);
	Py_DECREF(pArgs);
	Py_DECREF(pFunc);
	return pred;
}

void PythonWrapper::clean(){
	Py_DECREF(pModule);
}
