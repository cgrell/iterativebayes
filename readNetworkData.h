/*
 * readNetworkData.h
 *
 *  Created on: Jul 13, 2010
 *      Author: cgrell
 *
 *
*/

#ifndef HEADER_READNETWORKDATA_H
#define HEADER_READNETWORKDATA_H


#include <vector>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <map>

using std::string;
using std::pair;
using std::map;

class NetworkDataSet {
public:
	NetworkDataSet(string source);
	map<pair<string,string>,double > inputdataset;
	map<string, int> allvertices;
	string datasource;
	map<string,int>::size_type numvertices;
	void processInputGraph(string source);
	void printGeneList(map<string,int> genekeys);
};


#endif
