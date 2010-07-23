/*
 * readNetworkData.cpp
 *
 *  Created on: Jul 13, 2010
 *      Author: cgrell
 */
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <readNetworkData.h>

using namespace std;

vector<string> split(const string& s);
typedef string::size_type string_size;

NetworkDataSet::NetworkDataSet(string source) {
	datasource = source;
	processInputGraph(datasource);
	numvertices = allvertices.size();
}

void NetworkDataSet::processInputGraph(string source) {
	ifstream infile(source.c_str());
	string line;
	pair<string, string> connection;
	double defaultWeight = 0.0;
	vector<string> tokens;
	while (getline(infile, line)) {
		tokens = split(line);
		pair<string, string> connection;
		connection.first = tokens[0];
		connection.second = tokens[1];
		defaultWeight = strtod(tokens[2].c_str(), NULL);
		if (!(defaultWeight > 0)) defaultWeight = 1.0;
		inputdataset[connection] = defaultWeight;
		allvertices[tokens[0]] = 1;
		allvertices[tokens[1]] = 1;
		//\todo check if forward and backward pair are redundant
	}
}

vector<string> split(const string& s) {
	vector<string> ret;
	string_size i = 0;
	while (i != s.size()) {
		while (i != s.size() && isspace(s[i]))
			++i;
		string_size j = i;
		while (j != s.size() && !isspace(s[j]))
			++j;
		if (i != j) {
			ret.push_back(s.substr(i, j - i));
			i = j;
		}
	}
	return ret;
}

void NetworkDataSet::printGeneList(map<string,int> glist) {
	int counter = 0;
	for (map<string,int>::iterator i = glist.begin(); i != glist.end(); i++) {
	cout << i->first << endl;
	counter++;
	}
	cout << counter << endl;
}
