/*
 * readPathwayData.h
 *
 *  Created on: Jul 13, 2010
 *      Author: cgrell
 */

#ifndef HEADER_READPATHWAYDATA_H
#define HEADER_READPATHWAYDATA_H

#include <string>
#include <map>
#include <vector>
#include <readNetworkData.h>

using std::string;
using std::vector;

class Pathway {
public:
	Pathway(string name, const vector<string>& genelist, const vector<NetworkDataSet>& graphset);
	string pathwayName;
	map<string,int> path;
	map<string,int>::size_type genecount;
	map<string,int> pathInNet;
	map<string,int> pathNotInNet;
	void separateGroups(const map<string,int>& pathway, const vector<NetworkDataSet>& graphset);
	void printPathway(const map<string,int> genelist);
};


#endif
