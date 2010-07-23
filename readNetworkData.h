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

//Reads the data from all data sources, located in files listed on
//consecutive lines of $names_file. Each data source is assumed to
//be represented as an undirected weighted graph, with columns 0
//and 1 the vertices, and column 2 the weight.  If there is no column 2,
//then the weight is assumed to be 1. If there is an
//edge from i to j and from j to i, then the weight is taken to be
//the maximum of the weight for the edge from i to j and from j to
//i. The final data structure, @adj_list, is an array of directed
//graphs.  Value $adj_list[$source]{$in}{$out} is the weight on
//the edge from $in to $out in data source $source, and is
//undefined if there is no such edge.  The routine outputs a
//pointer to @adj_list, a pointer to %all_vertices, a list of
//genes/metagenes that appear in the data sources, a pointer to
//an array @sources which gives the data sources corresponding to $source,
//and a pointer to an array @avg_wts which gives the av. weights for each.


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
