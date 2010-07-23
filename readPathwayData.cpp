/*
 * readPathwayData.cpp
 *
 *  Created on: Jul 13, 2010
 *      Author: cgrell
 */

#include <string>
#include <map>
#include <vector>
#include <readPathwayData.h>
#include <readNetworkData.h>

using namespace std;

Pathway::Pathway(string name, const vector<string>& genelist, const vector<
		NetworkDataSet>& graphset) {
	typedef map<string,int>::size_type map_size;
	pathwayName = name;
	for (map_size i = 0; i < genelist.size(); i++) {
		path[genelist[i]] = 1;
	}
	genecount = genelist.size();
	separateGroups(path, graphset);
	//printPathway(pathNotInNet);
}

void Pathway::separateGroups(const map<string, int>& pathway, const vector<
		NetworkDataSet>& nets) {
	map<string, int>::const_iterator nextgene = pathway.begin();
	while (nextgene != pathway.end()) {
		vector<NetworkDataSet>::const_iterator iter = nets.begin();
		while (iter != nets.end()) {
			string gene = nextgene->first;
			if ((iter->allvertices.find(gene)) != iter->allvertices.end()) {
				pathInNet[gene] = 1;
			} else
				pathNotInNet[gene] = 1;
			++iter;
		}
		++nextgene;
	}
}

void Pathway::printPathway(const map<string, int> genelist) {
	map<string, int>::const_iterator nextgene = genelist.begin();
	while (nextgene != genelist.end()) {
		cout << nextgene->first << endl;
		++nextgene;
	}
}
