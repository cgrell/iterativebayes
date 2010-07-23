/*
 * setPriorProbability.cpp
 *
 *  Created on: Jul 13, 2010
 *      Author: cgrell
 */

#include <string>
#include <vector>
#include <map>
#include <readNetworkData.h>
#include <readPathwayData.h>
#include <setPriorProbability.h>
#include <cmath>

using namespace std;

Probability::Probability(const NetworkDataSet& network, const Pathway& path) {
	netsource = network.datasource;
	initializePrior(network, path, 0.1, 2.0);
	//printGeneProbability(geneodds);
}
void Probability::initializePrior(const NetworkDataSet& graph,
		const Pathway& pathway, double falsepos, double falseneg) {
	double prob = 0.0;
	typedef map<string, int>::size_type map_size;
	map_size numPathGenes = pathway.genecount;
	map_size numNetGenes = graph.numvertices;
	map<string, int>::const_iterator iter = graph.allvertices.begin();
	while (iter != graph.allvertices.end()) {
		string gene = iter->first;
		if ((pathway.path.find(gene)) != pathway.path.end()) {
			geneodds[gene] = log((1 - falsepos) / falsepos);
		} else {
			prob = (falseneg - 1 + falsepos) * numPathGenes
					/ (numNetGenes - numPathGenes);
			geneodds[gene] = log(prob/(1 - prob));
		}
		++iter;
	}
}
void Probability::printGeneProbability(const map<string, double> logodds) {
	size_t numGenes = logodds.size();
	cout << numGenes << endl;
	map<string, double>::const_iterator iter = logodds.begin();
	while (iter != logodds.end()) {
		cout << iter->first << "   " << iter->second << endl;
		++iter;
	}
}
