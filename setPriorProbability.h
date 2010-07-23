/*
 * setPriorProbability.h
 *
 *  Created on: Jul 13, 2010
 *      Author: cgrell
 */
#ifndef HEADER_SETPRIORPROBABILITY_H
#define HEADER_SETPRIORPROBABILITY_H

#include <string>
#include <vector>
#include <map>
#include <readNetworkData.h>
#include <readPathwayData.h>

class Probability {
public:
	Probability(const NetworkDataSet& network, const Pathway& path);
	string netsource;
	map<string,double> geneodds;
	void initializePrior(const NetworkDataSet& network, const Pathway& path, double falsepos, double falseneg);
	void printGeneProbability(const map<string,double> logodds);
};


#endif
