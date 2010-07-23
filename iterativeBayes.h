/*
 * iterativeBayes.h
 *
 *  Created on: Jul 13, 2010
 *      Author: cgrell
 */
#ifndef HEADER_ITERATIVEBAYES_H
#define HEADER_ITERATIVEBAYES_H

#include <readNetworkData.h>
#include <readPathwayData.h>
#include <setPriorProbability.h>
#include <bayesianUpdate.h>

class BayesData {
public:
	BayesData(const NetworkDataSet& net, const Pathway& pathway, Probability& odds, BayesUpdate& bayesscore);
	const NetworkDataSet * graphdata;
	const Pathway * genesubset;
	Probability * logodds;
	BayesUpdate * algorithms;
};





#endif
