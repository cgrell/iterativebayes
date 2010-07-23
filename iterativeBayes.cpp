/*
 * iterativeBayes.cpp
 *
 *  Created on: Jul 13, 2010
 *      Author: cgrell
 */

#include <readNetworkData.h>
#include <readPathwayData.h>
#include <setPriorProbability.h>
#include <bayesianUpdate.h>
#include <iterativeBayes.h>

using namespace std;

BayesData::BayesData(const NetworkDataSet& net, const Pathway& pathway, Probability& odds, BayesUpdate& bayesscore) {
		graphdata = &net;
		genesubset = &pathway;
		logodds = &odds;
		algorithms = &bayesscore;
}
