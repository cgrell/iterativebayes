/*
 * bayesianUpdate.h
 *
 *  Created on: Jul 13, 2010
 *      Author: cgrell
 */
#ifndef HEADER_BAYESIANUPDATE_H
#define HEADER_BAYESIANUPDATE_H

#include <string>
#include <vector>
#include <map>
#include <readNetworkData.h>
#include <readPathwayData.h>
#include <setPriorProbability.h>

using std::string;
using std::vector;

const int num_bins = 200;
const int num_graphs = 2;
const int beta = 30;

class BayesUpdate {
public:
	BayesUpdate(Probability& logodds,const NetworkDataSet& genenetwork,
			const Pathway& target);
	Probability *logodds_table;
	const NetworkDataSet *data_table;
	double totalWeight;
	double pathwayWeight;
	double averageWeight;
	vector<double> pos_smoothed;
	vector<double> neg_smoothed;
	map<string,double> genescore;
	map<string,double> posdistribution;
	map<string,double> negdistribution;
	double probabilitytoLogodds(double prob);
	double logoddstoProbability(double logodds);
	void scoreGeneNeighbors(const NetworkDataSet& genenetwork);
	void scoreGeneNeighbors();
	double calcAverageWeights(const NetworkDataSet& genenetwork);
	double computeBaselineProbability(map<string, double>& oddslist);
	double normalizeScores(map<string,double>& scores);
	void createHistogram(map<string, double>& oddslist);
	vector<double> exponentialSmoothing(map<string,double> dist, int res, int beta);
	void updateDistribution(map<string, double>& oddslist);
	void printGeneProbability(map<string, double>& oddslist);
};
#endif
