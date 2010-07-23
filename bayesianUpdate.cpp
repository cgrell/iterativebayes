/*
 * bayesianUpdate.cpp
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
#include <bayesianUpdate.h>
#include <cmath>
#include <algorithm>

using namespace std;

BayesUpdate::BayesUpdate(Probability& logodds, const NetworkDataSet& genenetwork,
		const Pathway& target) {
	logodds_table = &logodds;
	data_table = &genenetwork;
	scoreGeneNeighbors();
}

double BayesUpdate::probabilitytoLogodds(double prob) {
	double odds = prob / (1 - prob);
	return log(odds);
}

double BayesUpdate::logoddstoProbability(double logodds) {
	double odds = exp(logodds);
	return odds / (1 + odds);
}

void BayesUpdate::scoreGeneNeighbors(const NetworkDataSet& genenetwork) {
	pathwayWeight = 0.0;
	totalWeight = 0.0;
	averageWeight = 0.0;
	double prob = 0.0;
	map<string, int>::const_iterator geneindex =
			genenetwork.allvertices.begin();
	averageWeight = calcAverageWeights(genenetwork);
	double baseline = computeBaselineProbability(logodds_table->geneodds);
	double logodds = probabilitytoLogodds(baseline);
	while (geneindex != genenetwork.allvertices.end()) {
		string gene = geneindex->first;
		pair<string, string> connection;
		map<pair<string, string> , double>::const_iterator netindex =
				genenetwork.inputdataset.begin();
		while (netindex != genenetwork.inputdataset.end()) {
			connection = netindex->first;
			if (connection.first == gene) {
				totalWeight	+= netindex->second;
				double odds = logodds_table->geneodds[gene];
				prob = logoddstoProbability(odds);
				pathwayWeight += netindex->second*prob;
			}
			++netindex;
		}
		if ((averageWeight > 0) && (baseline > 0) && (baseline < 1)) {
			pathwayWeight /= averageWeight;
			totalWeight /= averageWeight;
			prob = (exp(totalWeight * logodds) - exp(pathwayWeight * logodds))
					/ (exp(totalWeight * logodds) - 1);
			if (prob == 0) {
				break;
			}
			genescore[gene] = -1 * log(prob);
		}
		++geneindex;
	}
	normalizeScores(genescore);
	return;
}

void BayesUpdate::scoreGeneNeighbors() {
	pathwayWeight = 0.0;
	totalWeight = 0.0;
	averageWeight = 0.0;
	double prob = 0.0;
	averageWeight = calcAverageWeights(*data_table);
	map<string, int>::const_iterator geneindex =
			data_table->allvertices.begin();
	double baseline = computeBaselineProbability(logodds_table->geneodds);
	double logodds = probabilitytoLogodds(baseline);
	while (geneindex != data_table->allvertices.end()) {
		string gene = geneindex->first;
		pair<string, string> connection;
		map<pair<string, string> , double>::const_iterator netindex =
				data_table->inputdataset.begin();
		while (netindex != data_table->inputdataset.end()) {
			connection = netindex->first;
			if (connection.first == gene) {
				totalWeight	+= netindex->second;
				double odds = logodds_table->geneodds[gene];
				prob = logoddstoProbability(odds);
				pathwayWeight += netindex->second*prob;
			}
			++netindex;
		}
		if ((averageWeight > 0) && (baseline > 0) && (baseline < 1)) {
			pathwayWeight /= averageWeight;
			totalWeight /= averageWeight;
			prob = (exp(totalWeight * logodds) - exp(pathwayWeight * logodds))
					/ (exp(totalWeight * logodds) - 1);
			if (prob == 0) {
				break;
			}
			genescore[gene] = -1 * log(prob);
		}
		++geneindex;
	}
	normalizeScores(genescore);
	return;
}

double BayesUpdate::calcAverageWeights(const NetworkDataSet& net) {
	double sumofweights = 0.000001;
	map<pair<string, string> , double>::const_iterator iter = net.inputdataset.begin();
	while (iter != net.inputdataset.end()) {
		double weight = iter->second;
		sumofweights += weight;
		++iter;
	}
	return sumofweights / net.inputdataset.size();
}

double BayesUpdate::computeBaselineProbability(map<string, double>& oddslist) {
	map<string, double>::const_iterator iter = oddslist.begin();
	double totalProbPathway = 0.000001;
	size_t numGenes = oddslist.size();
	while (iter != oddslist.end()) {
		string gene = iter->first;
		totalProbPathway += logoddstoProbability(oddslist[gene]);
		++iter;
	}
	return totalProbPathway / numGenes;
}

double BayesUpdate::normalizeScores(map<string,double>& scores) {
	double maxscore = 0.000001;
	map<string,double>::iterator iter = scores.begin();
	double score = 0.0;
	while (iter != scores.end()) {
		score = iter->second;
		maxscore = max(score,maxscore);
		++iter;
	}
	iter = scores.begin();
	while (iter != scores.end()) {
		score = iter->second;
		string gene = iter->first;
		scores[gene] = score /= maxscore;
		++iter;
	}
	return maxscore;
}

void BayesUpdate::createHistogram(map<string, double>& oddslist) {
	map<string,double>::iterator iter = oddslist.begin();
	double prob = 0.0;
	string gene;
	while (iter != oddslist.end()) {
		gene = iter->first;
		prob = logoddstoProbability(oddslist[gene]);
		posdistribution[gene] += prob;
		negdistribution[gene] += 1 - prob;
		++iter;
	}
	pos_smoothed = exponentialSmoothing(posdistribution, num_bins, beta);
	neg_smoothed = exponentialSmoothing(negdistribution, num_bins, beta);
	return;
}

 vector<double> BayesUpdate::exponentialSmoothing(map<string,double> dist, int bincount, int beta) {
	double weight;
	vector<double> smoothed_dist;
	vector<double> binned_dist;
	for (int i = 0; i < bincount; i++) {
		smoothed_dist.push_back(0.0);
		binned_dist.push_back(0.0);
	}
	map<string,double>::iterator iter = dist.begin();
	while (iter != dist.end()) {
		binned_dist[(int)floor(dist[iter->first]*bincount)] += dist[iter->first];
		++iter;
	}
	for (int i = 0; i < bincount; i++) {
		weight = (bincount / beta) * (2 - exp(-1 * beta * i / bincount) - exp(-1 * beta
						* ((bincount - i) / bincount)));
		if (weight == 0) break;
		for (int j = 0; j < i; j++) {
			smoothed_dist[j] += binned_dist[i] * bincount * exp(-1 * beta * i / bincount)
			* (exp(beta * (j + 1) / bincount) - exp(beta * j / bincount)) / (beta
					* weight);
		}
		for (int j = i; j < bincount; j++) {
			smoothed_dist[j] += binned_dist[i] * bincount * exp(beta * i / bincount)
			* (exp(-1 * beta * j / bincount)
					- exp(-1 * beta * (j + 1) / bincount)) / (beta * weight);
		}
	}
	double total_wt = 0.000001;
	for (int i = 0; i < bincount; i++) {
		total_wt += smoothed_dist[i];
	}
	if (total_wt > 0) {
		for (int i = 0; i < bincount; i++) {
			smoothed_dist[i] = smoothed_dist[i] / total_wt;
		}
	}
	return smoothed_dist;
}

void BayesUpdate::updateDistribution(map<string, double>& oddslist) {
	double logodds = 0.0;
	double pos_prob = 0.0;
	double neg_prob = 0.0;
	double score = 0.0;
	int bin = 0;
	map<string, double>::iterator iter = oddslist.begin();
	while (iter != oddslist.end()) {
		string gene = iter->first;
		logodds = oddslist[gene];
		for (int i = 0; i < num_graphs; i++) {
			score = genescore[gene];
			bin = (int) floor(num_bins * score);
			pos_prob = pos_smoothed[bin];
			neg_prob = neg_smoothed[bin];
			if ((pos_prob > 0) && (neg_prob > 0)) {
				logodds += (log(pos_prob) - log(neg_prob));
			} else if (pos_prob > 0) {
				logodds += (log(pos_prob) - log(0.000001));
			}
		}
		oddslist[gene] = logodds;
		++iter;
	}
	return;
}

void BayesUpdate::printGeneProbability(map<string, double>& oddslist) {
	double prob = 0.0;
	map<string, double>::iterator iter = oddslist.begin();
	while (iter != oddslist.end()) {
	string gene = iter->first;
	prob = logoddstoProbability(oddslist[gene]);
	cout << gene << "\t" << prob << endl;
	++iter;
	}
}
