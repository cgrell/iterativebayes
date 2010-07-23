/*
 * main.cpp
 *
 *  Created on: Jul 13, 2010
 *      Author: cgrell
 */

#include <iostream>
#include <ctime>
#include <cstdio>
#include <fstream>
#include <getopt.h>
#include <sstream>
#include <map>
#include <vector>
#include <readNetworkData.h>
#include <iterativeBayes.h>
#include <readPathwayData.h>
#include <setPriorProbability.h>
#include <bayesianUpdate.h>

using namespace std;

int iterations = 1;
vector<NetworkDataSet> bayesNetArray;
vector<BayesData> distByNet;

void readPathway(string pathwayname, const vector<string>& path, const vector<NetworkDataSet>& datasets) {
	vector<NetworkDataSet>::const_iterator iter = datasets.begin();
	while (iter != datasets.end()) {
	NetworkDataSet netName = *iter;
	Pathway *target = new Pathway(pathwayname,path,datasets);
	Probability *odds_table = new Probability(netName,*target);
	BayesUpdate *distribution = new BayesUpdate(*odds_table,netName,*target);
	BayesData *initialize = new BayesData(netName,*target,*odds_table,*distribution);
	distByNet.push_back(*initialize);
	++iter;
	}
	return;
}

void readNetworks(const vector<string>& netfiles, vector<NetworkDataSet>& dataset){
	vector<string>::const_iterator iter = netfiles.begin();
	while (iter != netfiles.end()) {
		string graphName = *iter;
		cout << graphName <<endl;
		NetworkDataSet *interactionGraph = new NetworkDataSet(graphName);
		dataset.push_back(*interactionGraph);
		++iter;
	}
	return;
}

//\todo need to add the case for synthesizing multiple networks into odds
void iterativeBayes(vector<BayesData>& dists, int iterations) {
	vector<BayesData>::const_iterator iter;
	for (int i = 0; i < iterations; i++) {
		iter = dists.begin();
		while (iter != dists.end()) {
			iter->algorithms->createHistogram(iter->logodds->geneodds);
			iter->algorithms->updateDistribution(iter->logodds->geneodds);
			iter->algorithms->printGeneProbability(iter->logodds->geneodds);
			++iter;
		}
	}
	return;
}

void print_usage(int signal)
{
  cerr << "iterativeBayes"
#ifdef VERSION
       << " -- " << VERSION
#endif
       << endl
       << "Usage:" << endl
       <<"  IterativeBayes [options] -p pathway.txt -d data.txt -i #iterations" << endl
       << "C++ program for finding new genes in a known pathway by iterative Bayes" << endl;
  exit(signal);
}

void die(string msg)
{
  cerr << msg << endl;
  exit(-1);
}

void setIterations(char * c) {
	  iterations = strtol(c, NULL ,10);
	  return;
}


int main(int argc, char *argv[])
{
	time_t start = time(NULL);

  const char* const short_options = "hp:b:c:d:i:o:v";
  const struct option long_options[] = {
    { "batch", 0, NULL, 'b' },
    { "config", 0, NULL, 'c' },
    { "pathway", 0, NULL, 'p' },
    { "data", 0, NULL, 'd' },
    { "iterations", 0, NULL, 'i'},
    { "output", 0, NULL, 'o' },
    { "help", 0, NULL, 'h' },
    { NULL, 0, NULL, 0 }
  };
  int next_options;

  string pathwayFilename;
  string batchPrefix;
  string configFile;
  string dataFileList;
  string actOutFile;
  string gene;
  string network;
  string pathwayname;
  vector<string> networkset;
  vector<string> pathwayset;


  // /////////////////////////////////////////////////
  // Read in command line options
  //
  do {
    next_options=getopt_long(argc,argv,short_options,long_options,NULL);
    switch (next_options){
    case 'h': print_usage(EXIT_SUCCESS); break;
    case 'p': pathwayFilename = optarg; break;
    case 'b': batchPrefix = optarg; break;
    case 'c': configFile = optarg; break;
    case 'd': dataFileList = optarg; break;
    case 'i': setIterations(optarg); break;
    case 'o': actOutFile = optarg; break;
    }
  } while (next_options != -1);


  cout << "pathway " << pathwayFilename << endl;
  cout << "data file list " << dataFileList << endl;
  cout << "number of iterations " << iterations << endl;

  ostream* outstream;
  ofstream outFileStream;
  if (actOutFile == "") {
    outstream = &cout;
  } else {
    outFileStream.open(actOutFile.c_str());
    if (!outFileStream.is_open()) {
      die("couldn't open output file");
    }
    outstream = &outFileStream;
  }

  // /////////////////////////////////////////////////
  // Verify that command line options are valid
  //

  if(dataFileList.length() == 0)
    {
      cerr << "Missing required arguments" << endl;
      print_usage(EXIT_FAILURE);
    }
  /*if(batchPrefix.length() == 0)
    {
      cerr << "Missing batch prefix" << endl;
      print_usage(EXIT_FAILURE);
    }*/
  /*if(configFile.length() == 0)
    {
      cerr << "Missing configuration file" << endl;
      print_usage(EXIT_FAILURE);
    }*/
/*
  // /////////////////////////////////////////////////
  // Load configuration
  RunConfiguration conf(configFile);
  if (conf.evidenceSize() == 0) {
    die("Must have at least one evidence file in configuration.");
  }
*/
  // /////////////////////////////////////////////////
  // Load pathway
  ifstream pathwayStream;
  pathwayStream.open(pathwayFilename.c_str());
  if (!pathwayStream.is_open()) {
    die("Could not open pathway stream");
  }
  //PathwayTab pathway = PathwayTab::create(pathwayStream, conf.pathwayProps());

  // /////////////////////////////////////////////////
  // Load network data file list
  ifstream dataStream;
  dataStream.open(dataFileList.c_str());
  if (!dataStream.is_open()) {
    die("Could not open data file stream");
  }

//\todo add capability for multiple pathways
pathwayStream >> pathwayname;
while (pathwayStream >> gene) pathwayset.push_back(gene);
while (getline(dataStream,network)) networkset.push_back(network);
readNetworks(networkset,bayesNetArray);
readPathway(pathwayname,pathwayset,bayesNetArray);
iterativeBayes(distByNet,iterations);
//\todo add delete for data structure

cout << pathwayname << endl;
for (vector<string>::iterator i = pathwayset.begin(); i != pathwayset.end(); i++) {
cout << *i << endl;
}
time_t timeinterval = (time(NULL) - start);
cout << timeinterval << " secs"<< endl;
return 0;
}
