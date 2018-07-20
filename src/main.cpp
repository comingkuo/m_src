#define Verbose 1

#include "Core.h"
#include "Core.cpp"
#include "dataManager/dataManager.h"

#include "unistd.h"
#include <stdio.h>
#include <stdlib.h>
#include "unitTest.h"
#include "algorithms/dimEst.h"
#include "algorithms/pageRank.h"
#include "tools/argParser.h"
#include "algorithms/maxAggregator.h"
#include "general.h"
#include "algorithms/SSSP.h"
#include "algorithms/ConnectedComponent.h"
#include "algorithms/Nstep.h"

using namespace std;

int main(int argc, char** argv) {

	argParser argp;
	MizanArgs myArgs = argp.parse(argc, argv);

	char ** inputBaseFile = argp.setInputPaths(myArgs.fs, myArgs.clusterSize,
			myArgs.graphName, myArgs.hdfsUserName, myArgs.partition);

#ifdef Verbose
	time_t begin_time = time(NULL);
#endif

	bool groupVoteToHalt;
	edgeStorage storageType;
	distType partType;

	int myWorkerID;

	if (myArgs.algorithm == 1) {// PageRank
		groupVoteToHalt = true;
		storageType = OutNbrStore;
		pageRank us(myArgs.superSteps);
		pageRankCombiner prc;

		Mizan<mLong, mDouble, mDouble, mLong> * mmk = new Mizan<mLong, mDouble,
				mDouble, mLong>(myArgs.communication, &us, storageType,
				inputBaseFile, myArgs.clusterSize, myArgs.fs, myArgs.migration, myArgs.threshold);

		mmk->registerMessageCombiner(&prc);

		mmk->setVoteToHalt(groupVoteToHalt);

		string output;
		output.append("/user/");
		output.append(myArgs.hdfsUserName.c_str());
		output.append("/m_run_output/");
		output.append(myArgs.graphName.c_str());
		mmk->setOutputPath(output.c_str());

		//User Defined aggregator
		char * maxAgg = "maxAggregator";
		maxAggregator maxi;
		mmk->registerAggregator(maxAgg, &maxi);

		mmk->run(argc, argv);

		myWorkerID = mmk->getPEID();

		delete mmk;

	}  else if (myArgs.algorithm == 2) { //SSSP
		groupVoteToHalt = false;
		storageType = OutNbrStore;
		SSSP sp(1,myArgs.superSteps);
		SSSPCombiner spc;


		Mizan<mLong, mLong, mLong, mLong> * mmk = new Mizan<mLong, mLong,
				mLong, mLong>(myArgs.communication, &sp, storageType,
				inputBaseFile, myArgs.clusterSize, myArgs.fs, myArgs.migration, myArgs.threshold);

		//mmk->registerMessageCombiner(&spc);

		mmk->setVoteToHalt(groupVoteToHalt);

		string output;
		output.append("/user/");
		output.append(myArgs.hdfsUserName.c_str());
		output.append("/m_run_output/");
		output.append(myArgs.graphName.c_str());
		mmk->setOutputPath(output.c_str());

		mmk->run(argc, argv);
		myWorkerID = mmk->getPEID();
		delete mmk;
	}
	else if (myArgs.algorithm == 3) { //Connected Component
		groupVoteToHalt = true;
		storageType = OutNbrStore;
		ConnectedComponent k(myArgs.superSteps);

		Mizan<mLong, mLong, mLong, mLong> * mmk = new Mizan<mLong, mLong,
			mLong, mLong>(myArgs.communication, &k, storageType,
				inputBaseFile, myArgs.clusterSize, myArgs.fs, myArgs.migration, myArgs.threshold);

		string output;
		output.append("/user/");
		output.append(myArgs.hdfsUserName.c_str());
		output.append("/m_run_output/");
		output.append(myArgs.graphName.c_str());
		mmk->setOutputPath(output.c_str());

		mmk->setVoteToHalt(groupVoteToHalt);
		mmk->run(argc, argv);
		myWorkerID = mmk->getPEID();
		delete mmk;

	}  else if (myArgs.algorithm == 4) { //N-steps
		groupVoteToHalt = true;
		storageType = OutNbrStore;
		Nstep NS(1, myArgs.superSteps);


		Mizan<mLong, mLong, mLong, mLong> * mmk = new Mizan<mLong, mLong,
			mLong, mLong>(myArgs.communication, &NS, storageType,
				inputBaseFile, myArgs.clusterSize, myArgs.fs, myArgs.migration, myArgs.threshold);

		mmk->setVoteToHalt(groupVoteToHalt);

		string output;
		output.append("/user/");
		output.append(myArgs.hdfsUserName.c_str());
		output.append("/m_run_output/");
		output.append(myArgs.graphName.c_str());
		mmk->setOutputPath(output.c_str());

		mmk->run(argc, argv);
		myWorkerID = mmk->getPEID();
		delete mmk;
	} else if (myArgs.algorithm == 5) { //Diameter Estimate
		groupVoteToHalt = true;
		storageType = InNbrStore;
		dimEst dE(myArgs.superSteps);

		Mizan<mLong, mLongArray, mLongArray, mLong> * mmk = new Mizan<mLong,
			mLongArray, mLongArray, mLong>(myArgs.communication, &dE, storageType, inputBaseFile, myArgs.clusterSize, myArgs.fs,
				myArgs.migration, myArgs.threshold);
		mmk->setVoteToHalt(groupVoteToHalt);

		string output;
		output.append("/user/");
		output.append(myArgs.hdfsUserName.c_str());
		output.append("/m_run_output/");
		output.append(myArgs.graphName.c_str());
		mmk->setOutputPath(output.c_str());

		mmk->run(argc, argv);
		myWorkerID = mmk->getPEID();
		delete mmk;
	}

#ifdef Verbose
	if (myWorkerID == 0) {
		std::cout << "-----TIME: Total Running Time = "
				<< float(time(NULL) - begin_time) << std::endl;
	}
#endif


	return 0;
}
