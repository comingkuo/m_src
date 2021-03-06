/*
 * Core.h
 *
 *  Created on: Mar 28, 2012
 *      Author: refops
 *  Fix on: Jan 31, 2017
 *		Author: Kuo
 *  Description: Main function of GPSer
 */

#ifndef MIZAN_H_
#define MIZAN_H_

#include "communication/dataStructures/general.h"
#include "communication/dataStructures/dht.h"
#include "dataManager/dataStructures/data/IdataType.h"
#include "dataManager/dataStructures/data/mIntTagDouble.h"
#include "dataManager/dataStructures/data/mDoubleArray.h"
#include "dataManager/dataStructures/data/mLongArray.h"
#include "dataManager/dataStructures/data/mLong.h"
#include "dataManager/dataStructures/data/mDouble.h"
#include "dataManager/dataStructures/data/mIntCharArrayPair.h"
#include "dataManager/dataStructures/data/mKCharArrayPair.h"
#include "dataManager/dataStructures/data/mArrayIntTagNK.h"
#include "dataManager/dataManager.h"
#include "communication/commManager.h"
#include "communication/sysComm.h"
#include "communication/thread_manager.h"
#include "computation/computeManager.h"
#include "IsuperStep.h"
#include "communication/dataStructures/multi-val-map.h"
#include "Icombiner.h"
#include "IAggregator.h"
#include "computation/systemWideInfo.h"
#include "general.h"
#include "dynamicPartitioner.h"

#include <boost/function.hpp>
#include <boost/bind.hpp>

template<class K, class V1, class M, class A> class DHT;
template<class K, class V1, class M, class A> class computeManager;
template<class K, class V1, class M, class A> class dataManager;
template<class K, class V1, class M, class A> class sysComm;
template<class K, class V1, class M, class A> class userComm;
template<class K, class V1, class M, class A> class comm_manager;
template<class K, class V1, class M, class A> class dynamicPartitioner;

template<class K, class V1, class M, class A>
class Mizan {
private:
	static const int defaultThreadCnt = 4;
	thread_manager tm;

	computeManager<K, V1, M, A> * cm;
	dataManager<K, V1, M, A> * dm;
	comm_manager<K, V1, M, A> commMang;
	DHT<K, V1, M, A> dht;
	IsuperStep<K, V1, M, A> * userST;
	Icombiner<K, V1, M, A> * userCombiner;
	dynamicPartitioner<K, V1, M, A> * dp;

	std::vector<subGraphSource *> subGraphContainer;



	systemDataPointer<K, V1, M, A> dataPtr;

	//Tracking system
	int myRank;
	int PECount;
	int subGCount;
	int mySubGCount;
	int superStepCounter;
	int prop = 5400;
  
    time_t sss = 0;//kuo time
    time_t sss2 = 0;

	typedef boost::function<void(void)> fun_t;
	std::map<SYS_CMDS, fun_t> sysGroupMessageContainer;

	typedef boost::function<void(int, char *)> fun_tv;
	std::map<SYS_CMDS, fun_tv> sysGroupMessageContainerValue;

	std::map<SYS_CMDS, int> sysGroupMessageCounter;
	boost::mutex sysCommandLock;

	boost::mutex messageSyncLock;
	bool messageSyncBool;

	std::map<int, long long> peSSResTime;
	std::map<int, long long> peVertexSumTime;
	std::map<int, long long> peSSResTimeWithDHT;
	//std::map<int, long long> peCommInGlobalCnt;
	std::map<int, long long> peCommInTotalCnt;
	std::map<int, long long> peCommInXTotalCnt;
	std::map<int, long long> peCommOutGlobalCnt;
	std::map<int, long long> peCommInGlobalCnt;
	std::map<int, long long> peCommInXGlobalCnt;
	//std::map<int, long long> peCommInLocalCnt;
	//std::map<int, long long> peCommOutLocalCnt;
	std::map<int, long long> peRemMem;
	std::map<int, long long> peIOCnt;

	int inThreadCnt;

	bool dynamicPart;
	int ssInterval;
	double globalZ;
	double vertexZ;

	int storageType;
	bool groupVoteToHalt;
	enableVTH enableVertices;

  bool terminate; // kuo 20170910

	bool doWriteToDisk;
	char * outputFilePath;
	bool steelMode;
	double myZ;

	bool withCombiner;

	char ** inputPaths;
	int pathSize;
	int thresholdBB;

public:

	Mizan(communicationType commType, IsuperStep<K, V1, M, A> * inST,
			int inStorageType, char ** inputPaths, int pathSize, fileSystem fs,
			migrationMode migrate, int thresholdB) {
		thresholdBB = thresholdB;
		superStepCounter = 0;
		userST = inST;
		commMang.setCommType(commType);
		dynamicPart = true;
		//globalZ = 0.95;
		//myZ = 1.27;
		myZ = 1.96;		//zscore threshold
		globalZ = myZ;
		vertexZ = myZ;
		stealEnabled = false;
		storageType = inStorageType;
		groupVoteToHalt = false;
		doWriteToDisk = false;
		steelMode = false;
		ssOutSync = 0;
		withCombiner = false;

		inThreadCnt = 4;

		this->inputPaths = inputPaths;
		this->pathSize = pathSize;
		for (int i = 0; i < pathSize; i++) {
			addSubGraphSourcePath(fs, inputPaths[i]);
		}
		setMigration(migrate);

		init();

	}
	void setMigration(migrationMode migrate) {
		configDynamicPart();
	}
	bool isDynamicPartEnabled() {
		return dynamicPart;
	}

	bool isSteelComputeMode() {
		return steelMode;
	}
	void setOutputPath(const char * inFilePath) {
		doWriteToDisk = true;
		outputFilePath = (char *) calloc(strlen(inFilePath) + 30, sizeof(char));
		strcpy(outputFilePath, inFilePath);
	}
	void setVoteToHalt(bool value) {
		groupVoteToHalt = value;
	}
	void configDynamicPart() {
		dynamicPart = true;
	}
	void recvGraphMutation(int size, char * data) {
		//cout << "PE" << myRank << " recvGraphMutation()" << std::endl;
		//std::cout << "PE" << myRank << " is getting recv Graph mutation" << std::endl;
		mArrayIntTagNK<K> dataObj;
		dataObj.byteDecode(size, data);
		graphMutateCommands gmc =
				(graphMutateCommands) dataObj.getTag().getValue();
		cm->dataMutate(gmc, dataObj.getArray()[0], dataObj.getArray()[1]);
		//std::cout << " I am at end of recvGraphMutation" << std::endl;
	}
	void recvVertexSize(int size, char * data) {
		//cout << "PE" << myRank << " recvVertexSize()" << std::endl;
		mArrayIntTagNK<K> myData;
		myData.byteDecode(size, data);
		K * array = myData.getArray();

		if (sysGroupMessageCounter[InitVertexCount] == 1) {
			dataPtr.sysInfo.setMinVertex(array[0]);
			dataPtr.sysInfo.setMaxVertex(array[1]);
		} else {
			dataPtr.sysInfo.updateMinVertex(array[0]);
			dataPtr.sysInfo.updateMaxVertex(array[1]);
		}
		dataPtr.sysInfo.incrementCounter("_Graph_Vertex_Size",
				myData.getTag().getValue());
		if (sysGroupMessageCounter[InitVertexCount] == PECount) {
			runUserInit();
		}
	}
	long getGlobalVertexSize() {
		return dataPtr.sysInfo.getCounter("_Graph_Vertex_Size");
	}
	void recvFinishInit() {
		//cout << "PE" << myRank << " recvFinishInit()" << std::endl;
		if (sysGroupMessageCounter.find(FinishInit)->second == PECount) {
			superStepCounter++;
			dataPtr.sc->BroadcastSysMessage(StartSS, INSTANT_PRIORITY);
		}
	}
	bool stealEnabled;
	int ssActualFinish;
	void recvEndofSSTerminate(int size, char * data) {
		//cout << "PE" << myRank << " recvEndofSSTerminate()" << std::endl;
		mLongArray value;
		value.byteDecode(size, data);
		mLong * array = value.getArray();

		peSSResTime[(int) (array[0].getValue())] = array[1].getValue();
		peVertexSumTime[(int) (array[0].getValue())] = array[2].getValue();
		//peCommOutLocalCnt[(int) (array[0].getValue())] = array[3].getValue();
		//peCommOutGlobalCnt[(int) (array[0].getValue())] = array[4].getValue();

		if (array[0].getValue() == myRank) {
			time_t stop_time_ss = time(NULL);
			ssActualFinish = stop_time_ss - cm->getSSBegin();
		}

		if (stealEnabled
				&& (sysGroupMessageCounter[Terminate]
						+ sysGroupMessageCounter[EndofSS] != PECount)
				&& (peSSResTime.find(myRank) != peSSResTime.end())) {
			stealVertex();
		}

		bool iWantToTerminate = false;
		bool next = false;
		enableVertices = noneEnable;

		if (sysGroupMessageCounter[Terminate] == PECount
				&& groupVoteToHalt == true) {
			iWantToTerminate = true;
			next = true;
		} else if (sysGroupMessageCounter[Terminate] == PECount
				&& groupVoteToHalt == false) {
			next = true;
			enableVertices = EnableWithMessages;

		} else if (sysGroupMessageCounter[Terminate]
				+ sysGroupMessageCounter[EndofSS] == PECount
				&& groupVoteToHalt == true) {
			next = true;
			enableVertices = EnableAll;
		} else if (sysGroupMessageCounter[Terminate]
				+ sysGroupMessageCounter[EndofSS] == PECount
				&& groupVoteToHalt == false) {
			next = true;
			enableVertices = EnableWithMessages;
		}
		if (next && iWantToTerminate) {
			terminateMizan();
		} else if (next && !iWantToTerminate) {

			sysGroupMessageCounter[Terminate] = 0;
			sysGroupMessageCounter[EndofSS] = 0;

			if (!stealEnabled) {
				//should it wait to process all incoming messages?
				endOfSSProc();
			} else {
				dataPtr.sc->BroadcastSysMessage(StealBarrier,
						AFTER_DATABUFFER_PRIORITY);
			}
		}
	}
	void terminateMizan() {
    if(groupVoteToHalt == false) {
		 //   cout << myRank << "  size!!!!" << peCommInTotalCnt.empty() << endl;
     for(int i = 0; i < peSSResTimeWithDHT.size(); i++) {
      //if(!peCommInTotalCnt.empty() || !peCommOutGlobalCnt.empty() || !peCommInXTotalCnt.empty()) {
      if(peCommInTotalCnt[i] != 0 || !peCommOutGlobalCnt[i] != 0|| !peCommInXTotalCnt[i]!= 0) {
         terminate = false;
		    //cout << myRank << "  terminating false!!!!!!!!!!!!!!!!!!" << endl;
      
      }
     }
    }
	//	cout << myRank << "KKK is terminating.."<< groupVoteToHalt  <<"  ter" << terminate<< endl;
    if(!terminate) {
		if (doWriteToDisk) {
			strcat(outputFilePath, "/part_");
			mInt rank(myRank);
			strcat(outputFilePath, rank.toString().c_str());
			dm->writeToDisk(outputFilePath);
		}
		dataPtr.sc->sendSelfExit();
		if (myRank == 0) {
			time_t appFinish = time(NULL);
			std::cout << "-----TIME: Total Running Time without IO = "
					<< (appFinish - cm->getAppStart()) << std::endl;

		}
    }
	}
	void recvStealBarrier() {
		if (sysGroupMessageCounter[StealBarrier] == PECount) {
			sysGroupMessageCounter[StealBarrier] = 0;
			endOfSSProc();
		}
	}
	void recvMigrateBarrier() {
		if (sysGroupMessageCounter[MigrateBarrier] == PECount) {
			sysGroupMessageCounter[MigrateBarrier] = 0;
			dataPtr.sc->BroadcastSysMessage(StartSS, AFTER_DATABUFFER_PRIORITY);
		}
	}

	void endOfSSProc() {
		messageSyncLock.lock();
		messageSyncBool = false;
		messageSyncLock.unlock();

		cm->gatherStats();
		cm->performDataMutations();
		bool subTerminate = cm->cleanUp(enableVertices);

 

    int subTerRank = (int)subTerminate;//kuo 20170923
    int rank_size = 10;
    int *subTers = new int[rank_size]();
    //MPI_Gather(&subTerRank, 1, MPI_INT, subTers, 1, MPI_INT, 0,MPI_COMM_WORLD);
    commMang.gatherMsg(subTers, subTerRank, 0);
    subTerRank=1;
    if (myRank == 0) {
      for (int i = 0; i < rank_size; i++) {
        if (subTers[i] == 0) {//if subTerminate = false
            subTerRank = 0;
            break;
        }
      }
    }
    commMang.sendBMsg(&subTerRank, 0);//kuo 20170923
    subTerminate = (bool) subTerRank;
    delete[] subTers;







		if (subTerminate == true && superStepCounter != 1) {
			terminateMizan();
		} else {
			mLong * array = new mLong[8];
			array[0].setValue(myRank);
			array[1].setValue(ssActualFinish);
			array[2].setValue(
					cm->getXSumInCommLocal() + cm->getXSumInCommGlobal());
			array[3].setValue(
					cm->getSumInCommLocal() + cm->getSumInCommGlobal());
			array[4].setValue(cm->getSumOutCommGlobal());
			array[5].setValue(cm->getSumInCommGlobal());
			array[6].setValue(cm->getXSumInCommGlobal());
			array[7].setValue(
					((long) (dm->getAvaliableSystemMemoryPercent() * 100)));
			mLongArray value(8, array);

      terminate = false; //kuo 20170910
      if(groupVoteToHalt == false) {// kuo 20170910
        for(int i = 0; i < peSSResTimeWithDHT.size(); i++) {
          if(peCommInTotalCnt[i] != 0 || !peCommOutGlobalCnt[i] != 0|| !peCommInXTotalCnt[i]!= 0) {
            terminate = true;
           }
        }
      }

			cm->storeInCommInfo();

			dataPtr.sc->BroadcastSysMessageValue(LateStats, value,
					AFTER_DATABUFFER_PRIORITY);
		}

	}
	void recvLateStats(int size, char * data) {
		//cout << "PE" << myRank << " recvLateStats()" << std::endl;
		mLongArray value;
		value.byteDecode(size, data);
		mLong * array = value.getArray();

		peSSResTimeWithDHT[(int) (array[0].getValue())] = array[1].getValue();
		peCommInXTotalCnt[(int) (array[0].getValue())] = array[2].getValue();
		peCommInTotalCnt[(int) (array[0].getValue())] = array[3].getValue();
		peCommOutGlobalCnt[(int) (array[0].getValue())] = array[4].getValue();
		peCommInGlobalCnt[(int) (array[0].getValue())] = array[5].getValue();
		peCommInXGlobalCnt[(int) (array[0].getValue())] = array[6].getValue();
		peRemMem[(int) (array[0].getValue())] = array[7].getValue();

		if (myRank == 0) {
			int theRank = (int) (array[0].getValue());
		  std::cout << " myRank  " <<myRank << "  " 				
			  << " -----Messages: Actual Finish = "
			  << peSSResTimeWithDHT[theRank] << " Global in comm = "
				<< peCommInXTotalCnt[theRank] << "/"
				<< peCommInTotalCnt[theRank] << " Global out Comm = "
				<< peCommOutGlobalCnt[theRank] << " Memory Rem = "
				<< peRemMem[theRank] << std::endl;
			cout.flush();
		}

		if (sysGroupMessageCounter[LateStats] == PECount) {
			sysGroupMessageCounter[LateStats] = 0;

			(tm.tp).schedule(
					boost::bind(&Mizan<K, V1, M, A>::bodyRecvEndOfSSTerminate,
							this));
		}
	}
	void bodyRecvEndOfSSTerminate() {
		//cout << "PE" << myRank << " bodyRecvEndOfSSTerminate()" << std::endl;
		//process dynamic partitioning

		long long maxTime = 0;
		long long minMem = 100;
		int maxTimeRank;
		int minMemRank;
		if (myRank == 0) {

			for (int i = 0; i < peSSResTimeWithDHT.size(); i++) {
				if (peSSResTimeWithDHT[i] > maxTime) {
					maxTime = peSSResTimeWithDHT[i];
					maxTimeRank = i;
				}
				if (peRemMem[i] < minMem) {
					minMem = peRemMem[i];
					minMemRank = i;
				}
			}

			std::cout << "PE" << maxTimeRank
					<< " -----Messages: Actual Finish = " << maxTime << " -----"
					<< std::endl;
			std::cout << "PE" << minMemRank << " -----Messages: Min Mem = "
					<< minMem << " -----" << std::endl;
		}


		if (dynamicPart) {
			//std::cout << "dynamicPart" << std::endl;
			//Move the original content of Soft Vertexs
			bool softStatus = moveSoftVertex();
			if (!softStatus) {
				dynamicPartition();
			}
		}
		steelMode = false;
		peSSResTime.clear();
		peVertexSumTime.clear();
		peSSResTimeWithDHT.clear();
		peCommInTotalCnt.clear();
		peCommInXTotalCnt.clear();
		peCommInGlobalCnt.clear();
		peCommInXGlobalCnt.clear();
		//peCommInLocalCnt.clear();
		//peCommInGlobalCnt.clear();
		peCommOutGlobalCnt.clear();
		//peCommOutLocalCnt.clear();
		peRemMem.clear();
		peIOCnt.clear();
		superStepCounter++;
		/*cout << "PE"<<myRank<<" is sending StartSS" << endl;
		 */

		if (dynamicPart) {
			dataPtr.sc->BroadcastSysMessage(MigrateBarrier,
					AFTER_DATABUFFER_PRIORITY);
		} else {
			dataPtr.sc->BroadcastSysMessage(StartSS, AFTER_DATABUFFER_PRIORITY);
		}
	}
	void recvStartSS() {
		//cout << "PE" << myRank << " recvStartSS()" << std::endl;
		//cout << "PE" << myRank << " [StartSS] = "
		//<< sysGroupMessageCounter[StartSS] << endl;
		//cout << "PE" << sc->getRank() << " Got out of SYNC message from " << dst.getValue() << ":" << srcRank << ":"<<*superStepCnt<< endl;
		if (sysGroupMessageCounter[StartSS] == 1 && cm->getCnt() > 0) {
			cout << "----- ERROR ----- PE" << myRank
					<< " Got out of SYNC forward messages = " << cm->getCnt()
					<< endl;
			//std::set<K> * bla = cm->getbla();

			/*typename std::set<K>::iterator it;
			 for (it = bla->begin(); it != bla->end(); it++) {
			 std::cout << "PE" << myRank << " wrong MessageID " << (mLong(*it)).getValue()
			 << " is recvd " << std::endl;
			 }*/
			cout.flush();
		}
		cm->setPendingMessageSync();

		messageSyncLock.lock();
		messageSyncBool = true;

		if (ssOutSync > 0) {
			std::cout << "---- ERROR ---- PE" << myRank
					<< " I got an out of Sync message between ENDSS and StartSS cnt = "
					<< ssOutSync << std::endl;
			cout.flush();
		}

		ssOutSync = 0;
		messageSyncLock.unlock();

		if (sysGroupMessageCounter[StartSS] == PECount) {
			sysGroupMessageCounter[StartSS] = 0;
			runSingleSuperStep();
		}
	}

	void dynamicPartition() {

		//cout << "PE" << myRank << " dynamicPartition()" << std::endl;
		//float meanTime = cm->getVerAveResTime();
		//float maxTime = cm->getVerMaxResTime();
		//int dstTime = dp->findPEPair(&peResTime);
		//bool myTestTime = dp->grubbsTest(timeDIff, &peResTime, dstTime,globalZ);

		//skip migration counting SS1


		double average = 0;
		int migrateNodes;
		bool migrationTest = dp->testForImbalance(&this->peSSResTimeWithDHT,
				average, thresholdBB);


		//remove extra used memory
		std::set<int> ignoreSet;
		for (int i = 0; i < peSSResTimeWithDHT.size(); i++) {
			if (this->peRemMem[i] < 6) {
				ignoreSet.insert(i);
			}
		}

		if (migrationTest) {

			int outMsgScore = dp->partitionMode(&this->peSSResTimeWithDHT,
					&peCommOutGlobalCnt);
			int inMsgScore = dp->partitionMode(&this->peSSResTimeWithDHT,
					&this->peCommInTotalCnt);// kuo modified 12/11/16


			if (outMsgScore > inMsgScore) {
				//Migrate for outMessages
				if (myRank == 0) {
					cout << "PE" << myRank << " outNetwork based migration "
							<< std::endl;
				}

				long long messageDiff = 0;
				float meanMessage = cm->getAveOutCommGlobal();

				int dstNetwork = dp->findPEPairLong(&peCommOutGlobalCnt,
						&ignoreSet, average, &peSSResTimeWithDHT);

				migrateNodes = (this->peSSResTimeWithDHT.at(myRank) - this->peSSResTimeWithDHT.at(dstNetwork)) * prop;
				bool myTestNetwork = dp->grubbsTestLong(messageDiff,
						&peCommOutGlobalCnt, dstNetwork, globalZ);
				if (myTestNetwork) {
					dp->findCandidateMessageOutGLDiff(vertexZ,
							messageDiff, dstNetwork, meanMessage, migrateNodes);
					hardMigrate(dstNetwork);
				}
			} else if (inMsgScore > outMsgScore) {
				//Migrate InMessages
				if (myRank == 8) {
					cout << "PE" << myRank << " InNetwork based migration "
							<< std::endl;
				}

				long long messageDiff = 0;
				float meanMessage = (float) peCommInTotalCnt[myRank] //(float) peCommInGlobalCnt[myRank]
				/ (float) dm->vertexSetSize();

				//int dstNetwork = dp->findPEPairLong(&peCommInGlobalCnt);
				int dstNetwork = dp->findPEPairLong(&peCommInTotalCnt,
						&ignoreSet, average, &peSSResTimeWithDHT);
					 migrateNodes = (this->peSSResTimeWithDHT.at(myRank) - this->peSSResTimeWithDHT.at(dstNetwork)) * prop;
     
				bool myTestNetwork = dp->grubbsTestLong(messageDiff,
						&peCommInTotalCnt, dstNetwork, globalZ);

				if (myTestNetwork) {
					dp->findCandidateMessageInComm(vertexZ,
							messageDiff, dstNetwork, meanMessage, migrateNodes);
					hardMigrate(dstNetwork);
				}
			} else if (outMsgScore == 0 && outMsgScore == inMsgScore) {
				if (myRank == 0) {
					cout << "PE" << myRank
							<< " vertex response based migration " << std::endl;
				}

				double timeDiff = 0;
				double meanMessage = cm->getVerAveResTime();

				int dstTime = dp->findPEPairLong(&peSSResTimeWithDHT,
						&ignoreSet, average, &peSSResTimeWithDHT);
        migrateNodes = (this->peSSResTimeWithDHT.at(myRank) - this->peSSResTimeWithDHT.at(dstTime)) * prop;
				bool myTestTime = dp->grubbsTestDouble(timeDiff,
					&peSSResTimeWithDHT, dstTime, globalZ);

				if (myTestTime) {
					dp->findCandidatePureExecTime(vertexZ, timeDiff,
							dstTime, meanMessage, migrateNodes);
					hardMigrate(dstTime);
				}
			} else { //add by kuo, mix outgoing and incoming
				if (myRank == 0) {
					cout << "PE" << myRank
						<< " vertex response based migration " << std::endl;
				}


				long long outMsgDiff = 0;
				long long inMsgDiff = 0;
				double timeDiff = 0;
				double meanMessage = cm->getVerAveResTime();
				double outMeanMessage = (double) cm->getAveOutCommGlobal(); 
				double inMeanMessage = (double)peCommInTotalCnt[myRank];
				double outMsgPer = 0.5;
				//kuo 

				int dstTime = dp->findPEPairLong(&peSSResTimeWithDHT,
					&ignoreSet, average, &peSSResTimeWithDHT);
        migrateNodes = (this->peSSResTimeWithDHT.at(myRank) - this->peSSResTimeWithDHT.at(dstTime)) * prop;


				bool myTestMix = dp->multiGrubbsTestLong(outMsgDiff,
					inMsgDiff, &peCommOutGlobalCnt,
					&peCommInTotalCnt, dstTime, globalZ, outMsgPer);


				if (myTestMix) {
					dp->findCandidateMix(vertexZ,
						outMsgDiff, inMsgDiff, dstTime, outMeanMessage, inMeanMessage, outMsgPer, migrateNodes);
					hardMigrate(dstTime);
				}
			}

		}

	}
	bool moveSoftVertex() {

		int globalSize = 0;
		if (dm->getSoftDynamicVertexSendSize() == 0) {
			return false;
		} else {
			for (int i = 0; i < dm->getSoftDynamicVertexSendSize(); i++) {
				if (i == 0) {
					cout << "PE" << myRank
							<< " entering VertexMigrate() of size = "
							<< dm->getSoftDynamicVertexSendSize() << endl;
				}
				mObject<K, V1, M> * kkk = dm->getSoftDynamicVertexSendObjByPos(
						i);
				int byteSize = 0;
				int dst = dm->getSoftDynamicVertexSendDst(kkk->getVertexID());
				char * vertex;
				vertex = kkk->byteEncode(byteSize, false);
				mCharArray verData(byteSize, vertex);

				dataPtr.sc->sendSysMessageValue(VertexMigrate, verData,
						INSTANT_PRIORITY, dst);

				globalSize = globalSize + byteSize;

			}

			dm->clearSoftDynamicVertexSend();
			cout << "PE" << myRank << " softMigration+1 move size = "
					<< globalSize << std::endl;
		}
		return true;
	}
	void recvSoftVertex(int size, char * data) {


		mObject<K, V1, M> * newVertex = new mObject<K, V1, M>();
		newVertex->byteDecode(size, data);
		newVertex->setMigrationMark();
		K vertexID = newVertex->getVertexID();
		if (!vertexExists(vertexID)) {
			addNewVertex(newVertex);
		}
	}
	void recvVertexMigrate(int size, char * data) {

		mObject<K, V1, M> * newVertex = new mObject<K, V1, M>();

		newVertex->byteDecode(size, data);

		K vertex = newVertex->getVertexID();


		if (dm->vertexExists(vertex)) {

			dm->applyUpdate(newVertex);
		} else {

			addNewVertex(newVertex);

		}
	}

	void softMigrate(int dst) {
		migrateFail.clear();
		int globalSize = 0;
		for (int i = 0; i < dm->getSoftDynamicVertexSendSize(); i++) {
			int size;
			mObject<K, V1, M> * transVertex =
					dm->getSoftDynamicVertexSendObjByPos(i);
			char * rawVertex = transVertex->byteEncode(size, false);
			mCharArray vertexContainer(size, rawVertex);
			bool msgSent = dataPtr.sc->sendSysMessageValue(SendSoftVertex,
					vertexContainer, AFTER_DATABUFFER_PRIORITY, dst);

			if (!msgSent) {
				transVertex->resetMigrationMark();
				migrateFail.insert(transVertex->getVertexID());

			} else {
				globalSize = globalSize + size;
			}

			globalSize = globalSize + size;
		}
		cout << "PE" << myRank << " softMigration move size = " << globalSize
				<< std::endl;
		dm->derigsterAllNBRSoftDynamic(&migrateFail);
		dm->removeMovedVertexs(&migrateFail);
	}
	std::set<K> migrateFail;
	void hardMigrate(int dst) {
		const clock_t start_hard_Migrate = clock();
		migrateFail.clear();
		int globalSize = 0;
		for (int i = 0; i < dm->getSoftDynamicVertexSendSize(); i++) {
			mObject<K, V1, M> * transVertex =
					dm->getSoftDynamicVertexSendObjByPos(i);

			int size = 0;
			char * rawVertex = transVertex->byteEncode(size, true);
			mCharArray vertexContainer(size, rawVertex);
			bool msgSent = dataPtr.sc->sendSysMessageValue(VertexMigrate,
					vertexContainer, AFTER_DATABUFFER_PRIORITY, dst);

			if (!msgSent) {
				transVertex->resetMigrationMark();
				migrateFail.insert(transVertex->getVertexID());

			} else {
				globalSize = globalSize + size;
			}
		}


		dm->derigsterAllNBRSoftDynamic(&migrateFail);
		dm->removeMovedVertexs(&migrateFail);
		dm->clearSoftDynamicVertexSend();
		const clock_t stop_hard_Migrate = clock();

		cout << "PE" << myRank << " hardMigration move size = " << globalSize
				<< " time = "
				<< ((double) (stop_hard_Migrate - start_hard_Migrate))
						/ ((double) CLOCKS_PER_SEC) << std::endl;
	}

	void init() {

		sysGroupMessageContainerValue[InitVertexCount] = boost::bind(
				&Mizan<K, V1, M, A>::recvVertexSize, this, _1, _2);
		sysGroupMessageContainer[FinishInit] = boost::bind(
				&Mizan<K, V1, M, A>::recvFinishInit, this);
		sysGroupMessageContainerValue[EndofSS] = boost::bind(
				&Mizan<K, V1, M, A>::recvEndofSSTerminate, this, _1, _2);
		sysGroupMessageContainerValue[Terminate] = boost::bind(
				&Mizan<K, V1, M, A>::recvEndofSSTerminate, this, _1, _2);
		sysGroupMessageContainer[StartSS] = boost::bind(
				&Mizan<K, V1, M, A>::recvStartSS, this);
		sysGroupMessageContainerValue[SendSoftVertex] = boost::bind(
				&Mizan<K, V1, M, A>::recvSoftVertex, this, _1, _2);
		sysGroupMessageContainerValue[VertexMigrate] = boost::bind(
				&Mizan<K, V1, M, A>::recvVertexMigrate, this, _1, _2);
		sysGroupMessageContainerValue[StealVertex] = boost::bind(
				&Mizan<K, V1, M, A>::processStealVertex, this, _1, _2);
		sysGroupMessageContainerValue[SendStolenVertex] = boost::bind(
				&Mizan<K, V1, M, A>::processStolenVertex, this, _1, _2);
		sysGroupMessageContainerValue[StolenVertexResult] = boost::bind(
				&Mizan<K, V1, M, A>::processStolenVertexResult, this, _1, _2);
		sysGroupMessageContainerValue[GraphMutation] = boost::bind(
				&Mizan<K, V1, M, A>::recvGraphMutation, this, _1, _2);
		sysGroupMessageContainerValue[LateStats] = boost::bind(
				&Mizan<K, V1, M, A>::recvLateStats, this, _1, _2);
		sysGroupMessageContainer[StealBarrier] = boost::bind(
				&Mizan<K, V1, M, A>::recvStealBarrier, this);
		sysGroupMessageContainerValue[Aggregator] = boost::bind(
				&Mizan<K, V1, M, A>::recvAggregator, this, _1, _2);
		sysGroupMessageContainer[MigrateBarrier] = boost::bind(
				&Mizan<K, V1, M, A>::recvMigrateBarrier, this);

		sysGroupMessageCounter.insert(make_pair(InitVertexCount, 0));
		sysGroupMessageCounter.insert(make_pair(FinishInit, 0));
		sysGroupMessageCounter.insert(make_pair(EndofSS, 0));
		sysGroupMessageCounter.insert(make_pair(Terminate, 0));
		sysGroupMessageCounter.insert(make_pair(StartSS, 0));
		sysGroupMessageCounter.insert(make_pair(SendSoftVertex, 0));
		sysGroupMessageCounter.insert(make_pair(VertexMigrate, 0));
		sysGroupMessageCounter.insert(make_pair(StealVertex, 0));
		sysGroupMessageCounter.insert(make_pair(SendStolenVertex, 0));
		sysGroupMessageCounter.insert(make_pair(StolenVertexResult, 0));
		sysGroupMessageCounter.insert(make_pair(GraphMutation, 0));
		sysGroupMessageCounter.insert(make_pair(LateStats, 0));
		sysGroupMessageCounter.insert(make_pair(StealBarrier, 0));
		sysGroupMessageCounter.insert(make_pair(Aggregator, 0));

		dataPtr.sysInfo.addCounter("_Graph_Vertex_Size", 0);
	}
	void recvAggregator(int size, char * data) {
		mKCharArrayPair<A> value;
		value.byteDecode(size, data);
		char * ptr = value.getValue().getValue();

		dataPtr.aggContainerLock.lock();

		typename std::map<char *, IAggregator<A> *>::iterator it;
		for (it = dataPtr.aggContainer.begin();
				it != dataPtr.aggContainer.end(); it++) {
			int length = strlen((*it).first);
			char * aggKey = (char*) calloc(length + 1, sizeof(char));
			strncpy(aggKey, (*it).first, length);
			if (strcmp(ptr, aggKey) == 0) {
				IAggregator<A> * usrAggClass = (*it).second;
				usrAggClass->aggregate(value.getTag());
			}

		}
		dataPtr.aggContainerLock.unlock();
	}
	void processSysCommand(SYS_CMDS message) {
		sysGroupMessageCounter[message]++;
		sysGroupMessageContainer.find(message)->second();
	}
	void processSysCommand(SYS_CMDS message, int size, char * data) {

		sysGroupMessageCounter[message]++;
		sysGroupMessageContainerValue.find(message)->second(size, data);
	}

	void run(int argc, char** argv) {

		tm.threadpool_init(inThreadCnt);

		commMang.init(argc, argv);

		myRank = commMang.getRank();

		PECount = commMang.getPsize();

		subGCount = PECount;
		mySubGCount = 1;

		dataPtr.sc = new sysComm<K, V1, M, A>(&commMang);

		dataPtr.uc = new userComm<K, V1, M, A>(&commMang);

		//sc = new sysComm<K, V, A>(&commMang);		dht.MapSysComm(dataPtr.sc);

		if (myRank == 0) {
			sss2 = time(NULL);
		}


		dm = new dataManager<K, V1, M, A>(subGraphContainer[myRank]->filePath,
				subGraphContainer[myRank]->fileSystemType, storageType);
		cm = new computeManager<K, V1, M, A>(dm, dataPtr.sc, dataPtr.uc, userST,
				&dataPtr.sysInfo, &superStepCounter, &(dataPtr.aggContainer),
				&(dataPtr.aggContainerLock), groupVoteToHalt, withCombiner);
		cm->setCombiner(userCombiner);
		dm->setComputeRank(myRank);

		dp = new dynamicPartitioner<K, V1, M, A>(dm, cm, myRank, myZ);

		commMang.MapInstances(this, &dht, &tm);

		(tm.tp).schedule(
				boost::bind(&comm_manager<K, V1, M, A>::listen_4cmds,
						&commMang));
		(tm.tp).schedule(
				boost::bind(&comm_manager<K, V1, M, A>::process_queueCmds,
						&commMang));
		(tm.tp).schedule(boost::bind(&computeManager<K, V1, M, A>::init, cm));
		wait();
	}
	void addSubGraphSourcePath(fileSystem inType, char * inputPath) {
		subGraphSource * subGSrc = new subGraphSource;
		subGSrc->fileSystemType = inType;
		subGSrc->filePath = (char *) calloc(strlen(inputPath) + 1,
				sizeof(char));
		strcpy(subGSrc->filePath, inputPath);
		subGraphContainer.push_back(subGSrc);
	}
	void runUserInit() {
		//superStep
		if (myRank == 0) {
			std::cout << "Starting user init.." << endl;
		}
		(tm.tp).schedule(
				boost::bind(&computeManager<K, V1, M, A>::userInit, cm));
	}
	void runSingleSuperStep() {
		//cout << "PE" << myRank << " runSingleSuperStep " << std::endl;
		//superStep
    float syzz = 0.0;
    if(superStepCounter %2 ==0 ) {
      sss2 = time(NULL);
      syzz = ((float) (sss2 - sss));
    } else {
      sss = time(NULL);
      syzz = ((float) (sss - sss2));
    }
		if (myRank == 0) {
			std::cout << " ---------- Starting Superstep " << superStepCounter
					<< " ----------Time:" << syzz << endl;
			cout.flush();
		}

		(tm.tp).schedule(
				boost::bind(&computeManager<K, V1, M, A>::superStep, cm));
		//superStepCounter++;
	}
	void registerMessageCombiner(Icombiner<K, V1, M, A> * combiner) {
		withCombiner = true;
		userCombiner = combiner;
	}
	void registerAggregator(char * aggName, IAggregator<A> * aggImp) {
		dataPtr.aggContainer[aggName] = aggImp;
	}
	void wait() {
		(tm.tp).wait();
	}
	computeStatus getStatus() {
	}

	bool vertexExists(K& vertex) {
		return dm->vertexExists(vertex);
	}
	void addNewVertex(K& vertex) {
		dataPtr.sc->InformOwnership(vertex);
		dm->addNewVertex(vertex);
	}
	void addNewVertex(mObject<K, V1, M> * verObj) {
		dataPtr.sc->InformOwnership(verObj->getVertexID());
		dm->addNewVertex(verObj);
	}

	int getPEID() {
		return myRank;
	}
	bool acceptGate() {
		bool answer;
		messageSyncLock.lock();

		answer = messageSyncBool;

		messageSyncLock.unlock();

		return answer;
	}
	void appendIncomeQueueAll(M& message) {
		if (acceptGate()) {
			cm->appendIncomeQueueAll(message);
		}
	}
	void appendIncomeQueueNbr(K& src, M& message, DATA_CMDS inOut) {
		if (acceptGate()) {
			cm->appendIncomeQueueNbr(src, message, inOut);
		}
	}
	int ssOutSync;
	void appendIncomeQueue(K& dst, M& message, int srcRank) {
		if (acceptGate()) {
			cm->appendIncomeQueue(dst, message, srcRank);
		} else {
			ssOutSync++;
		}
	}
	void appendLocalIncomeQueue(K& dst, M& message) {
		if (acceptGate()) {
			cm->appendLocalIncomeQueue(dst, message);
		}
	}
	void enableVertexSteal() {
		stealEnabled = true;
	}
	void stealVertex() {
		if (peSSResTime.size() != PECount || !peSSResTime.empty()) {
			int i = 0;
			for (; i < PECount; i++) {
				if (peSSResTime.find(i) == peSSResTime.end()) {
					/*cout << "PE" << myRank << " sending steal request to " << i
					 << endl;*/
					mInt id(myRank);
					dataPtr.sc->sendSysMessageValue(StealVertex, id,
							INSTANT_PRIORITY, i);
					//break;
				}
			}
		}
	}

	void processStealVertex(int size, char * data) {

		mInt id;
		id.byteDecode(size, data);
		mObject<K, V1, M> * vertexObj;

		bool successfulSend = false;

		while (!successfulSend && vertexObj != 0) {
			vertexObj = dm->iterateData(true);

			while (vertexObj != 0 && vertexObj->isMigrationMarked()) {
				dm->deleteFromStolen(vertexObj->getVertexID());
				vertexObj = dm->iterateData(true);
			}

			if (vertexObj != 0) {
				int outSize = 0;
				char * vertexChar = vertexObj->byteEncode(outSize, true);
				mIntCharArrayPair outVertexChar(myRank, outSize, vertexChar);

				successfulSend = dataPtr.sc->sendSysMessageValue(
						SendStolenVertex, outVertexChar, INSTANT_PRIORITY,
						id.getValue());

				free(vertexChar);

				if (!successfulSend) {
					cm->appendToRejectedMessageQueue(vertexObj);
					dm->deleteFromStolen(vertexObj->getVertexID());
				}
			}
		}
		//StealVertexBack

	}
	void processStolenVertex(int size, char * data) {
		//StolenVertexResult,
		steelMode = true;
		mIntCharArrayPair stolenVertexStorage; //calloc
		stolenVertexStorage.byteDecode(size, data);

		mObject<K, V1, M> stolenVertex;
		//const clock_t begin_time = clock();
		stolenVertex.byteDecode(stolenVertexStorage.getValue().getSize(),
				stolenVertexStorage.getValue().getValue());
		//float decodeTime = float(clock() - begin_time) / CLOCKS_PER_SEC;

		cm->singleVertexStep(&stolenVertex);
		int outSize = 0;
		char * vertexChar = stolenVertex.byteEncodeCompact(outSize); //calloc
		mCharArray outVertexChar(outSize, vertexChar); //calloc and encode calloc

		dataPtr.sc->sendSysMessageValue(StolenVertexResult, outVertexChar,
				INSTANT_PRIORITY, stolenVertexStorage.getTag().getValue());
		stealVertex();
	}
	void processStolenVertexResult(int size, char * data) {
		mObject<K, V1, M> * stolenVertex = new mObject<K, V1, M>();
		stolenVertex->byteDecodeCompact(size, data);
		dm->deleteFromStolen(stolenVertex->getVertexID());
		dm->applyUpdate(stolenVertex);
		/*cout << "PE" << myRank << " getting stolen vertex result "
		 << stolenVertex.getVertexID().getValue() << endl;*/
	}
	std::vector<edge<K, V1> *> * getInEdges(K vertex) {
		return dm->getInEdges(vertex);
	}
	std::vector<edge<K, V1> *> * getOutEdges(K vertex) {
		return dm->getOutEdges(vertex);
	}

	virtual ~Mizan() {
		delete (dp);
		delete (dataPtr.sc);
		delete (dataPtr.uc);
		delete (cm);
		delete (dm);
		if (doWriteToDisk) {
			free(outputFilePath);
		}
	}
};

#endif