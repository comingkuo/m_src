/*
 * dynamicPartitioner.h
 *
 *  Created on: Aug 7, 2012
 *      Author: refops
 */

#ifndef DYNAMICPARTITIONER_H_
#define DYNAMICPARTITIONER_H_

#include <vector>
#include<queue>
#include "dataManager/dataManager.h"
#include "computation/computeManager.h"

template<class K, class V1, class M, class A> class computeManager;
template<class K, class V1, class M, class A> class dataManager;
template<class K, class V1, class M, class A>




class dynamicPartitioner {
	dataManager<K, V1, M, A> * dataManagerPtr;
	computeManager<K, V1, M, A> * computeManagerPtr;
	int myRank;
	double imbalanceZ;

	//kuo 20170316
	struct heap_cmp_time {//need to fix!!!!
		bool operator() (mObject<K, V1, M>* a, mObject<K, V1, M>* b) {
			if (a->getSSResTime() > b->getSSResTime())
				return a > b;
			return b > a;
		}
	};

	struct heap_cmp_comm {
		bool operator() (mObject<K, V1, M>* a, mObject<K, V1, M>* b) {
			if (a->getInTotal() > b->getInTotal())
				return a > b;
			return b > a;
		}
	};

	struct heap_cmp_commOut {
		bool operator() (mObject<K, V1, M>* a, mObject<K, V1, M>* b) {
			if (a->getOutGlobal() > b->getOutGlobal())
				return a > b;
			return b > a;
		}
	};

	struct heap_cmp_mix {
		bool operator() (mObject<K, V1, M>* a, mObject<K, V1, M>* b) {
			if ((a->getOutGlobal() + a->getInTotal()) > (b->getOutGlobal() + b->getInTotal()))
				return a > b;
			return b > a;
		}
	};
	//kuo 20170316
public:
	dynamicPartitioner(dataManager<K, V1, M, A> * inDataManagerPtr,
			computeManager<K, V1, M, A> * inComputeManagerPtr, int inMyRank, double inZ) :
			dataManagerPtr(inDataManagerPtr), computeManagerPtr(
					inComputeManagerPtr), myRank(inMyRank) {
		imbalanceZ = inZ;
	}
	bool testForImbalance(std::map<int, long long> * timeMap, double &average, int thresholdB) {
		double timeAve = 0;
		long long maxTime = 0;
		for (int i = 0; i < timeMap->size(); i++) {
			timeAve = timeAve + timeMap->at(i);
			maxTime = max(maxTime,timeMap->at(i));
		}
		timeAve = timeAve / (double) timeMap->size();
		average = timeAve;//kuo step 1 x average

		double tmp = 0;

		//kuo 標準差 start
		double timeStdv = 0;
		for (int i = 0; i < timeMap->size(); i++) {
			tmp = pow(timeAve - timeMap->at(i), 2);
			timeStdv = timeStdv + tmp;
		}
		timeStdv = timeStdv / (double) timeMap->size();
		timeStdv = sqrt(timeStdv);
		//kuo 標準差 end

		//kuo (Coefficient of Varirance) * 100%

		double zTime;
		zTime = (timeStdv / timeAve) * 100;
//		for(int i=0;i<timeMap->size(); i++)
//			cout << "kuo line56 dynamicPartitioner.h RANKi=" << i << "  time:" << timeMap->at(i) << endl;

		if (zTime > thresholdB) { // 7 is imblanceZ
			return true;
		}
		return false;
		/*double zTime = 0;
		for (int i = 0; i < timeMap->size(); i++) {
			zTime = abs(
					(double) (((double) timeMap->at(i)) - maxTime)//timeAve)
							/ (double) timeStdv);
			if (zTime > imbalanceZ) {
				return true;
			}
		}
		return false;*/
	}
	int partitionMode(std::map<int, long long> * timeMap,
			std::map<int, long long> * networkMap) {

		//processing time
		//processing network

		vector<long long> timeVec;
		double timeAve = 0;

		vector<long long> networkVec;
		double networkAve = 0;

		for (int i = 0; i < timeMap->size(); i++) {
			timeAve = timeAve + timeMap->at(i);
			timeVec.push_back(timeMap->at(i));

			/*if (myRank == 0) {
			 cout << i << "-" << timeMap->at(i) << "-" << networkMap->at(i)
			 << std::endl;
			 }*/

			networkAve = networkAve + networkMap->at(i);
			networkVec.push_back(networkMap->at(i));
		}
		sort(timeVec.begin(), timeVec.end());
		timeAve = timeAve / timeVec.size();

		sort(networkVec.begin(), networkVec.end());
		networkAve = networkAve / networkVec.size();

		/*if (myRank == 0) {
		 cout << "timeAve = " << timeAve << " networkAve " << networkAve
		 << std::endl;
		 }*/

		double timeStdv = 0;

		double networkStdv = 0;

		double tmp;
		for (int i = 0; i < timeVec.size(); i++) {
			tmp = pow(timeAve - timeVec[i], 2);
			timeStdv = timeStdv + tmp;

			tmp = pow(networkAve - networkVec[i], 2);
			networkStdv = networkStdv + tmp;
		}
		timeStdv = timeStdv / ((double) (timeVec.size()));
		timeStdv = sqrt(timeStdv);

		networkStdv = networkStdv / ((double) (networkVec.size()));
		networkStdv = sqrt(networkStdv);

		//kuo testing correlation
		double sigmaXY=0, sigmaX=0, sigmaY=0, correlation=0;
		for (int i = 0; i < timeVec.size(); i++) {
			tmp = (timeMap->at(i) - timeAve) * (networkMap->at(i) - networkAve);
			sigmaXY += tmp;
		}
		for(int i=0;i<timeVec.size();i++)
			sigmaX += pow(timeVec[i] - timeAve, 2);
		sigmaX = sqrt(sigmaX);
		for(int i=0; i<networkVec.size(); i++)
			sigmaY += pow(networkVec[i] - networkAve, 2);
		sigmaY = sqrt(sigmaY);

		if (correlation < 0)
			cout << "kuo --- negative correlation!!!\\n";

		correlation = sigmaXY/(sigmaX*sigmaY);
		//cout << "my rank:" << myRank << " , inNout correlation value:" << correlation << endl;
	if (correlation >= 0.7) { // strong linear relationship
      return 3;
    } else if (correlation >= 0.5 && correlation <0.7) {//weak relationship
      return 2;
    } else if (correlation >= 0.3 && correlation < 0.5) {
      return 1;
    } else {
      return 0;
    }
	}

	int findPEPairLong(std::map<int, long long> * sample,//kuo-20170312
			std::set<int> * ignoreList, double average, std::map<int, long long> * timeCompare) {
		int tmp;
		vector<pair<int, long long>> timeRankPair;
		for (int i = 0; i < sample->size(); i++)
			timeRankPair.push_back(make_pair(i, sample->at(i)));// 將rank及migrate objective的配對放入

		if (timeRankPair.size() == 0)
			return myRank;

		sort(timeRankPair.begin(), timeRankPair.end(), //對migrate objective 排序
			boost::bind(&std::pair<int, long long>::second, _1) <
			boost::bind(&std::pair<int, long long>::second, _2));

		//20170520 kuo
		int tmpSize = timeRankPair.size() / 2;
		int tmp;
		for (int i = 0 ; i < tmpSize; i++) { //如果objective小的內有ignoreset成員，則去除
			tmp = timeRankPair[i].first;
			if (ignoreList->find(tmp) != ignoreList->end()) {
				timeRankPair.erase(timeRankPair.begin() + i);
				i--;
				tmpSize = timeRankPair.size() / 2;
			}
		}
		sort(timeRankPair.begin(), timeRankPair.end(), //針對去除ignore set的再次排序
			boost::bind(&std::pair<int, long long>::second, _1) <
			boost::bind(&std::pair<int, long long>::second, _2));
		//20170520 kuo


		int myPairPos = myRank;
		for (int i = 0; i < timeRankPair.size(); i++) {
			if (timeRankPair[i].first == myRank)
				myPairPos = timeRankPair[
					(timeRankPair.size() - 1) - i].first; //找出overloading 與 underloading配對
		}

		if (ignoreList->find(myPairPos) == ignoreList->end()) {
			  long long temp;
			  temp = abs(timeCompare->at(myRank) - timeCompare->at(myPairPos));
			  if(temp > 1)
				  return myPairPos;
		}
		return myRank;

	}
	void findCandidatePureExecTime(float vertexZ, double timeDiff, int dst, double mean, int migrateNodes) {
		dataManagerPtr->lockDataManager();

		double sumTime = 0;
		int countNodes = 0;
		double comm = 0;


		//kuo 20170316
		vector<mObject<K, V1, M>*> vertices;
		mObject<K, V1, M> * tmpObj;
		priority_queue<mObject<K, V1, M>*, vector<mObject<K, V1, M>*>, heap_cmp_time> verticeSorted;
		for (int i = 0; i < dataManagerPtr->vertexSetSize(); i++) {
			tmpObj = dataManagerPtr->getVertexObjByPos(i);
			if (!tmpObj->isHalted())
				verticeSorted.push(tmpObj);
		}
			
		while (!verticeSorted.empty())
		{
			tmpObj = verticeSorted.top();
			if ((timeDiff - sumTime - tmpObj->getSSResTime()) > -1) {
				countNodes++;
				sumTime = sumTime + tmpObj->getSSResTime();
				tmpObj->setMigrationMark();
				dataManagerPtr->copyVertexToSoftDynamicContainer(tmpObj,dst);
			} if (sumTime > timeDiff || countNodes > migrateNodes) {
				break;
			}
			
			verticeSorted.pop();
		}
		//kuo 20170316

		std::cout << "PE" << this->myRank << " want to transfer " << countNodes
			<< " nodes with TD = " << (timeDiff - sumTime) 
			 << " to "
			<< dst << " original timeDIff = " << timeDiff << std::endl;
		dataManagerPtr->unlockDataManager();
	}
	void findCandidateMix(float vertexZ,
			long long outDiff, long long inDiff, int dst, double outMean, double inMean, double outMsgPer, int migrateNodes) {
		dataManagerPtr->lockDataManager();

		//kuo 20170316
		std::cout << "outDiff = " << outDiff << ", inDiff = " << inDiff << std::endl;
		long long sumOutMsg = 0;
		long long sumInMsg = 0;
		int countNodes = 0;
		vector<mObject<K, V1, M>*> vertices;
		mObject<K, V1, M> * tmpObj;
		priority_queue<mObject<K, V1, M>*, vector<mObject<K, V1, M>*>, heap_cmp_mix> verticeSorted;
		for (int i = 0; i < dataManagerPtr->vertexSetSize(); i++) {
			tmpObj = dataManagerPtr->getVertexObjByPos(i);
			if (!tmpObj->isHalted()) {
				verticeSorted.push(tmpObj);
      }
		}

		while (!verticeSorted.empty())
		{
			tmpObj = verticeSorted.top();
			if ((inDiff - sumInMsg - tmpObj->getInTotal()) >(-1) || (outDiff - sumOutMsg - tmpObj->getOutGlobal()) > (-1)) {
				countNodes++;
				sumOutMsg = sumOutMsg + tmpObj->getOutGlobal();
				sumInMsg = sumInMsg + tmpObj->getInTotal();
				//sumInMsg = sumInMsg + tmpObj->getInGlobal();
				tmpObj->setMigrationMark();
				dataManagerPtr->copyVertexToSoftDynamicContainer(tmpObj,
					dst);
			} if (sumOutMsg > outDiff || sumInMsg > inDiff || countNodes > migrateNodes) {
				break;
			}
			verticeSorted.pop();
		} 
		//kuo 20170316

		std::cout << "PE" << this->myRank << " want to transfer " << countNodes
				<< " nodes with TD(out + in) = "
				<< (outDiff - sumOutMsg) << "+" << (inDiff - sumInMsg) << " to "
				<< dst << " original MixDIff = " << std::endl;
		dataManagerPtr->unlockDataManager();
	}
	void findCandidateMessageInComm(float vertexZ, long long diff,
			int dst, float mean, int migrateNodes) {
		dataManagerPtr->lockDataManager();

		long long sumComm = 0;
		int countNodes = 0;
		double comm = 0;
		//kuo 20170316
		vector<mObject<K, V1, M>*> vertices;
		mObject<K, V1, M> * tmpObj;
		priority_queue<mObject<K, V1, M>*, vector<mObject<K, V1, M>*>, heap_cmp_comm> verticeSorted;
		for (int i = 0; i < dataManagerPtr->vertexSetSize(); i++) {
			tmpObj = dataManagerPtr->getVertexObjByPos(i);
			if (!tmpObj->isHalted())
				verticeSorted.push(tmpObj);
		}

		while (!verticeSorted.empty())
		{
			tmpObj = verticeSorted.top();
			if ((diff - sumComm - tmpObj->getInTotal()) > -1) {
				countNodes++;
				sumComm = sumComm + tmpObj->getInTotal();
				tmpObj->setMigrationMark();
				dataManagerPtr->copyVertexToSoftDynamicContainer(tmpObj, dst);
			} if (sumComm > diff || countNodes > migrateNodes) {
				break;
			}
			verticeSorted.pop();
		}
		//kuo 20170316

		std::cout << "PE" << this->myRank << " want to transfer " << countNodes
				<< " nodes with TD = " << (diff - sumComm) << " to " << dst
				<< " original diff = " << diff << std::endl;

		dataManagerPtr->unlockDataManager();
	}
	void findCandidateMessageOutGLDiff(float vertexZ,
			long long diff, int dst, float mean, int migrateNodes) {
		dataManagerPtr->lockDataManager();

		long long sumComm = 0;
		int countNodes = 0;
		double comm = 0;

		//kuo 20170316
		vector<mObject<K, V1, M>*> vertices;
		mObject<K, V1, M> * tmpObj;
		priority_queue<mObject<K, V1, M>*, vector<mObject<K, V1, M>*>, heap_cmp_commOut> verticeSorted;
		for (int i = 0; i < dataManagerPtr->vertexSetSize(); i++) {
			tmpObj = dataManagerPtr->getVertexObjByPos(i);
			if (!tmpObj->isHalted())
				verticeSorted.push(tmpObj);
		}

		while (!verticeSorted.empty())
		{
			tmpObj = verticeSorted.top();
			if ((diff - sumComm - tmpObj->getOutGlobal()) > -1) {
				countNodes++;
				sumComm = sumComm + tmpObj->getOutGlobal();
				tmpObj->setMigrationMark();
				dataManagerPtr->copyVertexToSoftDynamicContainer(tmpObj, dst);
			} if (sumComm > diff || countNodes > migrateNodes) {
				break;
			}
			verticeSorted.pop();
		}
		//kuo 20170316
		std::cout << "PE" << this->myRank << " want to transfer " << countNodes
				<< " nodes with TD = " << (diff - sumComm) << " to " << dst
				<< " original diff = " << diff << std::endl;

		dataManagerPtr->unlockDataManager();
	}
	bool multiGrubbsTestLong(long long &diff1, long long &diff2,
			std::map<int, long long> * sample1,//out msg diff
			std::map<int, long long> * sample2, int dst, float globalZ, double outMsgPer) {

		float mean1 = 0;
		float mean2 = 0;
		//double outMsgPer = 0.5;

		for (int i = 0; i < sample1->size(); i++) {
			mean1 = mean1 + sample1->at(i);
			mean2 = mean2 + sample2->at(i);
		}
		//cout << "PE" << myRank << " mean = " << mean << " dst = " << dst
		//<< std::endl;
		if ((mean1 == 0 && mean2==0)|| dst == myRank) {
			return false;
		}
		mean1 = mean1 / (float) sample1->size();
		mean2 = mean2 / (float) sample1->size();
		//cout << "PE"<<myRank<<" mean = " << mean << " dst = " << dst<< std::endl;

		float stdv,outStdv = 0, inStdv=0;
		float tmp;
		for (int i = 0; i < sample1->size(); i++) {
			tmp = pow(mean1 - (float) sample1->at(i), 2);
			outStdv = outStdv + tmp;
			tmp = pow(mean2 - (float)sample2->at(i), 2);
			inStdv = inStdv + tmp;
		}
		outStdv = outStdv / (float) (sample1->size());
		outStdv = sqrt(outStdv);
		inStdv = inStdv / (float)(sample1->size());
		inStdv = sqrt(inStdv);
		stdv = (outStdv * outMsgPer) + (inStdv * (1.0 - outMsgPer));
		diff1 = sample1->at(myRank) - ((sample1->at(dst)==0)?mean1:sample1->at(dst)); //mean1;
		diff2 = sample2->at(myRank) - ((sample2->at(dst)==0)?mean2:sample2->at(dst)); //mean2;

		//kuo
		double mixDiff = (double)diff1 / outMsgPer;
		mixDiff += (double)diff2 / (1 - outMsgPer);
		//kuo
		float z = mixDiff / stdv;

		//cout << "PE" << myRank << " z = " << z << std::endl;
		//<< " diff = " << diff/2 << " sample->at(myRank) = " << sample->at(myRank) << " stdv = " << stdv << endl;

		diff1 = diff1 / 2;
		diff2 = diff2 / 2;
		if (z < globalZ) {
			return false;
		} else {
			return true;
		}

	}
	bool grubbsTestDouble(double &diff, std::map<int, long long> * sample,
		int dst, float globalZ) {
		float mean = 0;
		for (int i = 0; i < sample->size(); i++) {
			mean = mean + sample->at(i);
		}
		if (mean == 0 || dst == myRank) {
			return false;
		}
		mean = mean / (float)sample->size();
		//cout << "PE"<<myRank<<" mean = " << mean << " dst = " << dst<< std::endl;

		float stdv = 0;
		float tmp;
		for (int i = 0; i < sample->size(); i++) {
			tmp = pow(mean - (float)sample->at(i), 2);
			stdv = stdv + tmp;
		}
		stdv = stdv / (float)(sample->size());
		stdv = sqrt(stdv);
		diff = sample->at(myRank) - ((sample->at(dst) == 0) ? mean : sample->at(dst)); //mean; //
		float z = diff / stdv;

		cout << "kuo(dynamicPartitioner.h line546): Timediff: "<< diff << endl;
		//cout << "PE" << myRank << " z = " << z << std::endl;
		//<< " diff = " << diff/2 << " sample->at(myRank) = " << sample->at(myRank) << " stdv = " << stdv << endl;

		diff = diff / 2;
		if (z < globalZ) {
			return false;
		}
		else {
			return true;
		}

	}
	bool grubbsTestLong(long long &diff, std::map<int, long long> * sample,
			int dst, float globalZ) {
		float mean = 0;
		for (int i = 0; i < sample->size(); i++) {
			mean = mean + sample->at(i);
		}
		if (mean == 0 || dst == myRank) {
			return false;
		}
		mean = mean / (float) sample->size();
		//cout << "PE"<<myRank<<" mean = " << mean << " dst = " << dst<< std::endl;

		float stdv = 0;
		float tmp;
		for (int i = 0; i < sample->size(); i++) {
			tmp = pow(mean - (float) sample->at(i), 2);
			stdv = stdv + tmp;
		}
		stdv = stdv / (float) (sample->size());
		stdv = sqrt(stdv);
		diff = sample->at(myRank) - ((sample->at(dst)==0)?mean:sample->at(dst)); //mean; //
		float z = diff / stdv;

		//cout << "PE" << myRank << " z = " << z << std::endl;
		//<< " diff = " << diff/2 << " sample->at(myRank) = " << sample->at(myRank) << " stdv = " << stdv << endl;

		diff = diff / 2;
		if (z < globalZ) {
			return false;
		} else {
			return true;
		}

	}
	virtual ~dynamicPartitioner() {
	}
};

#endif /* DYNAMICPARTITIONER_H_ */
