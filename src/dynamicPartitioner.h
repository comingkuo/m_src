/*
 * dynamicPartitioner.h
 *
 *  Created on: Aug 7, 2012
 *      Author: refops
 */

#ifndef DYNAMICPARTITIONER_H_
#define DYNAMICPARTITIONER_H_

#include <vector>
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
public:
	dynamicPartitioner(dataManager<K, V1, M, A> * inDataManagerPtr,
			computeManager<K, V1, M, A> * inComputeManagerPtr, int inMyRank, double inZ) :
			dataManagerPtr(inDataManagerPtr), computeManagerPtr(
					inComputeManagerPtr), myRank(inMyRank) {
		imbalanceZ = inZ;
	}
	bool testForImbalance(std::map<int, long long> * timeMap, double &average) {
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

		if (zTime > 7) { // 7 is imblanceZ
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
		timeStdv = timeStdv / ((float) (timeVec.size()));
		timeStdv = sqrt(timeStdv);

		networkStdv = networkStdv / ((float) (networkVec.size()));
		networkStdv = sqrt(networkStdv);

		//kuo testing correlation
		double sigmaXY=0, sigmaX=0, sigmaY=0, correlation=0;
		for(int i=0;i< timeVec.size();i++) 
			sigmaXY += ((timeVec[i] - timeAve) * (networkVec[i] - networkAve));
		for(int i=0;i<timeVec.size();i++)
			sigmaX += pow(timeVec[i] - timeAve, 2);
		sigmaX = sqrt(sigmaX);
		for(int i=0; i<networkVec.size(); i++)
			sigmaY += pow(networkVec[i] - networkAve, 2);
		sigmaY = sqrt(sigmaY);

		if (correlation < 0)
			cout << "kuo --- negative correlation!!!\\n";

		correlation = sigmaXY/(sigmaX*sigmaY);
		if (correlation >= 0.7) { // strong linear relationship
			return 3;
		} else if (correlation >= 0.5) {//moderate relationship
			return 2;
		} else if (correlation >=0.3) {
			return 1;
		} else {
		return 0;
		}//kuo 尚未考慮負相關
	}

	int findPEPairLong(std::map<int, long long> * sample,
			std::set<int> * ignoreList, double average) {
		vector<long long> times;
		//cout << "sample->size() = " << sample->size() << std::endl;

		for (int i = 0; i < sample->size(); i++) {

			if (sample->at(i) < average
					&& ignoreList->find(i) == ignoreList->end()) {
				times.push_back(sample->at(i));
			}
		}

		if(times.size()==0){
			return myRank;
		}
		sort(times.begin(), times.end());

		int myTime = sample->at(myRank);
		int timePos = myRank;
		int error = 0;
		for (int i = 0; i < sample->size(); i++) {
			if (times[i] == myTime) {
				timePos = ((sample->size() - 1) - i);
				//cout << "PE" << myRank << " got timePos = " << timePos << std::endl;
				break;
			}
		}
		int pairPos = myRank;
		int myTime2 = times[timePos];
		for (int i = 0; i < sample->size(); i++) {
			if (sample->at(i) == myTime2) {
				pairPos = i;
				//cout << "PE" << myRank << " got pairPos = " << pairPos << std::endl;
				break;
			}
		}

		return pairPos;
	}
	void findCandidatePureExecTime(float vertexZ, double timeDiff, int dst, double mean) {
		dataManagerPtr->lockDataManager();
		double stdv = 0;
		double tmp;
		for (int i = 0; i < dataManagerPtr->vertexSetSize(); i++) {
			mObject<K, V1, M> * tmpObj = dataManagerPtr->getVertexObjByPos(i);
			tmp = pow(mean - tmpObj->getSSResTime(), 2);
			//cout << " tmpObj->getSSResTime() = " << tmpObj->getSSResTime() << std::endl;
			stdv = stdv + tmp;
		}
		stdv = stdv / ((double)(dataManagerPtr->vertexSetSize()));
		stdv = sqrt(stdv);
		double reqZ = vertexZ;
		//cout << "PE " << myRank << " reqZ = " << reqZ << " stdv = " << stdv << " mean " << mean << std::endl;
		double vz;
		double sumTime = 0;
		int countNodes = 0;
		double comm = 0;

		for (int i = 0; i < dataManagerPtr->vertexSetSize(); i++) {
			mObject<K, V1, M> * tmpObj = dataManagerPtr->getVertexObjByPos(i);
			if (!tmpObj->isHalted()) {
				vz = (tmpObj->getSSResTime() - mean) / stdv;
				//cout << "PE " << myRank << " vz = " << vz << std::endl;
				/*comm = 0.75 * ((double) tmpObj->getMessageCountGlobal())
				- ((double) tmpObj->getMessageCountLocal());*/
				//cout << " comm = " << comm << endl;
				if (reqZ < vz && (!tmpObj->isHalted())
					&& (timeDiff - sumTime - tmpObj->getSSResTime())
				>(-1)) {
					countNodes++;
					sumTime = sumTime + tmpObj->getSSResTime();
					//sumInMsg = sumInMsg + tmpObj->getInGlobal();
					dataManagerPtr->getVertexObjByPos(i)->setMigrationMark();
					dataManagerPtr->copyVertexToSoftDynamicContainer(tmpObj,
						dst);
				}
				else {
					dataManagerPtr->getVertexObjByPos(i)->resetMigrationMark();
				}
			}
			if (sumTime > timeDiff) {
				break;
			}
		}
		for (int i = 0; i < dataManagerPtr->vertexSetSize(); i++) {
			mObject<K, V1, M> * tmpObj = dataManagerPtr->getVertexObjByPos(i);
			if (!tmpObj->isHalted() && !tmpObj->isMigrationMarked()
				&& tmpObj->getSSResTime() > 0
				&& (timeDiff - sumTime - tmpObj->getSSResTime()) > (-1)) {

				countNodes++;
				sumTime = sumTime + tmpObj->getSSResTime();

				dataManagerPtr->getVertexObjByPos(i)->setMigrationMark();
				dataManagerPtr->copyVertexToSoftDynamicContainer(tmpObj, dst);
			}

			if (sumTime > timeDiff) {
				break;
			}
		}

		std::cout << "PE" << this->myRank << " want to transfer " << countNodes
			<< " nodes with TD = " << (timeDiff - sumTime) 
			 << " to "
			<< dst << " original timeDIff = " << timeDiff << std::endl;
		dataManagerPtr->unlockDataManager();
	}
	void findCandidateMix(float vertexZ,
			long long outDiff, long long inDiff, int dst, double outMean, double inMean, double outMsgPer) {
		dataManagerPtr->lockDataManager();
		//kuo
		//double outMsgPer = 0.5;
		double outStdv = 0;
		double inStdv = 0;
		//kuo
		double tmp;

		for (int i = 0; i < dataManagerPtr->vertexSetSize(); i++) {
			mObject<K, V1, M> * tmpObj = dataManagerPtr->getVertexObjByPos(i);
			tmp = pow(outMean - tmpObj->getOutGlobal(), 2); 
			outStdv = outStdv + tmp;
			tmp = pow(inMean - tmpObj->getInTotal(), 2); 
			inStdv = inStdv + tmp;
			//cout << " tmpObj->getSSResTime() = " << tmpObj->getSSResTime() << std::endl;
		}
		inStdv = inStdv / ((double)(dataManagerPtr->vertexSetSize()));
		inStdv = sqrt(inStdv);
		outStdv = outStdv / ((double)(dataManagerPtr->vertexSetSize()));
		outStdv = sqrt(outStdv);
		double reqZ = vertexZ;
		//cout << "PE " << myRank << " reqZ = " << reqZ << " stdv = " << stdv << " mean " << mean << std::endl;
		double vz;
		//kuo
		double outVz, inVz;
		//kuo


		long long sumOutMsg = 0;
		long long sumInMsg = 0;
		int countNodes = 0;
		double comm = 0;
		for (int i = 0; i < dataManagerPtr->vertexSetSize(); i++) {
			mObject<K, V1, M> * tmpObj = dataManagerPtr->getVertexObjByPos(i);
			if (!tmpObj->isHalted()) {
				outVz = (tmpObj->getOutGlobal() - outMean) / outStdv;
				inVz = (tmpObj->getInTotal() - inMean) / inStdv;

				vz = (outVz * outMsgPer) + (inVz * (1 - outMsgPer));
				//cout << "PE " << myRank << " Mixvz = " << vz << std::endl;
				/*comm = 0.75 * ((double) tmpObj->getMessageCountGlobal())
				 - ((double) tmpObj->getMessageCountLocal());*/
				//cout << " comm = " << comm << endl;
				if (reqZ < vz && (!tmpObj->isHalted()) &&
					((inDiff - sumInMsg - tmpObj->getInTotal()) >(-1) || (outDiff - sumOutMsg - tmpObj->getOutGlobal()) > (-1))) {
					countNodes++;
					sumOutMsg = sumOutMsg + tmpObj->getOutGlobal();
					sumInMsg = sumInMsg + tmpObj->getInTotal();
					//sumInMsg = sumInMsg + tmpObj->getInGlobal();
					dataManagerPtr->getVertexObjByPos(i)->setMigrationMark();
					dataManagerPtr->copyVertexToSoftDynamicContainer(tmpObj,
							dst);
				} else {
					dataManagerPtr->getVertexObjByPos(i)->resetMigrationMark();
				}
			}
			if (sumOutMsg > outDiff || sumInMsg > inDiff) {
				break;
			}
		}
		for (int i = 0; i < dataManagerPtr->vertexSetSize(); i++) {
			mObject<K, V1, M> * tmpObj = dataManagerPtr->getVertexObjByPos(i);
			if (!tmpObj->isHalted() && !tmpObj->isMigrationMarked()
					&& (tmpObj->getInTotal() > 0 || tmpObj->getOutGlobal() > 0)
				    && ((inDiff - sumInMsg - tmpObj->getInTotal()) >(-1) || (outDiff - sumOutMsg - tmpObj->getOutGlobal()) > (-1))) {

				countNodes++;
				sumOutMsg = sumOutMsg + tmpObj->getOutGlobal();
				sumInMsg = sumInMsg + tmpObj->getInTotal();
				//sumInMsg = sumInMsg + tmpObj->getInGlobal();
				dataManagerPtr->getVertexObjByPos(i)->setMigrationMark();
				dataManagerPtr->copyVertexToSoftDynamicContainer(tmpObj, dst);
			}

			if (sumOutMsg > outDiff || sumInMsg > inDiff) {
				break;
			}
		}
		std::cout << "PE" << this->myRank << " want to transfer " << countNodes
				<< " nodes with TD(out + in) = "
				<< (outDiff - sumOutMsg) << "+" << (inDiff - sumInMsg) << " to "
				<< dst << " original MixDIff = " << std::endl;
		dataManagerPtr->unlockDataManager();
	}
	void findCandidateMessageInComm(float vertexZ, long long diff,
			int dst, float mean) {
		dataManagerPtr->lockDataManager();

		float stdv = 0;
		float tmp;
		for (int i = 0; i < dataManagerPtr->vertexSetSize(); i++) {
			mObject<K, V1, M> * tmpObj = dataManagerPtr->getVertexObjByPos(i);
			tmp = pow(mean - tmpObj->getInTotal(), 2);
			//tmp = pow(mean - tmpObj->getInGlobal(), 2);
			stdv = stdv + tmp;
		}
		stdv = stdv / ((float) (dataManagerPtr->vertexSetSize()));
		stdv = sqrt(stdv);
		float reqZ = vertexZ;
		//cout << "PE " << myRank << " reqZ = " << reqZ << " stdv = " << stdv << std::endl;
		float vz;
		long long sumComm = 0;
		int countNodes = 0;
		double comm = 0;
		for (int i = 0; i < dataManagerPtr->vertexSetSize(); i++) {
			mObject<K, V1, M> * tmpObj = dataManagerPtr->getVertexObjByPos(i);
			if (!tmpObj->isHalted()) {
				vz = (tmpObj->getInTotal() - mean) / stdv;
				//vz = (tmpObj->getInGlobal() - mean) / stdv;
				//cout << "PE " << myRank << " vz = " << vz << std::endl;
				/*	comm = 0.75 * ((double) tmpObj->getMessageCountGlobal())
				 - ((double) tmpObj->getMessageCountLocal());*/
				//cout << " comm = " << comm << endl;
				if (reqZ < vz && (!tmpObj->isHalted())
						&& (diff - sumComm - tmpObj->getInTotal()) > (-1)) {
					countNodes++;
					sumComm = sumComm + tmpObj->getInTotal();
					//sumComm = sumComm + tmpObj->getInGlobal();
					dataManagerPtr->getVertexObjByPos(i)->setMigrationMark();
					dataManagerPtr->copyVertexToSoftDynamicContainer(tmpObj,
							dst);
				} else {
					dataManagerPtr->getVertexObjByPos(i)->resetMigrationMark();
				}
			}
			if (sumComm > diff) {
				break;
			}
		}
		for (int i = 0; i < dataManagerPtr->vertexSetSize(); i++) {
			mObject<K, V1, M> * tmpObj = dataManagerPtr->getVertexObjByPos(i);
			if (!tmpObj->isHalted() && !tmpObj->isMigrationMarked()
					&& tmpObj->getInTotal() > 0
					//&& tmpObj->getInGlobal() > 0
					&& (diff - sumComm - tmpObj->getInTotal()) > (-1)) {

				countNodes++;
				sumComm = sumComm + tmpObj->getInTotal();
				//sumComm = sumComm + tmpObj->getInGlobal();
				dataManagerPtr->getVertexObjByPos(i)->setMigrationMark();
				dataManagerPtr->copyVertexToSoftDynamicContainer(tmpObj, dst);
			}

			if (sumComm > diff) {
				break;
			}
		}
		std::cout << "PE" << this->myRank << " want to transfer " << countNodes
				<< " nodes with TD = " << (diff - sumComm) << " to " << dst
				<< " original diff = " << diff << std::endl;

		dataManagerPtr->unlockDataManager();
	}
	void findCandidateMessageOutGLDiff(float vertexZ,
			long long diff, int dst, float mean) {
		dataManagerPtr->lockDataManager();

		float stdv = 0;
		float tmp;
		for (int i = 0; i < dataManagerPtr->vertexSetSize(); i++) {
			mObject<K, V1, M> * tmpObj = dataManagerPtr->getVertexObjByPos(i);
			//tmp = pow(mean - tmpObj->getOutDiff(), 2);
			tmp = pow(mean - tmpObj->getOutGlobal(), 2);
			stdv = stdv + tmp;
		}
		stdv = stdv / ((float) (dataManagerPtr->vertexSetSize()));
		stdv = sqrt(stdv);
		float reqZ = vertexZ;
		//cout << "PE " << myRank << " reqZ = " << reqZ << " stdv = " << stdv << std::endl;
		float vz;
		long long sumComm = 0;
		int countNodes = 0;
		double comm = 0;
		for (int i = 0; i < dataManagerPtr->vertexSetSize(); i++) {
			mObject<K, V1, M> * tmpObj = dataManagerPtr->getVertexObjByPos(i);
			if (!tmpObj->isHalted()) {
				vz = (tmpObj->getOutGlobal() - mean) / stdv;
				//cout << "PE " << myRank << " vz = " << vz << std::endl;
				/*	comm = 0.75 * ((double) tmpObj->getMessageCountGlobal())
				 - ((double) tmpObj->getMessageCountLocal());*/
				//cout << " comm = " << comm << endl;
				if (reqZ < vz && (!tmpObj->isHalted())
						&& (diff - sumComm - tmpObj->getOutGlobal()) > (-1)) {
					countNodes++;
					sumComm = sumComm + tmpObj->getOutGlobal();
					dataManagerPtr->getVertexObjByPos(i)->setMigrationMark();
					dataManagerPtr->copyVertexToSoftDynamicContainer(tmpObj,
							dst);
				} else {
					dataManagerPtr->getVertexObjByPos(i)->resetMigrationMark();
				}
			}
			if (sumComm > diff) {
				break;
			}
		}
		for (int i = 0; i < dataManagerPtr->vertexSetSize(); i++) {
			mObject<K, V1, M> * tmpObj = dataManagerPtr->getVertexObjByPos(i);
			if (!tmpObj->isHalted() && !tmpObj->isMigrationMarked()
					&& tmpObj->getOutGlobal() > 0
					&& (diff - sumComm - tmpObj->getOutGlobal()) > (-1)) {

				countNodes++;
				sumComm = sumComm + tmpObj->getOutGlobal();
				dataManagerPtr->getVertexObjByPos(i)->setMigrationMark();
				dataManagerPtr->copyVertexToSoftDynamicContainer(tmpObj, dst);
			}

			if (sumComm > diff) {
				break;
			}
		}
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
