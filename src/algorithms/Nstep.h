#include "../IsuperStep.h"
#include "../Icombiner.h"
#include "../dataManager/dataStructures/data/mLong.h"

#define INF      mLong(LLONG_MAX)
/*
* Template types are <K, V1, M, A> where
*   K:  ID class
*   V1: vertex value class
*   M:  message value class
*   A:  aggregation class
*
* For SSSP, vertex and message values are both mLong
*/
class Nstep : public IsuperStep<mLong, mLong, mLong, mLong> {
private:
	mLong srcID;
  int maxSuperStep;

	bool isSrc(mLong id) {
		if (id.getValue() == srcID.getValue())
			return true;
		else
			return false;
	}

	mLong min(mLong a, mLong b) {
		return (a.getValue() < b.getValue()) ? a : b;
	}


public:
	/**
	* \param srcID The vertex ID of the source.
	*/
	Nstep(mLong sourceID, int maxSS) {
    srcID = sourceID;
    maxSuperStep = maxSS;
  }

	void initialize(userVertexObject<mLong, mLong, mLong, mLong> * data) {
		// start all vertices with INF distance
		data->setVertexValue(INF);

		for (int i = 0; i < data->getOutEdgeCount(); i++) {
			data->setOutEdgeValue(data->getOutEdgeID(i), mLong(1));
		}
	}

	void compute(messageIterator<mLong> * messages,
		userVertexObject<mLong, mLong, mLong, mLong> * data,
		messageManager<mLong, mLong, mLong, mLong> * comm) {

		// can use getValue() to convert mLong to long long
		mLong currDist = data->getVertexValue();

		// potential new minimum distance
		mLong newDist = isSrc(data->getVertexID()) ? mLong(0) : INF;

		while (messages->hasNext()) {
			newDist = min(newDist, messages->getNext());
		}

		// if new distance is smaller, notify out edges
		if (newDist.getValue() < currDist.getValue()) {
			data->setVertexValue(newDist);

			for (int i = 0; i < data->getOutEdgeCount(); i++) {
				// cout << "sending msg at ss=" << data->getCurrentSS() << " to id=" << data->getOutEdgeID(i).getValue() << endl;
				// (outEdgeValue is the value of an outgoing edge)
				comm->sendMessage(data->getOutEdgeID(i),
					mLong(newDist.getValue() + 1));
			}
		}

		if (data->getCurrentSS() > maxSuperStep) {
			data->voteToHalt();
    }
	}
};
