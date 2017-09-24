#include "../IsuperStep.h"
#include "../Icombiner.h"
#include "../dataManager/dataStructures/data/mLong.h"
#include "../dataManager/dataStructures/data/mDouble.h"
#include <time.h>

/*
* Template types are <K, V1, M, A> where
*   K:  ID class
*   V1: vertex value class
*   M:  message value class
*   A:  aggregation class
*
* For SSSP, vertex and message values are both mLong
*/
class Inference : public IsuperStep<mLong, mDouble, mDouble, mDouble> {
private:
	const double threshold = 0.3;
	const double ACTIVE_VALUE = 0.3;
	int maxSuperStep;

	double fRand()
	{
		srand(time(NULL));
		double f = (double)rand() / RAND_MAX;
		return f;
	}


public:

	Inference(int maxSS) : maxSuperStep(maxSS) {}

	void initialize(userVertexObject<mLong, mDouble, mDouble, mDouble> * data) {
		
		data->setVertexValue(fRand());
		double tmp;
		for (int i = 0; i < data->getOutEdgeCount(); i++) {
			tmp = fRand() / 10;
			data->setOutEdgeValue(data->getOutEdgeID(i), mDouble(tmp));
		}
	}

	void compute(messageIterator<mDouble> * messages,
		userVertexObject<mLong, mDouble, mDouble, mDouble> * data,
		messageManager<mLong, mDouble, mDouble, mDouble> * comm) {


		if (data->getCurrentSS() == 1) {
			if (data->getVertexValue().getValue() >= ACTIVE_VALUE) {
				for (int i = 0; i < data->getOutEdgeCount(); i++) {
					comm->sendMessage(data->getOutEdgeID(i), data->getOutEdgeValue(i));
				}
			}
			data->voteToHalt();
		}

		if (data->getVertexValue().getValue() < ACTIVE_VALUE) {
			double influence = 0;
			while (messages->hasNext()) {
				influence += messages->getNext().getValue();
			}
			influence += data->getVertexValue().getValue();
			if (influence > threshold) {
				for (int i = 0; i < data->getOutEdgeCount(); i++) {
					comm->sendMessage(data->getOutEdgeID(i), data->getOutEdgeValue(i));
				}
				data->setVertexValue(influence);
			}
		}

    //if(data->getCurrentSS() > maxSuperStep)
  		data->voteToHalt();
	}
};

