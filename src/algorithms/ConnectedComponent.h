#include "../IsuperStep.h"
#include "../Icombiner.h"
#include "../dataManager/dataStructures/data/mLong.h"

class ConnectedComponent : public IsuperStep<mLong, mLong, mLong, mLong> {
private:
  int a;
	mLong minVal(mLong A, mLong B) {
		return (A.getValue() <= B.getValue()) ? A : B;
	}

public:

	ConnectedComponent(int maxSS) {
    a = maxSS;
	}
	void initialize(userVertexObject<mLong, mLong, mLong, mLong> * data) {
		data->setVertexValue(data->getVertexID());
	}
	void compute(messageIterator<mLong> * messages,
		userVertexObject<mLong, mLong, mLong, mLong> * data,
		messageManager<mLong, mLong, mLong, mLong> * comm) {

    if (data->getCurrentSS() == 1) {
        for (int i = 0; i < data->getOutEdgeCount(); i++) {
                  comm->sendMessage(data->getOutEdgeID(i), data->getVertexID());
        }
    } else {
        mLong temp = data->getVertexValue();
        while (messages->hasNext()) {
          temp = minVal(temp, messages->getNext());
        }
//        std::cout << "temp:" << temp.getValue() << " ,data:" << data->getVertexValue().getValue() << std::endl;
        if (temp.getValue() < data->getVertexValue().getValue()) {
          data->setVertexValue(temp);
          for (int i = 0; i < data->getOutEdgeCount(); i++) {
            comm->sendMessage(data->getOutEdgeID(i), temp);
          }
        } else if (temp.getValue() > data->getVertexValue().getValue()) {
          for (int i = 0; i < data->getOutEdgeCount(); i++) {
            comm->sendMessage(data->getOutEdgeID(i), data->getVertexValue());
          }
        } else {
            data->voteToHalt();
        }
    }

	}
};
