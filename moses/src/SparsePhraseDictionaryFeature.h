#ifndef moses_SparsePhraseFeature_h
#define moses_SparsePhraseFeature_h

#include "FactorCollection.h"
#include "FeatureFunction.h"

namespace Moses
{

/**
  * Collection of sparse features attached to each phrase pair.
  **/
class SparsePhraseDictionaryFeature : public StatelessFeatureFunction {

public:
  SparsePhraseDictionaryFeature():
    StatelessFeatureFunction("stm", ScoreProducer::unlimited) {}

  void Evaluate(const TargetPhrase& cur_phrase,
                ScoreComponentCollection* accumulator) const;
  
  void EvaluateChart(
    const ChartHypothesis&,
    int /* featureID */,
    ScoreComponentCollection*) const {
		CHECK(0); // feature function not valid in chart decoder
	}

  // basic properties
	std::string GetScoreProducerWeightShortName(unsigned) const { return "stm"; }
	size_t GetNumInputScores() const { return 0; }

};


}


#endif