#include "util/check.hh"
#include <stdexcept>
#include <iostream>


#include "moses/ChartManager.h"
#include "moses/Hypothesis.h"
#include "moses/Manager.h"
#include "moses/StaticData.h"
#include "moses/TranslationModel/PhraseDictionaryDynSuffixArray.h"
#include "moses/TranslationSystem.h"
#include "moses/TreeInput.h"
#include "moses/LMList.h"
#include "moses/LM/ORLM.h"

#include <xmlrpc-c/base.hpp>
#include <xmlrpc-c/registry.hpp>
#include <xmlrpc-c/server_abyss.hpp>

using namespace Moses;
using namespace std;

typedef std::map<std::string, xmlrpc_c::value> params_t;

/** Find out which translation system to use */
const TranslationSystem& getTranslationSystem(params_t params)
{
  string system_id = TranslationSystem::DEFAULT;
  params_t::const_iterator pi = params.find("system");
  if (pi != params.end()) {
    system_id = xmlrpc_c::value_string(pi->second);
  }
  VERBOSE(1, "Using translation system " << system_id << endl;)
  return StaticData::Instance().GetTranslationSystem(system_id);
}

class DiscoveryManager: public xmlrpc_c::method {
public:
  DiscoveryManager() {
      // signature and help strings are documentation -- the client
      // can query this information with a system.methodSignature and
      // system.methodHelp RPC.
      this->_signature = "S:S";
      this->_help = "Gives information about this server instance";
  }
  void execute(xmlrpc_c::paramList const& paramList, xmlrpc_c::value *   const  retvalP) {
      // results
  	  map<string, xmlrpc_c::value> retData;
      
  	  // static conf data
      const StaticData &staticData = StaticData::Instance();
      
      // parameters to find
      const int nbParams = 18;
      const string parameters[nbParams] = {
      			"galateas-server-id", "galateas-server-name", "galateas-server-description",
      			"galateas-src-language", "galateas-tgt-language", 
      			"galateas-lex", "galateas-rev-lex", "galateas-lambdas"
      			"input-factors", "mapping", "ttable-file", "lmodel-file", "ttable-limit",
      			"weight-d", "weight-l", "weight-t", "weight-w", "n-best-list"
      			// TODO: add other knowns params
      };
      
      for(int i = 0; i < nbParams; i++) {
      	string paramName = parameters[i];
      	vector<string> paramValues = staticData.GetParam(paramName);
      	vector<xmlrpc_c::value> xmlParamValues;
      	// cerr << paramValues.size() << endl;
      	for(vector<string>::iterator it = paramValues.begin(); it != paramValues.end(); ++it) {
      		string value = *it;
      		xmlParamValues.push_back(xmlrpc_c::value_string(value));
      	}
		retData.insert(pair<string, xmlrpc_c::value_array>(paramName, xmlParamValues));	  
      }
      
      *retvalP = xmlrpc_c::value_struct(retData);
  }
};

class DownloadManager: public xmlrpc_c::method {
public:
  DownloadManager() {
  	  this->_signature = "S:S";
      this->_help = "Read a file from the server. The 'file' parameter must be the name of a moses.ini param.";
  }
  void execute(xmlrpc_c::paramList const& paramList, xmlrpc_c::value *   const  retvalP) {
  	  const map<std::string, xmlrpc_c::value> params = paramList.getStruct(0);
      paramList.verifyEnd(1);
      map<std::string, xmlrpc_c::value>::const_iterator si = params.find("file");
      if (si == params.end()) {
          throw xmlrpc_c::fault("Missing file parameter", xmlrpc_c::fault::CODE_PARSE);
      }
      string fileParameter((xmlrpc_c::value_string(si->second)));
        
      // results
  	  map<string, xmlrpc_c::value> retData;
      
  	  // static conf data
      const StaticData &staticData = StaticData::Instance();
      
      // reading file corresponding to the path in moses.ini
      readFileFromMosesParam(staticData, fileParameter, retData);
	  
      *retvalP = xmlrpc_c::value_struct(retData);
  }
  
  void readFileFromMosesParam(const StaticData& staticData, string paramName, map<string, xmlrpc_c::value>& retData) {
  	  vector<string> paramValues = staticData.GetParam(paramName);
      for(vector<string>::iterator it = paramValues.begin(); it != paramValues.end(); ++it) {
      	    string filePath = *it;
			ifstream file(filePath.c_str());
			if(file) {
			   stringstream buffer; 
			   buffer << file.rdbuf();
			   retData.insert(pair<string, xmlrpc_c::value_string>(paramName, xmlrpc_c::value_string(buffer.str())));	 
			   file.close();
			} else {
			   cerr << "No file found (or file not readable) " << filePath << endl;
			}
      }
  }
  
};

class UploadManager: public xmlrpc_c::method {
public:
  UploadManager() {
  	  this->_signature = "S:S";
      this->_help = "Write a file to the server. The 'dir' parameter must be the name of a moses.ini param and represents the destination directory. Mandatory. The 'fileName' parameter represents the destination file name. The 'fileContents' parameter must be the file contents as string or a base64 encoded representation of the file binary contents. Mandatory. The 'overwrite' parameter presence indicate the transfert will overwrite any existing file. Default false.\n";
  }
  void execute(xmlrpc_c::paramList const& paramList, xmlrpc_c::value *   const  retvalP) {
  	  const map<std::string, xmlrpc_c::value> params = paramList.getStruct(0);
      paramList.verifyEnd(1);
      
      // file name
      map<std::string, xmlrpc_c::value>::const_iterator si = params.find("fileName");
      if (si == params.end()) {
          throw xmlrpc_c::fault("Missing fileName parameter", xmlrpc_c::fault::CODE_PARSE);
      }
      string fileName((xmlrpc_c::value_string(si->second)));
      
      // file contents
      si = params.find("fileContents");
      if (si == params.end()) {
          throw xmlrpc_c::fault("Missing fileContents parameter", xmlrpc_c::fault::CODE_PARSE);
      }
      string fileContents((xmlrpc_c::value_string(si->second)));
      
      // destination directory
      si = params.find("dir");
      if (si == params.end()) {
          throw xmlrpc_c::fault("Missing dir parameter", xmlrpc_c::fault::CODE_PARSE);
      }
      string dirParameter((xmlrpc_c::value_string(si->second)));
        
      // overwrite option
      si = params.find("overwrite");
      bool overwrite = (si != params.end());
      
      // static conf data
      const StaticData &staticData = StaticData::Instance();
      
      // checking destination directory existence in moses.ini file
      vector<string> paramValues = staticData.GetParam(dirParameter);
	  if(paramValues.size() < 1){
	  	 throw xmlrpc_c::fault("Unknown dir parameter: " + dirParameter, xmlrpc_c::fault::CODE_PARSE);
	  } else {
	  	  // checking directory existence on the server
	  	  string dirPath = paramValues[0];
	  	  struct stat status;
	  	  if(stat(dirPath.c_str(),&status) == 0) { 
	  	  	// checking overwrite permissions
	  	  	string filePath = dirPath + "/" + fileName; 
	  	  	if(!overwrite && stat(filePath.c_str(),&status) == 0){
	  	  		throw xmlrpc_c::fault("File " + fileName + " already exists in " + dirPath + ". Use 'over' parameter to overwrite.", xmlrpc_c::fault::CODE_PARSE);
	  	  	}
	  	  	// writing contents
	  	  	ofstream myFile (filePath.c_str());
			myFile.write (fileContents.c_str(), fileContents.size());
			myFile.close();
	  	  } else { 
	  	  	throw xmlrpc_c::fault("Dir parameter doesn't exists on server: " + dirParameter, xmlrpc_c::fault::CODE_PARSE);
	  	  }
	  }
      
      *retvalP = xmlrpc_c::value_string("File uploaded to " + dirParameter + " directory.");
  }
    
};

class Updater: public xmlrpc_c::method
{
public:
  Updater() {
    // signature and help strings are documentation -- the client
    // can query this information with a system.methodSignature and
    // system.methodHelp RPC.
    this->_signature = "S:S";
    this->_help = "Updates stuff";
  }
  void
  execute(xmlrpc_c::paramList const& paramList,
          xmlrpc_c::value *   const  retvalP) {
    const params_t params = paramList.getStruct(0);
    breakOutParams(params);
    const TranslationSystem& system = getTranslationSystem(params);
    const PhraseDictionaryFeature* pdf = system.GetPhraseDictionaries()[0];
    PhraseDictionaryDynSuffixArray* pdsa = (PhraseDictionaryDynSuffixArray*) pdf->GetDictionary();
    cerr << "Inserting into address " << pdsa << endl;
    pdsa->insertSnt(source_, target_, alignment_);
    if(add2ORLM_) {       
      updateORLM();
    }
    cerr << "Done inserting\n";
    //PhraseDictionary* pdsa = (PhraseDictionary*) pdf->GetDictionary(*dummy);
    map<string, xmlrpc_c::value> retData;
    //*retvalP = xmlrpc_c::value_struct(retData);
    pdf = 0;
    pdsa = 0;
    *retvalP = xmlrpc_c::value_string("Phrase table updated");
  }
  string source_, target_, alignment_;
  bool bounded_, add2ORLM_;
  void updateORLM() {
    // TODO(level101): this belongs in the language model, not in moseserver.cpp
    vector<string> vl;
    map<vector<string>, int> ngSet;
    LMList lms = StaticData::Instance().GetLMList(); // get LM
    LMList::const_iterator lmIter = lms.begin();
    LanguageModelORLM* orlm = static_cast<LanguageModelORLM*>(static_cast<LMRefCount*>(*lmIter)->MosesServerCppShouldNotHaveLMCode());
    if(orlm == 0) {
      cerr << "WARNING: Unable to add target sentence to ORLM\n";
      return;
    }
    // break out new ngrams from sentence
    const int ngOrder(orlm->GetNGramOrder());
    const std::string sBOS = orlm->GetSentenceStart()->GetString();
    const std::string sEOS = orlm->GetSentenceEnd()->GetString();
    Utils::splitToStr(target_, vl, " ");
    // insert BOS and EOS 
    vl.insert(vl.begin(), sBOS); 
    vl.insert(vl.end(), sEOS);
    for(int j=0; j < vl.size(); ++j) {
      int i = (j<ngOrder) ? 0 : j-ngOrder+1;
      for(int t=j; t >= i; --t) {
        vector<string> ngVec;
        for(int s=t; s<=j; ++s) {
          ngVec.push_back(vl[s]);
          //cerr << vl[s] << " ";
        }
        ngSet[ngVec]++;
        //cerr << endl;
      }
    }
    // insert into LM in order from 1grams up (for LM well-formedness)
    cerr << "Inserting " << ngSet.size() << " ngrams into ORLM...\n";
    for(int i=1; i <= ngOrder; ++i) {
      iterate(ngSet, it) {
        if(it->first.size() == i)
          orlm->UpdateORLM(it->first, it->second);
      }
    }
  }
  void breakOutParams(const params_t& params) {
    params_t::const_iterator si = params.find("source");
    if(si == params.end())
      throw xmlrpc_c::fault("Missing source sentence", xmlrpc_c::fault::CODE_PARSE);
    source_ = xmlrpc_c::value_string(si->second);
    cerr << "source = " << source_ << endl;
    si = params.find("target");
    if(si == params.end())
      throw xmlrpc_c::fault("Missing target sentence", xmlrpc_c::fault::CODE_PARSE);
    target_ = xmlrpc_c::value_string(si->second);
    cerr << "target = " << target_ << endl;
    si = params.find("alignment");
    if(si == params.end())
      throw xmlrpc_c::fault("Missing alignment", xmlrpc_c::fault::CODE_PARSE);
    alignment_ = xmlrpc_c::value_string(si->second);
    cerr << "alignment = " << alignment_ << endl;
    si = params.find("bounded");
    bounded_ = (si != params.end());
    si = params.find("updateORLM");
    add2ORLM_ = (si != params.end());
  }
};

class Translator : public xmlrpc_c::method
{
public:
  Translator() {
    // signature and help strings are documentation -- the client
    // can query this information with a system.methodSignature and
    // system.methodHelp RPC.
    this->_signature = "S:S";
    this->_help = "Does translation";
  }

  void
  execute(xmlrpc_c::paramList const& paramList,
          xmlrpc_c::value *   const  retvalP) {

    const params_t params = paramList.getStruct(0);
    paramList.verifyEnd(1);
    params_t::const_iterator si = params.find("text");
    if (si == params.end()) {
      throw xmlrpc_c::fault(
        "Missing source text",
        xmlrpc_c::fault::CODE_PARSE);
    }
    const string source(
      (xmlrpc_c::value_string(si->second)));

    cerr << "Input: " << source << endl;
    si = params.find("align");
    bool addAlignInfo = (si != params.end());
    si = params.find("sg");
    bool addGraphInfo = (si != params.end());
    si = params.find("topt");
    bool addTopts = (si != params.end());
    si = params.find("report-all-factors");
    bool reportAllFactors = (si != params.end());

    // n-best translations params
    si = params.find("nBestSize");
    bool addNBest = (si != params.end());
    size_t nBestSize = xmlrpc_c::value_int(si->second);
    if(nBestSize <= 0){
    	throw xmlrpc_c::fault("nBestSize must be >= 1",xmlrpc_c::fault::CODE_PARSE);
    }
        
    // features param if we need to get all features scores
    si = params.find("features");
    bool getFeatures = (si != params.end()) ? xmlrpc_c::value_boolean(si->second) : xmlrpc_c::value_boolean(false);
    
    const StaticData &staticData = StaticData::Instance();

    if (addGraphInfo) {
      (const_cast<StaticData&>(staticData)).SetOutputSearchGraph(true);
    }

    const TranslationSystem& system = getTranslationSystem(params);
    stringstream out, graphInfo, transCollOpts;
    map<string, xmlrpc_c::value> retData;

    if (staticData.IsChart()) {
       TreeInput tinput; 
        const vector<FactorType> &inputFactorOrder =
          staticData.GetInputFactorOrder();
        stringstream in(source + "\n");
        tinput.Read(in,inputFactorOrder);
        ChartManager manager(tinput, &system);
        manager.ProcessSentence();
        const ChartHypothesis *hypo = manager.GetBestHypothesis();
        outputChartHypo(out,hypo);
    } else {
        Sentence sentence;
        const vector<FactorType> &inputFactorOrder =
          staticData.GetInputFactorOrder();
        stringstream in(source + "\n");
        sentence.Read(in,inputFactorOrder);
	size_t lineNumber = 0; // TODO: Include sentence request number here?
        Manager manager(lineNumber, sentence, staticData.GetSearchAlgorithm(), &system);
        manager.ProcessSentence();
        const Hypothesis* hypo = manager.GetBestHypothesis();

        vector<xmlrpc_c::value> alignInfo;
        outputHypo(out,hypo,addAlignInfo,alignInfo,reportAllFactors);

    map<string, xmlrpc_c::value> retData;
    pair<string, xmlrpc_c::value>
    text("text", xmlrpc_c::value_string(out.str()));
    cerr << "Output: " << out.str() << endl;
        if (addAlignInfo) {
          retData.insert(pair<string, xmlrpc_c::value>("align", xmlrpc_c::value_array(alignInfo)));
        }
    retData.insert(text);

    if(addNBest) {
		// getting n-best list
		TrellisPathList nBestList;
		manager.CalcNBest(nBestSize, nBestList, true);
		// cerr << "Total n-best size " << nBestList.GetSize() << endl;
		
		// building XML RPC n-best list
		vector<xmlrpc_c::value> nBestListData;
		
		// getting min(nBestIterator, nBestSize) results: 
		// i.e. nBestSize values, as asked by the user, or less, depending on the iterator size
		TrellisPathList::const_iterator iter;
		for (iter = nBestList.begin() ; iter != nBestList.end() && nBestListData.size() <= nBestSize; ++iter) {
		
			/*
			 * creating translation data as a map:
			 * - value = the translation
			 * - score = the translation score
			 * - features = map<string,vector<double>> // featureName/featureValues
			 * - alignInfos = vector<map<string,int>> // phrase align informations
			 *    example moses cmd output: 0-3=0-2 4=3
			 *    this output: {[src-start:0, src-end:3, tgt-start:0, tgt-end:3], [src-start:4, tgt-start:3]}
			 */
			map<string, xmlrpc_c::value> translationData;
			
			// getting a nbest value
			const TrellisPath &path = **iter;
			string pathString = path.GetSurfacePhrase().GetStringRep(staticData.GetOutputFactorOrder());
			translationData["value"] = xmlrpc_c::value_string(pathString);
			
			// getting the translation score
			translationData["score"] = xmlrpc_c::value_double(path.GetTotalScore());
			
			// getting features scores if necessary
			if(getFeatures) {
				map<string, xmlrpc_c::value> features;
				
				// distortion
				vector<xmlrpc_c::value> distortion;
				float d = path.GetScoreBreakdown().GetScoreForProducer((const Moses::ScoreProducer*)system.GetDistortionProducer());
				distortion.push_back(xmlrpc_c::value_double(d));
				features["d"] = xmlrpc_c::value_array(distortion);
				
				// language model
				vector<xmlrpc_c::value> languageModel;
				const LMList& lml = system.GetLanguageModels();
				if (lml.size() > 0) {
					LMList::const_iterator lmi = lml.begin();
					for (; lmi != lml.end(); ++lmi) {
						float score = path.GetScoreBreakdown().GetScoreForProducer(*lmi);
						languageModel.push_back(xmlrpc_c::value_double(score));
					}
				}
				features["lm"] = xmlrpc_c::value_array(languageModel);
				
				// translation model
				if (staticData.GetInputType() == SentenceInput){
					vector<xmlrpc_c::value> translationModel;
					vector<PhraseDictionaryFeature*> pds = system.GetPhraseDictionaries();
					if (pds.size() > 0) {
						vector<PhraseDictionaryFeature*>::iterator pdsIter;
						for(pdsIter = pds.begin(); pdsIter != pds.end() ; ++pdsIter) {
							vector<float> scores = path.GetScoreBreakdown().GetScoresForProducer(*pdsIter);
							for (size_t i = 0; i < scores.size(); ++i) {
								translationModel.push_back(xmlrpc_c::value_double(scores[i]));
							}
						}
					}
					features["tm"] = xmlrpc_c::value_array(translationModel);
				}
				
				// word penality
				vector<xmlrpc_c::value> wordPenality;
				float w = path.GetScoreBreakdown().GetScoreForProducer((const Moses::ScoreProducer*)system.GetWordPenaltyProducer());
				wordPenality.push_back(xmlrpc_c::value_double(w));
				features["w"] = xmlrpc_c::value_array(wordPenality);
				
				// adding features to the translation data
				translationData["features"] = xmlrpc_c::value_struct(features);
			}
			
			// phrase align information
			if(addAlignInfo){
				vector<xmlrpc_c::value> phraseAlignInfo;
				
				const std::vector<const Hypothesis *> &edges = path.GetEdges();
				for (int currEdge = (int)edges.size() - 2 ; currEdge >= 0 ; currEdge--) {
					const Hypothesis &edge = *edges[currEdge];
					const WordsRange &sourceRange = edge.GetCurrSourceWordsRange();
					WordsRange targetRange = path.GetTargetWordsRange(edge);
					map<string, xmlrpc_c::value> info;
					info["src-start"] = xmlrpc_c::value_int(sourceRange.GetStartPos());
					if (sourceRange.GetStartPos() < sourceRange.GetEndPos()) {
						info["src-end"] = xmlrpc_c::value_int(sourceRange.GetEndPos());
					}
					info["tgt-start"] = xmlrpc_c::value_int(targetRange.GetStartPos());
					if (targetRange.GetStartPos() < targetRange.GetEndPos()) {
						info["tgt-end"] = xmlrpc_c::value_int(targetRange.GetEndPos());
					}
					phraseAlignInfo.push_back(xmlrpc_c::value_struct(info));
				}
				
				translationData["phraseAlignInfo"] = xmlrpc_c::value_array(phraseAlignInfo);
			}
			
			// pushing data in the XMLRPC structure
			nBestListData.push_back(xmlrpc_c::value_struct(translationData));
		}
		
		retData.insert(pair<string, xmlrpc_c::value>("nbest", xmlrpc_c::value_array(nBestListData)));
	}

        if(addGraphInfo) {
          insertGraphInfo(manager,retData);
            (const_cast<StaticData&>(staticData)).SetOutputSearchGraph(false);
        }
        if (addTopts) {
          insertTranslationOptions(manager,retData);
        }
    }
    *retvalP = xmlrpc_c::value_struct(retData);
  }

  void outputHypo(ostream& out, const Hypothesis* hypo, bool addAlignmentInfo, vector<xmlrpc_c::value>& alignInfo, bool reportAllFactors = false) {
    if (hypo->GetPrevHypo() != NULL) {
      outputHypo(out,hypo->GetPrevHypo(),addAlignmentInfo, alignInfo, reportAllFactors);
      Phrase p = hypo->GetCurrTargetPhrase();
      if(reportAllFactors) {
        out << p << " ";
      } else {
        for (size_t pos = 0 ; pos < p.GetSize() ; pos++) {
          const Factor *factor = p.GetFactor(pos, 0);
          out << *factor << " ";
        }
      }

      if (addAlignmentInfo) {
        /**
         * Add the alignment info to the array. This is in target order and consists of
         *       (tgt-start, src-start, src-end) triples.
         **/
        map<string, xmlrpc_c::value> phraseAlignInfo;
        phraseAlignInfo["tgt-start"] = xmlrpc_c::value_int(hypo->GetCurrTargetWordsRange().GetStartPos());
        phraseAlignInfo["src-start"] = xmlrpc_c::value_int(hypo->GetCurrSourceWordsRange().GetStartPos());
        phraseAlignInfo["src-end"] = xmlrpc_c::value_int(hypo->GetCurrSourceWordsRange().GetEndPos());
        alignInfo.push_back(xmlrpc_c::value_struct(phraseAlignInfo));
      }
    }
  }

  void outputChartHypo(ostream& out, const ChartHypothesis* hypo) {
    Phrase outPhrase(20);
    hypo->CreateOutputPhrase(outPhrase);

    // delete 1st & last
    assert(outPhrase.GetSize() >= 2);
    outPhrase.RemoveWord(0);
    outPhrase.RemoveWord(outPhrase.GetSize() - 1);
    for (size_t pos = 0 ; pos < outPhrase.GetSize() ; pos++) {
      const Factor *factor = outPhrase.GetFactor(pos, 0);
      out << *factor << " ";
    }

  }

  void insertGraphInfo(Manager& manager, map<string, xmlrpc_c::value>& retData) {
    vector<xmlrpc_c::value> searchGraphXml;
    vector<SearchGraphNode> searchGraph;
    manager.GetSearchGraph(searchGraph);
    for (vector<SearchGraphNode>::const_iterator i = searchGraph.begin(); i != searchGraph.end(); ++i) {
      map<string, xmlrpc_c::value> searchGraphXmlNode;
      searchGraphXmlNode["forward"] = xmlrpc_c::value_double(i->forward);
      searchGraphXmlNode["fscore"] = xmlrpc_c::value_double(i->fscore);
      const Hypothesis* hypo = i->hypo;
      searchGraphXmlNode["hyp"] = xmlrpc_c::value_int(hypo->GetId());
      searchGraphXmlNode["stack"] = xmlrpc_c::value_int(hypo->GetWordsBitmap().GetNumWordsCovered());
      if (hypo->GetId() != 0) {
        const Hypothesis *prevHypo = hypo->GetPrevHypo();
        searchGraphXmlNode["back"] = xmlrpc_c::value_int(prevHypo->GetId());
        searchGraphXmlNode["score"] = xmlrpc_c::value_double(hypo->GetScore());
        searchGraphXmlNode["transition"] = xmlrpc_c::value_double(hypo->GetScore() - prevHypo->GetScore());
        if (i->recombinationHypo) {
          searchGraphXmlNode["recombined"] = xmlrpc_c::value_int(i->recombinationHypo->GetId());
        }
        searchGraphXmlNode["cover-start"] = xmlrpc_c::value_int(hypo->GetCurrSourceWordsRange().GetStartPos());
        searchGraphXmlNode["cover-end"] = xmlrpc_c::value_int(hypo->GetCurrSourceWordsRange().GetEndPos());
        searchGraphXmlNode["out"] =
          xmlrpc_c::value_string(hypo->GetCurrTargetPhrase().GetStringRep(StaticData::Instance().GetOutputFactorOrder()));
      }
      searchGraphXml.push_back(xmlrpc_c::value_struct(searchGraphXmlNode));
    }
    retData.insert(pair<string, xmlrpc_c::value>("sg", xmlrpc_c::value_array(searchGraphXml)));
  }

  void insertTranslationOptions(Manager& manager, map<string, xmlrpc_c::value>& retData) {
    const TranslationOptionCollection* toptsColl = manager.getSntTranslationOptions();
    vector<xmlrpc_c::value> toptsXml;
    for (size_t startPos = 0 ; startPos < toptsColl->GetSize() ; ++startPos) {
      size_t maxSize = toptsColl->GetSize() - startPos;
      size_t maxSizePhrase = StaticData::Instance().GetMaxPhraseLength();
      maxSize = std::min(maxSize, maxSizePhrase);

      for (size_t endPos = startPos ; endPos < startPos + maxSize ; ++endPos) {
        WordsRange range(startPos,endPos);
        const TranslationOptionList& fullList = toptsColl->GetTranslationOptionList(range);
        for (size_t i = 0; i < fullList.size(); i++) {
          const TranslationOption* topt = fullList.Get(i);
          map<string, xmlrpc_c::value> toptXml;
          toptXml["phrase"] = xmlrpc_c::value_string(topt->GetTargetPhrase().
                              GetStringRep(StaticData::Instance().GetOutputFactorOrder()));
          toptXml["fscore"] = xmlrpc_c::value_double(topt->GetFutureScore());
          toptXml["start"] =  xmlrpc_c::value_int(startPos);
          toptXml["end"] =  xmlrpc_c::value_int(endPos);
          vector<xmlrpc_c::value> scoresXml;
          const std::valarray<FValue> &scores = topt->GetScoreBreakdown().getCoreFeatures();
          for (size_t j = 0; j < scores.size(); ++j) {
            scoresXml.push_back(xmlrpc_c::value_double(scores[j]));
          }
          toptXml["scores"] = xmlrpc_c::value_array(scoresXml);
          toptsXml.push_back(xmlrpc_c::value_struct(toptXml));
        }
      }
    }
    retData.insert(pair<string, xmlrpc_c::value>("topt", xmlrpc_c::value_array(toptsXml)));
  }



};


int main(int argc, char** argv)
{

  //Extract port and log, send other args to moses
  char** mosesargv = new char*[argc+2];
  int mosesargc = 0;
  int port = 8080;
  const char* logfile = "/dev/null";
  bool isSerial = false;

  for (int i = 0; i < argc; ++i) {
    if (!strcmp(argv[i],"--server-port")) {
      ++i;
      if (i >= argc) {
        cerr << "Error: Missing argument to --server-port" << endl;
        exit(1);
      } else {
        port = atoi(argv[i]);
      }
    } else if (!strcmp(argv[i],"--server-log")) {
      ++i;
      if (i >= argc) {
        cerr << "Error: Missing argument to --server-log" << endl;
        exit(1);
      } else {
        logfile = argv[i];
      }
    } else if (!strcmp(argv[i], "--serial")) {
      cerr << "Running single-threaded server" << endl;
      isSerial = true;
    } else {
      mosesargv[mosesargc] = new char[strlen(argv[i])+1];
      strcpy(mosesargv[mosesargc],argv[i]);
      ++mosesargc;
    }
  }

  Parameter* params = new Parameter();
  if (!params->LoadParam(mosesargc,mosesargv)) {
    params->Explain();
    exit(1);
  }
  if (!StaticData::LoadDataStatic(params, argv[0])) {
    exit(1);
  }

  xmlrpc_c::registry myRegistry;

  xmlrpc_c::methodPtr const translator(new Translator);
  xmlrpc_c::methodPtr const updater(new Updater);

  myRegistry.addMethod("translate", translator);
  myRegistry.addMethod("updater", updater);
  // discovery method
  xmlrpc_c::methodPtr const discoveryManager(new DiscoveryManager);
  myRegistry.addMethod("discover", discoveryManager);
  
  // download & upload methods
  xmlrpc_c::methodPtr const downloadManager(new DownloadManager);
  myRegistry.addMethod("download", downloadManager);
  xmlrpc_c::methodPtr const uploadManager(new UploadManager);
  myRegistry.addMethod("upload", uploadManager);
  

  xmlrpc_c::serverAbyss myAbyssServer(
    myRegistry,
    port,              // TCP port on which to listen
    logfile
  );

  cerr << "Listening on port " << port << endl;
  if (isSerial) {
    while(1) {
      myAbyssServer.runOnce();
    }
  } else {
    myAbyssServer.run();
  }
  // xmlrpc_c::serverAbyss.run() never returns
  CHECK(false);
  return 0;
}
