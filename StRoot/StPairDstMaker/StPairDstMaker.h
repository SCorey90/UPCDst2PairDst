#ifndef ST_PAIR_DST_MAKER_H
#define ST_PAIR_DST_MAKER_H

#include "StMaker.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "FemtoPairFormat.h"
#include "../StUpcDst/StUPCEvent.h"

class StPairDstMaker : public StMaker {
public:
    StPairDstMaker(const char* name = "StPairDstMaker");
    virtual ~StPairDstMaker();

    virtual Int_t Init();
    virtual Int_t Make();
    virtual Int_t Finish();

    bool eventSelection(StUPCEvent* evt);

    void SetInputFileList(const char* fileList);
    void SetInputFile(const char* file) { fInputFile = file; }
    void SetOutputFile(const char* file) { fOutputFile = file; }
    //void SetTriggerId(Int_t triggerId) { fTriggerId = triggerId; }
    void SetTriggerIds(const std::vector<Int_t>& triggerIds) { fTriggerIds = triggerIds; }

private:
    TString fInputFile;
    TString fOutputFile;
    TFile* fOutFile;
    TTree* fTree;
    FemtoPair fFemtoPair;

    TChain* fChain;
    StUPCEvent* fUpcEvt;
    //Int_t fTriggerId;  // Desired trigger ID
    std::vector<Int_t> fTriggerIds;  // Desired trigger IDs

    void ResetFemtoPair();

    ClassDef(StPairDstMaker, 1)
};

#endif
