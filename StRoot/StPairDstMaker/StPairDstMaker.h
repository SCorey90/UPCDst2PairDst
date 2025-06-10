#ifndef ST_PAIR_DST_MAKER_H
#define ST_PAIR_DST_MAKER_H

#include "Rtypes.h"
#include "StMaker.h"
#include "TString.h"
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "FemtoPairFormat.h"
#include "../StUpcDst/StUPCEvent.h"
#include <map>

class StPairDstMaker : public StMaker {
public:
    StPairDstMaker(const char* name = "StPairDstMaker");
    virtual ~StPairDstMaker();

    virtual Int_t Init();
    virtual Int_t Make();
    virtual Int_t Finish();

    bool eventSelection(StUPCEvent* evt);
    bool trackSelection(StUPCTrack* track1, StUPCTrack* track2);

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

    // ints
    Long64_t fTotalEvents;
    Long64_t fPassTrigger;
    Long64_t fPassNTracks;
    Long64_t fPassVertex;
    Long64_t fPassPID;

    // Histograms
    TH1F* hEventCounter;

    TH1F* hTriggerId;
    TH1F* hNPrimTracksPreCut;
    TH1F* hNPrimTracks;
    TH1F* hNPrimVertices;
    TH1F* hChiPiPi;
    TH1F* hDcaXY1;
    TH1F* hDcaXY2;
    TH1F* hNHitsFit1;
    TH1F* hNHitsFit2;
    TH1F* hNHitsDedx1;
    TH1F* hNHitsDedx2;
    TH1F* hVertexR;

    TChain* fChain;
    StUPCEvent* fUpcEvt;
    //Int_t fTriggerId;  // Desired trigger ID
    std::vector<Int_t> fTriggerIds;  // Desired trigger IDs

    void ResetFemtoPair();

    ClassDef(StPairDstMaker, 1)
};

#endif
