#include "StPairDstMaker.h"
#include "../StUpcDst/StUPCEvent.h"
#include "../StUpcDst/StUPCTrack.h"
#include "../StUpcDst/StUPCVertex.h"
#include "TLorentzVector.h"
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include <iostream>
#include <utility>
#include <fstream>

ClassImp(StPairDstMaker)

StPairDstMaker::StPairDstMaker(const char* name) : StMaker(name), fOutFile(nullptr), fTree(nullptr), fChain(new TChain("mUPCTree")), fUpcEvt(nullptr) {
}

StPairDstMaker::~StPairDstMaker() {
}

Int_t StPairDstMaker::Init() {
    fChain->SetBranchAddress("mUPCEvent", &fUpcEvt);

    fOutFile = new TFile(fOutputFile, "RECREATE");
    fTree = new TTree("PairDst", "PairDst");
    fTree->Branch("Pairs", &fFemtoPair);

    return kStOK;
}

bool StPairDstMaker::eventSelection(StUPCEvent* evt) {
    if (!evt) return false;

    int nTracks = evt->getNPrimTracks();
    if (nTracks != 2) return false;

    int nVertices = evt->getNPrimVertices();
    if (nVertices != 1) return false;

    bool triggerPassed = false;
    for (int i = 0; i < fTriggerIds.size(); ++i) {
        if (evt->isTrigger(fTriggerIds[i])) {
            triggerPassed = true;
            break;
        }
    }
    if (!triggerPassed) return false;

    StUPCTrack* track1 = evt->getTrack(0);
    StUPCTrack* track2 = evt->getTrack(1);

    if (!track1 || !track2) return false;

    double nSigmaPi1 = track1->getNSigmasTPCPion();
    double nSigmaPi2 = track2->getNSigmasTPCPion();
    double chipipi2 = nSigmaPi1 * nSigmaPi1 + nSigmaPi2 * nSigmaPi2;

    if (chipipi2 > 20) return false;

    double dcaXY1 = track1->getDcaXY();
    double dcaXY2 = track2->getDcaXY();

    if (dcaXY1 > 3 && dcaXY2 > 3) return false;

    int nHitsFit1 = track1->getNhitsFit();
    int nHitsFit2 = track2->getNhitsFit();
    int nHitsDEDx1 = track1->getNhitsDEdx();
    int nHitsDEDx2 = track2->getNhitsDEdx();

    if (nHitsFit1 < 8 || nHitsFit2 < 8 || nHitsDEDx1 < 5 || nHitsDEDx2 < 5) return false;

    return true;
}

Int_t StPairDstMaker::Make() {
    Long64_t nEntries = fChain->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        fChain->GetEntry(i);

        // Print progress every 1000 entries
        if (i % 1000 == 0) {
            std::cout << "Processing entry " << i << " / " << nEntries << std::endl;
        }

        // Check if the event passes the selection criteria
        if (!eventSelection(fUpcEvt)) {
            continue;
        }

        // Get the two tracks
        StUPCTrack* track1 = fUpcEvt->getTrack(0);
        StUPCTrack* track2 = fUpcEvt->getTrack(1);

        // Ensure track1 is positive and track2 is negative if they have opposite charges
        if (track1->getCharge() * track2->getCharge() < 0) {
            if (track1->getCharge() < 0) {
                StUPCTrack* temp = track1;
                track1 = track2;
                track2 = temp;
            }
        }

        // Fill FemtoPair with track information
        ResetFemtoPair();

        fFemtoPair.mVertexZ = fUpcEvt->getVertex(0)->getPosZ();
        fFemtoPair.mGRefMult = fUpcEvt->getNGlobTracks();
        fFemtoPair.mZDCEast = fUpcEvt->getZDCUnAttEast();
        fFemtoPair.mZDCWest = fUpcEvt->getZDCUnAttWest();
        fFemtoPair.mChargeSum = track1->getCharge() + track2->getCharge();

        fFemtoPair.d1_mPt = track1->getPt();
        fFemtoPair.d1_mEta = track1->getEta();
        fFemtoPair.d1_mPhi = track1->getPhi();
        fFemtoPair.d1_mNSigmaPion = track1->getNSigmasTPCPion();
        fFemtoPair.d1_mNSigmaKaon = track1->getNSigmasTPCKaon();
        fFemtoPair.d1_mNSigmaProton = track1->getNSigmasTPCProton();
        fFemtoPair.d1_mNSigmaElectron = track1->getNSigmasTPCElectron();
        fFemtoPair.d1_mNHitsFit = track1->getNhitsFit();
        fFemtoPair.d1_mNHitsMax = track1->getNhitsMax();
        fFemtoPair.d1_mNHitsDedx = track1->getNhitsDEdx();
        fFemtoPair.d1_mDCA = track1->getDcaXY();
        fFemtoPair.d1_mTof = track1->getTofTime();
        double tofpathlength1 = track1->getTofPathLength();
        fFemtoPair.d1_mMatchFlag = (tofpathlength1 > 1) ? 1 : 0;
        fFemtoPair.d1_mLength = tofpathlength1;

        fFemtoPair.d2_mPt = track2->getPt();
        fFemtoPair.d2_mEta = track2->getEta();
        fFemtoPair.d2_mPhi = track2->getPhi();
        fFemtoPair.d2_mNSigmaPion = track2->getNSigmasTPCPion();
        fFemtoPair.d2_mNSigmaKaon = track2->getNSigmasTPCKaon();
        fFemtoPair.d2_mNSigmaProton = track2->getNSigmasTPCProton();
        fFemtoPair.d2_mNSigmaElectron = track2->getNSigmasTPCElectron();
        fFemtoPair.d2_mNHitsFit = track2->getNhitsFit();
        fFemtoPair.d2_mNHitsMax = track2->getNhitsMax();
        fFemtoPair.d2_mNHitsDedx = track2->getNhitsDEdx();
        fFemtoPair.d2_mDCA = track2->getDcaXY();
        fFemtoPair.d2_mTof = track2->getTofTime();
        double tofpathlength2 = track2->getTofPathLength();
        fFemtoPair.d2_mMatchFlag = (tofpathlength2 > 1) ? 1 : 0;
        fFemtoPair.d2_mLength = tofpathlength2;

        TLorentzVector lv1, lv2, lv;
        lv1.SetPtEtaPhiM(track1->getPt(), track1->getEta(), track1->getPhi(), 0.13957039);
        lv2.SetPtEtaPhiM(track2->getPt(), track2->getEta(), track2->getPhi(), 0.13957039);
        lv = lv1 + lv2;

        fFemtoPair.mPt = lv.Pt();
        fFemtoPair.mEta = lv.Eta();
        fFemtoPair.mPhi = lv.Phi();
        fFemtoPair.mMass = lv.M();
        fFemtoPair.mRapidity = lv.Rapidity();

        fTree->Fill();
    }
    return kStOK;
}

Int_t StPairDstMaker::Finish() {
    fOutFile->cd();
    fTree->Write();
    fOutFile->Close();
    return kStOK;
}

void StPairDstMaker::SetInputFileList(const char* fileList) {
    std::ifstream infile(fileList);
    std::string line;
    while (std::getline(infile, line)) {
        fChain->Add(line.c_str());
    }
}

void StPairDstMaker::ResetFemtoPair() {
    fFemtoPair.reset();
}
