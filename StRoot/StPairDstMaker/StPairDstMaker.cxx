#include "StPairDstMaker.h"
#include "../StUpcDst/StUPCEvent.h"
#include "../StUpcDst/StUPCTrack.h"
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include <iostream>

ClassImp(StPairDstMaker)

StPairDstMaker::StPairDstMaker(const char* name) : StMaker(name), fOutFile(nullptr), fTree(nullptr), fChain(nullptr), fUpcEvt(nullptr) {
}

StPairDstMaker::~StPairDstMaker() {
}

Int_t StPairDstMaker::Init() {
    fChain = new TChain("mUPCTree");
    fChain->Add(fInputFile);
    fChain->SetBranchAddress("mUPCEvent", &fUpcEvt);

    fOutFile = new TFile(fOutputFile, "RECREATE");
    fTree = new TTree("FemtoPairTree", "FemtoPairTree");
    fTree->Branch("FemtoPair", &fFemtoPair);

    return kStOK;
}

Int_t StPairDstMaker::Make() {
    Long64_t nEntries = fChain->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        fChain->GetEntry(i);

        int nTracks = fUpcEvt->getNumberOfTracks();
        for (int j = 0; j < nTracks; ++j) {
            StUPCTrack* track1 = fUpcEvt->getTrack(j);
            if (!track1) continue;

            for (int k = j + 1; k < nTracks; ++k) {
                StUPCTrack* track2 = fUpcEvt->getTrack(k);
                if (!track2) continue;

                // Fill FemtoPair with track information
                ResetFemtoPair();

                fFemtoPair.mVertexZ = fUpcEvt->getVertex(0)->getPosZ();
                fFemtoPair.mGRefMult = fUpcEvt->getRefMult();
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

                TLorentzVector lv1, lv2, lv;
                lv1.SetPtEtaPhiM(track1->getPt(), track1->getEta(), track1->getPhi(), 0.);
                lv2.SetPtEtaPhiM(track2->getPt(), track2->getEta(), track2->getPhi(), 0.493);
                lv = lv1 + lv2;

                fFemtoPair.mPt = lv.Pt();
                fFemtoPair.mEta = lv.Eta();
                fFemtoPair.mPhi = lv.Phi();
                fFemtoPair.mMass = lv.M();
                fFemtoPair.mRapidity = lv.Rapidity();

                fTree->Fill();
            }
        }
    }
    return kStOK;
}

Int_t StPairDstMaker::Finish() {
    fOutFile->cd();
    fTree->Write();
    fOutFile->Close();
    return kStOK;
}

void StPairDstMaker::ResetFemtoPair() {
    fFemtoPair.reset();
}
