#include "StPairDstMaker.h"
#include "../StUpcDst/StUPCEvent.h"
#include "../StUpcDst/StUPCTrack.h"
#include "../StUpcDst/StUPCVertex.h"
#include "TLorentzVector.h"
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TH1F.h"
#include <iostream>
#include <map>
#include <utility>

ClassImp(StPairDstMaker)

StPairDstMaker::StPairDstMaker(const char* name) : StMaker(name), fOutFile(nullptr), fTree(nullptr), fChain(new TChain("mUPCTree")), fUpcEvt(nullptr), fTriggerIds({-1}), fTotalEvents(0), fPassTrigger(0), fPassNTracks(0), fPassVertex(0), fPassPID(0), fOutputFile("PairDst.root") {
}

StPairDstMaker::~StPairDstMaker() {
}

Int_t StPairDstMaker::Init() {
    fChain->SetBranchAddress("mUPCEvent", &fUpcEvt);

    fOutFile = new TFile(fOutputFile, "RECREATE");
    fTree = new TTree("PairDst", "PairDst");
    fTree->Branch("Pairs", &fFemtoPair);

    // Create histograms
    hEventCounter = new TH1F("hEventCounter", "Event Selection", 4, 0, 4);
    hEventCounter->GetXaxis()->SetBinLabel(1, "Total Events");
    hEventCounter->GetXaxis()->SetBinLabel(2, "Trigger");
    hEventCounter->GetXaxis()->SetBinLabel(3, "NTracks=2");
    hEventCounter->GetXaxis()->SetBinLabel(4, "PID");

    hTriggerId = new TH1F("hTriggerId", "Trigger ID", 1000, 450000, 451000);
    hZDCWest = new TH1F("hZDCWest", "; ZDC ADC; Yield", 200, 0, 1200);
    hZDCEast = new TH1F("hZDCEast", "; ZDC ADC; Yield", 200, 0, 1200);
    hZDCWestvsEast = new TH2F("hZDCWestvsEast", "; ZDC West ADC; ZDC East ADC; Yield", 200, 0, 1200, 200, 0, 1200);
    hNPrimTracksPreCut = new TH1F("hNPrimTracksPreCut", "Number of Primary Tracks", 10, 0, 10);
    hNPrimTracks = new TH1F("hNPrimTracks", "Number of Primary Tracks", 10, 0, 10);
    hNPrimVertices = new TH1F("hNPrimVertices", "Number of Primary Vertices", 10, 0, 10);
    hChiPiPi = new TH1F("hChiPiPi", "Chi^2 of Pion Hypothesis", 100, 0, 50);
    hDcaXY1 = new TH1F("hDcaXY1", "DCA XY of Track 1", 100, 0, 10);
    hDcaXY2 = new TH1F("hDcaXY2", "DCA XY of Track 2", 100, 0, 10);
    hNHitsFit1 = new TH1F("hNHitsFit1", "Number of Hits Fit for Track 1", 50, 0, 50);
    hNHitsFit2 = new TH1F("hNHitsFit2", "Number of Hits Fit for Track 2", 50, 0, 50);
    hNHitsDedx1 = new TH1F("hNHitsDedx1", "Number of Hits dE/dx for Track 1", 50, 0, 50);
    hNHitsDedx2 = new TH1F("hNHitsDedx2", "Number of Hits dE/dx for Track 2", 50, 0, 50);
    hVertexR = new TH1F("hVertexR", "Vertex R", 100, 0, 10);

    return kStOK;
}

bool StPairDstMaker::manualUPCmain(StUPCEvent* evt) {
    if (!evt) return false;

    // Check if the event has BBC activity
    if (evt->getBBCSmallEast() > 50 || evt->getBBCSmallWest() > 50) {
        return false;
    }

    // Check if 2<= TOF Hits <= 6
    if (evt->getTOFMultiplicity() < 2 || evt->getTOFMultiplicity() > 6) {
        return false;
    }

    return true;
}

bool StPairDstMaker::eventSelection(StUPCEvent *evt){
    if (!evt) return false;

    fTotalEvents++; //iterate total events

    int nTracks = evt->getNPrimTracks();
    int nVertices = evt->getNPrimVertices();
    bool isTriggered = false;
    // bool isTriggered = true;

    // Check trigger
    for (int i = 0; i < fTriggerIds.size(); ++i) {
        if (evt->isTrigger(fTriggerIds[i])) {
            isTriggered = true;
            break;
        }
    }

    if (!isTriggered) return false;
    if (!manualUPCmain(evt)) return false;
    for (int i = 0; i < fTriggerIds.size(); ++i) {
        if (evt->isTrigger(fTriggerIds[i])) {
            hTriggerId->Fill(fTriggerIds[i]);
        }
    }

    fPassTrigger++; //iterate triggered events when one passes the trigger

    hNPrimVertices->Fill(nVertices);
    hNPrimTracksPreCut->Fill(nTracks);

    hZDCEast->Fill(evt->getZDCUnAttEast());
    hZDCWest->Fill(evt->getZDCUnAttWest());
    hZDCWestvsEast->Fill(evt->getZDCUnAttWest(), evt->getZDCUnAttEast());

    return true;
}

bool StPairDstMaker::trackSelection(StUPCTrack* track1, StUPCTrack* track2) {
    if (!track1 || !track2) return false;

    double nSigmaPi1 = track1->getNSigmasTPCPion();
    double nSigmaPi2 = track2->getNSigmasTPCPion();
    double chipipi2 = nSigmaPi1 * nSigmaPi1 + nSigmaPi2 * nSigmaPi2;

    double dcaXY1 = track1->getDcaXY();
    double dcaXY2 = track2->getDcaXY();

    int nHitsFit1 = track1->getNhitsFit();
    int nHitsFit2 = track2->getNhitsFit();
    int nHitsDEDx1 = track1->getNhitsDEdx();
    int nHitsDEDx2 = track2->getNhitsDEdx();

    hChiPiPi->Fill(chipipi2);
    hDcaXY1->Fill(dcaXY1);
    hDcaXY2->Fill(dcaXY2);
    hNHitsFit1->Fill(nHitsFit1);
    hNHitsFit2->Fill(nHitsFit2);
    hNHitsDedx1->Fill(nHitsDEDx1);
    hNHitsDedx2->Fill(nHitsDEDx2);

    // if (chipipi2 > 20) return false;

    if (dcaXY1 > 3 && dcaXY2 > 3) return false;

    // if (nHitsFit1 < 8 || nHitsFit2 < 8 || nHitsDEDx1 < 5 || nHitsDEDx2 < 5) return false;

    fPassPID++;

    return true;
}

Int_t StPairDstMaker::Make() {
    Long64_t nEntries = fChain->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        fChain->GetEntry(i);

        int nTracks = fUpcEvt->getNPrimTracks();
        std::map<int, std::vector<StUPCTrack*>> trackByVertex;
        for (int j=0; j<nTracks; j++) {
            StUPCTrack* track = fUpcEvt->getTrack(j);
            if (!track) continue;
            StUPCVertex* vertex = track->getVertex();
            if (!vertex) continue;
            trackByVertex[vertex->getId()].push_back(track);
        }

        int nVertices = fUpcEvt->getNPrimVertices();

        if (!eventSelection(fUpcEvt)) {
            continue;
        }

        for (int j=0; j<nVertices; j++) {
            StUPCVertex* vertex = fUpcEvt->getVertex(j);
            if (!vertex) continue;
            if (trackByVertex.find(vertex->getId()) == trackByVertex.end()) continue;
            if (trackByVertex[vertex->getId()].size() != 2) continue;
            fPassNTracks++;

            StUPCTrack* track1 = trackByVertex[vertex->getId()][0];
            StUPCTrack* track2 = trackByVertex[vertex->getId()][1];

            if ( !trackSelection(track1, track2) ) {
                continue;
            }

            if (track1->getCharge() * track2->getCharge() < 0) {
                if (track1->getCharge() < 0) {
                    StUPCTrack* temp = track1;
                    track1 = track2;
                    track2 = temp;
                }
            }

            ResetFemtoPair();

            double nVertexX = fUpcEvt->getVertex(0)->getPosX();
            double nVertexY = fUpcEvt->getVertex(0)->getPosY();
            double nVertexR = sqrt(nVertexX * nVertexX + nVertexY * nVertexY);
            hVertexR->Fill(nVertexR);

            if (nVertexR > 2.0) continue; // Vertex R cut

            fFemtoPair.mRunID = fUpcEvt->getRunNumber();
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
            fFemtoPair.d1_mMatchFlag = (tofpathlength1 > 0) ? 1 : 0;
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
            fFemtoPair.d2_mMatchFlag = (tofpathlength2 > 0) ? 1 : 0;
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
    }

    hEventCounter->SetBinContent(1, fTotalEvents);
    hEventCounter->SetBinContent(2, fPassTrigger);
    hEventCounter->SetBinContent(3, fPassNTracks);
    hEventCounter->SetBinContent(4, fPassPID);

    return kStOK;
}

Int_t StPairDstMaker::Finish() {
    fOutFile->cd();
    fTree->Write();

    // Write histograms
    hEventCounter->Write();

    hTriggerId->Write();
    hZDCWest->Write();
    hZDCEast->Write();
    hZDCWestvsEast->Write();
    hNPrimTracksPreCut->Write();
    hNPrimTracks->Write();
    hNPrimVertices->Write();
    hChiPiPi->Write();
    hDcaXY1->Write();
    hDcaXY2->Write();
    hNHitsFit1->Write();
    hNHitsFit2->Write();
    hNHitsDedx1->Write();
    hNHitsDedx2->Write();
    hVertexR->Write();

    TTree* counterTree = new TTree("Counters", "Counters");
    counterTree->Branch("TotalEvents", &fTotalEvents, "TotalEvents/L");
    counterTree->Branch("PassTrigger", &fPassTrigger, "PassTrigger/L");
    counterTree->Branch("PassNTracks", &fPassNTracks, "PassNTracks/L");
    counterTree->Branch("PassVertex", &fPassVertex, "PassVertex/L");
    counterTree->Branch("PassPID", &fPassPID, "PassPID/L");
    counterTree->Fill();
    counterTree->Write();
    delete counterTree;

    fOutFile->Close();
    return kStOK;
}

void StPairDstMaker::SetInputFileList(const char* fileList) {
    std::ifstream infile(fileList);
    std::string line;
    int fileCount = 0;
    while (std::getline(infile, line)) {
        fChain->Add(line.c_str());
        fileCount++;
    }
    std::cout << "Number of files added to the chain: " << fileCount << std::endl;
}

void StPairDstMaker::ResetFemtoPair() {
    fFemtoPair.reset();
}
