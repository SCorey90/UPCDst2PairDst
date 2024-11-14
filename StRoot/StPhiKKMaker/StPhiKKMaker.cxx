#include <algorithm>
#include <cmath>
#include <exception>
#include <fstream>
#include <functional>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
// ROOT headers
#include "TCanvas.h"
#include "TChain.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2.h"
#include "TH2F.h"
#include "TH2I.h"
#include "TH3.h"
#include "TList.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TString.h"
#include "TTimer.h"
#include "TTree.h"
// #include "Math/Vector4Dfwd.h"
#include "Math/Vector4D.h"
#include "TLorentzVector.h"
#include "TVector3.h"
// picoDst headers
#include "../StUpcDst/StUPCEvent.h"
#include "../StUpcDst/StUPCTrack.h"
#include "../StUpcDst/StUPCVertex.h"
#include "TObjString.h"

#include "StPhiKKMaker.h"
#include "KaonPair.h"

using namespace std;
map<string, TH1 *> hists;
TTimer *timer;
TRandom3 *myRandom;
bool passX2(double x2ee, double x2pipi);
float XpipiMaxCut = 800.0;
float XppMaxCut = 800.0;
float XeeMaxCut = 800.0;
float XkkMaxCut = 800.0;
double getShift(double p, TH1D *histo, double pLLimit, double pHLimit);

TFile *outfile = nullptr;

TChain *read_filelist(TString infile) {
  TList nameList;
  char str[2000];
  fstream file_op(infile.Data(), ios::in);
  while (file_op >> str) {
    nameList.Add(new TObjString(str));
  }
  TChain *t = new TChain("mUPCTree");
  TIter next(&nameList);
  TObjString *fileNm;
  int iList = 0;
  while ((fileNm = (TObjString *)next())) {
    t->AddFile(fileNm->String());
    iList++;
  }
  cout << "Added " << iList << " files" << endl;
  return t;
}

void make_hists() {
  const char *EventCounterLabels[] = {"All Events", "Trigger", "nVtx<5",
                                      "nTrk<25","MultPerVertex <=2","Pass PID","ChargeSum 0","In Mass Range", "In Coh P_{T} Range"};
  hists["hEventCounter"] = new TH1I("hEventCounter", "Event Counter", 9, 0, 9);
  for (int i = 1; i <= 9; i++)
    hists["hEventCounter"]->GetXaxis()->SetBinLabel(i,
                                                    EventCounterLabels[i - 1]);
  hists["hRefMult"] = new TH1I("hRefMult", "RefMult", 1000, 0, 1000);
  hists["hAvgMultPerVertex"] =
      new TH1I("hAvgMultPerVertex", "AvgMultPerVertex", 1000, 0, 100);
  hists["hMultPerVertex"] =
      new TH1I("hMultPerVertex", "MultPerVertex", 100, 0, 100);
  hists["hNumVertices"] = new TH1I("hNumVertices", "NumVertices", 40, 0, 40);

  hists["hX2kk"] = new TH1I("hX2kk", "#chi^{2}_{KK}", 4000, 0, 400);

  hists["hUlsMassPtkk"] =
      new TH2F("hUlsMassPt", ";mass;pt", 600, 0.9, 1.5, 1000, 0, 1.0);
  hists["hUlsMasskk"] =
      new TH1D("hUlsMasskk", "hmassID;M_{ee} (GeV/c^{2})", 2000, 0.9, 2.9);
  hists["hUlsPtkk"] = new TH1D("hUlsPtkk", "hptID;p_{T} (GeV/c)", 1500, 0, 1.5);

  hists["hLsMassPtkk"] =
      new TH2F("hLsMassPt", ";mass;pt", 600, 0.9, 1.5, 1000, 0, 1.0);
  hists["hLsMasskk"] =
      new TH1D("hLsMasskk", "hmassID;M_{ee} (GeV/c^{2})", 2000, 0.9, 2.9);
  hists["hLsPtkk"] = new TH1D("hLsPtkk", "hptID;p_{T} (GeV/c)", 1500, 0, 1.5);


  for (auto kv : hists) {
    if (kv.first.find("Ls") != string::npos) {
      kv.second->SetLineColor(kRed);
    }
  }

  // Track QA
  // hists["trackFlag"] = new TH1F( "trackFlag", "trackFlag", 20, 0, 20 )
  string pres[] = {"track", "uls", "ls"};
  for (string prefix : pres) {
    hists[prefix + "Curvature"] =
        new TH1F((prefix + "Curvature").c_str(), (prefix + "Curvature").c_str(),
                 500, 0, 0.05);
    hists[prefix + "DipAngle"] =
        new TH1F((prefix + "DipAngle").c_str(), (prefix + "DipAngle").c_str(),
                 300, -1.5, 1.5);
    hists[prefix + "Phase"] =
        new TH1F((prefix + "Phase").c_str(), (prefix + "Phase").c_str(), 300,
                 -3.14159, 3.14159);
    hists[prefix + "OriginX"] =
        new TH1F((prefix + "OriginX").c_str(), (prefix + "OriginX").c_str(),
                 400, -200, 200);
    hists[prefix + "OriginY"] =
        new TH1F((prefix + "OriginY").c_str(), (prefix + "OriginY").c_str(),
                 400, -200, 200);
    hists[prefix + "OriginZ"] =
        new TH1F((prefix + "OriginZ").c_str(), (prefix + "OriginZ").c_str(),
                 400, -200, 200);
    hists[prefix + "Pt"] =
        new TH1F((prefix + "Pt").c_str(), (prefix + "Pt").c_str(), 1000, 0, 10);
    hists[prefix + "Eta"] = new TH1F((prefix + "Eta").c_str(),
                                     (prefix + "Eta").c_str(), 500, -2.5, 2.5);
    hists[prefix + "Phi"] =
        new TH1F((prefix + "Phi").c_str(), (prefix + "Phi").c_str(), 300,
                 -3.14159, 3.14159);
    hists[prefix + "DcaXY"] = new TH1F((prefix + "DcaXY").c_str(),
                                       (prefix + "DcaXY").c_str(), 200, 0, 10);
    hists[prefix + "DcaZ"] = new TH1F((prefix + "DcaZ").c_str(),
                                      (prefix + "DcaZ").c_str(), 160, -80, 80);
    hists[prefix + "Charge"] = new TH1F((prefix + "Charge").c_str(),
                                        (prefix + "Charge").c_str(), 4, -2, 2);
    hists[prefix + "Nhits"] = new TH1F((prefix + "Nhits").c_str(),
                                       (prefix + "Nhits").c_str(), 100, 0, 100);
    hists[prefix + "NhitsFit"] =
        new TH1F((prefix + "NhitsFit").c_str(), (prefix + "NhitsFit").c_str(),
                 100, 0, 100);
    hists[prefix + "NhitsMax"] =
        new TH1F((prefix + "NhitsMax").c_str(), (prefix + "NhitsMax").c_str(),
                 100, 0, 100);
    hists[prefix + "Chi2"] = new TH1F((prefix + "Chi2").c_str(),
                                      (prefix + "Chi2").c_str(), 200, 0, 20);
    hists[prefix + "NhitsDEdx"] =
        new TH1F((prefix + "NhitsDEdx").c_str(), (prefix + "NhitsDEdx").c_str(),
                 100, 0, 100);
    hists[prefix + "DEdxSignal"] =
        new TH1F((prefix + "DEdxSignal").c_str(),
                 (prefix + "DEdxSignal").c_str(), 800, 0, 1.6e-4);
    hists[prefix + "NSigmasTPCElectron"] =
        new TH1F((prefix + "NSigmasTPCElectron").c_str(),
                 (prefix + "NSigmasTPCElectron").c_str(), 200, -100, 100);
    hists[prefix + "NSigmasTPCPion"] =
        new TH1F((prefix + "NSigmasTPCPion").c_str(),
                 (prefix + "NSigmasTPCPion").c_str(), 200, -100, 100);
    hists[prefix + "NSigmasTPCKaon"] =
        new TH1F((prefix + "NSigmasTPCKaon").c_str(),
                 (prefix + "NSigmasTPCKaon").c_str(), 200, -100, 100);
    hists[prefix + "NSigmasTPCProton"] =
        new TH1F((prefix + "NSigmasTPCProton").c_str(),
                 (prefix + "NSigmasTPCProton").c_str(), 200, -100, 100);
  }
  // TH1I *hNhits = new TH1I("hNhits","Number Of Tpc Hits",80,0,80);
  // TH1D *hNhitsRatio = new TH1D("hNhitsRatio","Number Of Tpc Hits / Number Of
  // Possible Tpc Hits",100,0.,1.); TH1I *hNhitsDedx = new
  // TH1I("hNhitsDedx","Number Of Tpc Dedx Hits",80,0,80); TH1D *hPt = new
  // TH1D("hPt", "hPt", 100, 0, 10); TH1D *hEta = new TH1D("hEta", "#eta;
  // #eta",100, -5, 5);
}

void fillTrackQA(StUPCTrack *trk, string prefix = "track") {
  // hists[prefix+"Flag"] ->Fill( trk->getFlag() );
  hists[prefix + "Curvature"]->Fill(trk->getCurvature());
  hists[prefix + "DipAngle"]->Fill(trk->getDipAngle());
  hists[prefix + "Phase"]->Fill(trk->getPhase());
  hists[prefix + "OriginX"]->Fill(trk->getOrigin().X());
  hists[prefix + "OriginY"]->Fill(trk->getOrigin().Y());
  hists[prefix + "OriginZ"]->Fill(trk->getOrigin().Z());
  hists[prefix + "Pt"]->Fill(trk->getPt());
  hists[prefix + "Eta"]->Fill(trk->getEta());
  hists[prefix + "Phi"]->Fill(trk->getPhi());
  hists[prefix + "DcaXY"]->Fill(trk->getDcaXY());
  hists[prefix + "DcaZ"]->Fill(trk->getDcaZ());
  hists[prefix + "Charge"]->Fill(trk->getCharge());
  hists[prefix + "Nhits"]->Fill(trk->getNhits());
  hists[prefix + "NhitsFit"]->Fill(trk->getNhitsFit());
  hists[prefix + "NhitsMax"]->Fill(trk->getNhitsMax());
  hists[prefix + "Chi2"]->Fill(trk->getChi2());
  hists[prefix + "NhitsDEdx"]->Fill(trk->getNhitsDEdx());
  hists[prefix + "DEdxSignal"]->Fill(trk->getDEdxSignal());
  hists[prefix + "NSigmasTPCElectron"]->Fill(trk->getNSigmasTPCElectron());
  hists[prefix + "NSigmasTPCPion"]->Fill(trk->getNSigmasTPCPion());
  hists[prefix + "NSigmasTPCKaon"]->Fill(trk->getNSigmasTPCKaon());
  hists[prefix + "NSigmasTPCProton"]->Fill(trk->getNSigmasTPCProton());
}

TLorentzVector getTrackKinne(StUPCTrack *trk, double mass) {
  TLorentzVector outputPtr;
  double pt;
  double eta;
  double phi;
  trk->getPtEtaPhi(pt, eta, phi);
  outputPtr.SetPtEtaPhiM(pt, eta, phi, mass);
  return outputPtr;
}

void post() {
  outfile->cd();
  for (auto kv : hists) {
    kv.second->Write();
  }
  outfile->Close();
}

//_____________________________________________________Main
// Algorithm______________________________________________________
void make_ana(
    TString infile = "/work/data/UPCDstList.txt",
    TString outputFile = "output.root") {

  make_hists();

  // Read file list and setup TChain
  TChain *t = read_filelist(infile);
  // open output file
  outfile = TFile::Open(outputFile, "recreate");
  outfile->cd();

  // Create a map of TTree pointers for different multPerVertex values
  map<int, TTree*> trackPairTrees;
  map<int, Double_t> track1_pt, track1_eta, track1_phi, track1_NSigmaPion,
      track1_NSigmaKaon, track1_NSigmaProton, track1_NSigmaElectron, track1_charge;
  map<int, Double_t> track2_pt, track2_eta, track2_phi, track2_NSigmaPion,
      track2_NSigmaKaon, track2_NSigmaProton, track2_NSigmaElectron, track2_charge;
  map<int, Int_t> track1_NHitsdEdx, track2_NHitsdEdx;
  map<int, Int_t> track1_index, track2_index;

  std::vector<int> mPV_values = {2,3};
  // Initialize trees and branches for the required multPerVertex values
  for (int mpv : mPV_values) {
    TString treeName = TString::Format("two-track_information_mpv_%d", mpv);
    trackPairTrees[mpv] = new TTree(treeName, "TTrack");

    trackPairTrees[mpv]->Branch("track1pt", &track1_pt[mpv]);
    trackPairTrees[mpv]->Branch("track1eta", &track1_eta[mpv]);
    trackPairTrees[mpv]->Branch("track1phi", &track1_phi[mpv]);
    trackPairTrees[mpv]->Branch("track1NSigmaPion", &track1_NSigmaPion[mpv]);
    trackPairTrees[mpv]->Branch("track1NSigmaKaon", &track1_NSigmaKaon[mpv]);
    trackPairTrees[mpv]->Branch("track1NSigmaProton", &track1_NSigmaProton[mpv]);
    trackPairTrees[mpv]->Branch("track1NSigmaElectron", &track1_NSigmaElectron[mpv]);
    trackPairTrees[mpv]->Branch("track1charge", &track1_charge[mpv]);
    trackPairTrees[mpv]->Branch("track1_NHitsdEdx", &track1_NHitsdEdx[mpv]);
    trackPairTrees[mpv]->Branch("track1_index",&track1_index[mpv]);

    trackPairTrees[mpv]->Branch("track2pt", &track2_pt[mpv]);
    trackPairTrees[mpv]->Branch("track2eta", &track2_eta[mpv]);
    trackPairTrees[mpv]->Branch("track2phi", &track2_phi[mpv]);
    trackPairTrees[mpv]->Branch("track2NSigmaPion", &track2_NSigmaPion[mpv]);
    trackPairTrees[mpv]->Branch("track2NSigmaKaon", &track2_NSigmaKaon[mpv]);
    trackPairTrees[mpv]->Branch("track2NSigmaProton", &track2_NSigmaProton[mpv]);
    trackPairTrees[mpv]->Branch("track2NSigmaElectron", &track2_NSigmaElectron[mpv]);
    trackPairTrees[mpv]->Branch("track2charge", &track2_charge[mpv]);
    trackPairTrees[mpv]->Branch("track2_NHitsdEdx", &track2_NHitsdEdx[mpv]);
    trackPairTrees[mpv]->Branch("track2_index",&track2_index[mpv]);
  }

  // connect upc event to the chain
  static StUPCEvent *upcEvt = 0x0;
  t->SetBranchAddress("mUPCEvent", &upcEvt);
  // ask for number of events
  Long64_t nev = t->GetEntries();
  cout << "Number of events: " << nev << endl;

  // nev = 100000;
  // event loop
  TTimer timer;
  TRandom3 myRandom;
  for (Long64_t iev = 0; iev < nev; iev++) {
    if (iev % 1000 == 0) {
      long long tmp = (long long)timer.GetAbsTime();
      UInt_t seed = tmp / myRandom.Rndm();
      myRandom.SetSeed(seed);
    }
    t->GetEntry(iev);
    if (iev % 50000 == 0) {
      cout << " Procesing event " << iev / 1000 << "k ... " << endl;
      // cout << " Event has " << upcEvt->getNumberofTriggers() << endl;
    }
    hists["hEventCounter"]->Fill(0); // All events

    if (!upcEvt->isTrigger(700001))
      continue; // use only ZDC-minibias trigger (main trigger for 2019)
    hists["hEventCounter"]->Fill(1);

    int nTracks = upcEvt->getNumberOfTracks();
    int nVertices = upcEvt->getNumberOfVertices();

    hists["hNumVertices"]->Fill(nVertices);
    hists["hRefMult"]->Fill(nTracks);
    hists["hAvgMultPerVertex"]->Fill(((double)nTracks) / ((double)nVertices));

    map<int, int> multPerVertex;
    map<int, vector<StUPCTrack *>> trksForVertex;

    for (Int_t i = 0; i < nTracks; i++) {
      StUPCTrack *trk = upcEvt->getTrack(i);
      if (trk == NULL)
        continue;
      UInt_t vtxId = trk->getVertexId();
      multPerVertex[vtxId]++;
      trksForVertex[vtxId].push_back(trk);

      // all track QA
      fillTrackQA(trk, "track");
    }

    for (auto kv : multPerVertex) {
      hists["hMultPerVertex"]->Fill(kv.second);
    }

    /* if (nVertices >= 5) */
    /*   continue; */
    hists["hEventCounter"]->Fill(2);
    /* if (nTracks >= 25) */
    /*   continue; */
    hists["hEventCounter"]->Fill(3);

    TLorentzVector lv(0, 0, 0, 0);
    TLorentzVector lv1(0, 0, 0, 0);
    TLorentzVector lv2(0, 0, 0, 0);
    // make pairs in each vertex
    for (int iVtx = 0; iVtx < nVertices; iVtx++) {
      vector<StUPCTrack *> tracks = trksForVertex[iVtx];
      int mPV = multPerVertex[iVtx];
      int mPV_count = std::count(mPV_values.begin(),mPV_values.end(),mPV);
      if (mPV_count == 0)
        continue;
      hists["hEventCounter"]->Fill(4);
      for (int iTrk1 = 0; iTrk1 < tracks.size(); iTrk1++) {
        StUPCTrack *trk1 = tracks[iTrk1];

        for (int iTrk2 = iTrk1; iTrk2 < tracks.size(); iTrk2++) {
          if (iTrk1 == iTrk2)
            continue;
          StUPCTrack *trk2 = tracks[iTrk2];
          double x2kk = pow(trk1->getNSigmasTPCKaon(), 2) +
                        pow(trk2->getNSigmasTPCKaon(), 2);
          double x2pipi = pow(trk1->getNSigmasTPCPion(), 2) +
                          pow(trk2->getNSigmasTPCPion(), 2);
          double x2ee = pow(trk1->getNSigmasTPCElectron(), 2) +
                        pow(trk2->getNSigmasTPCElectron(), 2);
          double x2pp = pow(trk1->getNSigmasTPCProton(), 2) +
                        pow(trk2->getNSigmasTPCProton(), 2);

          hists["hX2kk"]->Fill(x2kk);

          TVector3 p1, p2;
          trk1->getMomentum(p1);
          trk2->getMomentum(p2);
          lv1.SetPtEtaPhiM(p1.Pt(), p1.Eta(), p1.Phi(), 0.493);
          lv2.SetPtEtaPhiM(p2.Pt(), p2.Eta(), p2.Phi(), 0.493);
          lv = lv1 + lv2;

          int chargeSum = trk1->getCharge() + trk2->getCharge();

          bool passPidSelection =
              (x2ee > 50) && (x2kk < 50) && (x2pipi > 50) && (x2pp > 50);

          bool passPairSelection = (x2kk < 800) && (lv1.Pt() < 1.) && (lv2.Pt() < 1.);

          if (passPairSelection) {
            track1_pt[mPV] = lv1.Pt();
            track1_eta[mPV] = lv1.Eta();
            track1_phi[mPV] = lv1.Phi();
            track1_NSigmaKaon[mPV] = trk1->getNSigmasTPCKaon();
            track1_NSigmaPion[mPV] = trk1->getNSigmasTPCPion();
            track1_NSigmaProton[mPV] = trk1->getNSigmasTPCProton();
            track1_NSigmaElectron[mPV] = trk1->getNSigmasTPCElectron();
            track1_charge[mPV] = trk1->getCharge();
            track1_NHitsdEdx[mPV] = trk1->getNhitsDEdx();
            track1_index[mPV] = iTrk1;
            track2_pt[mPV] = lv2.Pt();
            track2_eta[mPV] = lv2.Eta();
            track2_phi[mPV] = lv2.Phi();
            track2_NSigmaKaon[mPV] = trk2->getNSigmasTPCKaon();
            track2_NSigmaPion[mPV] = trk2->getNSigmasTPCPion();
            track2_NSigmaProton[mPV] = trk2->getNSigmasTPCProton();
            track2_NSigmaElectron[mPV] = trk2->getNSigmasTPCElectron();
            track2_charge[mPV] = trk2->getCharge();
            track2_NHitsdEdx[mPV] = trk2->getNhitsDEdx();
            track2_index[mPV] = iTrk2;
            trackPairTrees[mPV]->Fill();
          }
          if (passPidSelection) {
            hists["hEventCounter"]->Fill(5);
            if (chargeSum == 0) {
              hists["hEventCounter"]->Fill(6);
              hists["hUlsMasskk"]->Fill(lv.M());
              hists["hUlsPtkk"]->Fill(lv.Pt());
              hists["hUlsMassPtkk"]->Fill(lv.M(), lv.Pt());

              if (lv.M() > 1.0 && lv.M() < 1.04) {
                hists["hEventCounter"]->Fill(7);
                if (lv.Pt() < 0.3) {
                  hists["hEventCounter"]->Fill(8);
                }
              }
              fillTrackQA(trk1, "uls");
              fillTrackQA(trk2, "uls");
            } else if (chargeSum != 0) {
              hists["hLsMasskk"]->Fill(lv.M());
              hists["hLsPtkk"]->Fill(lv.Pt());
              hists["hLsMassPtkk"]->Fill(lv.M(), lv.Pt());

              fillTrackQA(trk1, "ls");
              fillTrackQA(trk2, "ls");

            } // like sign
          }
        }
      }
    }
  } // end of event loop

  // Write all trees
  for (auto& tree : trackPairTrees) {
    tree.second->Write();
  }

  post();
}
//////////////////////////////////////////////////////////////////////////////////////////////////////
bool passX2(double x2ee, double x2pipi) {
  if (x2ee < 10 && x2pipi > 20)
    return true;
  return false;
}
/*bool passX2( double x2ee, double x2pipi ){
         if ( x2ee < 10 && x2pipi >10 && (2 * x2ee + 10 ) < x2pipi && x2pipi <
XpipiMaxCut ) return true; return false;
}*/
/*bool passX2( double x2ee, double x2pipi ){
        if ( x2ee < XeeCut && XeeXpipi * x2ee < x2pipi && x2pipi < XpipiMaxCut )
                return true;
        return false;
}*/
//////////////////////////////////////////////////////////////////////////////////////////////////////
double getShift(double p, TH1D *histo, double pLLimit, double pHLimit) {
  if (p < pLLimit)
    return 0;
  if (p > pHLimit)
    p = pHLimit;
  int Bin = histo->GetXaxis()->FindBin(p + 1.e-8);
  if (Bin < 0)
    return 0;
  return histo->GetBinContent(Bin);
}

ClassImp( StPhiKKMaker )

StPhiKKMaker::StPhiKKMaker(){

}

StPhiKKMaker::~StPhiKKMaker(){

}

void StPhiKKMaker::Make(){
	make_ana("/work/data/UPCDstList.txt", "output.root");
}
