//===========================================================================
//Note: ( @Wei SHI, Jun 28, 2019 )
//This program may be terminated in the case of large number of events,
//      depending on running cluster limit settings on CPU time, etc
//To check settings:
//      Bash: ulimit -H -a, ulimit -S -a; tcsh: limit -h, limit
//===========================================================================

#include <iostream>
#include "Riostream.h"
#include <string>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <math.h>
#include "TTree.h"
#include "TBranch.h"
#include "TFrame.h"
#include "TCanvas.h"
#include "TPaveLabel.h"
#include "TPaveText.h"
#include "TFile.h"
#include "TString.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1.h"
#include "TF1.h"
#include "TLorentzVector.h"
#include "math.h"
#include "time.h"
#include "TRandom.h"
#include "TSpectrum.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include <algorithm>
#include "TLine.h"



int MassWindow(){   
    
// declaring the branches as floats 


  Float_t diMuonC_FittedVtx_m;
  Float_t diMuonF_FittedVtx_m;
  Int_t run;
  Int_t lumi;
  Int_t event;

  Bool_t is1GenMuHighPt;
  Bool_t is2GenMuHighPt;
  Bool_t is3GenMuLowPt;
  Bool_t is4GenMuLowPt;
  Float_t genMu0_pT; //leading mu
  Float_t genMu1_pT;
  Float_t genMu2_pT;
  Float_t genMu3_pT;
  Float_t genMu0_eta;
  Float_t genMu1_eta;
  Float_t genMu2_eta;
  Float_t genMu3_eta;
  Float_t genMu0_phi;
  Float_t genMu1_phi;
  Float_t genMu2_phi;
  Float_t genMu3_phi;
  Float_t genA0_Lxy; //A0:dark photon that contains the most energetic muon; redt: wrt detector
  Float_t genA1_Lxy;
  Float_t genA0_Lz;
  Float_t genA1_Lz;
  Float_t genA0Mu0_pt;
  Float_t genA0Mu1_pt;
  Float_t genA1Mu0_pt;
  Float_t genA1Mu1_pt;
  Float_t genA0Mu_dR;
  Float_t genA1Mu_dR;
  Float_t genDiMu0_M;
  Float_t genDiMu1_M;

  Int_t  nRecoMu;
  Bool_t is1SelMuHighPt;
  Bool_t is2SelMuHighPt;
  Bool_t is3SelMuLowPt;
  Bool_t is4SelMuLowPt;
  Float_t selMu0_pT; //leading mu
  Float_t selMu1_pT;
  Float_t selMu2_pT;
  Float_t selMu3_pT;
  Float_t selMu0_eta;
  Float_t selMu1_eta;
  Float_t selMu2_eta;
  Float_t selMu3_eta;
  Float_t selMu0_phi;
  Float_t selMu1_phi;
  Float_t selMu2_phi;
  Float_t selMu3_phi;
  Float_t selMu0_charge;
  Float_t selMu1_charge;
  Float_t selMu2_charge;
  Float_t selMu3_charge;
  Float_t muJetC_Mu0_pt;
  Float_t muJetC_Mu1_pt;
  Float_t muJetF_Mu0_pt;
  Float_t muJetF_Mu1_pt;
  Float_t muJetC_Mu0_eta;
  Float_t muJetC_Mu1_eta;
  Float_t muJetF_Mu0_eta;
  Float_t muJetF_Mu1_eta;
  Float_t muJetC_Mu0_phi;
  Float_t muJetC_Mu1_phi;
  Float_t muJetF_Mu0_phi;
  Float_t muJetF_Mu1_phi;
  Int_t muJetC_Mu0_matched_segs;
  Int_t muJetC_Mu1_matched_segs;
  Int_t muJetF_Mu0_matched_segs;
  Int_t muJetF_Mu1_matched_segs;

  Int_t   nMuJets;
  Bool_t  is2MuJets;
  Bool_t  is2DiMuons;
  Bool_t  dimuC_Mu0_SA;
  Bool_t  dimuC_Mu1_SA;
  Bool_t  dimuF_Mu0_SA;
  Bool_t  dimuF_Mu1_SA;
  Int_t   SAmu_nTrkWP1;
  Int_t   SAmu_nTrkWP2;
  Int_t   SAmu_nTrkWP3;
  Int_t   SAmu_nTrkNodzWP1;
  Int_t   SAmu_nTrkNodzWP2;
  Int_t   SAmu_nTrkNodzWP3;
  Float_t SAmu_TrkIsoWP1;
  Float_t SAmu_TrkIsoWP2;
  Float_t SAmu_TrkIsoWP3;
  Float_t SAmu_TrkIsoNodzWP1;
  Float_t SAmu_TrkIsoNodzWP2;
  Float_t SAmu_TrkIsoNodzWP3;
  Int_t   nSAMu;
  Int_t   dimuC_nSAMu;
  Int_t   dimuF_nSAMu;
  Int_t   nNonTrackerMu;
  Int_t   dimuC_nNonTrackerMu;
  Int_t   dimuF_nNonTrackerMu;
  Float_t diMuonC_FittedVtx_prob;
  Float_t diMuonF_FittedVtx_prob;
  Float_t diMuonC_FittedVtx_dR;
  Float_t diMuonF_FittedVtx_dR;
  Float_t diMuonC_FittedVtx_Lxy;
  Float_t diMuonC_FittedVtx_L;
  Float_t diMuonF_FittedVtx_Lxy;
  Float_t diMuonF_FittedVtx_L;
  Float_t massC;
  Float_t massF;
  Int_t   NPATJet;
  Int_t   NPATJetTightB;
  Int_t   NPATJetMediumB;
  Int_t   NPATJetLooseB;

  Float_t reco4mu_m;
  Float_t recoRePaired2muleading_m;
  Float_t recoRePaired2mutrailing_m;
  Float_t recoRePaired2muleading_dR;
  Float_t recoRePaired2mutrailing_dR;

  Bool_t  isVtxOK;
  Float_t diMuons_dz_FittedVtx;
  Float_t diMuonC_FittedVtx_dz;
  Float_t diMuonF_FittedVtx_dz;
  Bool_t  is2DiMuonsMassOK;

  //Start DEBUG: many HLT paths
  Bool_t  isSignalHLT_0_Fired;
  Bool_t  isSignalHLT_1_Fired;
  Bool_t  isSignalHLT_2_Fired;
  Bool_t  isSignalHLT_3_Fired;
  //End DEBUG: many HLT paths
  Bool_t  isSignalHLTFired;
  Bool_t  isSignalHLTL1Fired;
  Bool_t  isOrthogonalHLTFired;

  Float_t diMuonC_IsoTk_FittedVtx;
  Float_t diMuonF_IsoTk_FittedVtx;
  Float_t diMuonCMu0_IsoTk0p3_FittedVtx;
  Float_t diMuonCMu0_IsoTk0p4_FittedVtx;
  Float_t diMuonCMu0_IsoTk0p5_FittedVtx;
  Float_t diMuonCMu1_IsoTk0p3_FittedVtx;
  Float_t diMuonCMu1_IsoTk0p4_FittedVtx;
  Float_t diMuonCMu1_IsoTk0p5_FittedVtx;
  Float_t diMuonFMu0_IsoTk0p3_FittedVtx;
  Float_t diMuonFMu0_IsoTk0p4_FittedVtx;
  Float_t diMuonFMu0_IsoTk0p5_FittedVtx;
  Float_t diMuonFMu1_IsoTk0p3_FittedVtx;
  Float_t diMuonFMu1_IsoTk0p4_FittedVtx;
  Float_t diMuonFMu1_IsoTk0p5_FittedVtx;

  Int_t  diMuonC_m1_FittedVtx_hitpix_Phase1;
  Int_t  diMuonC_m2_FittedVtx_hitpix_Phase1;
  Int_t  diMuonF_m1_FittedVtx_hitpix_Phase1;
  Int_t  diMuonF_m2_FittedVtx_hitpix_Phase1;

  Bool_t orph_passOffLineSelPtEta;
  Bool_t orph_isSignalHLTFired;
  Bool_t orph_isVertexOK;
  Int_t orph_dimu_Mu0_hitpix_Phase1;
  Int_t orph_dimu_Mu1_hitpix_Phase1;
  Float_t orph_dimu_mass;
  Float_t orph_dimu_isoTk;
  Float_t orph_dimu_Mu0_isoTk0p3;
  Float_t orph_dimu_Mu0_isoTk0p4;
  Float_t orph_dimu_Mu0_isoTk0p5;
  Float_t orph_dimu_Mu1_isoTk0p3;
  Float_t orph_dimu_Mu1_isoTk0p4;
  Float_t orph_dimu_Mu1_isoTk0p5;
  Float_t orph_dimu_z;
  Float_t orph_isoTk;

   




TString masses[10] = {"5", "15", "20", "25", "30", "35", "40", "45", "50", "55"};

for (int ii = 0; ii < 10; ii++){


TString path = "MZD_200_";
TString root = ".root";

TString title1 = "MZ_{D}=200 GeV,";
TString title2 = "Mf_{D_{1}}=";
TString title3="Gev";

TString title;
title = title1 + title2 + masses[ii] + title3;


TString temp;

temp = path + masses[ii] + root;


TString saveMW1 = "Figures/MassW_200_";


TString pdf = ".pdf";

TString saveMW;

saveMW = saveMW1 + masses[ii] + pdf;



   TFile *myFile = new TFile(temp);








    //opening the root file 

	// TFile *myFile = new TFile("MZD_200_40.root");

	TCanvas *cnv1 = new TCanvas();


    // TPad *pad1 = new TPad("pad1","",0,0,1,1);

     TH2F *h2 = new TH2F(title,"",100,0,80,100,0,90);
     h2->SetMarkerColor(2);
     h2->SetMarkerStyle(7);
     h2->SetMarkerSize(4);
     h2->SetXTitle("diMuonC_FittedVtx_m [GeV]");
     h2->SetYTitle("diMuonF_FittedVtx_m [GeV]");



    gDirectory->cd("cutFlowAnalyzerPXBL4PXFL3;1");
    gDirectory->pwd();
    TTree *myTree1 = nullptr;
    gDirectory->GetObject("Events;1",myTree1);

    //Number of entries
    int N = myTree1->GetEntries();
    cout << "Number of entries: " << N << endl;

   
    
// declaring branch addresses 

  int nentries;//entries in main tree
  int mentries;//entries in orphan tree

  myTree1->SetBranchAddress("run",   &run);
  myTree1->SetBranchAddress("lumi",  &lumi);
  myTree1->SetBranchAddress("event", &event);

  myTree1->SetBranchAddress("is1GenMuHighPt",&is1GenMuHighPt);
  myTree1->SetBranchAddress("is2GenMuHighPt",&is2GenMuHighPt);
  myTree1->SetBranchAddress("is3GenMuLowPt", &is3GenMuLowPt);
  myTree1->SetBranchAddress("is4GenMuLowPt", &is4GenMuLowPt);
  myTree1->SetBranchAddress("genMu0_pT",     &genMu0_pT); //leading mu
  myTree1->SetBranchAddress("genMu1_pT",     &genMu1_pT);
  myTree1->SetBranchAddress("genMu2_pT",     &genMu2_pT);
  myTree1->SetBranchAddress("genMu3_pT",     &genMu3_pT);
  myTree1->SetBranchAddress("genMu0_eta",    &genMu0_eta);
  myTree1->SetBranchAddress("genMu1_eta",    &genMu1_eta);
  myTree1->SetBranchAddress("genMu2_eta",    &genMu2_eta);
  myTree1->SetBranchAddress("genMu3_eta",    &genMu3_eta);
  myTree1->SetBranchAddress("genMu0_phi",    &genMu0_phi);
  myTree1->SetBranchAddress("genMu1_phi",    &genMu1_phi);
  myTree1->SetBranchAddress("genMu2_phi",    &genMu2_phi);
  myTree1->SetBranchAddress("genMu3_phi",    &genMu3_phi);
  myTree1->SetBranchAddress("genA0_Lxy",     &genA0_Lxy);
  myTree1->SetBranchAddress("genA1_Lxy",     &genA1_Lxy);
  myTree1->SetBranchAddress("genA0_Lz",      &genA0_Lz);
  myTree1->SetBranchAddress("genA1_Lz",      &genA1_Lz);
  myTree1->SetBranchAddress("genA0Mu0_pt",   &genA0Mu0_pt);
  myTree1->SetBranchAddress("genA0Mu1_pt",   &genA0Mu1_pt);
  myTree1->SetBranchAddress("genA1Mu0_pt",   &genA1Mu0_pt);
  myTree1->SetBranchAddress("genA1Mu1_pt",   &genA1Mu1_pt);
  myTree1->SetBranchAddress("genA0Mu_dR",    &genA0Mu_dR);
  myTree1->SetBranchAddress("genA1Mu_dR",    &genA1Mu_dR);
  myTree1->SetBranchAddress("genDiMu0_M",    &genDiMu0_M);
  myTree1->SetBranchAddress("genDiMu1_M",    &genDiMu1_M);

  myTree1->SetBranchAddress("nRecoMu",       &nRecoMu);
  myTree1->SetBranchAddress("is1SelMuHighPt",&is1SelMuHighPt);
  myTree1->SetBranchAddress("is2SelMuHighPt",&is2SelMuHighPt);
  myTree1->SetBranchAddress("is3SelMuLowPt", &is3SelMuLowPt);
  myTree1->SetBranchAddress("is4SelMuLowPt", &is4SelMuLowPt);
  myTree1->SetBranchAddress("selMu0_pT",     &selMu0_pT); //leading mu
  myTree1->SetBranchAddress("selMu1_pT",     &selMu1_pT);
  myTree1->SetBranchAddress("selMu2_pT",     &selMu2_pT);
  myTree1->SetBranchAddress("selMu3_pT",     &selMu3_pT);
  myTree1->SetBranchAddress("selMu0_eta",    &selMu0_eta);
  myTree1->SetBranchAddress("selMu1_eta",    &selMu1_eta);
  myTree1->SetBranchAddress("selMu2_eta",    &selMu2_eta);
  myTree1->SetBranchAddress("selMu3_eta",    &selMu3_eta);
  myTree1->SetBranchAddress("selMu0_phi",    &selMu0_phi);
  myTree1->SetBranchAddress("selMu1_phi",    &selMu1_phi);
  myTree1->SetBranchAddress("selMu2_phi",    &selMu2_phi);
  myTree1->SetBranchAddress("selMu3_phi",    &selMu3_phi);
  myTree1->SetBranchAddress("selMu0_charge", &selMu0_charge);
  myTree1->SetBranchAddress("selMu1_charge", &selMu1_charge);
  myTree1->SetBranchAddress("selMu2_charge", &selMu2_charge);
  myTree1->SetBranchAddress("selMu3_charge", &selMu3_charge);
  myTree1->SetBranchAddress("muJetC_Mu0_pt", &muJetC_Mu0_pt);
  myTree1->SetBranchAddress("muJetC_Mu1_pt", &muJetC_Mu1_pt);
  myTree1->SetBranchAddress("muJetF_Mu0_pt", &muJetF_Mu0_pt);
  myTree1->SetBranchAddress("muJetF_Mu1_pt", &muJetF_Mu1_pt);
  myTree1->SetBranchAddress("muJetC_Mu0_eta", &muJetC_Mu0_eta);
  myTree1->SetBranchAddress("muJetC_Mu1_eta", &muJetC_Mu1_eta);
  myTree1->SetBranchAddress("muJetF_Mu0_eta", &muJetF_Mu0_eta);
  myTree1->SetBranchAddress("muJetF_Mu1_eta", &muJetF_Mu1_eta);
  myTree1->SetBranchAddress("muJetC_Mu0_phi", &muJetC_Mu0_phi);
  myTree1->SetBranchAddress("muJetC_Mu1_phi", &muJetC_Mu1_phi);
  myTree1->SetBranchAddress("muJetF_Mu0_phi", &muJetF_Mu0_phi);
  myTree1->SetBranchAddress("muJetF_Mu1_phi", &muJetF_Mu1_phi);
  myTree1->SetBranchAddress("muJetC_Mu0_matched_segs", &muJetC_Mu0_matched_segs);
  myTree1->SetBranchAddress("muJetC_Mu1_matched_segs", &muJetC_Mu1_matched_segs);
  myTree1->SetBranchAddress("muJetF_Mu0_matched_segs", &muJetF_Mu0_matched_segs);
  myTree1->SetBranchAddress("muJetF_Mu1_matched_segs", &muJetF_Mu1_matched_segs);

  myTree1->SetBranchAddress("isVertexOK", &isVtxOK);
  myTree1->SetBranchAddress("nMuJets",    &nMuJets);
  myTree1->SetBranchAddress("is2MuJets",  &is2MuJets);
  myTree1->SetBranchAddress("is2DiMuons", &is2DiMuons);
  myTree1->SetBranchAddress("nSAMu",      &nSAMu);
  myTree1->SetBranchAddress("dimuC_nSAMu",&dimuC_nSAMu);
  myTree1->SetBranchAddress("dimuF_nSAMu",&dimuF_nSAMu);
  myTree1->SetBranchAddress("nNonTrackerMu", &nNonTrackerMu);
  myTree1->SetBranchAddress("dimuC_nNonTrackerMu", &dimuC_nNonTrackerMu);
  myTree1->SetBranchAddress("dimuF_nNonTrackerMu", &dimuF_nNonTrackerMu);
  myTree1->SetBranchAddress("dimuC_Mu0_SA", &dimuC_Mu0_SA);
  myTree1->SetBranchAddress("dimuC_Mu1_SA", &dimuC_Mu1_SA);
  myTree1->SetBranchAddress("dimuF_Mu0_SA", &dimuF_Mu0_SA);
  myTree1->SetBranchAddress("dimuF_Mu1_SA", &dimuF_Mu1_SA);
  myTree1->SetBranchAddress("SAmu_nTrkWP1", &SAmu_nTrkWP1);
  myTree1->SetBranchAddress("SAmu_nTrkWP2", &SAmu_nTrkWP2);
  myTree1->SetBranchAddress("SAmu_nTrkWP3", &SAmu_nTrkWP3);
  myTree1->SetBranchAddress("SAmu_nTrkNodzWP1", &SAmu_nTrkNodzWP1);
  myTree1->SetBranchAddress("SAmu_nTrkNodzWP2", &SAmu_nTrkNodzWP2);
  myTree1->SetBranchAddress("SAmu_nTrkNodzWP3", &SAmu_nTrkNodzWP3);
  myTree1->SetBranchAddress("SAmu_TrkIsoWP1", &SAmu_TrkIsoWP1);
  myTree1->SetBranchAddress("SAmu_TrkIsoWP2", &SAmu_TrkIsoWP2);
  myTree1->SetBranchAddress("SAmu_TrkIsoWP3", &SAmu_TrkIsoWP3);
  myTree1->SetBranchAddress("SAmu_TrkIsoNodzWP1", &SAmu_TrkIsoNodzWP1);
  myTree1->SetBranchAddress("SAmu_TrkIsoNodzWP2", &SAmu_TrkIsoNodzWP2);
  myTree1->SetBranchAddress("SAmu_TrkIsoNodzWP3", &SAmu_TrkIsoNodzWP3);
  myTree1->SetBranchAddress("diMuonC_FittedVtx_prob",&diMuonC_FittedVtx_prob);
  myTree1->SetBranchAddress("diMuonF_FittedVtx_prob",&diMuonF_FittedVtx_prob);
  myTree1->SetBranchAddress("diMuonC_FittedVtx_dR",  &diMuonC_FittedVtx_dR);
  myTree1->SetBranchAddress("diMuonF_FittedVtx_dR",  &diMuonF_FittedVtx_dR);
  myTree1->SetBranchAddress("diMuonC_FittedVtx_Lxy", &diMuonC_FittedVtx_Lxy);
  myTree1->SetBranchAddress("diMuonC_FittedVtx_L",   &diMuonC_FittedVtx_L);
  myTree1->SetBranchAddress("diMuonF_FittedVtx_Lxy", &diMuonF_FittedVtx_Lxy);
  myTree1->SetBranchAddress("diMuonF_FittedVtx_L",   &diMuonF_FittedVtx_L);
  //Need to use fitted vertex mass, not the muon pair mass
  //In most cases, they are the close, but not necessarily in some cases
  myTree1->SetBranchAddress("diMuonC_FittedVtx_m", &massC);
  myTree1->SetBranchAddress("diMuonF_FittedVtx_m", &massF);
  myTree1->SetBranchAddress("NPATJet",       &NPATJet);
  myTree1->SetBranchAddress("NPATJetTightB", &NPATJetTightB);
  myTree1->SetBranchAddress("NPATJetMediumB",&NPATJetMediumB);
  myTree1->SetBranchAddress("NPATJetLooseB", &NPATJetLooseB);

  myTree1->SetBranchAddress("reco4mu_m",                  &reco4mu_m);
  myTree1->SetBranchAddress("recoRePaired2muleading_m",   &recoRePaired2muleading_m);
  myTree1->SetBranchAddress("recoRePaired2mutrailing_m",  &recoRePaired2mutrailing_m);
  myTree1->SetBranchAddress("recoRePaired2muleading_dR",  &recoRePaired2muleading_dR);
  myTree1->SetBranchAddress("recoRePaired2mutrailing_dR", &recoRePaired2mutrailing_dR);

  myTree1->SetBranchAddress("diMuonC_m1_FittedVtx_hitpix_Phase1", &diMuonC_m1_FittedVtx_hitpix_Phase1);
  myTree1->SetBranchAddress("diMuonC_m2_FittedVtx_hitpix_Phase1", &diMuonC_m2_FittedVtx_hitpix_Phase1);
  myTree1->SetBranchAddress("diMuonF_m1_FittedVtx_hitpix_Phase1", &diMuonF_m1_FittedVtx_hitpix_Phase1);
  myTree1->SetBranchAddress("diMuonF_m2_FittedVtx_hitpix_Phase1", &diMuonF_m2_FittedVtx_hitpix_Phase1);

  myTree1->SetBranchAddress("diMuonC_FittedVtx_dz",          &diMuonC_FittedVtx_dz);
  myTree1->SetBranchAddress("diMuonF_FittedVtx_dz",          &diMuonF_FittedVtx_dz);
  myTree1->SetBranchAddress("diMuons_dz_FittedVtx",          &diMuons_dz_FittedVtx);
  myTree1->SetBranchAddress("diMuonC_IsoTk_FittedVtx",       &diMuonC_IsoTk_FittedVtx);
  myTree1->SetBranchAddress("diMuonF_IsoTk_FittedVtx",       &diMuonF_IsoTk_FittedVtx);
  myTree1->SetBranchAddress("diMuonCMu0_IsoTk0p3_FittedVtx", &diMuonCMu0_IsoTk0p3_FittedVtx);
  myTree1->SetBranchAddress("diMuonCMu0_IsoTk0p4_FittedVtx", &diMuonCMu0_IsoTk0p4_FittedVtx);
  myTree1->SetBranchAddress("diMuonCMu0_IsoTk0p5_FittedVtx", &diMuonCMu0_IsoTk0p5_FittedVtx);
  myTree1->SetBranchAddress("diMuonCMu1_IsoTk0p3_FittedVtx", &diMuonCMu1_IsoTk0p3_FittedVtx);
  myTree1->SetBranchAddress("diMuonCMu1_IsoTk0p4_FittedVtx", &diMuonCMu1_IsoTk0p4_FittedVtx);
  myTree1->SetBranchAddress("diMuonCMu1_IsoTk0p5_FittedVtx", &diMuonCMu1_IsoTk0p5_FittedVtx);
  myTree1->SetBranchAddress("diMuonFMu0_IsoTk0p3_FittedVtx", &diMuonFMu0_IsoTk0p3_FittedVtx);
  myTree1->SetBranchAddress("diMuonFMu0_IsoTk0p4_FittedVtx", &diMuonFMu0_IsoTk0p4_FittedVtx);
  myTree1->SetBranchAddress("diMuonFMu0_IsoTk0p5_FittedVtx", &diMuonFMu0_IsoTk0p5_FittedVtx);
  myTree1->SetBranchAddress("diMuonFMu1_IsoTk0p3_FittedVtx", &diMuonFMu1_IsoTk0p3_FittedVtx);
  myTree1->SetBranchAddress("diMuonFMu1_IsoTk0p4_FittedVtx", &diMuonFMu1_IsoTk0p4_FittedVtx);
  myTree1->SetBranchAddress("diMuonFMu1_IsoTk0p5_FittedVtx", &diMuonFMu1_IsoTk0p5_FittedVtx);
  //Start Debug: many HLT paths
  myTree1->SetBranchAddress("isSignalHLT_0_Fired",           &isSignalHLT_0_Fired);
  myTree1->SetBranchAddress("isSignalHLT_1_Fired",           &isSignalHLT_1_Fired);
  myTree1->SetBranchAddress("isSignalHLT_2_Fired",           &isSignalHLT_2_Fired);
  myTree1->SetBranchAddress("isSignalHLT_3_Fired",           &isSignalHLT_3_Fired);
  //End Debug: many HLT paths
  myTree1->SetBranchAddress("isSignalHLTFired",              &isSignalHLTFired);
  myTree1->SetBranchAddress("isSignalHLTL1Fired",            &isSignalHLTL1Fired);
  myTree1->SetBranchAddress("isOrthogonalHLTFired",          &isOrthogonalHLTFired);
  myTree1->SetBranchAddress("is2DiMuonsMassOK_FittedVtx",    &is2DiMuonsMassOK);

  Float_t R0 = 10.0;
  Float_t P0 = 0.2;
  Float_t L0 = 0.1;
  Float_t C0 = 2.0;
  //f(dR)-Poly4
  Float_t p0 = 8.53647;
  Float_t p1 = -50.4571;
  Float_t p2 = 109.83;
  Float_t p3 = -92.7445;
  Float_t p4 = 36.8351;


  //Get branch from orphan-dimuon tree
  // o->SetBranchAddress("orph_passOffLineSelPtEta",    &orph_passOffLineSelPtEta); //offline pT, eta selection same as signal
  // o->SetBranchAddress("orph_isSignalHLTFired",       &orph_isSignalHLTFired);
  // o->SetBranchAddress("orph_isVertexOK",             &orph_isVertexOK);
  // o->SetBranchAddress("orph_dimu_Mu0_hitpix_Phase1", &orph_dimu_Mu0_hitpix_Phase1);
  // o->SetBranchAddress("orph_dimu_Mu1_hitpix_Phase1", &orph_dimu_Mu1_hitpix_Phase1);
  // o->SetBranchAddress("orph_dimu_mass",              &orph_dimu_mass);
  // o->SetBranchAddress("orph_dimu_z",                 &orph_dimu_z);

  // o->SetBranchAddress("orph_dimu_isoTk",        &orph_dimu_isoTk);
  // o->SetBranchAddress("orph_dimu_Mu0_isoTk0p3", &orph_dimu_Mu0_isoTk0p3);
  // o->SetBranchAddress("orph_dimu_Mu0_isoTk0p4", &orph_dimu_Mu0_isoTk0p4);
  // o->SetBranchAddress("orph_dimu_Mu0_isoTk0p5", &orph_dimu_Mu0_isoTk0p5);
  // o->SetBranchAddress("orph_dimu_Mu1_isoTk0p3", &orph_dimu_Mu1_isoTk0p3);
  // o->SetBranchAddress("orph_dimu_Mu1_isoTk0p4", &orph_dimu_Mu1_isoTk0p4);
  // o->SetBranchAddress("orph_dimu_Mu1_isoTk0p5", &orph_dimu_Mu1_isoTk0p5);
  // o->SetBranchAddress("orph_isoTk",             &orph_isoTk);




// looping over all events 

for(int ii = 0; ii < N; ii++){
    myTree1->GetEntry(ii);

   if(is1GenMuHighPt == true && is2GenMuHighPt == true && is3GenMuLowPt == true && is4GenMuLowPt == true){
   	if ( ( genA0_Lxy < 16.0 && fabs(genA0_Lz) < 51.6 ) && ( genA1_Lxy < 16.0 && fabs(genA1_Lz) < 51.6 ) ){
   		if ( is1SelMuHighPt== true && is2SelMuHighPt== true && is3SelMuLowPt== true && is4SelMuLowPt== true){
   			if ( isVtxOK ){
   				if ( is2DiMuons && nSAMu <= 1 && diMuonC_FittedVtx_prob > P0*(1 - dimuC_nSAMu)*exp( -( p0 + p1*(sqrt(diMuonC_FittedVtx_dR)) + p2*pow(sqrt(diMuonC_FittedVtx_dR), 2) + p3*pow(sqrt(diMuonC_FittedVtx_dR), 3) + p4*pow(sqrt(diMuonC_FittedVtx_dR), 4) )*pow(fabs(diMuonC_FittedVtx_Lxy/R0), C0) ) &&
               diMuonF_FittedVtx_prob > P0*(1 - dimuF_nSAMu)*exp( -( p0 + p1*(sqrt(diMuonF_FittedVtx_dR)) + p2*pow(sqrt(diMuonF_FittedVtx_dR), 2) + p3*pow(sqrt(diMuonF_FittedVtx_dR), 3) + p4*pow(sqrt(diMuonF_FittedVtx_dR), 4) )*pow(fabs(diMuonF_FittedVtx_Lxy/R0), C0) ) ) {
   		 			if ( ( diMuonC_m1_FittedVtx_hitpix_Phase1 == 1 || diMuonC_m2_FittedVtx_hitpix_Phase1 == 1 ) && ( diMuonF_m1_FittedVtx_hitpix_Phase1 == 1 || diMuonF_m2_FittedVtx_hitpix_Phase1 == 1 ) ) {
   		 				if ( (massC > 11 && massC < 60 && massF > 11 && massF < 60 && (recoRePaired2mutrailing_dR >= 0.2 || recoRePaired2mutrailing_m >= 3) ) ||
                          (massC > 0.2113 && massC < 9 && massF > 0.2113 && massF < 9) ) {
   		 					if ( diMuonC_IsoTk_FittedVtx < 2.3 && diMuonF_IsoTk_FittedVtx < 2.3 ) {
   		 						if ( isSignalHLTFired ){
   		 							if ( (massC > 0.2113 && massC < 9 && massF > 0.2113 && massF < 9) ||
                                (massC > 11 && massC < 60 && massF > 11 && massF < 60 && (nSAMu == 0 || ( nSAMu == 1 && (diMuonC_FittedVtx_Lxy > L0 || diMuonF_FittedVtx_Lxy > L0) && ( (dimuC_Mu0_SA==1 && muJetC_Mu0_matched_segs>=2) || (dimuC_Mu1_SA==1 && muJetC_Mu1_matched_segs>=2) || (dimuF_Mu0_SA==1 && muJetF_Mu0_matched_segs>=2) || (dimuF_Mu1_SA==1 && muJetF_Mu1_matched_segs>=2) ) ) ) ) ) {


   										h2->Fill(massC,massF);


   									}

   								}
   							}
   						}
    				}
    			}
    		}
    	}
    }
 }
}
  h2->Draw("");


  double mean[11] = {0.2560, 0.4012, 0.7003, 1.0000, 1.9990, 4.9980, 8.4920, 14.990, 24.980, 34.970, 57.930};
  double window[11] = {0.0438535344, 0.0336164980, 0.0359654996, 0.0428032000, 0.0753576680, 0.2037460866, 0.3539943472, 0.7041010200, 1.1972469080, 1.6131510600, 3.9866777100};


  double m2Small[11];
  double m2Large[11];
  for (int i = 0; i < 11; i++) {
    m2Small[i] = mean[i] - window[i];
    m2Large[i] = mean[i] + window[i];
  }
  TGraph* corridorDn = new TGraph(11, mean, m2Small);
  TGraph* corridorUp = new TGraph(11, mean, m2Large);
 

  corridorDn->SetLineColor(1); corridorDn->SetLineStyle(9); corridorDn->SetLineWidth(2);
  corridorUp->SetLineColor(1); corridorUp->SetLineStyle(9); corridorUp->SetLineWidth(2); 
  corridorDn->Draw("L Same"); corridorUp->Draw("L Same");


cnv1->SaveAs(saveMW);    

}
 return 0;

}











