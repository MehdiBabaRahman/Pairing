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



int Pt_eta(){   
    
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

   


    //opening the root file 


// TString masses[10] = {"5.root", "15.root", "20.root", "25.root", "30.root", "35.root", "40.root", "45.root", "50.root", "55.root"};
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


TString savePT1 = "Figures/pt_200_";

TString saveEta1 = "Figures/eta_200_";

TString pdf = ".pdf";

TString savePT;
TString saveEta;

savePT = savePT1 + masses[ii] + pdf;
saveEta = saveEta1 + masses[ii] + pdf;






   TFile *myFile = new TFile(temp);
     
  





	 // TFile *myFile = new TFile("root_files/root_files/MZD_125/MZD_125_MSD_5.root");



    TCanvas *cnv1 = new TCanvas();
    TH1F *h1 = new TH1F(title, "gen #mu p_{T}", 100, 0, 150);
    h1->SetLineColor(1);
    h1->GetYaxis()->SetRangeUser(0,250);

    //TCanvas *cnv2 = new TCanvas();
    TH1F *h2 = new TH1F("Subleading gen di#mu p_{T} ", "Sub-leading gen di#mu p_{T}", 100, 0, 150);
    h2->SetLineColor(2);
    h2->GetYaxis()->SetRangeUser(0,250);

   ///TCanvas *cnv3 = new TCanvas();
    TH1F *h3 = new TH1F("3^{rd} gen di#mu p_{T}", "3^{rd} gen di#mu p_{T}", 100, 0, 150);
    h3->SetLineColor(3);
    h3->GetYaxis()->SetRangeUser(0,250);

    //TCanvas *cnv4 = new TCanvas();
    TH1F *h4 = new TH1F("4^{th} gen di#mu p_{T}", "4^{th} gen di#mu p_{T}", 100, 0, 150);
    h4->SetLineColor(4);
    h4->GetYaxis()->SetRangeUser(0,250);


    TCanvas *cnv5 = new TCanvas();
    TH1F *h5 = new TH1F(title, "gen #mu #eta", 100, -5, 5);
    h5->SetLineColor(1);
    h5->GetYaxis()->SetRangeUser(0,80);

    //TCanvas *cnv6 = new TCanvas();
    TH1F *h6 = new TH1F("Subleading gen di#mu #eta", "Sub-leading gen di#mu #eta", 100, -5, 5);
    h6->SetLineColor(2);
    h6->GetYaxis()->SetRangeUser(0,80);

   // TCanvas *cnv7 = new TCanvas();
    TH1F *h7 = new TH1F("3^{rd} gen di#mu #eta", "3^{rd} gen di#mu #eta", 100, -5, 5);
    h7->SetLineColor(3);
    h7->GetYaxis()->SetRangeUser(0,80);

  //  TCanvas *cnv8 = new TCanvas();
    TH1F *h8 = new TH1F("4^{th} gen di#mu #eta", "4^{th} gen di#mu #eta", 100, -5, 5);
    h8->SetLineColor(4);
    h8->GetYaxis()->SetRangeUser(0,80);



     // TH2F *h1 = new TH2F("p_{T}, Leading #mu","",100,0,70,100,0,70);
     // h2->SetMarkerColor(2);
     // h2->SetMarkerStyle(7);
     // h2->SetMarkerSize(4);
     // h2->SetXTitle("diMuonC_FittedVtx_m [GeV]");
     // h2->SetYTitle("diMuonF_FittedVtx_m [GeV]");



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

   // if(is1GenMuHighPt ){h5->Fill(genMu0_eta);} //h1->Fill(genMu0_pT) && 
   //  if(is2GenMuHighPt ){h6->Fill(genMu1_eta);} // h2->Fill(genMu1_pT)
   //    if(is3GenMuLowPt){h7->Fill(genMu2_eta);} // h3->Fill(genMu2_pT
        if(is4GenMuLowPt){
        h1->Fill(genMu0_pT);
        h2->Fill(genMu1_pT);
        h3->Fill(genMu2_pT);
        h4->Fill(genMu3_pT);
        h5->Fill(genMu0_eta);
        h6->Fill(genMu1_eta);
        h7->Fill(genMu2_eta);
        h8->Fill(genMu3_eta);} 
        }

   TLegend *leg = new TLegend(0.7,0.7,0.9,0.9);
   leg->SetHeader("Muon rank","c");
   leg->AddEntry(h1,"Leading #mu","l");
   leg->AddEntry(h2,"Sub-leading #mu","l");
   leg->AddEntry(h3,"3^{rd} #mu","l");
   leg->AddEntry(h4,"4^{th} #mu","l");


   TLegend *leg1 = new TLegend(0.7,0.7,0.9,0.9);
   leg1->SetHeader("Muon rank","c");
   leg1->AddEntry(h1,"Leading #mu","l");
   leg1->AddEntry(h2,"Sub-leading #mu","l");
   leg1->AddEntry(h3,"3^{rd} #mu","l");
   leg1->AddEntry(h4,"4^{th} #mu","l");


   
   gStyle->SetOptStat(1000000001);

    cnv1->cd();
    h1->Draw();

    //cnv2->cd();
    h2->Draw("Same");

    //cnv3->cd(); 
    h3->Draw("Same");  
   

    //cnv4->cd();   
    h4->Draw("Same");
    leg->Draw("Same");



    cnv5->cd();
    h5->Draw();

  // cnv6->cd();
    h6->Draw("Same");

   // cnv7->cd(); 
    h7->Draw("Same");  
   

  //  cnv8->cd();   
    h8->Draw("Same");
    leg1->Draw("Same");



    cnv1->SaveAs(savePT);

    cnv5->SaveAs(saveEta);

}

     return 0;


}
