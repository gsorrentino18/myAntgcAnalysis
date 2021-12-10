#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <iostream>
#include <fstream>
#include <THStack.h>
#include "CMS_lumi.C"
#include <TStyle.h>

using namespace std;

void Validation(){

   Bool_t is2017 = 1;
   Bool_t is2018 = 0;

   Bool_t isEE = 0;

   Bool_t do_cumulative = 0; //cumulative does not work here

   const Int_t nfile=99;
   Int_t nfile_eff = 99;
   TFile *file[nfile];

   TString path = "/local/cms/user/gsorrent/antgc_analysis/BDTValidation/UL17/EBNtuples/merged/";
   TString savePath = "/local/cms/user/gsorrent/antgc_analysis/BDTValidation/UL17/EBPlots/";
   if (isEE) {
      path = "/local/cms/user/gsorrent/antgc_analysis/BDTValidation/UL17/EENtuples/merged/";
      savePath = "/local/cms/user/gsorrent/antgc_analysis/BDTValidation/UL17/EEPlots/";
   }

   if (is2018) {
      path = "/local/cms/user/gsorrent/antgc_analysis/BDTValidation/UL18/EBNtuples/merged/";
      savePath = "/local/cms/user/gsorrent/antgc_analysis/BDTValidation/UL18/EBPlots/";
      if (isEE) {
        path = "/local/cms/user/gsorrent/antgc_analysis/BDTValidation/UL18/EENtuples/merged/";
        savePath = "/local/cms/user/gsorrent/antgc_analysis/BDTValidation/UL18/EEPlots/";
      }
   }

   ofstream fraction_file;
   fraction_file.open(savePath+"bkg_fractions.txt");

   if (is2017) {

      nfile_eff=32;
      file[0] = new TFile(path+"data.root");
      file[1] = new TFile(path+"DYJetsToLLM50HT100to200TuneCP5PSweights13TeVmadgraphMLMpythia8RunIISummer20UL17MiniAODv2106X.root");
      file[2] = new TFile(path+"DYJetsToLLM50HT1200to2500TuneCP5PSweights13TeVmadgraphMLMpythia8RunIISummer20UL17MiniAODv2106X.root");
      file[3] = new TFile(path+"DYJetsToLLM50HT200to400TuneCP5PSweights13TeVmadgraphMLMpythia8RunIISummer20UL17MiniAODv2106X.root");
      file[4] = new TFile(path+"DYJetsToLLM50HT2500toInfTuneCP5PSweights13TeVmadgraphMLMpythia8RunIISummer20UL17MiniAODv2106X.root");
      file[5] = new TFile(path+"DYJetsToLLM50HT400to600TuneCP5PSweights13TeVmadgraphMLMpythia8RunIISummer20UL17MiniAODv2106X.root");
      file[6] = new TFile(path+"DYJetsToLLM50HT600to800TuneCP5PSweights13TeVmadgraphMLMpythia8RunIISummer20UL17MiniAODv2106X.root");
      file[7] = new TFile(path+"DYJetsToLLM50HT70to100TuneCP5PSweights13TeVmadgraphMLMpythia8RunIISummer20UL17MiniAODv2106X.root");
      file[8] = new TFile(path+"DYJetsToLLM50HT800to1200TuneCP5PSweights13TeVmadgraphMLMpythia8RunIISummer20UL17MiniAODv2106X.root");
      file[9] = new TFile(path+"GJetsHT200To400TuneCP513TeVmadgraphMLMpythia8RunIISummer19UL17MiniAODv2106X.root");
      file[10] = new TFile(path+"GJetsHT400To600TuneCP513TeVmadgraphMLMpythia8RunIISummer19UL17MiniAODv2106X.root");
      file[11] = new TFile(path+"GJetsHT600ToInfTuneCP513TeVmadgraphMLMpythia8RunIISummer19UL17MiniAODv2106X.root");
      file[12] = new TFile(path+"QCDHT1000to1500TuneCP5PSWeights13TeVmadgraphMLMpythia8RunIISummer19UL17MiniAODv2106X.root");
      file[13] = new TFile(path+"QCDHT1500to2000TuneCP5PSWeights13TeVmadgraphMLMpythia8RunIISummer19UL17MiniAODv2106X.root");
      file[14] = new TFile(path+"QCDHT2000toInfTuneCP5PSWeights13TeVmadgraphMLMpythia8RunIISummer19UL17MiniAODv2106X.root");
      file[15] = new TFile(path+"QCDHT200to300TuneCP5PSWeights13TeVmadgraphMLMpythia8RunIISummer19UL17MiniAODv2106X.root");
      file[16] = new TFile(path+"QCDHT300to500TuneCP5PSWeights13TeVmadgraphMLMpythia8RunIISummer19UL17MiniAODv2106X.root");
      file[17] = new TFile(path+"QCDHT500to700TuneCP5PSWeights13TeVmadgraphMLMpythia8RunIISummer19UL17MiniAODv2106X.root");
      file[18] = new TFile(path+"QCDHT700to1000TuneCP5PSWeights13TeVmadgraphMLMpythia8RunIISummer19UL17MiniAODv2106X.root");
      file[19] = new TFile(path+"WGToLNuG01J5fPtG130TuneCP513TeVamcatnloFXFXpythia8RunIISummer20UL17MiniAODv2106X.root");
      file[20] = new TFile(path+"TGJetsTuneCP513TeVamcatnlomadspinpythia8RunIISummer19UL17MiniAODv2106X.root");
      file[21] = new TFile(path+"TTGJetsTuneCP513TeVamcatnloFXFXmadspinpythia8RunIISummer20UL17MiniAODv2106X.root");
      file[22] = new TFile(path+"WWTuneCP513TeVpythia8RunIISummer20UL17MiniAOD106X.root");
      file[23] = new TFile(path+"WZTuneCP513TeVpythia8RunIISummer20UL17MiniAODv2106X.root");
      file[24] = new TFile(path+"ZZTuneCP513TeVpythia8RunIISummer20UL17MiniAODv2106X.root");
      file[25] = new TFile(path+"WJetsToLNuHT100To200TuneCP513TeVmadgraphMLMpythia8RunIISummer20UL17MiniAODv2106X.root");
      file[26] = new TFile(path+"WJetsToLNuHT1200To2500TuneCP513TeVmadgraphMLMpythia8RunIISummer20UL17MiniAODv2106X.root");
      file[27] = new TFile(path+"WJetsToLNuHT200To400TuneCP513TeVmadgraphMLMpythia8RunIISummer20UL17MiniAODv2106X.root");
      file[28] = new TFile(path+"WJetsToLNuHT400To600TuneCP513TeVmadgraphMLMpythia8RunIISummer20UL17MiniAODv2106X.root");
      file[29] = new TFile(path+"WJetsToLNuHT600To800TuneCP513TeVmadgraphMLMpythia8RunIISummer20UL17MiniAODv2106X.root");
      file[30] = new TFile(path+"WJetsToLNuHT800To1200TuneCP513TeVmadgraphMLMpythia8RunIISummer20UL17MiniAOD106X.root");
      file[31] = new TFile(path+"DiPhotonJetsMGG80toInfTuneCP513TeVamcatnloFXFXpythia8RunIISummer20UL17MiniAOD106X.root");

   }

   if (is2018) {

     nfile_eff=32;
     file[0] = new TFile(path+"data.root");
     file[1] = new TFile(path+"DYJetsToLLM50HT100to200TuneCP5PSweights13TeVmadgraphMLMpythia8RunIISummer20UL18MiniAOD106X.root");     
     file[2] = new TFile(path+"DYJetsToLLM50HT1200to2500TuneCP5PSweights13TeVmadgraphMLMpythia8RunIISummer20UL18MiniAOD106X.root");
     file[3] = new TFile(path+"DYJetsToLLM50HT200to400TuneCP5PSweights13TeVmadgraphMLMpythia8RunIISummer20UL18MiniAOD106X.root");
     file[4] = new TFile(path+"DYJetsToLLM50HT2500toInfTuneCP5PSweights13TeVmadgraphMLMpythia8RunIISummer20UL18MiniAODv2106X.root");
     file[5] = new TFile(path+"DYJetsToLLM50HT400to600TuneCP5PSweights13TeVmadgraphMLMpythia8RunIISummer20UL18MiniAODv2106X.root");
     file[6] = new TFile(path+"DYJetsToLLM50HT600to800TuneCP5PSweights13TeVmadgraphMLMpythia8RunIISummer20UL18MiniAOD106X.root");
     file[7] = new TFile(path+"DYJetsToLLM50HT70to100TuneCP5PSweights13TeVmadgraphMLMpythia8RunIISummer20UL18MiniAOD106X.root");
     file[8] = new TFile(path+"DYJetsToLLM50HT800to1200TuneCP5PSweights13TeVmadgraphMLMpythia8RunIISummer20UL18MiniAOD106X.root");
     file[9] = new TFile(path+"GJetsHT200To400TuneCP513TeVmadgraphMLMpythia8RunIISummer19UL18MiniAOD106X.root");
     file[10] = new TFile(path+"GJetsHT400To600TuneCP513TeVmadgraphMLMpythia8RunIISummer19UL18MiniAOD106X.root");
     file[11] = new TFile(path+"GJetsHT600ToInfTuneCP513TeVmadgraphMLMpythia8RunIISummer19UL18MiniAOD106X.root");
     file[12] = new TFile(path+"QCDHT1000to1500TuneCP5PSWeights13TeVmadgraphMLMpythia8RunIISummer19UL18MiniAOD106X.root");
     file[13] = new TFile(path+"QCDHT1500to2000TuneCP5PSWeights13TeVmadgraphMLMpythia8RunIISummer19UL18MiniAOD106X.root");
     file[14] = new TFile(path+"QCDHT2000toInfTuneCP5PSWeights13TeVmadgraphMLMpythia8RunIISummer19UL18MiniAOD106X.root");
     file[15] = new TFile(path+"QCDHT200to300TuneCP5PSWeights13TeVmadgraphMLMpythia8RunIISummer19UL18MiniAOD106X.root");
     file[16] = new TFile(path+"QCDHT300to500TuneCP5PSWeights13TeVmadgraphMLMpythia8RunIISummer19UL18MiniAOD106X.root");
     file[17] = new TFile(path+"QCDHT500to700TuneCP5PSWeights13TeVmadgraphMLMpythia8RunIISummer19UL18MiniAOD106X.root");
     file[18] = new TFile(path+"QCDHT700to1000TuneCP5PSWeights13TeVmadgraphMLMpythia8RunIISummer19UL18MiniAOD106X.root");
     file[19] = new TFile(path+"WGToLNuG01J5fPtG130TuneCP513TeVamcatnloFXFXpythia8RunIISummer20UL18MiniAODv2106X.root");
     file[20] = new TFile(path+"TGJetsTuneCP513TeVamcatnlomadspinpythia8RunIISummer19UL18MiniAOD106X.root");
     file[21] = new TFile(path+"TTGJetsTuneCP513TeVamcatnloFXFXmadspinpythia8RunIISummer20UL18MiniAOD106X.root");
     file[22] = new TFile(path+"WWTuneCP513TeVpythia8RunIISummer20UL18MiniAODv2106X.root");
     file[23] = new TFile(path+"WZTuneCP513TeVpythia8RunIISummer20UL18MiniAODv2106X.root");
     file[24] = new TFile(path+"ZZTuneCP513TeVpythia8RunIISummer20UL18MiniAOD106X.root");
     file[25] = new TFile(path+"WJetsToLNuHT100To200TuneCP513TeVmadgraphMLMpythia8RunIISummer20UL18MiniAOD106X.root");
     file[26] = new TFile(path+"WJetsToLNuHT1200To2500TuneCP513TeVmadgraphMLMpythia8RunIISummer20UL18MiniAOD106X.root");
     file[27] = new TFile(path+"WJetsToLNuHT200To400TuneCP513TeVmadgraphMLMpythia8RunIISummer20UL18MiniAOD106X.root");
     file[28] = new TFile(path+"WJetsToLNuHT400To600TuneCP513TeVmadgraphMLMpythia8RunIISummer20UL18MiniAOD106X.root");
     file[29] = new TFile(path+"WJetsToLNuHT600To800TuneCP513TeVmadgraphMLMpythia8RunIISummer20UL18MiniAOD106X.root");
     file[30] = new TFile(path+"WJetsToLNuHT70To100TuneCP513TeVmadgraphMLMpythia8RunIISummer20UL18MiniAOD106X.root");
     file[31] = new TFile(path+"WJetsToLNuHT800To1200TuneCP513TeVmadgraphMLMpythia8RunIISummer20UL18MiniAODv2106X.root");
     file[32] = new TFile(path+"DiPhotonJetsMGG80toInfTuneCP513TeVamcatnloFXFXpythia8RunIISummer20UL17MiniAOD106X.root");
   }

   Float_t puWeight;
   Float_t ph_PFECALClusIsoCorr;
   Float_t ph_PFHCALClusIsoCorr;
   Float_t ph_TrkIsoCorr;
   Float_t ph_HoverE;
   Float_t ph_BDTpred;
   Float_t ph_et;
   Float_t ph_eta;
   Float_t ph_phi;
   Float_t ph_E2x2Full5x5;
   Float_t ph_R9Full5x5;
   Float_t ph_S4Full5x5;
   Float_t ph_E2ndOESCrFull5x5;
   Float_t ph_E1x3OESCrFull5x5;
   Float_t ph_E2x5OESCrFull5x5;
   Float_t ph_E5x5OESCrFull5x5;
   Float_t ph_2x2OE3x3Full5x5;
   Float_t ph_EmaxOESCrFull5x5;
   Float_t ph_SigmaIEtaIEta;
   Float_t ph_SigmaIEtaIPhi;
   Float_t ph_SigmaIPhiIPhi;
   Float_t ph_SieieOSipipFull5x5;
   Float_t ph_EtaWidth;
   Float_t ph_PhiWidth;
   Float_t ph_EtaWOPhiWFull5x5;
   Float_t ph_ESEffSigmaRR;
   Float_t elePt;
   Float_t eleEta;
   Double_t zmass;
   Bool_t pass95;
   Int_t doubleCounting;
   Char_t ph_QualityBits;

   TTree *tree[nfile];
   for (Int_t i=0; i<nfile_eff;i++) {
     tree[i] = (TTree*)file[i]->Get("tnpPhoIDs/fitter_tree");
     tree[i] = (TTree*)file[i]->Get("tnpPhoIDs/fitter_tree");
     tree[i]->SetBranchAddress("puWeight",&puWeight);
     tree[i]->SetBranchAddress("ph_PFECALClusIsoCorr",&ph_PFECALClusIsoCorr);
     tree[i]->SetBranchAddress("ph_PFHCALClusIsoCorr",&ph_PFHCALClusIsoCorr);
     tree[i]->SetBranchAddress("ph_TkrIsoCorr",&ph_TrkIsoCorr);
     tree[i]->SetBranchAddress("ph_hoe",&ph_HoverE);
     tree[i]->SetBranchAddress("ph_BDTpred",&ph_BDTpred);
     tree[i]->SetBranchAddress("ph_et",&ph_et);
     tree[i]->SetBranchAddress("ph_eta",&ph_eta);
     tree[i]->SetBranchAddress("ph_phi",&ph_phi);
     tree[i]->SetBranchAddress("ph_E2x2Full5x5",&ph_E2x2Full5x5);
     tree[i]->SetBranchAddress("ph_R9Full5x5", &ph_R9Full5x5);
     tree[i]->SetBranchAddress("ph_S4Full5x5", &ph_S4Full5x5);
     tree[i]->SetBranchAddress("ph_EmaxOESCrFull5x5", &ph_EmaxOESCrFull5x5);
     tree[i]->SetBranchAddress("ph_E2ndOESCrFull5x5", &ph_E2ndOESCrFull5x5);
     tree[i]->SetBranchAddress("ph_E1x3OESCrFull5x5", &ph_E1x3OESCrFull5x5);
     tree[i]->SetBranchAddress("ph_E2x5OESCrFull5x5", &ph_E2x5OESCrFull5x5);
     tree[i]->SetBranchAddress("ph_E5x5OESCrFull5x5", &ph_E5x5OESCrFull5x5);
     tree[i]->SetBranchAddress("ph_2x2OE3x3Full5x5", &ph_2x2OE3x3Full5x5);
     tree[i]->SetBranchAddress("ph_EmaxOESCrFull5x5", &ph_EmaxOESCrFull5x5);
     tree[i]->SetBranchAddress("ph_sieie", &ph_SigmaIEtaIEta);
     tree[i]->SetBranchAddress("ph_sieip", &ph_SigmaIEtaIPhi);
     tree[i]->SetBranchAddress("ph_sipip", &ph_SigmaIPhiIPhi);
     tree[i]->SetBranchAddress("ph_sieieOsipip", &ph_SieieOSipipFull5x5);
     tree[i]->SetBranchAddress("ph_se", &ph_EtaWidth);
     tree[i]->SetBranchAddress("ph_sp", &ph_PhiWidth);
     tree[i]->SetBranchAddress("ph_seOsp", &ph_EtaWOPhiWFull5x5);
     tree[i]->SetBranchAddress("tag_Ele_pt", &elePt);
     tree[i]->SetBranchAddress("tag_Ele_eta", &eleEta);
     tree[i]->SetBranchAddress("pair_mass", &zmass);
     tree[i]->SetBranchAddress("pass95", &pass95);
     tree[i]->SetBranchAddress("doubleCounting", &doubleCounting);
     tree[i]->SetBranchAddress("ph_QualityBits", &ph_QualityBits);
   }

   const Int_t nhist=24;
   TH1F* hist[nfile][nhist];
   THStack *hs[nhist];

   gStyle->SetOptStat(false);

   TString sample[nfile];
   for (Int_t i=0; i<nfile_eff;i++) {
     TString num; 
     num.Form("%d",i); 
     sample[i] = i;
   }

   for (Int_t f=0; f<nfile_eff;f++) {
       hist[f][0] = new TH1F("phoBDT"+sample[f],"phoBDT"+sample[f],25,0,1);
       hist[f][1] = new TH1F("phoEt"+sample[f],"phoEt"+sample[f],20,100,900);
       hist[f][2] = new TH1F("phoEta"+sample[f],"phoEta"+sample[f],30,-1.5,1.5);
       hist[f][3] = new TH1F("ph_PFECALClusIsoCorr"+sample[f],"ph_PFECALClusIsoCorr"+sample[f],25,-10,190);
       hist[f][4] = new TH1F("ph_PFHCALClusIsoCorr"+sample[f],"ph_PFHCALClusIsoCorr"+sample[f],25,-40,160);
       hist[f][5] = new TH1F("ph_TkrIsoCorr"+sample[f],"ph_TkrIsoCorr"+sample[f],20,-5,115);
       hist[f][6] = new TH1F("ph_hoe"+sample[f],"ph_hoe"+sample[f],20,0,0.05);
       hist[f][7] = new TH1F("ph_sieie"+sample[f],"ph_sieie"+sample[f],25,0,0.025);
       hist[f][8] = new TH1F("ph_sieip"+sample[f],"ph_sieip"+sample[f],25,-0.000025,0.00025);
       hist[f][9] = new TH1F("ph_sipip"+sample[f],"ph_sipip"+sample[f],30,0.005,0.025);
       hist[f][10] = new TH1F("ph_2x2OE3x3Full5x5"+sample[f],"ph_2x2OE3x3Full5x5"+sample[f],30,0.7,1);
       hist[f][11] = new TH1F("ph_E1x3OESCrFull5x5"+sample[f],"ph_E1x3OESCrFull5x5"+sample[f],20,0.2,1);
       hist[f][12] = new TH1F("ph_E2ndOESCrFull5x5"+sample[f],"ph_E2ndOESCrFull5x5"+sample[f],25,0,0.5);
       hist[f][13] = new TH1F("ph_E2x5OESCrFull5x5"+sample[f],"ph_E2x5OESCrFull5x5"+sample[f],30,0.6,1.2);
       hist[f][14] = new TH1F("ph_EmaxOESCrFull5x5"+sample[f],"ph_EmaxOESCrFull5x5"+sample[f],20,0.1,0.9);
       hist[f][15] = new TH1F("ph_EtaWidth"+sample[f],"ph_EtaWidth"+sample[f],40,0,0.06);
       hist[f][16] = new TH1F("ph_PhiWidth"+sample[f],"ph_PhiWidth"+sample[f],40,0,0.06);
       hist[f][17] = new TH1F("ph_EtaWOPhiWFull5x5"+sample[f],"ph_EtaWOPhiWFull5x5"+sample[f],20,0,2);
       hist[f][18] = new TH1F("ph_R9Full5x5"+sample[f],"ph_R9Full5x5"+sample[f],40,0.4,1.2);
       hist[f][19] = new TH1F("ph_S4Full5x5"+sample[f],"ph_S4Full5x5"+sample[f],40,0.4,1.2);
       hist[f][20] = new TH1F("ph_SigmaIEtaIEta_o_ph_SigmaIPhiIPhi"+sample[f],"ph_SigmaIEtaIEta_o_ph_SigmaIPhiIPhi"+sample[f],35,0,1.4);
       hist[f][21] = new TH1F("zmass"+sample[f],"zmass"+sample[f],30,80,110);
       hist[f][22] = new TH1F("ele_pt"+sample[f],"ele_pt"+sample[f],35,0,350);
       hist[f][23] = new TH1F("ele_eta"+sample[f],"ele_eta"+sample[f],30,-2.5,2.5);
   }

   if (isEE) {
      for (Int_t f=0; f<nfile_eff;f++) {
        hist[f][0] = new TH1F("phoBDT"+sample[f],"phoBDT"+sample[f],25,0,1);
        hist[f][1] = new TH1F("phoEt"+sample[f],"phoEt"+sample[f],20,100,900);
        hist[f][2] = new TH1F("phoEta"+sample[f],"phoEta"+sample[f],50,-2.5,2.5);
        hist[f][3] = new TH1F("ph_PFECALClusIsoCorr"+sample[f],"ph_PFECALClusIsoCorr"+sample[f],20,-100,20);
        hist[f][4] = new TH1F("ph_PFHCALClusIsoCorr"+sample[f],"ph_PFHCALClusIsoCorr"+sample[f],20,-20,140);
        hist[f][5] = new TH1F("ph_TkrIsoCorr"+sample[f],"ph_TkrIsoCorr"+sample[f],20,-5,115);
        hist[f][6] = new TH1F("ph_hoe"+sample[f],"ph_hoe"+sample[f],20,0,0.05);
        hist[f][7] = new TH1F("ph_sieie"+sample[f],"ph_sieie"+sample[f],30,0.02,0.04);
        hist[f][8] = new TH1F("ph_sieip"+sample[f],"ph_sieip"+sample[f],40,-0.0004,0.0004);
        hist[f][9] = new TH1F("ph_sipip"+sample[f],"ph_sipip"+sample[f],40,0.01,0.05);
        hist[f][10] = new TH1F("ph_2x2OE3x3Full5x5"+sample[f],"ph_2x2OE3x3Full5x5"+sample[f],30,0.7,1);
        hist[f][11] = new TH1F("ph_E1x3OESCrFull5x5"+sample[f],"ph_E1x3OESCrFull5x5"+sample[f],20,0.2,1);
        hist[f][12] = new TH1F("ph_E2ndOESCrFull5x5"+sample[f],"ph_E2ndOESCrFull5x5"+sample[f],25,0,0.5);
        hist[f][13] = new TH1F("ph_E2x5OESCrFull5x5"+sample[f],"ph_E2x5OESCrFull5x5"+sample[f],30,0.6,1.2);
        hist[f][14] = new TH1F("ph_EmaxOESCrFull5x5"+sample[f],"ph_EmaxOESCrFull5x5"+sample[f],20,0.1,0.9);
        hist[f][15] = new TH1F("ph_EtaWidth"+sample[f],"ph_EtaWidth"+sample[f],40,0,0.06);
        hist[f][16] = new TH1F("ph_PhiWidth"+sample[f],"ph_PhiWidth"+sample[f],40,0,0.06);
        hist[f][17] = new TH1F("ph_EtaWOPhiWFull5x5"+sample[f],"ph_EtaWOPhiWFull5x5"+sample[f],20,0,2);
        hist[f][18] = new TH1F("ph_R9Full5x5"+sample[f],"ph_R9Full5x5"+sample[f],40,0.4,1.2);
        hist[f][19] = new TH1F("ph_S4Full5x5"+sample[f],"ph_S4Full5x5"+sample[f],40,0.4,1.2);
        hist[f][20] = new TH1F("ph_SigmaIEtaIEta_o_ph_SigmaIPhiIPhi"+sample[f],"ph_SigmaIEtaIEta_o_ph_SigmaIPhiIPhi"+sample[f],35,0,1.4);
        hist[f][21] = new TH1F("zmass"+sample[f],"zmass"+sample[f],30,80,110);
        hist[f][22] = new TH1F("ele_pt"+sample[f],"ele_pt"+sample[f],35,0,350);
        hist[f][23] = new TH1F("ele_eta"+sample[f],"ele_eta"+sample[f],30,-2.5,2.5); 
     }
   }

   Int_t nentries[nfile] = {0};
   for (Int_t i=0; i<nfile_eff;i++) {
     nentries[i] = (Int_t)tree[i]->GetEntries();
   }

   for (Int_t fil=0; fil<nfile_eff; fil++) {
    TH1F* cutFlowGenWeight = (TH1F*)file[fil]->Get("tnpPhoIDs/tnpPhoIDscutFlowGenWeight");
    Double_t sumGenWeight = cutFlowGenWeight->GetBinContent(1);

     for (Int_t i=0; i<nentries[fil]; i++) {
        tree[fil]->GetEntry(i);

        if (is2017) {
           if ((fil > 11 && fil <= 18) && (doubleCounting == 1)) continue; //remove double counting QCD
           //if ((fil > 25) && (doubleCounting == 1)) continue; //remove double counting W+jets

        }
        if (is2018) { 
           if ((fil > 11 && fil <= 18) && (doubleCounting == 1)) continue; //remove double counting QCD
        }

        if (ph_et < 225) continue;
        if (!isEE) {
          if (ph_TrkIsoCorr > 18) continue;
          if (ph_PFHCALClusIsoCorr > 10) continue;
        }

        //Fill histograms here 
        if (fil==0) {
          hist[fil][0]->Fill(ph_BDTpred);
          hist[fil][1]->Fill(ph_et);
          hist[fil][2]->Fill(ph_eta);
          hist[fil][3]->Fill(ph_PFECALClusIsoCorr);
          hist[fil][4]->Fill(ph_PFHCALClusIsoCorr);
          hist[fil][5]->Fill(ph_TrkIsoCorr);
          hist[fil][6]->Fill(ph_HoverE);
          hist[fil][7]->Fill(ph_SigmaIEtaIEta);
          hist[fil][8]->Fill(ph_SigmaIEtaIPhi);
          hist[fil][9]->Fill(ph_SigmaIPhiIPhi);
          hist[fil][10]->Fill(ph_2x2OE3x3Full5x5);
          hist[fil][11]->Fill(ph_E1x3OESCrFull5x5);
          hist[fil][12]->Fill(ph_E2ndOESCrFull5x5);
          hist[fil][13]->Fill(ph_E2x5OESCrFull5x5);
          hist[fil][14]->Fill(ph_EmaxOESCrFull5x5);
          hist[fil][15]->Fill(ph_EtaWidth);
          hist[fil][16]->Fill(ph_PhiWidth);
          hist[fil][17]->Fill(ph_EtaWOPhiWFull5x5);
          hist[fil][18]->Fill(ph_R9Full5x5);
          hist[fil][19]->Fill(ph_S4Full5x5);
          hist[fil][20]->Fill(ph_SigmaIEtaIEta/ph_SigmaIPhiIPhi);
          hist[fil][21]->Fill(zmass);
          hist[fil][22]->Fill(elePt);
          hist[fil][23]->Fill(eleEta);
        } else if (fil!=0) {
          hist[fil][0]->Fill(ph_BDTpred, puWeight/sumGenWeight);
          hist[fil][1]->Fill(ph_et, puWeight/sumGenWeight);
          hist[fil][2]->Fill(ph_eta, puWeight/sumGenWeight);
          hist[fil][3]->Fill(ph_PFECALClusIsoCorr, puWeight/sumGenWeight);
          hist[fil][4]->Fill(ph_PFHCALClusIsoCorr, puWeight/sumGenWeight);
          hist[fil][5]->Fill(ph_TrkIsoCorr, puWeight/sumGenWeight);
          hist[fil][6]->Fill(ph_HoverE, puWeight/sumGenWeight);
          hist[fil][7]->Fill(ph_SigmaIEtaIEta, puWeight/sumGenWeight);
          hist[fil][8]->Fill(ph_SigmaIEtaIPhi, puWeight/sumGenWeight);
          hist[fil][9]->Fill(ph_SigmaIPhiIPhi, puWeight/sumGenWeight);
          hist[fil][10]->Fill(ph_2x2OE3x3Full5x5, puWeight/sumGenWeight);
          hist[fil][11]->Fill(ph_E1x3OESCrFull5x5, puWeight/sumGenWeight);
          hist[fil][12]->Fill(ph_E2ndOESCrFull5x5, puWeight/sumGenWeight);
          hist[fil][13]->Fill(ph_E2x5OESCrFull5x5, puWeight/sumGenWeight);
          hist[fil][14]->Fill(ph_EmaxOESCrFull5x5, puWeight/sumGenWeight);
          hist[fil][15]->Fill(ph_EtaWidth, puWeight/sumGenWeight);
          hist[fil][16]->Fill(ph_PhiWidth, puWeight/sumGenWeight);
          hist[fil][17]->Fill(ph_EtaWOPhiWFull5x5, puWeight/sumGenWeight);
          hist[fil][18]->Fill(ph_R9Full5x5, puWeight/sumGenWeight);
          hist[fil][19]->Fill(ph_S4Full5x5, puWeight/sumGenWeight);
          hist[fil][20]->Fill(ph_SigmaIEtaIEta/ph_SigmaIPhiIPhi, puWeight/sumGenWeight);
          hist[fil][21]->Fill(zmass, puWeight/sumGenWeight);
          hist[fil][22]->Fill(elePt, puWeight/sumGenWeight);
          hist[fil][23]->Fill(eleEta, puWeight/sumGenWeight);
        }
     }
   }

   TString title[nhist] = {"ph_BDTpred", "ph_et", "ph_eta", "ph_PFECALClusIsoCorr", "ph_PFHCALClusIsoCorr", "ph_TkrIsoCorr", "ph_hoe", "ph_sieie", "ph_sieip", "ph_sipip", "ph_2x2OE3x3Full5x5", "ph_E1x3OESCrFull5x5", "ph_E2ndOESCrFull5x5", "ph_E2x5OESCrFull5x5", "ph_EmaxOESCrFull5x5", "ph_EtaWidth", "ph_PhiWidth", "ph_EtaWOPhiWFull5x5", "ph_R9Full5x5", "ph_S4Full5x5", "ph_SigmaIEtaIEta_o_ph_SigmaIPhiIPhi", "zmass", "ele_pt", "ele_eta"};
   TString title_axis[nhist] = {"BDT score", "p_{T}^{#gamma}", "#eta^{#gamma}", "PF ECALClusIsoCorr", "PF HCALClusIsoCorr", "TkrIsoCorr", "H/E", "#sigma_{i#eta i #eta}", "#sigma_{i#eta i#phi}","#sigma_{i#phi i #phi}", "E_{2x2}/E_{3x3}","E_{1x3}/E^{raw}_{SC}","E_{2}/E^{raw}_{SC}","E_{2x5}/E^{raw}_{SC}","E_{max}/E^{raw}_{SC}","#sigma_{#eta}", "#sigma_{#phi}", "#sigma_{#eta}/#sigma_{#phi}", "S4","R9", "#sigma_{i#eta i #eta}/#sigma_{i#phi i #phi}", "e#gamma invariant mass", "p_{T}^{e}", "#eta^{e}"};

   TH1F* h_ratio[nhist];
   TH1F* h_all[nhist];
   TH1F* h_allScaled[nhist];

   TH1F* h_dy[nhist];
   TH1F* h_gjets[nhist];
   TH1F* h_qcd[nhist];
   TH1F* h_wglnu[nhist];
   TH1F* h_tg[nhist];
   TH1F* h_ww[nhist];
   TH1F* h_wz[nhist];
   TH1F* h_zz[nhist];  
   TH1F* h_wjets[nhist];

   Double_t overflow[nfile][nhist]; 
   Double_t underflow[nfile][nhist];

   for (Int_t fil=0; fil<nfile_eff; fil++) {
      for (int i=0; i<nhist; i++) {
         overflow[fil][i] = hist[fil][i]->GetBinContent(hist[fil][i]->GetNbinsX()+1);
         underflow[fil][i] = hist[fil][i]->GetBinContent(0);
         hist[fil][i]->SetBinContent(1, hist[fil][i]->GetBinContent(1)+underflow[fil][i]);
         hist[fil][i]->SetBinContent(hist[fil][i]->GetNbinsX(), hist[fil][i]->GetBinContent(hist[fil][i]->GetNbinsX())+overflow[fil][i]);
      }
   }

   for (int i=0; i<nhist; i++) {

      if (is2017) {
         h_dy[i] = (TH1F*)hist[1][i]->Clone("h_dy");
         h_gjets[i] = (TH1F*)hist[9][i]->Clone("h_gjets");
         h_qcd[i] = (TH1F*)hist[12][i]->Clone("h_qcd");
         h_wglnu[i] = (TH1F*)hist[19][i]->Clone("h_wglnu");
         h_tg[i] = (TH1F*)hist[20][i]->Clone("h_tg");
         h_ww[i] = (TH1F*)hist[22][i]->Clone("h_ww");
         h_wz[i] = (TH1F*)hist[23][i]->Clone("h_wz");
         h_zz[i] = (TH1F*)hist[24][i]->Clone("h_zz");
         h_wjets[i] = (TH1F*)hist[25][i]->Clone("h_wjets");

         for (int j=1; j<nfile_eff; j++) {
           if (j>1 && j<=8) {h_dy[i]->Add(hist[j][i]); }
           if (j > 8 && j <= 11) {h_gjets[i]->Add(hist[j][i]);}//GJets 
           if (j > 12 && j <= 18) { h_qcd[i]->Add(hist[j][i]);}//QCD
           if (j > 20 && j <= 21) { h_tg[i]->Add(hist[j][i]);}//TG, TTG
           if (j > 25 && j <= 30) { h_wjets[i]->Add(hist[j][i]);}//WJets
         }
      }
      if (is2018) {
         h_dy[i] = (TH1F*)hist[1][i]->Clone("h_dy");
         h_gjets[i] = (TH1F*)hist[9][i]->Clone("h_gjets");
         h_qcd[i] = (TH1F*)hist[12][i]->Clone("h_qcd");
         h_wglnu[i] = (TH1F*)hist[19][i]->Clone("h_wglnu");
         h_tg[i] = (TH1F*)hist[20][i]->Clone("h_tg");
         h_ww[i] = (TH1F*)hist[22][i]->Clone("h_ww");
         h_wz[i] = (TH1F*)hist[23][i]->Clone("h_wz");
         h_zz[i] = (TH1F*)hist[24][i]->Clone("h_zz");
         h_wjets[i] = (TH1F*)hist[25][i]->Clone("h_wjets");

         for (int j=1; j<nfile_eff; j++) {
           if (j>1 && j<=8) {h_dy[i]->Add(hist[j][i]); }
           if (j > 8 && j <= 11) {h_gjets[i]->Add(hist[j][i]);}//GJets 
           if (j > 12 && j <= 18) { h_qcd[i]->Add(hist[j][i]);}//QCD
           if (j > 20 && j <= 21) { h_tg[i]->Add(hist[j][i]);}//TG, TTG
           if (j > 25 && j <= 31) { h_wjets[i]->Add(hist[j][i]);}//WJets
         } 
      }

      h_all[i] = (TH1F*)hist[1][i]->Clone("h_all");
      for (int j=2; j<nfile_eff; j++) { //sum signal+bkgs correctly normalized
        h_all[i]->Add(hist[j][i]);
      }

      Double_t scale = hist[0][i]->Integral()/h_all[i]->Integral(); //scale for data
      for (int j=1; j<nfile_eff; j++) {
        hist[j][i]->Sumw2(1); 
        hist[j][i]->Scale(scale);
      }

      h_allScaled[i] = (TH1F*)hist[1][i]->Clone("h_allScaled");
      for (int j=2; j<nfile_eff; j++) { //sum signal+bkgs correctly normalized
        h_allScaled[i]->Add(hist[j][i]);
      }
      h_ratio[i] = (TH1F*)hist[0][i]->Clone("h_ratio");
      h_ratio[i]->Sumw2(1);
      h_ratio[i]->Divide(h_allScaled[i]);
      TCanvas *c1 = new TCanvas("c1", "c1");
      c1->cd();

      TPad* pad1 = new TPad("pad1", "pad1", 0.0, 0.3, 1.0, 1.0);
      pad1->SetBottomMargin(0);
      pad1->Draw();
      pad1->cd();

      hs[i] = new THStack(title[i]+"_stack", "");

      hist[0][i]->SetMarkerStyle(20); //data
      hist[0][i]->SetTitle("");

      if (is2017) {
         for (int j=1; j<nfile_eff; j++) {
           hist[j][i]->SetTitle("");

           if (j <= 8) { //DY
             hist[j][i]->SetFillColor(kOrange-4);
             hist[j][i]->SetLineColor(kOrange-4);
           }
           h_dy[i]->SetFillColor(kOrange-4);
           h_dy[i]->SetLineColor(kOrange-4);
             
           if (j > 8 && j <= 11) { //GJets
           hist[j][i]->SetFillColor(kMagenta-8);
           hist[j][i]->SetLineColor(kMagenta-8);
           }

           if (j > 11 && j <= 18) { //QCD
           hist[j][i]->SetFillColor(kPink+2);
           hist[j][i]->SetLineColor(kPink+2);
           }

           if (j == 19) { //WGLnu
           hist[j][i]->SetFillColor(kCyan-2);
           hist[j][i]->SetLineColor(kCyan-2);
           }

           if (j > 19 && j <= 21) { //TG, TTG
           hist[j][i]->SetFillColor(kBlue-9);
           hist[j][i]->SetLineColor(kBlue-9);
           }

           if (j == 22) { //WW
           hist[j][i]->SetFillColor(kTeal-6); //kBlue-7
           hist[j][i]->SetLineColor(kTeal-6); //kBlue-7
           }

           if (j == 23) { //WZ
           hist[j][i]->SetFillColor(kSpring+2);
           hist[j][i]->SetLineColor(kSpring+2);
           }

           if (j == 24) { //ZZ
           hist[j][i]->SetFillColor(kSpring+4);
           hist[j][i]->SetLineColor(kSpring+4);
           }

           if (j > 24 && j < 31) { //Wjets
           hist[j][i]->SetFillColor(kAzure-5);
           hist[j][i]->SetLineColor(kAzure-5);
           }

           if (j == 31) { //Diphoton
           hist[j][i]->SetFillColor(kBlue-7);
           hist[j][i]->SetLineColor(kBlue-7);
           }

           if (j > 8)  hs[i]->Add(hist[j][i]);
         }
      }

      if (is2018) {
         for (int j=1; j<nfile_eff; j++) {
           hist[j][i]->SetTitle("");

           if (j <= 8) { //DY
             hist[j][i]->SetFillColor(kOrange-4);
             hist[j][i]->SetLineColor(kOrange-4);
           }
           h_dy[i]->SetFillColor(kOrange-4);
           h_dy[i]->SetLineColor(kOrange-4);

           if (j > 8 && j <= 11) { //GJets
           hist[j][i]->SetFillColor(kMagenta-8);
           hist[j][i]->SetLineColor(kMagenta-8);
           }

           if (j > 11 && j <= 18) { //QCD
           hist[j][i]->SetFillColor(kPink+2);
           hist[j][i]->SetLineColor(kPink+2);
           }

           if (j == 19) { //WGLnu
           hist[j][i]->SetFillColor(kCyan-2);
           hist[j][i]->SetLineColor(kCyan-2);
           }

           if (j > 19 && j <= 21) { //TG, TTG
           hist[j][i]->SetFillColor(kBlue-9);
           hist[j][i]->SetLineColor(kBlue-9);
           }

           if (j == 22) { //WW
           hist[j][i]->SetFillColor(kTeal-6); //kBlue-7
           hist[j][i]->SetLineColor(kTeal-6); //kBlue-7
           }
            
           if (j == 23) { //WZ
           hist[j][i]->SetFillColor(kSpring+2);
           hist[j][i]->SetLineColor(kSpring+2);
           }

           if (j == 24) { //ZZ
           hist[j][i]->SetFillColor(kSpring+4);
           hist[j][i]->SetLineColor(kSpring+4);
           }

           if (j > 24) { //Wjets
           hist[j][i]->SetFillColor(kAzure-5);
           hist[j][i]->SetLineColor(kAzure-5);
           }

           if (j > 8)  hs[i]->Add(hist[j][i]);
         }
      }

      hs[i]->Add(hist[1][i]);
      hs[i]->Add(hist[2][i]);
      hs[i]->Add(hist[3][i]);
      hs[i]->Add(hist[4][i]);
      hs[i]->Add(hist[5][i]);
      hs[i]->Add(hist[6][i]);
      hs[i]->Add(hist[7][i]);
      hs[i]->Add(hist[8][i]);

      hist[0][i]->Draw("ep");
      hs[i]->Draw("hist");
      h_allScaled[i]->SetFillStyle(3001);
      h_allScaled[i]->SetFillColor(15);
      h_allScaled[i]->GetYaxis()->SetTitle("Events");
      h_allScaled[i]->Draw("e2, SAME");
      hist[0][i]->Draw("ep, SAME");
      hs[i]->SetTitle("");
      hs[i]->GetYaxis()->SetTitle("Events");

      TLegend* leg = new TLegend();
      leg = new TLegend(0.68, 0.5, 0.9, 0.891);
      leg->SetBorderSize(0);
      leg->SetEntrySeparation(0.01);
      leg->SetFillColor(0);
      leg->SetFillStyle(0);

      if (is2017) {
         leg->AddEntry(hist[0][i], "Data", "ep");
         leg->AddEntry(h_allScaled[i], "MC stat", "f");
         leg->AddEntry(hist[1][i], "Z(#rightarrow ee)+jets", "f");
         leg->AddEntry(hist[11][i], "#gamma +jets", "f");
         leg->AddEntry(hist[18][i], "QCD", "f");
         leg->AddEntry(hist[21][i], "t#gamma, tt#gamma", "f");
         leg->AddEntry(hist[19][i], "W(#rightarrow l#nu)+#gamma", "f");
         leg->AddEntry(hist[22][i], "WW", "f");
         leg->AddEntry(hist[23][i], "WZ", "f");
         leg->AddEntry(hist[24][i], "ZZ", "f");
         leg->AddEntry(hist[25][i], "WJets", "f");
         leg->AddEntry(hist[31][i], "#gamma #gamma", "f");
      }

      if (is2018) {
         leg->AddEntry(hist[0][i], "Data", "ep");
         leg->AddEntry(h_allScaled[i], "MC stat", "f");
         leg->AddEntry(hist[1][i], "DYToEE", "f");
         leg->AddEntry(hist[11][i], "GJets", "f");
         leg->AddEntry(hist[18][i], "QCD", "f");
         leg->AddEntry(hist[21][i], "TTG,TG", "f");
         leg->AddEntry(hist[19][i], "WGToLNu", "f");
         leg->AddEntry(hist[22][i], "WW", "f");
         leg->AddEntry(hist[23][i], "WZ", "f");
         leg->AddEntry(hist[24][i], "ZZ", "f");
         leg->AddEntry(hist[25][i], "WJets", "f");
      }

      c1->Update();
      leg->Draw("SAME");

      pad1->SetLogy();

      Int_t year = 17;
      if (is2018) {year=18;}
      CMS_lumi(pad1, year, 0);

      pad1->Update();
      c1->Update();
      c1->cd();

      TPad* pad2 = new TPad("pad2", "pad2", 0.0, 0.0, 1.0, 0.3);
      pad2->SetTopMargin(0.01);
      pad2->SetBottomMargin(0.2);
      pad2->Draw();
      pad2->cd();

      h_ratio[i]->SetTitle("");
      h_ratio[i]->GetXaxis()->SetTitle(title_axis[i]);
      h_ratio[i]->GetXaxis()->SetTitleFont(42);
      h_ratio[i]->GetXaxis()->SetTitleSize(0.11);
      h_ratio[i]->GetXaxis()->SetTitleOffset(0.8);
      h_ratio[i]->GetXaxis()->SetLabelFont(42);
      h_ratio[i]->GetXaxis()->SetLabelSize(0.1);
      h_ratio[i]->GetYaxis()->SetTitle("Data/MC");
      h_ratio[i]->GetYaxis()->SetTitleSize(0.11);
      h_ratio[i]->GetYaxis()->SetTitleOffset(0.43);
      h_ratio[i]->GetYaxis()->SetLabelSize(0.1);
      h_ratio[i]->GetYaxis()->SetLabelOffset(0.01);
      h_ratio[i]->GetYaxis()->SetNdivisions(505);
      h_ratio[i]->GetYaxis()->SetRangeUser(0.5, 1.5);
      h_ratio[i]->SetMarkerStyle(20);

      TLine* line = new TLine(h_ratio[i]->GetXaxis()->GetXmin(), 1.0, h_ratio[i]->GetXaxis()->GetXmax(), 1.0);
      line->SetLineColor(kRed);
      line->SetLineWidth(1);
      h_ratio[i]->Draw("ep");
      line->Draw("SAME");
      c1->Update();
      c1->SaveAs(savePath+title[i]+"stack.png");
   } //nhist loop


   //background fraction in BDT
   const int num_bkg = 8;
   Double_t bkg_frac[num_bkg]={0};
   Double_t bkg[num_bkg]={0};
   Double_t dy=0;
   TString bkg_name[num_bkg] = {"gjets/dy", "qcd/dy", "wglnu/dy", "ttgjets+tgjets/dy", "ww/dy", "wz/dy", "zz/dy", "wjets/dy"};
   for (Int_t ibin=1; ibin<=h_dy[0]->GetNbinsX() ; ibin++) {
     dy += h_dy[0]->GetBinContent(ibin);
     bkg[0] += h_gjets[0]->GetBinContent(ibin);
     bkg[1] += h_qcd[0]->GetBinContent(ibin);
     bkg[2] += h_wglnu[0]->GetBinContent(ibin);
     bkg[3] += h_tg[0]->GetBinContent(ibin);
     bkg[4] += h_ww[0]->GetBinContent(ibin);
     bkg[5] += h_wz[0]->GetBinContent(ibin);
     bkg[6] += h_zz[0]->GetBinContent(ibin);
     bkg[7] += h_wjets[0]->GetBinContent(ibin);
   }
   for (Int_t ibkg=0; ibkg<num_bkg; ibkg++) {
     if (bkg[ibkg] < 0) bkg[ibkg] = 0;
     bkg_frac[ibkg] = bkg[ibkg]/dy;
     fraction_file << bkg_name[ibkg] << "  "<< bkg_frac[ibkg]*100. << endl;
   }
      

   if (do_cumulative) {

           TH1F* h_all = (TH1F*)h_dy[0]->Clone("h_all");
           h_all->Add(h_gjets[0]);
           h_all->Add(h_qcd[0]);
           h_all->Add(h_tg[0]);
           h_all->Add(h_wjets[0]);
           h_all->Add(h_ww[0]);
           h_all->Add(h_wz[0]);
           h_all->Add(h_zz[0]);
           h_all->Add(h_wglnu[0]);

           Double_t bsum_data = 0;
           Double_t bsum_dy = 0;
           Double_t bsum_gjets = 0;
           Double_t bsum_qcd = 0;
           Double_t bsum_tg = 0;
           Double_t bsum_wjets = 0;
           Double_t bsum_ww = 0;
           Double_t bsum_wz = 0;
           Double_t bsum_zz = 0;
           Double_t bsum_wglnu = 0;

           for (int ibin=1; ibin<=hist[0][0]->GetNbinsX(); ibin++) {
              bsum_data += hist[0][0]->GetBinContent(ibin);
              bsum_dy += h_dy[0]->GetBinContent(ibin);
              bsum_gjets += h_gjets[0]->GetBinContent(ibin);
              bsum_qcd += h_qcd[0]->GetBinContent(ibin);
              bsum_tg += h_tg[0]->GetBinContent(ibin);
              bsum_wjets += h_wjets[0]->GetBinContent(ibin);
              bsum_ww += h_ww[0]->GetBinContent(ibin);
              bsum_wz += h_wz[0]->GetBinContent(ibin);
              bsum_zz += h_zz[0]->GetBinContent(ibin);
              bsum_wglnu += h_wglnu[0]->GetBinContent(ibin);
              hist[0][0]->SetBinContent(ibin, bsum_data);
              h_dy[0]->SetBinContent(ibin, bsum_dy);
              h_gjets[0]->SetBinContent(ibin, bsum_gjets);
              h_qcd[0]->SetBinContent(ibin, bsum_qcd);
              h_tg[0]->SetBinContent(ibin, bsum_tg);
              h_wjets[0]->SetBinContent(ibin, bsum_wjets);
              h_ww[0]->SetBinContent(ibin, bsum_ww);
              h_wz[0]->SetBinContent(ibin, bsum_wz);
              h_zz[0]->SetBinContent(ibin, bsum_zz);
              h_wglnu[0]->SetBinContent(ibin, bsum_wglnu);
           }          
   
           Double_t scale = hist[0][0]->Integral()/h_all->Integral();
           h_all->Sumw2(1);
           h_all->Scale(scale);

           TH1F* h_ratio = (TH1F*)hist[0][0]->Clone("h_ratio");
           h_ratio->Sumw2(1);
           h_ratio->Divide(h_all);
           TCanvas *c1 = new TCanvas("c1", "c1");
           c1->cd();

           TPad* pad1 = new TPad("pad1", "pad1", 0.0, 0.3, 1.0, 1.0);
           pad1->SetBottomMargin(0);
           pad1->Draw();
           pad1->cd(); 
            
           hist[0][0]->Draw("ep");
           h_all->Draw("hist, SAME");
           hist[0][0]->Draw("ep, SAME");
           h_all->SetTitle("");
           h_all->GetYaxis()->SetTitle("Events");

           TLegend* leg = new TLegend();
           leg = new TLegend(0.78, 0.65, 0.9, 0.9);
           leg->SetBorderSize(0);
           leg->SetEntrySeparation(0.01);
           leg->SetFillColor(0);
           leg->SetFillStyle(0);

           leg->AddEntry(hist[0][0], "Data", "ep");
           leg->AddEntry(h_dy[0], "DYToEE", "f");
           leg->AddEntry(h_gjets[0], "GJets", "f");
           leg->AddEntry(h_qcd[0], "QCD", "f");
           leg->AddEntry(h_tg[0], "TTG, TG", "f");
           leg->AddEntry(h_wglnu[0], "WGToLNu", "f");
           leg->AddEntry(h_ww[0], "WW", "f");
           leg->AddEntry(h_wz[0], "WZ", "f");
           leg->AddEntry(h_zz[0], "ZZ", "f");
           leg->AddEntry(h_wjets[0], "WJets", "f");
           c1->Update();
           leg->Draw("SAME");

           pad1->SetLogy();

           Int_t year = 17;
           if (is2018) {year=18;}
           CMS_lumi(pad1, year, 0);

           pad1->Update();
           c1->Update();
           c1->cd(); 
           
           TPad* pad2 = new TPad("pad2", "pad2", 0.0, 0.0, 1.0, 0.3);
           pad2->SetTopMargin(0.01);
           pad2->SetBottomMargin(0.2);
           pad2->Draw();
           pad2->cd();

           h_ratio->SetTitle("");
           h_ratio->GetXaxis()->SetTitle(title_axis[0]);
           h_ratio->GetXaxis()->SetTitleFont(42);
           h_ratio->GetXaxis()->SetTitleSize(0.11);
           h_ratio->GetXaxis()->SetTitleOffset(0.8);
           h_ratio->GetXaxis()->SetLabelFont(42);
           h_ratio->GetXaxis()->SetLabelSize(0.1);
           h_ratio->GetYaxis()->SetTitle("Data/MC");
           h_ratio->GetYaxis()->SetTitleSize(0.11);
           h_ratio->GetYaxis()->SetTitleOffset(0.43);
           h_ratio->GetYaxis()->SetLabelSize(0.1);
           h_ratio->GetYaxis()->SetLabelOffset(0.01);
           h_ratio->GetYaxis()->SetNdivisions(505);
           h_ratio->GetYaxis()->SetRangeUser(0.5, 1.5);
           h_ratio->SetMarkerStyle(20);

           TLine* line = new TLine(h_ratio->GetXaxis()->GetXmin(), 1.0, h_ratio->GetXaxis()->GetXmax(), 1.0);
           line->SetLineColor(kRed);
           line->SetLineWidth(1);
           h_ratio->Draw("ep");
           line->Draw("SAME");
           c1->Update();
           c1->SaveAs(savePath+"BDTcumulative_stack.png");
   }




}
