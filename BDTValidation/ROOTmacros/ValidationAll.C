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
using std::setprecision;

double correctBDT(double x){

   double bdt0 = 0.94; ///place where shift is 0                                                                             
   double shift = ( pow(x,2) - pow(bdt0,2) )/( pow(x,2) + pow(bdt0,2) );

   return shift;
}

void ValidationAll(){

   Bool_t is2017 = 0;
   Bool_t is2018 = 1;
   Bool_t isEE = 0;

   Bool_t BDT_shift = 0;  

   Bool_t do_cumulative = 1;
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

   if (is2017) {

      nfile_eff=31;

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
   Float_t elePt;
   Float_t eleEta;
   Double_t zmass;
   Bool_t pass95;
   Int_t doubleCounting;
   Char_t ph_QualityBits;

   TTree *tree[nfile];
   for (Int_t i=0; i<nfile_eff;i++) {
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
       hist[f][18] = new TH1F("ph_R9Full5x5"+sample[f],"ph_R9Full5x5"+sample[f],20,0.4,1.2);
       hist[f][19] = new TH1F("ph_S4Full5x5"+sample[f],"ph_S4Full5x5"+sample[f],20,0.4,1.2);
       hist[f][20] = new TH1F("ph_SigmaIEtaIEta_o_ph_SigmaIPhiIPhi"+sample[f],"ph_SigmaIEtaIEta_o_ph_SigmaIPhiIPhi"+sample[f],35,0,1.4);
       hist[f][21] = new TH1F("zmass"+sample[f],"zmass"+sample[f],30,80,110);
       hist[f][22] = new TH1F("ele_pt"+sample[f],"ele_pt"+sample[f],35,0,350);
       hist[f][23] = new TH1F("ele_eta"+sample[f],"ele_eta"+sample[f],50,-2.5,2.5);
   }

   if (isEE) {
      for (Int_t f=0; f<nfile_eff;f++) {
        hist[f][0] = new TH1F("phoBDT"+sample[f],"phoBDT"+sample[f],25,0,1);
        hist[f][1] = new TH1F("phoEt"+sample[f],"phoEt"+sample[f],20,100,900);
        hist[f][2] = new TH1F("phoEta"+sample[f],"phoEta"+sample[f],50,-2.5,2.5);
        hist[f][3] = new TH1F("ph_PFECALClusIsoCorr"+sample[f],"ph_PFECALClusIsoCorr"+sample[f],30,-100,20);
        hist[f][4] = new TH1F("ph_PFHCALClusIsoCorr"+sample[f],"ph_PFHCALClusIsoCorr"+sample[f],20,-20,140);
        hist[f][5] = new TH1F("ph_TkrIsoCorr"+sample[f],"ph_TkrIsoCorr"+sample[f],20,-5,115);
        hist[f][6] = new TH1F("ph_hoe"+sample[f],"ph_hoe"+sample[f],20,0,0.05);
        hist[f][7] = new TH1F("ph_sieie"+sample[f],"ph_sieie"+sample[f],20,0.02,0.04);
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
        hist[f][18] = new TH1F("ph_R9Full5x5"+sample[f],"ph_R9Full5x5"+sample[f],20,0.4,1.2);
        hist[f][19] = new TH1F("ph_S4Full5x5"+sample[f],"ph_S4Full5x5"+sample[f],20,0.4,1.2);
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

        if ((ph_QualityBits == 0) || (ph_QualityBits == 2)) continue;

        if (is2017) {
           if ((fil > 11 && fil <= 18) && (doubleCounting == 1)) continue; //remove double counting QCD
           //if (fil > 25 && doubleCounting == 1) continue; //remove double counting W+jets
        }
        if (is2018) {
           if ((fil > 11 && fil <= 18) && (doubleCounting == 1)) continue; //remove double counting QCD
        }

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

          if (!BDT_shift) {
             hist[fil][0]->Fill(ph_BDTpred, puWeight/sumGenWeight);
          }
 
          if ((BDT_shift) && (fil<=8)) { //shifting only DY

             Double_t deltaShift = 0.05;
             Double_t shift = deltaShift * correctBDT(ph_BDTpred);
             Double_t newBDTscore = ph_BDTpred+shift;

             newBDTscore = max(0.,newBDTscore);
             newBDTscore = min(1.,newBDTscore);

             hist[fil][0]->Fill(newBDTscore, puWeight/sumGenWeight);

          } else if ((BDT_shift) && (fil>8)) { 
            hist[fil][0]->Fill(ph_BDTpred, puWeight/sumGenWeight);
          }

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
   TString title_axis[nhist] = {"BDT score", "p_{T}^{#gamma}", "#eta^{#gamma}", "PF ECALClusIsoCorr", "PF HCALClusIsoCorr", "TkrIsoCorr", "H/E", "#sigma_{i#eta i #eta}", "#sigma_{i#eta i#phi}","#sigma_{i#phi i #phi}", "E_{2x2}/E_{3x3}","E_{1x3}/E^{raw}_{SC}","E_{2}/E^{raw}_{SC}","E_{2x5}/E^{raw}_{SC}","E_{max}/E^{raw}_{SC}","#sigma_{#eta}", "#sigma_{#phi}", "#sigma_{#eta}/#sigma_{#phi}", "S4","R9", "#sigma_{i#eta i #eta}/#sigma_{i#phi i #phi}", "e#gamma invariant mass", "p_{T}^{e}", "#eta^{e}" };

   TH1F* h_ratio[nhist];
   TH1F* h_all[nhist];

   Double_t overflow[nfile][nhist];
   Double_t underflow[nfile][nhist];
   for (Int_t fil=0; fil<nfile_eff; fil++) {
      for (int i=0; i<nhist; i++) {
         overflow[fil][i] = hist[fil][i]->GetBinContent(hist[fil][i]->GetNbinsX()+1);
         underflow[fil][i] = hist[fil][i]->GetBinContent(0);
         hist[fil][i]->SetBinContent(1, hist[fil][i]->GetBinContent(1)+underflow[fil][i]);
         hist[fil][i]->SetBinContent(hist[fil][i]->GetNbinsX(), hist[fil][i]->GetBinContent(hist[fil][i]->GetNbinsX())+overflow[fil][i]);
         //if (fil==4) {std::cout << title[i] << ":  "<< setprecision(7) << hist[4][i]->Integral() << std::endl;}
      }
   }

   for (int i=0; i<nhist; i++) {

      h_all[i] = (TH1F*)hist[1][i]->Clone("h_all");
      for (int j=2; j<nfile_eff; j++) { //sum signal+bkgs correctly normalized
        h_all[i]->Add(hist[j][i]);
      }

      Double_t scale = hist[0][i]->Integral()/h_all[i]->Integral();
     
      h_all[i]->Sumw2(1);
      h_all[i]->Scale(scale);
      
      h_ratio[i] = (TH1F*)hist[0][i]->Clone("h_ratio");
      h_ratio[i]->Sumw2(1);
      h_ratio[i]->Divide(h_all[i]);
      TCanvas *c1 = new TCanvas("c1", "c1");
      c1->cd();

      TPad* pad1 = new TPad("pad1", "pad1", 0.0, 0.3, 1.0, 1.0);
      pad1->SetBottomMargin(0);
      pad1->Draw();
      pad1->cd();

      hs[i] = new THStack(title[i]+"_all", "");

      hist[0][i]->SetMarkerStyle(20); //data
      hist[0][i]->SetTitle("");

      h_all[i]->SetFillColor(kOrange);
      h_all[i]->SetLineColor(kOrange);

      hs[i]->Add(h_all[i]);
 
      hist[0][i]->Draw("ep");
      hs[i]->Draw("hist, SAME");
      hist[0][i]->Draw("ep, SAME");
      hs[i]->SetTitle("");
      hs[i]->GetYaxis()->SetTitle("Events");
     
      TLegend* leg = new TLegend();
      leg = new TLegend(0.78, 0.65, 0.9, 0.9);
      leg->SetBorderSize(0);
      leg->SetEntrySeparation(0.01);
      leg->SetFillColor(0);
      leg->SetFillStyle(0);

      leg->AddEntry(hist[0][i], "Data", "ep");
      leg->AddEntry(h_all[i], "all MC", "f");
      

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
      c1->SaveAs(savePath+title[i]+"_all.png");
   }

   if (do_cumulative) {

      TH1F* h_dy = (TH1F*)hist[1][0]->Clone("h_dy");

      /*for (int j=2; j<=8; j++) { //only signal
        h_dy->Add(hist[j][0]);
      }*/

      for (int j=2; j<nfile_eff; j++) { //signal and bkg
        h_dy->Add(hist[j][0]);
      }

      Double_t bsum_data = 0;
      Double_t bsum_mc = 0;
      std::cout << "bdt: "<< h_dy->GetBinContent(1) << " " << hist[0][0]->GetBinContent(1) << std::endl;
      for (int ibin=1; ibin<=hist[0][0]->GetNbinsX(); ibin++) {
         std::cout << bsum_mc << std::endl;
         std::cout << bsum_data << std::endl;
         bsum_data += hist[0][0]->GetBinContent(ibin);
         bsum_mc += h_dy->GetBinContent(ibin);
         hist[0][0]->SetBinContent(ibin, bsum_data);
         h_dy->SetBinContent(ibin, bsum_mc);
      }


      Double_t scale = hist[0][0]->Integral()/h_dy->Integral();
      h_dy->Sumw2(1);
      h_dy->Scale(scale);

      TH1F* h_ratio = (TH1F*)hist[0][0]->Clone("h_ratio");
      h_ratio->Sumw2(1);
      h_ratio->Divide(h_dy);
      TCanvas *c1 = new TCanvas("c1", "c1");
      c1->cd();

      TPad* pad1 = new TPad("pad1", "pad1", 0.0, 0.3, 1.0, 1.0);
      pad1->SetBottomMargin(0);
      pad1->Draw();
      pad1->cd();     

      hist[0][0]->SetMarkerStyle(20); //data
      hist[0][0]->SetTitle("");

      h_dy->SetFillColor(kOrange);
      h_dy->SetLineColor(kOrange);

      hist[0][0]->Draw("ep");
      h_dy->Draw("hist, SAME");
      hist[0][0]->Draw("ep, SAME");
      h_dy->SetTitle("");
      h_dy->GetYaxis()->SetTitle("Events");

      TLegend* leg = new TLegend();
      leg = new TLegend(0.78, 0.65, 0.9, 0.9);
      leg->SetBorderSize(0);
      leg->SetEntrySeparation(0.01);
      leg->SetFillColor(0);
      leg->SetFillStyle(0);

      leg->AddEntry(hist[0][0], "Data", "ep");
      leg->AddEntry(h_dy, "allMC", "f");

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
      c1->SaveAs(savePath+"BDTcumulative.png");
   }








}
