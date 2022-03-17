#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <iostream>
#include <fstream>
#include <THStack.h>
#include "CMS_lumi.C"
#include <TStyle.h>
#include "rebin.h"

using namespace std;

Float_t roundoff(float value) {
  Float_t pow_10 = pow(10.0f, 2);
  return round(value * pow_10) / pow_10;
}

//Float_t getSF(TH2F* h, Float_t eta, Float_t pt) { 
Float_t getSF(Float_t eta, Float_t pt) { 
  const int mpt = 9999;
  Float_t pt_min[mpt];
  Float_t pt_max[mpt];
  Float_t eta_min[mpt];
  Float_t eta_max[mpt];
  Float_t sf[mpt];
  Float_t sf_ = -9999;

  /*TAxis *xaxis = h->GetXaxis();
  TAxis *yaxis = h->GetYaxis();
  Int_t binx = xaxis->FindBin(pt);
  Int_t biny = yaxis->FindBin(eta);
  sf_ = h->GetBinContent(binx,biny);
  */
  ifstream in;
  in.open("SFs_UL17.txt");
  Int_t npt = 0;
  while (1) {
     in >> eta_min[npt] >> eta_max[npt] >> pt_min[npt] >> pt_max[npt] >> sf[npt];

     if ((eta > eta_min[npt]) && (eta < eta_max[npt]) && (pt > pt_min[npt]) && (pt < pt_max[npt])) {
       sf_ = sf[npt];
     }

     if (sf_ > 0 ) {
        break;
     }

     npt++;
     if (( ! in.good() )) {
       std::cout << "final break" << std::endl;
       break;
     }
  }
  if (sf_ < 0) std::cout << eta << " " << pt << "---> sf: " << sf_ << std::endl;

  return sf_;
}

Float_t eleFakeWeight(Float_t eta) {

   Float_t weight = 0.0504885 + -0.00999998*fabs(eta) + -0.0112097*eta*eta;
   return weight;

}

Float_t jetFakeWeight(Float_t pt) {

   Float_t weight = (1.5e-5*pt + (-6.0e-3) + 2.3)/(pt-198.2) ;
   return weight;
}
void plots(){

   Bool_t is2016_preVFP = 0;
   Bool_t is2016_postVFP = 0;

   Bool_t is2016 = 0;
   Bool_t is2017 = 1;
   Bool_t is2018 = 0;

   Bool_t isEE = 0;

   const Int_t nfile=99;
   Int_t nfile_eff = 99;
   TFile *file[nfile];

   TString path = "/local/cms/user/gsorrent/antgc_analysis/sigRegion/EBNtuples_datamc/merged/";
   TString path2 = "/local/cms/user/gsorrent/antgc_analysis/sigRegion/EBNtuples_datamc/merged/";
   TString savePath = "/local/cms/user/gsorrent/antgc_analysis/sigRegion/EBPlots/";

   if (isEE) {
      path = "/local/cms/user/gsorrent/antgc_analysis/BDTValidation/UL17/EENtuples_good/merged/";
      savePath = "/local/cms/user/gsorrent/antgc_analysis/BDTValidation/UL17/EEPlots/";
   }

   if (is2018) {
      path = "/local/cms/user/gsorrent/antgc_analysis/BDTValidation/UL18/EBNtuples_good/merged/";
      savePath = "/local/cms/user/gsorrent/antgc_analysis/BDTValidation/UL18/EBPlots/";
      if (isEE) {
        path = "/local/cms/user/gsorrent/antgc_analysis/BDTValidation/UL18/EENtuples_good/merged/";
        savePath = "/local/cms/user/gsorrent/antgc_analysis/BDTValidation/UL18/EEPlots/";
      }
   }

   if (is2016_preVFP) {
      path = "/local/cms/user/gsorrent/antgc_analysis/BDTValidation/UL16preVFP/EBNtuples/merged/";
      savePath = "/local/cms/user/gsorrent/antgc_analysis/BDTValidation/UL16preVFP/EBPlots/";
      if (isEE) {
        path = "/local/cms/user/gsorrent/antgc_analysis/BDTValidation/UL16preVFP/EENtuples/merged/";
        savePath = "/local/cms/user/gsorrent/antgc_analysis/BDTValidation/UL16preVFP/EEPlots/";
      }
   }

   if (is2016_postVFP) {
      path = "/local/cms/user/gsorrent/antgc_analysis/BDTValidation/UL16postVFP/EBNtuples/merged/";
      savePath = "/local/cms/user/gsorrent/antgc_analysis/BDTValidation/UL16postVFP/EBPlots/";
      if (isEE) {
        path = "/local/cms/user/gsorrent/antgc_analysis/BDTValidation/UL16postVFP/EENtuples/merged/";
        savePath = "/local/cms/user/gsorrent/antgc_analysis/BDTValidation/UL16postVFP/EEPlots/";
      }
   }

   if (is2016) {
      path = "/local/cms/user/gsorrent/antgc_analysis/BDTValidation/UL16preVFP/EBNtuples/merged/";
      path2 = "/local/cms/user/gsorrent/antgc_analysis/BDTValidation/UL16postVFP/EBNtuples/merged/";
      savePath = "/local/cms/user/gsorrent/antgc_analysis/BDTValidation/UL16/EBPlots/";
      if (isEE) {
        path = "/local/cms/user/gsorrent/antgc_analysis/BDTValidation/UL16preVFP/EENtuples/merged/";
        path2 = "/local/cms/user/gsorrent/antgc_analysis/BDTValidation/UL16postVFP/EENtuples/merged/";
        savePath = "/local/cms/user/gsorrent/antgc_analysis/BDTValidation/UL16/EEPlots/";
      }
   }

   ofstream fraction_file;
   fraction_file.open(savePath+"bkg_fractions.txt");

   if (is2017) {

      nfile_eff=24;
      file[0] = new TFile(path+"data.root");
      file[1] = new TFile(path+"ZNuNuGJetsMonoPhotonPtG130TuneCP513TeVamcatnloFXFXpythia8RunIISummer20UL17MiniAOD106X.root");
      file[2] = new TFile(path+"WGJetsMonoPhotonPtG130TuneCP513TeVmadgraphpythia8RunIISummer20UL17MiniAOD106X.root");
      file[3] = new TFile(path+"ZLLGJetsMonoPhotonPtG130TuneCP513TeVamcatnloFXFXpythia8RunIISummer20UL17MiniAOD106X.root");
      file[4] = new TFile(path+"DiPhotonJetsMGG80toInfTuneCP513TeVamcatnloFXFXpythia8RunIISummer20UL17MiniAOD106X.root");
      file[5] = new TFile(path+"GJetsHT200To400TuneCP513TeVmadgraphMLMpythia8RunIISummer19UL17MiniAODv2106X.root");
      file[6] = new TFile(path+"GJetsHT400To600TuneCP513TeVmadgraphMLMpythia8RunIISummer19UL17MiniAODv2106X.root");
      file[7] = new TFile(path+"GJetsHT600ToInfTuneCP513TeVmadgraphMLMpythia8RunIISummer19UL17MiniAODv2106X.root");
      file[8] = new TFile(path+"QCDHT1000to1500TuneCP5PSWeights13TeVmadgraphMLMpythia8RunIISummer19UL17MiniAODv2106X.root");
      file[9] = new TFile(path+"QCDHT1500to2000TuneCP5PSWeights13TeVmadgraphMLMpythia8RunIISummer19UL17MiniAODv2106X.root");
      file[10] = new TFile(path+"QCDHT2000toInfTuneCP5PSWeights13TeVmadgraphMLMpythia8RunIISummer19UL17MiniAODv2106X.root");
      file[11] = new TFile(path+"QCDHT200to300TuneCP5PSWeights13TeVmadgraphMLMpythia8RunIISummer19UL17MiniAODv2106X.root");
      file[12] = new TFile(path+"QCDHT300to500TuneCP5PSWeights13TeVmadgraphMLMpythia8RunIISummer19UL17MiniAODv2106X.root");
      file[13] = new TFile(path+"QCDHT500to700TuneCP5PSWeights13TeVmadgraphMLMpythia8RunIISummer19UL17MiniAODv2106X.root");
      file[14] = new TFile(path+"QCDHT700to1000TuneCP5PSWeights13TeVmadgraphMLMpythia8RunIISummer19UL17MiniAODv2106X.root");
      file[15] = new TFile(path+"TGJetsTuneCP513TeVamcatnlomadspinpythia8RunIISummer19UL17MiniAODv2106X.root");
      file[16] = new TFile(path+"TTGJetsTuneCP513TeVamcatnloFXFXmadspinpythia8RunIISummer20UL17MiniAODv2106X.root");
      file[17] = new TFile(path+"WWTuneCP513TeVpythia8RunIISummer20UL17MiniAOD106X.root");
      file[18] = new TFile(path+"WZTuneCP513TeVpythia8RunIISummer20UL17MiniAODv2106X.root");
      file[19] = new TFile(path+"ZZTuneCP513TeVpythia8RunIISummer20UL17MiniAODv2106X.root");
      file[20] = new TFile(path+"WToMuNuM100TuneCP513TeVpythia8RunIISummer20UL17MiniAOD106X.root");
      file[21] = new TFile(path+"WToTauNuM100TuneCP513TeVpythia8tauolaRunIISummer20UL17MiniAOD106X.root");
      file[22] = new TFile(path+"efakes.root");
      file[23] = new TFile(path+"jetfakes.root");
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
     file[30] = new TFile(path+"WJetsToLNuHT800To1200TuneCP513TeVmadgraphMLMpythia8RunIISummer20UL18MiniAODv2106X.root");
     file[31] = new TFile(path+"DiPhotonJetsMGG80toInfTuneCP513TeVamcatnloFXFXpythia8RunIISummer20UL18MiniAOD106X.root");
   }

   if (is2016_preVFP) {

     nfile_eff=24;
     file[0] = new TFile(path+"data.root");
     file[1] = new TFile(path+"DYJetsToEEM50massWgtFixTuneCP513TeVpowhegMiNNLOpythia8photosRunIISummer20UL16MiniAODAPV106XpreVFP.root");
     //file[1] = new TFile("/local/cms/user/gsorrent/antgc_analysis/BDTValidation/UL16preVFP/EBNtuples_good/merged/DYJetsToEEM50massWgtFixTuneCP513TeVpowhegMiNNLOpythia8photosRunIISummer20UL16MiniAODAPV106XpreVFP.root");
     file[2] = new TFile(path+"GJetsHT200To400TuneCP513TeVmadgraphMLMpythia8RunIISummer20UL16MiniAODAPV106XpreVFP.root");
     file[3] = new TFile(path+"GJetsHT400To600TuneCP513TeVmadgraphMLMpythia8RunIISummer20UL16MiniAODAPV106XpreVFP.root");
     file[4] = new TFile(path+"GJetsHT600ToInfTuneCP513TeVmadgraphMLMpythia8RunIISummer20UL16MiniAODAPV106XpreVFP.root");
     file[5] = new TFile(path+"QCDHT1000to1500TuneCP5PSWeights13TeVmadgraphMLMpythia8RunIISummer20UL16MiniAODAPV106XpreVFP.root");
     file[6] = new TFile(path+"QCDHT1500to2000TuneCP5PSWeights13TeVmadgraphMLMpythia8RunIISummer20UL16MiniAODAPV106XpreVFP.root");
     file[7] = new TFile(path+"QCDHT2000toInfTuneCP5PSWeights13TeVmadgraphMLMpythia8RunIISummer20UL16MiniAODAPV106XpreVFP.root");
     file[8] = new TFile(path+"QCDHT200to300TuneCP5PSWeights13TeVmadgraphMLMpythia8RunIISummer20UL16MiniAODAPV106XpreVFP.root");
     file[9] = new TFile(path+"QCDHT300to500TuneCP5PSWeights13TeVmadgraphMLMpythia8RunIISummer20UL16MiniAODAPV106XpreVFP.root");
     file[10] = new TFile(path+"QCDHT500to700TuneCP5PSWeights13TeVmadgraphMLMpythia8RunIISummer20UL16MiniAODAPV106XpreVFP.root");
     file[11] = new TFile(path+"QCDHT700to1000TuneCP5PSWeights13TeVmadgraphMLMpythia8RunIISummer20UL16MiniAODAPV106XpreVFP.root");
     file[12] = new TFile(path+"WGToLNuG01J5fPtG130TuneCP513TeVamcatnloFXFXpythia8RunIISummer20UL16MiniAODAPV106XpreVFP.root");
     file[13] = new TFile(path+"TTGJetsTuneCP513TeVamcatnloFXFXmadspinpythia8RunIISummer20UL16MiniAODAPV106XpreVFP.root");
     file[14] = new TFile(path+"WWTuneCP513TeVpythia8RunIISummer20UL16MiniAODAPV106XpreVFP.root");
     file[15] = new TFile(path+"WZTuneCP513TeVpythia8RunIISummer20UL16MiniAODAPV106XpreVFP.root");
     file[16] = new TFile(path+"ZZTuneCP513TeVpythia8RunIISummer20UL16MiniAODAPV106XpreVFP.root");
     file[17] = new TFile(path+"WJetsToLNuHT100To200TuneCP513TeVmadgraphMLMpythia8RunIISummer20UL16MiniAODAPV106XpreVFP.root");
     file[18] = new TFile(path+"WJetsToLNuHT200To400TuneCP513TeVmadgraphMLMpythia8RunIISummer20UL16MiniAODAPV106XpreVFP.root");
     file[19] = new TFile(path+"WJetsToLNuHT400To600TuneCP513TeVmadgraphMLMpythia8RunIISummer20UL16MiniAODAPV106XpreVFP.root");
     file[20] = new TFile(path+"WJetsToLNuHT600To800TuneCP513TeVmadgraphMLMpythia8RunIISummer20UL16MiniAODAPV106XpreVFP.root");
     file[21] = new TFile(path+"WJetsToLNuHT70To100TuneCP513TeVmadgraphMLMpythia8RunIISummer20UL16MiniAODAPV106XpreVFP.root");
     file[22] = new TFile(path+"WJetsToLNuHT800To1200TuneCP513TeVmadgraphMLMpythia8RunIISummer20UL16MiniAODAPV106XpreVFP.root");
     file[23] = new TFile(path+"DiPhotonJetsMGG80toInfTuneCP513TeVamcatnloFXFXpythia8RunIISummer20UL16MiniAODAPV106XpreVFP.root");
   }


   if (is2016) {

     nfile_eff=48;
     file[0] = new TFile(path+"alldata.root");
     file[1] = new TFile(path+"DYJetsToEEM50massWgtFixTuneCP513TeVpowhegMiNNLOpythia8photosRunIISummer20UL16MiniAODAPV106XpreVFP.root");
     file[2] = new TFile(path2+"DYJetsToEEM50massWgtFixTuneCP513TeVpowhegMiNNLOpythia8photosRunIISummer20UL16MiniAOD106X.root");
     file[3] = new TFile(path+"GJetsHT200To400TuneCP513TeVmadgraphMLMpythia8RunIISummer20UL16MiniAODAPV106XpreVFP.root");
     file[4] = new TFile(path+"GJetsHT400To600TuneCP513TeVmadgraphMLMpythia8RunIISummer20UL16MiniAODAPV106XpreVFP.root");
     file[5] = new TFile(path+"GJetsHT600ToInfTuneCP513TeVmadgraphMLMpythia8RunIISummer20UL16MiniAODAPV106XpreVFP.root");
     file[6] = new TFile(path2+"GJetsHT200To400TuneCP513TeVmadgraphMLMpythia8RunIISummer20UL16MiniAOD106X.root");
     file[7] = new TFile(path2+"GJetsHT400To600TuneCP513TeVmadgraphMLMpythia8RunIISummer20UL16MiniAOD106X.root");
     file[8] = new TFile(path2+"GJetsHT600ToInfTuneCP513TeVmadgraphMLMpythia8RunIISummer20UL16MiniAOD106X.root");
     file[9] = new TFile(path+"QCDHT1000to1500TuneCP5PSWeights13TeVmadgraphMLMpythia8RunIISummer20UL16MiniAODAPV106XpreVFP.root");
     file[10] = new TFile(path+"QCDHT1500to2000TuneCP5PSWeights13TeVmadgraphMLMpythia8RunIISummer20UL16MiniAODAPV106XpreVFP.root");
     file[11] = new TFile(path+"QCDHT2000toInfTuneCP5PSWeights13TeVmadgraphMLMpythia8RunIISummer20UL16MiniAODAPV106XpreVFP.root");
     file[12] = new TFile(path+"QCDHT200to300TuneCP5PSWeights13TeVmadgraphMLMpythia8RunIISummer20UL16MiniAODAPV106XpreVFP.root");
     file[13] = new TFile(path+"QCDHT300to500TuneCP5PSWeights13TeVmadgraphMLMpythia8RunIISummer20UL16MiniAODAPV106XpreVFP.root");
     file[14] = new TFile(path+"QCDHT500to700TuneCP5PSWeights13TeVmadgraphMLMpythia8RunIISummer20UL16MiniAODAPV106XpreVFP.root");
     file[15] = new TFile(path+"QCDHT700to1000TuneCP5PSWeights13TeVmadgraphMLMpythia8RunIISummer20UL16MiniAODAPV106XpreVFP.root");
     file[16] = new TFile(path2+"QCDHT1000to1500TuneCP5PSWeights13TeVmadgraphMLMpythia8RunIISummer20UL16MiniAODv2106X.root");
     file[17] = new TFile(path2+"QCDHT1500to2000TuneCP5PSWeights13TeVmadgraphMLMpythia8RunIISummer20UL16MiniAODv2106X.root");
     file[18] = new TFile(path2+"QCDHT200to300TuneCP5PSWeights13TeVmadgraphMLMpythia8RunIISummer20UL16MiniAOD106X.root");
     file[19] = new TFile(path2+"QCDHT300to500TuneCP5PSWeights13TeVmadgraphMLMpythia8RunIISummer20UL16MiniAOD106X.root");
     file[20] = new TFile(path2+"QCDHT500to700TuneCP5PSWeights13TeVmadgraphMLMpythia8RunIISummer20UL16MiniAOD106X.root");
     file[21] = new TFile(path2+"QCDHT700to1000TuneCP5PSWeights13TeVmadgraphMLMpythia8RunIISummer20UL16MiniAOD106X.root");
     file[22] = new TFile(path2+"QCDHT2000toInfTuneCP5PSWeights13TeVmadgraphMLMpythia8RunIISummer20UL16MiniAOD106X.root");
     file[23] = new TFile(path+"WGToLNuG01J5fPtG130TuneCP513TeVamcatnloFXFXpythia8RunIISummer20UL16MiniAODAPV106XpreVFP.root");
     //file[24] = new TFile(path2+"WGJetsMonoPhotonPtG130TuneCP513TeVmadgraphpythia8RunIISummer20UL16MiniAOD106X.root"); 
     file[24] = new TFile(path2+"WGToLNuG01J5fPtG130TuneCP513TeVamcatnloFXFXpythia8RunIISummer20UL16MiniAOD106X.root");
     file[25] = new TFile(path+"TTGJetsTuneCP513TeVamcatnloFXFXmadspinpythia8RunIISummer20UL16MiniAODAPV106XpreVFP.root");
     file[26] = new TFile(path2+"TTGJetsTuneCP513TeVamcatnloFXFXmadspinpythia8RunIISummer20UL16MiniAOD106X.root");
     file[27] = new TFile(path+"WWTuneCP513TeVpythia8RunIISummer20UL16MiniAODAPV106XpreVFP.root");
     file[28] = new TFile(path2+"WWTuneCP513TeVpythia8RunIISummer20UL16MiniAOD106X.root");
     file[29] = new TFile(path+"WZTuneCP513TeVpythia8RunIISummer20UL16MiniAODAPV106XpreVFP.root");
     file[30] = new TFile(path2+"WZTuneCP513TeVpythia8RunIISummer20UL16MiniAOD106X.root");
     file[31] = new TFile(path+"ZZTuneCP513TeVpythia8RunIISummer20UL16MiniAODAPV106XpreVFP.root");
     file[32] = new TFile(path2+"ZZTuneCP513TeVpythia8RunIISummer20UL16MiniAOD106X.root");
     file[33] = new TFile(path+"WJetsToLNuHT100To200TuneCP513TeVmadgraphMLMpythia8RunIISummer20UL16MiniAODAPV106XpreVFP.root");
     file[34] = new TFile(path+"WJetsToLNuHT200To400TuneCP513TeVmadgraphMLMpythia8RunIISummer20UL16MiniAODAPV106XpreVFP.root");
     file[35] = new TFile(path+"WJetsToLNuHT400To600TuneCP513TeVmadgraphMLMpythia8RunIISummer20UL16MiniAODAPV106XpreVFP.root");
     file[36] = new TFile(path+"WJetsToLNuHT600To800TuneCP513TeVmadgraphMLMpythia8RunIISummer20UL16MiniAODAPV106XpreVFP.root");
     file[37] = new TFile(path+"WJetsToLNuHT70To100TuneCP513TeVmadgraphMLMpythia8RunIISummer20UL16MiniAODAPV106XpreVFP.root");
     file[38] = new TFile(path+"WJetsToLNuHT800To1200TuneCP513TeVmadgraphMLMpythia8RunIISummer20UL16MiniAODAPV106XpreVFP.root");
     file[39] = new TFile(path2+"WJetsToLNuHT100To200TuneCP513TeVmadgraphMLMpythia8RunIISummer20UL16MiniAOD106X.root");
     file[40] = new TFile(path2+"WJetsToLNuHT1200To2500TuneCP513TeVmadgraphMLMpythia8RunIISummer20UL16MiniAOD106X.root");
     file[41] = new TFile(path2+"WJetsToLNuHT200To400TuneCP513TeVmadgraphMLMpythia8RunIISummer20UL16MiniAOD106X.root");
     file[42] = new TFile(path2+"WJetsToLNuHT400To600TuneCP513TeVmadgraphMLMpythia8RunIISummer20UL16MiniAOD106X.root");
     file[43] = new TFile(path2+"WJetsToLNuHT600To800TuneCP513TeVmadgraphMLMpythia8RunIISummer20UL16MiniAOD106X.root");
     file[44] = new TFile(path2+"WJetsToLNuHT70To100TuneCP513TeVmadgraphMLMpythia8RunIISummer20UL16MiniAOD106X.root");
     file[45] = new TFile(path2+"WJetsToLNuHT800To1200TuneCP513TeVmadgraphMLMpythia8RunIISummer20UL16MiniAOD106X.root");
     file[46] = new TFile(path+"DiPhotonJetsMGG80toInfTuneCP513TeVamcatnloFXFXpythia8RunIISummer20UL16MiniAODAPV106XpreVFP.root");
     file[47] = new TFile(path2+"DiPhotonJetsMGG80toInfTuneCP513TeVamcatnloFXFXpythia8RunIISummer20UL16MiniAOD106X.root");
   }

   if (is2016_postVFP) {
     
     nfile_eff=24; 
     file[0] = new TFile(path+"data.root");
     file[1] = new TFile(path+"DYJetsToEEM50massWgtFixTuneCP513TeVpowhegMiNNLOpythia8photosRunIISummer20UL16MiniAOD106X.root");
     file[2] = new TFile(path+"GJetsHT200To400TuneCP513TeVmadgraphMLMpythia8RunIISummer20UL16MiniAOD106X.root");
     file[3] = new TFile(path+"GJetsHT400To600TuneCP513TeVmadgraphMLMpythia8RunIISummer20UL16MiniAOD106X.root");
     file[4] = new TFile(path+"GJetsHT600ToInfTuneCP513TeVmadgraphMLMpythia8RunIISummer20UL16MiniAOD106X.root");
     file[5] = new TFile(path+"QCDHT1000to1500TuneCP5PSWeights13TeVmadgraphMLMpythia8RunIISummer20UL16MiniAODv2106X.root");
     file[6] = new TFile(path+"QCDHT1500to2000TuneCP5PSWeights13TeVmadgraphMLMpythia8RunIISummer20UL16MiniAODv2106X.root");
     file[7] = new TFile(path+"QCDHT200to300TuneCP5PSWeights13TeVmadgraphMLMpythia8RunIISummer20UL16MiniAOD106X.root");
     file[8] = new TFile(path+"QCDHT300to500TuneCP5PSWeights13TeVmadgraphMLMpythia8RunIISummer20UL16MiniAOD106X.root");
     file[9] = new TFile(path+"QCDHT500to700TuneCP5PSWeights13TeVmadgraphMLMpythia8RunIISummer20UL16MiniAOD106X.root");
     file[10] = new TFile(path+"QCDHT700to1000TuneCP5PSWeights13TeVmadgraphMLMpythia8RunIISummer20UL16MiniAOD106X.root");
     file[11] = new TFile(path+"QCDHT2000toInfTuneCP5PSWeights13TeVmadgraphMLMpythia8RunIISummer20UL16MiniAOD106X.root");
     file[12] = new TFile(path+"WGToLNuG01J5fPtG130TuneCP513TeVamcatnloFXFXpythia8RunIISummer20UL16MiniAOD106X.root");
     file[13] = new TFile(path+"TTGJetsTuneCP513TeVamcatnloFXFXmadspinpythia8RunIISummer20UL16MiniAOD106X.root");
     file[14] = new TFile(path+"WWTuneCP513TeVpythia8RunIISummer20UL16MiniAOD106X.root");
     file[15] = new TFile(path+"WZTuneCP513TeVpythia8RunIISummer20UL16MiniAOD106X.root");
     file[16] = new TFile(path+"ZZTuneCP513TeVpythia8RunIISummer20UL16MiniAOD106X.root");
     file[17] = new TFile(path+"WJetsToLNuHT100To200TuneCP513TeVmadgraphMLMpythia8RunIISummer20UL16MiniAOD106X.root");
     file[18] = new TFile(path+"WJetsToLNuHT200To400TuneCP513TeVmadgraphMLMpythia8RunIISummer20UL16MiniAOD106X.root");
     file[19] = new TFile(path+"WJetsToLNuHT400To600TuneCP513TeVmadgraphMLMpythia8RunIISummer20UL16MiniAOD106X.root");
     file[20] = new TFile(path+"WJetsToLNuHT600To800TuneCP513TeVmadgraphMLMpythia8RunIISummer20UL16MiniAOD106X.root");
     file[21] = new TFile(path+"WJetsToLNuHT70To100TuneCP513TeVmadgraphMLMpythia8RunIISummer20UL16MiniAOD106X.root");
     file[22] = new TFile(path+"WJetsToLNuHT800To1200TuneCP513TeVmadgraphMLMpythia8RunIISummer20UL16MiniAOD106X.root");
     file[23] = new TFile(path+"DiPhotonJetsMGG80toInfTuneCP513TeVamcatnloFXFXpythia8RunIISummer20UL16MiniAOD106X.root");
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

   Float_t  ph_SeedTime;
   Float_t  deltaPhiPUPPImetPho; 
   Float_t  phoPtOverPUPPImet;
   Float_t  puppiMET;
   Float_t  minDeltaPhiPUPPImetJet30;
   Float_t  nPUPPIJet30;
   Float_t  PUPPIJetHt30;

   Bool_t   hasElectron;
   Bool_t   hasMuon;
   Bool_t   eleVeto;
   Bool_t   muoVeto;
   Bool_t   tauVeto;
   Bool_t   haspixelseed;
   Bool_t   passjetfakes;

   Int_t    genPDGid;

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

     tree[i]->SetBranchAddress("ph_seedtime", &ph_SeedTime);
     tree[i]->SetBranchAddress("deltaPhiPUPPImetPho", &deltaPhiPUPPImetPho);
     tree[i]->SetBranchAddress("phoPtOverPUPPImet", &phoPtOverPUPPImet);
     tree[i]->SetBranchAddress("puppiMET", &puppiMET);
     tree[i]->SetBranchAddress("minDeltaPhiPUPPImetJet30", &minDeltaPhiPUPPImetJet30);

     tree[i]->SetBranchAddress("hasElectron", &hasElectron);
     tree[i]->SetBranchAddress("hasMuon", &hasMuon);
     tree[i]->SetBranchAddress("eleVeto", &eleVeto);
     tree[i]->SetBranchAddress("muoVeto", &muoVeto);
     tree[i]->SetBranchAddress("tauVeto", &tauVeto);
     tree[i]->SetBranchAddress("haspixelseed", &haspixelseed);
     tree[i]->SetBranchAddress("passjetfakes", &passjetfakes);
     tree[i]->SetBranchAddress("genPDGid", &genPDGid);

     tree[i]->SetBranchAddress("nPUPPIJet30", &nPUPPIJet30);
     tree[i]->SetBranchAddress("PUPPIJetHt30", &PUPPIJetHt30);
   }

   const Int_t nhist=4;
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
       hist[f][0] = new TH1F("phoEta"+sample[f],"phoEta"+sample[f],15,0,1.5);
       hist[f][1] = new TH1F("phoPhi"+sample[f],"phoPhi"+sample[f],10,-3.15,3.15);
       hist[f][2] = new TH1F("nPUPPIjet"+sample[f],"nPUPPIjet"+sample[f],6,0,6);
       hist[f][3] = new TH1F("PUPPIjetHT"+sample[f],"PUPPIjetHT"+sample[f],500,0,1000);
       rebin(hist[f][3]);
       /*hist[f][0] = new TH1F("phoBDT"+sample[f],"phoBDT"+sample[f],25,0,1);
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
       hist[f][23] = new TH1F("ele_eta"+sample[f],"ele_eta"+sample[f],30,-2.5,2.5);*/
   }

   if (isEE) {
      for (Int_t f=0; f<nfile_eff;f++) {
        /*hist[f][0] = new TH1F("phoBDT"+sample[f],"phoBDT"+sample[f],25,0,1);
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
        hist[f][23] = new TH1F("ele_eta"+sample[f],"ele_eta"+sample[f],30,-2.5,2.5);*/
     }
   }

   //Sumw
   for (Int_t j=0; j<nfile_eff; j++) {
      for (Int_t i=0; i<nhist; i++) {
         hist[j][i]->Sumw2(1);
      }
   }

   const int mpt = 9999;
   Float_t pt_min[mpt];
   Float_t pt_max[mpt];
   Float_t eta_min[mpt];
   Float_t eta_max[mpt];
   Float_t sf[mpt];

   Int_t nentries[nfile] = {0};
   for (Int_t i=0; i<nfile_eff;i++) {
     nentries[i] = (Int_t)tree[i]->GetEntries();
   }

   //TFile *filein = TFile::Open("/local/cms/user/gsorrent/antgc_analysis/ScaleFactorNtuples/root/UL17_EB.root");
   //TH2F* sfHist = (TH2F*)filein->Get("EGamma_SF2D");

   for (Int_t fil=0; fil<nfile_eff; fil++) {
    TH1F* cutFlowGenWeight = (TH1F*)file[fil]->Get("tnpPhoIDs/tnpPhoIDscutFlowGenWeight");
    Float_t sumGenWeight = cutFlowGenWeight->GetBinContent(1);
    //std::cout << "sumGenWeight:  " << cutFlowGenWeight->GetBinContent(1) << std::endl;
    std::cout << "Processing file " << fil << std::endl;
     for (Int_t i=0; i<nentries[fil]; i++) {
        tree[fil]->GetEntry(i);

        //Data
        if (fil == 0) {
          if (ph_et < 225) continue;
          if (fabs(ph_eta) > 1.4442) continue;
          if (!pass95) continue;
          if (ph_SigmaIEtaIEta < 0.001) continue;
          if (ph_SigmaIPhiIPhi < 0.001) continue;
          if (abs(ph_SeedTime) > 3.) continue;
          if (deltaPhiPUPPImetPho < 2.0) continue;
          if (phoPtOverPUPPImet > 1.4) continue;
          if (puppiMET < 160) continue;
          if (minDeltaPhiPUPPImetJet30 < 0.5) continue;
          if (hasElectron) continue;
          if (hasMuon) continue;
          if (haspixelseed) continue;
          if (eleVeto == 1) continue;
          if (muoVeto == 1) continue;
          if (tauVeto == 1) continue;
        }
       
        //MC
        if ((fil != 0) && (fil !=22) && (fil != 23)) {
          if (genPDGid != 22) continue;
          if (ph_et < 225) continue;
          if (fabs(ph_eta) > 1.4442) continue;
          if (!pass95) continue;
          if (ph_SigmaIEtaIEta < 0.001) continue;
          if (ph_SigmaIPhiIPhi < 0.001) continue;
          if (abs(ph_SeedTime) > 3.) continue;
          if (deltaPhiPUPPImetPho < 2.0) continue;
          if (phoPtOverPUPPImet > 1.4) continue;
          if (puppiMET < 160) continue;
          if (minDeltaPhiPUPPImetJet30 < 0.5) continue;
          if (hasElectron) continue;
          if (hasMuon) continue;
          if (haspixelseed) continue;
          if (eleVeto == 1) continue;
          if (muoVeto == 1) continue;
          if (tauVeto == 1) continue;
        } 

        //efakes
        if (fil == 22) {
          if (ph_et < 225) continue;
          if (fabs(ph_eta) > 1.4442) continue;
          if (!pass95) continue;
          if (ph_SigmaIEtaIEta < 0.001) continue;
          if (ph_SigmaIPhiIPhi < 0.001) continue;
          if (abs(ph_SeedTime) > 3.) continue;
          if (deltaPhiPUPPImetPho < 2.0) continue;
          if (phoPtOverPUPPImet > 1.4) continue;
          if (puppiMET < 160) continue;
          if (minDeltaPhiPUPPImetJet30 < 0.5) continue;
          if (hasElectron) continue;
          if (hasMuon) continue;
          if (!haspixelseed) continue;
          if (eleVeto == 1) continue;
          if (muoVeto == 1) continue;
          if (tauVeto == 1) continue;
        }

        //jetfakes
        if (fil == 23) {
          if (ph_et < 225) continue;
          if (fabs(ph_eta) > 1.4442) continue;
          if (pass95) continue;
          if (!passjetfakes) continue;
          if (ph_SigmaIEtaIEta < 0.001) continue;
          if (ph_SigmaIPhiIPhi < 0.001) continue;
          if (abs(ph_SeedTime) > 3.) continue;
          if (deltaPhiPUPPImetPho < 2.0) continue;
          if (phoPtOverPUPPImet > 1.4) continue;
          if (puppiMET < 160) continue;
          if (minDeltaPhiPUPPImetJet30 < 0.5) continue;
          if (hasElectron) continue;
          if (hasMuon) continue;
          if (eleVeto == 1) continue;
          if (muoVeto == 1) continue;
          if (tauVeto == 1) continue;
        }

        //Fill histograms here 
        if (fil==0) {
          hist[fil][0]->Fill(fabs(ph_eta));
          hist[fil][1]->Fill(ph_phi);
          hist[fil][2]->Fill(nPUPPIJet30);
          hist[fil][3]->Fill(PUPPIJetHt30);
        } else if ((fil != 0) && (fil !=22) && (fil != 23)) {
          hist[fil][0]->Fill(fabs(ph_eta), (puWeight/sumGenWeight)*getSF(ph_eta, ph_et));
          hist[fil][1]->Fill(ph_phi, (puWeight/sumGenWeight)*getSF(ph_eta, ph_et));
          hist[fil][2]->Fill(nPUPPIJet30, (puWeight/sumGenWeight)*getSF(ph_eta, ph_et));
          hist[fil][3]->Fill(PUPPIJetHt30, (puWeight/sumGenWeight)*getSF(ph_eta, ph_et));
        } else if (fil == 22) {
          hist[fil][0]->Fill(fabs(ph_eta), eleFakeWeight(ph_eta));
          hist[fil][1]->Fill(ph_phi, eleFakeWeight(ph_eta));
          hist[fil][2]->Fill(nPUPPIJet30, eleFakeWeight(ph_eta));
          hist[fil][3]->Fill(PUPPIJetHt30, eleFakeWeight(ph_eta));
        } else if (fil == 23) {
          hist[fil][0]->Fill(fabs(ph_eta), jetFakeWeight(ph_et));
          hist[fil][1]->Fill(ph_phi, jetFakeWeight(ph_et));
          hist[fil][2]->Fill(nPUPPIJet30, jetFakeWeight(ph_et));
          hist[fil][3]->Fill(PUPPIJetHt30, jetFakeWeight(ph_et));
        }
     }
   }

    TString title[nhist] = {"ph_eta", "ph_phi", "nPUPPIjet", "PUPPIjetHT"};
    TString title_axis[nhist] = {"#eta^{#gamma}", "#phi^{#gamma}", "nPUPPI jets", "H_{T}^{PUPPI30}"};

   /*TString title[nhist] = {"ph_BDTpred", "ph_et", "ph_eta", "ph_PFECALClusIsoCorr", "ph_PFHCALClusIsoCorr", "ph_TkrIsoCorr", "ph_hoe", "ph_sieie", "ph_sieip", "ph_sipip", "ph_2x2OE3x3Full5x5", "ph_E1x3OESCrFull5x5", "ph_E2ndOESCrFull5x5", "ph_E2x5OESCrFull5x5", "ph_EmaxOESCrFull5x5", "ph_EtaWidth", "ph_PhiWidth", "ph_EtaWOPhiWFull5x5", "ph_R9Full5x5", "ph_S4Full5x5", "ph_SigmaIEtaIEta_o_ph_SigmaIPhiIPhi", "zmass", "ele_pt", "ele_eta"};
   TString title_axis[nhist] = {"BDT score", "p_{T}^{#gamma}", "#eta^{#gamma}", "PF ECALClusIsoCorr", "PF HCALClusIsoCorr", "TkrIsoCorr", "H/E", "#sigma_{i#eta i #eta}", "#sigma_{i#eta i#phi}","#sigma_{i#phi i #phi}", "E_{2x2}/E_{3x3}","E_{1x3}/E^{raw}_{SC}","E_{2}/E^{raw}_{SC}","E_{2x5}/E^{raw}_{SC}","E_{max}/E^{raw}_{SC}","#sigma_{#eta}", "#sigma_{#phi}", "#sigma_{#eta}/#sigma_{#phi}", "S4","R9", "#sigma_{i#eta i #eta}/#sigma_{i#phi i #phi}", "e#gamma invariant mass", "p_{T}^{e}", "#eta^{e}"};*/

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

   Float_t overflow[nfile][nhist]; 
   Float_t underflow[nfile][nhist];

   for (Int_t fil=0; fil<nfile_eff; fil++) {
      for (int i=0; i<nhist; i++) {
         overflow[fil][i] = hist[fil][i]->GetBinContent(hist[fil][i]->GetNbinsX()+1);
         underflow[fil][i] = hist[fil][i]->GetBinContent(0);
         hist[fil][i]->SetBinContent(1, hist[fil][i]->GetBinContent(1)+underflow[fil][i]);
         hist[fil][i]->SetBinContent(hist[fil][i]->GetNbinsX(), hist[fil][i]->GetBinContent(hist[fil][i]->GetNbinsX())+overflow[fil][i]);
      }
   }

   Float_t data_ev[nhist] = {0};
   Float_t ZvvG_ev[nhist] = {0};
   Float_t WlvG_ev[nhist] = {0};
   Float_t ZllG_ev[nhist] = {0};
   Float_t GG_ev[nhist] = {0};
   Float_t Gjets_ev[nhist] = {0};
   Float_t QCD_ev[nhist] = {0};
   Float_t TTg_ev[nhist] = {0};
   Float_t VV_ev[nhist] = {0};
   Float_t Wlv_ev[nhist] = {0};
   Float_t Efake_ev[nhist] = {0};
   Float_t Jfake_ev[nhist] = {0};

   for (int i=0; i<nhist; i++) {

      h_all[i] = (TH1F*)hist[1][i]->Clone("h_all");
      for (int j=2; j<nfile_eff; j++) { //sum signal+bkgs correctly normalized   
        h_all[i]->Add(hist[j][i]);
      }

      Float_t scale = hist[0][i]->Integral()/h_all[i]->Integral(); //scale for data
      for (int j=1; j<nfile_eff; j++) {
        hist[j][i]->Scale(scale);
      }

      h_allScaled[i] = (TH1F*)hist[1][i]->Clone("h_allScaled");
      for (int j=2; j<nfile_eff; j++) { //sum signal+bkgs correctly normalized
        h_allScaled[i]->Add(hist[j][i]);
      }
      h_ratio[i] = (TH1F*)hist[0][i]->Clone("h_ratio");
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
      hist[0][i]->GetYaxis()->SetRangeUser(0.01, 1000000);

      for (Int_t ibin=1; ibin<=hist[0][i]->GetNbinsX() ; ibin++) {
        data_ev[i]+=hist[0][i]->GetBinContent(ibin);
      }  

      if (is2017 || is2018) {
         for (int j=1; j<nfile_eff; j++) {

           hist[j][i]->SetTitle("");
           hist[j][i]->GetYaxis()->SetRangeUser(0.01, 1000000);

           if (j == 1) { //Z(vv)gamma 
             hist[j][i]->SetFillColor(kOrange-4);
             hist[j][i]->SetLineColor(kOrange-4);
             for (Int_t ibin=1; ibin<=hist[j][i]->GetNbinsX() ; ibin++) {
                ZvvG_ev[i]+=hist[j][i]->GetBinContent(ibin);
             }
           }

           if (j == 2) { //W(lv)gamma
             hist[j][i]->SetFillColor(kMagenta-8);
             hist[j][i]->SetLineColor(kMagenta-8);
             for (Int_t ibin=1; ibin<=hist[j][i]->GetNbinsX() ; ibin++) {
                WlvG_ev[i]+=hist[j][i]->GetBinContent(ibin);
             }
           }

           if (j == 3) { //Z(ll)gamma
             hist[j][i]->SetFillColor(kPink+2);
             hist[j][i]->SetLineColor(kPink+2);
             for (Int_t ibin=1; ibin<=hist[j][i]->GetNbinsX() ; ibin++) {
                ZllG_ev[i]+=hist[j][i]->GetBinContent(ibin);
             }
           }

           if (j == 4) { //gg
             hist[j][i]->SetFillColor(kCyan-2);
             hist[j][i]->SetLineColor(kCyan-2);
             for (Int_t ibin=1; ibin<=hist[j][i]->GetNbinsX() ; ibin++) {
                GG_ev[i]+=hist[j][i]->GetBinContent(ibin);
             }
           }

           if (j > 4 && j <= 7) { //gjets
             hist[j][i]->SetFillColor(kBlue-9);
             hist[j][i]->SetLineColor(kBlue-9);
             for (Int_t ibin=1; ibin<=hist[j][i]->GetNbinsX() ; ibin++) {
                Gjets_ev[i]+=hist[j][i]->GetBinContent(ibin);
             }
           }

           if (j > 7 && j <= 14) { //QCD
             hist[j][i]->SetFillColor(kAzure+8); //kBlue-7
             hist[j][i]->SetLineColor(kAzure+8); //kBlue-7
             for (Int_t ibin=1; ibin<=hist[j][i]->GetNbinsX() ; ibin++) {
                QCD_ev[i]+=hist[j][i]->GetBinContent(ibin);
             }
           }

           if (j > 14 && j<= 16) { //tg, ttg
             hist[j][i]->SetFillColor(kSpring+2);
             hist[j][i]->SetLineColor(kSpring+2);
             for (Int_t ibin=1; ibin<=hist[j][i]->GetNbinsX() ; ibin++) {
                TTg_ev[i]+=hist[j][i]->GetBinContent(ibin);
             }
           }

           if (j > 16 && j<= 19) { //WW, WZ, ZZ
             hist[j][i]->SetFillColor(kSpring+4);
             hist[j][i]->SetLineColor(kSpring+4);
             for (Int_t ibin=1; ibin<=hist[j][i]->GetNbinsX() ; ibin++) {
                VV_ev[i]+=hist[j][i]->GetBinContent(ibin);
             }
           }

           if (j > 19 && j <=21 ) { //W(muv)(tauv)
             hist[j][i]->SetFillColor(kAzure-5);
             hist[j][i]->SetLineColor(kAzure-5);
             for (Int_t ibin=1; ibin<=hist[j][i]->GetNbinsX() ; ibin++) {
                Wlv_ev[i]+=hist[j][i]->GetBinContent(ibin);
             }
           }

           if (j == 22) { //e fakes
             hist[j][i]->SetFillColor(kBlue-7);
             hist[j][i]->SetLineColor(kBlue-7);
             for (Int_t ibin=1; ibin<=hist[j][i]->GetNbinsX() ; ibin++) {
                Efake_ev[i]+=hist[j][i]->GetBinContent(ibin);
             }
           }

           if (j == 23) { //jet fakes
             hist[j][i]->SetFillColor(kRed);
             hist[j][i]->SetLineColor(kRed);
             for (Int_t ibin=1; ibin<=hist[j][i]->GetNbinsX() ; ibin++) {
                Jfake_ev[i]+=hist[j][i]->GetBinContent(ibin);
             }
           }

           if (j>7 && j <=14)  hs[i]->Add(hist[j][i]);
         }

         hs[i]->Add(hist[20][i]);
         hs[i]->Add(hist[21][i]);
         hs[i]->Add(hist[3][i]);
         hs[i]->Add(hist[17][i]);
         hs[i]->Add(hist[18][i]);
         hs[i]->Add(hist[19][i]);
         hs[i]->Add(hist[15][i]);
         hs[i]->Add(hist[16][i]);
         hs[i]->Add(hist[23][i]);
         hs[i]->Add(hist[4][i]);
         hs[i]->Add(hist[5][i]);
         hs[i]->Add(hist[6][i]);
         hs[i]->Add(hist[7][i]);
         hs[i]->Add(hist[2][i]);
         hs[i]->Add(hist[22][i]);
         hs[i]->Add(hist[1][i]);
      }

      if (is2016_preVFP || is2016_postVFP) {
         for (int j=1; j<nfile_eff; j++) {
           hist[j][i]->SetTitle("");

           if (j == 1) { //DY
             hist[j][i]->SetFillColor(kOrange-4);
             hist[j][i]->SetLineColor(kOrange-4);
           }

           if (j > 1 && j <= 4) { //GJets
           hist[j][i]->SetFillColor(kMagenta-8);
           hist[j][i]->SetLineColor(kMagenta-8);
           }

           if (j > 4 && j <= 11) { //QCD
           hist[j][i]->SetFillColor(kPink+2);
           hist[j][i]->SetLineColor(kPink+2);
           }

           if (j == 12) { //WGLnu
           hist[j][i]->SetFillColor(kCyan-2);
           hist[j][i]->SetLineColor(kCyan-2);
           }

           if (j == 13) { //TG, TTG
           hist[j][i]->SetFillColor(kBlue-9);
           hist[j][i]->SetLineColor(kBlue-9);
           }

           if (j == 14) { //WW
           hist[j][i]->SetFillColor(kTeal-6); //kBlue-7
           hist[j][i]->SetLineColor(kTeal-6); //kBlue-7
           }

           if (j == 15) { //WZ
           hist[j][i]->SetFillColor(kSpring+2);
           hist[j][i]->SetLineColor(kSpring+2);
           }

           if (j == 16) { //ZZ
           hist[j][i]->SetFillColor(kSpring+4);
           hist[j][i]->SetLineColor(kSpring+4);
           }

           if (j > 16 && j <= 22) { //Wjets
           hist[j][i]->SetFillColor(kAzure-5);
           hist[j][i]->SetLineColor(kAzure-5);
           }

           if (j == 23) { //Diphoton
           hist[j][i]->SetFillColor(kBlue-7);
           hist[j][i]->SetLineColor(kBlue-7);
           }

           if (j > 1)  hs[i]->Add(hist[j][i]);
         }

         hs[i]->Add(hist[1][i]);
      }

      if (is2016) {
         for (int j=1; j<nfile_eff; j++) {
           hist[j][i]->SetTitle("");

           if (j ==1 || j==2) { //DY
             hist[j][i]->SetFillColor(kOrange-4);
             hist[j][i]->SetLineColor(kOrange-4);
           }

           if (j > 2 && j <= 8) { //GJets
           hist[j][i]->SetFillColor(kMagenta-8);
           hist[j][i]->SetLineColor(kMagenta-8);
           }

           if (j > 8 && j <= 22) { //QCD
           hist[j][i]->SetFillColor(kPink+2);
           hist[j][i]->SetLineColor(kPink+2);
           }

           if (j == 23 || j == 24) { //WGLnu
           hist[j][i]->SetFillColor(kCyan-2);
           hist[j][i]->SetLineColor(kCyan-2);
           }

           if (j == 25 || j == 26) { //TG, TTG
           hist[j][i]->SetFillColor(kBlue-9);
           hist[j][i]->SetLineColor(kBlue-9);
           }
    
           if (j == 27 || j==28) { //WW
           hist[j][i]->SetFillColor(kTeal-6); //kBlue-7
           hist[j][i]->SetLineColor(kTeal-6); //kBlue-7
           }

           if (j == 29 || j==30) { //WZ
           hist[j][i]->SetFillColor(kSpring+2);
           hist[j][i]->SetLineColor(kSpring+2);
           }

           if (j == 31 || j==32) { //ZZ
           hist[j][i]->SetFillColor(kSpring+4);
           hist[j][i]->SetLineColor(kSpring+4);
           }

           if (j > 33 && j <= 45) { //Wjets
           hist[j][i]->SetFillColor(kAzure-5);
           hist[j][i]->SetLineColor(kAzure-5);
           }

           if (j == 46 || j==47) { //Diphoton
           hist[j][i]->SetFillColor(kBlue-7);
           hist[j][i]->SetLineColor(kBlue-7);
           }

           if (j > 2)  hs[i]->Add(hist[j][i]);
         }

         hs[i]->Add(hist[1][i]);
         hs[i]->Add(hist[2][i]);
      
      }
      hist[0][i]->SetMinimum(0.1);
      hist[0][i]->SetMaximum(1000000);
      hist[0][i]->Draw("ep");
      hs[i]->Draw("hist");
      hs[i]->SetMinimum(0.1);
      hs[i]->SetMaximum(1000000);
      h_allScaled[i]->SetFillStyle(3001);
      h_allScaled[i]->SetFillColor(15);
      h_allScaled[i]->GetYaxis()->SetTitle("Events");
      h_allScaled[i]->Draw("e2, SAME");
      h_allScaled[i]->SetMinimum(0.1);
      h_allScaled[i]->SetMaximum(1000000);
      hist[0][i]->Draw("ep, SAME");
      hist[0][i]->SetMinimum(0.1);
      hist[0][i]->SetMaximum(1000000);
      hs[i]->SetTitle("");
      hs[i]->GetYaxis()->SetTitle("Events");

      hist[0][i]->Draw("ep");
      hs[i]->Draw("hist");
      h_allScaled[i]->Draw("e2, SAME");
      hist[0][i]->Draw("ep, SAME");

      TLegend* leg = new TLegend();
      leg = new TLegend(0.88, 0.5, 0.12, 0.891);
      leg->SetBorderSize(0);
      leg->SetEntrySeparation(0.01);
      leg->SetFillColor(0);
      leg->SetFillStyle(0);
      leg->SetNColumns(3);

      string data_string = "Data "+std::to_string(data_ev[0]);
      string ZvvG_string = "Z(#rightarrow #nu #nu) #gamma "+std::to_string(round(ZvvG_ev[0]));
      string WlvG_string = "W (#rightarrow l #nu) #gamma "+std::to_string(round(WlvG_ev[0]));
      string ZllG_string = "Z(#rightarrow ll) #gamma "+std::to_string(round(ZllG_ev[0]));
      string GG_string = "#gamma #gamma "+std::to_string(round(GG_ev[0]));
      string Gjets_string = "#gamma + jets "+std::to_string(round(Gjets_ev[0]));
      string QCD_string = "QCD "+std::to_string(round(QCD_ev[0]));
      string TTg_string = "t#gamma, tt#gamma "+std::to_string(round(TTg_ev[0]));
      string VV_string = "WW, WZ, ZZ "+std::to_string(round(VV_ev[0]));
      string Wlv__string = "W(#rightarrow #mu/#tau #nu) "+std::to_string(round(Wlv_ev[0]));
      string Efakes_string = "e fakes "+std::to_string(round(Efake_ev[0]));
      string Jfakes_string = "jet fakes "+std::to_string(round(Jfake_ev[0]));
     
      if (is2017 || is2018) {
         leg->AddEntry(hist[0][i], data_string.c_str(), "ep");
         leg->AddEntry(h_allScaled[i], "MC stat", "f");
         leg->AddEntry(hist[1][i], ZvvG_string.c_str(), "f");
         leg->AddEntry(hist[2][i], WlvG_string.c_str(), "f");
         leg->AddEntry(hist[3][i], ZllG_string.c_str(), "f");
         leg->AddEntry(hist[4][i], GG_string.c_str(),  "f");
         leg->AddEntry(hist[5][i], Gjets_string.c_str(), "f");
         leg->AddEntry(hist[8][i], QCD_string.c_str(), "f");
         leg->AddEntry(hist[15][i],TTg_string.c_str(), "f");
         leg->AddEntry(hist[17][i],VV_string.c_str(), "f");
         leg->AddEntry(hist[20][i],Wlv__string.c_str(), "f");
         leg->AddEntry(hist[22][i],Efakes_string.c_str(), "f");
         leg->AddEntry(hist[23][i],Jfakes_string.c_str(), "f");
      }

      if (is2016_preVFP || is2016_postVFP) {
         leg->AddEntry(hist[0][i], "Data", "ep");
         leg->AddEntry(h_allScaled[i], "MC stat", "f");
         leg->AddEntry(hist[1][i], "Z(#rightarrow ee)+jets", "f");
         leg->AddEntry(hist[2][i], "#gamma +jets", "f");
         leg->AddEntry(hist[5][i], "QCD", "f");
         leg->AddEntry(hist[13][i], "tt#gamma", "f");
         leg->AddEntry(hist[12][i], "W(#rightarrow l#nu)+#gamma", "f");
         leg->AddEntry(hist[14][i], "WW", "f");
         leg->AddEntry(hist[15][i], "WZ", "f");
         leg->AddEntry(hist[16][i], "ZZ", "f");
         leg->AddEntry(hist[17][i], "WJets", "f");
         leg->AddEntry(hist[23][i], "#gamma #gamma", "f");
      }

      if (is2016) {
         leg->AddEntry(hist[0][i], "Data", "ep");
         leg->AddEntry(h_allScaled[i], "MC stat", "f");
         leg->AddEntry(hist[1][i], "Z(#rightarrow ee)+jets", "f");
         leg->AddEntry(hist[3][i], "#gamma +jets", "f");
         leg->AddEntry(hist[9][i], "QCD", "f");
         leg->AddEntry(hist[25][i], "tt#gamma", "f");
         leg->AddEntry(hist[23][i], "W(#rightarrow l#nu)+#gamma", "f");
         leg->AddEntry(hist[27][i], "WW", "f");
         leg->AddEntry(hist[29][i], "WZ", "f");
         leg->AddEntry(hist[31][i], "ZZ", "f");
         leg->AddEntry(hist[32][i], "WJets", "f");
         leg->AddEntry(hist[46][i], "#gamma #gamma", "f");
      }

      c1->Update();
      leg->Draw("SAME");

      pad1->SetLogy();

      Int_t year = 17;
      if (is2018) {year=18;}
      if (is2016) {year=16;}
      if (is2016_preVFP) {year=161;}
      if (is2016_postVFP) {year=162;}
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


   //background fraction
   /*Double_t data_ev[nhist] = {0};
   Double_t ZvvG_ev[nhist] = {0};
   Double_t WlvG_ev[nhist] = {0};
   Double_t ZllG_ev[nhist] = {0};
   Double_t GG_ev[nhist] = {0};
   Double_t Gjets_ev[nhist] = {0};
   Double_t QCD_ev[nhist] = {0};
   Double_t TTg_ev[nhist] = {0};
   Double_t VV_ev[nhist] = {0};
   Double_t Wlv_ev[nhist] = {0};
   Double_t Efake_ev[nhist] = {0};
   Double_t Jfake_ev[nhist] = {0};

   for (int i=0; i<nhist; i++) {
    for (int j=0; j<nfile_eff; j++) {
       for (Int_t ibin=1; ibin<=hist[j][i]->GetNbinsX() ; ibin++) {
          if (j == 0) { //data
             data_ev[i]+=hist[j][i]->GetBinContent(ibin);
          }
          if (j == 1) { //Z(vv)gamma 
             ZvvG_ev[i]+=hist[j][i]->GetBinContent(ibin);
          }
          if (j == 2) { //W(lv)gamma
             WlvG_ev[i]+=hist[j][i]->GetBinContent(ibin);
          }
          if (j == 3) { //Z(ll)gamma
             ZllG_ev[i]+=hist[j][i]->GetBinContent(ibin);
          }
          if (j == 4) { //gg
             GG_ev[i]+=hist[j][i]->GetBinContent(ibin);
          }
          if (j > 4 && j <= 7) { //gjets
             Gjets_ev[i]+=hist[j][i]->GetBinContent(ibin);
          }
          if (j > 7 && j <= 14) { //QCD
             QCD_ev[i]+=hist[j][i]->GetBinContent(ibin);
          }
          if (j > 14 && j<= 16) { //tg, ttg
             TTg_ev[i]+=hist[j][i]->GetBinContent(ibin);
          }
          if (j > 16 && j<= 19) { //WW, WZ, ZZ
             VV_ev[i]+=hist[j][i]->GetBinContent(ibin);
          }
          if (j > 19 && j <=21 ) { //W(muv)(tauv)
             Wlv_ev[i]+=hist[j][i]->GetBinContent(ibin);
          }
          if (j == 22) { //e fakes
             Efake_ev[i]+=hist[j][i]->GetBinContent(ibin);
          }
          if (j == 23) { //jet fakes
             Jfake_ev[i]+=hist[j][i]->GetBinContent(ibin);
          }
       }
     }
   }*/

   std::cout << "ZvvG: " << ZvvG_ev[0] << std::endl;
   std::cout << "WlvG: " << WlvG_ev[0] << std::endl;
   std::cout << "ZllG: " << ZllG_ev[0] << std::endl;
   std::cout << "GG: " << GG_ev[0] << std::endl;
   std::cout << "Gjets: " << Gjets_ev[0] << std::endl;
   std::cout << "QCD: " << QCD_ev[0] << std::endl;
   std::cout << "TTg:  " << TTg_ev[0] << std::endl;
   std::cout << "VV: " << VV_ev[0] << std::endl;
   std::cout << "Wlv: " << Wlv_ev[0] << std::endl;
   std::cout << "Efake: " << Efake_ev[0] << std::endl;
   std::cout << "Jfake: " << Jfake_ev[0] << std::endl;
   
   /*const int num_bkg = 8;
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
   }*/
      


}
