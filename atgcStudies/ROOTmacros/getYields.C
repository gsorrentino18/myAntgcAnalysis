#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <iostream>
#include <fstream>
#include <THStack.h>
#include "CMS_lumi.C"
#include "rebin.h"
#include <TStyle.h>

using namespace std;

void getYields(){

   //bool DEBUG = true;
   bool DEBUG = false;
   
   const Int_t nfile=2;
   TFile *file[nfile];
   TString atgc = "Zggh34p0em4h4m5p0em7";
   //TString atgc = "Zggh30p0h45p0em7";
   //TString atgc = "Zggh34p0em4h40p0";
   
   file[0] = new TFile("/local/cms/user/gsorrent/antgc_analysis/atgcStudies/22Jun/ZNuNuGJetsMonoPhotonPtG130TuneCP513TeVmadgraph.root");
   file[1] = new TFile("/local/cms/user/gsorrent/antgc_analysis/atgcStudies/22Jun/ZnunugammaaTGC"+atgc+".root");

   Float_t         puWeight;
   Float_t         genMET;
   Float_t         ph_pt;
   Float_t         ph_eta;
   Float_t         ph_phi;
   Float_t         ph_e;
   Float_t         ph_et;

   TTree *tree[nfile];
   for (Int_t i=0; i<nfile;i++) {
     tree[i] = (TTree*)file[i]->Get("atgc/fitter_tree");
     tree[i]->SetBranchAddress("puWeight",&puWeight);
     tree[i]->SetBranchAddress("genMET",&genMET);
     tree[i]->SetBranchAddress("ph_pt",&ph_pt);
     tree[i]->SetBranchAddress("ph_eta",&ph_eta);
     tree[i]->SetBranchAddress("ph_phi", &ph_phi);
     tree[i]->SetBranchAddress("ph_e",&ph_e);
     tree[i]->SetBranchAddress("ph_et",&ph_et);
   }

   if (DEBUG) std::cout << "histos" << std::endl;

   const Int_t nhist=4;
   TH1F* hist[nfile][nhist];
   THStack *hs[nhist];

   gStyle->SetOptStat(false);

   hist[0][0]   = new TH1F("phoPtEB","phoPtEB",10,225,2000);
   hist[0][1]   = new TH1F("phoPtEE","phoPtEE",10,225,2000);
   hist[0][2]   = new TH1F ("phoEta","phoEta",5,1.566,2.5);
   hist[0][3]   = new TH1F ("phoPhi","phoPhi",5,0,3.14);

   hist[1][0]   = new TH1F("phoPt_atgcEB","phoPt_atgcEB",10,225,2000);
   hist[1][1]   = new TH1F("phoPt_atgcEE","phoPt_atgcEE",10,225,2000);
   hist[1][2]   = new TH1F("phoEta_atgc","phoEta_atgc",5,1.566,2.5);
   hist[1][3]   = new TH1F ("phoPhi_atgc","phoPhi_atgc",5,0,3.14);

   if (DEBUG) std::cout << "entries" << std::endl;

   Int_t nentries[nfile] = {0};
   for (Int_t i=0; i<nfile;i++) {
     nentries[i] = (Int_t)tree[i]->GetEntries();
   }

   for (Int_t fil=0; fil<nfile; fil++) {

     for (Int_t i=0; i<nentries[fil]; i++) {
        tree[fil]->GetEntry(i);

        if ( ph_pt < 225.) continue;
        //if ( genMET < 150.) continue;

        if (std::abs(ph_eta) < 1.4442) {
          hist[fil][0]->Fill(ph_pt, puWeight);
        } else if (std::abs(ph_eta) < 2.5 || std::abs(ph_eta) > 1.566) {
          hist[fil][1]->Fill(ph_pt, puWeight);
        }
        hist[fil][2]->Fill(ph_eta, puWeight);
     }
   }

   float yield[nfile][2];
   for (int i=0; i<nfile; i++) {
     for (int j=0; j<2; j++) {
       yield[i][j]=0;
     }
   }
   for (int fil=0; fil<nfile; fil++) {
     for (int pos=0; pos<2; pos++) {
       for (int i=0; i<=hist[fil][pos]->GetNbinsX(); i++) {
         yield[fil][pos]+=hist[fil][pos]->GetBinContent(i);
       }
     }
   }

   std::cout << "yield signal EB:  " << yield[0][0] << "  yield signal EE:  " << yield[0][1] << std::endl;
   std::cout << "yield atgc EB:  " << yield[1][0] << "  yield atgc EE:  " << yield[1][1] << std::endl;
   std::cout << "atgc/signal EB: " << yield[1][0]/yield[0][0] << "  atgc/signal EE: " << yield[1][1]/yield[0][1] << std::endl;

   TH1F* h_ratio[nhist];
   TH1F* h_total[nhist];

   TString title[nhist] = {"ph_ptEB","ph_ptEE", "ph_eta", "ph_phi"};
   TString title_axis[nhist] = {"p_{T}^{#gamma}", "p_{T}^{#gamma}", "#eta", "#phi"};

   for (int i=0; i<nhist; i++) {

      hist[0][i] = rebin(hist[0][i]);
      hist[1][i] = rebin(hist[1][i]);

      h_ratio[i] = (TH1F*)hist[1][i]->Clone("h_ratio"); //atgc
      h_ratio[i]->Sumw2(1);
      h_total[i] = (TH1F*)hist[0][i]->Clone("h_total"); //SM

      /*if (nfile > 2) {
         for (Int_t j=2; j<nfile; j++) {
           h_total[i]->Add(hist[j][i]);
         }
      }*/
      h_ratio[i]->Divide(h_total[i]);

      TCanvas *c1 = new TCanvas("c1", "c1");
      c1->cd();
      
      TPad* pad1 = new TPad("pad1", "pad1", 0.0, 0.3, 1.0, 1.0);
      pad1->SetBottomMargin(0);
      pad1->Draw();
      pad1->cd();

      hs[i] = new THStack(title[i]+"_stack", "");

      hist[0][i]->SetLineColor(kRed); //data
      hist[1][i]->SetLineColor(kBlue);

      for (Int_t j=1; j<nfile; j++) {
        hs[i]->Add(hist[j][i]);
      }

      TLegend* leg = new TLegend();
      leg = new TLegend(0.68, 0.5, 0.9, 0.891);
      leg->SetBorderSize(0);
      leg->SetEntrySeparation(0.01);
      leg->SetFillColor(0);
      leg->SetFillStyle(0);

      leg->AddEntry(hist[0][i], "SM", "l");
      leg->AddEntry(hist[1][i], atgc, "l");

      hist[0][i]->SetTitle("");
      hist[0][i]->Draw("hist");
      hs[i]->Draw("hist, SAME");
      leg->Draw("SAME");
      hs[i]->SetTitle("");
      hs[i]->GetYaxis()->SetTitle("Events");

      pad1->SetLogy();

      CMS_lumi(pad1, 17, 0);
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
      h_ratio[i]->GetYaxis()->SetTitle("aTGC/SM");
      h_ratio[i]->GetYaxis()->SetTitleSize(0.11);
      h_ratio[i]->GetYaxis()->SetTitleOffset(0.43);
      h_ratio[i]->GetYaxis()->SetLabelSize(0.1);
      h_ratio[i]->GetYaxis()->SetLabelOffset(0.01);
      h_ratio[i]->GetYaxis()->SetNdivisions(505);
      h_ratio[i]->GetYaxis()->SetRangeUser(0.5, 1.5);

      TLine* line = new TLine(h_ratio[i]->GetXaxis()->GetXmin(), 1.0, h_ratio[i]->GetXaxis()->GetXmax(), 1.0);
      line->SetLineColor(kRed);
      line->SetLineWidth(1);
      h_ratio[i]->Draw("ep");
      line->Draw("SAME");
      c1->Update();

      c1->SaveAs(title[i]+".png");
   }










}
