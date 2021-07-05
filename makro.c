#include <TH1.h>
#include <TF1.h>
#include <TFile.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>

#include "RooRandom.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooCBShape.h"
#include "RooBreitWigner.h"
#include "RooArgusBG.h"
#include "RooExponential.h"
#include "RooPolynomial.h"
#include "RooExtendPdf.h"
#include "RooHistPdf.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooPlot.h"
#include <TROOT.h>
#include "TLatex.h"
#include <TCanvas.h>
   using namespace RooFit;

void makro(std::string inFileName){
  
  TFile *inFile=TFile::Open(inFileName.c_str(), "READ");

  TTree* tree = (TTree*) inFile->Get("DecayTree");

  Double_t A_MM             = 0;
  Double_t A_M              = 0;
  Double_t A_PT             = 0;
  Double_t muplus_PIDmu     = 0;
  Double_t muminus_PIDmu    = 0;
  Double_t muplus_PT        = 0;
  Double_t muminus_PT       = 0;
  Double_t muplus_ProbNNmu  = 0;
  Double_t muminus_ProbNNmu = 0; 
  Double_t muminus_IP_OWNPV = 0;
  Double_t muplus_IP_OWNPV  = 0;

  tree->SetBranchAddress("A_MM",&A_MM);
  tree->SetBranchAddress("A_M",&A_M);
  tree->SetBranchAddress("A_PT",&A_PT);
  tree->SetBranchAddress("mup_PIDmu",&muplus_PIDmu);
  tree->SetBranchAddress("mun_PIDmu",&muminus_PIDmu);
  tree->SetBranchAddress("mup_PT",&muplus_PT);
  tree->SetBranchAddress("mun_PT",&muminus_PT);
  tree->SetBranchAddress("mup_IP_OWNPV",&muplus_IP_OWNPV);
  tree->SetBranchAddress("mun_IP_OWNPV",&muminus_IP_OWNPV);

  tree->SetBranchAddress("mup_ProbNNmu",&muplus_ProbNNmu);
  tree->SetBranchAddress("mun_ProbNNmu",&muminus_ProbNNmu);

  TH1D* h_upsilon = new TH1D("Upsilon","",200,8000,12000);
  TH1D* h_Z = new TH1D("Z","",100,10000,150000);

  Long64_t num = tree->GetEntries();
  for(Long64_t entry =0; entry < tree->GetEntries(); ++entry){
    if (entry%1000000 == 0) cout<<"Progress: "<<entry<<" of "<<num<<endl;
    tree->GetEntry(entry);
    if(muplus_PIDmu > 5 && muminus_PIDmu > 5 && muminus_ProbNNmu > 0.95 && muplus_ProbNNmu > 0.95){
      if(muplus_PT > 10 && muminus_PT > 10){
        //if(muplus_IP_OWNPV > 0.1 && muminus_IP_OWNPV > 0.1){
          h_upsilon->Fill(A_MM);
          h_Z->Fill(A_MM);
        //}
      }
    }
  }

  h_upsilon->Draw("hist");
  h_upsilon->GetXaxis()->SetTitle("#Upsilon mass [MeV]");
  h_upsilon->GetXaxis()->SetTitleOffset(1.5);
  h_upsilon->GetYaxis()->SetTitle("Candidates");

  h_Z->Draw("hist");
  h_Z->GetXaxis()->SetTitle("Z mass [MeV]");
  h_Z->GetXaxis()->SetTitleOffset(1.5);
  h_Z->GetYaxis()->SetTitle("Candidates");

  gStyle->SetOptStat(0);
  gPad->SetLogx();
  gPad->SetLogy();

  TFile *outHistFile = TFile::Open("Plot.root","RECREATE");
  outHistFile->cd();
  h_upsilon->Write();
  h_Z->Write();
  outHistFile->Close();

}

