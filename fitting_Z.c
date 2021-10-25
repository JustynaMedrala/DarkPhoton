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

void fitting_Z(){
	
  TFile* input_file = new TFile("result_high_mass_less.root", "OPEN");

  TH1* mass_data = dynamic_cast<TH1*>(input_file->Get("high mass") ); 

  TCanvas* fitted=new TCanvas("Z","Z",0,0,900,700);

  Double_t min_mass = 35000;
  Double_t max_mass = 120000;
  Double_t mass_mean = 90000;
  Double_t signal_min = 82000;
  Double_t signal_max = 95000;

  RooRealVar Photon_M("A_M", "A_M",mass_mean, min_mass, max_mass);
  RooDataHist data("data","data",Photon_M,mass_data);

  gPad->SetLeftMargin(0.13);

  RooRealVar mean("mean","mean", mass_mean, signal_min, signal_max);
  RooRealVar sigma("sigma","sigma", 700, 1e2, 2e4);
  RooRealVar width("width","width", 1000, 500, 1e5);
  RooVoigtian fit("v","signal pdf",Photon_M,mean,width,sigma);

  RooRealVar lambda("lambda", "slope", -8e-5, -0.1, 0.0);
  RooExponential expo("exponential", "exponential PDF", Photon_M, lambda);

  RooRealVar nsig("nsig","#signal events",2000,1e3,1e5);
  RooRealVar nbkg("nbkg","#background events",500,500,1e5);

  Photon_M.setRange("signal", signal_min, signal_max); 
  Photon_M.setRange("all", min_mass, max_mass);

  RooAddPdf sum("sum","sum" , RooArgList(expo, fit), RooArgList(nbkg, nsig));
  sum.fitTo(data, RooFit::Extended(), Range("all"));

  RooPlot* frame = Photon_M.frame(RooFit::Title("  "));

  data.plotOn(frame, Name("Data"));

  sum.plotOn(frame, Name("sum"));
  sum.plotOn(frame, Components("v"), LineColor(kRed), LineStyle(9), LineWidth(2), Range(20000, 120000), Name("fit"));
  sum.plotOn(frame, Components("exponential"), LineColor(kGreen), LineStyle(10), LineWidth(4), Name("bkg"));
  frame->GetXaxis()->SetTitle("m(#mu^{+}#mu^{-}) [MeV/c^{2}]");
  frame->GetYaxis()->SetTitle("Candidates/ 1400 [MeV/c^{2}]");
  frame->Draw();

  TLegend *leg1 = new TLegend(0.20, 0.65, 0.45, 0.89);
  leg1->SetFillColor(kWhite);
  leg1->SetLineColor(kWhite);
  leg1->AddEntry("Data","Data","LP");
  leg1->AddEntry("sum","Signal + background","LP");
  leg1->AddEntry("bkg","Background only", "LP");
  leg1->AddEntry("fit","Signal only", "LP");
  leg1->Draw();

  TLatex *t = new TLatex(0.22, 0.62,"pp #sqrt{s}=13 TeV");
  t->SetTextSize(0.03);
  t->SetNDC(kTRUE);
  t->Draw();

  fitted->SaveAs("Z.png");

  TH1 *hdata = data.createHistogram("hdata",Photon_M,Binning((max_mass-min_mass)/1400,min_mass, max_mass));
  TH1 *hbkg = expo.createHistogram("hbkg",Photon_M,Binning((max_mass-min_mass)/1400,min_mass, max_mass));

  TCanvas* no_bkg=new TCanvas("Without bkg","  ",0,0,900,700);
  gPad->SetLeftMargin(0.15);
  gStyle->SetOptStat(0);
  //gPad->SetLogy();
  //gPad->SetLogx();
  hbkg->Sumw2();
  hbkg->Scale(nbkg.getVal());
  hdata->Add(hbkg, -1);
  hdata->Draw("p");
  hdata->SetTitle("   ");
  hdata->GetXaxis()->SetTitle("m(#mu^{+}#mu^{-}) [MeV/c^{2}]");
  hdata->GetYaxis()->SetTitle("Data - Background events/1400 [MeV/c^{2}]");
  no_bkg->SaveAs("Z_without_bkg.png");

}
