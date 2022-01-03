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
#include "RooVoigtian.h"
#include <vector>
#include <cmath>

using namespace RooFit;
using namespace std;
using namespace TMath;

void fitting_upsilon(){
	
  TFile* input_file = new TFile("result_high_mass.root", "OPEN");

  TH1* mass_data = dynamic_cast<TH1*>(input_file->Get("high mass") ); 

                                //wartości używane dla 3 pików
  Double_t min_mass = 8000;     //1 - 8000  //2 - 9600  //3 - 10150
  Double_t max_mass = 12000;    //1 - 9800  //2 - 10180 //3 - 12000
  Double_t mass_mean_u1 = 9450;    //1 - 9450  //2 - 9900  //3 - 10250
  Double_t signal_min_u1 = 9340;   //1 - 9340  //2 - 9900  //3 - 10150
  Double_t signal_max_u1 = 9600;  //1 - 9600  //2 - 10100 //3 - 10450
  Double_t mass_mean_u2 = 9900;
  Double_t signal_min_u2 = 9900;
  Double_t signal_max_u2 = 10100;
  Double_t mass_mean_u3 = 10250;
  Double_t signal_min_u3 = 10150;
  Double_t signal_max_u3 = 10450;  

  RooRealVar Photon_M("A_M", "A_M",mass_mean_u1, min_mass, max_mass);
  RooDataHist data("data","data",Photon_M,mass_data);

  RooRealVar mean_u1("mean_u1","mean", mass_mean_u1, signal_min_u1, signal_max_u1);
  RooRealVar mean_u2("mean_u2","mean", mass_mean_u2, signal_min_u2, signal_max_u2);
  RooRealVar mean_u3("mean_u3","mean", mass_mean_u3, signal_min_u3, signal_max_u3);
  RooRealVar sigma_u1("sigma_u1","sigma", 30, 10, 50);
  RooRealVar sigma_u2("sigma_u2","sigma", 40, 10, 50);
  RooRealVar sigma_u3("sigma_u3","sigma", 40, 10, 50);
  RooRealVar width_u1("width_u1","width", 20, 5, 50);
  RooRealVar width_u2("width_u2","width", 10, 5, 50);
  RooRealVar width_u3("width_u3","width", 40, 5, 60);
  RooVoigtian fit_u1("#Upsilon_1","signal pdf",Photon_M,mean_u1,width_u1,sigma_u1);
  RooVoigtian fit_u2("#Upsilon_2","signal pdf",Photon_M,mean_u2,width_u2,sigma_u2);
  RooVoigtian fit_u3("#Upsilon_3","signal pdf",Photon_M,mean_u3,width_u3,sigma_u3);

  RooRealVar nsig_u1("nsig_u1","#signal events",63000,1e2,1e6);
  RooRealVar nsig_u2("nsig_u2","#signal events",15000,1e2,1e5);
  RooRealVar nsig_u3("nsig_u3","#signal events",750,100,1e5);
  RooRealVar nbkg("nbkg","#background events",93000,1e3,1e8);

  RooRealVar lambda("lambda", "slope", 0., -1., 0.1);
  RooExponential expo("exponential", "exponential PDF", Photon_M, lambda);

  Photon_M.setRange("signal_u1", signal_min_u1, signal_max_u1);
  Photon_M.setRange("signal_u2", signal_min_u2, signal_max_u2);
  Photon_M.setRange("signal_u3", signal_min_u3, signal_max_u3);
  Photon_M.setRange("all", min_mass, max_mass);

  RooAddPdf sum("sum","sum" , RooArgList(expo, fit_u1, fit_u2, fit_u3), RooArgList(nbkg, nsig_u1, nsig_u2, nsig_u3)); 

  //---------------------------------------------------//
  //FITOWANIE------------------------------------------//

  sum.fitTo(data, RooFit::Extended(), Range("all"));

  TCanvas* fitted=new TCanvas("Upsilon","Upsilon",0,0,900,700);
  gPad->SetLeftMargin(0.15);

  RooPlot* frame = Photon_M.frame(RooFit::Title("  "));

  data.plotOn(frame, Name("Data"), DataError(RooAbsData::SumW2), MarkerSize(0.7));

  sum.plotOn(frame, Components("sum"), LineWidth(2), Name("sum"));
  sum.plotOn(frame, Components("#Upsilon_1"), LineColor(kRed), LineStyle(9), LineWidth(2), Range(9000, 9800), Name("fit_u1"));
  sum.plotOn(frame, Components("#Upsilon_2"), LineColor(kRed), LineStyle(9), LineWidth(2), Range(9800, 10200), Name("fit_u2"));
  sum.plotOn(frame, Components("#Upsilon_3"), LineColor(kRed), LineStyle(9), LineWidth(2), Range(10100, 10800), Name("fit_u3"));
  sum.plotOn(frame, Components("exponential"), LineColor(kGreen), LineStyle(10), LineWidth(4), Name("bkg"));
  frame->GetXaxis()->SetTitle("m(#mu^{+}#mu^{-}) [MeV/c^{2}]");
  frame->GetYaxis()->SetTitle("Candidates/50 [MeV/c^{2}]");
  frame->Draw();

  //---------------------------------------------------//
  //TEST CHI KWADRAT DLA DANYCH------------------------//

  cout<<"Chi kwadrat: "<<frame->chiSquare("sum", "Data", 17)<<", liczba binów:"<<mass_data->GetNbinsX()<<", "<<(max_mass-min_mass)/mass_data->GetBinWidth(1)<<endl;

  gStyle->SetOptStat(0);

  TLegend *leg1 = new TLegend(0.65, 0.65, 0.89, 0.89);
  leg1->SetFillColor(kWhite);
  leg1->SetLineColor(kWhite);
  leg1->AddEntry("Data","Data","LP");
  leg1->AddEntry("sum","Signal + background","LP");
  leg1->AddEntry("bkg","Background only", "LP");
  leg1->AddEntry("fit_u1","Signal only", "LP");
  leg1->Draw();

  TLatex *t = new TLatex(0.66, 0.62,"pp #sqrt{s}=13 TeV");
  t->SetTextSize(0.03);
  t->SetNDC(kTRUE);
  t->Draw();

  fitted->SaveAs("upsilon_50.png");

  RooHist* hresid = frame->residHist();
  RooHist* hpull = frame->pullHist() ;
  RooPlot* frame2 = Photon_M.frame(Title("Residual Distribution"));
  frame2->addPlotable(hresid,"P");
  RooPlot* frame3 = Photon_M.frame(Title("   "));
  frame3->addPlotable(hpull,"P");

  TCanvas* function=new TCanvas("Function","  ",0,0,900,700);
  gPad->SetLeftMargin(0.15);

  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  pad1->SetBottomMargin(0);
  pad1->Draw();
  pad1->cd();
  frame->GetYaxis()->SetTitleOffset(1.6); 
  frame->Draw();
  TAxis *axis = frame->GetYaxis();
  axis->ChangeLabel(1, -1, 0.);
  axis->Draw();

  TLegend *leg2 = new TLegend(0.65, 0.75, 0.89, 0.89);
  leg2->SetFillColor(kWhite);
  leg2->SetLineColor(kWhite);
  leg2->AddEntry("Data","Data","LP");
  leg2->AddEntry("sum","Signal + background","LP");
  leg2->Draw();

  TLatex *t2 = new TLatex(0.66, 0.71,"pp #sqrt{s}=13 TeV");
  t2->SetTextSize(0.03);
  t2->SetNDC(kTRUE);
  t2->Draw();

  function->cd();
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.01, 1, 0.3);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.2);
  pad2->Draw();
  pad2->cd();
  frame3->GetYaxis()->SetTitle("Pull / 50 [MeV/c^{2}]");
  frame3->GetXaxis()->SetTitle("m(#mu^{+}#mu^{-}) [MeV/c^{2}]");
  frame3->GetYaxis()->SetTitleOffset(1.6); 
  frame3->GetXaxis()->SetTitleSize(0.08);
  frame3->GetXaxis()->SetLabelSize(0.08);
  frame3->GetYaxis()->SetTitleSize(0.08);
  frame3->GetYaxis()->SetLabelSize(0.08);
  frame3->GetYaxis()->SetTitleOffset(0.5);
  frame3->Draw();

  function->SaveAs("Pull_upsilon.png");
 
}
