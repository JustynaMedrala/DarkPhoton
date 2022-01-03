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
  RooRealVar mean_bw("mean_bw","mean", mass_mean, signal_min, signal_max);
  RooRealVar sigma("sigma","sigma", 700, 1e2, 2e4);
  RooRealVar sigma_bw("sigma_bw","sigma", 700, 1e2, 2e4);
  RooRealVar width("width","width", 1000, 500, 1e5);
  RooRealVar alpha("alpha_omega","alpha_omega", 1, 0., 5.);
  RooRealVar N("N_omega", "N_omega", 1.);
  RooCBShape fit("#Z","signal pdf CB",Photon_M,mean,sigma, alpha, N);
  RooBreitWigner fit_bw("#Z_bw","signal pdf omega",Photon_M,mean,sigma_bw);

  RooRealVar lambda("lambda", "slope", -8e-5, -0.1, 0.0);
  RooExponential expo("exponential", "exponential PDF", Photon_M, lambda);

  RooRealVar nsig("nsig","#signal events",2000,1e3,1e5);
  RooRealVar nsig_bw("nsig_bw","#signal events",2000,1e3,1e5);
  RooRealVar nbkg("nbkg","#background events",500,500,1e5);

  Photon_M.setRange("signal", signal_min, signal_max); 
  Photon_M.setRange("all", min_mass, max_mass);

  RooAddPdf sum("sum","sum" , RooArgList(expo, fit, fit_bw), RooArgList(nbkg, nsig, nsig_bw));

  //---------------------------------------------------//
  //FITOWANIE------------------------------------------//
  sum.fitTo(data, RooFit::Extended(), Range("all"));

  RooPlot* frame = Photon_M.frame(RooFit::Title("  "));

  data.plotOn(frame, Name("Data"), DataError(RooAbsData::SumW2));

  sum.plotOn(frame, Name("sum"));
  sum.plotOn(frame, Components(RooArgSet(fit, fit_bw)), LineColor(kRed), LineStyle(9), LineWidth(2), Range(20000, 120000), Name("fit"));
  sum.plotOn(frame, Components("exponential"), LineColor(kGreen), LineStyle(10), LineWidth(4), Name("bkg"));
  frame->GetXaxis()->SetTitle("m(#mu^{+}#mu^{-}) [MeV/c^{2}]");
  frame->GetYaxis()->SetTitle("Candidates/ 1400 [MeV/c^{2}]");
  frame->Draw();

  //---------------------------------------------------//
  //TEST CHI KWADRAT DLA DANYCH------------------------//

  cout<<"Chi kwadrat: "<<frame->chiSquare("sum", "Data", 8)<<endl;

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

  RooHist* hresid = frame->residHist();
  RooHist* hpull = frame->pullHist() ;
  RooPlot* frame2 = Photon_M.frame(Title("Residual Distribution")) ;
  frame2->addPlotable(hresid,"P");
  RooPlot* frame3 = Photon_M.frame(Title("   ")) ;
  frame3->addPlotable(hpull,"P");

  TCanvas* c = new TCanvas("rf109_chi2residpull","rf109_chi2residpull",900,300) ;
  c->Divide(3);
  c->cd(1); gPad->SetLeftMargin(0.15); frame->GetYaxis()->SetTitleOffset(1.6); frame->Draw();
  c->cd(2); gPad->SetLeftMargin(0.15); frame2->GetYaxis()->SetTitleOffset(1.6); frame2->Draw();
  c->cd(3); gPad->SetLeftMargin(0.15); frame3->GetYaxis()->SetTitleOffset(1.6); frame3->Draw();

  TCanvas* function=new TCanvas("Function","  ",0,0,900,700);
  gPad->SetLeftMargin(0.15);

  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  pad1->SetBottomMargin(0);
  pad1->Draw();
  pad1->cd();
  frame->GetYaxis()->SetTitleOffset(1.2); 
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
  frame3->GetYaxis()->SetTitle("Pull / 1400 [MeV/c^{2}]");
  frame3->GetXaxis()->SetTitle("m(#mu^{+}#mu^{-}) [MeV/c^{2}]");
  frame3->GetYaxis()->SetTitleOffset(1.6); 
  frame3->GetXaxis()->SetTitleSize(0.08);
  frame3->GetXaxis()->SetLabelSize(0.08);
  frame3->GetYaxis()->SetTitleSize(0.08);
  frame3->GetYaxis()->SetLabelSize(0.08);
  frame3->GetYaxis()->SetTitleOffset(0.5);
  frame3->Draw();

  function->SaveAs("Pull_Z.png");

}
