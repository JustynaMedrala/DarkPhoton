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
  Double_t max_mass = 11976;    //1 - 9800  //2 - 10180 //3 - 12000
  Double_t mass_mean_u1 = 9450;    //1 - 9450  //2 - 9900  //3 - 10250
  Double_t signal_min_u1 = 9340;   //1 - 9340  //2 - 9900  //3 - 10150
  Double_t signal_max_u1 = 9600;  //1 - 9600  //2 - 10100 //3 - 10450
  Double_t mass_mean_u2 = 9900;
  Double_t signal_min_u2 = 9900;
  Double_t signal_max_u2 = 10100;
  Double_t mass_mean_u3 = 10250;
  Double_t signal_min_u3 = 10150;
  Double_t signal_max_u3 = 10450;  
  Double_t mass_mean_z = 90000;
  Double_t signal_min_z = 82000;
  Double_t signal_max_z = 95000;

  RooRealVar Photon_M("A_M", "A_M",mass_mean_u1, min_mass, max_mass);
  RooDataHist data("data","data",Photon_M,mass_data);

  RooRealVar mean_u1("mean_u1","mean", mass_mean_u1, signal_min_u1, signal_max_u1);
  RooRealVar mean_u2("mean_u2","mean", mass_mean_u2, signal_min_u2, signal_max_u2);
  RooRealVar mean_u3("mean_u3","mean", mass_mean_u3, signal_min_u3, signal_max_u3);
  RooRealVar mean_z("mean_z","mean", mass_mean_z, signal_min_z, signal_max_z);
  RooRealVar sigma_u1("sigma_u1","sigma", 30, 10, 50);
  RooRealVar sigma_u2("sigma_u2","sigma", 40, 10, 50);
  RooRealVar sigma_u3("sigma_u3","sigma", 40, 10, 50);
  RooRealVar sigma_z("sigma_z","sigma", 70, 0.0001, 2e4);
  RooRealVar width_u1("width_u1","width", 20, 5, 50);
  RooRealVar width_u2("width_u2","width", 10, 5, 50);
  RooRealVar width_u3("width_u3","width", 40, 5, 60);
  RooRealVar width_z("width_z","width", 50, 10, 1e5);
  RooVoigtian fit_u1("#Upsilon_1","signal pdf",Photon_M,mean_u1,width_u1,sigma_u1);
  RooVoigtian fit_u2("#Upsilon_2","signal pdf",Photon_M,mean_u2,width_u2,sigma_u2);
  RooVoigtian fit_u3("#Upsilon_3","signal pdf",Photon_M,mean_u3,width_u3,sigma_u3);
  RooVoigtian fit_z("Z0","signal pdf",Photon_M,mean_z,width_z,sigma_z);

  RooRealVar nsig_u1("nsig_u1","#signal events",63000,1e2,1e6);
  RooRealVar nsig_u2("nsig_u2","#signal events",15000,1e2,1e5);
  RooRealVar nsig_u3("nsig_u3","#signal events",750,100,1e5);
  RooRealVar nsig_z("nsig_z","#signal events",20000,1e3,1e9);
  RooRealVar nbkg("nbkg","#background events",93000,1e3,1e8);

  RooRealVar lambda("lambda", "slope", 0., -1., 0.1);
  RooExponential expo("exponential", "exponential PDF", Photon_M, lambda);

  Photon_M.setRange("signal_u1", signal_min_u1, signal_max_u1);
  Photon_M.setRange("signal_u2", signal_min_u2, signal_max_u2);
  Photon_M.setRange("signal_u3", signal_min_u3, signal_max_u3);
  Photon_M.setRange("signal_z", signal_min_z, signal_max_z);
  Photon_M.setRange("all", min_mass, max_mass);

  RooAddPdf sum("sum","sum" , RooArgList(expo, fit_u1, fit_u2, fit_u3), RooArgList(nbkg, nsig_u1, nsig_u2, nsig_u3)); 
  sum.fitTo(data, RooFit::Extended(), Range("all"));

  TCanvas* fitted=new TCanvas("Upsilon","Upsilon",0,0,900,700);
  gPad->SetLeftMargin(0.15);
  //gPad->SetLogy();
  //gPad->SetLogx();

  RooPlot* frame = Photon_M.frame(RooFit::Title("  "));

  data.plotOn(frame, Name("Data"));

  sum.plotOn(frame, Components("sum"), LineWidth(3), Name("sum"));
  sum.plotOn(frame, Components("#Upsilon_1"), LineColor(kRed), LineStyle(9), LineWidth(2), Range(9000, 9800), Name("fit_u1"));
  sum.plotOn(frame, Components("#Upsilon_2"), LineColor(kRed), LineStyle(9), LineWidth(2), Range(9800, 10200), Name("fit_u2"));
  sum.plotOn(frame, Components("#Upsilon_3"), LineColor(kRed), LineStyle(9), LineWidth(2), Range(10100, 10800), Name("fit_u3"));
  sum.plotOn(frame, Components("Z0"), LineColor(kRed), LineStyle(9), LineWidth(2), Range(30000, 120000), Name("fit_z"));
  sum.plotOn(frame, Components("exponential"), LineColor(kGreen), LineStyle(10), LineWidth(4), Name("bkg"));
  frame->GetXaxis()->SetTitle("m(#mu^{+}#mu^{-}) [MeV/c^{2}]");
  frame->GetYaxis()->SetTitle("Candidates/56 [MeV/c^{2}]");
  frame->Draw();

  //mass_data->Draw("hist");
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

  fitted->SaveAs("upsilon.png");

  //TCanvas* no_bkg=new TCanvas("Without bkg","  ",0,0,900,700);
  TH1::SetDefaultSumw2();
  TH1 *hdata = data.createHistogram("hdata",Photon_M,Binning((max_mass-min_mass)/56,min_mass, max_mass));
  hdata->Sumw2();
  /*TH1 *hfit = expo.createHistogram("hfit",Photon_M,Binning((max_mass-min_mass)/56,min_mass, max_mass));
  TH1 *hu1 = fit_u1.createHistogram("hu1",Photon_M,Binning((max_mass-min_mass)/56,min_mass, max_mass));
  TH1 *hu2 = fit_u2.createHistogram("hu2",Photon_M,Binning((max_mass-min_mass)/56,min_mass, max_mass));
  TH1 *hu3 = fit_u3.createHistogram("hu3",Photon_M,Binning((max_mass-min_mass)/56,min_mass, max_mass));

  gPad->SetLeftMargin(0.15);
  gStyle->SetOptStat(0);
  //gPad->SetLogy();
  //gPad->SetLogx();
  hfit->Scale(nbkg.getVal());
  hu1->Sumw2(); 
  hu1->Scale(nsig_u1.getVal());
  hu2->Sumw2(); 
  hu2->Scale(nsig_u2.getVal());
  hu3->Sumw2(); 
  hu3->Scale(nsig_u3.getVal());
  hfit->Add(hu1);
  hfit->Add(hu2);
  hfit->Add(hu3);
  hdata->Sumw2();
  hfit->Sumw2();
  hfit->SetLineColor(kRed);
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  pad1->SetBottomMargin(0);
  pad1->Draw();
  pad1->cd();
  hdata->SetTitle("   ");
  hdata->SetMarkerStyle(20);
  hdata->SetMarkerSize(0.8);
  hfit->SetLineWidth(2);
  hdata->GetXaxis()->SetTitle("m(#mu^{+}#mu^{-}) [MeV/c^{2}]");
  hdata->GetYaxis()->SetTitle("Data events/56 [MeV/c^{2}]");
  hdata->Draw("EP");
  hfit->Draw("same l");
  TAxis *axis = hdata->GetYaxis();
  axis->ChangeLabel(1, -1, 0.);
  axis->Draw();

  TLegend *leg2 = new TLegend(0.65, 0.65, 0.89, 0.89);
  leg2->SetFillColor(kWhite);
  leg2->SetLineColor(kWhite);
  leg2->AddEntry(hdata,"Data","P");
  leg2->AddEntry(hfit,"Signal + background","L");
  leg2->Draw();

  no_bkg->cd();
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.01, 1, 0.3);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.2);
  pad2->Draw();
  pad2->cd();
  TH1F *hdiff = (TH1F*)hdata->Clone("hdiff");
  hdiff->Sumw2();
  hdiff->Add(hfit, -1);
  //hdiff->GetYaxis()->SetTitleSize(20);
  hdiff->GetXaxis()->SetTitleSize(0.08);
  hdiff->GetXaxis()->SetLabelSize(0.08);
  hdiff->GetYaxis()->SetTitleSize(0.08);
  hdiff->GetYaxis()->SetLabelSize(0.08);
  hdiff->GetYaxis()->SetTitleOffset(0.5);
  hdiff->GetYaxis()->SetTitle("Difference");
  hdiff->Draw();
  //auto rp1 = new TRatioPlot(hdata, hfit, "diff");
  //rp1->SetH1DrawOpt("PLC");
  //rp1->Draw();
  no_bkg->SaveAs("Ups_without_bkg.png");*/

  TCanvas* function=new TCanvas("Function","  ",0,0,900,700);

  RooArgList pars_bkg(*expo.getParameters(RooArgSet(Photon_M)));
  RooArgSet prodSet_bkg(expo); 
  prodSet_bkg.add(nbkg);
  RooProduct unNormPdf_bkg("fitted_Function_bkg", "fitted_Function_bkg", prodSet_bkg);
  TF1 * f1 = unNormPdf_bkg.asTF(RooArgList(Photon_M), pars_bkg, RooArgList(Photon_M));
  f1->SetName("fit_bkg");

  RooArgList pars_u1(*fit_u1.getParameters(RooArgSet(Photon_M)));
  RooArgSet prodSet_u1(fit_u1); 
  prodSet_u1.add(nsig_u1);
  RooProduct unNormPdf_u1("fitted_Function_u1", "fitted_Function_u1", prodSet_u1);
  TF1 * f2 = unNormPdf_u1.asTF(RooArgList(Photon_M), pars_u1, RooArgList(Photon_M));
  f2->SetName("fit_u1");

  RooArgList pars_u2(*fit_u2.getParameters(RooArgSet(Photon_M)));
  RooArgSet prodSet_u2(fit_u2); 
  prodSet_u2.add(nsig_u2);
  RooProduct unNormPdf_u2("fitted_Function_u2", "ffitted_Function_u2", prodSet_u2);
  TF1 * f3 = unNormPdf_u2.asTF(RooArgList(Photon_M), pars_u2, RooArgList(Photon_M));
  f3->SetName("fit_u2"); 

  RooArgList pars_u3(*fit_u3.getParameters(RooArgSet(Photon_M)));
  RooArgSet prodSet_u3(fit_u3); 
  prodSet_u3.add(nsig_u3);
  RooProduct unNormPdf_u3("fitted_Function_u3", "fitted_Function_u3", prodSet_u3);
  TF1 * f4 = unNormPdf_u3.asTF(RooArgList(Photon_M), pars_u3, RooArgList(Photon_M));
  f4->SetName("fit_u3");

  RooArgList pars(*sum.getParameters(RooArgSet(Photon_M)));
  RooArgSet prodSet(sum);
  RooProduct unNormPdf("fitted_Function", "fitted_Function", prodSet);
  TF1 * f5 = unNormPdf.asTF(RooArgList(Photon_M), pars, RooArgList(Photon_M));
  f5->SetName("sum");

  TF1NormSum *fnorm = new TF1NormSum(f1, f2, f3, f1->Integral(min_mass, max_mass), f2->Integral(min_mass, max_mass), f3->Integral(min_mass, max_mass));
  TF1 *f_sum = new TF1("fsum", *fnorm, min_mass, max_mass, fnorm->GetNpar());
  f_sum->SetParameters(fnorm->GetParameters().data());
  TF1NormSum *fnorm2 = new TF1NormSum(f_sum,f4,f_sum->Integral(min_mass, max_mass), f4->Integral(min_mass, max_mass));
  TF1 *f_all = new TF1("fall", *fnorm2, min_mass, max_mass, fnorm2->GetNpar());
  f_all->SetParameters(fnorm2->GetParameters().data());

  TH1::SetDefaultSumw2();

 for (int i=0;i<(max_mass-min_mass)/56;i++)
 {
   float val=hdata->GetBinContent(i); 
   hdata->SetBinError(i, sqrt(val));
 }

  TH1F *hdif = new TH1F("hdif", " ", (max_mass-min_mass)/56, min_mass, max_mass);
  TH1F *hdiv = new TH1F("hdiv", " ", (max_mass-min_mass)/56, min_mass, max_mass);
  TH1F *hdat = (TH1F*) hdata->Clone();

  hdat->Sumw2();
  hdif->Sumw2();

  hdata->Fit(f_all, "R");

  //f_all->Draw();

  //hdata->Draw("same");
  //f1->DrawCopy("same");
  //f1->DrawCopy();
  //f3->DrawCopy("same");
  //f4->DrawCopy("same");

 for (int i=0;i<(max_mass-min_mass)/56;i++)
 {
   float x=hdat->GetBinCenter(i);
   float y=f_all->Eval(x);
   hdif->Fill(x, y);
 }

  hdat->Add(hdif,-1);

 for (int i=0;i<(max_mass-min_mass)/56;i++)
 {
   float x=hdat->GetBinCenter(i);
   float y=hdat->GetBinContent(i);
   hdiv->Fill(x, sqrt(abs(y)));
 }

  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  pad1->SetBottomMargin(0);
  pad1->Draw();
  pad1->cd();
  hdata->SetTitle("   ");
  hdata->SetMarkerStyle(20);
  hdata->SetMarkerSize(0.8);
  hdata->GetXaxis()->SetTitle("m(#mu^{+}#mu^{-}) [MeV/c^{2}]");
  hdata->GetYaxis()->SetTitle("Data events/56 [MeV/c^{2}]");
  hdata->GetYaxis()->SetRangeUser(0, 200e3);
  hdata->Draw();
  TAxis *axis = hdata->GetYaxis();
  axis->ChangeLabel(1, -1, 0.);
  axis->Draw();

  TLegend *leg2 = new TLegend(0.65, 0.75, 0.89, 0.89);
  leg2->SetFillColor(kWhite);
  leg2->SetLineColor(kWhite);
  leg2->AddEntry(hdata,"Data","P");
  leg2->AddEntry(f_all,"Signal + background","L");
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
  hdiv->GetXaxis()->SetTitleSize(0.08);
  hdiv->GetXaxis()->SetLabelSize(0.08);
  hdiv->GetYaxis()->SetTitleSize(0.08);
  hdiv->GetYaxis()->SetLabelSize(0.08);
  hdiv->GetYaxis()->SetTitleOffset(0.5);
  hdiv->GetXaxis()->SetTitle("m(#mu^{+}#mu^{-}) [MeV/c^{2}]");
  hdiv->GetYaxis()->SetTitle("Pull");
  hdiv->SetTitle("  ");
  hdiv->Draw();

  function->SaveAs("func_upsilon.png");

  //input_file->Close();

}
