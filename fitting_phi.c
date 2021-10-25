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

using namespace RooFit;
using namespace RooStats;

void fitting_phi(){
	
  TFile* input_file = new TFile("result.root", "OPEN");

  TH1* mass_data = dynamic_cast<TH1*>(input_file->Get("ID_1") ); 

  TCanvas* fitted=new TCanvas("Phi","Phi",0,0,900,700);

  Double_t min_mass = 300;
  Double_t max_mass = 1200;
  Double_t mass_mean_phi = 1030;
  Double_t signal_min_phi = 980;
  Double_t signal_max_phi = 1070;

  Double_t mass_mean_omega = 790;
  Double_t signal_min_omega = 750;
  Double_t signal_max_omega = 810;

  Double_t mass_mean_eta = 550;
  Double_t signal_min_eta = 540;
  Double_t signal_max_eta = 575;

  RooRealVar Photon_M("A_M", "A_M", min_mass, max_mass);
  RooDataHist data("data","data",Photon_M,mass_data);

  gPad->SetLeftMargin(0.15);
  //gPad->SetLogy();
  //gPad->SetLogx();

  RooRealVar mean_phi("mean_phi","mean", mass_mean_phi, signal_min_phi, signal_max_phi);
  RooRealVar sigma_phi("sigma_phi","sigma", 5, 1, 15);
  RooRealVar width_phi("width_phi","width", 5, 1, 15);
  RooVoigtian fit_phi("#Phi","signal pdf phi",Photon_M,mean_phi,width_phi,sigma_phi);

  RooRealVar mean_omega("mean_omega","mean", mass_mean_omega, signal_min_omega, signal_max_omega);
  RooRealVar mean_bw_omega("mean_bw_omega","mean", mass_mean_omega, signal_min_omega, signal_max_omega);
  RooRealVar sigma_omega("sigma_omega","sigma", 5, 1, 15);
  RooRealVar sigma_bw_omega("width_omega","width", 5, 1, 20);
  RooRealVar alpha_omega("alpha_omega","alpha_omega", 1, 0., 5.);
  RooRealVar N_omega("N_omega", "N_omega", 1.);
  RooCBShape fit_omega("#Omega","signal pdf omega CB",Photon_M,mean_omega,sigma_omega, alpha_omega, N_omega);
  RooBreitWigner fit_bw_omega("#Omega_bw","signal pdf omega",Photon_M,mean_bw_omega,sigma_bw_omega);

  RooRealVar mean_eta("mean_eta","mean", mass_mean_eta, signal_min_eta, signal_max_eta);
  RooRealVar sigma_eta("sigma_eta","sigma", 5, 0.1, 15);
  RooRealVar width_eta("width_eta","width", 5, 0.1, 15);
  RooVoigtian fit_eta("#Eta","signal pdf eta",Photon_M,mean_eta,width_eta,sigma_eta);

  RooRealVar lambda1("lambda1", "slope1", 0., -1., 0.1);
  RooRealVar lambda2("lambda2", "slope2", 0., -1., 0.1);
  RooExponential expo1("exponential_1", "exponential1 PDF", Photon_M, lambda1);
  RooExponential expo2("exponential_2", "exponential2 PDF", Photon_M, lambda2);

  RooRealVar nsig_phi("nsig_phi","#signal events phi",1000,1e3,1e5);
  RooRealVar nsig_omega("nsig_omega","#signal events omega",1000,1e3,1e6);
  RooRealVar nsig_bw_omega("nsig_omega","#signal events omega BW",1000,1e3,1e6);
  RooRealVar nsig_eta("nsig_eta","#signal events eta",1000,1e2,1e5);
  RooRealVar nbkg1("nbkg1","#background events 1",5000,100,1e6);
  RooRealVar nbkg2("nbkg2","#background events 2",5000,100,1e6);

  Photon_M.setRange("all", min_mass, max_mass);
   
  RooAddPdf sum("sum","sum" , RooArgList(expo1, expo2, fit_phi, fit_bw_omega, fit_omega, fit_eta), RooArgList(nbkg1, nbkg2,nsig_phi,nsig_bw_omega, nsig_omega, nsig_eta)); 
  sum.fitTo(data, RooFit::Extended(), Range("all"));

  RooPlot* frame = Photon_M.frame(RooFit::Title("  "));

  data.plotOn(frame, Name("Data"));

  sum.plotOn(frame, Components("sum"), LineWidth(3), Name("sum"));
  sum.plotOn(frame, Components("#Phi"), LineColor(kRed), LineStyle(9), LineWidth(2), Range(900, 1800), Name("fit_phi"));
  sum.plotOn(frame, Components(RooArgSet(fit_omega, fit_bw_omega)), LineColor(kRed), LineStyle(9), LineWidth(2), Range(650, 900), Name("fit_omega"));
  sum.plotOn(frame, Components("#Eta"), LineColor(kRed), LineStyle(9), LineWidth(2), Range(450, 600), Name("fit_eta"));
  sum.plotOn(frame, Components(RooArgSet(expo1, expo2)), LineColor(kGreen), LineStyle(10), LineWidth(4), Name("bkg"));
  frame->GetXaxis()->SetTitle("m(#mu^{+}#mu^{-}) [MeV/c^{2}]");
  frame->GetYaxis()->SetTitle("Candidates/5 [MeV/c^{2}]");
  frame->Draw();

  TLegend *leg1 = new TLegend(0.18, 0.65, 0.43, 0.89);
  leg1->SetFillColor(kWhite);
  leg1->SetLineColor(kWhite);
  leg1->AddEntry("Data","Data","LP");
  leg1->AddEntry("sum","Signal + background","LP");
  leg1->AddEntry("bkg","Background only", "LP");
  leg1->AddEntry("fit_phi","Signal only", "LP");
  leg1->Draw();

  TLatex *t = new TLatex(0.19, 0.62,"pp #sqrt{s}=13 TeV");
  t->SetTextSize(0.03);
  t->SetNDC(kTRUE);
  t->Draw();

  fitted->SaveAs("Phi_Omega_Eta.png");

  //TCanvas* no_bkg=new TCanvas("Without bkg","  ",0,0,900,700);
  TH1::SetDefaultSumw2();
  TH1 *hdata = data.createHistogram("hdata",Photon_M,Binning((max_mass-min_mass)/5,min_mass, max_mass));
  /*TH1 *hfit = expo1.createHistogram("hfit",Photon_M,Binning((max_mass-min_mass)/5,min_mass, max_mass));
  TH1 *hbkg2 = expo2.createHistogram("hbkg2",Photon_M,Binning((max_mass-min_mass)/5,min_mass, max_mass));
  TH1 *hphi = fit_phi.createHistogram("hphi",Photon_M,Binning((max_mass-min_mass)/5,min_mass, max_mass));
  TH1 *homega = fit_omega.createHistogram("homega",Photon_M,Binning((max_mass-min_mass)/5,min_mass, max_mass));
  TH1 *homega_bw = fit_bw_omega.createHistogram("homega_bw",Photon_M,Binning((max_mass-min_mass)/5,min_mass, max_mass));
  TH1 *heta = fit_eta.createHistogram("heta",Photon_M,Binning((max_mass-min_mass)/5,min_mass, max_mass));

  gPad->SetLeftMargin(0.15);
  gStyle->SetOptStat(0);
  //gPad->SetLogy();
  //gPad->SetLogx();
  hfit->Scale(nbkg1.getVal());
  hbkg2->Sumw2();  
  hbkg2->Scale(nbkg2.getVal());
  hphi->Sumw2(); 
  hphi->Scale(nsig_phi.getVal());
  homega->Sumw2(); 
  homega->Scale(nsig_omega.getVal());
  homega_bw->Sumw2(); 
  homega_bw->Scale(nsig_bw_omega.getVal());
  heta->Sumw2(); 
  heta->Scale(nsig_eta.getVal());

  hfit->Add(hbkg2);
  hfit->Add(hphi);
  hfit->Add(homega);
  hfit->Add(homega_bw);
  hfit->Add(heta);
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
  hdata->GetYaxis()->SetTitle("Data events/7 [MeV/c^{2}]");
  hdata->Draw("EP");
  hfit->Draw("same l");
  TAxis *axis = hdata->GetYaxis();
  axis->ChangeLabel(1, -1, 0.);
  axis->Draw();

  TLegend *leg2 = new TLegend(0.15, 0.65, 0.40, 0.89);
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
  no_bkg->SaveAs("Phi_without_bkg.png");*/

  TCanvas* function=new TCanvas("Function","  ",0,0,900,700);

  gPad->SetLeftMargin(0.15);
  gStyle->SetOptStat(0);

  RooArgList pars_bkg1(*expo1.getParameters(RooArgSet(Photon_M)));
  RooArgSet prodSet_bkg1(expo1); 
  prodSet_bkg1.add(nbkg1);
  RooProduct unNormPdf_bkg1("fitted_Function_bkg_1", "fitted_Function_bkg_1", prodSet_bkg1);
  TF1 * f1 = unNormPdf_bkg1.asTF(RooArgList(Photon_M), pars_bkg1, RooArgList(Photon_M));
  f1->SetName("fit_bkg1");

  RooArgList pars_bkg2(*expo2.getParameters(RooArgSet(Photon_M)));
  RooArgSet prodSet_bkg2(expo2); 
  prodSet_bkg2.add(nbkg2);
  RooProduct unNormPdf_bkg2("fitted_Function_bkg_2", "fitted_Function_bkg_2", prodSet_bkg2);
  TF1 * f2 = unNormPdf_bkg2.asTF(RooArgList(Photon_M), pars_bkg2, RooArgList(Photon_M));
  f2->SetName("fit_bkg2");

  RooArgList pars_eta(*fit_eta.getParameters(RooArgSet(Photon_M)));
  RooArgSet prodSet_eta(fit_eta); 
  prodSet_eta.add(nsig_eta);
  RooProduct unNormPdf_eta("fitted_Function_eta", "fitted_Function_eta", prodSet_eta);
  TF1 * f3 = unNormPdf_eta.asTF(RooArgList(Photon_M), pars_eta, RooArgList(Photon_M));
  f3->SetName("fit_eta");

  RooArgList pars_omega(*fit_omega.getParameters(RooArgSet(Photon_M)));
  RooArgSet prodSet_omega(fit_omega); 
  prodSet_omega.add(nsig_omega);
  RooProduct unNormPdf_omega("fitted_Function_omega", "ffitted_Function_omega", prodSet_omega);
  TF1 * f4 = unNormPdf_omega.asTF(RooArgList(Photon_M), pars_omega, RooArgList(Photon_M));
  f4->SetName("fit_omega"); 

  RooArgList pars_bw_omega(*fit_bw_omega.getParameters(RooArgSet(Photon_M)));
  RooArgSet prodSet_bw_omega(fit_bw_omega); 
  prodSet_bw_omega.add(nsig_bw_omega);
  RooProduct unNormPdf_bw_omega("fitted_Function_bw_omega", "fitted_Function_bw_omega", prodSet_bw_omega);
  TF1 * f5 = unNormPdf_bw_omega.asTF(RooArgList(Photon_M), pars_bw_omega, RooArgList(Photon_M));
  f5->SetName("fit_bw_omega");

  RooArgList pars_phi(*fit_phi.getParameters(RooArgSet(Photon_M)));
  RooArgSet prodSet_phi(fit_phi);
  prodSet_phi.add(nsig_phi);
  RooProduct unNormPdf_phi("fitted_Function_phi", "fitted_Function_phi", prodSet_phi);
  TF1 * f6 = unNormPdf_phi.asTF(RooArgList(Photon_M), pars_phi, RooArgList(Photon_M));
  f6->SetName("fit_phi");

  RooArgList pars(*sum.getParameters(RooArgSet(Photon_M)));
  RooArgSet prodSet(sum);
  RooProduct unNormPdf("fitted_Function", "fitted_Function", prodSet);
  TF1 * f7 = unNormPdf.asTF(RooArgList(Photon_M), pars, RooArgList(Photon_M));
  f7->SetName("fit");

  for (int i = 0; i < f7->GetNpar(); ++i)
      cout<<f7->GetParName(i)<<endl;

  TF1NormSum *fnorm = new TF1NormSum(f1,f2,f3,f1->Integral(min_mass, max_mass),f2->Integral(min_mass, max_mass),f3->Integral(min_mass, max_mass));
  TF1 *f_sum = new TF1("fsum", *fnorm, min_mass, max_mass, fnorm->GetNpar());
  f_sum->SetParameters(fnorm->GetParameters().data());
  
  TF1NormSum *fnorm2 = new TF1NormSum(f_sum,f4,f5,f_sum->Integral(min_mass, max_mass),f4->Integral(min_mass, max_mass),f5->Integral(min_mass, max_mass));
  TF1 *f_sum2 = new TF1("fsum2", *fnorm2, min_mass, max_mass, fnorm2->GetNpar());
  f_sum2->SetParameters(fnorm2->GetParameters().data());

  f_sum2->SetParName(4,"Nomega_bw");
  f_sum2->SetParName(3,"Nomega");
  f_sum2->SetParName(2,"Neta");
  f_sum2->SetParName(1,"Nexpo2");
  f_sum2->SetParName(0,"Nexpo1");

  for (int i = 0; i < f_sum2->GetNpar(); ++i)
      f_sum2->SetParName(i,fnorm2->GetParName(i));

  TF1NormSum *fnorm3 = new TF1NormSum(f_sum2,f6,f_sum2->Integral(min_mass, max_mass),f6->Integral(min_mass, max_mass));
  TF1 *f_all = new TF1("fall", *fnorm3, min_mass, max_mass, fnorm3->GetNpar());
  f_all->SetParameters(fnorm3->GetParameters().data());

  f_all->SetParName(7,"Nphi");
  f_all->SetParName(6,"Nsum2");
  f_all->SetParName(5,"Nomega_bw");
  f_all->SetParName(4,"Nomega");
  f_all->SetParName(3,"Nsum");
  f_all->SetParName(2,"Neta");
  f_all->SetParName(1,"Nexpo2");
  f_all->SetParName(0,"Nexpo1");

  for (int i = 8; i < f_all->GetNpar(); ++i)
      f_all->SetParName(i,fnorm3->GetParName(i));

  TH1::SetDefaultSumw2();

  for (int i=0;i<(max_mass-min_mass)/5;i++)
  {
    float val=hdata->GetBinContent(i); 
    hdata->SetBinError(i, sqrt(val));
  }

  TH1F *hdif = new TH1F("hdif", " ", (max_mass-min_mass)/5, min_mass, max_mass);
  TH1F *hdiv = new TH1F("hdiv", " ", (max_mass-min_mass)/5, min_mass, max_mass);
  TH1F *hdat = (TH1F*) hdata->Clone();

  hdat->Sumw2();
  hdif->Sumw2();

  f_all->SetParLimits(0, 1e5, 1e7);
  f_all->SetParLimits(1, 1e5, 1e7);
  f_all->SetParLimits(2, 1.5e5, 1e6);
  //f_all->SetParLimits(3, 1e5, 1e7);
  f_all->SetParLimits(4, 1e5, 1e7);
  f_all->SetParLimits(5, 1e5, 1e7);
  //f_all->SetParLimits(6, 1e5, 1e7);
  f_all->SetParLimits(7, 1e5, 1e6);
  f_all->SetParLimits(8, lambda1.getVal(), lambda1.getVal());
  f_all->SetParLimits(9, lambda2.getVal(), lambda2.getVal());
  f_all->SetParLimits(10, mean_eta.getVal(), mean_eta.getVal());
  f_all->SetParLimits(11, sigma_eta.getVal(), sigma_eta.getVal());
  f_all->SetParLimits(12, width_eta.getVal(), width_eta.getVal());
  f_all->SetParLimits(13, N_omega.getVal(), N_omega.getVal());
  f_all->SetParLimits(14, alpha_omega.getVal(), alpha_omega.getVal());
  f_all->SetParLimits(15, mean_omega.getVal(), mean_omega.getVal());
  f_all->SetParLimits(16, sigma_omega.getVal(), sigma_omega.getVal());
  f_all->SetParLimits(17, mean_bw_omega.getVal(), mean_bw_omega.getVal());
  f_all->SetParLimits(18, sigma_bw_omega.getVal(), sigma_bw_omega.getVal());
  f_all->SetParLimits(19, mean_phi.getVal(), mean_phi.getVal());
  f_all->SetParLimits(20, sigma_phi.getVal(), sigma_phi.getVal());
  f_all->SetParLimits(21, width_phi.getVal(), width_phi.getVal());

  f7->SetParLimits(0, N_omega.getVal(), N_omega.getVal());
  f7->SetParLimits(1, alpha_omega.getVal(), alpha_omega.getVal());
  f7->SetParLimits(2, lambda1.getVal(), lambda1.getVal());
  f7->SetParLimits(3, lambda2.getVal(), lambda2.getVal());
  f7->SetParLimits(4, mean_bw_omega.getVal(), mean_bw_omega.getVal());
  f7->SetParLimits(5, mean_eta.getVal(), mean_eta.getVal());
  f7->SetParLimits(6, mean_omega.getVal(), mean_omega.getVal());
  f7->SetParLimits(7, mean_phi.getVal(), mean_phi.getVal());
  f7->SetParLimits(8, nbkg1.getVal(), nbkg1.getVal());
  f7->SetParLimits(9, nbkg2.getVal(), nbkg2.getVal());
  f7->SetParLimits(10, nsig_eta.getVal(), nsig_eta.getVal());
  f7->SetParLimits(11, nsig_omega.getVal(), nsig_omega.getVal());
  f7->SetParLimits(12, nsig_phi.getVal(), nsig_phi.getVal());
  f7->SetParLimits(13, sigma_eta.getVal(), sigma_eta.getVal());
  f7->SetParLimits(14, sigma_omega.getVal(), sigma_omega.getVal());
  f7->SetParLimits(15, sigma_phi.getVal(), sigma_phi.getVal());
  f7->SetParLimits(16, width_eta.getVal(), width_eta.getVal());
  f7->SetParLimits(17, sigma_bw_omega.getVal(), sigma_bw_omega.getVal());
  f7->SetParLimits(18, width_phi.getVal(), width_phi.getVal());

  hdata->Fit(f_all, "R");

  //f_all->Draw();

  //hdata->Draw("same");
  //f1->DrawCopy("same");
  //f1->DrawCopy();
  //f3->DrawCopy("same");
  //f4->DrawCopy("same");

  for (int i=0;i<(max_mass-min_mass)/5;i++)
  {
    float x=hdat->GetBinCenter(i);
    float y=f_all->Eval(x);
    hdif->Fill(x, y);
  }

  hdat->Add(hdif,-1);

 for (int i=0;i<(max_mass-min_mass)/5;i++)
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
  hdata->GetYaxis()->SetTitle("Data events/5 [MeV/c^{2}]");
  hdata->Draw();
  TAxis *axis = hdata->GetYaxis();
  axis->ChangeLabel(1, -1, 0.);
  axis->Draw();

  TLegend *leg2 = new TLegend(0.15, 0.75, 0.35, 0.89);
  leg2->SetFillColor(kWhite);
  leg2->SetLineColor(kWhite);
  leg2->AddEntry(hdata,"Data","P");
  leg2->AddEntry(f_all,"Signal + background","L");
  leg2->Draw();

  TLatex *t2 = new TLatex(0.17, 0.71,"pp #sqrt{s}=13 TeV");
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
  hdiv->GetYaxis()->SetTitle("Pull");
  hdiv->GetXaxis()->SetTitle("m(#mu^{+}#mu^{-}) [MeV/c^{2}]");
  hdiv->SetTitle("  ");
  hdiv->Draw();

  function->SaveAs("func_phi.png");

  //input_file->Close();


}
