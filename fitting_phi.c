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

  Double_t mass_mean_A = 410;
  Double_t signal_min_A = 400;
  Double_t signal_max_A = 420;

  RooRealVar Photon_M("A_M", "A_M", min_mass, max_mass);
  RooRealVar eta_M("A_MM", "A_MM", min_mass, 600);
  RooDataHist data("data","data",Photon_M,mass_data);
  RooDataHist data_eta("data_eta","data_eta",eta_M,mass_data);

  gPad->SetLeftMargin(0.15);

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
  RooBreitWigner fit_bw_omega("#Omega_bw","signal pdf omega",Photon_M,mean_omega,sigma_bw_omega);

  RooRealVar mean_eta("mean_eta","mean", mass_mean_eta, signal_min_eta, signal_max_eta);
  RooRealVar sigma_eta("sigma_eta","sigma", 5, 0.1, 15);
  RooRealVar width_eta("width_eta","width", 5, 0.1, 15);
  RooVoigtian fit_eta("#Eta","signal pdf eta",Photon_M,mean_eta,width_eta,sigma_eta);

  RooRealVar mean_A("mean_A","mean", mass_mean_A, signal_min_A, signal_max_A);
  RooRealVar sigma_A("sigma_A","sigma", 5, 0.1, 10);
  RooRealVar width_A("width_A","width", 5, 0.1, 10);
  RooVoigtian fit_A("#A","signal pdf A",Photon_M,mean_A,width_A,sigma_A);

  RooRealVar lambda1("lambda1", "slope1", 0., -1., 0.1);
  RooRealVar lambda2("lambda2", "slope2", 0., -1., 0.1);
  RooExponential expo1("exponential_1", "exponential1 PDF", Photon_M, lambda1);
  RooExponential expo2("exponential_2", "exponential2 PDF", Photon_M, lambda2);

  RooRealVar a("a", "a", -1., -1., 0.);
  RooRealVar b("b", "b", 1., 0., 100.);
  RooRealVar d("d", "d", 1., -100., 100.);
  RooPolynomial pol("pol", "pol", Photon_M, RooArgList(a, b));

  RooRealVar nsig_phi("nsig_phi","#signal events phi",1000,1e3,1e5);
  RooRealVar nsig_omega("nsig_omega","#signal events omega",1000,1e3,1e6);
  RooRealVar nsig_bw_omega("nsig_bw_omega","#signal events omega BW",1000,1e3,1e6);
  RooRealVar nsig_eta("nsig_eta","#signal events eta",100,1e2,1e5);
  RooRealVar nsig_A("nsig_A","#signal events A",100,1e2,1e5);
  RooRealVar nbkg1("nbkg1","#background events 1",5000,100,1e6);
  RooRealVar nbkg2("nbkg2","#background events 2",5000,100,1e6);

  Photon_M.setRange("all", min_mass, max_mass);
   
  RooAddPdf sum("sum","sum" , RooArgList(expo1, expo2, fit_phi, fit_bw_omega, fit_omega, fit_eta, fit_A), RooArgList(nbkg1, nbkg2,nsig_phi,nsig_bw_omega, nsig_omega, nsig_eta, nsig_A)); 

  //---------------------------------------------------//
  //FITOWANIE------------------------------------------//
  sum.fitTo(data, RooFit::Extended(), Range("all"));

  RooPlot* frame = Photon_M.frame(RooFit::Title("  "));

  data.plotOn(frame, Name("Data"), DataError(RooAbsData::SumW2));

  sum.plotOn(frame, Components("sum"), LineWidth(3), Name("sum"));
  sum.plotOn(frame, Components("#Phi"), LineColor(kRed), LineStyle(9), LineWidth(2), Range(900, 1800), Name("fit_phi"));
  sum.plotOn(frame, Components(RooArgSet(fit_omega, fit_bw_omega)), LineColor(kRed), LineStyle(9), LineWidth(2), Range(650, 900), Name("fit_omega"));
  sum.plotOn(frame, Components("#Eta"), LineColor(kRed), LineStyle(9), LineWidth(2), Range(450, 600), Name("fit_eta"));
  sum.plotOn(frame, Components("#A"), LineColor(kRed), LineStyle(9), LineWidth(2), Range(350, 450), Name("fit_A"));
  sum.plotOn(frame, Components(RooArgSet(expo1, expo2)), LineColor(kGreen), LineStyle(10), LineWidth(4), Name("bkg"));
  frame->GetXaxis()->SetTitle("m(#mu^{+}#mu^{-}) [MeV/c^{2}]");
  frame->GetYaxis()->SetTitle("Candidates/5 [MeV/c^{2}]");
  frame->Draw();

  //---------------------------------------------------//
  //TEST CHI KWADRAT DLA DANYCH------------------------//

  cout<<"Chi kwadrat/ stopień swobody: "<<frame->chiSquare("sum", "Data", 18)<<endl;

  RooChi2Var chi2 ("chi2", "chi2", sum, data);
  Double_t chi2_val = chi2.getVal();
  cout<<"Chi kwadrat: "<<chi2_val<<" stopnie swobody: "<<(max_mass-min_mass)/mass_data->GetBinWidth(1)-18<<endl;

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

  TLegend *leg2 = new TLegend(0.18, 0.65, 0.43, 0.89);
  leg2->SetFillColor(kWhite);
  leg2->SetLineColor(kWhite);
  leg2->AddEntry("Data","Data","LP");
  leg2->AddEntry("sum","Signal + background","LP");
  leg2->Draw();

  TLatex *t2 = new TLatex(0.19, 0.62,"pp #sqrt{s}=13 TeV");
  t2->SetTextSize(0.03);
  t2->SetNDC(kTRUE);
  t2->Draw();

  function->cd();
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.01, 1, 0.3);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.2);
  pad2->Draw();
  pad2->cd();
  frame3->GetYaxis()->SetTitle("Pull / 5 [MeV/c^{2}]");
  frame3->GetXaxis()->SetTitle("m(#mu^{+}#mu^{-}) [MeV/c^{2}]");
  frame3->GetYaxis()->SetTitleOffset(1.6); 
  frame3->GetXaxis()->SetTitleSize(0.08);
  frame3->GetXaxis()->SetLabelSize(0.08);
  frame3->GetYaxis()->SetTitleSize(0.08);
  frame3->GetYaxis()->SetLabelSize(0.08);
  frame3->GetYaxis()->SetTitleOffset(0.5);
  frame3->GetYaxis()->SetRangeUser(-8.5, 7);
  frame3->Draw();

  function->SaveAs("Pull_phi.png");

  TCanvas* c1 = new TCanvas("eta","eta",900,700);
  gPad->SetLeftMargin(0.12);
  eta_M.setRange("eta", min_mass, 600);
  RooRealVar nsig_eta_1("nsig_eta_1","#signal events eta",100,1e2,1e5);
  RooRealVar nsig_A_1("nsig_A_1","#signal events A",100,1e2,1e5);
  RooRealVar nbkg1_1("nbkg1_1","#background events 1",5000,100,1e6);
  RooRealVar nbkg2_1("nbkg2_1","#background events 2",5000,100,1e6);

  RooRealVar mean_eta_1("mean_eta_1","mean", mean_eta.getVal());
  RooRealVar sigma_eta_1("sigma_eta_1","sigma", sigma_eta.getVal());
  RooRealVar width_eta_1("width_eta_1","width", width_eta.getVal());
  RooVoigtian fit_eta_1("#Eta_1","signal pdf eta",eta_M,mean_eta_1,width_eta_1,sigma_eta_1);

  RooRealVar mean_A_1("mean_A_1","mean", mean_A.getVal());
  RooRealVar sigma_A_1("sigma_A_1","sigma", sigma_A.getVal());
  RooRealVar width_A_1("width_A_1","width", width_A.getVal());
  RooVoigtian fit_A_1("#A_1","signal pdf A",eta_M,mean_A_1,width_A_1,sigma_A_1);

  RooRealVar lambda1_1("lambda1_1", "slope1", lambda1.getVal());
  RooRealVar lambda2_1("lambda2_1", "slope2", lambda2.getVal());
  RooExponential expo1_1("exponential_1", "exponential1 PDF", eta_M, lambda1_1);
  RooExponential expo2_1("exponential_2", "exponential2 PDF", eta_M, lambda2_1);

  RooRealVar a_1("a1", "a1", a.getVal());
  RooRealVar b_1("b1", "b1", b.getVal());
  RooRealVar d_1("c1", "c1", d.getVal());
  RooPolynomial pol_1("pol1", "pol1", eta_M, RooArgList(a_1, b_1, d_1));


  RooAddPdf sum_eta("sum_eta","sum_eta" , RooArgList(expo1_1, expo2_1, fit_eta_1, fit_A_1), RooArgList(nbkg1_1, nbkg2_1, nsig_eta_1, nsig_A_1)); 
  sum_eta.fitTo(data_eta, RooFit::Extended(), Range("eta"));

  RooAddPdf bkg("bkg","bkg" , RooArgList(expo1_1, expo2_1), RooArgList(nbkg1_1, nbkg2_1));

  //---------------------------------------------------//
  //GENERACJA-----------------------------------------//

  eta_M.setBins((600-min_mass)/5);
  RooDataSet *data_ar = (RooDataSet*) sum_eta.generate(eta_M, 1.1e5);
  RooDataHist *data_h = (RooDataHist*) data_ar->binnedClone();
  RooDataSet *data_bkg = (RooDataSet*) bkg.generate(eta_M, 1.1e5);
  RooDataHist *data_h_bkg = (RooDataHist*) data_bkg->binnedClone();

  RooPlot* frame1 = eta_M.frame(RooFit::Title("  "));
  data_h_bkg->plotOn(frame1, Name("Data"));
  bkg.plotOn(frame1, Name("bkg"));
  bkg.plotOn(frame1, LineStyle(9), LineColor(kRed), Range(365, 435), Name("bkg_A"));
  bkg.plotOn(frame1, LineStyle(9), LineColor(kRed), Range(525, 595), Name("bkg_eta"));
  frame1->GetXaxis()->SetTitle("m(#mu^{+}#mu^{-}) [MeV/c^{2}]");
  frame1->GetYaxis()->SetTitle("Candidates/5 [MeV/c^{2}]");
  frame1->Draw();

  TLegend *leg3 = new TLegend(0.62, 0.75, 0.89, 0.89);
  leg3->SetFillColor(kWhite);
  leg3->SetLineColor(kWhite);
  leg3->AddEntry("Data","Data","LP");
  leg3->AddEntry("bkg","Background only","LP");
  leg3->Draw();

  //---------------------------------------------------//
  //TEST CHI KWADRAT DLA PSEUDO - DANYCH---------------//

  cout<<"Chi kwadrat/ stopień swobody dla A, tło: "<<frame1->chiSquare("bkg_A", "Data", 4)<<" p-value: "<<ROOT::Math::chisquared_cdf_c(frame1->chiSquare("bkg_A", "Data", 4), 1)<<endl;
  cout<<"Chi kwadrat/ stopień swobody dla eta, tło: "<<frame1->chiSquare("bkg_eta", "Data", 4)<<" p-value: "<<ROOT::Math::chisquared_cdf_c(frame1->chiSquare("bkg_eta", "Data", 4), 1)<<endl;

  c1->SaveAs("Eta.png");   

}
