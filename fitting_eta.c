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

void fitting_eta(){
	
  TFile* input_file = new TFile("result.root", "OPEN");

  TH1* mass_data = dynamic_cast<TH1*>(input_file->Get("ID_1") ); 

  TCanvas* fitted=new TCanvas("Phi","Phi",0,0,900,700);

  Double_t min_mass = 450;
  Double_t max_mass = 650;

  Double_t mass_mean_eta = 552;
  Double_t signal_min_eta = 540;
  Double_t signal_max_eta = 560;

  RooRealVar Photon_M("A_M", "A_M", min_mass, max_mass);
  RooDataHist data("data","data",Photon_M,mass_data);

  gPad->SetLeftMargin(0.15);

  RooRealVar mean_eta("mean_eta","mean", mass_mean_eta, signal_min_eta, signal_max_eta);
  RooRealVar sigma_eta("sigma_eta","sigma", 5, 1.8, 15);
  RooRealVar width_eta("width_eta","width", 5, 1.8, 15);
  RooVoigtian fit_eta("#Eta","signal pdf eta",Photon_M,mean_eta,width_eta,sigma_eta);

  RooRealVar lambda1("lambda1", "slope1", -0.001, -1., 0.);
  RooExponential expo1("exponential_1", "exponential1 PDF", Photon_M, lambda1);

  RooRealVar lambda2("lambda2", "slope2", -0.01, -1., 0.);
  RooExponential expo2("exponential_2", "exponential2 PDF", Photon_M, lambda2);

  RooRealVar a("a", "a", -1, -50., 50.);
  RooRealVar b("b", "b", 500., -500., 500.);
  RooRealVar d("d", "d", 0.01, -50., 50.);
  RooPolynomial pol("pol", "pol", Photon_M, RooArgList(a, b));

  RooRealVar nsig_eta("nsig_eta","#signal events eta",100,1e2,1e5);
  RooRealVar nbkg1("nbkg1","#background events 1",500,100,1e5);
  RooRealVar nbkg2("nbkg2","#background events 2",500,100,1e5);
  RooRealVar nbkg3("nbkg3","#background events 3",500,100,1e5);

  Photon_M.setRange("bkg1", min_mass, 550);
  Photon_M.setRange("bkg2", 555, max_mass);
  Photon_M.setRange("all", min_mass, max_mass);

  RooAddPdf bkg("bkg","bkg" , RooArgList(expo1, expo2, pol), RooArgList(nbkg1, nbkg2, nbkg3)); 
  RooAddPdf signal("signal","signal" , RooArgList(fit_eta), RooArgList(nsig_eta)); 
   
  RooAddPdf sum("sum","sum" , RooArgList(bkg, signal)); 

  //---------------------------------------------------//
  //FITOWANIE------------------------------------------//
  sum.fitTo(data, RooFit::Extended(), Range("all"));

  RooArgSet* comps = sum.getComponents();
  RooAbsArg* bkg_arg = comps->find("bkg");
  RooAddPdf bkg_comp("bkg_comp","bkg_comp" , *bkg_arg); 


  RooPlot* frame = Photon_M.frame(RooFit::Title("  "));

  data.plotOn(frame, Name("Data"), DataError(RooAbsData::Poisson));

  sum.plotOn(frame, Components("sum"), LineWidth(3), Name("sum"));
  sum.plotOn(frame, Components("#Eta"), LineColor(kRed), LineStyle(9), LineWidth(2), Range(450, 600), Name("fit_eta"));
  sum.plotOn(frame, Components(RooArgSet(expo1, expo2, pol)), LineColor(kGreen), LineStyle(10), LineWidth(4), Name("bkg"));
  frame->GetXaxis()->SetTitle("m(#mu^{+}#mu^{-}) [MeV/c^{2}]");
  frame->GetYaxis()->SetTitle("Candidates/5 [MeV/c^{2}]");
  frame->Draw();

  //---------------------------------------------------//
  //TEST CHI KWADRAT DLA DANYCH------------------------//

  cout<<endl;
  cout<<"Chi kwadrat/ stopień swobody: "<<frame->chiSquare("sum", "Data", 11)<<" p-value: "<<ROOT::Math::chisquared_cdf_c(frame->chiSquare("sum", "Data", 11), 1)<<endl;
  cout<<"Chi kwadrat/ stopień swobody: "<<frame->chiSquare("bkg", "Data", 7)<<" p-value: "<<ROOT::Math::chisquared_cdf_c(frame->chiSquare("bkg", "Data", 7), 1)<<endl;

  TLegend *leg1 = new TLegend(0.18, 0.75, 0.43, 0.89);
  leg1->SetFillColor(kWhite);
  leg1->SetLineColor(kWhite);
  leg1->AddEntry("Data","Data","LP");
  leg1->AddEntry("sum","Signal + background","LP");
  leg1->AddEntry("bkg","Background only", "LP");
  leg1->Draw();

  TLatex *t = new TLatex(0.19, 0.72,"pp #sqrt{s}=13 TeV");
  t->SetTextSize(0.03);
  t->SetNDC(kTRUE);
  t->Draw();

  fitted->SaveAs("Eta_1.png");

  TCanvas* c1 = new TCanvas("eta","eta",900,700);
  gPad->SetLeftMargin(0.12);

  //---------------------------------------------------//
  //GENERACJA-----------------------------------------//

  Photon_M.setBins((max_mass-min_mass)/5);
  RooDataSet *data_ar = (RooDataSet*) sum.generate(Photon_M, Extended());
  RooDataHist *data_h = (RooDataHist*) data_ar->binnedClone();
  RooDataSet *data_bkg = (RooDataSet*) bkg.generate(Photon_M, Extended());
  RooDataHist *data_h_bkg = (RooDataHist*) data_bkg->binnedClone();

  RooPlot* frame1 = Photon_M.frame(RooFit::Title("  "));
  data_h_bkg->plotOn(frame1, Name("Data"), DataError(RooAbsData::Poisson));
  bkg.plotOn(frame1, Name("bkg"));
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

  cout<<"Chi kwadrat/ stopień swobody dla eta, tło: "<<frame1->chiSquare("bkg", "Data", 7)<<" p-value: "<<ROOT::Math::chisquared_cdf_c(frame1->chiSquare("bkg", "Data", 7), 1)<<endl;

  c1->SaveAs("Eta_pol.png"); 


}
