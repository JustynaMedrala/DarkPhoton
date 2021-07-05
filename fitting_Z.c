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
	
  TFile* input_file = new TFile("Plot.root", "OPEN");

  TH1* mass_data = dynamic_cast<TH1*>(input_file->Get("Z") ); 

  TCanvas* fitted=new TCanvas("Z","Z",0,0,800,600);

  Double_t min_mass = 30000;
  Double_t max_mass = 120000;
  Double_t mass_mean = 90000;
  Double_t signal_min = 82000;
  Double_t signal_max = 95000;

  RooRealVar Photon_M("A_M", "A_M",mass_mean, min_mass, max_mass);
  RooDataHist data("data","data",Photon_M,mass_data);

  RooRealVar mean("mean","mean", mass_mean, signal_min, signal_max);//90e3, 9.3e3, 10.2e3);//90e3, 90e3, 100e3);
  RooRealVar sigma("sigma","sigma", 700, 1e2, 2e4);//70, 0.001, 200);//700, 200, 20000);
  RooRealVar width("width","width", 1000, 500, 1e5);//50, 1, 1000);//500, 100, 10000);
  RooVoigtian fit("v","signal pdf",Photon_M,mean,width,sigma);

  RooRealVar lambda("lambda", "slope", 0., -1., 0.5);
  RooExponential expo("exponential", "exponential PDF", Photon_M, lambda);

  RooRealVar nsig("nsig","#signal events",2000,1e3,1e4);
  RooRealVar nbkg("nbkg","#background events",500,500,1e3);

  Photon_M.setRange("signal", signal_min, signal_max); 
  Photon_M.setRange("all", min_mass, max_mass);


  RooAddPdf Zmodel("Zmodel","Zmodel" , RooArgList(fit), RooArgList(nsig));
  RooAddPdf sum("sum","sum" , RooArgList(expo, fit), RooArgList(nbkg, nsig)); 
  Zmodel.fitTo(data, RooFit::Extended(), Range("signal"));
  sum.fitTo(data, RooFit::Extended(), Range("all"));

  //RooRealVar Photon_M1("A_M", "A_M",8000,12000);
  //RooDataHist data1("data","data",Photon_M1,RooFit::Import(*h_A_MM));

  RooPlot* frame = Photon_M.frame(RooFit::Title("  "));
  //RooPlot* frame1 = Photon_M1.frame(RooFit::Title("#Upsilon"));

  data.plotOn(frame);
  //data1.plotOn(frame1);

  sum.plotOn(frame);
  frame->GetXaxis()->SetTitle("Z mass [MeV/c^{2}]");
  frame->GetYaxis()->SetTitle("Candidates");
  frame->Draw();

  //RooRealVar mean_upsilon("mean","mean", 9000, 7000, 11000);
  //RooRealVar sigma_upsilon("sigma","sigma", 50, 1, 1000);
  //RooRealVar width_upsilon("width","width", 50, 1, 100);
  //RooVoigtian fit_upsilon("v","signal pdf",Photon_M1,mean_phi,width_phi,sigma_phi);

  //fitted->cd(2);
  //fit_upsilon.fitTo(data1);
  //fit_upsilon.plotOn(frame1);
  //frame1->GetXaxis()->SetTitle("mass [MeV/c^{2}]");
  //frame1->GetYaxis()->SetTitle("Candidates");
  //frame1->GetYaxis()->SetRangeUser(100,250);
  //frame1->Draw();

  fitted->SaveAs("Z.png");

}
