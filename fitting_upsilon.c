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

void fitting_upsilon(){
	
  TFile* input_file = new TFile("Plot.root", "OPEN");

  TCanvas c; 
  c.SetTitle("Upsilon");

  TH1* mass_data = dynamic_cast<TH1*>(input_file->Get("Upsilon") ); 

  TCanvas* fitted=new TCanvas("Upsilon","Upsilon",0,0,800,600);
                                //wartości używane dla 3 pików
  Double_t min_mass = 10150;     //1 - 8000  //2 - 9600  //3 - 10150
  Double_t max_mass = 12000;    //1 - 9800  //2 - 10180 //3 - 12000
  Double_t mass_mean = 10250;    //1 - 9450  //2 - 9900  //3 - 10250
  Double_t signal_min = 10150;   //1 - 9380  //2 - 9900  //3 - 10150
  Double_t signal_max = 10450;  //1 - 9500  //2 - 10100 //3 - 10450

  RooRealVar Photon_M("A_M", "A_M",mass_mean, min_mass, max_mass);
  RooDataHist data("data","data",Photon_M,mass_data);

  RooRealVar mean("mean","mean", mass_mean, signal_min, signal_max);
  RooRealVar sigma("sigma","sigma", 70, 0.0001, 2e3);
  RooRealVar width("width","width", 50, 10, 1e3);
  RooVoigtian fit("v","signal pdf",Photon_M,mean,width,sigma);

  RooRealVar lambda("lambda", "slope", 0., -5., 1.);
  RooExponential expo("exponential", "exponential PDF", Photon_M, lambda);

  RooRealVar nsig("nsig","#signal events",2000,1e3,1e7);
  RooRealVar nbkg("nbkg","#background events",20,1,1e5);

  Photon_M.setRange("signal", signal_min, signal_max); 
  Photon_M.setRange("all", min_mass, max_mass);


  RooAddPdf Zmodel("Bmodel","Bmodel" , RooArgList(fit), RooArgList(nsig));
  RooAddPdf sum("sum","sum" , RooArgList(expo, fit), RooArgList(nbkg, nsig)); 
  Zmodel.fitTo(data, RooFit::Extended(), Range("signal"));
  sum.fitTo(data, RooFit::Extended(), Range("all"));

  RooPlot* frame = Photon_M.frame(RooFit::Title("  "));

  data.plotOn(frame);

  sum.plotOn(frame);
  frame->GetXaxis()->SetTitle("#Upsilon mass [MeV/c^{2}]");
  frame->GetYaxis()->SetTitle("Candidates");
  frame->Draw();

  fitted->SaveAs("upsilon_3.png");

}
