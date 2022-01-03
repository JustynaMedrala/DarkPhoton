//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jul  6 12:11:19 2021 by ROOT version 6.22/02
// from TTree DecayTree/DecayTree
// found on file: A2mumu_1.root
//////////////////////////////////////////////////////////

#ifndef selector_high_h
#define selector_high_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TH1F.h>
#include <TCanvas.h>

// Headers needed by this particular selector


class selector_high : public TSelector {
public :
   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain
   TH1F *h_ID1 = nullptr;
   TH1F *h_ID2 = nullptr;
   TH1F *h_ID3 = nullptr;
   TH1F *h_high_mass = nullptr;

   // Readers to access the data (delete the ones you do not need).
   TTreeReaderValue<Double_t> A_ENDVERTEX_X = {fReader, "A_ENDVERTEX_X"};
   TTreeReaderValue<Double_t> A_ENDVERTEX_Y = {fReader, "A_ENDVERTEX_Y"};
   TTreeReaderValue<Double_t> A_ENDVERTEX_Z = {fReader, "A_ENDVERTEX_Z"};
   TTreeReaderValue<Double_t> A_ENDVERTEX_XERR = {fReader, "A_ENDVERTEX_XERR"};
   TTreeReaderValue<Double_t> A_ENDVERTEX_YERR = {fReader, "A_ENDVERTEX_YERR"};
   TTreeReaderValue<Double_t> A_ENDVERTEX_ZERR = {fReader, "A_ENDVERTEX_ZERR"};
   TTreeReaderValue<Double_t> A_ENDVERTEX_CHI2 = {fReader, "A_ENDVERTEX_CHI2"};
   TTreeReaderValue<Int_t> A_ENDVERTEX_NDOF = {fReader, "A_ENDVERTEX_NDOF"};
   TTreeReaderValue<Double_t> A_OWNPV_X = {fReader, "A_OWNPV_X"};
   TTreeReaderValue<Double_t> A_OWNPV_Y = {fReader, "A_OWNPV_Y"};
   TTreeReaderValue<Double_t> A_OWNPV_Z = {fReader, "A_OWNPV_Z"};
   TTreeReaderValue<Double_t> A_OWNPV_XERR = {fReader, "A_OWNPV_XERR"};
   TTreeReaderValue<Double_t> A_OWNPV_YERR = {fReader, "A_OWNPV_YERR"};
   TTreeReaderValue<Double_t> A_OWNPV_ZERR = {fReader, "A_OWNPV_ZERR"};
   TTreeReaderValue<Double_t> A_OWNPV_CHI2 = {fReader, "A_OWNPV_CHI2"};
   TTreeReaderValue<Int_t> A_OWNPV_NDOF = {fReader, "A_OWNPV_NDOF"};
   TTreeReaderValue<Double_t> A_IP_OWNPV = {fReader, "A_IP_OWNPV"};
   TTreeReaderValue<Double_t> A_IPCHI2_OWNPV = {fReader, "A_IPCHI2_OWNPV"};
   TTreeReaderValue<Double_t> A_FD_OWNPV = {fReader, "A_FD_OWNPV"};
   TTreeReaderValue<Double_t> A_FDCHI2_OWNPV = {fReader, "A_FDCHI2_OWNPV"};
   TTreeReaderValue<Double_t> A_DIRA_OWNPV = {fReader, "A_DIRA_OWNPV"};
   TTreeReaderValue<Double_t> A_P = {fReader, "A_P"};
   TTreeReaderValue<Double_t> A_PT = {fReader, "A_PT"};
   TTreeReaderValue<Double_t> A_PE = {fReader, "A_PE"};
   TTreeReaderValue<Double_t> A_PX = {fReader, "A_PX"};
   TTreeReaderValue<Double_t> A_PY = {fReader, "A_PY"};
   TTreeReaderValue<Double_t> A_PZ = {fReader, "A_PZ"};
   TTreeReaderValue<Double_t> A_MM = {fReader, "A_MM"};
   TTreeReaderValue<Double_t> A_MMERR = {fReader, "A_MMERR"};
   TTreeReaderValue<Double_t> A_M = {fReader, "A_M"};
   TTreeReaderValue<Int_t> A_ID = {fReader, "A_ID"};
   TTreeReaderValue<Double_t> mup_MC12TuneV2_ProbNNe = {fReader, "mup_MC12TuneV2_ProbNNe"};
   TTreeReaderValue<Double_t> mup_MC12TuneV2_ProbNNmu = {fReader, "mup_MC12TuneV2_ProbNNmu"};
   TTreeReaderValue<Double_t> mup_MC12TuneV2_ProbNNpi = {fReader, "mup_MC12TuneV2_ProbNNpi"};
   TTreeReaderValue<Double_t> mup_MC12TuneV2_ProbNNk = {fReader, "mup_MC12TuneV2_ProbNNk"};
   TTreeReaderValue<Double_t> mup_MC12TuneV2_ProbNNp = {fReader, "mup_MC12TuneV2_ProbNNp"};
   TTreeReaderValue<Double_t> mup_MC12TuneV2_ProbNNghost = {fReader, "mup_MC12TuneV2_ProbNNghost"};
   TTreeReaderValue<Double_t> mup_MC12TuneV3_ProbNNe = {fReader, "mup_MC12TuneV3_ProbNNe"};
   TTreeReaderValue<Double_t> mup_MC12TuneV3_ProbNNmu = {fReader, "mup_MC12TuneV3_ProbNNmu"};
   TTreeReaderValue<Double_t> mup_MC12TuneV3_ProbNNpi = {fReader, "mup_MC12TuneV3_ProbNNpi"};
   TTreeReaderValue<Double_t> mup_MC12TuneV3_ProbNNk = {fReader, "mup_MC12TuneV3_ProbNNk"};
   TTreeReaderValue<Double_t> mup_MC12TuneV3_ProbNNp = {fReader, "mup_MC12TuneV3_ProbNNp"};
   TTreeReaderValue<Double_t> mup_MC12TuneV3_ProbNNghost = {fReader, "mup_MC12TuneV3_ProbNNghost"};
   TTreeReaderValue<Double_t> mup_MC12TuneV4_ProbNNe = {fReader, "mup_MC12TuneV4_ProbNNe"};
   TTreeReaderValue<Double_t> mup_MC12TuneV4_ProbNNmu = {fReader, "mup_MC12TuneV4_ProbNNmu"};
   TTreeReaderValue<Double_t> mup_MC12TuneV4_ProbNNpi = {fReader, "mup_MC12TuneV4_ProbNNpi"};
   TTreeReaderValue<Double_t> mup_MC12TuneV4_ProbNNk = {fReader, "mup_MC12TuneV4_ProbNNk"};
   TTreeReaderValue<Double_t> mup_MC12TuneV4_ProbNNp = {fReader, "mup_MC12TuneV4_ProbNNp"};
   TTreeReaderValue<Double_t> mup_MC12TuneV4_ProbNNghost = {fReader, "mup_MC12TuneV4_ProbNNghost"};
   TTreeReaderValue<Double_t> mup_MC15TuneV1_ProbNNe = {fReader, "mup_MC15TuneV1_ProbNNe"};
   TTreeReaderValue<Double_t> mup_MC15TuneV1_ProbNNmu = {fReader, "mup_MC15TuneV1_ProbNNmu"};
   TTreeReaderValue<Double_t> mup_MC15TuneV1_ProbNNpi = {fReader, "mup_MC15TuneV1_ProbNNpi"};
   TTreeReaderValue<Double_t> mup_MC15TuneV1_ProbNNk = {fReader, "mup_MC15TuneV1_ProbNNk"};
   TTreeReaderValue<Double_t> mup_MC15TuneV1_ProbNNp = {fReader, "mup_MC15TuneV1_ProbNNp"};
   TTreeReaderValue<Double_t> mup_MC15TuneV1_ProbNNghost = {fReader, "mup_MC15TuneV1_ProbNNghost"};
   TTreeReaderValue<Double_t> mup_OWNPV_X = {fReader, "mup_OWNPV_X"};
   TTreeReaderValue<Double_t> mup_OWNPV_Y = {fReader, "mup_OWNPV_Y"};
   TTreeReaderValue<Double_t> mup_OWNPV_Z = {fReader, "mup_OWNPV_Z"};
   TTreeReaderValue<Double_t> mup_OWNPV_XERR = {fReader, "mup_OWNPV_XERR"};
   TTreeReaderValue<Double_t> mup_OWNPV_YERR = {fReader, "mup_OWNPV_YERR"};
   TTreeReaderValue<Double_t> mup_OWNPV_ZERR = {fReader, "mup_OWNPV_ZERR"};
   TTreeReaderValue<Double_t> mup_OWNPV_CHI2 = {fReader, "mup_OWNPV_CHI2"};
   TTreeReaderValue<Int_t> mup_OWNPV_NDOF = {fReader, "mup_OWNPV_NDOF"};
   TTreeReaderValue<Double_t> mup_IP_OWNPV = {fReader, "mup_IP_OWNPV"};
   TTreeReaderValue<Double_t> mup_IPCHI2_OWNPV = {fReader, "mup_IPCHI2_OWNPV"};
   TTreeReaderValue<Double_t> mup_ORIVX_X = {fReader, "mup_ORIVX_X"};
   TTreeReaderValue<Double_t> mup_ORIVX_Y = {fReader, "mup_ORIVX_Y"};
   TTreeReaderValue<Double_t> mup_ORIVX_Z = {fReader, "mup_ORIVX_Z"};
   TTreeReaderValue<Double_t> mup_ORIVX_XERR = {fReader, "mup_ORIVX_XERR"};
   TTreeReaderValue<Double_t> mup_ORIVX_YERR = {fReader, "mup_ORIVX_YERR"};
   TTreeReaderValue<Double_t> mup_ORIVX_ZERR = {fReader, "mup_ORIVX_ZERR"};
   TTreeReaderValue<Double_t> mup_ORIVX_CHI2 = {fReader, "mup_ORIVX_CHI2"};
   TTreeReaderValue<Int_t> mup_ORIVX_NDOF = {fReader, "mup_ORIVX_NDOF"};
   TTreeReaderValue<Double_t> mup_P = {fReader, "mup_P"};
   TTreeReaderValue<Double_t> mup_PT = {fReader, "mup_PT"};
   TTreeReaderValue<Double_t> mup_PE = {fReader, "mup_PE"};
   TTreeReaderValue<Double_t> mup_PX = {fReader, "mup_PX"};
   TTreeReaderValue<Double_t> mup_PY = {fReader, "mup_PY"};
   TTreeReaderValue<Double_t> mup_PZ = {fReader, "mup_PZ"};
   TTreeReaderValue<Double_t> mup_M = {fReader, "mup_M"};
   TTreeReaderValue<Int_t> mup_ID = {fReader, "mup_ID"};
   TTreeReaderValue<Double_t> mup_PIDe = {fReader, "mup_PIDe"};
   TTreeReaderValue<Double_t> mup_PIDmu = {fReader, "mup_PIDmu"};
   TTreeReaderValue<Double_t> mup_PIDK = {fReader, "mup_PIDK"};
   TTreeReaderValue<Double_t> mup_PIDp = {fReader, "mup_PIDp"};
   TTreeReaderValue<Double_t> mup_PIDd = {fReader, "mup_PIDd"};
   TTreeReaderValue<Double_t> mup_ProbNNe = {fReader, "mup_ProbNNe"};
   TTreeReaderValue<Double_t> mup_ProbNNk = {fReader, "mup_ProbNNk"};
   TTreeReaderValue<Double_t> mup_ProbNNp = {fReader, "mup_ProbNNp"};
   TTreeReaderValue<Double_t> mup_ProbNNpi = {fReader, "mup_ProbNNpi"};
   TTreeReaderValue<Double_t> mup_ProbNNmu = {fReader, "mup_ProbNNmu"};
   TTreeReaderValue<Double_t> mup_ProbNNd = {fReader, "mup_ProbNNd"};
   TTreeReaderValue<Double_t> mup_ProbNNghost = {fReader, "mup_ProbNNghost"};
   TTreeReaderValue<Bool_t> mup_hasMuon = {fReader, "mup_hasMuon"};
   TTreeReaderValue<Bool_t> mup_isMuon = {fReader, "mup_isMuon"};
   TTreeReaderValue<Bool_t> mup_hasRich = {fReader, "mup_hasRich"};
   TTreeReaderValue<Bool_t> mup_UsedRichAerogel = {fReader, "mup_UsedRichAerogel"};
   TTreeReaderValue<Bool_t> mup_UsedRich1Gas = {fReader, "mup_UsedRich1Gas"};
   TTreeReaderValue<Bool_t> mup_UsedRich2Gas = {fReader, "mup_UsedRich2Gas"};
   TTreeReaderValue<Bool_t> mup_RichAboveElThres = {fReader, "mup_RichAboveElThres"};
   TTreeReaderValue<Bool_t> mup_RichAboveMuThres = {fReader, "mup_RichAboveMuThres"};
   TTreeReaderValue<Bool_t> mup_RichAbovePiThres = {fReader, "mup_RichAbovePiThres"};
   TTreeReaderValue<Bool_t> mup_RichAboveKaThres = {fReader, "mup_RichAboveKaThres"};
   TTreeReaderValue<Bool_t> mup_RichAbovePrThres = {fReader, "mup_RichAbovePrThres"};
   TTreeReaderValue<Bool_t> mup_hasCalo = {fReader, "mup_hasCalo"};
   TTreeReaderValue<Double_t> mun_MC12TuneV2_ProbNNe = {fReader, "mun_MC12TuneV2_ProbNNe"};
   TTreeReaderValue<Double_t> mun_MC12TuneV2_ProbNNmu = {fReader, "mun_MC12TuneV2_ProbNNmu"};
   TTreeReaderValue<Double_t> mun_MC12TuneV2_ProbNNpi = {fReader, "mun_MC12TuneV2_ProbNNpi"};
   TTreeReaderValue<Double_t> mun_MC12TuneV2_ProbNNk = {fReader, "mun_MC12TuneV2_ProbNNk"};
   TTreeReaderValue<Double_t> mun_MC12TuneV2_ProbNNp = {fReader, "mun_MC12TuneV2_ProbNNp"};
   TTreeReaderValue<Double_t> mun_MC12TuneV2_ProbNNghost = {fReader, "mun_MC12TuneV2_ProbNNghost"};
   TTreeReaderValue<Double_t> mun_MC12TuneV3_ProbNNe = {fReader, "mun_MC12TuneV3_ProbNNe"};
   TTreeReaderValue<Double_t> mun_MC12TuneV3_ProbNNmu = {fReader, "mun_MC12TuneV3_ProbNNmu"};
   TTreeReaderValue<Double_t> mun_MC12TuneV3_ProbNNpi = {fReader, "mun_MC12TuneV3_ProbNNpi"};
   TTreeReaderValue<Double_t> mun_MC12TuneV3_ProbNNk = {fReader, "mun_MC12TuneV3_ProbNNk"};
   TTreeReaderValue<Double_t> mun_MC12TuneV3_ProbNNp = {fReader, "mun_MC12TuneV3_ProbNNp"};
   TTreeReaderValue<Double_t> mun_MC12TuneV3_ProbNNghost = {fReader, "mun_MC12TuneV3_ProbNNghost"};
   TTreeReaderValue<Double_t> mun_MC12TuneV4_ProbNNe = {fReader, "mun_MC12TuneV4_ProbNNe"};
   TTreeReaderValue<Double_t> mun_MC12TuneV4_ProbNNmu = {fReader, "mun_MC12TuneV4_ProbNNmu"};
   TTreeReaderValue<Double_t> mun_MC12TuneV4_ProbNNpi = {fReader, "mun_MC12TuneV4_ProbNNpi"};
   TTreeReaderValue<Double_t> mun_MC12TuneV4_ProbNNk = {fReader, "mun_MC12TuneV4_ProbNNk"};
   TTreeReaderValue<Double_t> mun_MC12TuneV4_ProbNNp = {fReader, "mun_MC12TuneV4_ProbNNp"};
   TTreeReaderValue<Double_t> mun_MC12TuneV4_ProbNNghost = {fReader, "mun_MC12TuneV4_ProbNNghost"};
   TTreeReaderValue<Double_t> mun_MC15TuneV1_ProbNNe = {fReader, "mun_MC15TuneV1_ProbNNe"};
   TTreeReaderValue<Double_t> mun_MC15TuneV1_ProbNNmu = {fReader, "mun_MC15TuneV1_ProbNNmu"};
   TTreeReaderValue<Double_t> mun_MC15TuneV1_ProbNNpi = {fReader, "mun_MC15TuneV1_ProbNNpi"};
   TTreeReaderValue<Double_t> mun_MC15TuneV1_ProbNNk = {fReader, "mun_MC15TuneV1_ProbNNk"};
   TTreeReaderValue<Double_t> mun_MC15TuneV1_ProbNNp = {fReader, "mun_MC15TuneV1_ProbNNp"};
   TTreeReaderValue<Double_t> mun_MC15TuneV1_ProbNNghost = {fReader, "mun_MC15TuneV1_ProbNNghost"};
   TTreeReaderValue<Double_t> mun_OWNPV_X = {fReader, "mun_OWNPV_X"};
   TTreeReaderValue<Double_t> mun_OWNPV_Y = {fReader, "mun_OWNPV_Y"};
   TTreeReaderValue<Double_t> mun_OWNPV_Z = {fReader, "mun_OWNPV_Z"};
   TTreeReaderValue<Double_t> mun_OWNPV_XERR = {fReader, "mun_OWNPV_XERR"};
   TTreeReaderValue<Double_t> mun_OWNPV_YERR = {fReader, "mun_OWNPV_YERR"};
   TTreeReaderValue<Double_t> mun_OWNPV_ZERR = {fReader, "mun_OWNPV_ZERR"};
   TTreeReaderValue<Double_t> mun_OWNPV_CHI2 = {fReader, "mun_OWNPV_CHI2"};
   TTreeReaderValue<Int_t> mun_OWNPV_NDOF = {fReader, "mun_OWNPV_NDOF"};
   TTreeReaderValue<Double_t> mun_IP_OWNPV = {fReader, "mun_IP_OWNPV"};
   TTreeReaderValue<Double_t> mun_IPCHI2_OWNPV = {fReader, "mun_IPCHI2_OWNPV"};
   TTreeReaderValue<Double_t> mun_ORIVX_X = {fReader, "mun_ORIVX_X"};
   TTreeReaderValue<Double_t> mun_ORIVX_Y = {fReader, "mun_ORIVX_Y"};
   TTreeReaderValue<Double_t> mun_ORIVX_Z = {fReader, "mun_ORIVX_Z"};
   TTreeReaderValue<Double_t> mun_ORIVX_XERR = {fReader, "mun_ORIVX_XERR"};
   TTreeReaderValue<Double_t> mun_ORIVX_YERR = {fReader, "mun_ORIVX_YERR"};
   TTreeReaderValue<Double_t> mun_ORIVX_ZERR = {fReader, "mun_ORIVX_ZERR"};
   TTreeReaderValue<Double_t> mun_ORIVX_CHI2 = {fReader, "mun_ORIVX_CHI2"};
   TTreeReaderValue<Int_t> mun_ORIVX_NDOF = {fReader, "mun_ORIVX_NDOF"};
   TTreeReaderValue<Double_t> mun_P = {fReader, "mun_P"};
   TTreeReaderValue<Double_t> mun_PT = {fReader, "mun_PT"};
   TTreeReaderValue<Double_t> mun_PE = {fReader, "mun_PE"};
   TTreeReaderValue<Double_t> mun_PX = {fReader, "mun_PX"};
   TTreeReaderValue<Double_t> mun_PY = {fReader, "mun_PY"};
   TTreeReaderValue<Double_t> mun_PZ = {fReader, "mun_PZ"};
   TTreeReaderValue<Double_t> mun_M = {fReader, "mun_M"};
   TTreeReaderValue<Int_t> mun_ID = {fReader, "mun_ID"};
   TTreeReaderValue<Double_t> mun_PIDe = {fReader, "mun_PIDe"};
   TTreeReaderValue<Double_t> mun_PIDmu = {fReader, "mun_PIDmu"};
   TTreeReaderValue<Double_t> mun_PIDK = {fReader, "mun_PIDK"};
   TTreeReaderValue<Double_t> mun_PIDp = {fReader, "mun_PIDp"};
   TTreeReaderValue<Double_t> mun_PIDd = {fReader, "mun_PIDd"};
   TTreeReaderValue<Double_t> mun_ProbNNe = {fReader, "mun_ProbNNe"};
   TTreeReaderValue<Double_t> mun_ProbNNk = {fReader, "mun_ProbNNk"};
   TTreeReaderValue<Double_t> mun_ProbNNp = {fReader, "mun_ProbNNp"};
   TTreeReaderValue<Double_t> mun_ProbNNpi = {fReader, "mun_ProbNNpi"};
   TTreeReaderValue<Double_t> mun_ProbNNmu = {fReader, "mun_ProbNNmu"};
   TTreeReaderValue<Double_t> mun_ProbNNd = {fReader, "mun_ProbNNd"};
   TTreeReaderValue<Double_t> mun_ProbNNghost = {fReader, "mun_ProbNNghost"};
   TTreeReaderValue<Bool_t> mun_hasMuon = {fReader, "mun_hasMuon"};
   TTreeReaderValue<Bool_t> mun_isMuon = {fReader, "mun_isMuon"};
   TTreeReaderValue<Bool_t> mun_hasRich = {fReader, "mun_hasRich"};
   TTreeReaderValue<Bool_t> mun_UsedRichAerogel = {fReader, "mun_UsedRichAerogel"};
   TTreeReaderValue<Bool_t> mun_UsedRich1Gas = {fReader, "mun_UsedRich1Gas"};
   TTreeReaderValue<Bool_t> mun_UsedRich2Gas = {fReader, "mun_UsedRich2Gas"};
   TTreeReaderValue<Bool_t> mun_RichAboveElThres = {fReader, "mun_RichAboveElThres"};
   TTreeReaderValue<Bool_t> mun_RichAboveMuThres = {fReader, "mun_RichAboveMuThres"};
   TTreeReaderValue<Bool_t> mun_RichAbovePiThres = {fReader, "mun_RichAbovePiThres"};
   TTreeReaderValue<Bool_t> mun_RichAboveKaThres = {fReader, "mun_RichAboveKaThres"};
   TTreeReaderValue<Bool_t> mun_RichAbovePrThres = {fReader, "mun_RichAbovePrThres"};
   TTreeReaderValue<Bool_t> mun_hasCalo = {fReader, "mun_hasCalo"};
   TTreeReaderValue<UInt_t> nCandidate = {fReader, "nCandidate"};
   TTreeReaderValue<ULong64_t> totCandidates = {fReader, "totCandidates"};
   TTreeReaderValue<ULong64_t> EventInSequence = {fReader, "EventInSequence"};
   TTreeReaderValue<UInt_t> runNumber = {fReader, "runNumber"};
   TTreeReaderValue<ULong64_t> eventNumber = {fReader, "eventNumber"};
   TTreeReaderValue<UInt_t> BCID = {fReader, "BCID"};
   TTreeReaderValue<Int_t> BCType = {fReader, "BCType"};
   TTreeReaderValue<UInt_t> OdinTCK = {fReader, "OdinTCK"};
   TTreeReaderValue<UInt_t> L0DUTCK = {fReader, "L0DUTCK"};
   TTreeReaderValue<UInt_t> HLT1TCK = {fReader, "HLT1TCK"};
   TTreeReaderValue<UInt_t> HLT2TCK = {fReader, "HLT2TCK"};
   TTreeReaderValue<ULong64_t> GpsTime = {fReader, "GpsTime"};
   TTreeReaderValue<Short_t> Polarity = {fReader, "Polarity"};


   selector_high(TTree * /*tree*/ =0) { }
   virtual ~selector_high() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   ClassDef(selector_high,0);

};

#endif

#ifdef selector_high_cxx
void selector_high::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
}

Bool_t selector_high::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef selector_high_cxx
