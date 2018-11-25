
#include <cstdlib>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <map>
#include <string>
#include <iomanip>
#include <stdio.h>

#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TChain.h"
#include "TBranch.h"
#include "TLorentzVector.h"
#include "TSpectrum.h"
#include "TVirtualFFT.h"
#include "TRandom.h"


char target[50] = "H1426";
char method[50] = "MSCW";

double Elev_cut_lower = 55;
double Elev_cut_upper = 85;

double MSCW_cut_lower = -0.5;
double MSCW_cut_upper = 0.6;
double MSCL_cut_lower = -0.5;
double MSCL_cut_upper = 0.6;
double Depth_cut_lower = 7;
double Depth_cut_upper = 13;


string  filename;
double TelElevation = 0;
double ErecS = 0;
double EChi2S = 0;
double MSCW = 0;
double MSCL = 0;
double SlantDepth = 0;

//const int N_energy_bins = 1;
//double energy_bins[N_energy_bins+1] = {pow(10,2.0),pow(10,2.5)};
const int N_energy_bins = 23;
double energy_bins[N_energy_bins+1] = {pow(10,2.1),pow(10,2.2),pow(10,2.3),pow(10,2.4),pow(10,2.5),pow(10,2.6),pow(10,2.7),pow(10,2.8),pow(10,2.9),pow(10,3.0),pow(10,3.1),pow(10,3.2),pow(10,3.3),pow(10,3.4),pow(10,3.5),pow(10,3.6),pow(10,3.7),pow(10,3.8),pow(10,3.9),pow(10,4.0),pow(10,4.1),pow(10,4.2),pow(10,4.3),pow(10,4.4)};

const int N_bins_for_deconv = int(pow(2,10));

vector<int> GetRunList(string source) {
        vector<int> list;
        if (source=="H1426") {
                list.push_back(84601);
                list.push_back(84602);
                list.push_back(84603);
                list.push_back(84604);
                list.push_back(84795);
                list.push_back(85012);
                list.push_back(85013);
                list.push_back(85014);
                list.push_back(85015);
                list.push_back(85121);
                list.push_back(85122);
                list.push_back(85123);
                list.push_back(85604);
                list.push_back(85605);
                list.push_back(85831);
                list.push_back(85832);
                list.push_back(85833);
                list.push_back(85834);
                list.push_back(85999);
                list.push_back(86000);
                list.push_back(86001);
                list.push_back(86002);
                list.push_back(86187);
                list.push_back(86188);
                list.push_back(86189);
                list.push_back(86690);
                list.push_back(89102);
                list.push_back(89293);
                list.push_back(89372);
                list.push_back(89373);
                list.push_back(89374);
                list.push_back(89712);
                list.push_back(89713);
                list.push_back(89783);
                list.push_back(89784);
                list.push_back(89924);
                list.push_back(89925);
                list.push_back(89968);
                list.push_back(89969);
                list.push_back(89990);
                list.push_back(90019);
                list.push_back(90044);
                list.push_back(90119);
                list.push_back(90120);
                list.push_back(90239);
                list.push_back(90240);
                list.push_back(90241);
                list.push_back(90242);
                list.push_back(90279);
                list.push_back(90281);
                list.push_back(90322);
                list.push_back(90323);
                list.push_back(90325);
                list.push_back(90341);
                list.push_back(90370);
                list.push_back(90371);
                list.push_back(90419);
                list.push_back(90448);
                list.push_back(90449);
                list.push_back(90549);
                list.push_back(90550);
                list.push_back(90551);
                list.push_back(90584);
        }
        if (source=="Crab") {
                list.push_back(85232);
                list.push_back(85244);
                list.push_back(85245);
                list.push_back(85329);
                list.push_back(85373);
                list.push_back(85374);
                list.push_back(85375);
                list.push_back(85386);
                list.push_back(85387);
                list.push_back(85414);
                list.push_back(85415);
                list.push_back(85648);
                list.push_back(85649);
                list.push_back(88690);
                list.push_back(88714);
                list.push_back(88715);
                list.push_back(88716);
                list.push_back(88719);
                list.push_back(88720);
                list.push_back(88742);
                list.push_back(88746);
                list.push_back(88747);
                list.push_back(88772);
                list.push_back(88796);
                list.push_back(86931);
                list.push_back(86932);
                list.push_back(86952);
                list.push_back(86953);
                list.push_back(86972);
                list.push_back(86993);
                list.push_back(86994);
                list.push_back(87054);
                list.push_back(87097);
                list.push_back(87098);
                list.push_back(87165);
                list.push_back(87166);
                list.push_back(87190);
                list.push_back(87214);
                list.push_back(87215);
                list.push_back(87220);
                list.push_back(87221);
                list.push_back(87375);
                list.push_back(87376);
                list.push_back(87395);
                list.push_back(87399);
                list.push_back(87400);
                list.push_back(87401);
                list.push_back(87441);
                list.push_back(87510);
                list.push_back(87511);
                list.push_back(87534);
                list.push_back(87535);
                list.push_back(87575);
                list.push_back(87576);
                list.push_back(87605);
                list.push_back(87606);
                list.push_back(87658);
                list.push_back(87716);
                list.push_back(87717);
                list.push_back(87733);
                list.push_back(87970);
                list.push_back(87971);
                list.push_back(88026);
                list.push_back(88098);
                list.push_back(88099);
                list.push_back(88129);
                list.push_back(88130);
                list.push_back(88132);
                list.push_back(88133);
                list.push_back(88134);
                list.push_back(88140);
                list.push_back(88141);
                list.push_back(88148);
                list.push_back(88149);
                list.push_back(88153);
                list.push_back(88225);
                list.push_back(88254);
                list.push_back(88269);
                list.push_back(88272);
                list.push_back(88298);
                list.push_back(88299);
                list.push_back(88314);
                list.push_back(88315);
                list.push_back(88316);
                list.push_back(88317);
                list.push_back(88344);
                list.push_back(88363);
                list.push_back(88364);
                list.push_back(88408);
                list.push_back(88433);
                list.push_back(88434);
                list.push_back(88454);
                list.push_back(88455);
                list.push_back(88456);
                list.push_back(88457);
                list.push_back(88458);
                list.push_back(88459);
                list.push_back(88489);
                list.push_back(88490);
                list.push_back(88499);
                list.push_back(88500);
                list.push_back(88509);
                list.push_back(88521);
                list.push_back(88539);
                list.push_back(88561);
                list.push_back(88562);
                list.push_back(88569);
                list.push_back(88570);
                list.push_back(88576);
                list.push_back(88577);
                list.push_back(88590);
                list.push_back(88633);
                list.push_back(88649);
                list.push_back(88666);
                list.push_back(88673);
                list.push_back(88710);
                list.push_back(88711);
                list.push_back(88712);
                list.push_back(88713);
                list.push_back(88721);
                list.push_back(88722);
                list.push_back(88723);
                list.push_back(88724);
                list.push_back(88725);
                list.push_back(88737);
                list.push_back(88743);
                list.push_back(88748);
                list.push_back(88749);
                list.push_back(88750);
                list.push_back(88751);
                list.push_back(88752);
                list.push_back(88822);
                list.push_back(88842);
                list.push_back(88882);
                list.push_back(88883);
                list.push_back(88884);
                list.push_back(88907);
                list.push_back(88908);
                list.push_back(88909);
                list.push_back(88942);
                list.push_back(88943);
                list.push_back(88944);
                list.push_back(89029);
                list.push_back(89077);
                list.push_back(89082);
                list.push_back(89107);
                list.push_back(89138);
                list.push_back(89139);
                list.push_back(89166);
                list.push_back(89167);
                list.push_back(89234);
                list.push_back(89394);
                list.push_back(89395);
                list.push_back(89507);
        }
        return list;
}
bool SignalSelection() {
    //if (EChi2S>0.2) return false;
    if (MSCL<MSCL_cut_lower) return false;
    if (MSCL>MSCL_cut_upper) return false;
    if (SlantDepth*100./37.<Depth_cut_lower) return false;
    if (SlantDepth*100./37.>Depth_cut_upper) return false;
    return true;
}
bool AuxSignalSelection() {
    //if (EChi2S>0.2) return false;
    if (MSCL<MSCL_cut_lower) return false;
    if (MSCL>MSCL_cut_upper) return false;
    if (SlantDepth*100./37.>Depth_cut_lower && SlantDepth*100./37.<Depth_cut_upper) return false;
    return true;
}
bool ControlSelection() {
    //if (EChi2S>0.2) return false;
    if (MSCL<MSCL_cut_upper*1.5) return false;
    if (MSCL>MSCL_cut_upper*3.) return false;
    if (SlantDepth*100./37.<Depth_cut_lower) return false;
    if (SlantDepth*100./37.>Depth_cut_upper) return false;
    return true;
}
bool AuxControlSelection() {
    //if (EChi2S>0.2) return false;
    if (MSCL<MSCL_cut_upper*1.5) return false;
    if (MSCL>MSCL_cut_upper*3.) return false;
    if (SlantDepth*100./37.>Depth_cut_lower && SlantDepth*100./37.<Depth_cut_upper) return false;
    return true;
}

void Deconvolution(TH1* Hist_source, TH1* Hist_response, TH1* Hist_Deconv) {
        Double_t *source = new Double_t[N_bins_for_deconv];
        Double_t *response = new Double_t[N_bins_for_deconv];
        for (int i=0;i<N_bins_for_deconv;i++) {
            source[i] = max(Hist_source->GetBinContent(i+1),0.);
            if (i<N_bins_for_deconv/2) response[i] = max(Hist_response->GetBinContent(i+1+N_bins_for_deconv/2),0.);
            else response[i] = 0;
        }
        TSpectrum sp;
        sp.Deconvolution(source,response,N_bins_for_deconv,50,1,1);
        for (int i=0;i<N_bins_for_deconv;i++) {
            Hist_Deconv->SetBinContent(i+1,source[i]);
        }
}

void DeconvolutionMethodForExtendedSources() {
        TRandom rnd;
        vector<TH1D> Hist_CrabOff_SR_ErecS;
        vector<TH1D> Hist_CrabOff_SR_MSCW;
        vector<TH1D> Hist_CrabOff_CR_MSCW;
        vector<TH1D> Hist_CrabOff_Deconv_MSCW;
        vector<TH1D> Hist_CrabOff_BkgTemp_MSCW;
        vector<TH1D> Hist_CrabOff_Bkg_MSCW;
        vector<TH1D> Hist_Target_SR_ErecS;
        vector<TH1D> Hist_Target_SR_MSCW;
        vector<TH1D> Hist_Target_CR_MSCW;
        vector<TH1D> Hist_Target_ASR_MSCW;
        vector<TH1D> Hist_Target_ACR_MSCW;
        vector<TH1D> Hist_Target_BkgTemp_MSCW;
        vector<TH1D> Hist_Target_Bkg_MSCW;
        vector<TH1D> Hist_Target_ABkgTemp_MSCW;
        vector<TH1D> Hist_Target_ABkg_MSCW;
        vector<TH1D> Hist_Target_Deconv_MSCW;
        for (int e=0;e<N_energy_bins;e++) {
            char e_low[50];
            sprintf(e_low, "%i", int(energy_bins[e]));
            char e_up[50];
            sprintf(e_up, "%i", int(energy_bins[e+1]));
            Hist_CrabOff_SR_ErecS.push_back(TH1D("Hist_CrabOff_SR_ErecS_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_energy_bins,energy_bins));
            Hist_CrabOff_SR_MSCW.push_back(TH1D("Hist_CrabOff_SR_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-100,50));
            Hist_CrabOff_CR_MSCW.push_back(TH1D("Hist_CrabOff_CR_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-100,50));
            Hist_CrabOff_Deconv_MSCW.push_back(TH1D("Hist_CrabOff_Deconv_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-100,50));
            Hist_CrabOff_BkgTemp_MSCW.push_back(TH1D("Hist_CrabOff_BkgTemp_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-100,50));
            Hist_CrabOff_Bkg_MSCW.push_back(TH1D("Hist_CrabOff_Bkg_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-100,50));
            Hist_Target_SR_ErecS.push_back(TH1D("Hist_Target_SR_ErecS_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_energy_bins,energy_bins));
            Hist_Target_SR_MSCW.push_back(TH1D("Hist_Target_SR_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-100,50));
            Hist_Target_CR_MSCW.push_back(TH1D("Hist_Target_CR_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-100,50));
            Hist_Target_ASR_MSCW.push_back(TH1D("Hist_Target_ASR_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-100,50));
            Hist_Target_ACR_MSCW.push_back(TH1D("Hist_Target_ACR_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-100,50));
            Hist_Target_BkgTemp_MSCW.push_back(TH1D("Hist_Target_BkgTemp_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-100,50));
            Hist_Target_Bkg_MSCW.push_back(TH1D("Hist_Target_Bkg_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-100,50));
            Hist_Target_ABkgTemp_MSCW.push_back(TH1D("Hist_Target_ABkgTemp_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-100,50));
            Hist_Target_ABkg_MSCW.push_back(TH1D("Hist_Target_ABkg_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-100,50));
            Hist_Target_Deconv_MSCW.push_back(TH1D("Hist_Target_Deconv_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-100,50));
        }
        vector<int> Crab_runlist = GetRunList("Crab");
        for (int run=0;run<Crab_runlist.size();run++) {
                char run_number[50];
                sprintf(run_number, "%i", int(Crab_runlist.at(run)));
                filename = TString("$VERITAS_USER_DATA_DIR/Crab_V6_Moderate-TMVA-BDT.RB."+TString(run_number)+".root");
                TFile*  input_file = TFile::Open(filename.c_str());
                std::cout << "Open " << filename << std::endl;
                TTree* pointing_tree = (TTree*) input_file->Get("run_"+TString(run_number)+"/stereo/pointingDataReduced");
                pointing_tree->SetBranchAddress("TelElevation",&TelElevation);
                pointing_tree->GetEntry(0);
                std::cout << TelElevation << std::endl;
                if (TelElevation<Elev_cut_lower) continue;
                if (TelElevation>Elev_cut_upper) continue;

                TTree* CrabOff_tree = (TTree*) input_file->Get("run_"+TString(run_number)+"/stereo/data_off");
                CrabOff_tree->SetBranchAddress("ErecS",&ErecS);
                CrabOff_tree->SetBranchAddress("EChi2S",&EChi2S);
                CrabOff_tree->SetBranchAddress("MSCW",&MSCW);
                CrabOff_tree->SetBranchAddress("MSCL",&MSCL);
                CrabOff_tree->SetBranchAddress("SlantDepth",&SlantDepth);
                for (int entry=0;entry<CrabOff_tree->GetEntries();entry++) {
                        CrabOff_tree->GetEntry(entry);
                        int e = Hist_CrabOff_SR_ErecS.at(0).FindBin(ErecS*1000.)-1;
                        if (e>=N_energy_bins) continue;
                        if (e<0) continue;
                        if (SignalSelection()) Hist_CrabOff_SR_MSCW.at(e).Fill(MSCW);
                        if (ControlSelection()) Hist_CrabOff_CR_MSCW.at(e).Fill(MSCW);
                        //MSCW = rnd.Gaus (1.0, 3.0);
                        //Hist_CrabOff_SR_MSCW.at(e).Fill(MSCW);
                        //MSCW = rnd.Gaus (2.0, 2.0);
                        //Hist_CrabOff_CR_MSCW.at(e).Fill(MSCW);
                }
                input_file->Close();

        }

        vector<int> Target_runlist = GetRunList(target);
        for (int run=0;run<Target_runlist.size();run++) {
                char run_number[50];
                sprintf(run_number, "%i", int(Target_runlist[run]));
                filename = TString("$VERITAS_USER_DATA_DIR/H1426_V6_Moderate-TMVA-BDT.RB."+TString(run_number)+".root");
                TFile*  input_file = TFile::Open(filename.c_str());
                std::cout << "Open " << filename << std::endl;
                TTree* pointing_tree = (TTree*) input_file->Get("run_"+TString(run_number)+"/stereo/pointingDataReduced");
                pointing_tree->SetBranchAddress("TelElevation",&TelElevation);
                pointing_tree->GetEntry(0);
                std::cout << TelElevation << std::endl;
                if (TelElevation<Elev_cut_lower) continue;
                if (TelElevation>Elev_cut_upper) continue;

                TTree* Target_tree = (TTree*) input_file->Get("run_"+TString(run_number)+"/stereo/data_on");
                Target_tree->SetBranchAddress("ErecS",&ErecS);
                Target_tree->SetBranchAddress("EChi2S",&EChi2S);
                Target_tree->SetBranchAddress("MSCW",&MSCW);
                Target_tree->SetBranchAddress("MSCL",&MSCL);
                Target_tree->SetBranchAddress("SlantDepth",&SlantDepth);
                for (int entry=0;entry<Target_tree->GetEntries();entry++) {
                        Target_tree->GetEntry(entry);
                        int e = Hist_Target_SR_ErecS.at(0).FindBin(ErecS*1000.)-1;
                        if (e>=N_energy_bins) continue;
                        if (e<0) continue;
                        if (SignalSelection()) Hist_Target_SR_MSCW.at(e).Fill(MSCW);
                        if (ControlSelection()) Hist_Target_CR_MSCW.at(e).Fill(MSCW);
                        if (AuxSignalSelection()) Hist_Target_ASR_MSCW.at(e).Fill(MSCW);
                        if (AuxControlSelection()) Hist_Target_ACR_MSCW.at(e).Fill(MSCW);
                }
                input_file->Close();

        }

        for (int e=0;e<N_energy_bins;e++) {
                Deconvolution(&Hist_CrabOff_CR_MSCW.at(e), &Hist_CrabOff_SR_MSCW.at(e), &Hist_CrabOff_Deconv_MSCW.at(e));
                Deconvolution(&Hist_CrabOff_CR_MSCW.at(e), &Hist_CrabOff_Deconv_MSCW.at(e), &Hist_CrabOff_BkgTemp_MSCW.at(e));
                Deconvolution(&Hist_Target_ACR_MSCW.at(e), &Hist_Target_ASR_MSCW.at(e), &Hist_Target_Deconv_MSCW.at(e));
                Deconvolution(&Hist_Target_ACR_MSCW.at(e), &Hist_Target_Deconv_MSCW.at(e), &Hist_Target_ABkgTemp_MSCW.at(e));
                Deconvolution(&Hist_Target_CR_MSCW.at(e), &Hist_Target_Deconv_MSCW.at(e), &Hist_Target_BkgTemp_MSCW.at(e));
        }

        for (int e=0;e<N_energy_bins;e++) {
                double offset = Hist_CrabOff_SR_MSCW.at(e).GetMean()-Hist_CrabOff_BkgTemp_MSCW.at(e).GetMean();
                for (int i=0;i<N_bins_for_deconv;i++) {
                        int b = Hist_CrabOff_Bkg_MSCW.at(e).FindBin(Hist_CrabOff_Bkg_MSCW.at(e).GetBinCenter(i+1)-offset);
                        Hist_CrabOff_Bkg_MSCW.at(e).SetBinContent(i+1,Hist_CrabOff_BkgTemp_MSCW.at(e).GetBinContent(b));
                        Hist_CrabOff_Bkg_MSCW.at(e).SetBinError(i+1,Hist_CrabOff_BkgTemp_MSCW.at(e).GetBinError(b));
                }
                int norm_bin_low_crab = Hist_CrabOff_Bkg_MSCW.at(e).FindBin(3*MSCW_cut_upper);
                int norm_bin_up_crab = Hist_CrabOff_Bkg_MSCW.at(e).FindBin(100*MSCW_cut_upper);
                double scale_crab = Hist_CrabOff_SR_MSCW.at(e).Integral(norm_bin_low_crab,norm_bin_up_crab)/Hist_CrabOff_Bkg_MSCW.at(e).Integral(norm_bin_low_crab,norm_bin_up_crab);
                Hist_CrabOff_Bkg_MSCW.at(e).Scale(scale_crab);

                offset = Hist_Target_ASR_MSCW.at(e).GetMean()-Hist_Target_ABkgTemp_MSCW.at(e).GetMean();
                for (int i=0;i<N_bins_for_deconv;i++) {
                        int b = Hist_Target_Bkg_MSCW.at(e).FindBin(Hist_Target_Bkg_MSCW.at(e).GetBinCenter(i+1)-offset);
                        Hist_Target_ABkg_MSCW.at(e).SetBinContent(i+1,Hist_Target_ABkgTemp_MSCW.at(e).GetBinContent(b));
                        Hist_Target_ABkg_MSCW.at(e).SetBinError(i+1,Hist_Target_ABkgTemp_MSCW.at(e).GetBinError(b));
                        Hist_Target_Bkg_MSCW.at(e).SetBinContent(i+1,Hist_Target_BkgTemp_MSCW.at(e).GetBinContent(b));
                        Hist_Target_Bkg_MSCW.at(e).SetBinError(i+1,Hist_Target_BkgTemp_MSCW.at(e).GetBinError(b));
                }
                int norm_bin_low_target = Hist_Target_Bkg_MSCW.at(e).FindBin(2*MSCW_cut_upper);
                int norm_bin_up_target = Hist_Target_Bkg_MSCW.at(e).FindBin(100*MSCW_cut_upper);
                double scale_target = Hist_Target_SR_MSCW.at(e).Integral(norm_bin_low_target,norm_bin_up_target)/Hist_Target_Bkg_MSCW.at(e).Integral(norm_bin_low_target,norm_bin_up_target);
                Hist_Target_Bkg_MSCW.at(e).Scale(scale_target);
        }
        TFile OutputFile("output/Histograms_"+TString(target)+"_"+TString(method)+".root","recreate"); 
        for (int e=0;e<N_energy_bins;e++) {
                Hist_CrabOff_SR_MSCW.at(e).Write();
                Hist_CrabOff_CR_MSCW.at(e).Write();
                Hist_CrabOff_Deconv_MSCW.at(e).Write();
                Hist_CrabOff_Bkg_MSCW.at(e).Write();
                Hist_Target_SR_MSCW.at(e).Write();
                Hist_Target_CR_MSCW.at(e).Write();
                Hist_Target_Bkg_MSCW.at(e).Write();
        }
        OutputFile.Close();

}
