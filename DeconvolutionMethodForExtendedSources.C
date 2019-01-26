
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
#include "TH2.h"
#include "TGraph.h"
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

#include "GetRunList.h"

char target[50] = "";
char Region[50] = "SR";
//char Region[50] = "VR";

double Elev_cut_lower = 0;
double Elev_cut_upper = 0;

double Target_Elev_cut_lower = 0;
double Target_Elev_cut_upper = 0;

double MSCW_cut_lower = -0.5;
double MSCW_cut_upper = 0.5;
double MSCL_signal_cut_lower = -0.5;
double MSCL_signal_cut_upper = 0.5;
double MSCL_control_cut_lower = 0.8;
double MSCL_control_cut_upper = 1.2;

//double MSCW_cut_lower = 0.8;
//double MSCW_cut_upper = 1.2;
//double MSCL_signal_cut_lower = -0.5+1.;
//double MSCL_signal_cut_upper = 0.5+1.;
//double MSCL_control_cut_lower = 2.0;
//double MSCL_control_cut_upper = 3.0;

double Norm_Lower = MSCW_cut_upper;
double Norm_Upper = 20.0;

string  filename;
double TelElevation = 0;
double TelAzimuth = 0;
double ErecS = 0;
double EChi2S = 0;
double MSCW = 0;
double MSCL = 0;
double SlantDepth = 0;
float EmissionHeight = 0;
float EmissionHeightChi2 = 0;
double theta2 = 0;
double ra_sky = 0;
double dec_sky = 0;

//const int N_energy_bins = 1;
//double energy_bins[N_energy_bins+1] = {1500,2000};
//const int N_energy_bins = 10;
//double energy_bins[N_energy_bins+1] = {150,200,250,300,400,600,1000,1500,2000,4000,10000};
//const int N_energy_bins = 7;
//double energy_bins[N_energy_bins+1] = {300,400,600,1000,1500,2000,4000,10000};
const int N_energy_bins = 4;
double energy_bins[N_energy_bins+1] = {1000,1500,2000,4000,10000};

int N_bins_for_deconv = 480;

bool FoV() {
    //if (theta2<4.0) return true;
    //if (theta2>2.0) return true;
    //if (theta2<0.5) return true;
    //if (theta2<1.0) return true;
    //if (theta2>2.0 && theta2<4.0) return true;
    if (theta2>1.0 && theta2<3.0) return true;
    return false;
}
bool RingFoV() {
    if (theta2<0.5) return false;
    if (FoV()) return false;
    //if (theta2<1.0) return false;
    //if (theta2>4.0) return false;
    return true;
}
bool QualitySelection() {
    if (EmissionHeightChi2/EmissionHeight>0.2) return false;
    return true;
}
bool SignalSelectionMSCL() {
    if (MSCW<MSCW_cut_lower) return false;
    if (MSCW>MSCW_cut_upper) return false;
    return true;
}
bool SignalSelectionMSCW() {
    if (MSCL<MSCL_signal_cut_lower) return false;
    if (MSCL>MSCL_signal_cut_upper) return false;
    return true;
}
bool ControlSelectionMSCW() {
    if (MSCL<MSCL_control_cut_lower) return false;
    if (MSCL>MSCL_control_cut_upper) return false;
    return true;
}
double background(Double_t *x, Double_t *par) {
    return (1./par[0])*(x[0]+par[1])*exp(-pow((x[0]+par[1])/par[0],2));
}
double GetChi2(TH1* Hist_SR, TH1* Hist_Bkg, bool includeSR) {
    double chi2_temp = 0.;
    double mean = Hist_SR->GetMean();
    double rms = Hist_SR->GetRMS();
    for (int i=0;i<Hist_SR->GetNbinsX();i++) {
        //if (Hist_SR->GetBinCenter(i+1)<mean-3.*rms) continue;
        //if (Hist_SR->GetBinCenter(i+1)>mean+3.*rms) continue;
        if (Hist_SR->GetBinCenter(i+1)<mean-4.*rms) continue;
        if (Hist_SR->GetBinCenter(i+1)>mean+4.*rms) continue;
        double bkg = Hist_Bkg->GetBinContent(i+1);
        double data = Hist_SR->GetBinContent(i+1);
        double bkg_err = Hist_Bkg->GetBinError(i+1);
        double data_err = Hist_SR->GetBinError(i+1);
        //if (data_err==0) data_err = 0.1;
        if (!includeSR && Hist_Bkg->GetBinCenter(i+1)>MSCW_cut_lower && Hist_Bkg->GetBinCenter(i+1)<MSCW_cut_upper) {
            continue;
        }
        if ((data_err*data_err+bkg_err*bkg_err)==0) continue;
        chi2_temp += pow(bkg-data,2)/(data_err*data_err+bkg_err*bkg_err);
    }
    chi2_temp = 1./chi2_temp;
    return chi2_temp;
}
double ShiftAndNormalize(TH1* Hist_SR, TH1* Hist_BkgTemp, TH1* Hist_Bkg, double shift_begin, bool doShift) {
    double shift_fit = 0;
    double scale_fit = 0;
    double chi2 = 0.;
    double mean = Hist_SR->GetMean();
    double rms = Hist_SR->GetRMS();
    int norm_bin_low = 0;
    int norm_bin_up = 0;
    double norm = Hist_BkgTemp->Integral();
    if (doShift) {
        for (int fit=0;fit<200;fit++) {
                double shift = shift_begin-10.0+double(fit)*0.1;
                for (int i=0;i<Hist_SR->GetNbinsX();i++) {
                        int b = Hist_SR->FindBin(Hist_SR->GetBinCenter(i+1)-shift);
                        Hist_Bkg->SetBinContent(i+1,Hist_BkgTemp->GetBinContent(b));
                        Hist_Bkg->SetBinError(i+1,Hist_BkgTemp->GetBinError(b));
                }
                Hist_Bkg->Scale(norm/Hist_Bkg->Integral());
                norm_bin_low = Hist_SR->FindBin(mean-8.*rms);
                norm_bin_up = Hist_SR->FindBin(MSCW_cut_lower);
                double SR_area1 = Hist_SR->Integral(norm_bin_low,norm_bin_up);
                double Bkg_area1 = Hist_Bkg->Integral(norm_bin_low,norm_bin_up);
                norm_bin_low = Hist_SR->FindBin(MSCW_cut_upper);
                norm_bin_up = Hist_SR->FindBin(mean+8.*rms);
                double SR_area2 = Hist_SR->Integral(norm_bin_low,norm_bin_up);
                double Bkg_area2 = Hist_Bkg->Integral(norm_bin_low,norm_bin_up);
                double scale_begin = (SR_area1+SR_area2)/(Bkg_area1+Bkg_area2);
                for (int ds=0;ds<20;ds++) {
                    double scale = scale_begin -0.2*scale_begin +0.4*scale_begin*double(ds)/20.;
                    Hist_Bkg->Scale(norm/Hist_Bkg->Integral());
                    Hist_Bkg->Scale(scale);
                    double chi2_temp = GetChi2(Hist_SR, Hist_Bkg,false);
                    if (chi2<chi2_temp) {
                        chi2 = chi2_temp;
                        shift_fit = shift;
                        scale_fit = scale;
                    } 
                }
                //    Hist_Bkg->Scale(scale_begin);
                //    double chi2_temp = GetChi2(Hist_SR, Hist_Bkg,false);
                //    if (chi2<chi2_temp) {
                //        chi2 = chi2_temp;
                //        shift_fit = shift;
                //        scale_fit = scale_begin;
                //    } 
        }
    }
    else {
        shift_fit = shift_begin;
    }
    for (int i=0;i<Hist_SR->GetNbinsX();i++) {
            int b = Hist_Bkg->FindBin(Hist_Bkg->GetBinCenter(i+1)-shift_fit);
            Hist_Bkg->SetBinContent(i+1,Hist_BkgTemp->GetBinContent(b));
            Hist_Bkg->SetBinError(i+1,Hist_BkgTemp->GetBinError(b));
    }
    //norm_bin_low = Hist_SR->FindBin(mean-8.*rms);
    //norm_bin_up = Hist_SR->FindBin(MSCW_cut_lower);
    //double SR_area1 = Hist_SR->Integral(norm_bin_low,norm_bin_up);
    //double Bkg_area1 = Hist_Bkg->Integral(norm_bin_low,norm_bin_up);
    //norm_bin_low = Hist_SR->FindBin(MSCW_cut_upper);
    //norm_bin_up = Hist_SR->FindBin(mean+8.*rms);
    //double SR_area2 = Hist_SR->Integral(norm_bin_low,norm_bin_up);
    //double Bkg_area2 = Hist_Bkg->Integral(norm_bin_low,norm_bin_up);
    //double scale = (SR_area1+SR_area2)/(Bkg_area1+Bkg_area2);
    Hist_Bkg->Scale(norm/Hist_Bkg->Integral());
    Hist_Bkg->Scale(scale_fit);
    return shift_fit;
}

void Deconvolution(TH1* Hist_source, TH1* Hist_response, TH1* Hist_Deconv, int n_iteration) {
        const int N_bins = Hist_source->GetNbinsX();
        Double_t *source = new Double_t[N_bins];
        Double_t *response = new Double_t[N_bins];
        for (int i=0;i<N_bins;i++) {
            source[i] = max(Hist_source->GetBinContent(i+1),0.);
            if (i<N_bins/2) response[i] = max(Hist_response->GetBinContent(i+1+N_bins/2),0.);
            else response[i] = 0;
        }
        TSpectrum sp;
        sp.Deconvolution(source,response,N_bins,n_iteration,1,100000); // new best option
        for (int i=0;i<N_bins;i++) {
            Hist_Deconv->SetBinContent(i+1,max(source[i],0.));
        }
}

double FindNIteration(TH1* Hist_SR, TH1* Hist_CR, TH1* Hist_Bkg, TH1* Hist_BkgTemp, TH1* Hist_Deconv, double rms, double mean, bool includeSR) {
    double n_iter_final = 1;
    double chi2_best = 0.;
    chi2_best = 0.;
    TF1 *func = new TF1("func", "gaus", -50, 50);
    func->SetParameters(10.,1.5,0.5);
    for (int delta_n_iter = 1;delta_n_iter<=100;delta_n_iter++) {
          int n_iter = delta_n_iter;
          double offset_begin = 0;
          double chi2 = 0;
          func->SetParameters(10.,mean,rms);
          Hist_Deconv->Reset();
          Hist_Deconv->FillRandom("func",Hist_SR->Integral()*100);
          Deconvolution(Hist_CR,Hist_Deconv,Hist_BkgTemp,n_iter);
          offset_begin = Hist_SR->GetMean()-Hist_BkgTemp->GetMean();
          offset_begin = ShiftAndNormalize(Hist_SR,Hist_BkgTemp,Hist_Bkg,offset_begin,true);
          chi2 = GetChi2(Hist_SR,Hist_Bkg,includeSR);
          if (chi2_best<chi2) {
              chi2_best = chi2;
              n_iter_final = n_iter;
          } 
    }
    return n_iter_final;
}
double FindRMS(TH1* Hist_SR, TH1* Hist_CR, TH1* Hist_Bkg, TH1* Hist_BkgTemp, TH1* Hist_Deconv, double rms_begin, double mean, double n_iter, bool includeSR) {
    double chi2_best = 0.;
    double rms_final = 0.;
    TF1 *func = new TF1("func", "gaus", -50, 50);
    func->SetParameters(10.,1.5,0.5);
    for (int n_rms = 0; n_rms<=50;n_rms++) {
        double offset_begin = 0;
        double chi2 = 0;
        double rms = rms_begin-0.5*rms_begin+double(n_rms)*1.5*rms_begin/50.;
        func->SetParameters(10.,mean,rms);
        Hist_Deconv->Reset();
        Hist_Deconv->FillRandom("func",Hist_SR->Integral()*100);
        Deconvolution(Hist_CR,Hist_Deconv,Hist_BkgTemp,n_iter);
        offset_begin = Hist_SR->GetMean()-Hist_BkgTemp->GetMean();
        offset_begin = ShiftAndNormalize(Hist_SR,Hist_BkgTemp,Hist_Bkg,offset_begin,true);
        chi2 = GetChi2(Hist_SR, Hist_Bkg,includeSR);
        if (chi2_best<chi2) {
            chi2_best = chi2;
            rms_final = rms;
        } 
    }
    return rms_final;
}
void Convolution(TH1D* Hist_source, TH1D* Hist_response, TH1D* Hist_Conv) {
        if (Hist_source->GetNbinsX()!=Hist_response->GetNbinsX()) 
            std::cout << "Convolving histograms with different bins!!!" << std::endl;
        if (Hist_source->Integral()==0) std::cout << "Hist_source->Integral()==0!!!" << std::endl;
        if (Hist_response->Integral()==0) std::cout << "Hist_response->Integral()==0!!!" << std::endl;
        for (int i=1;i<=Hist_source->GetNbinsX();i++) {
            double conv = 0;
            for (int j=1;j<=Hist_source->GetNbinsX();j++) {
                if (i>j)
                    conv += Hist_source->GetBinContent(j)*Hist_response->GetBinContent(i-j);
                else
                    conv += Hist_source->GetBinContent(j)*Hist_response->GetBinContent(i-j+N_bins_for_deconv);
            }
            Hist_Conv->SetBinContent(i,conv);
        }
}

void DeconvolutionMethodForExtendedSources(string target_data, double elev_lower, double elev_upper) {

        //TH1::SetDefaultSumw2();
        sprintf(target, "%s", target_data.c_str());
        Target_Elev_cut_lower = elev_lower;
        Target_Elev_cut_upper = elev_upper;
        Elev_cut_lower = elev_lower;
        Elev_cut_upper = elev_upper;

        if (TString(Region)=="VR") {
                MSCW_cut_lower = 1.5;
                MSCW_cut_upper = 2.5;
                MSCL_signal_cut_lower = 1.5;
                MSCL_signal_cut_upper = 2.5;
                MSCL_control_cut_lower = 2.5;
                MSCL_control_cut_upper = 3.5;
        }

        TRandom rnd;
        TH2D Hist_Dark_TelElevAzim("Hist_Dark_TelElevAzim","",90,0,90,360,0,360);
        TH2D Hist_Target_TelElevAzim("Hist_Target_TelElevAzim","",90,0,90,360,0,360);
        TH1D Hist_Target_ON_MSCW_Alpha("Hist_Target_ON_MSCW_Alpha","",100,0,10);
        TH1D Hist_Target_OFF_MSCW_Alpha("Hist_Target_OFF_MSCW_Alpha","",100,0,10);
        vector<TH1D> Hist_Dark_SR_ErecS;
        vector<TH1D> Hist_Dark_SR_MSCW;
        vector<TH1D> Hist_Dark_SR_MSCL;
        vector<TH1D> Hist_Dark_SRB_MSCW;
        vector<TH1D> Hist_Dark_CR_MSCW;
        vector<TH1D> Hist_Dark_ASR1_MSCW;
        vector<TH1D> Hist_Dark_ACR1_MSCW;
        vector<TH1D> Hist_Dark_ASR2_MSCW;
        vector<TH1D> Hist_Dark_ACR2_MSCW;
        vector<TH1D> Hist_Dark_BkgTemp_MSCW;
        vector<TH1D> Hist_Dark_Bkg_MSCW;
        vector<TH1D> Hist_Dark_ABkg1Temp_MSCW;
        vector<TH1D> Hist_Dark_ABkg1_MSCW;
        vector<TH1D> Hist_Dark_ABkg2Temp_MSCW;
        vector<TH1D> Hist_Dark_ABkg2_MSCW;
        vector<TH1D> Hist_Dark_Deconv_MSCW;
        vector<TH1D> Hist_Dark_TrueDeconv_MSCW;
        vector<TH1D> Hist_Dark_Elec_MSCW;
        vector<TH1D> Hist_Target_SR_ErecS;
        vector<TH1D> Hist_Target_SR_MSCW;
        vector<TH1D> Hist_Target_SR_MSCL;
        vector<TH1D> Hist_Target_SRB_MSCW;
        vector<TH1D> Hist_Target_CR_MSCW;
        vector<TH1D> Hist_Target_ASR1_MSCW;
        vector<TH1D> Hist_Target_ACR1_MSCW;
        vector<TH1D> Hist_Target_ASR2_MSCW;
        vector<TH1D> Hist_Target_ACR2_MSCW;
        vector<TH1D> Hist_Target_BkgTemp_MSCW;
        vector<TH1D> Hist_Target_Bkg_MSCW;
        vector<TH1D> Hist_Target_ABkg1Temp_MSCW;
        vector<TH1D> Hist_Target_ABkg1_MSCW;
        vector<TH1D> Hist_Target_ABkg2Temp_MSCW;
        vector<TH1D> Hist_Target_ABkg2_MSCW;
        vector<TH1D> Hist_Target_Deconv_MSCW;
        vector<TH1D> Hist_Target_TrueDeconv_MSCW;
        vector<TH1D> Hist_Target_Elec_MSCW;
        vector<TH1D> Hist_Target_Ring_MSCW;
        vector<TH1D> Hist_Target_Ring_MSCL;
        vector<TH1D> Hist_Target_TotalBkg_MSCW;
        vector<TH2D> Hist_Target_SR_SkyMap;
        vector<TH2D> Hist_Target_CR_SkyMap;
        vector<TH2D> Hist_Target_Bkg_SkyMap;
        vector<TH1D> Hist_Target_SR_Theta2;
        vector<TH1D> Hist_Target_CR_Theta2;
        vector<TH1D> Hist_Target_Bkg_Theta2;
        for (int e=0;e<N_energy_bins;e++) {
            char e_low[50];
            sprintf(e_low, "%i", int(energy_bins[e]));
            char e_up[50];
            sprintf(e_up, "%i", int(energy_bins[e+1]));
            N_bins_for_deconv = 480;
            Hist_Dark_SR_ErecS.push_back(TH1D("Hist_Dark_SR_ErecS_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_energy_bins,energy_bins));
            Hist_Dark_SR_MSCW.push_back(TH1D("Hist_Dark_SR_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-15,15));
            Hist_Dark_SR_MSCL.push_back(TH1D("Hist_Dark_SR_MSCL_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-15,15));
            Hist_Dark_SRB_MSCW.push_back(TH1D("Hist_Dark_SRB_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-15,15));
            Hist_Dark_CR_MSCW.push_back(TH1D("Hist_Dark_CR_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-15,15));
            Hist_Dark_BkgTemp_MSCW.push_back(TH1D("Hist_Dark_BkgTemp_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-15,15));
            Hist_Dark_Bkg_MSCW.push_back(TH1D("Hist_Dark_Bkg_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-15,15));
            Hist_Dark_Deconv_MSCW.push_back(TH1D("Hist_Dark_Deconv_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-15,15));
            Hist_Dark_TrueDeconv_MSCW.push_back(TH1D("Hist_Dark_TrueDeconv_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-15,15));
            Hist_Target_SR_ErecS.push_back(TH1D("Hist_Target_SR_ErecS_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_energy_bins,energy_bins));
            Hist_Target_SR_MSCW.push_back(TH1D("Hist_Target_SR_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-15,15));
            Hist_Target_SR_MSCL.push_back(TH1D("Hist_Target_SR_MSCL_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-15,15));
            Hist_Target_SRB_MSCW.push_back(TH1D("Hist_Target_SRB_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-15,15));
            Hist_Target_CR_MSCW.push_back(TH1D("Hist_Target_CR_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-15,15));
            Hist_Target_BkgTemp_MSCW.push_back(TH1D("Hist_Target_BkgTemp_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-15,15));
            Hist_Target_Bkg_MSCW.push_back(TH1D("Hist_Target_Bkg_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-15,15));
            Hist_Target_Deconv_MSCW.push_back(TH1D("Hist_Target_Deconv_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-15,15));
            Hist_Target_TrueDeconv_MSCW.push_back(TH1D("Hist_Target_TrueDeconv_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-15,15));
            Hist_Target_Ring_MSCW.push_back(TH1D("Hist_Target_Ring_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-15,15));
            Hist_Target_Ring_MSCL.push_back(TH1D("Hist_Target_Ring_MSCL_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-15,15));
        }

        std::cout << "Getting dark runs... " << std::endl;
        vector<int> Dark_runlist = GetRunList("BrandonValidation");
        //vector<int> Dark_runlist = GetRunList("Crab");
        //vector<int> Dark_runlist = GetRunList("PKS1424");
        for (int run=0;run<Dark_runlist.size();run++) {
                char run_number[50];
                sprintf(run_number, "%i", int(Dark_runlist[run]));
                //filename = TString("$VERITAS_USER_DATA_DIR/Crab_V6_Moderate-TMVA-BDT.RB."+TString(run_number)+".root");
                //filename = TString("$VERITAS_USER_DATA_DIR/PKS1424_V6_Moderate-TMVA-BDT.RB."+TString(run_number)+".root");
                filename = TString("/veritas/userspace/brandon/VERITAS/Background/anasum/"+TString(run_number)+".anasum.root");
                TFile*  input_file = TFile::Open(filename.c_str());
                TTree* pointing_tree = (TTree*) input_file->Get("run_"+TString(run_number)+"/stereo/pointingDataReduced");
                pointing_tree->SetBranchAddress("TelElevation",&TelElevation);
                pointing_tree->SetBranchAddress("TelAzimuth",&TelAzimuth);
                pointing_tree->GetEntry(0);
                if (TelElevation<Elev_cut_lower) continue;
                if (TelElevation>Elev_cut_upper) continue;
                //if (TelElevation<55.) continue;
                //if (TelElevation>85.) continue;
                for (int entry=0;entry<pointing_tree->GetEntries();entry++) {
                        pointing_tree->GetEntry(entry);
                        Hist_Dark_TelElevAzim.Fill(TelElevation,TelAzimuth);
                }

                //TTree* Dark_tree = (TTree*) input_file->Get("run_"+TString(run_number)+"/stereo/data_off");
                TTree* Dark_tree = (TTree*) input_file->Get("run_"+TString(run_number)+"/stereo/data_on");
                Dark_tree->SetBranchAddress("theta2",&theta2);
                Dark_tree->SetBranchAddress("ErecS",&ErecS);
                Dark_tree->SetBranchAddress("EChi2S",&EChi2S);
                Dark_tree->SetBranchAddress("MSCW",&MSCW);
                Dark_tree->SetBranchAddress("MSCL",&MSCL);
                //Dark_tree->SetBranchAddress("MWR",&MSCW);
                //Dark_tree->SetBranchAddress("MLR",&MSCL);
                for (int entry=0;entry<Dark_tree->GetEntries();entry++) {
                        Dark_tree->GetEntry(entry);
                        int e = Hist_Dark_SR_ErecS.at(0).FindBin(ErecS*1000.)-1;
                        if (e>=N_energy_bins) continue;
                        if (e<0) continue;
                        if (!QualitySelection()) continue;
                        if (FoV()) {
                            if (SignalSelectionMSCW()) Hist_Dark_SR_MSCW.at(e).Fill(MSCW);
                            if (SignalSelectionMSCL()) Hist_Dark_SR_MSCL.at(e).Fill(MSCL);
                            if (ControlSelectionMSCW()) Hist_Dark_CR_MSCW.at(e).Fill(MSCW);
                        }
                }
                input_file->Close();

        }

        std::cout << "Getting target runs... " << std::endl;
        vector<int> Target_runlist = GetRunList(target);
        char observation[50];
        sprintf(observation, "%s", target);
        if (TString(target)=="2ndCrab") sprintf(observation, "%s", "Crab");
        for (int run=0;run<Target_runlist.size();run++) {
                char run_number[50];
                sprintf(run_number, "%i", int(Target_runlist[run]));
                filename = TString("$VERITAS_USER_DATA_DIR/"+TString(observation)+"_V6_Moderate-TMVA-BDT.RB."+TString(run_number)+".root");
                if (TString(observation)=="Coma" || TString(observation)=="BrandonValidation") {
                  filename = TString("/veritas/userspace/brandon/VERITAS/Background/anasum/"+TString(run_number)+".anasum.root");
                }
                TFile*  input_file = TFile::Open(filename.c_str());
                TTree* pointing_tree = (TTree*) input_file->Get("run_"+TString(run_number)+"/stereo/pointingDataReduced");
                pointing_tree->SetBranchAddress("TelElevation",&TelElevation);
                pointing_tree->SetBranchAddress("TelAzimuth",&TelAzimuth);
                pointing_tree->GetEntry(0);
                if (TelElevation<Target_Elev_cut_lower) continue;
                if (TelElevation>Target_Elev_cut_upper) continue;
                for (int entry=0;entry<pointing_tree->GetEntries();entry++) {
                        pointing_tree->GetEntry(entry);
                        Hist_Target_TelElevAzim.Fill(TelElevation,TelAzimuth);
                }

                TString root_file = "run_"+TString(run_number)+"/stereo/data_on";
                TTree* Target_tree = (TTree*) input_file->Get(root_file);
                Target_tree->SetBranchAddress("theta2",&theta2);
                Target_tree->SetBranchAddress("ErecS",&ErecS);
                Target_tree->SetBranchAddress("EChi2S",&EChi2S);
                Target_tree->SetBranchAddress("MSCW",&MSCW);
                Target_tree->SetBranchAddress("MSCL",&MSCL);
                //Target_tree->SetBranchAddress("MWR",&MSCW);
                //Target_tree->SetBranchAddress("MLR",&MSCL);
                for (int entry=0;entry<Target_tree->GetEntries();entry++) {
                        Target_tree->GetEntry(entry);
                        int e = Hist_Target_SR_ErecS.at(0).FindBin(ErecS*1000.)-1;
                        if (e>=N_energy_bins) continue;
                        if (e<0) continue;
                        if (!QualitySelection()) continue;
                        if (MSCL<MSCL_signal_cut_upper) {
                            if (theta2<0.1) {
                                Hist_Target_ON_MSCW_Alpha.Fill(MSCW);
                            }
                            else {
                                Hist_Target_OFF_MSCW_Alpha.Fill(MSCW);
                            }
                        }
                        if (FoV()) {
                                if (SignalSelectionMSCW()) Hist_Target_SR_MSCW.at(e).Fill(MSCW);
                                if (SignalSelectionMSCL()) Hist_Target_SR_MSCL.at(e).Fill(MSCL);
                                if (ControlSelectionMSCW()) Hist_Target_CR_MSCW.at(e).Fill(MSCW);
                        }
                        else if (RingFoV()) {
                                if (SignalSelectionMSCW()) Hist_Target_Ring_MSCW.at(e).Fill(MSCW);
                                if (SignalSelectionMSCL()) Hist_Target_Ring_MSCL.at(e).Fill(MSCL);
                        }
                }
                input_file->Close();

        }

        std::cout << "Rebin histograms... " << std::endl;
        for (int e=0;e<N_energy_bins;e++) {
            int n_rebins = 0;
            double rms = Hist_Target_SR_MSCW.at(e).GetRMS();
            for (int n_rebins=0;n_rebins<10;n_rebins++) {
                if (Hist_Target_SR_MSCW.at(e).GetNbinsX()<=120) break;
                double hist_integral = 0;
                double hist_error = 0;
                for (int b=0;b<Hist_Target_SR_MSCW.at(e).GetNbinsX();b++) {
                    hist_integral += Hist_Target_SR_MSCW.at(e).GetBinContent(b+1);
                    hist_error += Hist_Target_SR_MSCW.at(e).GetBinError(b+1);
                }
                double bin_width = Hist_Target_SR_MSCW.at(e).GetBinWidth(1);
                if (hist_error/hist_integral>0.1) {
                //if (bin_width/rms<0.2 || hist_error/hist_integral>0.2) {
                        Hist_Dark_SR_MSCW.at(e).Rebin(2);
                        Hist_Dark_SR_MSCL.at(e).Rebin(2);
                        Hist_Dark_SRB_MSCW.at(e).Rebin(2);
                        Hist_Dark_CR_MSCW.at(e).Rebin(2);
                        Hist_Dark_Bkg_MSCW.at(e).Rebin(2);
                        Hist_Dark_BkgTemp_MSCW.at(e).Rebin(2);
                        Hist_Dark_Deconv_MSCW.at(e).Rebin(2);
                        Hist_Dark_TrueDeconv_MSCW.at(e).Rebin(2);
                        Hist_Target_SR_MSCW.at(e).Rebin(2);
                        Hist_Target_SR_MSCL.at(e).Rebin(2);
                        Hist_Target_SRB_MSCW.at(e).Rebin(2);
                        Hist_Target_CR_MSCW.at(e).Rebin(2);
                        Hist_Target_Bkg_MSCW.at(e).Rebin(2);
                        Hist_Target_BkgTemp_MSCW.at(e).Rebin(2);
                        Hist_Target_Deconv_MSCW.at(e).Rebin(2);
                        Hist_Target_TrueDeconv_MSCW.at(e).Rebin(2);
                        Hist_Target_Ring_MSCW.at(e).Rebin(2);
                        Hist_Target_Ring_MSCL.at(e).Rebin(2);
                }
                else break;
            }
        }

        vector<int> N_iter;
        vector<double> N_rms;
        vector<double> N_mean;
        vector<double> N_shift;
        for (int e=0;e<N_energy_bins;e++) {
                N_iter.push_back(30);
                N_rms.push_back(1);
                N_mean.push_back(1);
                N_shift.push_back(0);
        }

        TF1 *myfunc = new TF1("myfunc", "gaus", -50, 50);
        myfunc->SetParameters(10.,1.5,0.5);

        std::cout << "Dark run deconvolution... " << std::endl;
        // Dark deconvolution method
        for (int e=0;e<N_energy_bins;e++) {
            std::cout << "=================================================================" << std::endl;
            std::cout << "Dark, e " << energy_bins[e] << ", N_iter.at(e) = " << N_iter.at(e) << std::endl;


            Deconvolution(&Hist_Dark_CR_MSCW.at(e),&Hist_Dark_SR_MSCW.at(e),&Hist_Dark_Deconv_MSCW.at(e),N_iter.at(e));
            double rms_begin = Hist_Dark_Deconv_MSCW.at(e).GetRMS();
            double mean_begin = Hist_Dark_Deconv_MSCW.at(e).GetMean();
            std::cout << "Dark, e " << energy_bins[e] << ", rms_begin = " << rms_begin << std::endl;

            N_iter.at(e) = FindNIteration(&Hist_Dark_SR_MSCW.at(e),&Hist_Dark_CR_MSCW.at(e),&Hist_Dark_Bkg_MSCW.at(e),&Hist_Dark_BkgTemp_MSCW.at(e),&Hist_Dark_Deconv_MSCW.at(e),rms_begin,mean_begin,true);
            std::cout << "Dark, e " << energy_bins[e] << ", N_iter.at(e) = " << N_iter.at(e) << std::endl;
            N_rms.at(e) = FindRMS(&Hist_Dark_SR_MSCW.at(e),&Hist_Dark_CR_MSCW.at(e),&Hist_Dark_Bkg_MSCW.at(e),&Hist_Dark_BkgTemp_MSCW.at(e),&Hist_Dark_Deconv_MSCW.at(e),rms_begin,mean_begin,N_iter.at(e),true);
            std::cout << "Dark, e " << energy_bins[e] << ", N_rms.at(e) = " << N_rms.at(e) << std::endl;

            N_iter.at(e) = FindNIteration(&Hist_Dark_SR_MSCW.at(e),&Hist_Dark_CR_MSCW.at(e),&Hist_Dark_Bkg_MSCW.at(e),&Hist_Dark_BkgTemp_MSCW.at(e),&Hist_Dark_Deconv_MSCW.at(e),N_rms.at(e),mean_begin,true);
            std::cout << "Dark, e " << energy_bins[e] << ", N_iter.at(e) = " << N_iter.at(e) << std::endl;
            N_rms.at(e) = FindRMS(&Hist_Dark_SR_MSCW.at(e),&Hist_Dark_CR_MSCW.at(e),&Hist_Dark_Bkg_MSCW.at(e),&Hist_Dark_BkgTemp_MSCW.at(e),&Hist_Dark_Deconv_MSCW.at(e),N_rms.at(e),mean_begin,N_iter.at(e),true);
            std::cout << "Dark, e " << energy_bins[e] << ", N_rms.at(e) = " << N_rms.at(e) << std::endl;

            N_iter.at(e) = FindNIteration(&Hist_Dark_SR_MSCW.at(e),&Hist_Dark_CR_MSCW.at(e),&Hist_Dark_Bkg_MSCW.at(e),&Hist_Dark_BkgTemp_MSCW.at(e),&Hist_Dark_Deconv_MSCW.at(e),N_rms.at(e),mean_begin,true);
            std::cout << "Dark, e " << energy_bins[e] << ", N_iter.at(e) = " << N_iter.at(e) << std::endl;
            N_rms.at(e) = FindRMS(&Hist_Dark_SR_MSCW.at(e),&Hist_Dark_CR_MSCW.at(e),&Hist_Dark_Bkg_MSCW.at(e),&Hist_Dark_BkgTemp_MSCW.at(e),&Hist_Dark_Deconv_MSCW.at(e),N_rms.at(e),mean_begin,N_iter.at(e),true);
            std::cout << "Dark, e " << energy_bins[e] << ", N_rms.at(e) = " << N_rms.at(e) << std::endl;

            N_iter.at(e) = FindNIteration(&Hist_Dark_SR_MSCW.at(e),&Hist_Dark_CR_MSCW.at(e),&Hist_Dark_Bkg_MSCW.at(e),&Hist_Dark_BkgTemp_MSCW.at(e),&Hist_Dark_Deconv_MSCW.at(e),N_rms.at(e),mean_begin,true);
            std::cout << "Dark, e " << energy_bins[e] << ", N_iter.at(e) = " << N_iter.at(e) << std::endl;
            N_rms.at(e) = FindRMS(&Hist_Dark_SR_MSCW.at(e),&Hist_Dark_CR_MSCW.at(e),&Hist_Dark_Bkg_MSCW.at(e),&Hist_Dark_BkgTemp_MSCW.at(e),&Hist_Dark_Deconv_MSCW.at(e),N_rms.at(e),mean_begin,N_iter.at(e),true);
            std::cout << "Dark, e " << energy_bins[e] << ", N_rms.at(e) = " << N_rms.at(e) << std::endl;

            N_iter.at(e) = FindNIteration(&Hist_Dark_SR_MSCW.at(e),&Hist_Dark_CR_MSCW.at(e),&Hist_Dark_Bkg_MSCW.at(e),&Hist_Dark_BkgTemp_MSCW.at(e),&Hist_Dark_Deconv_MSCW.at(e),N_rms.at(e),mean_begin,true);
            std::cout << "Dark, e " << energy_bins[e] << ", N_iter.at(e) = " << N_iter.at(e) << std::endl;
            N_rms.at(e) = FindRMS(&Hist_Dark_SR_MSCW.at(e),&Hist_Dark_CR_MSCW.at(e),&Hist_Dark_Bkg_MSCW.at(e),&Hist_Dark_BkgTemp_MSCW.at(e),&Hist_Dark_Deconv_MSCW.at(e),N_rms.at(e),mean_begin,N_iter.at(e),true);
            std::cout << "Dark, e " << energy_bins[e] << ", N_rms.at(e) = " << N_rms.at(e) << std::endl;

            N_mean.at(e) = mean_begin;

            double rms = N_rms.at(e);
            double mean = N_mean.at(e);
            myfunc->SetParameters(10.,mean,rms);
            Hist_Dark_Deconv_MSCW.at(e).Reset();
            Hist_Dark_Deconv_MSCW.at(e).FillRandom("myfunc",Hist_Dark_SR_MSCW.at(e).Integral()*100);
            Deconvolution(&Hist_Dark_CR_MSCW.at(e),&Hist_Dark_Deconv_MSCW.at(e),&Hist_Dark_BkgTemp_MSCW.at(e),N_iter.at(e));
            Deconvolution(&Hist_Dark_CR_MSCW.at(e),&Hist_Dark_SR_MSCW.at(e),&Hist_Dark_TrueDeconv_MSCW.at(e),N_iter.at(e));
            std::cout << "Dark, e " << energy_bins[e] << ", final mean = " << Hist_Dark_Deconv_MSCW.at(e).GetMean() << std::endl;
            std::cout << "Dark, e " << energy_bins[e] << ", final rms = " << Hist_Dark_Deconv_MSCW.at(e).GetRMS() << std::endl;
            std::cout << "Dark, e " << energy_bins[e] << ", true Mean = " << Hist_Dark_TrueDeconv_MSCW.at(e).GetMean() << std::endl;
            std::cout << "Dark, e " << energy_bins[e] << ", true RMS = " << Hist_Dark_TrueDeconv_MSCW.at(e).GetRMS() << std::endl;
            double offset_begin = 0;
            offset_begin = Hist_Dark_SR_MSCW.at(e).GetMean()-Hist_Dark_BkgTemp_MSCW.at(e).GetMean();
            std::cout << "offset_begin = " << offset_begin << std::endl;
            std::cout << "Dark, e " << energy_bins[e] << ", offset_begin = " << offset_begin << std::endl;
            offset_begin = ShiftAndNormalize(&Hist_Dark_SR_MSCW.at(e),&Hist_Dark_BkgTemp_MSCW.at(e),&Hist_Dark_Bkg_MSCW.at(e),offset_begin,true);
            N_shift.at(e) = offset_begin;
            std::cout << "offset_begin = " << offset_begin << std::endl;
        }

        std::cout << "Target run deconvolution... " << std::endl;
        // Target deconvolution method
        for (int e=0;e<N_energy_bins;e++) {
            std::cout << "=================================================================" << std::endl;
            std::cout << "Target, e " << energy_bins[e] << ", N_iter.at(e) = " << N_iter.at(e) << std::endl;


            double rms_begin = N_rms.at(e);
            double mean_begin = N_mean.at(e);

            //N_iter.at(e) = FindNIteration(&Hist_Target_SR_MSCW.at(e),&Hist_Target_CR_MSCW.at(e),&Hist_Target_Bkg_MSCW.at(e),&Hist_Target_BkgTemp_MSCW.at(e),&Hist_Target_Deconv_MSCW.at(e),rms_begin,mean_begin,false);
            std::cout << "Target, e " << energy_bins[e] << ", N_iter.at(e) = " << N_iter.at(e) << std::endl;
            N_rms.at(e) = FindRMS(&Hist_Target_SR_MSCW.at(e),&Hist_Target_CR_MSCW.at(e),&Hist_Target_Bkg_MSCW.at(e),&Hist_Target_BkgTemp_MSCW.at(e),&Hist_Target_Deconv_MSCW.at(e),rms_begin,mean_begin,N_iter.at(e),false);
            std::cout << "Target, e " << energy_bins[e] << ", N_rms.at(e) = " << N_rms.at(e) << std::endl;

            //N_iter.at(e) = FindNIteration(&Hist_Target_SR_MSCW.at(e),&Hist_Target_CR_MSCW.at(e),&Hist_Target_Bkg_MSCW.at(e),&Hist_Target_BkgTemp_MSCW.at(e),&Hist_Target_Deconv_MSCW.at(e),N_rms.at(e),mean_begin,false);
            std::cout << "Target, e " << energy_bins[e] << ", N_iter.at(e) = " << N_iter.at(e) << std::endl;
            N_rms.at(e) = FindRMS(&Hist_Target_SR_MSCW.at(e),&Hist_Target_CR_MSCW.at(e),&Hist_Target_Bkg_MSCW.at(e),&Hist_Target_BkgTemp_MSCW.at(e),&Hist_Target_Deconv_MSCW.at(e),N_rms.at(e),mean_begin,N_iter.at(e),false);
            std::cout << "Target, e " << energy_bins[e] << ", N_rms.at(e) = " << N_rms.at(e) << std::endl;

            //N_iter.at(e) = FindNIteration(&Hist_Target_SR_MSCW.at(e),&Hist_Target_CR_MSCW.at(e),&Hist_Target_Bkg_MSCW.at(e),&Hist_Target_BkgTemp_MSCW.at(e),&Hist_Target_Deconv_MSCW.at(e),N_rms.at(e),mean_begin,false);
            std::cout << "Target, e " << energy_bins[e] << ", N_iter.at(e) = " << N_iter.at(e) << std::endl;
            N_rms.at(e) = FindRMS(&Hist_Target_SR_MSCW.at(e),&Hist_Target_CR_MSCW.at(e),&Hist_Target_Bkg_MSCW.at(e),&Hist_Target_BkgTemp_MSCW.at(e),&Hist_Target_Deconv_MSCW.at(e),N_rms.at(e),mean_begin,N_iter.at(e),false);
            std::cout << "Target, e " << energy_bins[e] << ", N_rms.at(e) = " << N_rms.at(e) << std::endl;

            //N_iter.at(e) = FindNIteration(&Hist_Target_SR_MSCW.at(e),&Hist_Target_CR_MSCW.at(e),&Hist_Target_Bkg_MSCW.at(e),&Hist_Target_BkgTemp_MSCW.at(e),&Hist_Target_Deconv_MSCW.at(e),N_rms.at(e),mean_begin,false);
            std::cout << "Target, e " << energy_bins[e] << ", N_iter.at(e) = " << N_iter.at(e) << std::endl;
            N_rms.at(e) = FindRMS(&Hist_Target_SR_MSCW.at(e),&Hist_Target_CR_MSCW.at(e),&Hist_Target_Bkg_MSCW.at(e),&Hist_Target_BkgTemp_MSCW.at(e),&Hist_Target_Deconv_MSCW.at(e),N_rms.at(e),mean_begin,N_iter.at(e),false);
            std::cout << "Target, e " << energy_bins[e] << ", N_rms.at(e) = " << N_rms.at(e) << std::endl;

            //N_iter.at(e) = FindNIteration(&Hist_Target_SR_MSCW.at(e),&Hist_Target_CR_MSCW.at(e),&Hist_Target_Bkg_MSCW.at(e),&Hist_Target_BkgTemp_MSCW.at(e),&Hist_Target_Deconv_MSCW.at(e),N_rms.at(e),mean_begin,false);
            std::cout << "Target, e " << energy_bins[e] << ", N_iter.at(e) = " << N_iter.at(e) << std::endl;
            N_rms.at(e) = FindRMS(&Hist_Target_SR_MSCW.at(e),&Hist_Target_CR_MSCW.at(e),&Hist_Target_Bkg_MSCW.at(e),&Hist_Target_BkgTemp_MSCW.at(e),&Hist_Target_Deconv_MSCW.at(e),N_rms.at(e),mean_begin,N_iter.at(e),false);
            std::cout << "Target, e " << energy_bins[e] << ", N_rms.at(e) = " << N_rms.at(e) << std::endl;

            N_mean.at(e) = mean_begin;

            double rms = N_rms.at(e);
            double mean = N_mean.at(e);
            myfunc->SetParameters(10.,mean,rms);
            Hist_Target_Deconv_MSCW.at(e).Reset();
            Hist_Target_Deconv_MSCW.at(e).FillRandom("myfunc",Hist_Target_SR_MSCW.at(e).Integral()*100);
            Deconvolution(&Hist_Target_CR_MSCW.at(e),&Hist_Target_Deconv_MSCW.at(e),&Hist_Target_BkgTemp_MSCW.at(e),N_iter.at(e));
            Deconvolution(&Hist_Target_CR_MSCW.at(e),&Hist_Target_SR_MSCW.at(e),&Hist_Target_TrueDeconv_MSCW.at(e),N_iter.at(e));
            std::cout << "Target, e " << energy_bins[e] << ", final mean = " << Hist_Target_Deconv_MSCW.at(e).GetMean() << std::endl;
            std::cout << "Target, e " << energy_bins[e] << ", final rms = " << Hist_Target_Deconv_MSCW.at(e).GetRMS() << std::endl;
            std::cout << "Target, e " << energy_bins[e] << ", true Mean = " << Hist_Target_TrueDeconv_MSCW.at(e).GetMean() << std::endl;
            std::cout << "Target, e " << energy_bins[e] << ", true RMS = " << Hist_Target_TrueDeconv_MSCW.at(e).GetRMS() << std::endl;
            double offset_begin = 0;
            offset_begin = Hist_Target_SR_MSCW.at(e).GetMean()-Hist_Target_BkgTemp_MSCW.at(e).GetMean();
            std::cout << "offset_begin = " << offset_begin << std::endl;
            std::cout << "Target, e " << energy_bins[e] << ", offset_begin = " << offset_begin << std::endl;
            offset_begin = ShiftAndNormalize(&Hist_Target_SR_MSCW.at(e),&Hist_Target_BkgTemp_MSCW.at(e),&Hist_Target_Bkg_MSCW.at(e),offset_begin,true);
            N_shift.at(e) = offset_begin;
            std::cout << "offset_begin = " << offset_begin << std::endl;
        }


        // simple dark template method
        std::cout << "Simple dark template... " << std::endl;
        for (int e=0;e<N_energy_bins;e++) {
                int norm_bin_low_target = Hist_Target_Bkg_MSCW.at(e).FindBin(Norm_Lower);
                int norm_bin_up_target = Hist_Target_Bkg_MSCW.at(e).FindBin(Norm_Upper);
                double scale_target = Hist_Target_SR_MSCW.at(e).Integral(norm_bin_low_target,norm_bin_up_target)/Hist_Dark_SR_MSCW.at(e).Integral(norm_bin_low_target,norm_bin_up_target);
                Hist_Dark_SR_MSCW.at(e).Scale(scale_target);
                Hist_Dark_SR_MSCL.at(e).Scale(scale_target);
                Hist_Dark_SRB_MSCW.at(e).Scale(scale_target);
                Hist_Dark_Bkg_MSCW.at(e).Scale(scale_target);
        }

        // Get ring bkg
        for (int e=0;e<N_energy_bins;e++) {
                int norm_bin_low_ring = Hist_Target_Ring_MSCW.at(e).FindBin(Norm_Lower);
                int norm_bin_up_ring = Hist_Target_Ring_MSCW.at(e).FindBin(Norm_Upper);
                double scale_ring = Hist_Target_SR_MSCW.at(e).Integral(norm_bin_low_ring,norm_bin_up_ring)/Hist_Target_Ring_MSCW.at(e).Integral(norm_bin_low_ring,norm_bin_up_ring);
                if (!(scale_ring>0)) scale_ring = 0;
                if (!(Hist_Target_Ring_MSCW.at(e).Integral(norm_bin_low_ring,norm_bin_up_ring)>0)) scale_ring = 0;
                Hist_Target_Ring_MSCW.at(e).Scale(scale_ring);
                Hist_Target_Ring_MSCL.at(e).Scale(scale_ring);
        }


        int norm_bin_low_ring = Hist_Target_ON_MSCW_Alpha.FindBin(Norm_Lower);
        int norm_bin_up_ring = Hist_Target_ON_MSCW_Alpha.FindBin(Norm_Upper);
        double scale_ring = Hist_Target_ON_MSCW_Alpha.Integral(norm_bin_low_ring,norm_bin_up_ring)/Hist_Target_OFF_MSCW_Alpha.Integral(norm_bin_low_ring,norm_bin_up_ring);

        TFile OutputFile("output/Deconvolution_"+TString(target)+"_Elev"+std::to_string(int(Target_Elev_cut_lower))+"to"+std::to_string(int(Target_Elev_cut_upper))+"_"+TString(Region)+".root","recreate"); 
        TTree InfoTree("InfoTree","info tree");
        InfoTree.Branch("MSCW_cut_lower",&MSCW_cut_lower,"MSCW_cut_lower/D");
        InfoTree.Branch("MSCW_cut_upper",&MSCW_cut_upper,"MSCW_cut_upper/D");
        InfoTree.Branch("MSCL_signal_cut_lower",&MSCL_signal_cut_lower,"MSCL_signal_cut_lower/D");
        InfoTree.Branch("MSCL_signal_cut_upper",&MSCL_signal_cut_upper,"MSCL_signal_cut_upper/D");
        InfoTree.Branch("MSCL_control_cut_lower",&MSCL_control_cut_lower,"MSCL_control_cut_lower/D");
        InfoTree.Branch("MSCL_control_cut_upper",&MSCL_control_cut_upper,"MSCL_control_cut_upper/D");
        InfoTree.Branch("Elev_cut_lower",&Elev_cut_lower,"Elev_cut_lower/D");
        InfoTree.Branch("Elev_cut_upper",&Elev_cut_upper,"Elev_cut_upper/D");
        InfoTree.Branch("Target_Elev_cut_lower",&Target_Elev_cut_lower,"Target_Elev_cut_lower/D");
        InfoTree.Branch("Target_Elev_cut_upper",&Target_Elev_cut_upper,"Target_Elev_cut_upper/D");
        InfoTree.Fill();
        InfoTree.Write();
        Hist_Target_TelElevAzim.Write();
        Hist_Dark_TelElevAzim.Write();
        for (int e=0;e<N_energy_bins;e++) {
                Hist_Dark_SR_MSCW.at(e).Write();
                Hist_Dark_SR_MSCL.at(e).Write();
                Hist_Dark_SRB_MSCW.at(e).Write();
                Hist_Dark_CR_MSCW.at(e).Write();
                Hist_Dark_Bkg_MSCW.at(e).Write();
                Hist_Dark_Deconv_MSCW.at(e).Write();
                Hist_Dark_TrueDeconv_MSCW.at(e).Write();
                Hist_Target_SR_MSCW.at(e).Write();
                Hist_Target_SR_MSCL.at(e).Write();
                Hist_Target_SRB_MSCW.at(e).Write();
                Hist_Target_CR_MSCW.at(e).Write();
                Hist_Target_Bkg_MSCW.at(e).Write();
                Hist_Target_Deconv_MSCW.at(e).Write();
                Hist_Target_TrueDeconv_MSCW.at(e).Write();
                Hist_Target_Ring_MSCW.at(e).Write();
                Hist_Target_Ring_MSCL.at(e).Write();
        }
        OutputFile.Close();

}
