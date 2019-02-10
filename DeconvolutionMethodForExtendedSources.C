
#include <cstdlib>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <map>
#include <string>
#include <iomanip>
#include <stdio.h>
#include <utility>

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

double Azim_cut_lower = 0;
double Azim_cut_upper = 0;
double Target_Azim_cut_lower = 0;
double Target_Azim_cut_upper = 0;

double Theta2_cut_lower = 0;
double Theta2_cut_upper = 0;

// SR
double MSCW_cut_lower = -1.0;
double MSCW_cut_upper = 1.0;
double MSCL_signal_cut_lower = -0.5;
double MSCL_signal_cut_upper = 0.5;
double MSCL_control1_cut_lower = 0.5;
double MSCL_control1_cut_upper = 1.0;
double MSCL_control2_cut_lower = 1.0;
double MSCL_control2_cut_upper = 1.55;

// VR
//double MSCW_cut_lower = -0.5;
//double MSCW_cut_upper = 1.0;
//double MSCL_signal_cut_lower = 0.5;
//double MSCL_signal_cut_upper = 1.0;
//double MSCL_control_cut_lower = 1.0;
//double MSCL_control_cut_upper = 1.5;

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
double TelRAJ2000 = 0;
double TelDecJ2000 = 0;
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

const int N_energy_bins = 11;
double energy_bins[N_energy_bins+1] =     {400,500,600,700,800,1000,1200,1500,2000,3000,5000,10000};
int number_runs_included[N_energy_bins] = {4  ,4  ,8  ,8  ,32 ,32  ,32  ,32  ,32  ,64  ,64};
//const int N_energy_bins = 8;
//double energy_bins[N_energy_bins+1] =     {600,800,1000,1200,1500,2000,3000,5000,10000};
//int number_runs_included[N_energy_bins] = {8  ,16 ,16  ,16  ,16  ,32  ,64  ,64};
//const int N_energy_bins = 2;
//double energy_bins[N_energy_bins+1] =     {3000,5000,10000};
//int number_runs_included[N_energy_bins] = {16  ,64};

int N_bins_for_deconv = 480;

bool FoV() {
    if (theta2>10.0) return false;
    if (theta2<Theta2_cut_lower) return false;
    if (theta2>Theta2_cut_upper) return false;
    return true;
}
bool RingFoV() {
    if (theta2<0.5) return false;
    if (FoV()) return false;
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
bool ControlSelection1MSCW() {
    if (MSCL<MSCL_control1_cut_lower) return false;
    if (MSCL>MSCL_control1_cut_upper) return false;
    return true;
}
bool ControlSelection2MSCW() {
    MSCL_control2_cut_upper = 1.55;
    if (ErecS*1000.>=1500.) MSCL_control2_cut_upper = 1.55;
    if (MSCL<MSCL_control2_cut_lower) return false;
    if (MSCL>MSCL_control2_cut_upper) return false;
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
        //if (Hist_SR->GetBinCenter(i+1)<mean-5.*rms) continue;
        //if (Hist_SR->GetBinCenter(i+1)>mean+5.*rms) continue;
        //if (Hist_SR->GetBinCenter(i+1)<mean-6.*rms) continue;
        //if (Hist_SR->GetBinCenter(i+1)>mean+6.*rms) continue;
        //if (Hist_SR->GetBinCenter(i+1)<mean-8.*rms) continue;
        //if (Hist_SR->GetBinCenter(i+1)>mean+8.*rms) continue;
        double bkg = Hist_Bkg->GetBinContent(i+1);
        double data = Hist_SR->GetBinContent(i+1);
        double bkg_err = Hist_Bkg->GetBinError(i+1);
        double data_err = Hist_SR->GetBinError(i+1);
        //if (Hist_Bkg->GetBinCenter(i+1)<MSCW_cut_lower)
        //{
        //    data = 0.;
        //    data_err = 0.;
        //}
        if (!includeSR && Hist_Bkg->GetBinCenter(i+1)>MSCW_cut_lower && Hist_Bkg->GetBinCenter(i+1)<MSCW_cut_upper) {
            continue;
        }
        if ((data_err*data_err+bkg_err*bkg_err)==0) continue;
        chi2_temp += pow(bkg-data,2)/(data_err*data_err+bkg_err*bkg_err);
    }
    chi2_temp = 1./chi2_temp;
    return chi2_temp;
}
double ShiftAndNormalize(TH1* Hist_SR, TH1* Hist_BkgTemp, TH1* Hist_Bkg, double shift_begin, bool doShift, bool includeSR) {
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
                    double chi2_temp = GetChi2(Hist_SR, Hist_Bkg,includeSR);
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

double FindNIteration(TH1* Hist_SR, TH1* Hist_CR, TH1* Hist_Bkg, TH1* Hist_BkgTemp, TH1* Hist_Deconv, double rms, double mean, int n_iter_begin, bool includeSR) {
    double n_iter_final = 1;
    double chi2_best = 0.;
    chi2_best = 0.;
    TF1 *func = new TF1("func", "gaus", -50, 50);
    func->SetParameters(10.,1.5,0.5);
    for (int delta_n_iter = 1;delta_n_iter<=20;delta_n_iter++) {
          int n_iter = n_iter_begin-10+delta_n_iter;
          if (n_iter<=0) n_iter = 1;
          double offset_begin = 0;
          double chi2 = 0;
          func->SetParameters(10.,mean,rms);
          Hist_Deconv->Reset();
          Hist_Deconv->FillRandom("func",Hist_SR->Integral()*100);
          Deconvolution(Hist_CR,Hist_Deconv,Hist_BkgTemp,n_iter);
          offset_begin = Hist_SR->GetMean()-Hist_BkgTemp->GetMean();
          offset_begin = ShiftAndNormalize(Hist_SR,Hist_BkgTemp,Hist_Bkg,offset_begin,true,includeSR);
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
        double rms = rms_begin-0.1*rms_begin+double(n_rms)*0.2*rms_begin/50.;
        func->SetParameters(10.,mean,rms);
        Hist_Deconv->Reset();
        Hist_Deconv->FillRandom("func",Hist_SR->Integral()*100);
        Deconvolution(Hist_CR,Hist_Deconv,Hist_BkgTemp,n_iter);
        offset_begin = Hist_SR->GetMean()-Hist_BkgTemp->GetMean();
        offset_begin = ShiftAndNormalize(Hist_SR,Hist_BkgTemp,Hist_Bkg,offset_begin,true,includeSR);
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
vector< std::pair <double,int> > SortList(vector< std::pair <double,int> > list_match)
{
        vector< std::pair <double,int> > sorted_list_match;
        double best_match = 100000.;
        int best_run = 0;
        for (int run1=0;run1<list_match.size();run1++)
        {
            best_match = 100000.;
            best_run = 0;
            for (int run2=0;run2<list_match.size();run2++)
            {
                if (best_match>list_match.at(run2).first)
                {
                    bool is_it_already_in_list = false;
                    for (int run3=0;run3<sorted_list_match.size();run3++)
                    {
                        if (sorted_list_match.at(run3).second==list_match.at(run2).second) is_it_already_in_list = true;
                    }
                    if (is_it_already_in_list) continue;
                    best_match = list_match.at(run2).first;
                    best_run = list_match.at(run2).second;
                }
            }
            if (best_run==0) continue;
            sorted_list_match.push_back(std::make_pair(best_match,best_run));
        }
        return sorted_list_match;
}
vector<int> FindRunSublist(string source, vector<int> Target_runlist, vector<int> *already_used_runs, int number_runs_included)
{
        std::cout << "Getting sublist runs from " << source << std::endl;
        vector<int> sublist;
        vector< std::pair <double,int> > list_match;
        double elev_primary = 0;
        double azim_primary = 0;
        double elev_this = 0;
        double azim_this = 0;
        double ra_primary = 0;
        double dec_primary = 0;
        double ra_this = 0;
        double dec_this = 0;
        char observation[50];
        sprintf(observation, "%s", source.c_str());
        if (TString(source)=="2ndCrab") sprintf(observation, "%s", "Crab");
        for (int run=0;run<Target_runlist.size();run++) {
                char run_number[50];
                bool skip_this_run = false;
                for (int old_run=0;old_run<already_used_runs->size();old_run++)
                {
                    if (Target_runlist[run]==already_used_runs->at(old_run)) 
                    {
                        skip_this_run = true;
                    }
                }
                if (skip_this_run) continue;
                if (sublist.size()==0) 
                {
                    sprintf(run_number, "%i", int(Target_runlist[run]));
                    filename = TString("$VERITAS_USER_DATA_DIR/"+TString(observation)+"_V6_Moderate-TMVA-BDT.RB."+TString(run_number)+".root");
                    if (TString(observation)=="Coma" || TString(observation)=="BrandonValidation") {
                      filename = TString("/veritas/userspace/brandon/VERITAS/Background/anasum/"+TString(run_number)+".anasum.root");
                    }
                    TFile*  input_file = TFile::Open(filename.c_str());
                    TTree* pointing_tree = (TTree*) input_file->Get("run_"+TString(run_number)+"/stereo/pointingDataReduced");
                    pointing_tree->SetBranchAddress("TelElevation",&TelElevation);
                    pointing_tree->SetBranchAddress("TelAzimuth",&TelAzimuth);
                    pointing_tree->SetBranchAddress("TelRAJ2000",&TelRAJ2000);
                    pointing_tree->SetBranchAddress("TelDecJ2000",&TelDecJ2000);
                    pointing_tree->GetEntry(0);
                    //if (TelRAJ2000>3.0 && TelRAJ2000<4.0 && TelDecJ2000>0.2 && TelDecJ2000<0.6) 
                    //{
                    //    input_file->Close();
                    //    continue;
                    //}
                    if (TelElevation<Target_Elev_cut_lower) 
                    {
                        input_file->Close();
                        continue;
                    }
                    if (TelElevation>Target_Elev_cut_upper)
                    {
                        input_file->Close();
                        continue;
                    }
                    if (TelAzimuth<Target_Azim_cut_lower)
                    {
                        input_file->Close();
                        continue;
                    }
                    if (TelAzimuth>Target_Azim_cut_upper)
                    {
                        input_file->Close();
                        continue;
                    }
                    elev_primary = TelElevation;
                    azim_primary = TelAzimuth;
                    ra_primary = TelRAJ2000;
                    dec_primary = TelDecJ2000;
                    sublist.push_back(Target_runlist[run]);
                    already_used_runs->push_back(Target_runlist[run]);
                    input_file->Close();
                    continue;
                }
                sprintf(run_number, "%i", int(Target_runlist[run]));
                filename = TString("$VERITAS_USER_DATA_DIR/"+TString(observation)+"_V6_Moderate-TMVA-BDT.RB."+TString(run_number)+".root");
                if (TString(observation)=="Coma" || TString(observation)=="BrandonValidation") {
                  filename = TString("/veritas/userspace/brandon/VERITAS/Background/anasum/"+TString(run_number)+".anasum.root");
                }
                TFile*  input_file = TFile::Open(filename.c_str());
                TTree* pointing_tree = (TTree*) input_file->Get("run_"+TString(run_number)+"/stereo/pointingDataReduced");
                pointing_tree->SetBranchAddress("TelElevation",&TelElevation);
                pointing_tree->SetBranchAddress("TelAzimuth",&TelAzimuth);
                pointing_tree->SetBranchAddress("TelRAJ2000",&TelRAJ2000);
                pointing_tree->SetBranchAddress("TelDecJ2000",&TelDecJ2000);
                pointing_tree->GetEntry(0);
                //if (TelRAJ2000>3.0 && TelRAJ2000<4.0 && TelDecJ2000>0.2 && TelDecJ2000<0.6) 
                //{
                //    input_file->Close();
                //    continue;
                //}
                if (TelElevation<Target_Elev_cut_lower) 
                {
                    input_file->Close();
                    continue;
                }
                if (TelElevation>Target_Elev_cut_upper)
                {
                    input_file->Close();
                    continue;
                }
                if (TelAzimuth<Target_Azim_cut_lower)
                {
                    input_file->Close();
                    continue;
                }
                if (TelAzimuth>Target_Azim_cut_upper)
                {
                    input_file->Close();
                    continue;
                }
                elev_this = TelElevation;
                azim_this = TelAzimuth;
                ra_this = TelRAJ2000;
                dec_this = TelDecJ2000;
                //double match = pow(elev_primary-elev_this,2)+pow(azim_primary-azim_this,2);
                double match = pow(elev_primary-elev_this,2)+pow(azim_primary-azim_this,2)/9.;
                list_match.push_back(std::make_pair(match,Target_runlist[run]));
                input_file->Close();
        }
        list_match = SortList(list_match);
        for (int run=0;run<list_match.size();run++)
        {
                if (run>=number_runs_included-1) continue;
                sublist.push_back(list_match.at(run).second);
                already_used_runs->push_back(list_match.at(run).second);
                std::cout << "list_match.at(run).first = " << list_match.at(run).first << std::endl;
                std::cout << "list_match.at(run).second = " << list_match.at(run).second << std::endl;
        }
        return sublist;
}
void DeconvolutionMethodForExtendedSources(string target_data, double elev_lower, double elev_upper, double azim_lower, double azim_upper, double theta2_cut_lower_input, double theta2_cut_upper_input) {

        //TH1::SetDefaultSumw2();
        sprintf(target, "%s", target_data.c_str());
        Target_Elev_cut_lower = elev_lower;
        Target_Elev_cut_upper = elev_upper;
        Elev_cut_lower = elev_lower;
        Elev_cut_upper = elev_upper;
        Target_Azim_cut_lower = azim_lower;
        Target_Azim_cut_upper = azim_upper;
        Azim_cut_lower = azim_lower;
        Azim_cut_upper = azim_upper;
        Theta2_cut_lower = theta2_cut_lower_input;
        Theta2_cut_upper = theta2_cut_upper_input;

        //if (TString(Region)=="VR") {
        //        MSCW_cut_lower = 1.5;
        //        MSCW_cut_upper = 2.5;
        //        MSCL_signal_cut_lower = 1.5;
        //        MSCL_signal_cut_upper = 2.5;
        //        MSCL_control_cut_lower = 2.5;
        //        MSCL_control_cut_upper = 3.5;
        //}

        TRandom rnd;
        TH2D Hist_Dark_TelElevAzim("Hist_Dark_TelElevAzim","",18,0,90,18,0,360);
        TH2D Hist_Dark_TelRaDec("Hist_Dark_TelRaDec","",100,0,5,100,-1,1);
        TH2D Hist_Target_TelElevAzim("Hist_Target_TelElevAzim","",18,0,90,18,0,360);
        TH2D Hist_Target_TelRaDec("Hist_Target_TelRaDec","",100,0,5,100,-1,1);
        TH2D Hist_Target_TelRaDec_AfterCut("Hist_Target_TelRaDec_AfterCut","",100,0,5,100,-1,1);
        TH1D Hist_Target_ON_MSCW_Alpha("Hist_Target_ON_MSCW_Alpha","",100,0,10);
        TH1D Hist_Target_OFF_MSCW_Alpha("Hist_Target_OFF_MSCW_Alpha","",100,0,10);
        vector<int> Run_sublist;
        vector<TH1D> Hist_Dark_SR_ErecS;
        vector<TH1D> Hist_Dark_SR_MSCW;
        vector<TH1D> Hist_Dark_SR_MSCL;
        vector<TH1D> Hist_Dark_CR_MSCW;
        vector<TH1D> Hist_Dark_BkgTemp_MSCW;
        vector<TH1D> Hist_Dark_Bkg_MSCW;
        vector<TH1D> Hist_Dark_Deconv_MSCW;
        vector<TH1D> Hist_Dark_TrueDeconv_MSCW;
        vector<TH1D> Hist_Dark_Elec_MSCW;
        vector<TH1D> Hist_Target_SR_ErecS;
        vector<TH1D> Hist_Target_SR_theta2;
        vector<TH1D> Hist_Target_Bkg_theta2;
        vector<TH1D> Hist_Target_Bkg_theta2_Sum;
        vector<TH1D> Hist_Target_CR_theta2;
        vector<TH2D> Hist_Target_SR_RaDec;
        vector<TH2D> Hist_Target_Bkg_RaDec;
        vector<TH2D> Hist_Target_Bkg_RaDec_Sum;
        vector<TH2D> Hist_Target_CR_RaDec;
        vector<TH1D> Hist_Target_SR_MSCW;
        vector<TH1D> Hist_Target_SR_MSCW_Sum;
        vector<TH1D> Hist_Target_SR_MSCL;
        vector<TH1D> Hist_Target_CR1_MSCW;
        vector<TH1D> Hist_Target_CR1_MSCW_Sum;
        vector<TH1D> Hist_Target_CR2_MSCW;
        vector<TH1D> Hist_Target_CR2_MSCW_Sum;
        vector<TH1D> Hist_Target_BkgTemp_MSCW;
        vector<TH1D> Hist_Target_Bkg_MSCW;
        vector<TH1D> Hist_Target_BkgCR1_MSCW;
        vector<TH1D> Hist_Target_Bkg_MSCW_Sum;
        vector<TH1D> Hist_Target_BkgCR1_MSCW_Sum;
        vector<TH1D> Hist_Target_Deconv_MSCW;
        vector<TH1D> Hist_Target_TrueDeconv_MSCW;
        vector<TH1D> Hist_Target_Ring_MSCW;
        vector<TH1D> Hist_Target_Ring_MSCL;
        for (int e=0;e<N_energy_bins;e++) {
            char e_low[50];
            sprintf(e_low, "%i", int(energy_bins[e]));
            char e_up[50];
            sprintf(e_up, "%i", int(energy_bins[e+1]));
            if (energy_bins[e]>=100.) N_bins_for_deconv = 960;
            if (energy_bins[e]>=1500.) N_bins_for_deconv = 480;
            Hist_Dark_SR_ErecS.push_back(TH1D("Hist_Dark_SR_ErecS_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_energy_bins,energy_bins));
            Hist_Dark_SR_MSCW.push_back(TH1D("Hist_Dark_SR_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-30,30));
            Hist_Dark_SR_MSCL.push_back(TH1D("Hist_Dark_SR_MSCL_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-30,30));
            Hist_Dark_CR_MSCW.push_back(TH1D("Hist_Dark_CR_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-30,30));
            Hist_Dark_BkgTemp_MSCW.push_back(TH1D("Hist_Dark_BkgTemp_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-30,30));
            Hist_Dark_Bkg_MSCW.push_back(TH1D("Hist_Dark_Bkg_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-30,30));
            Hist_Dark_Deconv_MSCW.push_back(TH1D("Hist_Dark_Deconv_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-30,30));
            Hist_Dark_TrueDeconv_MSCW.push_back(TH1D("Hist_Dark_TrueDeconv_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-30,30));
            Hist_Target_SR_ErecS.push_back(TH1D("Hist_Target_SR_ErecS_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_energy_bins,energy_bins));
            Hist_Target_SR_theta2.push_back(TH1D("Hist_Target_SR_theta2_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",70,0,7));
            Hist_Target_Bkg_theta2.push_back(TH1D("Hist_Target_Bkg_theta2_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",70,0,7));
            Hist_Target_Bkg_theta2_Sum.push_back(TH1D("Hist_Target_Bkg_theta2_Sum_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",70,0,7));
            Hist_Target_CR_theta2.push_back(TH1D("Hist_Target_CR_theta2_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",70,0,7));
            Hist_Target_SR_RaDec.push_back(TH2D("Hist_Target_SR_RaDec_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",400,-2,2,400,-2,2));
            Hist_Target_Bkg_RaDec.push_back(TH2D("Hist_Target_Bkg_RaDec_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",400,-2,2,400,-2,2));
            Hist_Target_Bkg_RaDec_Sum.push_back(TH2D("Hist_Target_Bkg_RaDec_Sum_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",400,-2,2,400,-2,2));
            Hist_Target_CR_RaDec.push_back(TH2D("Hist_Target_CR_RaDec_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",400,-2,2,400,-2,2));
            Hist_Target_SR_MSCW.push_back(TH1D("Hist_Target_SR_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-30,30));
            Hist_Target_SR_MSCW_Sum.push_back(TH1D("Hist_Target_SR_MSCW_Sum_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-30,30));
            Hist_Target_SR_MSCL.push_back(TH1D("Hist_Target_SR_MSCL_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-30,30));
            Hist_Target_CR1_MSCW.push_back(TH1D("Hist_Target_CR1_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-30,30));
            Hist_Target_CR1_MSCW_Sum.push_back(TH1D("Hist_Target_CR1_MSCW_Sum_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-30,30));
            Hist_Target_CR2_MSCW.push_back(TH1D("Hist_Target_CR2_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-30,30));
            Hist_Target_CR2_MSCW_Sum.push_back(TH1D("Hist_Target_CR2_MSCW_Sum_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-30,30));
            Hist_Target_BkgTemp_MSCW.push_back(TH1D("Hist_Target_BkgTemp_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-30,30));
            Hist_Target_Bkg_MSCW.push_back(TH1D("Hist_Target_Bkg_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-30,30));
            Hist_Target_BkgCR1_MSCW.push_back(TH1D("Hist_Target_BkgCR1_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-30,30));
            Hist_Target_Bkg_MSCW_Sum.push_back(TH1D("Hist_Target_Bkg_MSCW_Sum_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-30,30));
            Hist_Target_BkgCR1_MSCW_Sum.push_back(TH1D("Hist_Target_BkgCR1_MSCW_Sum_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-30,30));
            Hist_Target_Deconv_MSCW.push_back(TH1D("Hist_Target_Deconv_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-30,30));
            Hist_Target_TrueDeconv_MSCW.push_back(TH1D("Hist_Target_TrueDeconv_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-30,30));
            Hist_Target_Ring_MSCW.push_back(TH1D("Hist_Target_Ring_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-30,30));
            Hist_Target_Ring_MSCL.push_back(TH1D("Hist_Target_Ring_MSCL_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-30,30));
        }

        std::cout << "Getting dark runs... " << std::endl;
        //vector<int> Dark_runlist = GetRunList("BrandonValidation");
        //vector<int> Dark_runlist = GetRunList("Crab");
        vector<int> Dark_runlist = GetRunList("PKS1424");
        for (int run=0;run<Dark_runlist.size();run++) {
                char run_number[50];
                sprintf(run_number, "%i", int(Dark_runlist[run]));
                //filename = TString("$VERITAS_USER_DATA_DIR/Crab_V6_Moderate-TMVA-BDT.RB."+TString(run_number)+".root");
                filename = TString("$VERITAS_USER_DATA_DIR/PKS1424_V6_Moderate-TMVA-BDT.RB."+TString(run_number)+".root");
                //filename = TString("/veritas/userspace/brandon/VERITAS/Background/anasum/"+TString(run_number)+".anasum.root");
                TFile*  input_file = TFile::Open(filename.c_str());
                TTree* pointing_tree = (TTree*) input_file->Get("run_"+TString(run_number)+"/stereo/pointingDataReduced");
                pointing_tree->SetBranchAddress("TelElevation",&TelElevation);
                pointing_tree->SetBranchAddress("TelAzimuth",&TelAzimuth);
                pointing_tree->SetBranchAddress("TelRAJ2000",&TelRAJ2000);
                pointing_tree->SetBranchAddress("TelDecJ2000",&TelDecJ2000);
                pointing_tree->GetEntry(0);
                for (int entry=0;entry<pointing_tree->GetEntries();entry++) {
                        pointing_tree->GetEntry(entry);
                        Hist_Dark_TelElevAzim.Fill(TelElevation,TelAzimuth);
                        Hist_Dark_TelRaDec.Fill(TelRAJ2000,TelDecJ2000);
                }
                //if (TelElevation<Elev_cut_lower) continue;
                //if (TelElevation>Elev_cut_upper) continue;
                //if (TelAzimuth<Azim_cut_lower) continue;
                //if (TelAzimuth>Azim_cut_upper) continue;

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
                            if (ControlSelection1MSCW()) Hist_Dark_CR_MSCW.at(e).Fill(MSCW);
                        }
                }
                input_file->Close();

        }

        std::cout << "Getting target runs... " << std::endl;
        char observation[50];
        sprintf(observation, "%s", target);
        if (TString(target)=="2ndCrab") sprintf(observation, "%s", "Crab");
        vector<int> already_used_runs;
        vector<int> Target_runlist = GetRunList(target);

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

        for (int e=0;e<N_energy_bins;e++)
        {
            std::cout << "=================================================================" << std::endl;
            std::cout << "Target, e " << energy_bins[e] << std::endl;
            vector<int> Sublist;
            bool DoNextSublist = true;
            already_used_runs.clear();
            while (DoNextSublist) 
            {
                Hist_Target_SR_MSCW.at(e).Reset();
                Hist_Target_Bkg_MSCW.at(e).Reset();
                Hist_Target_BkgCR1_MSCW.at(e).Reset();
                Hist_Target_CR1_MSCW.at(e).Reset();
                Hist_Target_CR2_MSCW.at(e).Reset();
                Hist_Target_CR_theta2.at(e).Reset();
                Hist_Target_Bkg_theta2.at(e).Reset();
                Hist_Target_CR_RaDec.at(e).Reset();
                Hist_Target_Bkg_RaDec.at(e).Reset();
                Sublist = FindRunSublist(target,Target_runlist,&already_used_runs,number_runs_included[e]);
                std::cout << "Sublist.size() = " << Sublist.size() << std::endl;
                if (Sublist.size()==0) 
                {
                    DoNextSublist = false;
                    continue;
                }
                for (int run=0;run<Sublist.size();run++)
                {
                    char run_number[50];
                    sprintf(run_number, "%i", int(Sublist[run]));
                    std::cout << "Reading run " << run_number << std::endl;
                    filename = TString("$VERITAS_USER_DATA_DIR/"+TString(observation)+"_V6_Moderate-TMVA-BDT.RB."+TString(run_number)+".root");
                    if (TString(observation)=="Coma" || TString(observation)=="BrandonValidation") {
                      filename = TString("/veritas/userspace/brandon/VERITAS/Background/anasum/"+TString(run_number)+".anasum.root");
                    }
                    TFile*  input_file = TFile::Open(filename.c_str());
                    TTree* pointing_tree = (TTree*) input_file->Get("run_"+TString(run_number)+"/stereo/pointingDataReduced");
                    pointing_tree->SetBranchAddress("TelElevation",&TelElevation);
                    pointing_tree->SetBranchAddress("TelAzimuth",&TelAzimuth);
                    pointing_tree->SetBranchAddress("TelRAJ2000",&TelRAJ2000);
                    pointing_tree->SetBranchAddress("TelDecJ2000",&TelDecJ2000);
                    for (int entry=0;entry<pointing_tree->GetEntries();entry++) {
                            pointing_tree->GetEntry(entry);
                            Hist_Target_TelElevAzim.Fill(TelElevation,TelAzimuth);
                            Hist_Target_TelRaDec.Fill(TelRAJ2000,TelDecJ2000);
                    }
                    pointing_tree->GetEntry(0);
                    //if (TelRAJ2000>3.0 && TelRAJ2000<4.0 && TelDecJ2000>0.2 && TelDecJ2000<0.6) 
                    //{
                    //    input_file->Close();
                    //    continue;
                    //}
                    if (TelElevation<Target_Elev_cut_lower) 
                    {
                        input_file->Close();
                        continue;
                    }
                    if (TelElevation>Target_Elev_cut_upper)
                    {
                        input_file->Close();
                        continue;
                    }
                    if (TelAzimuth<Target_Azim_cut_lower)
                    {
                        input_file->Close();
                        continue;
                    }
                    if (TelAzimuth>Target_Azim_cut_upper)
                    {
                        input_file->Close();
                        continue;
                    }
                    for (int entry=0;entry<pointing_tree->GetEntries();entry++) {
                            pointing_tree->GetEntry(entry);
                            Hist_Target_TelRaDec_AfterCut.Fill(TelRAJ2000,TelDecJ2000);
                    }
                    TString root_file = "run_"+TString(run_number)+"/stereo/data_on";
                    TTree* Target_tree = (TTree*) input_file->Get(root_file);
                    Target_tree->SetBranchAddress("theta2",&theta2);
                    Target_tree->SetBranchAddress("ra",&ra_sky);
                    Target_tree->SetBranchAddress("dec",&dec_sky);
                    Target_tree->SetBranchAddress("ErecS",&ErecS);
                    Target_tree->SetBranchAddress("EChi2S",&EChi2S);
                    Target_tree->SetBranchAddress("MSCW",&MSCW);
                    Target_tree->SetBranchAddress("MSCL",&MSCL);
                    //Target_tree->SetBranchAddress("MWR",&MSCW);
                    //Target_tree->SetBranchAddress("MLR",&MSCL);
                    for (int entry=0;entry<Target_tree->GetEntries();entry++) {
                            Target_tree->GetEntry(entry);
                            int energy = Hist_Target_SR_ErecS.at(0).FindBin(ErecS*1000.)-1;
                            if (e!=energy) continue;
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
                                    if (SignalSelectionMSCW()&&SignalSelectionMSCL()) 
                                    {
                                        Hist_Target_SR_theta2.at(e).Fill(theta2);
                                        Hist_Target_SR_RaDec.at(e).Fill(ra_sky,dec_sky);
                                    }
                                    if (MSCL>MSCL_signal_cut_upper*1.0 || MSCW>MSCW_cut_upper*1.0) 
                                    {
                                        if (MSCL<MSCL_signal_cut_upper*2.0 && MSCW<MSCW_cut_upper*2.0) 
                                        {
                                            Hist_Target_CR_theta2.at(e).Fill(theta2);
                                            Hist_Target_CR_RaDec.at(e).Fill(ra_sky,dec_sky);
                                        }
                                    }
                                    if (SignalSelectionMSCW()) Hist_Target_SR_MSCW.at(e).Fill(MSCW);
                                    if (SignalSelectionMSCL()) Hist_Target_SR_MSCL.at(e).Fill(MSCL);
                                    if (ControlSelection1MSCW()) Hist_Target_CR1_MSCW.at(e).Fill(MSCW);
                                    if (ControlSelection2MSCW()) Hist_Target_CR2_MSCW.at(e).Fill(MSCW);
                            }
                            else if (RingFoV()) {
                                    if (SignalSelectionMSCW()) Hist_Target_Ring_MSCW.at(e).Fill(MSCW);
                                    if (SignalSelectionMSCL()) Hist_Target_Ring_MSCL.at(e).Fill(MSCL);
                            }
                    }
                    input_file->Close();
                }
                std::cout << "=================================================================" << std::endl;
                std::cout << "Target, e " << energy_bins[e] << ", N_iter.at(e) = " << N_iter.at(e) << std::endl;
                std::cout << "Hist_Target_SR_MSCW.at(e).Integral() = " << Hist_Target_SR_MSCW.at(e).Integral() << std::endl;
                std::cout << "Hist_Target_SR_MSCW.at(e).GetMean() = " << Hist_Target_SR_MSCW.at(e).GetMean() << std::endl;
                std::cout << "Hist_Target_SR_MSCW.at(e).GetRMS() = " << Hist_Target_SR_MSCW.at(e).GetRMS() << std::endl;

                Deconvolution(&Hist_Target_CR2_MSCW.at(e),&Hist_Target_CR1_MSCW.at(e),&Hist_Target_Deconv_MSCW.at(e),N_iter.at(e));
                double rms_begin = Hist_Target_Deconv_MSCW.at(e).GetRMS();
                double mean_begin = Hist_Target_Deconv_MSCW.at(e).GetMean();

                N_iter.at(e) = FindNIteration(&Hist_Target_CR1_MSCW.at(e),&Hist_Target_CR2_MSCW.at(e),&Hist_Target_BkgCR1_MSCW.at(e),&Hist_Target_BkgTemp_MSCW.at(e),&Hist_Target_Deconv_MSCW.at(e),rms_begin,mean_begin,20,true);
                std::cout << "Target, e " << energy_bins[e] << ", N_iter.at(e) = " << N_iter.at(e) << std::endl;
                N_rms.at(e) = FindRMS(&Hist_Target_CR1_MSCW.at(e),&Hist_Target_CR2_MSCW.at(e),&Hist_Target_BkgCR1_MSCW.at(e),&Hist_Target_BkgTemp_MSCW.at(e),&Hist_Target_Deconv_MSCW.at(e),rms_begin,mean_begin,N_iter.at(e),true);
                std::cout << "Target, e " << energy_bins[e] << ", N_rms.at(e) = " << N_rms.at(e) << std::endl;

                N_iter.at(e) = FindNIteration(&Hist_Target_CR1_MSCW.at(e),&Hist_Target_CR2_MSCW.at(e),&Hist_Target_BkgCR1_MSCW.at(e),&Hist_Target_BkgTemp_MSCW.at(e),&Hist_Target_Deconv_MSCW.at(e),N_rms.at(e),mean_begin,N_iter.at(e),true);
                std::cout << "Target, e " << energy_bins[e] << ", N_iter.at(e) = " << N_iter.at(e) << std::endl;
                N_rms.at(e) = FindRMS(&Hist_Target_CR1_MSCW.at(e),&Hist_Target_CR2_MSCW.at(e),&Hist_Target_BkgCR1_MSCW.at(e),&Hist_Target_BkgTemp_MSCW.at(e),&Hist_Target_Deconv_MSCW.at(e),N_rms.at(e),mean_begin,N_iter.at(e),true);
                std::cout << "Target, e " << energy_bins[e] << ", N_rms.at(e) = " << N_rms.at(e) << std::endl;

                N_iter.at(e) = FindNIteration(&Hist_Target_CR1_MSCW.at(e),&Hist_Target_CR2_MSCW.at(e),&Hist_Target_BkgCR1_MSCW.at(e),&Hist_Target_BkgTemp_MSCW.at(e),&Hist_Target_Deconv_MSCW.at(e),N_rms.at(e),mean_begin,N_iter.at(e),true);
                std::cout << "Target, e " << energy_bins[e] << ", N_iter.at(e) = " << N_iter.at(e) << std::endl;
                N_rms.at(e) = FindRMS(&Hist_Target_CR1_MSCW.at(e),&Hist_Target_CR2_MSCW.at(e),&Hist_Target_BkgCR1_MSCW.at(e),&Hist_Target_BkgTemp_MSCW.at(e),&Hist_Target_Deconv_MSCW.at(e),N_rms.at(e),mean_begin,N_iter.at(e),true);
                std::cout << "Target, e " << energy_bins[e] << ", N_rms.at(e) = " << N_rms.at(e) << std::endl;


                myfunc->SetParameters(10.,mean_begin,N_rms.at(e));
                Hist_Target_Deconv_MSCW.at(e).Reset();
                Hist_Target_Deconv_MSCW.at(e).FillRandom("myfunc",Hist_Target_SR_MSCW.at(e).Integral()*100);
                Deconvolution(&Hist_Target_CR2_MSCW.at(e),&Hist_Target_Deconv_MSCW.at(e),&Hist_Target_BkgTemp_MSCW.at(e),N_iter.at(e));
                double offset_begin = 0;
                offset_begin = Hist_Target_CR1_MSCW.at(e).GetMean()-Hist_Target_BkgTemp_MSCW.at(e).GetMean();
                offset_begin = ShiftAndNormalize(&Hist_Target_CR1_MSCW.at(e),&Hist_Target_BkgTemp_MSCW.at(e),&Hist_Target_BkgCR1_MSCW.at(e),offset_begin,true,true);

                //Deconvolution(&Hist_Target_CR1_MSCW.at(e),&Hist_Target_SR_MSCW.at(e),&Hist_Target_Deconv_MSCW.at(e),N_iter.at(e));
                //double rms_begin = Hist_Target_Deconv_MSCW.at(e).GetRMS();
                //double mean_begin = Hist_Target_Deconv_MSCW.at(e).GetMean();

                //N_iter.at(e) = FindNIteration(&Hist_Target_SR_MSCW.at(e),&Hist_Target_CR1_MSCW.at(e),&Hist_Target_Bkg_MSCW.at(e),&Hist_Target_BkgTemp_MSCW.at(e),&Hist_Target_Deconv_MSCW.at(e),rms_begin,mean_begin,20,false);
                //std::cout << "Target, e " << energy_bins[e] << ", N_iter.at(e) = " << N_iter.at(e) << std::endl;
                //N_rms.at(e) = FindRMS(&Hist_Target_SR_MSCW.at(e),&Hist_Target_CR1_MSCW.at(e),&Hist_Target_Bkg_MSCW.at(e),&Hist_Target_BkgTemp_MSCW.at(e),&Hist_Target_Deconv_MSCW.at(e),rms_begin,mean_begin,N_iter.at(e),false);
                //std::cout << "Target, e " << energy_bins[e] << ", N_rms.at(e) = " << N_rms.at(e) << std::endl;

                //N_iter.at(e) = FindNIteration(&Hist_Target_SR_MSCW.at(e),&Hist_Target_CR1_MSCW.at(e),&Hist_Target_Bkg_MSCW.at(e),&Hist_Target_BkgTemp_MSCW.at(e),&Hist_Target_Deconv_MSCW.at(e),N_rms.at(e),mean_begin,N_iter.at(e),false);
                //std::cout << "Target, e " << energy_bins[e] << ", N_iter.at(e) = " << N_iter.at(e) << std::endl;
                N_rms.at(e) = FindRMS(&Hist_Target_SR_MSCW.at(e),&Hist_Target_CR1_MSCW.at(e),&Hist_Target_Bkg_MSCW.at(e),&Hist_Target_BkgTemp_MSCW.at(e),&Hist_Target_Deconv_MSCW.at(e),N_rms.at(e),mean_begin,N_iter.at(e),false);
                std::cout << "Target, e " << energy_bins[e] << ", N_rms.at(e) = " << N_rms.at(e) << std::endl;

                N_mean.at(e) = mean_begin;

                double rms = N_rms.at(e);
                double mean = N_mean.at(e);
                myfunc->SetParameters(10.,mean,rms);
                Hist_Target_Deconv_MSCW.at(e).Reset();
                Hist_Target_Deconv_MSCW.at(e).FillRandom("myfunc",Hist_Target_SR_MSCW.at(e).Integral()*100);

                Deconvolution(&Hist_Target_CR1_MSCW.at(e),&Hist_Target_Deconv_MSCW.at(e),&Hist_Target_BkgTemp_MSCW.at(e),N_iter.at(e));
                Deconvolution(&Hist_Target_CR1_MSCW.at(e),&Hist_Target_SR_MSCW.at(e),&Hist_Target_TrueDeconv_MSCW.at(e),N_iter.at(e));
                std::cout << "Target, e " << energy_bins[e] << ", final rms = " << Hist_Target_Deconv_MSCW.at(e).GetRMS() << std::endl;
                std::cout << "Target, e " << energy_bins[e] << ", true RMS = " << Hist_Target_TrueDeconv_MSCW.at(e).GetRMS() << std::endl;
                offset_begin = Hist_Target_SR_MSCW.at(e).GetMean()-Hist_Target_BkgTemp_MSCW.at(e).GetMean();
                std::cout << "offset_begin = " << offset_begin << std::endl;
                std::cout << "Target, e " << energy_bins[e] << ", offset_begin = " << offset_begin << std::endl;
                offset_begin = ShiftAndNormalize(&Hist_Target_SR_MSCW.at(e),&Hist_Target_BkgTemp_MSCW.at(e),&Hist_Target_Bkg_MSCW.at(e),offset_begin,true,false);
                N_shift.at(e) = offset_begin;

                // predict bkg theta2 distribution
                int norm_bin_low_target = Hist_Target_Bkg_MSCW.at(e).FindBin(MSCW_cut_lower);
                int norm_bin_up_target = Hist_Target_Bkg_MSCW.at(e).FindBin(MSCW_cut_upper);
                double total_bkg = (double) Hist_Target_Bkg_MSCW.at(e).Integral(norm_bin_low_target,norm_bin_up_target);
                double old_integral = (double) Hist_Target_CR_theta2.at(e).Integral();
                double total_bkg_err = 0.;
                for (int bin=norm_bin_low_target;bin<=norm_bin_up_target;bin++)
                {
                    total_bkg_err += Hist_Target_Bkg_MSCW.at(e).GetBinError(bin); // remember this error is correlated.
                }
                double scale_sky = total_bkg/old_integral;
                double scale_sky_err = total_bkg_err/old_integral;
                for (int bin=1;bin<=Hist_Target_Bkg_theta2.at(e).GetNbinsX();bin++)
                {
                    double bin_content = Hist_Target_CR_theta2.at(e).GetBinContent(bin)*scale_sky;
                    double bin_error_1 = Hist_Target_CR_theta2.at(e).GetBinError(bin)*scale_sky;
                    double bin_error_2 = Hist_Target_CR_theta2.at(e).GetBinContent(bin)*scale_sky_err;
                    Hist_Target_Bkg_theta2.at(e).SetBinContent(bin,bin_content);
                    Hist_Target_Bkg_theta2.at(e).SetBinError(bin,pow(bin_error_1*bin_error_1+bin_error_2*bin_error_2,0.5));
                }
                for (int binx=1;binx<=Hist_Target_Bkg_RaDec.at(e).GetNbinsX();binx++)
                {
                    for (int biny=1;biny<=Hist_Target_Bkg_RaDec.at(e).GetNbinsY();biny++)
                    {
                        double bin_content = Hist_Target_CR_RaDec.at(e).GetBinContent(binx,biny)*scale_sky;
                        double bin_error_1 = Hist_Target_CR_RaDec.at(e).GetBinError(binx,biny)*scale_sky;
                        double bin_error_2 = Hist_Target_CR_RaDec.at(e).GetBinContent(binx,biny)*scale_sky_err;
                        Hist_Target_Bkg_RaDec.at(e).SetBinContent(binx,biny,bin_content);
                        Hist_Target_Bkg_RaDec.at(e).SetBinError(binx,biny,pow(bin_error_1*bin_error_1+bin_error_2*bin_error_2,0.5));
                    }
                }

                Hist_Target_SR_MSCW_Sum.at(e).Add(&Hist_Target_SR_MSCW.at(e));
                Hist_Target_CR1_MSCW_Sum.at(e).Add(&Hist_Target_CR1_MSCW.at(e));
                Hist_Target_CR2_MSCW_Sum.at(e).Add(&Hist_Target_CR2_MSCW.at(e));
                Hist_Target_Bkg_MSCW_Sum.at(e).Add(&Hist_Target_Bkg_MSCW.at(e));
                Hist_Target_BkgCR1_MSCW_Sum.at(e).Add(&Hist_Target_BkgCR1_MSCW.at(e));
                Hist_Target_Bkg_theta2_Sum.at(e).Add(&Hist_Target_Bkg_theta2.at(e));
                Hist_Target_Bkg_RaDec_Sum.at(e).Add(&Hist_Target_Bkg_RaDec.at(e));

            }
        }

/*
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
                        Hist_Dark_SR_MSCW.at(e).Rebin(2);
                        Hist_Dark_SR_MSCL.at(e).Rebin(2);
                        Hist_Dark_CR_MSCW.at(e).Rebin(2);
                        Hist_Dark_Bkg_MSCW.at(e).Rebin(2);
                        Hist_Dark_BkgTemp_MSCW.at(e).Rebin(2);
                        Hist_Dark_Deconv_MSCW.at(e).Rebin(2);
                        Hist_Dark_TrueDeconv_MSCW.at(e).Rebin(2);
                        Hist_Target_SR_MSCW.at(e).Rebin(2);
                        Hist_Target_SR_MSCW_Sum.at(e).Rebin(2);
                        Hist_Target_SR_MSCL.at(e).Rebin(2);
                        Hist_Target_CR1_MSCW.at(e).Rebin(2);
                        Hist_Target_Bkg_MSCW.at(e).Rebin(2);
                        Hist_Target_Bkg_MSCW_Sum.at(e).Rebin(2);
                        Hist_Target_BkgTemp_MSCW.at(e).Rebin(2);
                        Hist_Target_Deconv_MSCW.at(e).Rebin(2);
                        Hist_Target_TrueDeconv_MSCW.at(e).Rebin(2);
                        Hist_Target_Ring_MSCW.at(e).Rebin(2);
                        Hist_Target_Ring_MSCL.at(e).Rebin(2);
                }
                else break;
            }
        }
*/

        // simple dark template method
        std::cout << "Simple dark template... " << std::endl;
        for (int e=0;e<N_energy_bins;e++) {
                int norm_bin_low_target = Hist_Target_Bkg_MSCW_Sum.at(e).FindBin(Norm_Lower);
                int norm_bin_up_target = Hist_Target_Bkg_MSCW_Sum.at(e).FindBin(Norm_Upper);
                double scale_target = Hist_Target_SR_MSCW_Sum.at(e).Integral(norm_bin_low_target,norm_bin_up_target)/Hist_Dark_SR_MSCW.at(e).Integral(norm_bin_low_target,norm_bin_up_target);
                Hist_Dark_SR_MSCW.at(e).Scale(scale_target);
                Hist_Dark_SR_MSCL.at(e).Scale(scale_target);
                Hist_Dark_Bkg_MSCW.at(e).Scale(scale_target);
        }

        // Get ring bkg
        for (int e=0;e<N_energy_bins;e++) {
                int norm_bin_low_ring = Hist_Target_Ring_MSCW.at(e).FindBin(Norm_Lower);
                int norm_bin_up_ring = Hist_Target_Ring_MSCW.at(e).FindBin(Norm_Upper);
                double scale_ring = Hist_Target_SR_MSCW_Sum.at(e).Integral(norm_bin_low_ring,norm_bin_up_ring)/Hist_Target_Ring_MSCW.at(e).Integral(norm_bin_low_ring,norm_bin_up_ring);
                if (!(scale_ring>0)) scale_ring = 0;
                if (!(Hist_Target_Ring_MSCW.at(e).Integral(norm_bin_low_ring,norm_bin_up_ring)>0)) scale_ring = 0;
                Hist_Target_Ring_MSCW.at(e).Scale(scale_ring);
                Hist_Target_Ring_MSCL.at(e).Scale(scale_ring);
        }


        int norm_bin_low_ring = Hist_Target_ON_MSCW_Alpha.FindBin(Norm_Lower);
        int norm_bin_up_ring = Hist_Target_ON_MSCW_Alpha.FindBin(Norm_Upper);
        double scale_ring = Hist_Target_ON_MSCW_Alpha.Integral(norm_bin_low_ring,norm_bin_up_ring)/Hist_Target_OFF_MSCW_Alpha.Integral(norm_bin_low_ring,norm_bin_up_ring);

        TFile OutputFile("output/Deconvolution_"+TString(target)+"_Elev"+std::to_string(int(Target_Elev_cut_lower))+"to"+std::to_string(int(Target_Elev_cut_upper))+"_Azim"+std::to_string(int(Target_Azim_cut_lower))+"to"+std::to_string(int(Target_Azim_cut_upper))+"_Theta2"+std::to_string(int(Theta2_cut_lower))+"to"+std::to_string(int(Theta2_cut_upper))+"_"+TString(Region)+".root","recreate"); 
        TTree InfoTree("InfoTree","info tree");
        InfoTree.Branch("MSCW_cut_lower",&MSCW_cut_lower,"MSCW_cut_lower/D");
        InfoTree.Branch("MSCW_cut_upper",&MSCW_cut_upper,"MSCW_cut_upper/D");
        InfoTree.Branch("MSCL_signal_cut_lower",&MSCL_signal_cut_lower,"MSCL_signal_cut_lower/D");
        InfoTree.Branch("MSCL_signal_cut_upper",&MSCL_signal_cut_upper,"MSCL_signal_cut_upper/D");
        InfoTree.Branch("MSCL_control_cut_lower",&MSCL_control1_cut_lower,"MSCL_control_cut_lower/D");
        InfoTree.Branch("MSCL_control_cut_upper",&MSCL_control1_cut_upper,"MSCL_control_cut_upper/D");
        InfoTree.Branch("Elev_cut_lower",&Elev_cut_lower,"Elev_cut_lower/D");
        InfoTree.Branch("Elev_cut_upper",&Elev_cut_upper,"Elev_cut_upper/D");
        InfoTree.Branch("Target_Elev_cut_lower",&Target_Elev_cut_lower,"Target_Elev_cut_lower/D");
        InfoTree.Branch("Target_Elev_cut_upper",&Target_Elev_cut_upper,"Target_Elev_cut_upper/D");
        InfoTree.Branch("Azim_cut_lower",&Azim_cut_lower,"Azim_cut_lower/D");
        InfoTree.Branch("Azim_cut_upper",&Azim_cut_upper,"Azim_cut_upper/D");
        InfoTree.Branch("Target_Azim_cut_lower",&Target_Azim_cut_lower,"Target_Azim_cut_lower/D");
        InfoTree.Branch("Target_Azim_cut_upper",&Target_Azim_cut_upper,"Target_Azim_cut_upper/D");
        InfoTree.Fill();
        InfoTree.Write();
        Hist_Target_TelElevAzim.Write();
        Hist_Target_TelRaDec.Write();
        Hist_Target_TelRaDec_AfterCut.Write();
        Hist_Dark_TelElevAzim.Write();
        Hist_Dark_TelRaDec.Write();
        for (int e=0;e<N_energy_bins;e++) {
                Hist_Dark_SR_MSCW.at(e).Write();
                Hist_Dark_SR_MSCL.at(e).Write();
                Hist_Dark_CR_MSCW.at(e).Write();
                Hist_Dark_Bkg_MSCW.at(e).Write();
                Hist_Dark_Deconv_MSCW.at(e).Write();
                Hist_Dark_TrueDeconv_MSCW.at(e).Write();
                Hist_Target_SR_theta2.at(e).Write();
                Hist_Target_Bkg_theta2_Sum.at(e).Write();
                Hist_Target_CR_theta2.at(e).Write();
                Hist_Target_SR_RaDec.at(e).Write();
                Hist_Target_Bkg_RaDec_Sum.at(e).Write();
                Hist_Target_CR_RaDec.at(e).Write();
                Hist_Target_SR_MSCW_Sum.at(e).Write();
                Hist_Target_SR_MSCL.at(e).Write();
                Hist_Target_CR1_MSCW_Sum.at(e).Write();
                Hist_Target_CR1_MSCW.at(e).Write();
                Hist_Target_CR2_MSCW_Sum.at(e).Write();
                Hist_Target_CR2_MSCW.at(e).Write();
                Hist_Target_Bkg_MSCW_Sum.at(e).Write();
                Hist_Target_BkgCR1_MSCW_Sum.at(e).Write();
                Hist_Target_Deconv_MSCW.at(e).Write();
                Hist_Target_TrueDeconv_MSCW.at(e).Write();
                Hist_Target_Ring_MSCW.at(e).Write();
                Hist_Target_Ring_MSCL.at(e).Write();
        }
        OutputFile.Close();

}
