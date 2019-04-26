
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

#include <algorithm>
#include <functional>

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

// VEGAS
//#define VEGAS
//#include "aclicPreProcCommands.h"
//#include <Riostream.h>
//#include <VAShowerData.h>

// VEGAS
#ifndef VEGAS
bool UseVegas =false;
#else
bool UseVegas = true;
#endif

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
double Theta2_upper_limit = 0;

#ifndef VEGAS
// EVDISP
double MSCW_cut_lower = -0.8;  // set this to -1.5 for V5
double MSCW_cut_blind = 1.5;
double MSCW_cut_upper = 1.5;
const int Number_of_SR = 6;
double MSCL_signal_cut_lower[Number_of_SR] = {0.75,0.50,0.25,0.00,-0.25,-0.50};
double MSCL_signal_cut_upper[Number_of_SR] = {1.00,0.75,0.50,0.25, 0.00,-0.25};
const int Number_of_CR = 2;
double MSCL_control_cut_lower[Number_of_CR] = {1.25,1.00};
double MSCL_control_cut_upper[Number_of_CR] = {1.50,1.25};
//double MSCW_cut_lower = -0.8;  // set this to -1.5 for V5
//double MSCW_cut_blind = 1.5;
//double MSCW_cut_upper = 1.5;
//const int Number_of_SR = 1;
//double MSCL_signal_cut_lower[Number_of_SR] = {0.75};
//double MSCL_signal_cut_upper[Number_of_SR] = {1.00};
//const int Number_of_CR = 2;
//double MSCL_control_cut_lower[Number_of_CR] = {1.25,1.00};
//double MSCL_control_cut_upper[Number_of_CR] = {1.50,1.25};
#else
// VEGAS
//double MSCW_cut_lower = 0.7;
//double MSCW_cut_blind = 1.3;
//double MSCW_cut_upper = 1.3;
//const int Number_of_SR = 12;
//double MSCL_signal_cut_lower[Number_of_SR] = {1.25,1.20,1.15,1.10,1.05,1.00,0.95,0.90,0.85,0.80,0.75,0.70};
//double MSCL_signal_cut_upper[Number_of_SR] = {1.30,1.25,1.20,1.15,1.10,1.05,1.00,0.95,0.90,0.85,0.80,0.75};
//const int Number_of_CR = 2;
//double MSCL_control_cut_lower[Number_of_CR] = {1.35,1.30};
//double MSCL_control_cut_upper[Number_of_CR] = {1.40,1.35};
double MSCW_cut_lower = 0.7;
double MSCW_cut_blind = 1.3;
double MSCW_cut_upper = 1.3;
const int Number_of_SR = 6;
double MSCL_signal_cut_lower[Number_of_SR] = {1.20,1.10,1.00,0.90,0.80,0.70};
double MSCL_signal_cut_upper[Number_of_SR] = {1.30,1.20,1.10,1.00,0.90,0.80};
const int Number_of_CR = 2;
double MSCL_control_cut_lower[Number_of_CR] = {1.40,1.30};
double MSCL_control_cut_upper[Number_of_CR] = {1.50,1.40};
#endif


double Norm_Lower = MSCW_cut_blind;
double Norm_Upper = 20.0;

string  filename;
double TelElevation = 0;
double TelAzimuth = 0;
double TelRAJ2000 = 0;
double TelDecJ2000 = 0;
double ErecS = 0;
double EChi2S = 0;
double MSCW = 0;
int NImages = 0;
double MSCL = 0;
double Time = 0;
double SlantDepth = 0;
float EmissionHeight = 0;
float EmissionHeightChi2 = 0;
double Xoff = 0;
double Yoff = 0;
double R2off = 0;
double theta2 = 0;
double ra_sky = 0;
double dec_sky = 0;
vector<int> used_runs;
vector<double> energy_vec;
double exposure_hours = 0.;
vector<double> blinded_cr_rms;
vector<double> blinded_sr_rms;
vector<double> initial_kernel_rms;
vector<double> final_kernel_rms;
vector<double> unblinded_kernel_rms;
vector<double> initial_shift;
vector<double> blinded_shift;
vector<double> unblinded_shift;

//const int N_energy_bins = 16;
//double energy_bins[N_energy_bins+1] =     {200,237,282,335,398,473,562,663,794,937,1122,1585,2239,3162,4467,6310,8913};
//int number_runs_included[N_energy_bins] = {1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,2   ,2   ,4   ,4   ,8   ,8};
const int N_energy_bins = 11;
double energy_bins[N_energy_bins+1] =     {200,282,398,562,794,1122,1585,2239,3162,4467,6310,8913};
int number_runs_included[N_energy_bins] = {99 ,99 ,99 ,99 ,99 ,99  ,99  ,99  ,99  ,99  ,99};
//const int N_energy_bins = 6;
//double energy_bins[N_energy_bins+1] =     {1122,1585,2239,3162,4467,6310,8913};
//int number_runs_included[N_energy_bins] = {4   ,4   ,8   ,8   ,16   ,64};
//const int N_energy_bins = 1;
//double energy_bins[N_energy_bins+1] =     {562,794};
//int number_runs_included[N_energy_bins] = {99};
//const int N_energy_bins = 1;
//double energy_bins[N_energy_bins+1] =     {6310,8913};
//int number_runs_included[N_energy_bins] = {8};
//const int N_energy_bins = 1;
//double energy_bins[N_energy_bins+1] =     {282,398};
//int number_runs_included[N_energy_bins] = {1};

int N_bins_for_deconv = 480;
double MSCW_plot_lower = -30.;
double MSCW_plot_upper = 30.;

TH1D Hist_Target_TelElev("Hist_Target_TelElev","",18,0,90);
TH1D Hist_Dark_TelElev("Hist_Dark_TelElev","",18,0,90);
TH2D Hist_Target_TelElevAzim("Hist_Target_TelElevAzim","",18,0,90,18,0,360);
TH2D Hist_Target_TelRaDec("Hist_Target_TelRaDec","",100,0,5,100,-1,1);

bool FoV() {
    //if (theta2<0.2) return false;
    //if (R2off<Theta2_cut_lower) return false;
    //if (R2off>Theta2_cut_upper) return false;
    if (theta2<Theta2_cut_lower) return false;
    if (theta2>Theta2_cut_upper) return false;
    return true;
}
bool RingFoV() {
    if (theta2<Theta2_cut_upper) return false;
    return true;
}
bool SelectNImages(int Nmin, int Nmax)
{
    if (NImages<Nmin) return false;
    if (NImages>Nmax) return false;
    return true;
}
bool SignalSelectionMSCL() {
    if (MSCW<MSCW_cut_lower) return false;
    if (MSCW>MSCW_cut_upper) return false;
    return true;
}
bool SignalSelectionMSCW(int whichSR) {
    if (MSCL<MSCL_signal_cut_lower[whichSR]) return false;
    if (MSCL>MSCL_signal_cut_upper[whichSR]) return false;
    return true;
}
bool ControlSelectionMSCW(int whichCR) {
    if (MSCL<MSCL_control_cut_lower[whichCR]) return false;
    if (MSCL>MSCL_control_cut_upper[whichCR]) return false;
    return true;
}
bool ControlSelectionTheta2(int whichSR)
{
    double deltaL = MSCL_signal_cut_upper[0]-MSCL_signal_cut_lower[0];
    if (MSCL<MSCL_signal_cut_upper[0]+(whichSR)*deltaL) return false;
    if (MSCL>MSCL_signal_cut_upper[0]+(whichSR+1)*deltaL) return false;
    if (MSCW<MSCW_cut_blind*1.0) return false;
    if (MSCW>MSCW_cut_blind*3.0) return false;
    return true;
}
double background(Double_t *x, Double_t *par) {
    return (1./par[0])*(x[0]+par[1])*exp(-pow((x[0]+par[1])/par[0],2));
}
double Kernel(Double_t *x, Double_t *par) {
    double xx =x[0];
    return exp(-0.5*pow((xx)/(par[0]),2));
}
double GetLinearUncWithRange(TH1* Hist_Bkg, double lower_end, double upper_end) {
    double total_err = 0.;
    for (int i=0;i<Hist_Bkg->GetNbinsX();i++) {
        double bkg = Hist_Bkg->GetBinContent(i+1);
        double bkg_err = Hist_Bkg->GetBinError(i+1);
        if (Hist_Bkg->GetBinCenter(i+1)<lower_end || Hist_Bkg->GetBinCenter(i+1)>upper_end) {
            continue;
        }
        total_err += bkg_err;
    }
    return total_err;
}
double GetBlindedRMS(TH1* Hist_SR)
{
    TH1D Hist_SR_blinded = TH1D("Hist_SR_blinded","",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper);
    for (int i=0;i<Hist_SR->GetNbinsX();i++) {
        if (Hist_SR->GetBinCenter(i+1)<MSCW_cut_blind) {
            continue;
        }
        Hist_SR_blinded.SetBinContent(i+1,Hist_SR->GetBinContent(i+1));
        Hist_SR_blinded.SetBinError(i+1,Hist_SR->GetBinError(i+1));
    }
    return Hist_SR_blinded.GetRMS();
}
double GetChi2WithRange(TH1* Hist_SR, TH1* Hist_Bkg, double lower_end, double upper_end) {
    double chi2_temp = 0.;
    double sign_chi2_temp = 0.;
    double mean = Hist_SR->GetMean();
    double rms = Hist_SR->GetRMS();
    double nbins = 0.;
    for (int i=0;i<Hist_SR->GetNbinsX();i++) {
        double bkg = Hist_Bkg->GetBinContent(i+1);
        double data = Hist_SR->GetBinContent(i+1);
        double bkg_err = Hist_Bkg->GetBinError(i+1);
        double data_err = Hist_SR->GetBinError(i+1);
        if (Hist_Bkg->GetBinCenter(i+1)<lower_end || Hist_Bkg->GetBinCenter(i+1)>upper_end) {
            continue;
        }
        if ((data_err*data_err+bkg_err*bkg_err)==0) continue;
        chi2_temp += pow(bkg-data,2)/(data_err*data_err+bkg_err*bkg_err);
    }
    chi2_temp = 1./(chi2_temp);
    return chi2_temp;
}
double GetChi2(TH1* Hist_CR, TH1* Hist_SR, TH1* Hist_Bkg, bool includeSR) {
    double inflation = 1.;
    if (!includeSR)
    {
        if (Hist_CR->GetRMS()<Hist_Bkg->GetRMS()) inflation = inflation*10000.;
        if (Hist_CR->GetMean()<Hist_Bkg->GetMean()) inflation = inflation*10000.;
        if (Hist_SR->GetRMS()>Hist_Bkg->GetRMS()) inflation = inflation*10000.;
        if (Hist_SR->GetMean()>Hist_Bkg->GetMean()) inflation = inflation*10000.;
    }
    double chi2_temp = 0.;
    double sign_chi2_temp = 0.;
    double mean = Hist_SR->GetMean();
    double rms = Hist_SR->GetRMS();
    double nbins = 0.;
    double sign_nbins = 0.;
    for (int i=0;i<Hist_SR->GetNbinsX();i++) {
        double bkg = Hist_Bkg->GetBinContent(i+1);
        double data = Hist_SR->GetBinContent(i+1);
        double bkg_err = Hist_Bkg->GetBinError(i+1);
        double data_err = Hist_SR->GetBinError(i+1);
        if (!includeSR && Hist_Bkg->GetBinCenter(i+1)>MSCW_cut_lower && Hist_Bkg->GetBinCenter(i+1)<MSCW_cut_blind) {
            continue;
        }
        if ((data_err*data_err+bkg_err*bkg_err)==0) continue;
        chi2_temp += inflation*pow(bkg-data,2)/(data_err*data_err+bkg_err*bkg_err);
        nbins += 1.;
        if (Hist_Bkg->GetBinCenter(i+1)>(MSCW_cut_blind+rms)) continue;
        sign_chi2_temp += (bkg-data)/pow(data_err*data_err+bkg_err*bkg_err,0.5);
        sign_nbins += 1.;
    }
    //chi2_temp = 1./(chi2_temp/nbins+pow(sign_chi2_temp/sign_nbins,2));
    chi2_temp = 1./(chi2_temp);
    return chi2_temp;
}
void AddBkgStatistics(TH1* Hist_Bkg)
{
    for (int i=0;i<Hist_Bkg->GetNbinsX();i++)
    {
        double stat_err = pow(Hist_Bkg->GetBinContent(i+1),0.5);
        double old_err = Hist_Bkg->GetBinError(i+1);
        Hist_Bkg->SetBinError(i+1,pow(stat_err*stat_err+old_err*old_err,0.5));
    }
}
void AddSystematics3(TH1* Hist_Bkg)
{
    double rel_syst = 0.;
    double area_blind = 0.;
    double area_unblind = 0.;
    for (int i=0;i<Hist_Bkg->GetNbinsX();i++)
    {
        if (Hist_Bkg->GetBinCenter(i+1)>MSCW_cut_blind) {
            area_unblind += pow(Hist_Bkg->GetBinContent(i+1),2);
        }
        else
        {
            area_blind += pow(Hist_Bkg->GetBinContent(i+1),2);
        }
    }
    rel_syst = pow(area_blind/(area_blind+area_unblind),0.5);
    for (int i=0;i<Hist_Bkg->GetNbinsX();i++)
    {
        if (Hist_Bkg->GetBinCenter(i+1)>MSCW_cut_blind) {
            continue;
        }
        double old_err = Hist_Bkg->GetBinError(i+1);
        double new_err = Hist_Bkg->GetBinContent(i+1)*rel_syst;
        if (Hist_Bkg->GetBinCenter(i+1)>MSCW_cut_lower && Hist_Bkg->GetBinCenter(i+1)<MSCW_cut_blind)
        {
            Hist_Bkg->SetBinError(i+1,pow(old_err*old_err+new_err*new_err,0.5));
        }
    }
}
void AddSystematics2(TH1* Hist_CR,TH1* Hist_CRBkg, TH1* Hist_Bkg)
{
    double mean = Hist_CR->GetMean();
    double rms = Hist_CR->GetRMS();
    for (int i=0;i<Hist_CR->GetNbinsX();i++)
    {
        if (Hist_CR->GetBinCenter(i+1)>MSCW_cut_blind) {
            continue;
        }
        double rel_syst = 0.;
        if (Hist_CRBkg->GetBinContent(i+1)!=0)
        {
            rel_syst = abs(Hist_CR->GetBinContent(i+1)-Hist_CRBkg->GetBinContent(i+1))/Hist_CRBkg->GetBinContent(i+1);
        }
        else
        {
            rel_syst = 1.;
        }
        double old_content = Hist_Bkg->GetBinContent(i+1);
        double new_content = old_content+old_content*rel_syst;
        double old_err = Hist_Bkg->GetBinError(i+1);
        double new_err = Hist_Bkg->GetBinContent(i+1)*rel_syst;
        if (Hist_Bkg->GetBinCenter(i+1)>MSCW_cut_lower && Hist_Bkg->GetBinCenter(i+1)<MSCW_cut_blind)
        {
            //Hist_Bkg->SetBinContent(i+1,new_content);
            Hist_Bkg->SetBinError(i+1,pow(old_err*old_err+new_err*new_err,0.5));
        }
    }
}
void AddSystematics(TH1* Hist_SR,TH1* Hist_Bkg)
{
    double mean = Hist_SR->GetMean();
    double rms = Hist_SR->GetRMS();
    double Nbins = 0.;
    double rel_syst = 0.;
    for (int i=0;i<Hist_SR->GetNbinsX();i++)
    {
        if (Hist_SR->GetBinCenter(i+1)>mean+5.*rms) continue;
        if (Hist_SR->GetBinCenter(i+1)<MSCW_cut_blind) {
            continue;
        }
        Nbins += 1.;
        if (Hist_Bkg->GetBinContent(i+1)!=0)
        {
            rel_syst += abs(Hist_SR->GetBinContent(i+1)-Hist_Bkg->GetBinContent(i+1))/Hist_Bkg->GetBinError(i+1);
        }
        else
        {
            rel_syst += 1.;
        }
    }
    rel_syst = rel_syst/Nbins;
    for (int i=0;i<Hist_SR->GetNbinsX();i++)
    {
        double old_err = Hist_Bkg->GetBinError(i+1);
        double new_err = Hist_Bkg->GetBinError(i+1)*rel_syst;
        if (Hist_Bkg->GetBinCenter(i+1)>MSCW_cut_lower && Hist_Bkg->GetBinCenter(i+1)<MSCW_cut_blind)
        {
            Hist_Bkg->SetBinError(i+1,pow(old_err*old_err+new_err*new_err,0.5));
        }
    }
}
void MakeBkgPrevious(TH1* Hist_SR,TH1* Hist_Bkg,TH1* Hist_Previous)
{
    for (int i=0;i<Hist_SR->GetNbinsX();i++)
    {
        if (Hist_Bkg->GetBinCenter(i+1)>MSCW_cut_lower && Hist_Bkg->GetBinCenter(i+1)<MSCW_cut_blind)
        {
            Hist_Previous->SetBinContent(i+1,Hist_Bkg->GetBinContent(i+1));
            Hist_Previous->SetBinError(i+1,Hist_Bkg->GetBinError(i+1));
        }
        else
        {
            Hist_Previous->SetBinContent(i+1,Hist_SR->GetBinContent(i+1));
            Hist_Previous->SetBinError(i+1,Hist_SR->GetBinError(i+1));
        }
        //Hist_Previous->SetBinContent(i+1,Hist_Bkg->GetBinContent(i+1));
        //Hist_Previous->SetBinError(i+1,Hist_Bkg->GetBinError(i+1));
    }
}
double ConvergeFunction(double x, double threshold, double amplitude)
{
    //return 1.;
    //
    //if (x-(threshold)>=0.) return 1.;
    //return 0.;
    //return 1./(1.+exp(-1.*(x-(threshold))/amplitude));

    //if (x-(threshold+2.*amplitude)>=0.) return 1.;
    //if (x-(threshold)<0.) return 0.;
    //return 0.5+0.5*(x-(threshold+amplitude))/amplitude;

    if (x-(threshold+amplitude)>=0.) return 1.;
    if (x-(threshold)<0.) return 0.;
    return pow((x-threshold)/amplitude,2);

    //if (x-(threshold+2.*amplitude)>=0.) return 1.;
    //if (x-(threshold)<0.) return 0.;
    //return exp((x-(threshold+2.*amplitude))/amplitude);

}
double FindConvergeThreshold(TH1* Hist_SR, TH1* Hist_Bkg, TH1* Hist_Bkg_Temp, double amplitude)
{
    double threshold = 0.;
    double chi2_best = 0.;
    double norm_best = 0.;
    double mean = Hist_SR->GetMean();
    double rms = Hist_SR->GetRMS();

    for (int th=0;th<50;th++)
    {
#ifdef VEGAS
        double try_threshold = MSCW_cut_lower + (MSCW_cut_blind-MSCW_cut_lower)*double(th)/50.;
#else
        double try_threshold = MSCW_cut_lower + (MSCW_cut_blind-MSCW_cut_lower)*double(th)/50.;
#endif
        double chi2 = 0.;
        double norm = 0.;
        for (int i=0;i<Hist_Bkg->GetNbinsX();i++)
        {
            double old_content = Hist_Bkg->GetBinContent(i+1);
            double new_content = old_content*ConvergeFunction(Hist_Bkg->GetBinCenter(i+1),try_threshold,amplitude);
            double old_error = Hist_Bkg->GetBinError(i+1);
            double new_error = old_error*ConvergeFunction(Hist_Bkg->GetBinCenter(i+1),try_threshold,amplitude);
            Hist_Bkg_Temp->SetBinContent(i+1,new_content);
        }
        //chi2 = GetChi2WithRange(Hist_SR, Hist_Bkg_Temp,MSCW_cut_lower,MSCW_cut_blind);
        //if (chi2_best<chi2) {
        //    chi2_best = chi2;
        //    threshold = try_threshold;
        //} 
        double data_counts = 0.;
        for (int i=0;i<Hist_Bkg->GetNbinsX();i++)
        {
            if (Hist_Bkg->GetBinCenter(i+1)<MSCW_cut_lower) continue;
            if (Hist_Bkg->GetBinCenter(i+1)>MSCW_cut_blind) continue;
            data_counts += Hist_SR->GetBinContent(i+1);
        }
        double bkg_counts = 0.;
        for (int i=0;i<Hist_Bkg->GetNbinsX();i++)
        {
            if (Hist_Bkg->GetBinCenter(i+1)<MSCW_cut_lower) continue;
            if (Hist_Bkg->GetBinCenter(i+1)>MSCW_cut_blind) continue;
            bkg_counts += Hist_Bkg_Temp->GetBinContent(i+1);
        }
        if (bkg_counts-data_counts==0) continue;
        norm = 1./pow(bkg_counts-data_counts,2);
        if (norm_best<norm) {
            norm_best = norm;
            threshold = try_threshold;
        } 
    }
    return threshold;
}
double FindConvergeAmplitude(TH1* Hist_SR, TH1* Hist_Bkg, TH1* Hist_Bkg_Temp, double threshold, double init_amplitude)
{
    double amplitude = 1.0;
    double chi2_best = 0.;
    double norm_best = 0.;
    double mean = Hist_SR->GetMean();
    double rms = Hist_SR->GetRMS();

    double data_counts = 0.;
    double nbin_counts = 0.;
    for (int i=0;i<Hist_Bkg->GetNbinsX();i++)
    {
        if (Hist_Bkg->GetBinCenter(i+1)<MSCW_cut_lower) continue;
        if (Hist_Bkg->GetBinCenter(i+1)>mean) continue;
        data_counts += Hist_SR->GetBinContent(i+1);
        nbin_counts += 1;
    }
    std::cout << "data statisitcs used in finding converge amplitude = " << data_counts << std::endl;
    std::cout << "N bins used in finding converge amplitude = " << nbin_counts << std::endl;
    std::cout << "Better to keep n_data/n_bin > 10. (if it is <10, merge more runs)" << std::endl;
    for (int amp=0;amp<1000;amp++)
    {
        //double try_amplitude = 0. + 2.0*double(amp)/1000.;
        double try_amplitude = 0. + (mean-MSCW_cut_lower)*double(amp)/1000.;
        double norm = 0.;
        double chi2 = 0.;
        for (int i=0;i<Hist_Bkg->GetNbinsX();i++)
        {
            double old_content = Hist_Bkg->GetBinContent(i+1);
            double new_content = old_content*ConvergeFunction(Hist_Bkg->GetBinCenter(i+1),threshold,try_amplitude);
            double old_error = Hist_Bkg->GetBinError(i+1);
            double new_error = old_error*ConvergeFunction(Hist_Bkg->GetBinCenter(i+1),threshold,try_amplitude);
            Hist_Bkg_Temp->SetBinContent(i+1,new_content);
        }
        chi2 = GetChi2WithRange(Hist_SR, Hist_Bkg_Temp,MSCW_cut_lower,mean);
        if (chi2_best<chi2) {
            chi2_best = chi2;
            amplitude = try_amplitude;
        } 
        //double data_counts = 0.;
        //for (int i=0;i<Hist_Bkg->GetNbinsX();i++)
        //{
        //    if (Hist_Bkg->GetBinCenter(i+1)<MSCW_cut_lower) continue;
        //    if (Hist_Bkg->GetBinCenter(i+1)>MSCW_cut_blind) continue;
        //    data_counts += Hist_SR->GetBinContent(i+1);
        //}
        //double bkg_counts = 0.;
        //for (int i=0;i<Hist_Bkg->GetNbinsX();i++)
        //{
        //    if (Hist_Bkg->GetBinCenter(i+1)<MSCW_cut_lower) continue;
        //    if (Hist_Bkg->GetBinCenter(i+1)>MSCW_cut_blind) continue;
        //    bkg_counts += Hist_Bkg_Temp->GetBinContent(i+1);
        //}
        //if (bkg_counts-data_counts==0) continue;
        //norm = 1./pow(bkg_counts-data_counts,2);
        //if (norm_best<norm) {
        //    norm_best = norm;
        //    amplitude = try_amplitude;
        //} 
    }
    return amplitude;
}
std::pair <double,double> FindConverge(TH1* Hist_SR, TH1* Hist_Bkg, TH1* Hist_Bkg_Temp)
{
    double threshold = 0.;
    double amplitude = 1.0;
    double init_amplitude = 1.0;
    double chi2_best = 0.;
    double norm_best = 0.;
    double mean = Hist_Bkg->GetMean();
    double rms = Hist_Bkg->GetRMS();
    if (UseVegas) 
    {
        //init_amplitude = 10.0;
        init_amplitude = 0.1;
    }
    else
    {
        //init_amplitude = 2.0;
        init_amplitude = 0.5;
    }
    // try a few iterations if possible
    //threshold = FindConvergeThreshold(Hist_SR,Hist_Bkg,Hist_Bkg_Temp,init_amplitude);
    //amplitude = FindConvergeAmplitude(Hist_SR,Hist_Bkg,Hist_Bkg_Temp,threshold,init_amplitude);
    //threshold = FindConvergeThreshold(Hist_SR,Hist_Bkg,Hist_Bkg_Temp,amplitude);
    //amplitude = FindConvergeAmplitude(Hist_SR,Hist_Bkg,Hist_Bkg_Temp,threshold,amplitude);
    //threshold = FindConvergeThreshold(Hist_SR,Hist_Bkg,Hist_Bkg_Temp,amplitude);
    //amplitude = FindConvergeAmplitude(Hist_SR,Hist_Bkg,Hist_Bkg_Temp,threshold,amplitude);
    //threshold = mean-1.7*rms;
    threshold = MSCW_cut_lower;
    amplitude = FindConvergeAmplitude(Hist_SR,Hist_Bkg,Hist_Bkg_Temp,threshold,init_amplitude);
    std::cout << "found threshold = " << threshold << ", amplitude = " << amplitude << std::endl;
    return std::make_pair(threshold,amplitude);
}
void Converge(TH1* Hist_Bkg, double threshold, double amplitude)
{
    double mean = Hist_Bkg->GetMean();
    double rms = Hist_Bkg->GetRMS();
    //double threshold = mean-rms_threshold*rms;
    std::cout << "using threshold = " << threshold << ", amplitude = " << amplitude << std::endl;
    for (int i=0;i<Hist_Bkg->GetNbinsX();i++)
    {
        double old_content = Hist_Bkg->GetBinContent(i+1);
        double new_content = old_content*ConvergeFunction(Hist_Bkg->GetBinCenter(i+1),threshold,amplitude);
        double old_error = Hist_Bkg->GetBinError(i+1);
        double new_error = old_error*ConvergeFunction(Hist_Bkg->GetBinCenter(i+1),threshold,amplitude);
        Hist_Bkg->SetBinContent(i+1,new_content);
    }
}
std::pair <double,double> ShiftAndNormalize(TH1* Hist_CR, TH1* Hist_SR, TH1* Hist_BkgTemp, TH1* Hist_Bkg, double scaled_shift_begin, bool doShift, bool includeSR) {
    double shift_fit = 0;
    double unblinded_shift_fit = 0;
    double scale_fit = 0;
    double chi2 = 0.;
    double unblinded_chi2 = 0.;
    double mean = Hist_SR->GetMean();
    double rms = Hist_SR->GetRMS();
    int norm_bin_low = 0;
    int norm_bin_up = 0;
    double norm = Hist_BkgTemp->Integral();
    Hist_Bkg->Reset();
    double shift_begin = scaled_shift_begin*Hist_BkgTemp->GetRMS();
    double scaled_shift_final = 0.;
    double scaled_unblinded_shift_final = 0.;
    if (doShift) {
        for (int fit=0;fit<1000;fit++) {
                if (Hist_CR->GetRMS()==0.)
                {
                    continue;
                }
                double shift = shift_begin-1.0*Hist_BkgTemp->GetRMS()+2.0*Hist_BkgTemp->GetRMS()*double(fit)*0.001;
                double inflation;
                //if ((shift-shift_begin)/Hist_BkgTemp->GetRMS()-0.13<0.)
                //{
                //    inflation = exp(-0.5*pow(((shift-shift_begin)/Hist_BkgTemp->GetRMS()-0.13)/0.26,2));
                //}
                //else
                //{
                //    inflation = 1.;
                //}
                if ((shift-shift_begin)/Hist_BkgTemp->GetRMS()-0.<0.)
                {
                    inflation = exp(-0.5*pow(((shift-shift_begin)/Hist_BkgTemp->GetRMS()-0.)/0.26,2));
                }
                else
                {
                    inflation = 1.;
                }
                if (includeSR) inflation = 1.;
                //if (shift<0) continue;
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
                norm_bin_low = Hist_SR->FindBin(MSCW_cut_blind);
                norm_bin_up = Hist_SR->FindBin(mean+8.*rms);
                double SR_area2 = Hist_SR->Integral(norm_bin_low,norm_bin_up);
                double Bkg_area2 = Hist_Bkg->Integral(norm_bin_low,norm_bin_up);
                double scale_begin = (SR_area1+SR_area2)/(Bkg_area1+Bkg_area2);
                //double scale_begin = (SR_area2)/(Bkg_area2);
                double scale = scale_begin;
                Hist_Bkg->Scale(norm/Hist_Bkg->Integral());
                Hist_Bkg->Scale(scale);
                double chi2_temp = inflation*GetChi2(Hist_CR, Hist_SR, Hist_Bkg,includeSR);
                double unblinded_chi2_temp = GetChi2(Hist_CR, Hist_SR, Hist_Bkg,true);
                if (chi2<chi2_temp && Hist_Bkg->Integral()!=0) {
                    chi2 = chi2_temp;
                    shift_fit = shift;
                    scale_fit = scale;
                } 
                if (unblinded_chi2<unblinded_chi2_temp && Hist_Bkg->Integral()!=0) {
                    unblinded_chi2 = unblinded_chi2_temp;
                    unblinded_shift_fit = shift;
                } 
        }
    }
    else {
        shift_fit = shift_begin-0.13*Hist_BkgTemp->GetRMS();
        unblinded_shift_fit = shift_begin-0.13*Hist_BkgTemp->GetRMS();
    }
    for (int i=0;i<Hist_SR->GetNbinsX();i++) {
            int b = Hist_Bkg->FindBin(Hist_Bkg->GetBinCenter(i+1)-shift_fit);
            Hist_Bkg->SetBinContent(i+1,Hist_BkgTemp->GetBinContent(b));
            Hist_Bkg->SetBinError(i+1,Hist_BkgTemp->GetBinError(b));
    }
    if (Hist_Bkg->Integral()==0) 
    {
        std::cout << "Hist_Bkg->Integral() = 0!!!" << std::endl;
        Hist_Bkg->Scale(0.);
    }
    else 
    {
        Hist_Bkg->Scale(norm/Hist_Bkg->Integral());
    }
    Hist_Bkg->Scale(scale_fit);
    scaled_shift_final = shift_fit/Hist_Bkg->GetRMS();
    scaled_unblinded_shift_final = unblinded_shift_fit/Hist_Bkg->GetRMS();
    return std::make_pair(scaled_shift_final,scaled_unblinded_shift_final);
}

void Deconvolution(TH1* Hist_source, TH1* Hist_response, TH1* Hist_Deconv, int n_iteration) {
        Hist_Deconv->Reset();
        const int N_bins = Hist_source->GetNbinsX();
        Double_t *source = new Double_t[N_bins];
        Double_t *response = new Double_t[N_bins];
        for (int i=0;i<N_bins;i++) {
            source[i] = max(Hist_source->GetBinContent(i+1),0.);
            //if (Hist_source->GetBinCenter(i+1)<MSCW_cut_lower) source[i] = 0.;
            if (i<N_bins/2) 
            {
                response[i] = max(Hist_response->GetBinContent(i+1+N_bins/2),0.);
                //if (Hist_source->GetBinCenter(i+1+N_bins/2)<MSCW_cut_lower) source[i] = 0.;
            }
            else response[i] = 0;
        }
        TSpectrum sp;
        sp.Deconvolution(source,response,N_bins,n_iteration,1,100000); // new best option
        for (int i=0;i<N_bins;i++) {
            Hist_Deconv->SetBinContent(i+1,max(source[i],0.));
        }
}

double FindNIteration(TH1* Hist_SR, TH1* Hist_CR, TH1* Hist_Bkg, TH1* Hist_BkgTemp, TH1* Hist_Deconv, double scaled_shift, double rms, double mean, int n_iter_begin, bool includeSR) {
    if (!includeSR) return n_iter_begin;
    if (n_iter_begin==1) return 1;
    double n_iter_final = n_iter_begin;
    double chi2_best = 0.;
    chi2_best = 0.;
    TF1 *func = new TF1("func",Kernel,-50.,50.,1);
    func->SetParameter(0,0.5);
    int nbins = Hist_SR->GetNbinsX();
    for (int delta_n_iter = 0;delta_n_iter<=10;delta_n_iter++) {
          int n_iter = n_iter_begin;
          n_iter = n_iter_begin-delta_n_iter;
          if (n_iter<=1) {
              n_iter = 1;
              continue;
          }
          double chi2 = 0;
          func->SetParameter(0,rms);
          Hist_Deconv->Reset();
          Hist_Deconv->FillRandom("func",1000000);
          Deconvolution(Hist_CR,Hist_Deconv,Hist_BkgTemp,n_iter);
          std::pair <double,double> offset = ShiftAndNormalize(Hist_CR,Hist_SR,Hist_BkgTemp,Hist_Bkg,scaled_shift,true,includeSR);
          chi2 = GetChi2(Hist_CR,Hist_SR,Hist_Bkg,includeSR);
          if (chi2_best<chi2) {
              chi2_best = chi2;
              n_iter_final = n_iter;
          } 
    }
    return n_iter_final;
}
std::pair <double,double> FindRMS(TH1* Hist_SR, TH1* Hist_CR, TH1* Hist_Bkg, TH1* Hist_BkgTemp, TH1* Hist_Deconv, double scaled_shift, double rms_begin, double mean, double n_iter, bool includeSR) {
    double chi2_best = 0.;
    double rms_final = rms_begin;
    double unblinded_chi2_best = 0.;
    double unblinded_rms_final = rms_begin;
    double blinded_sr_rms = GetBlindedRMS(Hist_SR);
    double blinded_cr_rms = GetBlindedRMS(Hist_CR);
    double estimated_kernel_rms = 0;
    if (blinded_cr_rms>blinded_sr_rms)
    {
        estimated_kernel_rms = pow(pow(blinded_cr_rms,2)-pow(blinded_sr_rms,2),0.5);
    }
    TF1 *func = new TF1("func",Kernel,-50.,50.,1);
    func->SetParameter(0,0.5);
    for (int n_rms = 0; n_rms<=50;n_rms++) {
        double chi2 = 0;
        double unblinded_chi2 = 0;
        double rms = rms_begin;
        rms = rms_begin-0.5*rms_begin+double(n_rms)*2.5*rms_begin/50.;
        if (!includeSR && estimated_kernel_rms!=0.)
        {
            double rms_diff = (rms-estimated_kernel_rms)/estimated_kernel_rms;
            if (rms_diff<0.) continue;
            if (rms_diff>1.) continue;
        }
        double inflation;
        if ((rms-rms_begin)/rms_begin-0.<0.)
        {
            inflation = exp(-0.5*pow(((rms-rms_begin)/rms_begin-0.)/0.44,2));
        }
        else
        {
            inflation = 1.;
        }
        if (includeSR) inflation = 1.;
        func->SetParameter(0,rms);
        Hist_Deconv->Reset();
        //Hist_Deconv->FillRandom("func",Hist_SR->Integral()*100);
        //Hist_Deconv->FillRandom("func",Hist_SR->Integral()*10);
        Hist_Deconv->FillRandom("func",1000000);
        Deconvolution(Hist_CR,Hist_Deconv,Hist_BkgTemp,n_iter);
        std::pair <double,double> offset = ShiftAndNormalize(Hist_CR,Hist_SR,Hist_BkgTemp,Hist_Bkg,scaled_shift,true,includeSR);
        chi2 = inflation*GetChi2(Hist_CR, Hist_SR, Hist_Bkg,includeSR);
        unblinded_chi2 = GetChi2(Hist_CR, Hist_SR, Hist_Bkg,true);
        if (chi2_best<chi2) {
            chi2_best = chi2;
            rms_final = rms;
        } 
        if (unblinded_chi2_best<unblinded_chi2) {
            unblinded_chi2_best = unblinded_chi2;
            unblinded_rms_final = rms;
        } 
    }
    return std::make_pair(rms_final,unblinded_rms_final);
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

double PointingDistribution(string file_name,int run, bool isDark, bool fillHist)
{
    char run_number[50];
    sprintf(run_number, "%i", int(run));
    TFile*  input_file = TFile::Open(file_name.c_str());
    TTree* pointing_tree = nullptr;
    if (!UseVegas) 
    {
        pointing_tree = (TTree*) input_file->Get("run_"+TString(run_number)+"/stereo/pointingDataReduced");
        pointing_tree->SetBranchAddress("TelElevation",&TelElevation);
        pointing_tree->SetBranchAddress("TelAzimuth",&TelAzimuth);
        pointing_tree->SetBranchAddress("TelRAJ2000",&TelRAJ2000);
        pointing_tree->SetBranchAddress("TelDecJ2000",&TelDecJ2000);
        double total_entries = (double)pointing_tree->GetEntries();
        pointing_tree->GetEntry(int(total_entries/2.));
        if (fillHist)
        {
                if (isDark)
                {
                    Hist_Dark_TelElev.Fill(TelElevation);
                }
                else
                {
                    Hist_Target_TelElev.Fill(TelElevation);
                }
        }
    }
    input_file->Close();
    return TelElevation;
}
bool PointingSelection(string file_name,int run, bool fillHist, double Elev_cut_lower, double Elev_cut_upper, double Azim_cut_lower, double Azim_cut_upper)
{
    char run_number[50];
    sprintf(run_number, "%i", int(run));
    TFile*  input_file = TFile::Open(file_name.c_str());
    TTree* pointing_tree = nullptr;
    if (!UseVegas) 
    {
        pointing_tree = (TTree*) input_file->Get("run_"+TString(run_number)+"/stereo/pointingDataReduced");
        pointing_tree->SetBranchAddress("TelElevation",&TelElevation);
        pointing_tree->SetBranchAddress("TelAzimuth",&TelAzimuth);
        pointing_tree->SetBranchAddress("TelRAJ2000",&TelRAJ2000);
        pointing_tree->SetBranchAddress("TelDecJ2000",&TelDecJ2000);
        double total_entries = (double)pointing_tree->GetEntries();
        pointing_tree->GetEntry(int(total_entries/2.));
        if (TelElevation<Elev_cut_lower) 
        {
            input_file->Close();
            return false;
        }
        if (TelElevation>Elev_cut_upper)
        {
            input_file->Close();
            return false;
        }
        if (TelAzimuth<Azim_cut_lower)
        {
            input_file->Close();
            return false;
        }
        if (TelAzimuth>Azim_cut_upper)
        {
            input_file->Close();
            return false;
        }
        if (fillHist)
        {
            //for (int entry=0;entry<pointing_tree->GetEntries();entry++) {
            //        pointing_tree->GetEntry(entry);
                    Hist_Target_TelElevAzim.Fill(TelElevation,TelAzimuth);
                    Hist_Target_TelRaDec.Fill(TelRAJ2000,TelDecJ2000);
            //}
        }
    }
    else
    {
        pointing_tree = (TTree*) input_file->Get("ShowerEvents/ShowerEventsTree");
        // VEGAS
#ifdef VEGAS
	VAShowerData* sh = nullptr;
	pointing_tree->SetBranchAddress("S", &sh);
        //for (int entry=0;entry<pointing_tree->GetEntries();entry++) {
        //for (int entry=0;entry<1;entry++) {
        //        pointing_tree->GetEntry(entry);
        //        TelElevation = sh->fArrayTrackingElevation_Deg;
        //        TelAzimuth = sh->fArrayTrackingAzimuth_Deg;
        //        Hist_Target_TelElevAzim.Fill(TelElevation,TelAzimuth);
        //        Hist_Target_TelRaDec.Fill(TelRAJ2000,TelDecJ2000);
        //}
        double total_entries = (double)pointing_tree->GetEntries();
        pointing_tree->GetEntry(int(total_entries/2.));
        TelElevation = sh->fArrayTrackingElevation_Deg;
        TelAzimuth = sh->fArrayTrackingAzimuth_Deg;
        if (TelElevation<Elev_cut_lower) 
        {
            input_file->Close();
            return false;
        }
        if (TelElevation>Elev_cut_upper)
        {
            input_file->Close();
            return false;
        }
        if (TelAzimuth<Azim_cut_lower)
        {
            input_file->Close();
            return false;
        }
        if (TelAzimuth>Azim_cut_upper)
        {
            input_file->Close();
            return false;
        }
        //for (int entry=0;entry<pointing_tree->GetEntries();entry++) {
        //        pointing_tree->GetEntry(entry);
        //        TelElevation = sh->fArrayTrackingElevation_Deg;
        //        TelAzimuth = sh->fArrayTrackingAzimuth_Deg;
        //        Hist_Target_TelRaDec_AfterCut.Fill(TelRAJ2000,TelDecJ2000);
        //}
#endif
    }
    input_file->Close();
    return true;
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
vector<int> SortRunListByElevation(string source,vector<int> Target_runlist)
{
        char observation[50];
        vector< pair <double,int> > list_elev_run;
        sprintf(observation, "%s", source.c_str());
        if (TString(source)=="CrabA") sprintf(observation, "%s", "Crab");
        if (TString(source)=="CrabB") sprintf(observation, "%s", "Crab");
        if (TString(source)=="Segue1AV6") sprintf(observation, "%s", "Segue1V6");
        if (TString(source)=="Segue1BV6") sprintf(observation, "%s", "Segue1V6");
        for (int run=0;run<Target_runlist.size();run++) {
            char run_number[50];
            sprintf(run_number, "%i", int(Target_runlist[run]));
            filename = TString("$VERITAS_USER_DATA_DIR/"+TString(observation)+"_V6_Moderate-TMVA-BDT.RB."+TString(run_number)+".root");
            if (TString(observation)=="Coma" || TString(observation)=="BrandonValidation") {
              filename = TString("/veritas/userspace/brandon/VERITAS/Background/anasum/"+TString(run_number)+".anasum.root");
            }
            if (TString(observation)=="VA_Segue1" || TString(observation)=="VA_Geminga")
            {
              filename = TString("/veritas/upload/Geminga_Iowa/"+TString(run_number)+".stage4.root");
            }
            TFile*  input_file = TFile::Open(filename.c_str());
            TTree* pointing_tree = nullptr;
            if (!UseVegas) 
            {
                pointing_tree = (TTree*) input_file->Get("run_"+TString(run_number)+"/stereo/pointingDataReduced");
                pointing_tree->SetBranchAddress("TelElevation",&TelElevation);
                pointing_tree->SetBranchAddress("TelAzimuth",&TelAzimuth);
                pointing_tree->SetBranchAddress("TelRAJ2000",&TelRAJ2000);
                pointing_tree->SetBranchAddress("TelDecJ2000",&TelDecJ2000);
                double total_entries = (double)pointing_tree->GetEntries();
                pointing_tree->GetEntry(int(total_entries/2.));
            }
            else
            {
                pointing_tree = (TTree*) input_file->Get("ShowerEvents/ShowerEventsTree");
                // VEGAS
#ifdef VEGAS
	        VAShowerData* sh = nullptr;
	        pointing_tree->SetBranchAddress("S", &sh);
                double total_entries = (double)pointing_tree->GetEntries();
                pointing_tree->GetEntry(int(total_entries/2.));
                TelElevation = sh->fArrayTrackingElevation_Deg;
                TelAzimuth = sh->fArrayTrackingAzimuth_Deg;
#endif
            }
            double elev_primary = TelElevation;
            double azim_primary = TelAzimuth;
            list_elev_run.push_back(make_pair(elev_primary,Target_runlist[run]));
            input_file->Close();
        }
        std::sort(list_elev_run.begin(), list_elev_run.end());
        vector<int> new_list;
        for (int run=0;run<Target_runlist.size();run++)
        {
            new_list.push_back(list_elev_run.at(run).second);
        }
        return new_list;

}
vector<vector<int>> FindRunSublist(string source, vector<int> Target_runlist)
{
        std::cout << "Getting sublist runs from " << source << std::endl;
        vector<vector<int>> list;
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
        if (TString(source)=="CrabA") sprintf(observation, "%s", "Crab");
        if (TString(source)=="CrabB") sprintf(observation, "%s", "Crab");
        if (TString(source)=="Segue1AV6") sprintf(observation, "%s", "Segue1V6");
        if (TString(source)=="Segue1BV6") sprintf(observation, "%s", "Segue1V6");
        double delta_elev = 5.;
        double delta_azim = 20.;
        //double delta_elev = 10.;
        //double delta_azim = 20.;
        for (int elev=0;elev<(90./delta_elev);elev++)
        {
            double elev_low = delta_elev*((double)elev);
            double elev_up = delta_elev*((double)elev+1);
            if (elev_up<=Target_Elev_cut_lower) continue;
            if (elev_low>=Target_Elev_cut_upper) continue;

            for (int azim=0;azim<(360./delta_azim);azim++)
            {
                double azim_low = delta_azim*((double)azim);
                double azim_up = delta_azim*((double)azim+1);
                if (azim_up<=Target_Azim_cut_lower) continue;
                if (azim_low>=Target_Azim_cut_upper) continue;
                std::cout << "run elev. " << elev_low << "-" << elev_up << ", azim. " << azim_low << "-" << azim_up << std::endl;
                vector<int> sublist;
                for (int run=0;run<Target_runlist.size();run++) {
                    char run_number[50];
                    sprintf(run_number, "%i", int(Target_runlist[run]));
                    filename = TString("$VERITAS_USER_DATA_DIR/"+TString(observation)+"_V6_Moderate-TMVA-BDT.RB."+TString(run_number)+".root");
                    if (TString(observation)=="Coma" || TString(observation)=="BrandonValidation") {
                      filename = TString("/veritas/userspace/brandon/VERITAS/Background/anasum/"+TString(run_number)+".anasum.root");
                    }
                    if (TString(observation)=="VA_Segue1" || TString(observation)=="VA_Geminga")
                    {
                      filename = TString("/veritas/upload/Geminga_Iowa/"+TString(run_number)+".stage4.root");
                    }
                    TFile*  input_file = TFile::Open(filename.c_str());
                    TTree* pointing_tree = nullptr;
                    if (!UseVegas) 
                    {
                        pointing_tree = (TTree*) input_file->Get("run_"+TString(run_number)+"/stereo/pointingDataReduced");
                        pointing_tree->SetBranchAddress("TelElevation",&TelElevation);
                        pointing_tree->SetBranchAddress("TelAzimuth",&TelAzimuth);
                        pointing_tree->SetBranchAddress("TelRAJ2000",&TelRAJ2000);
                        pointing_tree->SetBranchAddress("TelDecJ2000",&TelDecJ2000);
                        double total_entries = (double)pointing_tree->GetEntries();
                        pointing_tree->GetEntry(int(total_entries/2.));
                        if (TelElevation<elev_low) 
                        {
                            input_file->Close();
                            continue;
                        }
                        if (TelElevation>elev_up)
                        {
                            input_file->Close();
                            continue;
                        }
                        if (TelAzimuth<azim_low)
                        {
                            input_file->Close();
                            continue;
                        }
                        if (TelAzimuth>azim_up)
                        {
                            input_file->Close();
                            continue;
                        }
                    }
                    else
                    {
                        pointing_tree = (TTree*) input_file->Get("ShowerEvents/ShowerEventsTree");
                        // VEGAS
#ifdef VEGAS
	                VAShowerData* sh = nullptr;
	                pointing_tree->SetBranchAddress("S", &sh);
                        double total_entries = (double)pointing_tree->GetEntries();
                        pointing_tree->GetEntry(int(total_entries/2.));
                        TelElevation = sh->fArrayTrackingElevation_Deg;
                        TelAzimuth = sh->fArrayTrackingAzimuth_Deg;
#endif
                        if (TelElevation<elev_low) 
                        {
                            input_file->Close();
                            continue;
                        }
                        if (TelElevation>elev_up)
                        {
                            input_file->Close();
                            continue;
                        }
                        if (TelAzimuth<azim_low)
                        {
                            input_file->Close();
                            continue;
                        }
                        if (TelAzimuth>azim_up)
                        {
                            input_file->Close();
                            continue;
                        }
                    }
                    sublist.push_back(Target_runlist[run]);
                    input_file->Close();
                    continue;
                }
                if (!(sublist.size()==0)) list.push_back(sublist);
            }
        }
        return list;
}
void DeconvolutionMethodForExtendedSources(string target_data, int NTelMin, int NTelMax, double elev_lower, double elev_upper, double azim_lower, double azim_upper, double theta2_cut_lower_input, double theta2_cut_upper_input, double MSCW_cut_blind_input, double MSCW_cut_upper_input, double MSCW_cut_lower_input,bool DoConverge) 
{

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
        MSCW_cut_lower = MSCW_cut_lower_input;
        MSCW_cut_upper = MSCW_cut_upper_input;
        MSCW_cut_blind = MSCW_cut_blind_input;

#ifndef VEGAS
        Theta2_upper_limit = 10.;
#else
        Theta2_upper_limit = 20.;
#endif

        TRandom rnd;
        TH2D Hist_Target_TelRaDec_AfterCut("Hist_Target_TelRaDec_AfterCut","",100,0,5,100,-1,1);
        TH1D Hist_Target_ON_MSCW_Alpha("Hist_Target_ON_MSCW_Alpha","",100,0,10);
        TH1D Hist_Target_OFF_MSCW_Alpha("Hist_Target_OFF_MSCW_Alpha","",100,0,10);
        vector<int> Run_sublist;
        vector<TH1D> Hist_Target_SR_ErecS;
        vector<vector<TH1D>> Hist_Target_SR_theta2;
        vector<vector<TH1D>> Hist_Dark_SR_theta2;
        vector<vector<TH1D>> Hist_DarkLZA_SR_theta2;
        vector<vector<TH1D>> Hist_DarkSZA_SR_theta2;
        vector<vector<TH1D>> Hist_Dark_CR_theta2;
        vector<vector<TH1D>> Hist_DarkLZA_CR_theta2;
        vector<vector<TH1D>> Hist_DarkSZA_CR_theta2;
        vector<vector<TH1D>> Hist_Target_CR_theta2;
        vector<vector<TH1D>> Hist_TargetLZA_CR_theta2;
        vector<vector<TH1D>> Hist_TargetSZA_CR_theta2;
        vector<vector<TH2D>> Hist_Target_SR_RaDec;
        vector<vector<TH2D>> Hist_Target_CR_RaDec;
        vector<vector<TH1D>> Hist_Target_SR_MSCW;
        vector<vector<TH1D>> Hist_Target_SR_MSCW_SumRuns;
        vector<vector<TH1D>> Hist_Target_CR_MSCW;
        vector<vector<TH1D>> Hist_Target_CR_MSCW_SumRuns;
        vector<TH1D> Hist_Target_SR_MSCW_SumRuns_SumSRs;
        vector<TH1D> Hist_Target_SR_MSCL;
        vector<TH1D> Hist_Target_BkgTemp_MSCW;
        vector<TH1D> Hist_Target_BkgPrevious_MSCW;
        vector<vector<TH1D>> Hist_Target_BkgSR_MSCW;
        vector<vector<TH1D>> Hist_Target_BkgSR_MSCW_AllCR;
        vector<double> Weight_Target_BkgSR_MSCW_AllCR;
        vector<vector<TH1D>> Hist_Target_BkgCR_MSCW;
        vector<vector<TH1D>> Hist_Target_BkgSR_MSCW_SumRuns;
        vector<vector<TH1D>> Hist_Target_BkgCR_MSCW_SumRuns;
        vector<TH1D> Hist_Target_BkgSR_MSCW_SumRuns_SumSRs;
        vector<TH1D> Hist_Target_Deconv_MSCW;
        vector<TH2D> Hist_Target_MSCLW;
        for (int e=0;e<=N_energy_bins;e++) {
            energy_vec.push_back(energy_bins[e]);
        }
        for (int e=0;e<N_energy_bins;e++) {
            char e_low[50];
            sprintf(e_low, "%i", int(energy_bins[e]));
            char e_up[50];
            sprintf(e_up, "%i", int(energy_bins[e+1]));
            //if (energy_bins[e]>=100.) N_bins_for_deconv = 1920;
            //if (energy_bins[e]>=350.) N_bins_for_deconv = 960;
            //if (energy_bins[e]>=1000.) N_bins_for_deconv = 480;
            //if (energy_bins[e]>=2000.) N_bins_for_deconv = 240;
            //if (energy_bins[e]>=6000.) N_bins_for_deconv = 120;
            if (energy_bins[e]>=200.) N_bins_for_deconv = 960;
            if (energy_bins[e]>=800.) N_bins_for_deconv = 480;
            if (energy_bins[e]>=3200.) N_bins_for_deconv = 240;
            if (UseVegas)
            {
              if (energy_bins[e]>=100.) N_bins_for_deconv = 960;
              if (energy_bins[e]>=1000.) N_bins_for_deconv = 960;
              if (energy_bins[e]>=3000.) N_bins_for_deconv = 960*2;
            }
            Hist_Target_SR_ErecS.push_back(TH1D("Hist_Target_SR_ErecS_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_energy_bins,energy_bins));
            Hist_Target_SR_MSCL.push_back(TH1D("Hist_Target_SR_MSCL_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
            Hist_Target_BkgTemp_MSCW.push_back(TH1D("Hist_Target_BkgTemp_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
            Hist_Target_BkgPrevious_MSCW.push_back(TH1D("Hist_Target_BkgPrevious_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
            Hist_Target_Deconv_MSCW.push_back(TH1D("Hist_Target_Deconv_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
            if (UseVegas)
            {
            Hist_Target_MSCLW.push_back(TH2D("Hist_Target_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",Number_of_SR+Number_of_CR,MSCL_signal_cut_lower[Number_of_SR-1],MSCL_control_cut_upper[0],100,0,5));
            }
            else
            {
            Hist_Target_MSCLW.push_back(TH2D("Hist_Target_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",Number_of_SR+Number_of_CR,MSCL_signal_cut_lower[Number_of_SR-1],MSCL_control_cut_upper[0],100,-1,15));
            }
            vector<TH1D> Hist_Target_ThisE_SR_MSCW;
            vector<TH1D> Hist_Target_ThisE_SR_MSCW_SumRuns;
            vector<TH1D> Hist_Target_ThisE_BkgSR_MSCW;
            vector<TH1D> Hist_Target_ThisE_BkgSR_MSCW_AllCR;
            vector<TH1D> Hist_Target_ThisE_BkgSR_MSCW_SumRuns;
            vector<TH1D> Hist_Target_ThisE_CR_MSCW;
            vector<TH1D> Hist_Target_ThisE_CR_MSCW_SumRuns;
            vector<TH1D> Hist_Target_ThisE_BkgCR_MSCW;
            vector<TH1D> Hist_Target_ThisE_BkgCR_MSCW_SumRuns;
            vector<TH1D> Hist_Target_ThisE_SR_theta2;
            vector<TH1D> Hist_Dark_ThisE_SR_theta2;
            vector<TH1D> Hist_DarkLZA_ThisE_SR_theta2;
            vector<TH1D> Hist_DarkSZA_ThisE_SR_theta2;
            vector<TH1D> Hist_Dark_ThisE_CR_theta2;
            vector<TH1D> Hist_DarkLZA_ThisE_CR_theta2;
            vector<TH1D> Hist_DarkSZA_ThisE_CR_theta2;
            vector<TH1D> Hist_Target_ThisE_CR_theta2;
            vector<TH1D> Hist_TargetLZA_ThisE_CR_theta2;
            vector<TH1D> Hist_TargetSZA_ThisE_CR_theta2;
            vector<TH2D> Hist_Target_ThisE_SR_RaDec;
            vector<TH2D> Hist_Target_ThisE_CR_RaDec;
            for (int s=0;s<Number_of_SR;s++)
            {
                char nsr[50];
                sprintf(nsr, "%i", int(s+1));
                Hist_Target_ThisE_SR_MSCW.push_back(TH1D("Hist_Target_SR"+TString(nsr)+"_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
                Hist_Target_ThisE_SR_MSCW_SumRuns.push_back(TH1D("Hist_Target_SR"+TString(nsr)+"_MSCW_SumRuns_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
                Hist_Target_ThisE_BkgSR_MSCW.push_back(TH1D("Hist_Target_BkgSR"+TString(nsr)+"_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
                Hist_Target_ThisE_BkgSR_MSCW_AllCR.push_back(TH1D("Hist_Target_BkgSR"+TString(nsr)+"_MSCW_AllCR_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
                Hist_Target_ThisE_BkgSR_MSCW_SumRuns.push_back(TH1D("Hist_Target_BkgSR"+TString(nsr)+"_MSCW_SumRuns_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
                Hist_Target_ThisE_SR_theta2.push_back(TH1D("Hist_Target_SR"+TString(nsr)+"_theta2_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",1024,0,Theta2_upper_limit));
                Hist_Dark_ThisE_SR_theta2.push_back(TH1D("Hist_Dark_SR"+TString(nsr)+"_theta2_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",20,0,Theta2_upper_limit));
                Hist_DarkLZA_ThisE_SR_theta2.push_back(TH1D("Hist_DarkLZA_SR"+TString(nsr)+"_theta2_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",20,0,Theta2_upper_limit));
                Hist_DarkSZA_ThisE_SR_theta2.push_back(TH1D("Hist_DarkSZA_SR"+TString(nsr)+"_theta2_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",20,0,Theta2_upper_limit));
                Hist_Dark_ThisE_CR_theta2.push_back(TH1D("Hist_Dark_CR"+TString(nsr)+"_theta2_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",20,0,Theta2_upper_limit));
                Hist_DarkLZA_ThisE_CR_theta2.push_back(TH1D("Hist_DarkLZA_CR"+TString(nsr)+"_theta2_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",20,0,Theta2_upper_limit));
                Hist_DarkSZA_ThisE_CR_theta2.push_back(TH1D("Hist_DarkSZA_CR"+TString(nsr)+"_theta2_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",20,0,Theta2_upper_limit));
                Hist_Target_ThisE_CR_theta2.push_back(TH1D("Hist_Target_CR"+TString(nsr)+"_theta2_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",1024,0,Theta2_upper_limit));
                Hist_TargetLZA_ThisE_CR_theta2.push_back(TH1D("Hist_TargetLZA_CR"+TString(nsr)+"_theta2_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",1024,0,Theta2_upper_limit));
                Hist_TargetSZA_ThisE_CR_theta2.push_back(TH1D("Hist_TargetSZA_CR"+TString(nsr)+"_theta2_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",1024,0,Theta2_upper_limit));
                Hist_Target_ThisE_SR_RaDec.push_back(TH2D("Hist_Target_SR"+TString(nsr)+"_RaDec_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",100,-2,2,100,-2,2));
                Hist_Target_ThisE_CR_RaDec.push_back(TH2D("Hist_Target_CR"+TString(nsr)+"_RaDec_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",100,-2,2,100,-2,2));
            }
            Hist_Target_SR_MSCW.push_back(Hist_Target_ThisE_SR_MSCW);
            Hist_Target_SR_MSCW_SumRuns.push_back(Hist_Target_ThisE_SR_MSCW_SumRuns);
            Hist_Target_BkgSR_MSCW.push_back(Hist_Target_ThisE_BkgSR_MSCW);
            Hist_Target_BkgSR_MSCW_AllCR.push_back(Hist_Target_ThisE_BkgSR_MSCW_AllCR);
            Hist_Target_BkgSR_MSCW_SumRuns.push_back(Hist_Target_ThisE_BkgSR_MSCW_SumRuns);
            Hist_Target_SR_theta2.push_back(Hist_Target_ThisE_SR_theta2);
            Hist_Dark_SR_theta2.push_back(Hist_Dark_ThisE_SR_theta2);
            Hist_DarkLZA_SR_theta2.push_back(Hist_DarkLZA_ThisE_SR_theta2);
            Hist_DarkSZA_SR_theta2.push_back(Hist_DarkSZA_ThisE_SR_theta2);
            Hist_Dark_CR_theta2.push_back(Hist_Dark_ThisE_CR_theta2);
            Hist_DarkLZA_CR_theta2.push_back(Hist_DarkLZA_ThisE_CR_theta2);
            Hist_DarkSZA_CR_theta2.push_back(Hist_DarkSZA_ThisE_CR_theta2);
            Hist_Target_CR_theta2.push_back(Hist_Target_ThisE_CR_theta2);
            Hist_TargetLZA_CR_theta2.push_back(Hist_TargetLZA_ThisE_CR_theta2);
            Hist_TargetSZA_CR_theta2.push_back(Hist_TargetSZA_ThisE_CR_theta2);
            Hist_Target_SR_RaDec.push_back(Hist_Target_ThisE_SR_RaDec);
            Hist_Target_CR_RaDec.push_back(Hist_Target_ThisE_CR_RaDec);
            Hist_Target_SR_MSCW_SumRuns_SumSRs.push_back(TH1D("Hist_Target_SR_MSCW_SumRuns_SumSRs_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
            Hist_Target_BkgSR_MSCW_SumRuns_SumSRs.push_back(TH1D("Hist_Target_BkgSR_MSCW_SumRuns_SumSRs_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
            for (int s=0;s<Number_of_CR;s++)
            {
                char nCR[50];
                sprintf(nCR, "%i", int(s+1));
                Hist_Target_ThisE_CR_MSCW.push_back(TH1D("Hist_Target_CR"+TString(nCR)+"_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
                Hist_Target_ThisE_CR_MSCW_SumRuns.push_back(TH1D("Hist_Target_CR"+TString(nCR)+"_MSCW_SumRuns_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
                Hist_Target_ThisE_BkgCR_MSCW.push_back(TH1D("Hist_Target_BkgCR"+TString(nCR)+"_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
                Hist_Target_ThisE_BkgCR_MSCW_SumRuns.push_back(TH1D("Hist_Target_BkgCR"+TString(nCR)+"_MSCW_SumRuns_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
            }
            Hist_Target_CR_MSCW.push_back(Hist_Target_ThisE_CR_MSCW);
            Hist_Target_CR_MSCW_SumRuns.push_back(Hist_Target_ThisE_CR_MSCW_SumRuns);
            Hist_Target_BkgCR_MSCW.push_back(Hist_Target_ThisE_BkgCR_MSCW);
            Hist_Target_BkgCR_MSCW_SumRuns.push_back(Hist_Target_ThisE_BkgCR_MSCW_SumRuns);
        }


        std::cout << "Getting dark pointing distributions... " << std::endl;
        vector<int> Dark_runlist = GetRunList("DarkField");
        char Dark_observation[50];
        sprintf(Dark_observation, "%s", "DarkField");
        if (TString(target)=="Segue1V5" || TString(target)=="MGRO_J1908_V5")
        {
            Dark_runlist = GetRunList("Segue1V5");
            sprintf(Dark_observation, "%s", "Segue1V5");
        }
        if (UseVegas) 
        {
            Dark_runlist = GetRunList("VA_DarkSegue1");
            sprintf(Dark_observation, "%s", "VA_DarkSegue1");
        }
        for (int run=0;run<Dark_runlist.size();run++)
        {
            char run_number[50];
            sprintf(run_number, "%i", int(Dark_runlist[run]));
            //std::cout << "Reading run " << run_number << std::endl;
            filename = TString("$VERITAS_USER_DATA_DIR/"+TString(Dark_observation)+"_V6_Moderate-TMVA-BDT.RB."+TString(run_number)+".root");
            if (TString(Dark_observation)=="VA_DarkSegue1")
            {
              filename = TString("/veritas/upload/Geminga_Iowa/"+TString(run_number)+".stage4.root");
            }

            double tel_elev = PointingDistribution(filename,int(Dark_runlist[run]), true, true);
        }

        std::cout << "Getting target pointing distributions... " << std::endl;
        char observation[50];
        sprintf(observation, "%s", target);
        if (TString(target)=="CrabA") sprintf(observation, "%s", "Crab");
        if (TString(target)=="CrabB") sprintf(observation, "%s", "Crab");
        if (TString(target)=="Segue1AV6") sprintf(observation, "%s", "Segue1V6");
        if (TString(target)=="Segue1BV6") sprintf(observation, "%s", "Segue1V6");
        vector<int> Target_runlist = GetRunList(target);
        for (int run=0;run<Target_runlist.size();run++)
        {
            char run_number[50];
            sprintf(run_number, "%i", int(Target_runlist[run]));
            //std::cout << "Reading run " << run_number << std::endl;
            filename = TString("$VERITAS_USER_DATA_DIR/"+TString(observation)+"_V6_Moderate-TMVA-BDT.RB."+TString(run_number)+".root");
            if (TString(observation)=="VA_Segue1" || TString(observation)=="VA_Geminga")
            {
              filename = TString("/veritas/upload/Geminga_Iowa/"+TString(run_number)+".stage4.root");
            }

            double tel_elev = PointingDistribution(filename,int(Target_runlist[run]), false, true);
        }

        std::cout << "Getting dark runs... " << std::endl;
        //for (int e=0;e<N_energy_bins;e++)
        //{
        //    std::cout << "=================================================================" << std::endl;
        //    std::cout << "Dark, e " << energy_bins[e] << std::endl;
        //    for (int s=0;s<Number_of_SR;s++)
        //    {
        //        Hist_Dark_SR_theta2.at(e).at(s).Reset();
        //        Hist_DarkLZA_SR_theta2.at(e).at(s).Reset();
        //        Hist_DarkSZA_SR_theta2.at(e).at(s).Reset();
        //        Hist_Dark_CR_theta2.at(e).at(s).Reset();
        //        Hist_DarkLZA_CR_theta2.at(e).at(s).Reset();
        //        Hist_DarkSZA_CR_theta2.at(e).at(s).Reset();
        //    }
            for (int run=0;run<Dark_runlist.size();run++)
            {
                char run_number[50];
                sprintf(run_number, "%i", int(Dark_runlist[run]));
                //std::cout << "Reading run " << run_number << std::endl;
                filename = TString("$VERITAS_USER_DATA_DIR/"+TString(Dark_observation)+"_V6_Moderate-TMVA-BDT.RB."+TString(run_number)+".root");
                if (TString(Dark_observation)=="VA_DarkSegue1")
                {
                  filename = TString("/veritas/upload/Geminga_Iowa/"+TString(run_number)+".stage4.root");
                }

                //if (!PointingSelection(filename,int(Dark_runlist[run]),false,Target_Elev_cut_lower,Target_Elev_cut_upper,Target_Azim_cut_lower,Target_Azim_cut_upper)) continue;
                bool nominal_dark = false;
                bool LZA_dark = false;
                bool SZA_dark = false;
                if (PointingSelection(filename,int(Dark_runlist[run]),false,50,90,0,360)) nominal_dark = true;
                if (PointingSelection(filename,int(Dark_runlist[run]),false,50,70,0,360)) LZA_dark = true;
                if (PointingSelection(filename,int(Dark_runlist[run]),false,70,90,0,360)) SZA_dark = true;
                double tel_elev = PointingDistribution(filename,int(Dark_runlist[run]), true, false);

                TFile*  input_file = TFile::Open(filename.c_str());
                TString root_file = "run_"+TString(run_number)+"/stereo/data_on";
                if (UseVegas) root_file = "ShowerEvents/ShowerEventsTree";
                TTree* Dark_tree = (TTree*) input_file->Get(root_file);
                // VEGAS
#ifdef VEGAS
	        VAShowerData* sh = nullptr;
	            Dark_tree->SetBranchAddress("S", &sh);
#else
                    Dark_tree->SetBranchAddress("Xoff",&Xoff);
                    Dark_tree->SetBranchAddress("Yoff",&Yoff);
                    Dark_tree->SetBranchAddress("theta2",&theta2);
                    Dark_tree->SetBranchAddress("ra",&ra_sky);
                    Dark_tree->SetBranchAddress("dec",&dec_sky);
                    Dark_tree->SetBranchAddress("ErecS",&ErecS);
                    Dark_tree->SetBranchAddress("EChi2S",&EChi2S);
                    Dark_tree->SetBranchAddress("MSCW",&MSCW);
                    Dark_tree->SetBranchAddress("MSCL",&MSCL);
                    Dark_tree->SetBranchAddress("NImages",&NImages);
#endif
                for (int entry=0;entry<Dark_tree->GetEntries();entry++) {
                    theta2 = 0;
                    ra_sky = 0;
                    dec_sky = 0;
                    ErecS = 0;
                    EChi2S = 0;
                    NImages = 0;
                    MSCW = 0;
                    MSCL = 0;
                    Dark_tree->GetEntry(entry);
                    R2off = Xoff*Xoff+Yoff*Yoff;
                    //theta2 = ra_sky*ra_sky+dec_sky*dec_sky;
                    // VEGAS
#ifdef VEGAS
                    if (UseVegas)
                    {
                        if (sh->fAccurateShowerMaxHeight_KM<7.) continue;
                        ErecS = sh->fEnergy_GeV/1000.;
                        theta2 = sh->fTheta2_Deg2;
                        //MSCW = sh->fMRSW;
                        //MSCL = sh->fMRSL;
                        MSCW = sh->fMSW;
                        MSCL = sh->fMSL;
                        ra_sky = sh->fDirectionRA_J2000_Rad*180./TMath::Pi()-151.77;
                        dec_sky = sh->fDirectionDec_J2000_Rad*180./TMath::Pi()-16.08194;
                    }
#endif
                    int energy = Hist_Target_SR_ErecS.at(0).FindBin(ErecS*1000.)-1;
                    if (energy<0) continue;
                    if (energy>=N_energy_bins) continue;
                    int e = energy;
                    //if (e!=energy) continue;
                    if (!SelectNImages(NTelMin,NTelMax)) continue;
                    if (FoV()) {
                            if (SignalSelectionMSCL()) 
                            {
                                for (int s=0;s<Number_of_SR;s++)
                                {
                                    if (SignalSelectionMSCW(s))
                                    {
                                        if (nominal_dark)
                                        {
                                            int elev_bin = Hist_Dark_TelElev.FindBin(tel_elev);
                                            double n_dark = Hist_Dark_TelElev.GetBinContent(elev_bin);
                                            double n_target = Hist_Target_TelElev.GetBinContent(elev_bin);
                                            double weight = 0;
                                            if (n_dark!=0) weight = n_target/n_dark;
                                            Hist_Dark_SR_theta2.at(e).at(s).Fill(R2off,weight);
                                        }
                                        if (LZA_dark)
                                        {
                                            Hist_DarkLZA_SR_theta2.at(e).at(s).Fill(R2off);
                                        }
                                        if (SZA_dark)
                                        {
                                            Hist_DarkSZA_SR_theta2.at(e).at(s).Fill(R2off);
                                        }
                                    }
                                }
                            }
                            for (int s=0;s<Number_of_SR;s++)
                            {
                                if (ControlSelectionTheta2(s)) 
                                {
                                        if (nominal_dark)
                                        {
                                            int elev_bin = Hist_Dark_TelElev.FindBin(tel_elev);
                                            double n_dark = Hist_Dark_TelElev.GetBinContent(elev_bin);
                                            double n_target = Hist_Target_TelElev.GetBinContent(elev_bin);
                                            double weight = 0;
                                            if (n_dark!=0) weight = n_target/n_dark;
                                            Hist_Dark_CR_theta2.at(e).at(s).Fill(R2off,weight);
                                        }
                                        if (LZA_dark)
                                        {
                                            Hist_DarkLZA_CR_theta2.at(e).at(s).Fill(R2off);
                                        }
                                        if (SZA_dark)
                                        {
                                            Hist_DarkSZA_CR_theta2.at(e).at(s).Fill(R2off);
                                        }
                                }
                            }
                    }
                }
                input_file->Close();
                //break;  // for test purpose
            }
        //}


        std::cout << "Getting target runs... " << std::endl;

        vector<int> N_iter;
        for (int e=0;e<N_energy_bins;e++) {
                N_iter.push_back(30);
        }

        TF1 *myfunc = new TF1("myfunc",Kernel,-50.,50.,1);
        myfunc->SetParameter(0,0.5);

        for (int s=0;s<Number_of_SR;s++)
        {
                Weight_Target_BkgSR_MSCW_AllCR.push_back(0.);
        }

        vector<vector<int>> Sublist;
        Target_runlist = SortRunListByElevation(target,Target_runlist);
        Sublist = FindRunSublist(target,Target_runlist);
        for (int e=0;e<N_energy_bins;e++)
        {
            std::cout << "=================================================================" << std::endl;
            std::cout << "Target, e " << energy_bins[e] << std::endl;
            std::cout << "Sublist.size() = " << Sublist.size() << std::endl;
            for (int subrun=0;subrun<Sublist.size();subrun++)
            {
                std::cout << subrun << "/" << Sublist.size() << " completed." << std::endl;
                for (int s=0;s<Number_of_SR;s++)
                {
                    Hist_Target_SR_MSCW.at(e).at(s).Reset();
                    Hist_Target_BkgSR_MSCW.at(e).at(s).Reset();
                    Hist_Target_BkgSR_MSCW_AllCR.at(e).at(s).Reset();
                }
                for (int s=0;s<Number_of_CR;s++)
                {
                    Hist_Target_CR_MSCW.at(e).at(s).Reset();
                    Hist_Target_BkgCR_MSCW.at(e).at(s).Reset();
                }
                if (Sublist.at(subrun).size()==0) continue;
                for (int run=0;run<Sublist.at(subrun).size();run++)
                {
                    char run_number[50];
                    sprintf(run_number, "%i", int(Sublist.at(subrun)[run]));
                    //std::cout << "Reading run " << run_number << std::endl;
                    filename = TString("$VERITAS_USER_DATA_DIR/"+TString(observation)+"_V6_Moderate-TMVA-BDT.RB."+TString(run_number)+".root");
                    if (TString(observation)=="Coma" || TString(observation)=="BrandonValidation") {
                      filename = TString("/veritas/userspace/brandon/VERITAS/Background/anasum/"+TString(run_number)+".anasum.root");
                    }
                    if (TString(observation)=="VA_Segue1" || TString(observation)=="VA_Geminga")
                    {
                      filename = TString("/veritas/upload/Geminga_Iowa/"+TString(run_number)+".stage4.root");
                    }

                    if (!PointingSelection(filename,int(Sublist.at(subrun)[run]),true,Target_Elev_cut_lower,Target_Elev_cut_upper,Target_Azim_cut_lower,Target_Azim_cut_upper)) continue;

                    TFile*  input_file = TFile::Open(filename.c_str());
                    TString root_file = "run_"+TString(run_number)+"/stereo/data_on";
                    if (UseVegas) root_file = "ShowerEvents/ShowerEventsTree";
                    TTree* Target_tree = (TTree*) input_file->Get(root_file);
                    // VEGAS
#ifdef VEGAS
	            VAShowerData* sh = nullptr;
	                Target_tree->SetBranchAddress("S", &sh);
#else
                        Target_tree->SetBranchAddress("Xoff",&Xoff);
                        Target_tree->SetBranchAddress("Yoff",&Yoff);
                        Target_tree->SetBranchAddress("theta2",&theta2);
                        Target_tree->SetBranchAddress("ra",&ra_sky);
                        Target_tree->SetBranchAddress("dec",&dec_sky);
                        Target_tree->SetBranchAddress("ErecS",&ErecS);
                        Target_tree->SetBranchAddress("EChi2S",&EChi2S);
                        Target_tree->SetBranchAddress("MSCW",&MSCW);
                        Target_tree->SetBranchAddress("MSCL",&MSCL);
                        Target_tree->SetBranchAddress("Time",&Time);
                        Target_tree->SetBranchAddress("NImages",&NImages);
                        //Target_tree->SetBranchAddress("MWR",&MSCW);
                        //Target_tree->SetBranchAddress("MLR",&MSCL);
#endif
                    Target_tree->GetEntry(0);
                    R2off = Xoff*Xoff+Yoff*Yoff;
                    double time_0 = Time;
                    Target_tree->GetEntry(Target_tree->GetEntries()-1);
                    R2off = Xoff*Xoff+Yoff*Yoff;
                    double time_1 = Time;
                    if (e==0) exposure_hours += (time_1-time_0)/3600.;
                    if (e==0) used_runs.push_back(int(Sublist.at(subrun)[run]));
                    for (int entry=0;entry<Target_tree->GetEntries();entry++) {
                        theta2 = 0;
                        ra_sky = 0;
                        dec_sky = 0;
                        ErecS = 0;
                        EChi2S = 0;
                        NImages = 0;
                        MSCW = 0;
                        MSCL = 0;
                        Target_tree->GetEntry(entry);
                        R2off = Xoff*Xoff+Yoff*Yoff;
                        //theta2 = ra_sky*ra_sky+dec_sky*dec_sky;
                        // VEGAS
#ifdef VEGAS
                        if (UseVegas)
                        {
                            if (sh->fAccurateShowerMaxHeight_KM<7.) continue;
                            ErecS = sh->fEnergy_GeV/1000.;
                            theta2 = sh->fTheta2_Deg2;
                            //MSCW = sh->fMRSW;
                            //MSCL = sh->fMRSL;
                            MSCW = sh->fMSW;
                            MSCL = sh->fMSL;
                            ra_sky = sh->fDirectionRA_J2000_Rad*180./TMath::Pi()-151.77;
                            dec_sky = sh->fDirectionDec_J2000_Rad*180./TMath::Pi()-16.08194;
                        }
#endif
                        int energy = Hist_Target_SR_ErecS.at(0).FindBin(ErecS*1000.)-1;
                        if (e!=energy) continue;
                        if (MSCL<MSCL_signal_cut_upper[0]) {
                            if (theta2<0.1) {
                                Hist_Target_ON_MSCW_Alpha.Fill(MSCW);
                            }
                            else {
                                Hist_Target_OFF_MSCW_Alpha.Fill(MSCW);
                            }
                        }
                        if (!SelectNImages(NTelMin,NTelMax)) continue;
                        if (FoV()) {
                                if (SignalSelectionMSCL()) 
                                {
                                    for (int s=0;s<Number_of_SR;s++)
                                    {
                                        if (SignalSelectionMSCW(s))
                                        {
                                            Hist_Target_SR_theta2.at(e).at(s).Fill(theta2);
                                            Hist_Target_SR_RaDec.at(e).at(s).Fill(ra_sky,dec_sky);
                                        }
                                    }
                                }
                                for (int s=0;s<Number_of_SR;s++)
                                {
                                    if (ControlSelectionTheta2(s)) 
                                    {
                                        double weight = 0.;
                                        double dark_cr = (double) Hist_Dark_CR_theta2.at(e).at(s).GetBinContent(Hist_Dark_CR_theta2.at(e).at(s).GetXaxis()->FindBin(R2off));
                                        double dark_sr = (double) Hist_Dark_SR_theta2.at(e).at(s).GetBinContent(Hist_Dark_SR_theta2.at(e).at(s).GetXaxis()->FindBin(R2off));
                                        if (dark_cr==0.) weight = 0.;
                                        else weight = dark_sr/dark_cr;
                                        Hist_Target_CR_theta2.at(e).at(s).Fill(theta2,weight);
                                        Hist_Target_CR_RaDec.at(e).at(s).Fill(ra_sky,dec_sky,weight);

                                        double dark_lza_cr = (double) Hist_DarkLZA_CR_theta2.at(e).at(s).GetBinContent(Hist_DarkLZA_CR_theta2.at(e).at(s).GetXaxis()->FindBin(R2off));
                                        double dark_lza_sr = (double) Hist_DarkLZA_SR_theta2.at(e).at(s).GetBinContent(Hist_DarkLZA_SR_theta2.at(e).at(s).GetXaxis()->FindBin(R2off));
                                        if (dark_lza_cr==0.) weight = 0.;
                                        else weight = dark_lza_sr/dark_lza_cr;
                                        Hist_TargetLZA_CR_theta2.at(e).at(s).Fill(theta2,weight);

                                        double dark_sza_cr = (double) Hist_DarkSZA_CR_theta2.at(e).at(s).GetBinContent(Hist_DarkSZA_CR_theta2.at(e).at(s).GetXaxis()->FindBin(R2off));
                                        double dark_sza_sr = (double) Hist_DarkSZA_SR_theta2.at(e).at(s).GetBinContent(Hist_DarkSZA_SR_theta2.at(e).at(s).GetXaxis()->FindBin(R2off));
                                        if (dark_sza_cr==0.) weight = 0.;
                                        else weight = dark_sza_sr/dark_sza_cr;
                                        Hist_TargetSZA_CR_theta2.at(e).at(s).Fill(theta2,weight);
                                    }
                                }
                                for (int s=0;s<Number_of_SR;s++)
                                {
                                    if (SignalSelectionMSCW(s)) Hist_Target_SR_MSCW.at(e).at(s).Fill(MSCW);
                                }
                                for (int s=0;s<Number_of_CR;s++)
                                {
                                    if (ControlSelectionMSCW(s)) Hist_Target_CR_MSCW.at(e).at(s).Fill(MSCW);
                                }
                                Hist_Target_SR_MSCL.at(e).Fill(MSCL);
                                Hist_Target_MSCLW.at(e).Fill(MSCL,MSCW);
                        }
                    }
                    input_file->Close();
                }
                std::cout << "=================================================================" << std::endl;
                std::cout << "Target, e " << energy_bins[e] << ", N_iter.at(e) = " << N_iter.at(e) << std::endl;
                std::cout << "Hist_Target_SR_MSCW.at(e).Integral() = " << Hist_Target_SR_MSCW.at(e).at(0).Integral() << std::endl;
                std::cout << "Hist_Target_SR_MSCW.at(e).GetMean() = " << Hist_Target_SR_MSCW.at(e).at(0).GetMean() << std::endl;
                std::cout << "Hist_Target_SR_MSCW.at(e).GetRMS() = " << Hist_Target_SR_MSCW.at(e).at(0).GetRMS() << std::endl;

                std::pair <double,double> offset;
                std::pair <double,double> kernel_rms;
                for (int c1=0;c1<Number_of_CR-1;c1++)
                {
                    Deconvolution(&Hist_Target_CR_MSCW.at(e).at(c1),&Hist_Target_CR_MSCW.at(e).at(Number_of_CR-1),&Hist_Target_Deconv_MSCW.at(e),N_iter.at(e));
                    double rms_begin = Hist_Target_Deconv_MSCW.at(e).GetRMS();
                    double mean_begin = Hist_Target_Deconv_MSCW.at(e).GetMean();

                    std::cout << "Target, e " << energy_bins[e] << ", running CR " << c1 << std::endl;
                    std::cout << "Target, e " << energy_bins[e] << ", mean_begin = " << mean_begin << std::endl;
                    std::cout << "Target, e " << energy_bins[e] << ", rms_begin = " << rms_begin << std::endl;

                    N_iter.at(e) = FindNIteration(&Hist_Target_CR_MSCW.at(e).at(Number_of_CR-1),&Hist_Target_CR_MSCW.at(e).at(c1),&Hist_Target_BkgCR_MSCW.at(e).at(Number_of_CR-1),&Hist_Target_BkgTemp_MSCW.at(e),&Hist_Target_Deconv_MSCW.at(e),0.,rms_begin,mean_begin,10,true);
                    std::cout << "Target, e " << energy_bins[e] << ", N_iter.at(e) = " << N_iter.at(e) << std::endl;
                    kernel_rms = FindRMS(&Hist_Target_CR_MSCW.at(e).at(Number_of_CR-1),&Hist_Target_CR_MSCW.at(e).at(c1),&Hist_Target_BkgCR_MSCW.at(e).at(Number_of_CR-1),&Hist_Target_BkgTemp_MSCW.at(e),&Hist_Target_Deconv_MSCW.at(e),0.,rms_begin,mean_begin,N_iter.at(e),true);
                    std::cout << "Target, e " << energy_bins[e] << ", CR_Kernel_RMS = " << kernel_rms.first << std::endl;
                    std::cout << "Target, e " << energy_bins[e] << ", unblinded CR_Kernel_RMS = " << kernel_rms.second << std::endl;

                    //N_iter.at(e) = FindNIteration(&Hist_Target_CR_MSCW.at(e).at(Number_of_CR-1),&Hist_Target_CR_MSCW.at(e).at(c1),&Hist_Target_BkgCR_MSCW.at(e).at(Number_of_CR-1),&Hist_Target_BkgTemp_MSCW.at(e),&Hist_Target_Deconv_MSCW.at(e),kernel_rms.first,mean_begin,N_iter.at(e),true);
                    std::cout << "Target, e " << energy_bins[e] << ", N_iter.at(e) = " << N_iter.at(e) << std::endl;
                    kernel_rms = FindRMS(&Hist_Target_CR_MSCW.at(e).at(Number_of_CR-1),&Hist_Target_CR_MSCW.at(e).at(c1),&Hist_Target_BkgCR_MSCW.at(e).at(Number_of_CR-1),&Hist_Target_BkgTemp_MSCW.at(e),&Hist_Target_Deconv_MSCW.at(e),0.,kernel_rms.first,mean_begin,N_iter.at(e),true);
                    std::cout << "Target, e " << energy_bins[e] << ", CR_blinded_RMS = " << GetBlindedRMS(&Hist_Target_CR_MSCW.at(e).at(c1)) << std::endl;
                    std::cout << "Target, e " << energy_bins[e] << ", SR_blinded_RMS = " << GetBlindedRMS(&Hist_Target_CR_MSCW.at(e).at(Number_of_CR-1)) << std::endl;
                    std::cout << "Target, e " << energy_bins[e] << ", CR_Kernel_RMS = " << kernel_rms.first << std::endl;
                    std::cout << "Target, e " << energy_bins[e] << ", unblinded CR_Kernel_RMS = " << kernel_rms.second << std::endl;

                    double initial_RMS = kernel_rms.first;



                    // estimate unblinded CR bkg
                    std::cout << "Target, e " << energy_bins[e] << ", running CR " << c1 << " to CR " << Number_of_CR-1 << std::endl;
                    myfunc->SetParameter(0,initial_RMS);
                    Hist_Target_Deconv_MSCW.at(e).Reset();
                    //Hist_Target_Deconv_MSCW.at(e).FillRandom("myfunc",Hist_Target_CR_MSCW.at(e).at(Number_of_CR-1).Integral()*100);
                    //Hist_Target_Deconv_MSCW.at(e).FillRandom("myfunc",Hist_Target_CR_MSCW.at(e).at(Number_of_CR-1).Integral()*10);
                    Hist_Target_Deconv_MSCW.at(e).FillRandom("myfunc",1000000);
                    Deconvolution(&Hist_Target_CR_MSCW.at(e).at(c1),&Hist_Target_Deconv_MSCW.at(e),&Hist_Target_BkgTemp_MSCW.at(e),N_iter.at(e));
                    offset = ShiftAndNormalize(&Hist_Target_CR_MSCW.at(e).at(c1),&Hist_Target_CR_MSCW.at(e).at(Number_of_CR-1),&Hist_Target_BkgTemp_MSCW.at(e),&Hist_Target_BkgCR_MSCW.at(e).at(Number_of_CR-1),0.,true,true);
                    double scaled_shift = offset.second;
                    std::cout << "CR blinded offset = " << offset.first << std::endl;
                    std::cout << "CR unblinded offset = " << offset.second << std::endl;
                    std::cout << "CR Mean = " << Hist_Target_BkgCR_MSCW.at(e).at(Number_of_CR-1).GetMean() << std::endl;
                    std::cout << "CR RMS = " << Hist_Target_BkgCR_MSCW.at(e).at(Number_of_CR-1).GetRMS() << std::endl;
                    std::pair <double,double> converge = FindConverge(&Hist_Target_CR_MSCW.at(e).at(Number_of_CR-1),&Hist_Target_BkgCR_MSCW.at(e).at(Number_of_CR-1),&Hist_Target_BkgTemp_MSCW.at(e));
                    if (DoConverge) Converge(&Hist_Target_BkgCR_MSCW.at(e).at(Number_of_CR-1),converge.first,converge.second);

                    Hist_Target_CR_MSCW_SumRuns.at(e).at(c1).Add(&Hist_Target_CR_MSCW.at(e).at(Number_of_CR-1));
                    Hist_Target_BkgCR_MSCW_SumRuns.at(e).at(c1).Add(&Hist_Target_BkgCR_MSCW.at(e).at(Number_of_CR-1));

                    // estimate SR1 bkg
                    MakeBkgPrevious(&Hist_Target_CR_MSCW.at(e).at(Number_of_CR-1),&Hist_Target_BkgCR_MSCW.at(e).at(Number_of_CR-1),&Hist_Target_BkgPrevious_MSCW.at(e));
                    std::cout << "Target, e " << energy_bins[e] << ", running CR " << c1 << " to SR " << 0 << std::endl;
                    int SR1_Niter = N_iter.at(e);
                    SR1_Niter = FindNIteration(&Hist_Target_SR_MSCW.at(e).at(0),&Hist_Target_BkgPrevious_MSCW.at(e),&Hist_Target_BkgSR_MSCW.at(e).at(0),&Hist_Target_BkgTemp_MSCW.at(e),&Hist_Target_Deconv_MSCW.at(e),scaled_shift,initial_RMS,mean_begin,N_iter.at(e),false);
                    std::cout << "Target, e " << energy_bins[e] << ", SR2_Niter = " << SR1_Niter << std::endl;
                    initial_kernel_rms.push_back(initial_RMS);
                    kernel_rms = FindRMS(&Hist_Target_SR_MSCW.at(e).at(0),&Hist_Target_BkgPrevious_MSCW.at(e),&Hist_Target_BkgSR_MSCW.at(e).at(0),&Hist_Target_BkgTemp_MSCW.at(e),&Hist_Target_Deconv_MSCW.at(e),scaled_shift,initial_RMS,mean_begin,SR1_Niter,false);
                    std::cout << "Target, e " << energy_bins[e] << ", CR_blinded_RMS = " << GetBlindedRMS(&Hist_Target_CR_MSCW.at(e).at(Number_of_CR-1)) << std::endl;
                    std::cout << "Target, e " << energy_bins[e] << ", SR_blinded_RMS = " << GetBlindedRMS(&Hist_Target_SR_MSCW.at(e).at(0)) << std::endl;
                    std::cout << "Target, e " << energy_bins[e] << ", SR0_Kernel_RMS = " << kernel_rms.first << std::endl;
                    std::cout << "Target, e " << energy_bins[e] << ", unblinded SR0_Kernel_RMS = " << kernel_rms.second << std::endl;
                    N_iter.at(e) = SR1_Niter;
                    final_kernel_rms.push_back(kernel_rms.first);
                    unblinded_kernel_rms.push_back(kernel_rms.second);
                    blinded_cr_rms.push_back(GetBlindedRMS(&Hist_Target_CR_MSCW.at(e).at(Number_of_CR-1)));
                    blinded_sr_rms.push_back(GetBlindedRMS(&Hist_Target_SR_MSCW.at(e).at(0)));

                    myfunc->SetParameter(0,kernel_rms.first);
                    Hist_Target_Deconv_MSCW.at(e).Reset();
                    //Hist_Target_Deconv_MSCW.at(e).FillRandom("myfunc",Hist_Target_SR_MSCW.at(e).at(0).Integral()*100);
                    //Hist_Target_Deconv_MSCW.at(e).FillRandom("myfunc",Hist_Target_SR_MSCW.at(e).at(0).Integral()*10);
                    Hist_Target_Deconv_MSCW.at(e).FillRandom("myfunc",1000000);
                    
                    Deconvolution(&Hist_Target_BkgPrevious_MSCW.at(e),&Hist_Target_Deconv_MSCW.at(e),&Hist_Target_BkgTemp_MSCW.at(e),SR1_Niter);
                    //Deconvolution(&Hist_Target_CR_MSCW.at(e).at(Number_of_CR-1),&Hist_Target_Deconv_MSCW.at(e),&Hist_Target_BkgTemp_MSCW.at(e),SR1_Niter);
                    std::cout << "Target, e " << energy_bins[e] << ", final rms = " << Hist_Target_Deconv_MSCW.at(e).GetRMS() << std::endl;
                    offset = ShiftAndNormalize(&Hist_Target_BkgPrevious_MSCW.at(e),&Hist_Target_SR_MSCW.at(e).at(0),&Hist_Target_BkgTemp_MSCW.at(e),&Hist_Target_BkgSR_MSCW.at(e).at(0),scaled_shift,true,false);
                    std::cout << "SR0 blinded offset = " << offset.first << std::endl;
                    std::cout << "SR0 unblinded offset = " << offset.second << std::endl;
                    std::cout << "SR0 Mean = " << Hist_Target_BkgSR_MSCW.at(e).at(0).GetMean() << std::endl;
                    std::cout << "SR0 RMS = " << Hist_Target_BkgSR_MSCW.at(e).at(0).GetRMS() << std::endl;
                    initial_shift.push_back(scaled_shift);
                    blinded_shift.push_back(offset.first);
                    unblinded_shift.push_back(offset.second);
                    if (DoConverge) Converge(&Hist_Target_BkgSR_MSCW.at(e).at(0),converge.first,converge.second);
                    //AddBkgStatistics(&Hist_Target_BkgSR_MSCW.at(e).at(0));
                    //AddSystematics(&Hist_Target_SR_MSCW.at(e).at(0),&Hist_Target_BkgSR_MSCW.at(e).at(0));
                    //AddSystematics2(&Hist_Target_CR_MSCW.at(e).at(Number_of_CR-1),&Hist_Target_BkgCR_MSCW.at(e).at(Number_of_CR-1),&Hist_Target_BkgSR_MSCW.at(e).at(0));
                    //AddSystematics3(&Hist_Target_BkgSR_MSCW.at(e).at(0));

                    double weight = GetLinearUncWithRange(&Hist_Target_BkgSR_MSCW.at(e).at(0),MSCW_cut_lower,Hist_Target_BkgSR_MSCW.at(e).at(0).GetMean());
                    std::cout << "weight = " << weight << std::endl; 
                    if (weight!=0) 
                    {
                        weight = 1./weight;
                        Weight_Target_BkgSR_MSCW_AllCR.at(0) += weight;
                        Hist_Target_BkgSR_MSCW_AllCR.at(e).at(0).Add(&Hist_Target_BkgSR_MSCW.at(e).at(0),weight);
                    }

                    // estimate SR bkg
                    for (int s=1;s<Number_of_SR;s++)
                    {
                        MakeBkgPrevious(&Hist_Target_SR_MSCW.at(e).at(s-1),&Hist_Target_BkgSR_MSCW.at(e).at(s-1),&Hist_Target_BkgPrevious_MSCW.at(e));
                        std::cout << "Target, e " << energy_bins[e] << ", running CR " << c1 << " to SR " << s << std::endl;
                        int SR_Niter = N_iter.at(e);
                        SR_Niter = FindNIteration(&Hist_Target_SR_MSCW.at(e).at(s),&Hist_Target_BkgPrevious_MSCW.at(e),&Hist_Target_BkgSR_MSCW.at(e).at(s),&Hist_Target_BkgTemp_MSCW.at(e),&Hist_Target_Deconv_MSCW.at(e),scaled_shift,initial_RMS,mean_begin,N_iter.at(e),false);
                        std::cout << "Target, e " << energy_bins[e] << ", SR2_Niter = " << SR_Niter << std::endl;
                        initial_kernel_rms.push_back(initial_RMS);
                        kernel_rms = FindRMS(&Hist_Target_SR_MSCW.at(e).at(s),&Hist_Target_BkgPrevious_MSCW.at(e),&Hist_Target_BkgSR_MSCW.at(e).at(s),&Hist_Target_BkgTemp_MSCW.at(e),&Hist_Target_Deconv_MSCW.at(e),scaled_shift,initial_RMS,mean_begin,SR_Niter,false);
                        std::cout << "Target, e " << energy_bins[e] << ", CR_blinded_RMS = " << GetBlindedRMS(&Hist_Target_SR_MSCW.at(e).at(s-1)) << std::endl;
                        std::cout << "Target, e " << energy_bins[e] << ", SR_blinded_RMS = " << GetBlindedRMS(&Hist_Target_SR_MSCW.at(e).at(s)) << std::endl;
                        std::cout << "Target, e " << energy_bins[e] << ", SR_Kernel_RMS = " << kernel_rms.first << std::endl;
                        std::cout << "Target, e " << energy_bins[e] << ", unblinded SR_Kernel_RMS = " << kernel_rms.second << std::endl;
                        N_iter.at(e) = SR1_Niter;
                        final_kernel_rms.push_back(kernel_rms.first);
                        unblinded_kernel_rms.push_back(kernel_rms.second);
                        blinded_cr_rms.push_back(GetBlindedRMS(&Hist_Target_SR_MSCW.at(e).at(s-1)));
                        blinded_sr_rms.push_back(GetBlindedRMS(&Hist_Target_SR_MSCW.at(e).at(s)));

                        myfunc->SetParameter(0,kernel_rms.first);
                        Hist_Target_Deconv_MSCW.at(e).Reset();
                        //Hist_Target_Deconv_MSCW.at(e).FillRandom("myfunc",Hist_Target_SR_MSCW.at(e).at(0).Integral()*100);
                        //Hist_Target_Deconv_MSCW.at(e).FillRandom("myfunc",Hist_Target_SR_MSCW.at(e).at(0).Integral()*10);
                        Hist_Target_Deconv_MSCW.at(e).FillRandom("myfunc",1000000);
                        
                        Deconvolution(&Hist_Target_BkgPrevious_MSCW.at(e),&Hist_Target_Deconv_MSCW.at(e),&Hist_Target_BkgTemp_MSCW.at(e),SR_Niter);
                        //Deconvolution(&Hist_Target_BkgSR_MSCW.at(e).at(s-1),&Hist_Target_Deconv_MSCW.at(e),&Hist_Target_BkgTemp_MSCW.at(e),SR_Niter);
                        std::cout << "Target, e " << energy_bins[e] << ", final rms = " << Hist_Target_Deconv_MSCW.at(e).GetRMS() << std::endl;
                        offset = ShiftAndNormalize(&Hist_Target_BkgPrevious_MSCW.at(e),&Hist_Target_SR_MSCW.at(e).at(s),&Hist_Target_BkgTemp_MSCW.at(e),&Hist_Target_BkgSR_MSCW.at(e).at(s),scaled_shift,true,false);
                        std::cout << "SR" << s << " blinded offset = " << offset.first << std::endl;
                        std::cout << "SR" << s << " unblinded offset = " << offset.second << std::endl;
                        std::cout << "SR" << s << " Mean = " << Hist_Target_BkgSR_MSCW.at(e).at(s).GetMean() << std::endl;
                        std::cout << "SR" << s << " RMS = " << Hist_Target_BkgSR_MSCW.at(e).at(s).GetRMS() << std::endl;
                        initial_shift.push_back(scaled_shift);
                        blinded_shift.push_back(offset.first);
                        unblinded_shift.push_back(offset.second);
                        if (DoConverge) Converge(&Hist_Target_BkgSR_MSCW.at(e).at(s),converge.first,converge.second);
                        //AddBkgStatistics(&Hist_Target_BkgSR_MSCW.at(e).at(s));
                        //AddSystematics(&Hist_Target_SR_MSCW.at(e).at(s),&Hist_Target_BkgSR_MSCW.at(e).at(s));
                        //AddSystematics2(&Hist_Target_CR_MSCW.at(e).at(Number_of_CR-1),&Hist_Target_BkgCR_MSCW.at(e).at(Number_of_CR-1),&Hist_Target_BkgSR_MSCW.at(e).at(s));
                        //AddSystematics3(&Hist_Target_BkgSR_MSCW.at(e).at(s));

                        double weight = GetLinearUncWithRange(&Hist_Target_BkgSR_MSCW.at(e).at(s),MSCW_cut_lower,Hist_Target_BkgSR_MSCW.at(e).at(s).GetMean());
                        std::cout << "weight = " << weight << std::endl; 
                        if (weight!=0) 
                        {
                            weight = 1./weight;
                            Weight_Target_BkgSR_MSCW_AllCR.at(s) += weight;
                            Hist_Target_BkgSR_MSCW_AllCR.at(e).at(s).Add(&Hist_Target_BkgSR_MSCW.at(e).at(s),weight);
                        }

                    }
                }
                for (int s=0;s<Number_of_SR;s++)
                {
                    if (Weight_Target_BkgSR_MSCW_AllCR.at(s)!=0) 
                    {
                        Hist_Target_BkgSR_MSCW_AllCR.at(e).at(s).Scale(1./Weight_Target_BkgSR_MSCW_AllCR.at(s));
                    }
                    Weight_Target_BkgSR_MSCW_AllCR.at(s) = 0.;
                }
                for (int s=0;s<Number_of_SR;s++)
                {
                    Hist_Target_SR_MSCW_SumRuns.at(e).at(s).Add(&Hist_Target_SR_MSCW.at(e).at(s));
                    Hist_Target_SR_MSCW_SumRuns_SumSRs.at(e).Add(&Hist_Target_SR_MSCW.at(e).at(s));
                    Hist_Target_BkgSR_MSCW_SumRuns.at(e).at(s).Add(&Hist_Target_BkgSR_MSCW_AllCR.at(e).at(s)); // this is summing over sublist, errors are not correlated.
                }
                for (int run=0;run<Sublist.at(subrun).size();run++)
                {
                    std::cout << "run " << int(Sublist.at(subrun)[run]) << std::endl;
                }
                for (int s=0;s<Number_of_SR;s++)
                {
                    int norm_bin_low = Hist_Target_BkgSR_MSCW_SumRuns.at(e).at(s).FindBin(MSCW_cut_lower);
                    int norm_bin_up = Hist_Target_BkgSR_MSCW_SumRuns.at(e).at(s).FindBin(MSCW_cut_blind);
                    std::cout << "Hist_Target_BkgSR_MSCW_SumRuns.at(e).at("<< s<< ").Integral() = " << Hist_Target_BkgSR_MSCW_SumRuns.at(e).at(s).Integral(norm_bin_low,norm_bin_up) << std::endl;
                    std::cout << "Hist_Target_SR_MSCW_SumRuns.at(e).at("<< s<< ").Integral() = " << Hist_Target_SR_MSCW_SumRuns.at(e).at(s).Integral(norm_bin_low,norm_bin_up) << std::endl;
                }
            }
        }
        for (int e=0;e<N_energy_bins;e++)
        {
            for (int s=0;s<Number_of_SR;s++)
            {
                Hist_Target_BkgSR_MSCW_SumRuns_SumSRs.at(e).Add(&Hist_Target_BkgSR_MSCW_SumRuns.at(e).at(s));
            }
        }



        int norm_bin_low_ring = Hist_Target_ON_MSCW_Alpha.FindBin(Norm_Lower);
        int norm_bin_up_ring = Hist_Target_ON_MSCW_Alpha.FindBin(Norm_Upper);
        double scale_ring = Hist_Target_ON_MSCW_Alpha.Integral(norm_bin_low_ring,norm_bin_up_ring)/Hist_Target_OFF_MSCW_Alpha.Integral(norm_bin_low_ring,norm_bin_up_ring);

        int Number_of_CR_new = Number_of_CR;
        int Number_of_SR_new = Number_of_SR;
        TString ConvergeOrNot = "";
        if (!DoConverge) ConvergeOrNot = "_NoConverge";
        TFile OutputFile("output_Apr23/Deconvolution_"+TString(target)+"_Ntel"+std::to_string(NTelMin)+"to"+std::to_string(NTelMax)+"_Elev"+std::to_string(int(Target_Elev_cut_lower))+"to"+std::to_string(int(Target_Elev_cut_upper))+"_Azim"+std::to_string(int(Target_Azim_cut_lower))+"to"+std::to_string(int(Target_Azim_cut_upper))+"_Theta2"+std::to_string(int(10.*Theta2_cut_lower))+"to"+std::to_string(int(10.*Theta2_cut_upper))+"_MSCWCut"+std::to_string(int(10.*MSCW_cut_upper))+"_MSCWBlind"+std::to_string(int(10.*MSCW_cut_blind))+ConvergeOrNot+".root","recreate");
        TTree InfoTree("InfoTree","info tree");
        InfoTree.Branch("Number_of_CR",&Number_of_CR_new,"Number_of_CR/I");
        InfoTree.Branch("Number_of_SR",&Number_of_SR_new,"Number_of_SR/I");
        InfoTree.Branch("MSCW_cut_lower",&MSCW_cut_lower,"MSCW_cut_lower/D");
        InfoTree.Branch("MSCW_cut_upper",&MSCW_cut_upper,"MSCW_cut_upper/D");
        InfoTree.Branch("MSCW_cut_blind",&MSCW_cut_blind,"MSCW_cut_blind/D");
        InfoTree.Branch("MSCL_cut_lower",&MSCL_signal_cut_lower[Number_of_SR-1],"MSCL_cut_lower/D");
        InfoTree.Branch("MSCL_cut_upper",&MSCL_signal_cut_upper[0],"MSCL_cut_upper/D");
        InfoTree.Branch("Theta2_upper_limit",&Theta2_upper_limit,"Theta2_upper_limit/D");
        InfoTree.Branch("Elev_cut_lower",&Elev_cut_lower,"Elev_cut_lower/D");
        InfoTree.Branch("Elev_cut_upper",&Elev_cut_upper,"Elev_cut_upper/D");
        InfoTree.Branch("Target_Elev_cut_lower",&Target_Elev_cut_lower,"Target_Elev_cut_lower/D");
        InfoTree.Branch("Target_Elev_cut_upper",&Target_Elev_cut_upper,"Target_Elev_cut_upper/D");
        InfoTree.Branch("Azim_cut_lower",&Azim_cut_lower,"Azim_cut_lower/D");
        InfoTree.Branch("Azim_cut_upper",&Azim_cut_upper,"Azim_cut_upper/D");
        InfoTree.Branch("Target_Azim_cut_lower",&Target_Azim_cut_lower,"Target_Azim_cut_lower/D");
        InfoTree.Branch("Target_Azim_cut_upper",&Target_Azim_cut_upper,"Target_Azim_cut_upper/D");
        InfoTree.Branch("used_runs","std::vector<int>",&used_runs);
        InfoTree.Branch("energy_vec","std::vector<double>",&energy_vec);
        InfoTree.Branch("exposure_hours",&exposure_hours,"exposure_hours/D");
        InfoTree.Branch("unblinded_shift","std::vector<double>",&unblinded_shift);
        InfoTree.Branch("blinded_shift","std::vector<double>",&blinded_shift);
        InfoTree.Branch("initial_shift","std::vector<double>",&initial_shift);
        InfoTree.Branch("initial_kernel_rms","std::vector<double>",&initial_kernel_rms);
        InfoTree.Branch("final_kernel_rms","std::vector<double>",&final_kernel_rms);
        InfoTree.Branch("unblinded_kernel_rms","std::vector<double>",&unblinded_kernel_rms);
        InfoTree.Branch("blinded_cr_rms","std::vector<double>",&blinded_cr_rms);
        InfoTree.Branch("blinded_sr_rms","std::vector<double>",&blinded_sr_rms);
        InfoTree.Fill();
        InfoTree.Write();
        Hist_Target_TelElevAzim.Write();
        Hist_Target_TelRaDec.Write();
        Hist_Target_TelRaDec_AfterCut.Write();
        for (int e=0;e<N_energy_bins;e++) {
                Hist_Target_SR_MSCL.at(e).Write();
                Hist_Target_Deconv_MSCW.at(e).Write();
                Hist_Target_MSCLW.at(e).Write();
                Hist_Target_SR_MSCW_SumRuns_SumSRs.at(e).Write();
                Hist_Target_BkgSR_MSCW_SumRuns_SumSRs.at(e).Write();
                for (int s=0;s<Number_of_SR;s++)
                {
                    Hist_Target_SR_MSCW_SumRuns.at(e).at(s).Write();
                    Hist_Target_BkgSR_MSCW_SumRuns.at(e).at(s).Write();
                    Hist_Target_SR_theta2.at(e).at(s).Write();
                    Hist_Dark_SR_theta2.at(e).at(s).Write();
                    Hist_DarkLZA_SR_theta2.at(e).at(s).Write();
                    Hist_DarkSZA_SR_theta2.at(e).at(s).Write();
                    Hist_Dark_CR_theta2.at(e).at(s).Write();
                    Hist_DarkLZA_CR_theta2.at(e).at(s).Write();
                    Hist_DarkSZA_CR_theta2.at(e).at(s).Write();
                    Hist_Target_CR_theta2.at(e).at(s).Write();
                    Hist_TargetLZA_CR_theta2.at(e).at(s).Write();
                    Hist_TargetSZA_CR_theta2.at(e).at(s).Write();
                    Hist_Target_SR_RaDec.at(e).at(s).Write();
                    Hist_Target_CR_RaDec.at(e).at(s).Write();
                }
                for (int s=0;s<Number_of_CR;s++)
                {
                    Hist_Target_CR_MSCW_SumRuns.at(e).at(s).Write();
                    Hist_Target_BkgCR_MSCW_SumRuns.at(e).at(s).Write();
                }
        }
        OutputFile.Close();

}
