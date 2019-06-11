
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
#include "TProfile.h"
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
#include "TSpline.h"

#include "GetRunList.h"

ClassImp(TSplinePoly);
ClassImp(TSplinePoly3);
ClassImp(TSplinePoly5);
ClassImp(TSpline3);
ClassImp(TSpline5);
ClassImp(TSpline);


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
double MSCW_cut_lower = -1.0;  // set this to -1.5 for V5
double MSCW_cut_blind = 1.5;
double MSCW_cut_upper = 1.5;
const int Number_of_SR = 6;
double MSCL_signal_cut_lower[Number_of_SR] = {0.75,0.50,0.25,0.00,-0.25,-0.50};
double MSCL_signal_cut_upper[Number_of_SR] = {1.00,0.75,0.50,0.25, 0.00,-0.25};
const int Number_of_CR = 6;
double MSCL_control_cut_lower[Number_of_CR] = {2.25,2.00,1.75,1.50,1.25,1.00};
double MSCL_control_cut_upper[Number_of_CR] = {2.50,2.25,2.00,1.75,1.50,1.25};
//const int Number_of_CR = 8;
//double MSCL_control_cut_lower[Number_of_CR] = {2.75,2.50,2.25,2.00,1.75,1.50,1.25,1.00};
//double MSCL_control_cut_upper[Number_of_CR] = {3.00,2.75,2.50,2.25,2.00,1.75,1.50,1.25};
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
double current_sr = 0.;
std::pair<double,std::pair<double,double>> estimated_parameters;
double estimated_amplitude = 0.;
double estimated_mean = 0.;
double estimated_mean_err = 0.;
double estimated_rms_previous = 0.;
double estimated_kernel_rms = 0.;
double estimated_rms = 0.;
double estimated_rms_err = 0.;
double dark_initial_mean_delta = 0.;
double dark_initial_mean = 0.;
double dark_current_mean = 0.;
double target_initial_mean_delta = 0.;
double target_initial_stat = 0.;
double target_initial_mean = 0.;
double target_current_mean = 0.;
double dark_initial_rms_delta = 0.;
double dark_initial_rms = 0.;
double dark_current_rms = 0.;
double target_initial_rms_delta = 0.;
double target_initial_rms = 0.;
double target_current_rms = 0.;
double dark_initial_amp_delta = 0.;
double dark_initial_amp = 0.;
double dark_current_amp = 0.;
double target_initial_amp_delta = 0.;
double target_initial_amp = 0.;
double target_current_amp = 0.;
double bkg_mean = 0.;
double bkg_rms = 0.;


std::pair <double,double> converge;
vector<int> used_runs;
vector<double> cosmic_electron;
double exposure_hours = 0.;
double delta_mscw_mean = 0.;
vector<std::pair<double,double>> dark_converge;
vector<vector<double>> dark_kernel_rms;
vector<vector<double>> dark_kernel_shift;
vector<vector<double>> target_kernel_rms;
vector<vector<double>> target_kernel_shift;

const int N_energy_bins = 18;
double energy_bins[N_energy_bins+1] = {200,237,282,335,398,473,562,667,794,943,1122,1332,1585,1882,2239,3162,4467,6310,8913};
bool use_this_energy_bin[N_energy_bins] = {false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false};
double electron_flux[N_energy_bins] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
double electron_flux_err[N_energy_bins] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
double electron_count[N_energy_bins] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
double electron_count_err[N_energy_bins] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

//const int N_energy_bins = 11;
//double energy_bins[N_energy_bins+1] = {200,282,398,562,794,1122,1585,2239,3162,4467,6310,8913};
//bool use_this_energy_bin[N_energy_bins] = {false,false,false,false,false,false,false,false,false,false,false};
//double electron_flux[N_energy_bins] = {0,0,0,0,0,0,0,0,0,0,0};
//double electron_flux_err[N_energy_bins] = {0,0,0,0,0,0,0,0,0,0,0};
//double electron_count[N_energy_bins] = {0,0,0,0,0,0,0,0,0,0,0};
//double electron_count_err[N_energy_bins] = {0,0,0,0,0,0,0,0,0,0,0};

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
TObject* getEffAreaHistogram( TFile* fAnasumDataFile, int runnumber)
{
        double iSlizeY = -9999;
        string dirname = "energyHistograms";
        string hisname = "herecEffectiveArea_on";
	if( !fAnasumDataFile )
	{
		return 0;
	}
	
	char dx[600];
	if( runnumber > 1 )
	{
		sprintf( dx, "run_%d/stereo/%s", runnumber, dirname.c_str() );
	}
	else
	{
		if( runnumber == 0 )
		{
			sprintf( dx, "total/stereo/%s", dirname.c_str() );
		}
		else if( runnumber == 1 )
		{
			sprintf( dx, "total_%d/stereo/%s", runnumber, dirname.c_str() );
		}
		else
		{
			sprintf( dx, "total_%d/stereo/%s", -1 * runnumber, dirname.c_str() );
		}
	}
	
	fAnasumDataFile->cd( dx );
	TDirectory* iDir = gDirectory;
	if( !iDir )
	{
		return 0;
	}
	
	TObject* h = ( TObject* )iDir->Get( hisname.c_str() );
	
	if( h && iSlizeY < -9998. )
	{
		return h->Clone();
	}
	else if( h )
	{
		string iClassName = h->ClassName();
		if( iClassName.find( "TH2" ) != string::npos )
		{
			TH2* i_h2 = ( TH2* )h;
			string iN = hisname + "px";
			TH1* i_h = ( TH1* )i_h2->ProjectionX( iN.c_str(), i_h2->GetYaxis()->FindBin( iSlizeY ), i_h2->GetYaxis()->FindBin( iSlizeY ) );
			return i_h->Clone();
		}
	}
	
	
	return 0;
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
std::pair<double,double> GetRangeMeanRMS(TH1* Hist_SR,double low,double up)
{
    TH1D Hist_SR_blinded = TH1D("Hist_SR_blinded","",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper);
    for (int i=0;i<Hist_SR->GetNbinsX();i++) {
        if (Hist_SR->GetBinCenter(i+1)<low || Hist_SR->GetBinCenter(i+1)>up) {
                Hist_SR_blinded.SetBinContent(i+1,0.);
                Hist_SR_blinded.SetBinError(i+1,0.);
        }
        else
        {
                Hist_SR_blinded.SetBinContent(i+1,Hist_SR->GetBinContent(i+1));
                Hist_SR_blinded.SetBinError(i+1,Hist_SR->GetBinError(i+1));
        }
    }
    return std::make_pair(Hist_SR_blinded.GetMean(),Hist_SR_blinded.GetRMS());
}
double GetBlindedRMS(TH1* Hist_SR)
{
    TH1D Hist_SR_blinded = TH1D("Hist_SR_blinded","",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper);
    for (int i=0;i<Hist_SR->GetNbinsX();i++) {
        if (Hist_SR->GetBinCenter(i+1)<MSCW_cut_blind) {
                Hist_SR_blinded.SetBinContent(i+1,0.);
                Hist_SR_blinded.SetBinError(i+1,0.);
        }
        else
        {
                Hist_SR_blinded.SetBinContent(i+1,Hist_SR->GetBinContent(i+1));
                Hist_SR_blinded.SetBinError(i+1,Hist_SR->GetBinError(i+1));
        }
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
double GetChi2(TH1* Hist_SR, TH1* Hist_Bkg, bool includeSR, int chi2_type) {
    if (Hist_Bkg->Integral()==0) return 0.;
    double chi2_temp = 0.;
    double chi2_mean = 0.;
    double chi2_rms = 0.;
    double sign_chi2_temp = 0.;
    double mean = Hist_SR->GetMean();
    double rms = Hist_SR->GetRMS();
    double nbins = 0.;
    double sign_nbins = 0.;
    double inflation_mean = 1.;
    double inflation_rms = 1.;
    int norm_bin_low = Hist_SR->FindBin(MSCW_cut_lower);
    int norm_bin_blind = Hist_SR->FindBin(MSCW_cut_blind);
    int norm_bin_up = Hist_SR->FindBin(30.);
    double weight_blinded = double(Hist_SR->Integral(norm_bin_low,norm_bin_blind));
    double weight_unblinded = double(Hist_SR->Integral(norm_bin_blind,norm_bin_up));
    double weight_ratio = weight_blinded/(weight_blinded+weight_unblinded);
    if (!includeSR) 
    {
        //inflation_rms = exp(-0.5*pow((Hist_Bkg->GetRMS()-estimated_rms)/(0.25*estimated_rms_err*weight_ratio),2));
        inflation_rms = exp(-0.5*pow((Hist_Bkg->GetRMS()-estimated_rms)/(0.5*(MSCW_cut_blind-MSCW_cut_lower)),2));
        //inflation_mean = exp(-0.5*pow((Hist_Bkg->GetMean()-estimated_mean)/(0.25*estimated_mean_err),2));
        chi2_mean = pow(Hist_Bkg->GetMean()-estimated_mean,2)/(estimated_mean_err*estimated_mean_err);
        chi2_rms = pow((Hist_Bkg->GetRMS()-estimated_rms)/(estimated_rms_err),2);
    }
    for (int i=0;i<Hist_SR->GetNbinsX();i++) {
        double bkg = Hist_Bkg->GetBinContent(i+1);
        double data = Hist_SR->GetBinContent(i+1);
        double bkg_err = Hist_Bkg->GetBinError(i+1);
        double data_err = Hist_SR->GetBinError(i+1);
        if (!includeSR && Hist_Bkg->GetBinCenter(i+1)>MSCW_cut_lower && Hist_Bkg->GetBinCenter(i+1)<MSCW_cut_blind) {
            continue;
        }
        if ((data_err*data_err)==0) data_err = 1.;
        chi2_temp += pow(bkg-data,2)/(data_err*data_err);
        nbins += 1.;
    }
    chi2_temp = chi2_temp/nbins;
    double chi2_total = 1.;
    if (includeSR || chi2_type==1) 
    {
        chi2_total = chi2_temp;
    }
    else
    {
        //chi2_total = chi2_rms;
        //chi2_total += chi2_temp;
        chi2_total = chi2_temp;
    }
    chi2_total = 1./chi2_total;
    //chi2_total = inflation_rms*1./chi2_total;
    return chi2_total;
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
double FitConstantFunction(Double_t *x, Double_t *par) {
    double xx =x[0];
    return par[0];
}
double FitLinearFunction(Double_t *x, Double_t *par) {
    double xx =x[0];
    return par[0]*xx+par[1];
}
double FitPowerLawFunction(Double_t *x, Double_t *par) {
    double xx =x[0];
    return par[0]*pow(xx,par[1]);
}
double FitFunction(Double_t *x, Double_t *par) {
    double xx =x[0];
    //if (xx-par[1]<0.) return 0.;
    //return par[0]*pow((xx-par[1])/(par[2]),2);
    if (xx-par[1]<0.) return 0.;
    return par[0]*pow((xx-par[1])/(par[2]),1)*exp(-0.5*pow((xx-par[3])/(par[4]),2));
}
void MakeSmoothSplineFunction(TH1* Hist_SR)
{
    double bin_width = Hist_SR->GetBinLowEdge(2)-Hist_SR->GetBinLowEdge(1);
    double nbins_in_RMS = Hist_SR->GetRMS()/bin_width;
    int knot_size = int(nbins_in_RMS/1.);
    if (knot_size<1) knot_size = 1;
    else if (knot_size<2) knot_size = 1;
    else if (knot_size<4) knot_size = 2;
    else if (knot_size<8) knot_size = 4;
    else if (knot_size<16) knot_size = 8;
    else knot_size = 16;
    //std::cout << "Smooth spline bin size = " << knot_size << std::endl;
    TH1D Hist_SR_rebin = TH1D("Hist_SR_rebin","",Hist_SR->GetNbinsX(),MSCW_plot_lower,MSCW_plot_upper);
    Hist_SR_rebin.Reset();
    Hist_SR_rebin.Add(Hist_SR);
    Hist_SR_rebin.Rebin(knot_size);
    TSpline3 spline = TSpline3(&Hist_SR_rebin);
    for (int i=0;i<Hist_SR->GetNbinsX();i++) {
        double x = Hist_SR->GetBinCenter(i+1);
        Hist_SR->SetBinContent(i+1,spline.Eval(x)/double(knot_size));
    }
}
void MakeBkgPrevious(TH1* Hist_SR,TH1* Hist_Bkg,TH1* Hist_Previous, bool includeSR)
{
    for (int i=0;i<Hist_SR->GetNbinsX();i++)
    {
        if (Hist_Bkg->GetBinCenter(i+1)<MSCW_cut_blind && !includeSR)
        //if (!includeSR)
        {
            Hist_Previous->SetBinContent(i+1,Hist_Bkg->GetBinContent(i+1));
            Hist_Previous->SetBinError(i+1,Hist_Bkg->GetBinError(i+1));
        }
        else
        {
            Hist_Previous->SetBinContent(i+1,Hist_SR->GetBinContent(i+1));
            Hist_Previous->SetBinError(i+1,Hist_SR->GetBinError(i+1));
        }
    }
    //MakeSmoothSplineFunction(Hist_Previous);
    if (Hist_Previous->Integral()==0) 
    {
        std::cout << "Hist_Previous->Integral() = 0!!!" << std::endl;
    }
}
double ConvergeFunction(double x, double endpoint_0, double endpoint_1, int type)
{
    if (x-(endpoint_1)>=0.) return 1.;
    if (x-(endpoint_0)<0.) return 0.;
    return pow((x-endpoint_0)/(endpoint_1-endpoint_0),1);
    //return 1./(1.+exp(-(x-(endpoint_0+endpoint_1)/2.)/(endpoint_1-endpoint_0)));

}
//void Converge(TH1* Hist_Bkg, double endpoint_0_width, double endpoint_1_width)
void Converge(TH1* Hist_Bkg, double endpoint_0, double endpoint_1)
{
    //double endpoint_0 = estimated_mean-endpoint_0_width*estimated_rms;
    //double endpoint_1 = estimated_mean-endpoint_1_width*estimated_rms;
    //std::cout << "estimated_mean = " << estimated_mean << std::endl;
    //std::cout << "estimated_rms = " << estimated_rms << std::endl;
    //std::cout << "endpoint_0_width = " << endpoint_0_width << std::endl;
    //std::cout << "endpoint_1_width = " << endpoint_1_width << std::endl;
    //std::cout << "endpoint_0 = " << endpoint_0 << std::endl;
    //std::cout << "endpoint_1 = " << endpoint_1 << std::endl;
    for (int i=0;i<Hist_Bkg->GetNbinsX();i++)
    {
        double old_content = Hist_Bkg->GetBinContent(i+1);
        double new_content = old_content*ConvergeFunction(Hist_Bkg->GetBinCenter(i+1),endpoint_0,endpoint_1,0);
        double old_error = Hist_Bkg->GetBinError(i+1);
        double new_error = old_error;
        Hist_Bkg->SetBinContent(i+1,new_content);
        Hist_Bkg->SetBinError(i+1,new_error);
    }
}
std::pair <double,double> FindConvergeEndpoints(TH1* Hist_SR,TH1* Hist_Bkg,double init_endpoint_0,double init_endpoint_1)
{
    double endpoint_0_best = 0.;
    double endpoint_1_best = 0.;
    double chi2_best = 0.;
    for (int ep0=0;ep0<50;ep0++) 
    {
        double tmp_endpoint_0 = -1.0-1.0+2.0*double(ep0)/50.;
        for (int ep1=0;ep1<50;ep1++) 
        {
            double tmp_endpoint_1 = 1.0-1.0+2.0*double(ep1)/50.;
            TH1D Hist_Bkg_Tmp = TH1D("Hist_Bkg_Tmp","",Hist_SR->GetNbinsX(),MSCW_plot_lower,MSCW_plot_upper);
            Hist_Bkg_Tmp.Reset();
            Hist_Bkg_Tmp.Add(Hist_Bkg);
            Converge(&Hist_Bkg_Tmp,tmp_endpoint_0,tmp_endpoint_1);
            double chi2 = GetChi2(Hist_SR,&Hist_Bkg_Tmp,true,1);
            if (chi2_best<chi2) {
                chi2_best = chi2;
                endpoint_0_best = tmp_endpoint_0;
                endpoint_1_best = tmp_endpoint_1;
            } 
        }
    }
    return std::make_pair(endpoint_0_best,endpoint_1_best);

}
void FindAverage(TH1* Hist_CR)
{
    TF1 *func = new TF1("func",FitConstantFunction,MSCL_signal_cut_lower[Number_of_SR-1],MSCL_control_cut_upper[0],1);
    func->SetParameter(0,-1.);
    Hist_CR->Fit("func","","",MSCL_control_cut_lower[Number_of_CR-1],MSCL_control_cut_upper[0]);
    for (int s=0;s<Number_of_SR;s++)
    {
        int bin = Hist_CR->FindBin(0.5*(MSCL_signal_cut_lower[s]+MSCL_signal_cut_upper[s]));
        Hist_CR->SetBinContent(bin,func->Eval(0.5*(MSCL_signal_cut_lower[s]+MSCL_signal_cut_upper[s])));
    }
}
void FindSRMean(TH1* Hist_CR, double energy)
{
    TF1 *func = new TF1("func",FitLinearFunction,MSCL_signal_cut_lower[Number_of_SR-1],MSCL_control_cut_upper[0],2);
    func->SetParameter(0,2.);
    func->SetParameter(1,2.);
    Hist_CR->Fit("func","","",MSCL_control_cut_lower[Number_of_CR-1],MSCL_control_cut_upper[0]);
    for (int s=0;s<Number_of_SR;s++)
    {
        int bin = Hist_CR->FindBin(0.5*(MSCL_signal_cut_lower[s]+MSCL_signal_cut_upper[s]));
        Hist_CR->SetBinContent(bin,func->Eval(0.5*(MSCL_signal_cut_lower[s]+MSCL_signal_cut_upper[s])));
        double width = 5.;
        if (energy>500.) width = 1.0;
        if (energy>700.) width = 0.5;
        if (energy>1000.) width = 0.2;
        if (energy>1500.) width = 0.1;
        //if (energy>1800.) width = 5.;
        Hist_CR->SetBinError(bin,width*(MSCW_cut_blind-MSCW_cut_lower));
    }
}
std::pair <bool,std::pair <double,double>> ShiftAndNormalize(TH1* Hist_SR, TH1* Hist_BkgTemp, TH1* Hist_Bkg, bool doShift, bool includeSR, int chi2_type) {
    includeSR = false;
    double shift_fit = 0;
    double unblinded_shift_fit = 0;
    double scale_fit = 0;
    double chi2 = 0.;
    double unblinded_chi2 = 0.;
    double mean = Hist_SR->GetMean();
    double rms = Hist_SR->GetRMS();
    double norm = Hist_BkgTemp->Integral();
    double shift_begin = 0.;
    if (!includeSR) 
    {
      shift_begin = -Hist_BkgTemp->GetMean()+estimated_mean; 
    }
    Hist_Bkg->Reset();
    int norm_bin_low = Hist_SR->FindBin(MSCW_cut_lower);
    int norm_bin_blind = Hist_SR->FindBin(MSCW_cut_blind);
    int norm_bin_up = Hist_SR->FindBin(30.);
    double weight_blinded = double(Hist_SR->Integral(norm_bin_low,norm_bin_blind));
    double weight_unblinded = double(Hist_SR->Integral(norm_bin_blind,norm_bin_up));
    double weight_ratio = weight_unblinded/(weight_blinded+weight_unblinded);
    if (doShift) {
        for (int fit=0;fit<100;fit++) {
                double shift = shift_begin-5.0*estimated_mean_err+10.0*estimated_mean_err*double(fit)*0.01;
                if (!includeSR) shift = shift_begin-1.0*estimated_mean_err*weight_ratio*double(fit)*0.01;
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
                double chi2_temp = GetChi2(Hist_SR, Hist_Bkg,includeSR,1);
                double unblinded_chi2_temp = GetChi2(Hist_SR, Hist_Bkg,true,1);
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
        for (int i=0;i<Hist_SR->GetNbinsX();i++) {
                int b = Hist_SR->FindBin(Hist_SR->GetBinCenter(i+1)-shift_begin);
                Hist_Bkg->SetBinContent(i+1,Hist_BkgTemp->GetBinContent(b));
                Hist_Bkg->SetBinError(i+1,Hist_BkgTemp->GetBinError(b));
        }
        Hist_Bkg->Scale(norm/Hist_Bkg->Integral());
        norm_bin_low = Hist_SR->FindBin(-2.);
        norm_bin_up = Hist_SR->FindBin(MSCW_cut_lower);
        double SR_area1 = Hist_SR->Integral(norm_bin_low,norm_bin_up);
        double Bkg_area1 = Hist_Bkg->Integral(norm_bin_low,norm_bin_up);
        norm_bin_low = Hist_SR->FindBin(MSCW_cut_blind);
        norm_bin_up = Hist_SR->FindBin(30.);
        double SR_area2 = Hist_SR->Integral(norm_bin_low,norm_bin_up);
        double Bkg_area2 = Hist_Bkg->Integral(norm_bin_low,norm_bin_up);
        double scale_begin = (SR_area1+SR_area2)/(Bkg_area1+Bkg_area2);
        scale_fit = scale_begin;
        for (int fit=0;fit<1000;fit++)
        {
            double scale = scale_begin-(double(fit))/1000.*0.1*scale_begin;
            Hist_Bkg->Scale(norm/Hist_Bkg->Integral());
            Hist_Bkg->Scale(scale);
            double chi2_temp = GetChi2(Hist_SR, Hist_Bkg,includeSR,1);
            if (chi2<chi2_temp && Hist_Bkg->Integral()!=0) {
                chi2 = chi2_temp;
                scale_fit = scale;
            } 
        }
        shift_fit = shift_begin;
        unblinded_shift_fit = shift_begin;
    }
    for (int i=0;i<Hist_SR->GetNbinsX();i++) {
            int b = Hist_Bkg->FindBin(Hist_Bkg->GetBinCenter(i+1)-shift_fit);
            Hist_Bkg->SetBinContent(i+1,Hist_BkgTemp->GetBinContent(b));
            Hist_Bkg->SetBinError(i+1,Hist_BkgTemp->GetBinError(b));
    }
    if (Hist_BkgTemp->Integral()==0) 
    {
        std::cout << "Hist_BkgTemp->Integral() = 0!!!" << std::endl;
    }
    if (Hist_Bkg->Integral()==0) 
    {
        std::cout << "Hist_Bkg->Integral() = 0!!!" << std::endl;
        for (int i=0;i<Hist_SR->GetNbinsX();i++) {
                Hist_Bkg->SetBinContent(i+1,0.);
                Hist_Bkg->SetBinError(i+1,0.);
        }
    }
    else 
    {
        Hist_Bkg->Scale(norm/Hist_Bkg->Integral());
    }
    Hist_Bkg->Scale(scale_fit);
    return std::make_pair(doShift,std::make_pair(shift_fit,unblinded_shift_fit));
}
bool IsReasonableResult(TH1* Hist_Bkg)
{
    double total_bkg = 0.;
    double total_err = 0.;
    for (int i=0;i<Hist_Bkg->GetNbinsX();i++)
    {
        total_bkg += Hist_Bkg->GetBinContent(i);
        total_err += Hist_Bkg->GetBinError(i);
    }
    if (total_bkg<=0.) return false;
    if (total_err/total_bkg>10.) return false;
    return true;
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
            Hist_Deconv->SetBinError(i+1,pow(max(source[i],0.),0.5));
        }
}

double FindNIteration(TH1* Hist_SR, TH1* Hist_CR, TH1* Hist_Bkg, TH1* Hist_BkgTemp, TH1* Hist_Deconv, double rms, double mean, int n_iter_begin, bool includeSR, int chi2_type) {
    return 10;
    //if (!includeSR) return n_iter_begin;
    //if (n_iter_begin==1) return 1;
    double n_iter_final = n_iter_begin;
    double chi2_best = 0.;
    chi2_best = 0.;
    TF1 *func = new TF1("func",Kernel,-50.,50.,1);
    func->SetParameter(0,0.5);
    int nbins = Hist_SR->GetNbinsX();
    for (int delta_n_iter = 0;delta_n_iter<10;delta_n_iter++) {
          int n_iter = n_iter_begin;
          //if (includeSR) n_iter = n_iter_begin+5-delta_n_iter;
          //else n_iter = n_iter_begin-delta_n_iter;
          n_iter = n_iter_begin+5-delta_n_iter;
          if (!includeSR) n_iter = n_iter_begin-delta_n_iter;
          if (n_iter<1) {
              continue;
          }
          double chi2 = 0;
          func->SetParameter(0,rms);
          Hist_Deconv->Reset();
          //Hist_Deconv->FillRandom("func",1000000);
          for (int b=0;b<Hist_Deconv->GetNbinsX();b++)
          {
              double content = func->Eval(Hist_Deconv->GetBinCenter(b+1));
              Hist_Deconv->SetBinContent(b+1,content);
              Hist_Deconv->SetBinError(b+1,pow(content,0.5));
          }
          Deconvolution(Hist_CR,Hist_Deconv,Hist_BkgTemp,n_iter);
          if (!IsReasonableResult(Hist_BkgTemp)) continue;
          std::pair <bool,std::pair <double,double>> offset = ShiftAndNormalize(Hist_SR,Hist_BkgTemp,Hist_Bkg,true,includeSR,chi2_type);
          chi2 = GetChi2(Hist_SR,Hist_Bkg,includeSR,1);
          if (chi2_best<chi2) {
              chi2_best = chi2;
              n_iter_final = n_iter;
          } 
    }
    return n_iter_final;
}
double FindRMS(TH1* Hist_SR, TH1* Hist_CR, TH1* Hist_Bkg, TH1* Hist_BkgTemp, TH1* Hist_Deconv, double rms_begin, double mean, double n_iter, bool includeSR, int chi2_type, double endpoint_0, double endpoint_1) {
    bool DoShift = true;
    double chi2_best = 0.;
    double rms_final = rms_begin;
    double unblinded_chi2_best = 0.;
    double unblinded_rms_final = rms_begin;
    double this_blinded_sr_rms = GetBlindedRMS(Hist_SR);
    double this_blinded_cr_rms = GetBlindedRMS(Hist_CR);
    int final_n_rms = 0;
    TF1 *func = new TF1("func",Kernel,-50.,50.,1);
    func->SetParameter(0,0.5);
    for (int n_rms = 0; n_rms<=40;n_rms++) {
        double chi2 = 0;
        double unblinded_chi2 = 0;
        double rms = rms_begin;
        rms = 0.1+double(n_rms+1)*2.0/40.;
        if (!includeSR) rms = rms_begin*(0.5+double(n_rms)*1.0/40.);
        func->SetParameter(0,rms);
        Hist_Deconv->Reset();
        for (int b=0;b<Hist_Deconv->GetNbinsX();b++)
        {
            double content = func->Eval(Hist_Deconv->GetBinCenter(b+1));
            Hist_Deconv->SetBinContent(b+1,content);
            Hist_Deconv->SetBinError(b+1,pow(content,0.5));
        }
        Deconvolution(Hist_CR,Hist_Deconv,Hist_BkgTemp,n_iter);
        if (!IsReasonableResult(Hist_BkgTemp)) continue;
        std::pair <bool,std::pair <double,double>> offset = ShiftAndNormalize(Hist_SR,Hist_BkgTemp,Hist_Bkg,DoShift,false,chi2_type);
        Converge(Hist_Bkg,endpoint_0,endpoint_1);
        chi2 = GetChi2(Hist_SR, Hist_Bkg,includeSR,2);
        unblinded_chi2 = GetChi2(Hist_SR, Hist_Bkg,true,2);
        if (chi2_best<chi2) {
            chi2_best = chi2;
            rms_final = rms;
            final_n_rms = n_rms;
        } 
        if (unblinded_chi2_best<unblinded_chi2) {
            unblinded_chi2_best = unblinded_chi2;
            unblinded_rms_final = rms;
        } 
    }
    //std::cout << "estimated_rms_previous = " << estimated_rms_previous << std::endl;
    //std::cout << "estimated_rms = " << estimated_rms << std::endl;
    //std::cout << "final_kernel_rms = " << rms_final << std::endl;
    //if (final_n_rms==0) std::cout << "kernel RMS solution attemps to move outside lower bound!!" << std::endl;
    //if (final_n_rms==40) std::cout << "kernel RMS solution attemps to move outside upper bound!!" << std::endl;
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

std::pair <double,std::pair <double,double>> PredictNextLayerHadron(TH1* Hist_ElectronMC, TH1* Hist_SR, TH1* Hist_SR_Previous,TH1* Hist_Bkg_Previous, TH1* Hist_Bkg, double energy, std::pair <double,std::pair <double,double>> parameters, bool useOldSR, bool isDark, bool doConverge)
{
    Hist_Bkg->Reset();
    TH1D Hist_Bkg_Privous_Adapt = TH1D("Hist_Bkg_Privous_Adapt","",Hist_SR->GetNbinsX(),MSCW_plot_lower,MSCW_plot_upper);
    if (!isDark)
    {
        MakeBkgPrevious(Hist_SR_Previous,Hist_Bkg_Previous,&Hist_Bkg_Privous_Adapt,useOldSR);
    }
    else
    {
        Hist_Bkg_Privous_Adapt.Reset();
        Hist_Bkg_Privous_Adapt.Add(Hist_SR_Previous);
    }

    TH1D Hist_SR_Temp = TH1D("Hist_SR_Temp","",Hist_SR->GetNbinsX(),MSCW_plot_lower,MSCW_plot_upper);
    Hist_SR_Temp.Reset();
    for (int b=0;b<Hist_SR_Temp.GetNbinsX();b++)
    {
        double content = 0;
        content = Hist_SR->GetBinContent(b+1)-Hist_ElectronMC->GetBinContent(b+1);
        double error = Hist_SR->GetBinError(b+1);
        Hist_SR_Temp.SetBinContent(b+1,content);
        Hist_SR_Temp.SetBinError(b+1,error);
    }
    std::cout << "Hist_SR_Temp.Integral() = " << Hist_SR_Temp.Integral() << std::endl;
    TH1D Hist_Bkg_Temp = TH1D("Hist_Bkg_Temp","",Hist_SR->GetNbinsX(),MSCW_plot_lower,MSCW_plot_upper);
    TH1D Hist_Kernel = TH1D("Hist_Kernel","",Hist_SR->GetNbinsX(),MSCW_plot_lower,MSCW_plot_upper);
    int SR_Niter=10;
    double kernel_rms = parameters.first;
    double endpoint_0 = parameters.second.first;
    double endpoint_1 = parameters.second.second;
    if (isDark) kernel_rms = FindRMS(&Hist_SR_Temp,&Hist_Bkg_Privous_Adapt,Hist_Bkg,&Hist_Bkg_Temp,&Hist_Kernel,kernel_rms,0,SR_Niter,true,1,endpoint_0,endpoint_1);
    else kernel_rms = FindRMS(&Hist_SR_Temp,&Hist_Bkg_Privous_Adapt,Hist_Bkg,&Hist_Bkg_Temp,&Hist_Kernel,kernel_rms,0,SR_Niter,false,1,endpoint_0,endpoint_1);
    TF1 *myfunc = new TF1("myfunc",Kernel,-50.,50.,1);
    myfunc->SetParameter(0,kernel_rms);
    Hist_Kernel.Reset();
    for (int b=0;b<Hist_Kernel.GetNbinsX();b++)
    {
        double content = myfunc->Eval(Hist_Kernel.GetBinCenter(b+1));
        Hist_Kernel.SetBinContent(b+1,content);
        Hist_Kernel.SetBinError(b+1,pow(content,0.5));
    }
                    
    Deconvolution(&Hist_Bkg_Privous_Adapt,&Hist_Kernel,&Hist_Bkg_Temp,SR_Niter);
    std::pair <bool,std::pair <double,double>> offset = ShiftAndNormalize(&Hist_SR_Temp,&Hist_Bkg_Temp,Hist_Bkg,true,false,2);
    MakeSmoothSplineFunction(Hist_Bkg);

    std::pair <double,double> converge;
    converge.first = endpoint_0;
    converge.second = endpoint_1;
    if (isDark)
    {
        converge = FindConvergeEndpoints(&Hist_SR_Temp,Hist_Bkg,endpoint_0,endpoint_1);
    }
    if (doConverge){ 
        std::cout << "Apply converge, endpoint 0 = " << converge.first << ", endpoint 1 = " << converge.second << std::endl;
        Converge(Hist_Bkg,converge.first,converge.second);
    }
    parameters.first = kernel_rms;
    parameters.second.first = converge.first;
    parameters.second.second = converge.second;
    return parameters;
}
void PredictNextLayer(TH1* Hist_GammaMC, TH1* Hist_ElectronMC, TH1* Hist_SR, TH1* Hist_SR_Previous,TH1* Hist_Bkg_Previous, TH1* Hist_Bkg, TH1* Hist_DarkSR, TH1* Hist_DarkSR_Previous, double energy, std::pair <double,std::pair <double,double>> parameters, bool useOldSR, bool doConverge)
{

    parameters = PredictNextLayerHadron(Hist_ElectronMC,Hist_DarkSR,Hist_DarkSR_Previous,Hist_Bkg_Previous,Hist_Bkg,energy,parameters,useOldSR,true,doConverge);
    std::cout << "Dark field found kernel RMS = " << parameters.first << std::endl;
    std::cout << "Dark field found endpoint 0 = " << parameters.second.first << std::endl;
    std::cout << "Dark field found endpoint 1 = " << parameters.second.second << std::endl;
    parameters = PredictNextLayerHadron(Hist_ElectronMC,Hist_SR,Hist_SR_Previous,Hist_Bkg_Previous,Hist_Bkg,energy,parameters,useOldSR,false,doConverge);
    TH1D Hist_MC_Temp2 = TH1D("Hist_MC_Temp2","",Hist_SR->GetNbinsX(),MSCW_plot_lower,MSCW_plot_upper);
    int norm_bin_low = Hist_SR->FindBin(MSCW_cut_lower);
    int norm_bin_up = Hist_SR->FindBin(MSCW_cut_blind);
    double old_integral = Hist_GammaMC->Integral(norm_bin_low,norm_bin_up);
    double new_integral = Hist_SR->Integral(norm_bin_low,norm_bin_up)-Hist_Bkg->Integral(norm_bin_low,norm_bin_up)-Hist_ElectronMC->Integral(norm_bin_low,norm_bin_up);
    if (new_integral<0)
    {
        Hist_GammaMC->Scale(0);
    }
    else
    {
        Hist_GammaMC->Scale(0);
        //double scale = new_integral/old_integral;
        //Hist_GammaMC->Scale(scale);
    }

}

double PointingDistribution(string file_name,int run, bool isDark, bool fillHist)
{
    if (run>100000) return 0.;
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
    if (run>100000) return true;
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
vector<vector<int>> FindRunSublist(string source, vector<int> Target_runlist, double energy)
{
        std::cout << "Getting sublist runs from " << source << std::endl;
        vector<vector<int>> list;
        vector<int> totallist_runnumber;
        vector<double> totallist_elevation;
        vector<double> totallist_azimuth;
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
        double delta_elev = 1.0;
        double delta_azim = 2.0;
        //double delta_elev = 0.1;
        //double delta_azim = 0.2;
        if (energy>300.)
        {
            delta_elev = 5.0;
            delta_azim = 360.;
        }
        if (energy>500.)
        {
            delta_elev = 10.;
            delta_azim = 360.;
        }
        if (energy>700.)
        {
            delta_elev = 10.;
            delta_azim = 360.;
        }
        if (energy>1000.)
        {
            delta_elev = 10.;
            delta_azim = 360.;
        }
        if (energy>1500.)
        {
            delta_elev = 20.;
            delta_azim = 360.;
        }
        if (energy>2000.)
        {
            delta_elev = 40.;
            delta_azim = 360.;
        }
        if (energy>3000.)
        {
            delta_elev = 40.;
            delta_azim = 360.;
        }
        if (energy>4000.)
        {
            delta_elev = 40.;
            delta_azim = 360.;
        }
        if (energy>6000.)
        {
            delta_elev = 40.;
            delta_azim = 360.;
        }
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
            totallist_runnumber.push_back(Target_runlist[run]);
            totallist_elevation.push_back(TelElevation);
            totallist_azimuth.push_back(TelAzimuth);
            input_file->Close();
            continue;
        }
        for (int elev=0;elev<(90./delta_elev);elev++)
        {
            double elev_low = delta_elev*((double)elev)+10.;
            double elev_up = delta_elev*((double)elev+1)+10.;
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
                for (int run=0;run<totallist_runnumber.size();run++) {
                    TelElevation = totallist_elevation.at(run);
                    TelAzimuth = totallist_azimuth.at(run);
                    if (TelElevation<elev_low) 
                    {
                        continue;
                    }
                    if (TelElevation>elev_up)
                    {
                        continue;
                    }
                    if (TelAzimuth<azim_low)
                    {
                        continue;
                    }
                    if (TelAzimuth>azim_up)
                    {
                        continue;
                    }
                    sublist.push_back(totallist_runnumber.at(run));
                    continue;
                }
                if (!(sublist.size()==0)) list.push_back(sublist);
            }
        }
        return list;
}
std::pair <double,double> GetMcGillElectronFlux(double energy)
{
    double McGill_energy[11] = {335,376,422,473,531,596,668,749,840,995,1253};
    double McGill_energy_error[11] = {0,0,0,0,0,0,0,0,0,0,0};
    double McGill_electron_flux[11] = {314,210,149,97.6,71.7,44.8,35.2,23.3,13.8,6.4,2.69};
    double McGill_electron_error[11] = {9.7,8.5,5.7,4.7,3.1,2.5,1.4,1.2,0.72,0.47,0.27};
    TH1D Hist_Flux = TH1D("Hist_Flux","",N_energy_bins,energy_bins);
    for (int e=0;e<N_energy_bins;e++)
    {
        for (int e2=0;e2<10;e2++)
        {
            if (McGill_energy[e2]<energy_bins[e] && McGill_energy[e2+1]>energy_bins[e])
            {
                double delta_x = McGill_energy[e2+1]-McGill_energy[e2];
                double delta_y = McGill_electron_flux[e2+1]-McGill_electron_flux[e2];
                double delta_yerr = McGill_electron_error[e2+1]-McGill_electron_error[e2];
                double new_delta_x = energy_bins[e]-McGill_energy[e2];
                double new_delta_y = delta_y/delta_x*new_delta_x;
                double new_delta_yerr = delta_yerr/delta_x*new_delta_x;
                Hist_Flux.SetBinContent(e+1,McGill_electron_flux[e2]+new_delta_y);
                Hist_Flux.SetBinError(e+1,McGill_electron_error[e2]+new_delta_yerr);
            }
        }
    }
    TF1 *func = new TF1("func",FitPowerLawFunction,energy_bins[0],energy_bins[N_energy_bins-1],2);
    func->SetParameter(0,314);
    func->SetParameter(1,-2.);
    Hist_Flux.Fit("func","","",energy_bins[0],energy_bins[N_energy_bins-1]);
    return std::make_pair(func->Eval(energy),0.);
}
void MLDeconvolutionMethodForExtendedSources(string target_data, int NTelMin, int NTelMax, double elev_lower, double elev_upper, double azim_lower, double azim_upper, double theta2_cut_lower_input, double theta2_cut_upper_input, double MSCW_cut_blind_input, double MSCW_cut_upper_input, double MSCW_cut_lower_input,bool DoConverge, int run_energy_bin, TString WriteType) 
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

        use_this_energy_bin[run_energy_bin] = true;
        //for (int e=0;e<N_energy_bins;e++)
        //{
        //    //if (energy_bins[e]<280) use_this_energy_bin[e] = false;
        //    //if (energy_bins[e]>300) use_this_energy_bin[e] = false;
        //    //if (energy_bins[e]<500) use_this_energy_bin[e] = false;
        //    //if (energy_bins[e]>600) use_this_energy_bin[e] = false;
        //    //if (energy_bins[e]<1000) use_this_energy_bin[e] = false;
        //    //if (energy_bins[e]>1200) use_this_energy_bin[e] = false;
        //    //if (energy_bins[e]<1500) use_this_energy_bin[e] = false;
        //    //if (energy_bins[e]>1800) use_this_energy_bin[e] = false;
        //    if (energy_bins[e]<1800) use_this_energy_bin[e] = false;
        //    if (energy_bins[e]>2000) use_this_energy_bin[e] = false;
        //    //if (energy_bins[e]<2000) use_this_energy_bin[e] = false;
        //    //if (energy_bins[e]>3000) use_this_energy_bin[e] = false;
        //}


#ifndef VEGAS
        Theta2_upper_limit = 10.;
#else
        Theta2_upper_limit = 20.;
#endif

        TRandom rnd;
        TH2D Hist_Target_TelRaDec_AfterCut("Hist_Target_TelRaDec_AfterCut","",100,0,5,100,-1,1);
        TH1D Hist_Target_ON_MSCW_Alpha("Hist_Target_ON_MSCW_Alpha","",100,0,10);
        TH1D Hist_Target_OFF_MSCW_Alpha("Hist_Target_OFF_MSCW_Alpha","",100,0,10);
        TProfile Hist_Measured_Electron_Flux("Hist_Measured_Electron_Flux","",N_energy_bins,energy_bins,-10000.,1000000.);
        TProfile Hist_Target_EndPoint_0("Hist_Target_EndPoint_0","",N_energy_bins,energy_bins,-5.,5.);
        TProfile Hist_Target_EndPoint_1("Hist_Target_EndPoint_1","",N_energy_bins,energy_bins,-5.,5.);
        TH1D Hist_EffAreaTime("Hist_EffAreaTime","",N_energy_bins,energy_bins);
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
        vector<vector<TH1D>> Hist_ElectronMC_SR_MSCW_SumRuns;
        vector<vector<TH1D>> Hist_ElectronMC_CR_MSCW_SumRuns;
        vector<vector<TH1D>> Hist_Scaled_ElectronMC_SR_MSCW_SumRuns;
        vector<vector<TH1D>> Hist_Scaled_ElectronMC_SR_MSCW;
        vector<vector<TH1D>> Hist_Scaled_ElectronMC_CR_MSCW_SumRuns;
        vector<vector<TH1D>> Hist_Scaled_ElectronMC_CR_MSCW;
        vector<vector<TH1D>> Hist_GammaMC_SR_MSCW_SumRuns;
        vector<vector<TH1D>> Hist_GammaMC_CR_MSCW_SumRuns;
        vector<vector<TH1D>> Hist_Scaled_GammaMC_SR_MSCW_SumRuns;
        vector<vector<TH1D>> Hist_Scaled_GammaMC_SR_MSCW;
        vector<vector<TH1D>> Hist_Scaled_GammaMC_CR_MSCW_SumRuns;
        vector<vector<TH1D>> Hist_Scaled_GammaMC_CR_MSCW;
        vector<vector<TH1D>> Hist_Target_CR_MSCW;
        vector<vector<TH1D>> Hist_Target_CR_MSCW_SumRuns;
        vector<TH1D> Hist_Target_Excess_EachRun;
        vector<TH1D> Hist_Target_SR_MSCW_SumRuns_SumSRs;
        vector<TH1D> Hist_Target_CR_MSCW_SumSRs;
        vector<TH1D> Hist_Target_BkgCR_MSCW_SumSRs;
        vector<TH1D> Hist_ElectronMC_SR_MSCW_SumRuns_SumSRs;
        vector<TH1D> Hist_GammaMC_SR_MSCW_SumRuns_SumSRs;
        vector<TH1D> Hist_Target_SR_MSCL;
        vector<TH1D> Hist_Target_Amplitude;
        vector<TH1D> Hist_Target_Mean;
        vector<TH1D> Hist_Target_RMS;
        vector<TH1D> Hist_Target_BkgTemp_MSCW;
        vector<TH1D> Hist_Target_BkgTemp2_MSCW;
        vector<TH1D> Hist_Target_BkgTemp3_MSCW;
        vector<vector<TH1D>> Hist_Target_BkgSR_MSCW;
        vector<double> Weight_Target_BkgSR_MSCW_AllCR;
        vector<vector<TH1D>> Hist_Target_BkgCR_MSCW;
        vector<vector<TH1D>> Hist_Target_BkgSR_MSCW_SumRuns;
        vector<vector<TH1D>> Hist_Target_BkgCR_MSCW_SumRuns;
        vector<TH1D> Hist_Target_BkgSR_MSCW_SumRuns_SumSRs;
        vector<TH1D> Hist_Target_Deconv_MSCW;
        vector<vector<TH1D>> Hist_Dark_SR_MSCW;
        vector<vector<TH1D>> Hist_Dark_SR_MSCW_SumRuns;
        vector<vector<TH1D>> Hist_Dark_CR_MSCW;
        vector<vector<TH1D>> Hist_Dark_CR_MSCW_SumRuns;
        vector<TH1D> Hist_Dark_SR_MSCW_SumRuns_SumSRs;
        vector<TH1D> Hist_Dark_SR_MSCL;
        vector<TH1D> Hist_Dark_BkgTemp_MSCW;
        vector<TH1D> Hist_Dark_BkgPrevious_MSCW;
        vector<vector<TH1D>> Hist_Dark_BkgSR_MSCW;
        vector<vector<TH1D>> Hist_Dark_BkgSR_MSCW_AllCR;
        vector<double> Weight_Dark_BkgSR_MSCW_AllCR;
        vector<vector<TH1D>> Hist_Dark_BkgCR_MSCW;
        vector<vector<TH1D>> Hist_Dark_BkgSR_MSCW_SumRuns;
        vector<vector<TH1D>> Hist_Dark_BkgCR_MSCW_SumRuns;
        vector<TH1D> Hist_Dark_BkgSR_MSCW_SumRuns_SumSRs;
        vector<TH1D> Hist_Dark_Deconv_MSCW;
        vector<TH2D> Hist_Target_MSCLW;
        for (int e=0;e<N_energy_bins;e++) {
            char e_low[50];
            sprintf(e_low, "%i", int(energy_bins[e]));
            char e_up[50];
            sprintf(e_up, "%i", int(energy_bins[e+1]));
            //N_bins_for_deconv = 960;
            if (energy_bins[e]>=200.) N_bins_for_deconv = 960;
            if (energy_bins[e]>=500.) N_bins_for_deconv = 480;
            if (energy_bins[e]>=3200.) N_bins_for_deconv = 240;
            if (UseVegas)
            {
              if (energy_bins[e]>=100.) N_bins_for_deconv = 960;
              if (energy_bins[e]>=1000.) N_bins_for_deconv = 960;
              if (energy_bins[e]>=3000.) N_bins_for_deconv = 960*2;
            }
            Hist_Target_SR_ErecS.push_back(TH1D("Hist_Target_SR_ErecS_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_energy_bins,energy_bins));
            Hist_Target_SR_MSCL.push_back(TH1D("Hist_Target_SR_MSCL_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
            Hist_Target_Amplitude.push_back(TH1D("Hist_Target_Amplitude_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",Number_of_SR+Number_of_CR,MSCL_signal_cut_lower[Number_of_SR-1],MSCL_control_cut_upper[0]));
            Hist_Target_Mean.push_back(TH1D("Hist_Target_Mean_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",Number_of_SR+Number_of_CR,MSCL_signal_cut_lower[Number_of_SR-1],MSCL_control_cut_upper[0]));
            Hist_Target_RMS.push_back(TH1D("Hist_Target_RMS_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",Number_of_SR+Number_of_CR,MSCL_signal_cut_lower[Number_of_SR-1],MSCL_control_cut_upper[0]));
            Hist_Target_BkgTemp_MSCW.push_back(TH1D("Hist_Target_BkgTemp_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
            Hist_Target_BkgTemp2_MSCW.push_back(TH1D("Hist_Target_BkgTemp2_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
            Hist_Target_BkgTemp3_MSCW.push_back(TH1D("Hist_Target_BkgTemp3_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
            Hist_Target_Deconv_MSCW.push_back(TH1D("Hist_Target_Deconv_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
            Hist_Dark_SR_MSCL.push_back(TH1D("Hist_Dark_SR_MSCL_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
            Hist_Dark_BkgTemp_MSCW.push_back(TH1D("Hist_Dark_BkgTemp_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
            Hist_Dark_BkgPrevious_MSCW.push_back(TH1D("Hist_Dark_BkgPrevious_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
            Hist_Dark_Deconv_MSCW.push_back(TH1D("Hist_Dark_Deconv_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
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
            vector<TH1D> Hist_ElectronMC_ThisE_SR_MSCW_SumRuns;
            vector<TH1D> Hist_ElectronMC_ThisE_CR_MSCW_SumRuns;
            vector<TH1D> Hist_Scaled_ElectronMC_ThisE_SR_MSCW_SumRuns;
            vector<TH1D> Hist_Scaled_ElectronMC_ThisE_CR_MSCW_SumRuns;
            vector<TH1D> Hist_Scaled_ElectronMC_ThisE_SR_MSCW;
            vector<TH1D> Hist_Scaled_ElectronMC_ThisE_CR_MSCW;
            vector<TH1D> Hist_GammaMC_ThisE_SR_MSCW_SumRuns;
            vector<TH1D> Hist_GammaMC_ThisE_CR_MSCW_SumRuns;
            vector<TH1D> Hist_Scaled_GammaMC_ThisE_SR_MSCW_SumRuns;
            vector<TH1D> Hist_Scaled_GammaMC_ThisE_CR_MSCW_SumRuns;
            vector<TH1D> Hist_Scaled_GammaMC_ThisE_SR_MSCW;
            vector<TH1D> Hist_Scaled_GammaMC_ThisE_CR_MSCW;
            vector<TH1D> Hist_Target_ThisE_BkgSR_MSCW;
            vector<TH1D> Hist_Target_ThisE_BkgSR_MSCW_AllCR;
            vector<TH1D> Hist_Target_ThisE_BkgSR_MSCW_SumRuns;
            vector<TH1D> Hist_Target_ThisE_CR_MSCW;
            vector<TH1D> Hist_Target_ThisE_CR_MSCW_SumRuns;
            vector<TH1D> Hist_Target_ThisE_BkgCR_MSCW;
            vector<TH1D> Hist_Target_ThisE_BkgCR_MSCW_SumRuns;
            vector<TH1D> Hist_Dark_ThisE_SR_MSCW;
            vector<TH1D> Hist_Dark_ThisE_SR_MSCW_SumRuns;
            vector<TH1D> Hist_Dark_ThisE_BkgSR_MSCW;
            vector<TH1D> Hist_Dark_ThisE_BkgSR_MSCW_AllCR;
            vector<TH1D> Hist_Dark_ThisE_BkgSR_MSCW_SumRuns;
            vector<TH1D> Hist_Dark_ThisE_CR_MSCW;
            vector<TH1D> Hist_Dark_ThisE_CR_MSCW_SumRuns;
            vector<TH1D> Hist_Dark_ThisE_BkgCR_MSCW;
            vector<TH1D> Hist_Dark_ThisE_BkgCR_MSCW_SumRuns;
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
                Hist_ElectronMC_ThisE_SR_MSCW_SumRuns.push_back(TH1D("Hist_ElectronMC_SR"+TString(nsr)+"_MSCW_SumRuns_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
                Hist_Scaled_ElectronMC_ThisE_SR_MSCW_SumRuns.push_back(TH1D("Hist_Scaled_ElectronMC_SR"+TString(nsr)+"_MSCW_SumRuns_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
                Hist_Scaled_ElectronMC_ThisE_SR_MSCW.push_back(TH1D("Hist_Scaled_ElectronMC_SR"+TString(nsr)+"_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
                Hist_GammaMC_ThisE_SR_MSCW_SumRuns.push_back(TH1D("Hist_GammaMC_SR"+TString(nsr)+"_MSCW_SumRuns_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
                Hist_Scaled_GammaMC_ThisE_SR_MSCW_SumRuns.push_back(TH1D("Hist_Scaled_GammaMC_SR"+TString(nsr)+"_MSCW_SumRuns_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
                Hist_Scaled_GammaMC_ThisE_SR_MSCW.push_back(TH1D("Hist_Scaled_GammaMC_SR"+TString(nsr)+"_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
                Hist_Target_ThisE_BkgSR_MSCW.push_back(TH1D("Hist_Target_BkgSR"+TString(nsr)+"_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
                Hist_Target_ThisE_BkgSR_MSCW_AllCR.push_back(TH1D("Hist_Target_BkgSR"+TString(nsr)+"_MSCW_AllCR_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
                Hist_Target_ThisE_BkgSR_MSCW_SumRuns.push_back(TH1D("Hist_Target_BkgSR"+TString(nsr)+"_MSCW_SumRuns_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
                Hist_Dark_ThisE_SR_MSCW.push_back(TH1D("Hist_Dark_SR"+TString(nsr)+"_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
                Hist_Dark_ThisE_SR_MSCW_SumRuns.push_back(TH1D("Hist_Dark_SR"+TString(nsr)+"_MSCW_SumRuns_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
                Hist_Dark_ThisE_BkgSR_MSCW.push_back(TH1D("Hist_Dark_BkgSR"+TString(nsr)+"_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
                Hist_Dark_ThisE_BkgSR_MSCW_AllCR.push_back(TH1D("Hist_Dark_BkgSR"+TString(nsr)+"_MSCW_AllCR_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
                Hist_Dark_ThisE_BkgSR_MSCW_SumRuns.push_back(TH1D("Hist_Dark_BkgSR"+TString(nsr)+"_MSCW_SumRuns_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
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
            Hist_ElectronMC_SR_MSCW_SumRuns.push_back(Hist_ElectronMC_ThisE_SR_MSCW_SumRuns);
            Hist_Scaled_ElectronMC_SR_MSCW_SumRuns.push_back(Hist_Scaled_ElectronMC_ThisE_SR_MSCW_SumRuns);
            Hist_Scaled_ElectronMC_SR_MSCW.push_back(Hist_Scaled_ElectronMC_ThisE_SR_MSCW);
            Hist_GammaMC_SR_MSCW_SumRuns.push_back(Hist_GammaMC_ThisE_SR_MSCW_SumRuns);
            Hist_Scaled_GammaMC_SR_MSCW_SumRuns.push_back(Hist_Scaled_GammaMC_ThisE_SR_MSCW_SumRuns);
            Hist_Scaled_GammaMC_SR_MSCW.push_back(Hist_Scaled_GammaMC_ThisE_SR_MSCW);
            Hist_Target_BkgSR_MSCW.push_back(Hist_Target_ThisE_BkgSR_MSCW);
            Hist_Target_BkgSR_MSCW_SumRuns.push_back(Hist_Target_ThisE_BkgSR_MSCW_SumRuns);
            Hist_Dark_SR_MSCW.push_back(Hist_Dark_ThisE_SR_MSCW);
            Hist_Dark_SR_MSCW_SumRuns.push_back(Hist_Dark_ThisE_SR_MSCW_SumRuns);
            Hist_Dark_BkgSR_MSCW.push_back(Hist_Dark_ThisE_BkgSR_MSCW);
            Hist_Dark_BkgSR_MSCW_AllCR.push_back(Hist_Dark_ThisE_BkgSR_MSCW_AllCR);
            Hist_Dark_BkgSR_MSCW_SumRuns.push_back(Hist_Dark_ThisE_BkgSR_MSCW_SumRuns);
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
            Hist_Target_Excess_EachRun.push_back(TH1D("Hist_Target_Excess_EachRun_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
            Hist_Target_SR_MSCW_SumRuns_SumSRs.push_back(TH1D("Hist_Target_SR_MSCW_SumRuns_SumSRs_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
            Hist_Target_CR_MSCW_SumSRs.push_back(TH1D("Hist_Target_CR_MSCW_SumSRs_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
            Hist_Target_BkgCR_MSCW_SumSRs.push_back(TH1D("Hist_Target_BkgCR_MSCW_SumSRs_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
            Hist_ElectronMC_SR_MSCW_SumRuns_SumSRs.push_back(TH1D("Hist_ElectronMC_SR_MSCW_SumRuns_SumSRs_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
            Hist_GammaMC_SR_MSCW_SumRuns_SumSRs.push_back(TH1D("Hist_GammaMC_SR_MSCW_SumRuns_SumSRs_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
            Hist_Target_BkgSR_MSCW_SumRuns_SumSRs.push_back(TH1D("Hist_Target_BkgSR_MSCW_SumRuns_SumSRs_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
            Hist_Dark_SR_MSCW_SumRuns_SumSRs.push_back(TH1D("Hist_Dark_SR_MSCW_SumRuns_SumSRs_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
            Hist_Dark_BkgSR_MSCW_SumRuns_SumSRs.push_back(TH1D("Hist_Dark_BkgSR_MSCW_SumRuns_SumSRs_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
            for (int s=0;s<Number_of_CR;s++)
            {
                char nCR[50];
                sprintf(nCR, "%i", int(s+1));
                Hist_ElectronMC_ThisE_CR_MSCW_SumRuns.push_back(TH1D("Hist_ElectronMC_CR"+TString(nCR)+"_MSCW_SumRuns_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
                Hist_Scaled_ElectronMC_ThisE_CR_MSCW_SumRuns.push_back(TH1D("Hist_Scaled_ElectronMC_CR"+TString(nCR)+"_MSCW_SumRuns_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
                Hist_Scaled_ElectronMC_ThisE_CR_MSCW.push_back(TH1D("Hist_Scaled_ElectronMC_CR"+TString(nCR)+"_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
                Hist_GammaMC_ThisE_CR_MSCW_SumRuns.push_back(TH1D("Hist_GammaMC_CR"+TString(nCR)+"_MSCW_SumRuns_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
                Hist_Scaled_GammaMC_ThisE_CR_MSCW_SumRuns.push_back(TH1D("Hist_Scaled_GammaMC_CR"+TString(nCR)+"_MSCW_SumRuns_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
                Hist_Scaled_GammaMC_ThisE_CR_MSCW.push_back(TH1D("Hist_Scaled_GammaMC_CR"+TString(nCR)+"_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
                Hist_Target_ThisE_CR_MSCW.push_back(TH1D("Hist_Target_CR"+TString(nCR)+"_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
                Hist_Target_ThisE_CR_MSCW_SumRuns.push_back(TH1D("Hist_Target_CR"+TString(nCR)+"_MSCW_SumRuns_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
                Hist_Target_ThisE_BkgCR_MSCW.push_back(TH1D("Hist_Target_BkgCR"+TString(nCR)+"_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
                Hist_Target_ThisE_BkgCR_MSCW_SumRuns.push_back(TH1D("Hist_Target_BkgCR"+TString(nCR)+"_MSCW_SumRuns_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
                Hist_Dark_ThisE_CR_MSCW.push_back(TH1D("Hist_Dark_CR"+TString(nCR)+"_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
                Hist_Dark_ThisE_CR_MSCW_SumRuns.push_back(TH1D("Hist_Dark_CR"+TString(nCR)+"_MSCW_SumRuns_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
                Hist_Dark_ThisE_BkgCR_MSCW.push_back(TH1D("Hist_Dark_BkgCR"+TString(nCR)+"_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
                Hist_Dark_ThisE_BkgCR_MSCW_SumRuns.push_back(TH1D("Hist_Dark_BkgCR"+TString(nCR)+"_MSCW_SumRuns_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
            }
            Hist_ElectronMC_CR_MSCW_SumRuns.push_back(Hist_ElectronMC_ThisE_CR_MSCW_SumRuns);
            Hist_Scaled_ElectronMC_CR_MSCW_SumRuns.push_back(Hist_Scaled_ElectronMC_ThisE_CR_MSCW_SumRuns);
            Hist_Scaled_ElectronMC_CR_MSCW.push_back(Hist_Scaled_ElectronMC_ThisE_CR_MSCW);
            Hist_GammaMC_CR_MSCW_SumRuns.push_back(Hist_GammaMC_ThisE_CR_MSCW_SumRuns);
            Hist_Scaled_GammaMC_CR_MSCW_SumRuns.push_back(Hist_Scaled_GammaMC_ThisE_CR_MSCW_SumRuns);
            Hist_Scaled_GammaMC_CR_MSCW.push_back(Hist_Scaled_GammaMC_ThisE_CR_MSCW);
            Hist_Target_CR_MSCW.push_back(Hist_Target_ThisE_CR_MSCW);
            Hist_Target_CR_MSCW_SumRuns.push_back(Hist_Target_ThisE_CR_MSCW_SumRuns);
            Hist_Target_BkgCR_MSCW.push_back(Hist_Target_ThisE_BkgCR_MSCW);
            Hist_Target_BkgCR_MSCW_SumRuns.push_back(Hist_Target_ThisE_BkgCR_MSCW_SumRuns);
            Hist_Dark_CR_MSCW.push_back(Hist_Dark_ThisE_CR_MSCW);
            Hist_Dark_CR_MSCW_SumRuns.push_back(Hist_Dark_ThisE_CR_MSCW_SumRuns);
            Hist_Dark_BkgCR_MSCW.push_back(Hist_Dark_ThisE_BkgCR_MSCW);
            Hist_Dark_BkgCR_MSCW_SumRuns.push_back(Hist_Dark_ThisE_BkgCR_MSCW_SumRuns);
        }


        std::cout << "Getting dark pointing distributions... " << std::endl;
        vector<int> Dark_runlist = GetRunList("DarkField");
        if (TString(target)=="Proton") Dark_runlist = GetRunList("Proton");
        char Dark_observation[50];
        sprintf(Dark_observation, "%s", "DarkField");
        if (TString(target)=="Proton") sprintf(Dark_observation, "%s", "Proton");
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

        vector<int> N_iter;
        for (int e=0;e<N_energy_bins;e++) {
                N_iter.push_back(30);
        }

        TF1 *myfunc = new TF1("myfunc",Kernel,-50.,50.,1);
        myfunc->SetParameter(0,0.5);

        for (int e=0;e<N_energy_bins;e++)
        {
            vector<double> thisE_dark_kernel_rms;
            vector<double> thisE_dark_kernel_shift;
            vector<double> thisE_target_kernel_rms;
            vector<double> thisE_target_kernel_shift;
            for (int s=0;s<Number_of_SR;s++)
            {
                thisE_dark_kernel_rms.push_back(0.);
                thisE_dark_kernel_shift.push_back(0.);
                thisE_target_kernel_rms.push_back(0.);
                thisE_target_kernel_shift.push_back(0.);
            }
            dark_converge.push_back(std::make_pair(0.,0.));
            dark_kernel_rms.push_back(thisE_dark_kernel_rms);
            dark_kernel_shift.push_back(thisE_dark_kernel_shift);
            target_kernel_rms.push_back(thisE_target_kernel_rms);
            target_kernel_shift.push_back(thisE_target_kernel_shift);
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
                        for (int s=0;s<Number_of_SR;s++)
                        {
                            if (SignalSelectionMSCW(s)) Hist_Dark_SR_MSCW.at(e).at(s).Fill(MSCW);
                        }
                        for (int s=0;s<Number_of_CR;s++)
                        {
                            if (ControlSelectionMSCW(s)) Hist_Dark_CR_MSCW.at(e).at(s).Fill(MSCW);
                        }
                }
            }
            input_file->Close();
            //break;  // for test purpose
        }

        std::cout << "Getting MC runs... " << std::endl;
        vector<int> MC_runlist = GetRunList("Photon");
        for (int run=0;run<MC_runlist.size();run++)
        {
            char run_number[50];
            sprintf(run_number, "%i", int(MC_runlist[run]));
            //std::cout << "Reading run " << run_number << std::endl;
            sprintf(Dark_observation, "%s", "Photon");
            filename = TString("$VERITAS_USER_DATA_DIR/"+TString(Dark_observation)+"_V6_Moderate-TMVA-BDT.RB."+TString(run_number)+".root");
            TFile*  input_file = TFile::Open(filename.c_str());
            TString root_file = "run_"+TString(run_number)+"/stereo/data_on";
            TTree* MC_tree = (TTree*) input_file->Get(root_file);
            MC_tree->SetBranchAddress("ErecS",&ErecS);
            MC_tree->SetBranchAddress("EChi2S",&EChi2S);
            MC_tree->SetBranchAddress("MSCW",&MSCW);
            MC_tree->SetBranchAddress("MSCL",&MSCL);
            MC_tree->SetBranchAddress("NImages",&NImages);
            for (int entry=0;entry<MC_tree->GetEntries();entry++) {
                ErecS = 0;
                EChi2S = 0;
                NImages = 0;
                MSCW = 0;
                MSCL = 0;
                MC_tree->GetEntry(entry);
                int energy = Hist_Target_SR_ErecS.at(0).FindBin(ErecS*1000.)-1;
                if (energy<0) continue;
                if (energy>=N_energy_bins) continue;
                int e = energy;
                if (!SelectNImages(NTelMin,NTelMax)) continue;
                for (int s=0;s<Number_of_SR;s++)
                {
                    if (SignalSelectionMSCW(s)) 
                    {
                        Hist_ElectronMC_SR_MSCW_SumRuns.at(e).at(s).Fill(MSCW+0.2);
                        Hist_GammaMC_SR_MSCW_SumRuns.at(e).at(s).Fill(MSCW+0.2);
                    }
                }
                for (int s=0;s<Number_of_CR;s++)
                {
                    if (ControlSelectionMSCW(s)) 
                    {
                        Hist_ElectronMC_CR_MSCW_SumRuns.at(e).at(s).Fill(MSCW+0.2);
                        Hist_GammaMC_CR_MSCW_SumRuns.at(e).at(s).Fill(MSCW+0.2);
                    }
                }
            }
            input_file->Close();
        }

        std::cout << "Getting target runs... " << std::endl;

        for (int s=0;s<Number_of_SR;s++)
        {
                Weight_Target_BkgSR_MSCW_AllCR.push_back(0.);
        }


        MSCW_cut_blind = MSCW_cut_blind_input;
        vector<vector<int>> Sublist;
        //Target_runlist = SortRunListByElevation(target,Target_runlist);
        for (int e=0;e<N_energy_bins;e++)
        {
            if (!use_this_energy_bin[e]) continue;
            if (TString(target)=="Proton")
            {
                Sublist.clear();
                Sublist.push_back(Target_runlist);
            }
            else
            {
                Sublist = FindRunSublist(target,Target_runlist,energy_bins[e]);
            }
            std::cout << "=================================================================" << std::endl;
            std::cout << "Target, e " << energy_bins[e] << std::endl;
            std::cout << "Sublist.size() = " << Sublist.size() << std::endl;
            Hist_Target_SR_MSCW_SumRuns_SumSRs.at(e).Reset();
            Hist_Target_BkgSR_MSCW_SumRuns_SumSRs.at(e).Reset();
            for (int s=0;s<Number_of_SR;s++)
            {
                Hist_Target_SR_MSCW_SumRuns.at(e).at(s).Reset();
                Hist_Target_BkgSR_MSCW_SumRuns.at(e).at(s).Reset();
            }
            for (int s=0;s<Number_of_CR;s++)
            {
                Hist_Target_CR_MSCW_SumRuns.at(e).at(s).Reset();
                Hist_Target_BkgCR_MSCW_SumRuns.at(e).at(s).Reset();
            }
            for (int subrun=0;subrun<Sublist.size();subrun++)
            {
                std::cout << subrun << "/" << Sublist.size() << " completed." << std::endl;
                for (int s=0;s<Number_of_SR;s++)
                {
                    Hist_Target_SR_MSCW.at(e).at(s).Reset();
                    Hist_Target_BkgSR_MSCW.at(e).at(s).Reset();
                }
                for (int s=0;s<Number_of_CR;s++)
                {
                    Hist_Target_CR_MSCW.at(e).at(s).Reset();
                    Hist_Target_BkgCR_MSCW.at(e).at(s).Reset();
                }
                if (Sublist.at(subrun).size()==0) continue;
                Hist_EffAreaTime.Reset();
                electron_count[e] = 0.;
                electron_count_err[e] = 0.;
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
		    TH1* i_hEffAreaP = ( TH1* )getEffAreaHistogram(input_file,Sublist.at(subrun)[run]);
                    double eff_area = i_hEffAreaP->GetBinContent( i_hEffAreaP->FindBin( log10(0.5*(energy_bins[e]+energy_bins[e+1])/1000.)));
                    std::cout << "eff_area = " << eff_area << std::endl; 
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
                    std::pair <double,double> mcgillflux = GetMcGillElectronFlux((energy_bins[e+1]+energy_bins[e])/2.);
                    electron_flux[e] = mcgillflux.first;
                    electron_flux_err[e] = mcgillflux.second;
                    std::cout << "electron_flux[e] = " << electron_flux[e] << std::endl;
                    double expected_electrons = 1e-12*10000.*electron_flux[e]*eff_area*(time_1-time_0)*(energy_bins[e+1]-energy_bins[e])/1000.;
                    double expected_electrons_err = 1e-12*10000.*electron_flux_err[e]*eff_area*(time_1-time_0)*(energy_bins[e+1]-energy_bins[e])/1000.;
                    electron_count[e] += expected_electrons;
                    electron_count_err[e] += expected_electrons_err; // correlated error, add in linear
                    Hist_EffAreaTime.Fill(energy_bins[e],eff_area*(time_1-time_0));
                    bool run_already_used = false;
                    for (int this_run=0;this_run<used_runs.size();this_run++)
                    {
                        if (used_runs.at(this_run)==Sublist.at(subrun)[run]) run_already_used = true;
                    }
                    if (!run_already_used) exposure_hours += (time_1-time_0)/3600.;
                    if (!run_already_used) used_runs.push_back(int(Sublist.at(subrun)[run]));
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
                double gamma_integral = 0.;
                double hadron_integral = 0.;
                for (int s=0;s<Number_of_SR;s++)
                {
                    gamma_integral += Hist_GammaMC_SR_MSCW_SumRuns.at(e).at(s).Integral();
                    hadron_integral += Hist_Target_SR_MSCW.at(e).at(s).Integral();
                }
                double fake_signal_scale = 0.05*hadron_integral/gamma_integral;
                for (int s=0;s<Number_of_SR;s++)
                {
                    for (int bin=0;bin<Hist_GammaMC_SR_MSCW_SumRuns.at(e).at(s).GetNbinsX();bin++)
                    {
                        double old_content = Hist_GammaMC_SR_MSCW_SumRuns.at(e).at(s).GetBinContent(bin+1);
                        double old_error = Hist_GammaMC_SR_MSCW_SumRuns.at(e).at(s).GetBinError(bin+1);
                        double new_content = old_content*fake_signal_scale;
                        Hist_Scaled_GammaMC_SR_MSCW.at(e).at(s).SetBinContent(bin+1,new_content);
                        Hist_Scaled_GammaMC_SR_MSCW.at(e).at(s).SetBinError(bin+1,pow(new_content,0.5));
                    }
                }
                for (int s=0;s<Number_of_SR;s++)
                {
                    if (TString(target)=="Proton" && theta2_cut_lower_input==0)
                    //if (theta2_cut_lower_input==0)
                    {
                        Hist_Target_SR_MSCW.at(e).at(s).Add(&Hist_Scaled_GammaMC_SR_MSCW.at(e).at(s));
                    }
                    Hist_Target_SR_MSCW_SumRuns.at(e).at(s).Add(&Hist_Target_SR_MSCW.at(e).at(s));
                    Hist_Target_SR_MSCW_SumRuns_SumSRs.at(e).Add(&Hist_Target_SR_MSCW.at(e).at(s));
                }
                for (int c=0;c<Number_of_CR;c++)
                {
                    Hist_Target_CR_MSCW_SumRuns.at(e).at(c).Add(&Hist_Target_CR_MSCW.at(e).at(c));
                }

                std::cout << "Normalizing MC histograms..." << std::endl;
                TH1D Hist_ElectronMC_Temp = TH1D("Hist_ElectronMC_Temp","",Hist_Target_SR_MSCW_SumRuns_SumSRs.at(e).GetNbinsX(),MSCW_plot_lower,MSCW_plot_upper);
                int norm_bin_low = Hist_Target_SR_MSCW_SumRuns_SumSRs.at(e).FindBin(MSCW_cut_lower);
                int norm_bin_up = Hist_Target_SR_MSCW_SumRuns_SumSRs.at(e).FindBin(MSCW_cut_blind);
                std::cout << "Energy "  << energy_bins[e] << std::endl;
                std::cout << "predicted electrons = " << electron_count[e] << std::endl;
                double old_integral = 0.;
                for (int s=0;s<Number_of_SR;s++)
                {
                    old_integral += double(Hist_ElectronMC_SR_MSCW_SumRuns.at(e).at(s).Integral());
                }
                for (int s=0;s<Number_of_CR;s++)
                {
                    old_integral += double(Hist_ElectronMC_CR_MSCW_SumRuns.at(e).at(s).Integral());
                }
                double scale = electron_count[e]/old_integral;
                double scale_err = electron_count_err[e]/old_integral;
                if (TString(target)=="Proton") scale = 0.;
                if (TString(target)=="Proton") scale_err = 0.;
                for (int s=0;s<Number_of_SR;s++)
                {
                    for (int bin=0;bin<Hist_ElectronMC_SR_MSCW_SumRuns.at(e).at(s).GetNbinsX();bin++)
                    {
                        double old_content = Hist_ElectronMC_SR_MSCW_SumRuns.at(e).at(s).GetBinContent(bin+1);
                        double old_error = Hist_ElectronMC_SR_MSCW_SumRuns.at(e).at(s).GetBinError(bin+1);
                        double new_content = old_content*scale;
                        double new_error = old_error*old_error*scale*scale+(old_content*scale_err)*(old_content*scale_err);
                        if (new_error>0) new_error = pow(new_error,0.5);
                        Hist_Scaled_ElectronMC_SR_MSCW.at(e).at(s).SetBinContent(bin+1,new_content);
                        //Hist_Scaled_ElectronMC_SR_MSCW.at(e).at(s).SetBinError(bin+1,pow(new_content,0.5));
                        Hist_Scaled_ElectronMC_SR_MSCW.at(e).at(s).SetBinError(bin+1,0.3*new_content);
                    }
                    Hist_Scaled_ElectronMC_SR_MSCW_SumRuns.at(e).at(s).Add(&Hist_Scaled_ElectronMC_SR_MSCW.at(e).at(s));
                    Hist_ElectronMC_SR_MSCW_SumRuns_SumSRs.at(e).Add(&Hist_Scaled_ElectronMC_SR_MSCW.at(e).at(s));
                }
                for (int s=0;s<Number_of_CR;s++)
                {
                    for (int bin=0;bin<Hist_ElectronMC_CR_MSCW_SumRuns.at(e).at(s).GetNbinsX();bin++)
                    {
                        double old_content = Hist_ElectronMC_CR_MSCW_SumRuns.at(e).at(s).GetBinContent(bin+1);
                        double old_error = Hist_ElectronMC_CR_MSCW_SumRuns.at(e).at(s).GetBinError(bin+1);
                        double new_content = old_content*scale;
                        double new_error = old_error*old_error*scale*scale+(old_content*scale_err)*(old_content*scale_err);
                        if (new_error>0) new_error = pow(new_error,0.5);
                        Hist_Scaled_ElectronMC_CR_MSCW.at(e).at(s).SetBinContent(bin+1,new_content);
                        Hist_Scaled_ElectronMC_CR_MSCW.at(e).at(s).SetBinError(bin+1,pow(new_content,0.5));
                    }
                    Hist_Scaled_ElectronMC_CR_MSCW_SumRuns.at(e).at(s).Add(&Hist_Scaled_ElectronMC_CR_MSCW.at(e).at(s));
                }


                std::cout << "=================================================================" << std::endl;
                N_iter.at(e) = 10;
                std::cout << "Target, e " << energy_bins[e] << ", N_iter.at(e) = " << N_iter.at(e) << std::endl;
                std::cout << "Hist_Target_SR_MSCW.at(e).Integral() = " << Hist_Target_SR_MSCW.at(e).at(0).Integral() << std::endl;
                std::cout << "Hist_Target_SR_MSCW.at(e).GetMean() = " << Hist_Target_SR_MSCW.at(e).at(0).GetMean() << std::endl;
                std::cout << "Hist_Target_SR_MSCW.at(e).GetRMS() = " << Hist_Target_SR_MSCW.at(e).at(0).GetRMS() << std::endl;

                std::pair <bool,std::pair <double,double>> offset;
                std::pair <double,double> kernel_rms;
                Hist_Target_Mean.at(e).Reset();
                Hist_Target_RMS.at(e).Reset();
                for (int c=0;c<Number_of_CR;c++)
                {
                    int bin = Hist_Target_Mean.at(e).FindBin(0.5*(MSCL_control_cut_lower[c]+MSCL_control_cut_upper[c]));
                    Hist_Target_Mean.at(e).SetBinContent(bin,Hist_Target_CR_MSCW.at(e).at(c).GetMean());
                    Hist_Target_RMS.at(e).SetBinContent(bin,Hist_Target_CR_MSCW.at(e).at(c).GetRMS());
                    Hist_Target_Mean.at(e).SetBinError(bin,MSCW_cut_blind-MSCW_cut_lower);
                    Hist_Target_RMS.at(e).SetBinError(bin,MSCW_cut_blind-MSCW_cut_lower);
                }
                FindSRMean(&Hist_Target_Mean.at(e),energy_bins[e]);
                FindSRMean(&Hist_Target_RMS.at(e),energy_bins[e]);

                double initial_RMS = 0;
                double initial_shift = 0;
                double rms_begin = 0;
                double mean_begin = 0;

                Hist_Target_CR_MSCW_SumSRs.at(e).Reset();
                Hist_Target_BkgCR_MSCW_SumSRs.at(e).Reset();
                for (int c=0;c<Number_of_CR;c++)
                {
                    Hist_Target_CR_MSCW_SumSRs.at(e).Add(&Hist_Target_CR_MSCW.at(e).at(c));
                }

                bool useOldSR = false; bool doConverge = true;


                // estimate unblinded CR bkg
                for (int c=0;c<Number_of_CR;c++)
                {
                    Hist_Target_BkgCR_MSCW.at(e).at(c).Reset();
                }
                for (int c=1;c<Number_of_CR;c++)
                {

                    std::cout << "Target, e " << energy_bins[e] << ", running CR " << c << " to CR " << c << std::endl;

                    int bin = Hist_Target_Mean.at(e).FindBin(0.5*(MSCL_control_cut_lower[c]+MSCL_control_cut_upper[c]));
                    estimated_mean = Hist_Target_Mean.at(e).GetBinContent(bin);
                    estimated_rms = Hist_Target_RMS.at(e).GetBinContent(bin);
                    bin = Hist_Target_Mean.at(e).FindBin(0.5*(MSCL_control_cut_lower[c-1]+MSCL_control_cut_upper[c-1]));
                    estimated_rms_previous = Hist_Target_RMS.at(e).GetBinContent(bin);
                    bin = Hist_Target_Mean.at(e).FindBin(0.5*(MSCL_control_cut_lower[c]+MSCL_control_cut_upper[c]));

                    useOldSR = true; doConverge = true;
                    if (c>1) useOldSR = false;
                    Hist_Scaled_GammaMC_CR_MSCW.at(e).at(c).Reset();
                    PredictNextLayer(&Hist_Scaled_GammaMC_CR_MSCW.at(e).at(c),&Hist_Scaled_ElectronMC_CR_MSCW.at(e).at(c),&Hist_Target_CR_MSCW.at(e).at(c),&Hist_Target_CR_MSCW.at(e).at(c-1),&Hist_Target_BkgCR_MSCW.at(e).at(c-1),&Hist_Target_BkgCR_MSCW.at(e).at(c),&Hist_Dark_CR_MSCW.at(e).at(c),&Hist_Dark_CR_MSCW.at(e).at(c-1),energy_bins[e],estimated_parameters,useOldSR,doConverge);

                    Hist_Target_BkgCR_MSCW_SumRuns.at(e).at(c).Add(&Hist_Target_BkgCR_MSCW.at(e).at(c));

                }

                std::cout << "Target, e " << energy_bins[e] << ", running CR " << Number_of_CR-1 << " to SR " << 0 << std::endl;
                Hist_Target_BkgSR_MSCW.at(e).at(0).Reset();

                int bin;
                bin = Hist_Target_Mean.at(e).FindBin(0.5*(MSCL_signal_cut_lower[0]+MSCL_signal_cut_upper[0]));
                estimated_mean_err = Hist_Target_RMS.at(e).GetBinContent(bin)/2.;
                estimated_amplitude = Hist_Target_Amplitude.at(e).GetBinContent(bin);
                estimated_mean = Hist_Target_Mean.at(e).GetBinContent(bin);
                estimated_rms = Hist_Target_RMS.at(e).GetBinContent(bin);
                estimated_rms_err = Hist_Target_RMS.at(e).GetBinError(bin);
                bin = Hist_Target_Mean.at(e).FindBin(0.5*(MSCL_control_cut_lower[Number_of_CR-1]+MSCL_control_cut_upper[Number_of_CR-1]));
                estimated_rms_previous = Hist_Target_RMS.at(e).GetBinContent(bin);
                bin = Hist_Target_Mean.at(e).FindBin(0.5*(MSCL_signal_cut_lower[0]+MSCL_signal_cut_upper[0]));
                if (TString(target)=="Proton")
                {
                    estimated_mean = Hist_Target_SR_MSCW.at(e).at(0).GetMean();
                    estimated_rms = Hist_Target_SR_MSCW.at(e).at(0).GetRMS();
                }

                useOldSR = true; doConverge = true;
                PredictNextLayer(&Hist_Scaled_GammaMC_SR_MSCW.at(e).at(0),&Hist_Scaled_ElectronMC_SR_MSCW.at(e).at(0),&Hist_Target_SR_MSCW.at(e).at(0),&Hist_Target_CR_MSCW.at(e).at(Number_of_CR-1),&Hist_Target_BkgCR_MSCW.at(e).at(Number_of_CR-1),&Hist_Target_BkgSR_MSCW.at(e).at(0),&Hist_Dark_SR_MSCW.at(e).at(0),&Hist_Dark_CR_MSCW.at(e).at(Number_of_CR-1),energy_bins[e],estimated_parameters,useOldSR,doConverge);

                //AddBkgStatistics(&Hist_Target_BkgSR_MSCW.at(e).at(0));
                //AddSystematics(&Hist_Target_SR_MSCW.at(e).at(0),&Hist_Target_BkgSR_MSCW.at(e).at(0));
                //AddSystematics2(&Hist_Target_CR_MSCW.at(e).at(Number_of_CR-1),&Hist_Target_BkgCR_MSCW.at(e).at(Number_of_CR-1),&Hist_Target_BkgSR_MSCW.at(e).at(0));
                //AddSystematics3(&Hist_Target_BkgSR_MSCW.at(e).at(0));


                // estimate SR bkg
                for (int s=1;s<Number_of_SR;s++)
                {

                    std::cout << "===============================================================" << std::endl;
                    std::cout << "Target, e " << energy_bins[e] << ", running CR " << Number_of_CR-1 << " to SR " << s << std::endl;
                    Hist_Target_BkgSR_MSCW.at(e).at(s).Reset();

                    bin = Hist_Target_Mean.at(e).FindBin(0.5*(MSCL_signal_cut_lower[s]+MSCL_signal_cut_upper[s]));
                    estimated_amplitude = Hist_Target_Amplitude.at(e).GetBinContent(bin);
                    estimated_mean = Hist_Target_Mean.at(e).GetBinContent(bin);
                    estimated_rms = Hist_Target_RMS.at(e).GetBinContent(bin);
                    estimated_rms_err = Hist_Target_RMS.at(e).GetBinError(bin);
                    bin = Hist_Target_Mean.at(e).FindBin(0.5*(MSCL_signal_cut_lower[s-1]+MSCL_signal_cut_upper[s-1]));
                    estimated_rms_previous = Hist_Target_RMS.at(e).GetBinContent(bin);
                    bin = Hist_Target_Mean.at(e).FindBin(0.5*(MSCL_signal_cut_lower[s]+MSCL_signal_cut_upper[s]));
                    if (TString(target)=="Proton")
                    {
                        estimated_mean = Hist_Target_SR_MSCW.at(e).at(s).GetMean();
                        estimated_rms = Hist_Target_SR_MSCW.at(e).at(s).GetRMS();
                        estimated_rms_previous = Hist_Target_SR_MSCW.at(e).at(s-1).GetRMS();
                    }

                    useOldSR = false; doConverge = true;
                    PredictNextLayer(&Hist_Scaled_GammaMC_SR_MSCW.at(e).at(s),&Hist_Scaled_ElectronMC_SR_MSCW.at(e).at(s),&Hist_Target_SR_MSCW.at(e).at(s),&Hist_Target_SR_MSCW.at(e).at(s-1),&Hist_Target_BkgSR_MSCW.at(e).at(s-1),&Hist_Target_BkgSR_MSCW.at(e).at(s),&Hist_Dark_SR_MSCW.at(e).at(s),&Hist_Dark_SR_MSCW.at(e).at(s-1),energy_bins[e],estimated_parameters,useOldSR,doConverge);

                    //AddBkgStatistics(&Hist_Target_BkgSR_MSCW.at(e).at(s));
                    //AddSystematics(&Hist_Target_SR_MSCW.at(e).at(s),&Hist_Target_BkgSR_MSCW.at(e).at(s));
                    //AddSystematics2(&Hist_Target_CR_MSCW.at(e).at(Number_of_CR-1),&Hist_Target_BkgCR_MSCW.at(e).at(Number_of_CR-1),&Hist_Target_BkgSR_MSCW.at(e).at(s));
                    //AddSystematics3(&Hist_Target_BkgSR_MSCW.at(e).at(s));


                }


                // here we calculate the measured e/gamma flux
                Hist_Target_Excess_EachRun.at(e).Reset();
                for (int s=0;s<Number_of_SR;s++)
                {
                    Hist_Target_Excess_EachRun.at(e).Add(&Hist_Target_SR_MSCW.at(e).at(s));
                    Hist_Target_Excess_EachRun.at(e).Add(&Hist_Target_BkgSR_MSCW.at(e).at(s),-1.);
                }
                norm_bin_low = Hist_Target_Excess_EachRun.at(e).FindBin(MSCW_cut_lower);
                norm_bin_up = Hist_Target_Excess_EachRun.at(e).FindBin(MSCW_cut_blind);
                double excess_this_run = Hist_Target_Excess_EachRun.at(e).Integral(norm_bin_low,norm_bin_up);
                double effarea_time = Hist_EffAreaTime.GetBinContent(Hist_EffAreaTime.FindBin(energy_bins[e]));
                double flux_this_run = excess_this_run*1000./(effarea_time*(energy_bins[e+1]-energy_bins[e])*1e-12*10000.);
                std::cout << "excess_this_run = " << excess_this_run << std::endl;
                std::cout << "flux_this_run = " << flux_this_run << std::endl;
                Hist_Measured_Electron_Flux.Fill(energy_bins[e],flux_this_run,effarea_time);
                
                for (int s=0;s<Number_of_SR;s++)
                {
                    Hist_Target_BkgSR_MSCW_SumRuns.at(e).at(s).Add(&Hist_Target_BkgSR_MSCW.at(e).at(s)); // this is summing over sublist, errors are not correlated.
                    Hist_Scaled_GammaMC_SR_MSCW_SumRuns.at(e).at(s).Add(&Hist_Scaled_GammaMC_SR_MSCW.at(e).at(s));
                }
                for (int run=0;run<Sublist.at(subrun).size();run++)
                {
                    std::cout << "run " << int(Sublist.at(subrun)[run]) << std::endl;
                }
                for (int s=1;s<Number_of_CR;s++)
                {
                    int norm_bin_low = Hist_Target_BkgCR_MSCW_SumRuns.at(e).at(s).FindBin(MSCW_cut_lower);
                    int norm_bin_up = Hist_Target_BkgCR_MSCW_SumRuns.at(e).at(s).FindBin(MSCW_cut_blind);
                    std::cout << "Hist_Target_BkgCR_MSCW_SumRuns.at(e).at("<< s<< ").Integral() = " << Hist_Target_BkgCR_MSCW_SumRuns.at(e).at(s).Integral(norm_bin_low,norm_bin_up) << std::endl;
                    std::cout << "Hist_Target_CR_MSCW_SumRuns.at(e).at("<< s<< ").Integral() = " << Hist_Target_CR_MSCW_SumRuns.at(e).at(s).Integral(norm_bin_low,norm_bin_up) << std::endl;
                }
                for (int s=0;s<Number_of_SR;s++)
                {
                    int norm_bin_low = Hist_Target_BkgSR_MSCW_SumRuns.at(e).at(s).FindBin(MSCW_cut_lower);
                    int norm_bin_up = Hist_Target_BkgSR_MSCW_SumRuns.at(e).at(s).FindBin(MSCW_cut_blind);
                    std::cout << "Hist_Target_BkgSR_MSCW_SumRuns.at(e).at("<< s<< ").Integral() = " << Hist_Target_BkgSR_MSCW_SumRuns.at(e).at(s).Integral(norm_bin_low,norm_bin_up) << std::endl;
                    std::cout << "Hist_Target_SR_MSCW_SumRuns.at(e).at("<< s<< ").Integral() = " << Hist_Target_SR_MSCW_SumRuns.at(e).at(s).Integral(norm_bin_low,norm_bin_up) << std::endl;
                }
                std::cout << "electron_count[e] = " << electron_count[e] << std::endl;
            }
        }

        for (int e=0;e<N_energy_bins;e++)
        {
            Hist_Target_BkgSR_MSCW_SumRuns_SumSRs.at(e).Reset();
            for (int s=0;s<Number_of_SR;s++)
            {
                std::cout << "=================================================" << std::endl;
                std::cout << "Energy "  << energy_bins[e] << std::endl;
                std::cout << "SR "  << s << std::endl;
                std::cout << "dark_kernel_rms = " << dark_kernel_rms.at(e).at(s) << std::endl;
                std::cout << "target_kernel_rms = " << target_kernel_rms.at(e).at(s) << std::endl;
                std::cout << "dark_kernel_shift = " << dark_kernel_shift.at(e).at(s) << std::endl;
                std::cout << "target_kernel_shift = " << target_kernel_shift.at(e).at(s) << std::endl;
                Hist_Target_BkgSR_MSCW_SumRuns_SumSRs.at(e).Add(&Hist_Target_BkgSR_MSCW_SumRuns.at(e).at(s));
            }
        }

        std::cout << "=================================================" << std::endl;
        for (int e=0;e<=N_energy_bins;e++) {
            cosmic_electron.push_back(electron_count[e]);
            std::cout << "// Energy "  << energy_bins[e] << std::endl;
            std::cout << "electron_flux[" << e << "] = " << Hist_Measured_Electron_Flux.GetBinContent(Hist_Measured_Electron_Flux.FindBin(energy_bins[e]))  << ";" << std::endl;;
            std::cout << "electron_flux_err[" << e << "] = " << Hist_Measured_Electron_Flux.GetBinError(Hist_Measured_Electron_Flux.FindBin(energy_bins[e]))  << ";" << std::endl;;
        }
        for (int e=0;e<=N_energy_bins;e++) {
            std::cout << "if (energy>=" << energy_bins[e] << "&&energy<" << energy_bins[e+1] << ")";
            std::cout << " return std::make_pair(" << Hist_Target_EndPoint_0.GetBinContent(e+1) << "," << Hist_Target_EndPoint_1.GetBinContent(e+1) << ");" << std::endl;
        }

        
        int norm_bin_low_ring = Hist_Target_ON_MSCW_Alpha.FindBin(Norm_Lower);
        int norm_bin_up_ring = Hist_Target_ON_MSCW_Alpha.FindBin(Norm_Upper);
        double scale_ring = Hist_Target_ON_MSCW_Alpha.Integral(norm_bin_low_ring,norm_bin_up_ring)/Hist_Target_OFF_MSCW_Alpha.Integral(norm_bin_low_ring,norm_bin_up_ring);

        int Number_of_CR_new = Number_of_CR;
        int Number_of_SR_new = Number_of_SR;
        TString ConvergeOrNot = "";
        if (!DoConverge) ConvergeOrNot = "_NoConverge";
        TFile OutputFile("output_May25/Deconvolution_"+TString(target)+"_Ntel"+std::to_string(NTelMin)+"to"+std::to_string(NTelMax)+"_Elev"+std::to_string(int(Target_Elev_cut_lower))+"to"+std::to_string(int(Target_Elev_cut_upper))+"_Azim"+std::to_string(int(Target_Azim_cut_lower))+"to"+std::to_string(int(Target_Azim_cut_upper))+"_Theta2"+std::to_string(int(10.*Theta2_cut_lower))+"to"+std::to_string(int(10.*Theta2_cut_upper))+"_MSCWCut"+std::to_string(int(10.*MSCW_cut_upper))+"_MSCWBlind"+std::to_string(int(10.*MSCW_cut_blind))+ConvergeOrNot+".root",WriteType);
        if (WriteType=="recreate")
        {
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
          InfoTree.Branch("cosmic_electron","std::vector<double>",&cosmic_electron);
          InfoTree.Branch("exposure_hours",&exposure_hours,"exposure_hours/D");
          InfoTree.Fill();
          InfoTree.Write();
          Hist_Target_TelElevAzim.Write();
          Hist_Target_TelRaDec.Write();
          Hist_Target_TelRaDec_AfterCut.Write();
          Hist_Measured_Electron_Flux.Write();
          Hist_Target_EndPoint_0.Write();
          Hist_Target_EndPoint_1.Write();
        }
        for (int e=0;e<N_energy_bins;e++) {
                if (e!=run_energy_bin) continue;
                Hist_Target_SR_MSCL.at(e).Write();
                Hist_Target_Amplitude.at(e).Write();
                Hist_Target_Mean.at(e).Write();
                Hist_Target_RMS.at(e).Write();
                Hist_Target_Deconv_MSCW.at(e).Write();
                Hist_Target_MSCLW.at(e).Write();
                Hist_Target_SR_MSCW_SumRuns_SumSRs.at(e).Write();
                Hist_ElectronMC_SR_MSCW_SumRuns_SumSRs.at(e).Write();
                Hist_GammaMC_SR_MSCW_SumRuns_SumSRs.at(e).Write();
                Hist_Target_BkgSR_MSCW_SumRuns_SumSRs.at(e).Write();
                for (int s=0;s<Number_of_SR;s++)
                {
                    Hist_Scaled_ElectronMC_SR_MSCW_SumRuns.at(e).at(s).Write();
                    Hist_Scaled_GammaMC_SR_MSCW_SumRuns.at(e).at(s).Write();
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
                    Hist_Scaled_ElectronMC_CR_MSCW_SumRuns.at(e).at(s).Write();
                    Hist_Scaled_GammaMC_CR_MSCW_SumRuns.at(e).at(s).Write();
                    Hist_Target_CR_MSCW_SumRuns.at(e).at(s).Write();
                    Hist_Target_BkgCR_MSCW_SumRuns.at(e).at(s).Write();
                }
        }
        OutputFile.Close();

}
