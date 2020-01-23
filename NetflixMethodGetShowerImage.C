//
// 
// This code reads output from EventDisplay and creates a root file which stores histograms of MSCW and MSCL.
// What histograms are created by this code:
//      1. 2D histograms of data events in MSCL and MSCW parameter space (as function of energy).
//      2. 2D histograms of MC electron events (normalized to electron flux) in MSCL and MSCW parameter space.
//      3. 1D theta2 histograms of data events (gamma-ray region and CR region).
//      4. 2D skymap histograms of data events (gamma-ray region and CR region).
// Author: Ruo Yu Shang
//
//
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
#include <math.h> 

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
#include "TChain.h"
#include "TBranch.h"

#include "/home/rshang/EventDisplay/EVNDISP-480e/inc/VEvndispRunParameter.h"

#include "GetRunList.h"

#include <complex>
#include "/home/rshang/Eigen/eigen-eigen-323c052e1731/Eigen/Dense"
#include "/home/rshang/Eigen/eigen-eigen-323c052e1731/Eigen/StdVector"
using namespace Eigen;
//using Eigen::MatrixXd;

double MSCW_cut_lower = -1.0;
double MSCW_cut_blind = 1.0;
double MSCW_cut_upper = 1.0;
double MSCL_cut_lower = -1.0;
double MSCL_cut_blind = 1.0;
double MSCL_cut_upper = 1.0;

const int N_energy_bins = 1;
double energy_bins[N_energy_bins+1] = {100,1e4};
int N_bins_for_deconv_at_E[N_energy_bins] = {30};
//const int N_energy_bins = 11;
//double energy_bins[N_energy_bins+1] = {200,237,282,335,398,473,562,794,1122,2239,4467,8913};
//int N_bins_for_deconv_at_E[N_energy_bins] = {40,40,40,40,40,40,40,40,40,40,40};
//const int N_energy_bins = 5;
//double energy_bins[N_energy_bins+1] = {pow(10,2.0),pow(10,2.2),pow(10,2.4),pow(10,2.6),pow(10,2.8),pow(10,4.0)};
//int N_bins_for_deconv_at_E[N_energy_bins] = {40,40,40,40,40};
//const int N_energy_bins = 12;
//double energy_bins[N_energy_bins+1] = {pow(10,2.0),pow(10,2.1),pow(10,2.2),pow(10,2.3),pow(10,2.4),pow(10,2.5),pow(10,2.6),pow(10,2.7),pow(10,2.8),pow(10,3.0),pow(10,3.2),pow(10,3.6),pow(10,4.0)};
//int N_bins_for_deconv_at_E[N_energy_bins] = {40,40,40,40,40,40,40,40,40,40,40,40};

const int N_energy_fine_bins = 20;
double energy_fine_bins[N_energy_fine_bins+1] = {pow(10,2.0),pow(10,2.1),pow(10,2.2),pow(10,2.3),pow(10,2.4),pow(10,2.5),pow(10,2.6),pow(10,2.7),pow(10,2.8),pow(10,2.9),pow(10,3.0),pow(10,3.1),pow(10,3.2),pow(10,3.3),pow(10,3.4),pow(10,3.5),pow(10,3.6),pow(10,3.7),pow(10,3.8),pow(10,3.9),pow(10,4.0)};
double gamma_flux[N_energy_fine_bins] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
double gamma_count[N_energy_fine_bins] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
double raw_gamma_count[N_energy_fine_bins] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

int N_bins_for_deconv = 30;
double MSCW_plot_lower = -1.;
double MSCW_plot_upper = 2.;
double MSCL_plot_lower = -1.;
double MSCL_plot_upper = 2.;
double Theta2_cut_lower = 0;
double Theta2_cut_upper = 0;
double Theta2_upper_limit = 10;
double TelElev_upper = 90.;
double TelElev_lower = 60.;
char target[50] = "";

double TelElevation = 0;
double TelAzimuth = 0;
double TelRAJ2000 = 0;
double TelDecJ2000 = 0;
double ErecS = 0;
double EChi2S = 0;
double MSCW = 0;
int NImages = 0;
double Xcore = 0.;
double Ycore = 0.;
double SizeSecondMax = 0;
double MSCL = 0;
double Time = 0;
double Shower_Ze = 0;
double Shower_Az = 0;
double SlantDepth = 0;
float EmissionHeight = 0;
float EmissionHeightChi2 = 0;
double Xoff = 0;
double Yoff = 0;
double Xoff_derot = 0;
double Yoff_derot = 0;
double R2off = 0;
double Phioff = 0;
double theta2 = 0;
double ra_sky = 0;
double dec_sky = 0;
double exposure_hours = 0.;
double exposure_hours_dark = 0.;
double NSB_avg = 0.;
double NSB_avg_dark = 0.;
double NSB_diff_dark = 0.;
double Elev_diff_dark = 0.;
double NSB_low_cut = 0.;
double NSB_high_cut = 100.;
double mean_tele_point_ra = 0.;
double mean_tele_point_dec = 0.;
double mean_tele_point_l = 0.;
double mean_tele_point_b = 0.;
double roi_ra = 0.;
double roi_dec = 0.;
double roi_radius = 0.;
double DarkRun_theta2_lower = 0.;
double DarkRun_theta2_upper = 10.;

vector<pair<double,double>> ON_pointing_radec_for_DarkRun;

pair<double,double> GetSourceRaDec(TString source_name)
{
    double Source_RA = 0.;
    double Source_Dec = 0.;
    if (source_name=="SgrAV6")
    {
            Source_RA = 266.415;
                Source_Dec = -29.006;
    }
    if (source_name=="GemingaV6")
    {
            Source_RA = 98.117;
                Source_Dec = 17.367;
    }
    if (source_name=="GemingaV5")
    {
            Source_RA = 98.117;
                Source_Dec = 17.367;
    }
    if (source_name=="PKS1441V6")
    {
            Source_RA = 220.987;
                Source_Dec = 25.029;
    }
    if (source_name=="MS1221V6")
    {
            Source_RA = 186.101;
                Source_Dec = 24.607;
    }
    if (source_name=="S3_1227_V6")
    {
            Source_RA = 187.559;
                Source_Dec = 25.302;
    }
    if (source_name=="Crab")
    {
            Source_RA = 83.633;
                Source_Dec = 22.014;
    }
    if (source_name=="CrabV5")
    {
            Source_RA = 83.633;
                Source_Dec = 22.014;
    }
    if (source_name=="Mrk421")
    {
            Source_RA = 166.079;
                Source_Dec = 38.195;
    }
    if (source_name=="H1426")
    {
            Source_RA = 217.136;
                Source_Dec = 42.673;
    }
    if (source_name=="1ES0229")
    {
            Source_RA = 38.222;
                Source_Dec = 20.273;
    }
    if (source_name=="PKS1424")
    {
            Source_RA = 216.750;
                Source_Dec = 23.783;
    }
    if (source_name=="3C264")
    {
            Source_RA = 176.271;
                Source_Dec = 19.606;
    }
    if (source_name=="OJ287V6")
    {
            Source_RA = 133.705;
                Source_Dec = 20.100;
    }
    if (source_name=="RBS0413V6")
    {
            Source_RA = 49.946;
                Source_Dec = 18.762;
    }
    if (source_name=="PG1553V6")
    {
            Source_RA = 238.936;
                Source_Dec = 11.195;
    }
    if (source_name=="Segue1V6")
    {
            Source_RA = 151.767;
                Source_Dec = 16.082;
    }
    if (source_name=="1ES1011V6")
    {
            Source_RA = 153.767;
                Source_Dec = 49.434;
    }
    if (source_name=="NGC1275V6")
    {
            Source_RA = 49.950;
                Source_Dec = 41.512;
    }
    if (source_name=="1ES0647V6")
    {
            Source_RA = 102.694;
                Source_Dec = 25.050;
    }
    if (source_name=="1ES1440V6")
    {
            Source_RA = 220.701;
                Source_Dec = 12.011;
    }
    if (source_name=="1ES1741V6")
    {
            Source_RA = 266.005;
                Source_Dec = 19.546;
    }
    if (source_name=="IC443HotSpot")
    {
            Source_RA = 94.511;
                Source_Dec = 22.660;
    }
    if (source_name=="RGBJ0710")
    {
            Source_RA = 107.610;
                Source_Dec = 59.150;
    }
    if (source_name=="CasA")
    {
            Source_RA = 350.808;
                Source_Dec = 58.807;
    }
    if (source_name=="M82")
    {
            Source_RA = 148.970;
                Source_Dec = 69.679;
    }
    if (source_name=="G079")
    {
            Source_RA = 308.119;
                Source_Dec = 40.328;
    }
    if (source_name=="WComaeV6")
    {
            Source_RA = 185.382;
                Source_Dec = 28.233;
    }
    if (source_name=="1ES1218V6")
    {
            Source_RA = 185.360;
                Source_Dec = 30.191;
    }
    if (source_name=="MGRO_J1908_V6")
    {
            Source_RA = 286.975;
                Source_Dec = 6.269;
    }
    if (source_name=="MGRO_J1908_V5")
    {
            Source_RA = 286.975;
                Source_Dec = 6.269;
    }
    if (source_name=="Segue1V5")
    {
            Source_RA = 151.767;
                Source_Dec = 16.082;
    }
    if (source_name=="IC443HotSpotV5")
    {
            Source_RA = 94.511;
                Source_Dec = 22.660;
    }
    if (source_name=="Proton")
    {
            Source_RA = 0.;
                Source_Dec = 0.;
    }
    if (source_name=="Proton_NSB750")
    {
            Source_RA = 0.;
                Source_Dec = 0.;
    }
    if (source_name=="Everything")
    {
            Source_RA = 0.;
                Source_Dec = 0.;
    }
    return std::make_pair(Source_RA,Source_Dec);
}
pair<double,double> ConvertRaDecToGalactic(double Ra, double Dec)
{
    double delta = Dec*M_PI/180.;
    double delta_G = 27.12825*M_PI/180.;
    double alpha = Ra*M_PI/180.;
    double alpha_G = 192.85948*M_PI/180.;
    double l_NCP = 122.93192*M_PI/180.;
    double sin_b = sin(delta)*sin(delta_G)+cos(delta)*cos(delta_G)*cos(alpha-alpha_G);
    double cos_b = cos(asin(sin_b));
    double sin_l_NCP_m_l = cos(delta)*sin(alpha-alpha_G)/cos_b;
    double b = (asin(sin_b))*180./M_PI;
    double l = (l_NCP-asin(sin_l_NCP_m_l))*180./M_PI;
    double b_round = floor(b*pow(10,2))/pow(10,2);
    double l_round = floor(l*pow(10,2))/pow(10,2);
    return std::make_pair(l_round,b_round);
}
pair<double,double> GetRunRaDec(string file_name, int run)
{
    char run_number[50];
    sprintf(run_number, "%i", int(run));
    TFile*  input_file = TFile::Open(file_name.c_str());
    TTree* pointing_tree = nullptr;
    pointing_tree = (TTree*) input_file->Get("run_"+TString(run_number)+"/stereo/pointingDataReduced");
    pointing_tree->SetBranchAddress("TelElevation",&TelElevation);
    pointing_tree->SetBranchAddress("TelAzimuth",&TelAzimuth);
    pointing_tree->SetBranchAddress("TelRAJ2000",&TelRAJ2000);
    pointing_tree->SetBranchAddress("TelDecJ2000",&TelDecJ2000);
    double total_entries = (double)pointing_tree->GetEntries();
    pointing_tree->GetEntry(int(total_entries/2.));
    double TelRAJ2000_tmp = TelRAJ2000*180./M_PI;
    double TelDecJ2000_tmp = TelDecJ2000*180./M_PI;
    input_file->Close();
    return std::make_pair(TelRAJ2000_tmp,TelDecJ2000_tmp);
}
pair<double,double> GetRunElevAzim(string file_name, int run)
{
    char run_number[50];
    sprintf(run_number, "%i", int(run));
    TFile*  input_file = TFile::Open(file_name.c_str());
    TTree* pointing_tree = nullptr;
    pointing_tree = (TTree*) input_file->Get("run_"+TString(run_number)+"/stereo/pointingDataReduced");
    pointing_tree->SetBranchAddress("TelElevation",&TelElevation);
    pointing_tree->SetBranchAddress("TelAzimuth",&TelAzimuth);
    pointing_tree->SetBranchAddress("TelRAJ2000",&TelRAJ2000);
    pointing_tree->SetBranchAddress("TelDecJ2000",&TelDecJ2000);
    double total_entries = (double)pointing_tree->GetEntries();
    pointing_tree->GetEntry(int(total_entries/2.));
    double TelElevation_tmp = TelElevation;
    double TelAzimuth_tmp = TelAzimuth;
    input_file->Close();
    return std::make_pair(TelElevation_tmp,TelAzimuth_tmp);
}
bool PointingSelection(string file_name,int run, double Elev_cut_lower, double Elev_cut_upper, double Azim_cut_lower, double Azim_cut_upper)
{
    if (run>100000) return true;
    char run_number[50];
    sprintf(run_number, "%i", int(run));
    TFile*  input_file = TFile::Open(file_name.c_str());

    //VEvndispRunParameter* fPar = 0;
    //fPar = ( VEvndispRunParameter* )input_file->Get( "run_"+TString(run_number)+"/stereo/runparameterV2" );
    //std::cout << "NSB scale = " << fPar->fNSBscale << std::endl;

    TTree* pointing_tree = nullptr;
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
    input_file->Close();
    return true;
}

double GetRunNSB(int run_number)
{
    string line;
    char delimiter = ',';
    string acc_runnumber = "";
    string acc_nsb = "";
    int nth_line = 0;
    int nth_delimiter = 0;
    std::string::size_type sz;
    double NSB = 0.;

    ifstream myfile ("/home/rshang/EventDisplay/NewBkgMethodExtendedSource/allrunsdatabaseinfo.txt");
    if (myfile.is_open())
    {
        while ( getline(myfile,line) )
        {
            acc_runnumber = "";
            acc_nsb = "";
            nth_delimiter = 0;
            for(int i = 0; i < line.size(); i++)
            {
                if(line[i] == delimiter && nth_line!=0)
                {
                    //if (nth_delimiter==1 && std::stoi(acc_runnumber,nullptr,10)==run_number) 
                    //{
                    //    cout << "run = " << acc_runnumber << '\n';
                    //}
                    if (nth_delimiter==17 && std::stoi(acc_runnumber,nullptr,10)==run_number) 
                    {
                        //cout << "NSB = " << std::stod(acc_nsb,&sz) << '\n';
                        NSB = std::stod(acc_nsb,&sz);
                    }
                    nth_delimiter += 1;
                }
                else if (nth_delimiter==1)
                {
                    acc_runnumber += line[i];
                }
                else if (nth_delimiter==17)
                {
                    acc_nsb += line[i];
                }
            }
            nth_line += 1;
        }
        myfile.close();
    }
    else cout << "Unable to open file"; 

    return NSB;
}

vector<pair<string,int>> SelectONRunList(vector<pair<string,int>> Data_runlist, double Elev_cut_lower, double Elev_cut_upper, double Azim_cut_lower, double Azim_cut_upper)
{
    vector<pair<string,int>> new_list;
    for (int run=0;run<Data_runlist.size();run++)
    {
        char run_number[50];
        char Data_observation[50];
        sprintf(run_number, "%i", int(Data_runlist[run].second));
        sprintf(Data_observation, "%s", Data_runlist[run].first.c_str());
        double NSB_thisrun = GetRunNSB(int(Data_runlist[run].second));
        if (TString(Data_observation).Contains("NSB075")) NSB_thisrun = 3.30;
        if (TString(Data_observation).Contains("NSB100")) NSB_thisrun = 3.84;
        if (TString(Data_observation).Contains("NSB150")) NSB_thisrun = 4.69;
        if (TString(Data_observation).Contains("NSB200")) NSB_thisrun = 5.41;
        if (TString(Data_observation).Contains("NSB250")) NSB_thisrun = 6.00;
        if (TString(Data_observation).Contains("NSB325")) NSB_thisrun = 6.83;
        if (TString(Data_observation).Contains("NSB425")) NSB_thisrun = 7.73;
        if (TString(Data_observation).Contains("NSB550")) NSB_thisrun = 8.75;
        if (TString(Data_observation).Contains("NSB750")) NSB_thisrun = 10.18;
        if (TString(Data_observation).Contains("NSB1000")) NSB_thisrun = 11.77;
        if (NSB_thisrun < NSB_low_cut) continue;
        if (NSB_thisrun > NSB_high_cut) continue;
        string filename;
        filename = TString("$VERITAS_USER_DATA_DIR/"+TString(Data_observation)+"_V6_Moderate-TMVA-BDT.RB."+TString(run_number)+".root");

        if (!PointingSelection(filename,int(Data_runlist[run].second),Elev_cut_lower,Elev_cut_upper,Azim_cut_lower,Azim_cut_upper)) continue;
        new_list.push_back(std::make_pair(Data_runlist[run].first,Data_runlist[run].second));

    }
    return new_list;
}

vector<pair<string,int>> SelectOFFRunList(vector<pair<string,int>> ON_runlist, vector<pair<string,int>> OFF_runlist)
{

    vector<pair<double,double>> ON_pointing;
    vector<pair<double,double>> ON_pointing_radec;
    vector<double> ON_NSB;
    for (int on_run=0;on_run<ON_runlist.size();on_run++)
    {
        char ON_runnumber[50];
        char ON_observation[50];
        sprintf(ON_runnumber, "%i", int(ON_runlist[on_run].second));
        sprintf(ON_observation, "%s", ON_runlist[on_run].first.c_str());
        string ON_filename;
        ON_filename = TString("$VERITAS_USER_DATA_DIR/"+TString(ON_observation)+"_V6_Moderate-TMVA-BDT.RB."+TString(ON_runnumber)+".root");
        if (TString(ON_observation).Contains("Proton")) ON_pointing.push_back(std::make_pair(70,0));
        else ON_pointing.push_back(GetRunElevAzim(ON_filename,int(ON_runlist[on_run].second)));
        if (TString(ON_observation).Contains("Proton")) ON_pointing_radec.push_back(std::make_pair(0,0));
        else ON_pointing_radec.push_back(GetRunRaDec(ON_filename,int(ON_runlist[on_run].second)));
        double NSB_thisrun = GetRunNSB(int(ON_runlist[on_run].second));
        if (TString(ON_observation).Contains("NSB075")) NSB_thisrun = 3.30;
        if (TString(ON_observation).Contains("NSB100")) NSB_thisrun = 3.84;
        if (TString(ON_observation).Contains("NSB150")) NSB_thisrun = 4.69;
        if (TString(ON_observation).Contains("NSB200")) NSB_thisrun = 5.41;
        if (TString(ON_observation).Contains("NSB250")) NSB_thisrun = 6.00;
        if (TString(ON_observation).Contains("NSB325")) NSB_thisrun = 6.83;
        if (TString(ON_observation).Contains("NSB425")) NSB_thisrun = 7.73;
        if (TString(ON_observation).Contains("NSB550")) NSB_thisrun = 8.75;
        if (TString(ON_observation).Contains("NSB750")) NSB_thisrun = 10.18;
        if (TString(ON_observation).Contains("NSB1000")) NSB_thisrun = 11.77;
        ON_NSB.push_back(NSB_thisrun);
    }

    vector<pair<double,double>> OFF_pointing;
    vector<double> OFF_NSB;
    for (int off_run=0;off_run<OFF_runlist.size();off_run++)
    {
        char OFF_runnumber[50];
        char OFF_observation[50];
        sprintf(OFF_runnumber, "%i", int(OFF_runlist[off_run].second));
        sprintf(OFF_observation, "%s", OFF_runlist[off_run].first.c_str());
        string OFF_filename;
        OFF_filename = TString("$VERITAS_USER_DATA_DIR/"+TString(OFF_observation)+"_V6_Moderate-TMVA-BDT.RB."+TString(OFF_runnumber)+".root");
        if (TString(OFF_observation).Contains("Proton")) OFF_pointing.push_back(std::make_pair(70,0));
        else OFF_pointing.push_back(GetRunElevAzim(OFF_filename,int(OFF_runlist[off_run].second)));
        double NSB_thisrun = GetRunNSB(int(OFF_runlist[off_run].second));
        if (TString(OFF_observation).Contains("NSB075")) NSB_thisrun = 3.30;
        if (TString(OFF_observation).Contains("NSB100")) NSB_thisrun = 3.84;
        if (TString(OFF_observation).Contains("NSB150")) NSB_thisrun = 4.69;
        if (TString(OFF_observation).Contains("NSB200")) NSB_thisrun = 5.41;
        if (TString(OFF_observation).Contains("NSB250")) NSB_thisrun = 6.00;
        if (TString(OFF_observation).Contains("NSB325")) NSB_thisrun = 6.83;
        if (TString(OFF_observation).Contains("NSB425")) NSB_thisrun = 7.73;
        if (TString(OFF_observation).Contains("NSB550")) NSB_thisrun = 8.75;
        if (TString(OFF_observation).Contains("NSB750")) NSB_thisrun = 10.18;
        if (TString(OFF_observation).Contains("NSB1000")) NSB_thisrun = 11.77;
        OFF_NSB.push_back(NSB_thisrun);
    }

    vector<pair<string,int>> new_list;
    vector<pair<double,double>> ON_pointing_radec_new;
    for (int n=0;n<1;n++)
    {
        for (int on_run=0;on_run<ON_runlist.size();on_run++)
        {

            pair<string,int> best_match;
            pair<double,double> best_pointing;
            double best_chi2 = 10000.;
            for (int off_run=0;off_run<OFF_runlist.size();off_run++)
            {
                if (ON_runlist[on_run].first.find("Proton")==std::string::npos)
                {
                    if (ON_runlist[on_run].first.compare(OFF_runlist[off_run].first) == 0) continue;
                }
                bool already_used_run = false;
                for (int new_run=0;new_run<new_list.size();new_run++)
                {
                    if (int(new_list[new_run].second)==int(OFF_runlist[off_run].second)) already_used_run = true;
                }
                if (already_used_run) continue;
                double chi2 = pow(ON_pointing[on_run].first+Elev_diff_dark-OFF_pointing[off_run].first,2);
                chi2 += 10.*pow(ON_NSB[on_run]+NSB_diff_dark-OFF_NSB[off_run],2);
                if (best_chi2>chi2)
                {
                    best_chi2 = chi2;
                    best_match = OFF_runlist[off_run];
                    best_pointing = OFF_pointing[off_run];
                }
            }
            if (best_chi2<10000.) 
            {
                new_list.push_back(best_match);
                ON_pointing_radec_new.push_back(ON_pointing_radec[on_run]);
                std::cout << "add run:" << std::endl;
                std::cout << best_match.first << " " << best_match.second << std::endl;
                std::cout << best_pointing.first << " " << best_pointing.second << std::endl;
            }
        }
    }
    ON_pointing_radec_for_DarkRun = ON_pointing_radec_new;
    return new_list;
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
double FitPowerLawFunction(Double_t *x, Double_t *par) {
    double xx =x[0];
    return par[0]*pow(xx,par[1]);
}
double GetCrabFlux(double energy_gev)
{
    //double flux = 1e-7*3.75*pow(energy_gev/1000.,-2.2);
    double flux = 3.75*pow(10,-7)*pow(energy_gev/1000.,-2.467-0.16*log(energy_gev/1000.));
    return flux;
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
    return std::make_pair(1e-12*10000.*func->Eval(energy),1e-12*10000.*func->Eval(energy));
    //
    //TGraph *func = new TGraph(&Hist_Flux);
    //return std::make_pair(1.*func->Eval(energy),1.*func->Eval(energy));
}
bool GammaFoV() {
    double x = ra_sky-mean_tele_point_ra;
    double y = dec_sky-(mean_tele_point_dec-0.3);
    double size1 = 0.5;
    double diff1 = pow(x*x+y*y-size1,3)-2.*(x*x*y*y*y);
    double size2 = 0.9;
    double diff2 = pow(x*x+y*y-size2,3)-2.*(x*x*y*y*y);
    if (diff1<0) return false;
    if (diff2>0) return false;
    return true;
}
bool DarkFoV() {
    //if (R2off<Theta2_cut_lower) return false;
    //if (R2off>Theta2_cut_upper) return false;
    if (theta2<0.3) return false;
    //if (theta2>10.) return false;
    //if (R2off<DarkRun_theta2_lower) return false;
    //if (R2off>DarkRun_theta2_upper) return false;
    return true;
}
bool FoV() {
    //if (R2off<Theta2_cut_lower) return false;
    //if (R2off>Theta2_cut_upper) return false;
    if (theta2<Theta2_cut_lower) return false;
    //if (theta2>Theta2_cut_upper) return false;
    return true;
}
bool RoIFoV() {
    double x = ra_sky-roi_ra;
    double y = dec_sky-roi_dec;
    double radius = pow(x*x+y*y,0.5);
    if (radius>roi_radius) return false;
    return true;
}
bool RingFoV() {
    if (theta2<Theta2_cut_upper) return false;
    if (theta2>Theta2_cut_upper*5) return false;
    return true;
}
bool SelectNImages(int Nmin, int Nmax)
{
    if (NImages<Nmin) return false;
    if (NImages>Nmax) return false;
    return true;
}
bool SignalSelectionTheta2()
{
    if (MSCW>MSCW_cut_blind) return false;
    if (MSCW<MSCW_cut_lower) return false;
    if (MSCL>MSCL_cut_blind) return false;
    if (MSCL<MSCL_cut_lower) return false;
    return true;
}
bool ControlSelectionTheta2()
{
    if (SignalSelectionTheta2()) return false;
    if (MSCL>MSCL_cut_blind+1.0) return false;
    if (MSCW>MSCW_cut_blind+1.0) return false;
    //if (MSCL<1.0 && MSCW<1.0) return false;
    //if (MSCL>3.0) return false;
    //if (MSCW>2.0) return false;
    return true;
}
void NetflixMethodGetShowerImage(string target_data, double PercentCrab, double tel_elev_lower_input, double tel_elev_upper_input, double Elev_diff, double NSB_diff, bool isON, double MSCW_cut_input, double MSCL_cut_input, double Theta2_cut_input)
{

    TH1::SetDefaultSumw2();
    sprintf(target, "%s", target_data.c_str());
    MSCW_cut_blind = MSCW_cut_input;
    MSCL_cut_blind = MSCL_cut_input;
    TelElev_lower = tel_elev_lower_input;
    TelElev_upper = tel_elev_upper_input;
    Theta2_cut_lower = 0.3;
    //Theta2_cut_upper = Theta2_cut_input;
    Theta2_cut_upper = 1e4;
    if (isON) Theta2_cut_lower = 0.;
    else Theta2_cut_upper = max(Theta2_cut_input,1.0);
    TString file_tag;
    if (isON) file_tag = "ON";
    else file_tag = "OFF";
    Elev_diff_dark = Elev_diff;
    NSB_diff_dark = NSB_diff;

    vector<pair<string,int>> PhotonMC_runlist = GetRunList("Photon");
    vector<pair<string,int>> PhotonData_runlist = GetRunList("Crab");
    if (TString(target).Contains("Mrk421")) PhotonData_runlist = GetRunList("Crab");
    if (TString(target).Contains("Crab")) PhotonData_runlist = GetRunList("Mrk421");
    if (TString(target).Contains("V5")) PhotonData_runlist = GetRunList("CrabV5");
    PhotonData_runlist = SelectONRunList(PhotonData_runlist,TelElev_lower,TelElev_upper,0,360);

    // Get a list of target observation runs
    vector<pair<string,int>> Data_runlist_init = GetRunList(target);
    vector<pair<string,int>> Data_runlist;
    if (!TString(target).Contains("Proton")) Data_runlist = SelectONRunList(Data_runlist_init,TelElev_lower,TelElev_upper,0,360);
    else Data_runlist = Data_runlist_init;
    std::cout << "Data_runlist size = " << Data_runlist.size() << std::endl;
    if (Data_runlist.size()==0) return;

    // Get a list of dark observation runs
    vector<pair<string,int>> Dark_runlist_init = GetRunList("Everything");
    if (TString(target).Contains("V5")) Dark_runlist_init = GetRunList("EverythingV5");
    if (TString(target).Contains("Proton")) Dark_runlist_init = GetRunList("EverythingProton");
    vector<pair<string,int>> Dark_runlist;
    std::cout << "initial Dark_runlist size = " << Dark_runlist_init.size() << std::endl;
    Dark_runlist = SelectOFFRunList(Data_runlist, Dark_runlist_init);
    std::cout << "final Dark_runlist size = " << Dark_runlist.size() << std::endl;
    //vector<pair<string,int>> Dark_runlist = GetRunList("Proton");

    mean_tele_point_ra = 0.;
    mean_tele_point_dec = 0.;
    pair<double,double> source_ra_dec = GetSourceRaDec(TString(target));
    mean_tele_point_ra = source_ra_dec.first;
    mean_tele_point_dec = source_ra_dec.second;

    roi_ra = mean_tele_point_ra;
    roi_dec = mean_tele_point_dec;
    roi_radius = 0.5;
    if (TString(target).Contains("Geminga")) 
    {
        roi_radius = 1.0;
    }
    if (TString(target).Contains("Segue1")) 
    {
        roi_ra = 151.8;
        roi_dec = 16.7;
        roi_radius = 0.5;
    }

    TH1D Hist_ErecS = TH1D("Hist_ErecS","",N_energy_bins,energy_bins);
    TH1D Hist_ErecS_fine = TH1D("Hist_ErecS_fine","",N_energy_fine_bins,energy_fine_bins);
    TH1D Hist_EffArea = TH1D("Hist_EffArea","",N_energy_fine_bins,energy_fine_bins);
    TH1D Hist_Dark_NSB = TH1D("Hist_Dark_NSB","",20,4,14);
    TH1D Hist_Data_NSB = TH1D("Hist_Data_NSB","",20,4,14);
    TH2D Hist_Dark_ShowerDirection = TH2D("Hist_Dark_ShowerDirection","",180,0,360,90,0,90);
    TH2D Hist_Data_ShowerDirection = TH2D("Hist_Data_ShowerDirection","",180,0,360,90,0,90);
    TH2D Hist_Data_MSCLW_incl = TH2D("Hist_Data_MSCLW_incl","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper);
    vector<TH2D> Hist_GammaDark_MSCLW;
    vector<TH2D> Hist_GammaMC_MSCLW;
    vector<TH2D> Hist_GammaData_MSCLW;
    vector<TH2D> Hist_GammaDataON_MSCLW;
    vector<TH2D> Hist_GammaDataOFF_MSCLW;
    vector<TH2D> Hist_Data_MSCLW;
    vector<TH2D> Hist_TrueBkgd_MSCLW;
    vector<TH2D> Hist_Dark_MSCLW;
    vector<TH1D> Hist_TrueBkgd_SR_Energy;
    vector<TH1D> Hist_Data_SR_Energy;
    vector<TH1D> Hist_Data_CR_Energy;
    vector<TH1D> Hist_Data_SR_RoI_Energy;
    vector<TH1D> Hist_Data_CR_RoI_Energy;
    vector<TH1D> Hist_Dark_SR_Energy;
    vector<TH1D> Hist_Dark_CR_Energy;
    vector<TH1D> Hist_TrueBkgd_SR_Skymap_Theta2;
    vector<TH1D> Hist_Data_SR_Skymap_Theta2;
    vector<TH1D> Hist_Data_CR_Skymap_Theta2;
    vector<TH1D> Hist_Data_CR_Skymap_Theta2_Raw;
    vector<TH1D> Hist_Data_CR_CameraFoV_Theta2;
    vector<TH1D> Hist_Data_CR_CameraFoV_Theta2_Raw;
    vector<TH2D> Hist_Data_SR_Skymap;
    vector<TH2D> Hist_Data_CR_Skymap;
    vector<TH2D> Hist_Data_CR_Skymap_Raw;
    vector<TH2D> Hist_Data_SR_Skymap_Galactic;
    vector<TH2D> Hist_Data_CR_Skymap_Galactic;
    vector<TH2D> Hist_Data_SR_CameraFoV;
    vector<TH2D> Hist_Data_CR_CameraFoV;
    vector<TH2D> Hist_Data_CR_CameraFoV_Raw;
    vector<TH2D> Hist_Dark_SR_CameraFoV;
    vector<TH2D> Hist_Dark_CR_CameraFoV;
    vector<TH1D> Hist_Dark_SR_Theta2;
    vector<TH1D> Hist_Dark_CR_Theta2;
    for (int e=0;e<N_energy_bins;e++) 
    {

        N_bins_for_deconv = N_bins_for_deconv_at_E[e];

        char e_low[50];
        sprintf(e_low, "%i", int(energy_bins[e]));
        char e_up[50];
        sprintf(e_up, "%i", int(energy_bins[e+1]));
        Hist_GammaDark_MSCLW.push_back(TH2D("Hist_GammaDark_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_GammaDark_MSCLW.at(e).SetBinErrorOption(TH1::kPoisson);
        Hist_GammaMC_MSCLW.push_back(TH2D("Hist_GammaMC_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_GammaMC_MSCLW.at(e).SetBinErrorOption(TH1::kPoisson);
        Hist_GammaData_MSCLW.push_back(TH2D("Hist_GammaData_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_GammaData_MSCLW.at(e).SetBinErrorOption(TH1::kPoisson);
        Hist_GammaDataON_MSCLW.push_back(TH2D("Hist_GammaDataON_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_GammaDataON_MSCLW.at(e).SetBinErrorOption(TH1::kPoisson);
        Hist_GammaDataOFF_MSCLW.push_back(TH2D("Hist_GammaDataOFF_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_GammaDataOFF_MSCLW.at(e).SetBinErrorOption(TH1::kPoisson);
        Hist_Data_MSCLW.push_back(TH2D("Hist_Data_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_Data_MSCLW.at(e).SetBinErrorOption(TH1::kPoisson);
        Hist_TrueBkgd_MSCLW.push_back(TH2D("Hist_TrueBkgd_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_TrueBkgd_MSCLW.at(e).SetBinErrorOption(TH1::kPoisson);
        Hist_Dark_MSCLW.push_back(TH2D("Hist_Dark_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_Dark_MSCLW.at(e).SetBinErrorOption(TH1::kPoisson);
        Hist_TrueBkgd_SR_Energy.push_back(TH1D("Hist_TrueBkgd_SR_Energy_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_energy_fine_bins,energy_fine_bins));
        Hist_Data_SR_Energy.push_back(TH1D("Hist_Data_SR_Energy_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_energy_fine_bins,energy_fine_bins));
        Hist_Data_CR_Energy.push_back(TH1D("Hist_Data_CR_Energy_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_energy_fine_bins,energy_fine_bins));
        Hist_Data_SR_RoI_Energy.push_back(TH1D("Hist_Data_SR_RoI_Energy_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_energy_fine_bins,energy_fine_bins));
        Hist_Data_CR_RoI_Energy.push_back(TH1D("Hist_Data_CR_RoI_Energy_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_energy_fine_bins,energy_fine_bins));
        Hist_Dark_SR_Energy.push_back(TH1D("Hist_Dark_SR_Energy_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_energy_fine_bins,energy_fine_bins));
        Hist_Dark_CR_Energy.push_back(TH1D("Hist_Dark_CR_Energy_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_energy_fine_bins,energy_fine_bins));
        Hist_TrueBkgd_SR_Skymap_Theta2.push_back(TH1D("Hist_TrueBkgd_SR_Skymap_Theta2_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",50,0,10));
        Hist_Data_SR_Skymap_Theta2.push_back(TH1D("Hist_Data_SR_Skymap_Theta2_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",50,0,10));
        Hist_Data_CR_Skymap_Theta2.push_back(TH1D("Hist_Data_CR_Skymap_Theta2_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",50,0,10));
        Hist_Data_CR_Skymap_Theta2_Raw.push_back(TH1D("Hist_Data_CR_Skymap_Theta2_Raw_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",50,0,10));
        Hist_Data_CR_CameraFoV_Theta2.push_back(TH1D("Hist_Data_CR_CameraFoV_Theta2_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",50,0,10));
        Hist_Data_CR_CameraFoV_Theta2_Raw.push_back(TH1D("Hist_Data_CR_CameraFoV_Theta2_Raw_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",50,0,10));
        Hist_Data_SR_Skymap.push_back(TH2D("Hist_Data_SR_Skymap_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",150,mean_tele_point_ra-3,mean_tele_point_ra+3,150,mean_tele_point_dec-3,mean_tele_point_dec+3));
        Hist_Data_CR_Skymap.push_back(TH2D("Hist_Data_CR_Skymap_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",150,mean_tele_point_ra-3,mean_tele_point_ra+3,150,mean_tele_point_dec-3,mean_tele_point_dec+3));
        Hist_Data_CR_Skymap_Raw.push_back(TH2D("Hist_Data_CR_Skymap_Raw_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",150,mean_tele_point_ra-3,mean_tele_point_ra+3,150,mean_tele_point_dec-3,mean_tele_point_dec+3));

        pair<double,double> tele_point_l_b = ConvertRaDecToGalactic(mean_tele_point_ra, mean_tele_point_dec);
        mean_tele_point_l = tele_point_l_b.first;
        mean_tele_point_b = tele_point_l_b.second;
        Hist_Data_SR_Skymap_Galactic.push_back(TH2D("Hist_Data_SR_Skymap_Galactic_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",150,tele_point_l_b.first-3,tele_point_l_b.first+3,150,tele_point_l_b.second-3,tele_point_l_b.second+3));
        Hist_Data_CR_Skymap_Galactic.push_back(TH2D("Hist_Data_CR_Skymap_Galactic_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",150,tele_point_l_b.first-3,tele_point_l_b.first+3,150,tele_point_l_b.second-3,tele_point_l_b.second+3));

        Hist_Data_SR_CameraFoV.push_back(TH2D("Hist_Data_SR_CameraFoV_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",150,-3,3,150,-3,3));
        Hist_Data_CR_CameraFoV.push_back(TH2D("Hist_Data_CR_CameraFoV_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",150,-3,3,150,-3,3));
        Hist_Data_CR_CameraFoV_Raw.push_back(TH2D("Hist_Data_CR_CameraFoV_Raw_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",150,-3,3,150,-3,3));
        Hist_Dark_SR_CameraFoV.push_back(TH2D("Hist_Dark_SR_CameraFoV_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",20,0,5,6,0,2*M_PI));
        Hist_Dark_CR_CameraFoV.push_back(TH2D("Hist_Dark_CR_CameraFoV_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",20,0,5,6,0,2*M_PI));
        Hist_Dark_SR_Theta2.push_back(TH1D("Hist_Dark_SR_Theta2_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",20,0,10));
        Hist_Dark_CR_Theta2.push_back(TH1D("Hist_Dark_CR_Theta2_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",20,0,10));
    }

    //std::cout << "Prepare dark run samples..." << std::endl;
    for (int run=0;run<Dark_runlist.size();run++)
    {
        char run_number[50];
        char Dark_observation[50];
        sprintf(run_number, "%i", int(Dark_runlist[run].second));
        sprintf(Dark_observation, "%s", Dark_runlist[run].first.c_str());
        //std::cout << "Use dark run " << Dark_runlist[run].first << std::endl;
        string filename;
        filename = TString("$VERITAS_USER_DATA_DIR/"+TString(Dark_observation)+"_V6_Moderate-TMVA-BDT.RB."+TString(run_number)+".root");

        pair<double,double> onrun_tele_point_ra_dec = std::make_pair(0,0);
        if (!TString(target).Contains("Proton")) onrun_tele_point_ra_dec = ON_pointing_radec_for_DarkRun[run];
        pair<double,double> source_ra_dec = GetSourceRaDec(TString(target));
        double darkrun_theta2_off = pow(onrun_tele_point_ra_dec.first-source_ra_dec.first,2)+pow(onrun_tele_point_ra_dec.second-source_ra_dec.second,2);
        DarkRun_theta2_lower = max(0.,darkrun_theta2_off-Theta2_cut_upper);
        DarkRun_theta2_upper = darkrun_theta2_off+Theta2_cut_upper;

        TFile*  input_file = TFile::Open(filename.c_str());
	TH1* i_hEffAreaP = ( TH1* )getEffAreaHistogram(input_file,Dark_runlist[run].second);
        TString root_file = "run_"+TString(run_number)+"/stereo/data_on";
        TTree* Dark_tree = (TTree*) input_file->Get(root_file);
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
        Dark_tree->SetBranchAddress("Xcore",&Xcore);
        Dark_tree->SetBranchAddress("Ycore",&Ycore);
        Dark_tree->SetBranchAddress("SizeSecondMax",&SizeSecondMax);
        Dark_tree->SetBranchAddress("Time",&Time);
        Dark_tree->SetBranchAddress("Shower_Ze",&Shower_Ze);
        Dark_tree->SetBranchAddress("Shower_Az",&Shower_Az);

        Dark_tree->GetEntry(0);
        double time_0 = Time;
        Dark_tree->GetEntry(Dark_tree->GetEntries()-1);
        double time_1 = Time;
        exposure_hours_dark += (time_1-time_0)/3600.;
        double NSB_thisrun = GetRunNSB(int(Dark_runlist[run].second));
        if (TString(Dark_observation).Contains("NSB075")) NSB_thisrun = 3.30;
        if (TString(Dark_observation).Contains("NSB100")) NSB_thisrun = 3.84;
        if (TString(Dark_observation).Contains("NSB150")) NSB_thisrun = 4.69;
        if (TString(Dark_observation).Contains("NSB200")) NSB_thisrun = 5.41;
        if (TString(Dark_observation).Contains("NSB250")) NSB_thisrun = 6.00;
        if (TString(Dark_observation).Contains("NSB325")) NSB_thisrun = 6.83;
        if (TString(Dark_observation).Contains("NSB425")) NSB_thisrun = 7.73;
        if (TString(Dark_observation).Contains("NSB550")) NSB_thisrun = 8.75;
        if (TString(Dark_observation).Contains("NSB750")) NSB_thisrun = 10.18;
        if (TString(Dark_observation).Contains("NSB1000")) NSB_thisrun = 11.77;
        NSB_avg_dark += (time_1-time_0)/3600.*NSB_thisrun;
        Hist_Dark_NSB.Fill(NSB_thisrun);

        for (int entry=0;entry<Dark_tree->GetEntries();entry++) 
        {
            ErecS = 0;
            EChi2S = 0;
            NImages = 0;
            Xcore = 0;
            Ycore = 0;
            SizeSecondMax = 0;
            MSCW = 0;
            MSCL = 0;
            R2off = 0;
            Dark_tree->GetEntry(entry);
            R2off = Xoff*Xoff+Yoff*Yoff;
            Phioff = atan2(Yoff,Xoff)+M_PI;
            int energy = Hist_ErecS.FindBin(ErecS*1000.)-1;
            int energy_fine = Hist_ErecS_fine.FindBin(ErecS*1000.)-1;
            if (energy<0) continue;
            if (energy>=N_energy_bins) continue;
            int e = energy;
            if (!SelectNImages(3,4)) continue;
            if (SizeSecondMax<600.) continue;
            if (pow(Xcore*Xcore+Ycore*Ycore,0.5)>350) continue;
            if (R2off>4.) continue;
            if (Dark_runlist[run].first.find("Proton")!=std::string::npos && NSB_diff_dark==0.) 
            {
                //if (Data_runlist[0].first.find("NSB")==std::string::npos && abs(Shower_Az-180.)<90.) continue;
                if (abs(Shower_Az-180.)<90.) continue;
            }
            Hist_Dark_ShowerDirection.Fill(Shower_Az,Shower_Ze);
            if (DarkFoV() || Dark_runlist[run].first.find("Proton")!=std::string::npos)
            {
                Hist_Dark_MSCLW.at(e).Fill(MSCL,MSCW);
                if (SignalSelectionTheta2())
                {
                    Hist_Dark_SR_CameraFoV.at(e).Fill(R2off,Phioff);
                    Hist_Dark_SR_Theta2.at(e).Fill(Xoff*Xoff+Yoff*Yoff);
                    Hist_Dark_SR_Energy.at(e).Fill(ErecS*1000.);
                }
                if (ControlSelectionTheta2())
                {
                    Hist_Dark_CR_CameraFoV.at(e).Fill(R2off,Phioff);
                    Hist_Dark_CR_Theta2.at(e).Fill(Xoff*Xoff+Yoff*Yoff);
                    Hist_Dark_CR_Energy.at(e).Fill(ErecS*1000.);
                }
            }
        }
        input_file->Close();
    }



    //std::cout << "Prepare ON run samples..." << std::endl;
    for (int run=0;run<Data_runlist.size();run++)
    {
        char run_number[50];
        char Data_observation[50];
        sprintf(run_number, "%i", int(Data_runlist[run].second));
        sprintf(Data_observation, "%s", Data_runlist[run].first.c_str());
        string filename;
        filename = TString("$VERITAS_USER_DATA_DIR/"+TString(Data_observation)+"_V6_Moderate-TMVA-BDT.RB."+TString(run_number)+".root");

        //std::cout << "Get telescope pointing RA and Dec..." << std::endl;
        pair<double,double> tele_point_ra_dec = std::make_pair(0,0);
        if (!TString(target).Contains("Proton")) tele_point_ra_dec = GetRunRaDec(filename,int(Data_runlist[run].second));

        TFile*  input_file = TFile::Open(filename.c_str());
        TString root_file = "run_"+TString(run_number)+"/stereo/data_on";
        TTree* Data_tree = (TTree*) input_file->Get(root_file);
        Data_tree->SetBranchAddress("Xoff",&Xoff);
        Data_tree->SetBranchAddress("Yoff",&Yoff);
        Data_tree->SetBranchAddress("theta2",&theta2);
        Data_tree->SetBranchAddress("Xoff_derot",&Xoff_derot);
        Data_tree->SetBranchAddress("Yoff_derot",&Yoff_derot);
        Data_tree->SetBranchAddress("ErecS",&ErecS);
        Data_tree->SetBranchAddress("EChi2S",&EChi2S);
        Data_tree->SetBranchAddress("MSCW",&MSCW);
        Data_tree->SetBranchAddress("MSCL",&MSCL);
        Data_tree->SetBranchAddress("NImages",&NImages);
        Data_tree->SetBranchAddress("Xcore",&Xcore);
        Data_tree->SetBranchAddress("Ycore",&Ycore);
        Data_tree->SetBranchAddress("SizeSecondMax",&SizeSecondMax);
        Data_tree->SetBranchAddress("Time",&Time);
        Data_tree->SetBranchAddress("Shower_Ze",&Shower_Ze);
        Data_tree->SetBranchAddress("Shower_Az",&Shower_Az);

        // Get effective area and livetime and determine the cosmic electron counts for this run.
        //std::cout << "Get effective area and livetime..." << std::endl;
	TH1* i_hEffAreaP = ( TH1* )getEffAreaHistogram(input_file,Data_runlist[run].second);
        Data_tree->GetEntry(0);
        double time_0 = Time;
        Data_tree->GetEntry(Data_tree->GetEntries()-1);
        double time_1 = Time;
        exposure_hours += (time_1-time_0)/3600.;
        double NSB_thisrun = GetRunNSB(int(Data_runlist[run].second));
        if (TString(Data_observation).Contains("NSB075")) NSB_thisrun = 3.30;
        if (TString(Data_observation).Contains("NSB100")) NSB_thisrun = 3.84;
        if (TString(Data_observation).Contains("NSB150")) NSB_thisrun = 4.69;
        if (TString(Data_observation).Contains("NSB200")) NSB_thisrun = 5.41;
        if (TString(Data_observation).Contains("NSB250")) NSB_thisrun = 6.00;
        if (TString(Data_observation).Contains("NSB325")) NSB_thisrun = 6.83;
        if (TString(Data_observation).Contains("NSB425")) NSB_thisrun = 7.73;
        if (TString(Data_observation).Contains("NSB550")) NSB_thisrun = 8.75;
        if (TString(Data_observation).Contains("NSB750")) NSB_thisrun = 10.18;
        if (TString(Data_observation).Contains("NSB1000")) NSB_thisrun = 11.77;
        NSB_avg += (time_1-time_0)/3600.*NSB_thisrun;
        Hist_Data_NSB.Fill(NSB_thisrun);

        for (int e=0;e<N_energy_fine_bins;e++) 
        {
            double eff_area = i_hEffAreaP->GetBinContent( i_hEffAreaP->FindBin( log10(0.5*(energy_fine_bins[e]+energy_fine_bins[e+1])/1000.)));
            //std::pair <double,double> mcgillflux = GetMcGillElectronFlux((energy_fine_bins[e+1]+energy_fine_bins[e])/2.);
            //gamma_flux[e] = mcgillflux.first;
            gamma_flux[e] = PercentCrab/100.*GetCrabFlux((energy_fine_bins[e+1]+energy_fine_bins[e])/2.);
            //std::cout << "energy = " << (energy_fine_bins[e+1]+energy_fine_bins[e])/2. << std::endl;
            //std::cout << "gamma_flux[e] = " << gamma_flux[e] << std::endl;
            double expected_electrons = gamma_flux[e]*eff_area*(time_1-time_0)*(energy_fine_bins[e+1]-energy_fine_bins[e])/1000.;
            //std::cout << "expected_electrons = " << expected_electrons << std::endl;
            gamma_count[e] += expected_electrons; // this is used to normalize MC electron template.

            Hist_EffArea.SetBinContent(e+1,Hist_EffArea.GetBinContent(e+1)+eff_area*(time_1-time_0));

        }
        //std::cout << "Finished getting effective area and livetime..." << std::endl;


        for (int entry=0;entry<Data_tree->GetEntries();entry++) 
        {
            ErecS = 0;
            EChi2S = 0;
            NImages = 0;
            Xcore = 0;
            Ycore = 0;
            SizeSecondMax = 0;
            MSCW = 0;
            MSCL = 0;
            R2off = 0;
            Data_tree->GetEntry(entry);
            R2off = Xoff*Xoff+Yoff*Yoff;
            Phioff = atan2(Yoff,Xoff)+M_PI;
            ra_sky = tele_point_ra_dec.first+Xoff_derot;
            dec_sky = tele_point_ra_dec.second+Yoff_derot;
            if (TString(target).Contains("Everything"))
            {
                ra_sky = Xoff;
                dec_sky = Yoff;
            }
            pair<double,double> evt_l_b = ConvertRaDecToGalactic(ra_sky,dec_sky);
            int energy = Hist_ErecS.FindBin(ErecS*1000.)-1;
            int energy_fine = Hist_ErecS_fine.FindBin(ErecS*1000.)-1;
            if (energy<0) continue;
            if (energy>=N_energy_bins) continue;
            int e = energy;
            if (!SelectNImages(3,4)) continue;
            if (SizeSecondMax<600.) continue;
            if (pow(Xcore*Xcore+Ycore*Ycore,0.5)>350) continue;
            if (R2off>4.) continue;
            if (Data_runlist[run].first.find("Proton")!=std::string::npos && NSB_diff_dark==0.) 
            {
                //if (Data_runlist[run].first.find("NSB")==std::string::npos && abs(Shower_Az-180.)>90.) continue;
                if (abs(Shower_Az-180.)>90.) continue;
            }
            Hist_Data_ShowerDirection.Fill(Shower_Az,Shower_Ze);
            if (FoV() || Data_runlist[run].first.find("Proton")!=std::string::npos)
            {
                Hist_Data_MSCLW_incl.Fill(MSCL,MSCW);
                Hist_Data_MSCLW.at(e).Fill(MSCL,MSCW);
                Hist_TrueBkgd_MSCLW.at(e).Fill(MSCL,MSCW);
            }
            if (SignalSelectionTheta2())
            {
                if (FoV() || Data_runlist[run].first.find("Proton")!=std::string::npos)
                {
                    Hist_TrueBkgd_SR_Energy.at(e).Fill(ErecS*1000.);
                    Hist_Data_SR_Energy.at(e).Fill(ErecS*1000.);
                    if (RoIFoV()) Hist_Data_SR_RoI_Energy.at(e).Fill(ErecS*1000.);
                    Hist_TrueBkgd_SR_Skymap_Theta2.at(e).Fill(theta2);
                    Hist_Data_SR_Skymap_Theta2.at(e).Fill(theta2);
                    Hist_Data_SR_Skymap.at(e).Fill(ra_sky,dec_sky);
                    Hist_Data_SR_CameraFoV.at(e).Fill(Xoff,Yoff);
                    Hist_Data_SR_Skymap_Galactic.at(e).Fill(evt_l_b.first,evt_l_b.second);
                }
            }
            if (ControlSelectionTheta2())
            {
                if (FoV() || Data_runlist[run].first.find("Proton")!=std::string::npos)
                {
                    int binx = Hist_Dark_SR_CameraFoV.at(e).GetXaxis()->FindBin(R2off);
                    int biny = Hist_Dark_SR_CameraFoV.at(e).GetYaxis()->FindBin(Phioff);
                    double dark_cr_content = Hist_Dark_CR_CameraFoV.at(e).GetBinContent(binx,biny);
                    double dark_sr_content = Hist_Dark_SR_CameraFoV.at(e).GetBinContent(binx,biny);
                    double weight = 0.;
                    if (dark_cr_content>0.) weight = dark_sr_content/dark_cr_content;
                    int bin_e = Hist_Dark_SR_Energy.at(e).FindBin(ErecS*1000.);
                    double dark_cr_content_e = Hist_Dark_CR_Energy.at(e).GetBinContent(bin_e);
                    double dark_sr_content_e = Hist_Dark_SR_Energy.at(e).GetBinContent(bin_e);
                    double weight_e = 0.;
                    if (dark_cr_content_e>0.) weight_e = dark_sr_content_e/dark_cr_content_e;
                    Hist_Data_CR_Energy.at(e).Fill(ErecS*1000.,weight_e);
                    if (RoIFoV()) Hist_Data_CR_RoI_Energy.at(e).Fill(ErecS*1000.,weight_e);
                    Hist_Data_CR_Skymap_Theta2.at(e).Fill(theta2,weight);
                    Hist_Data_CR_Skymap_Theta2_Raw.at(e).Fill(theta2,1.);
                    Hist_Data_CR_CameraFoV_Theta2.at(e).Fill(R2off,weight);
                    Hist_Data_CR_CameraFoV_Theta2_Raw.at(e).Fill(R2off,1.);
                    Hist_Data_CR_Skymap.at(e).Fill(ra_sky,dec_sky,weight);
                    Hist_Data_CR_Skymap_Raw.at(e).Fill(ra_sky,dec_sky,1.);
                    Hist_Data_CR_Skymap_Galactic.at(e).Fill(evt_l_b.first,evt_l_b.second,weight);
                    Hist_Data_CR_CameraFoV.at(e).Fill(Xoff,Yoff,weight);
                    Hist_Data_CR_CameraFoV_Raw.at(e).Fill(Xoff,Yoff,1.);
                }
            }
        }
        input_file->Close();
    }


    //std::cout << "Prepare photon MC samples..." << std::endl;
    double n_photon = 0.;
    for (int run=0;run<PhotonMC_runlist.size();run++)
    {
        char run_number[50];
        char Data_observation[50];
        sprintf(run_number, "%i", int(PhotonMC_runlist[run].second));
        sprintf(Data_observation, "%s", PhotonMC_runlist[run].first.c_str());
        string filename;
        filename = TString("$VERITAS_USER_DATA_DIR/"+TString(Data_observation)+"_V6_Moderate-TMVA-BDT.RB."+TString(run_number)+".root");
        TFile*  input_file = TFile::Open(filename.c_str());
        TString root_file = "run_"+TString(run_number)+"/stereo/data_on";
        TTree* Data_tree = (TTree*) input_file->Get(root_file);
        Data_tree->SetBranchAddress("Xoff",&Xoff);
        Data_tree->SetBranchAddress("Yoff",&Yoff);
        Data_tree->SetBranchAddress("theta2",&theta2);
        Data_tree->SetBranchAddress("ra",&ra_sky);
        Data_tree->SetBranchAddress("dec",&dec_sky);
        Data_tree->SetBranchAddress("ErecS",&ErecS);
        Data_tree->SetBranchAddress("EChi2S",&EChi2S);
        Data_tree->SetBranchAddress("MSCW",&MSCW);
        Data_tree->SetBranchAddress("MSCL",&MSCL);
        Data_tree->SetBranchAddress("NImages",&NImages);
        Data_tree->SetBranchAddress("Xcore",&Xcore);
        Data_tree->SetBranchAddress("Ycore",&Ycore);
        Data_tree->SetBranchAddress("SizeSecondMax",&SizeSecondMax);
        Data_tree->SetBranchAddress("Time",&Time);
        Data_tree->SetBranchAddress("Shower_Ze",&Shower_Ze);
        Data_tree->SetBranchAddress("Shower_Az",&Shower_Az);

        for (int entry=0;entry<Data_tree->GetEntries();entry++) 
        {
            ErecS = 0;
            EChi2S = 0;
            NImages = 0;
            Xcore = 0;
            Ycore = 0;
            SizeSecondMax = 0;
            MSCW = 0;
            MSCL = 0;
            R2off = 0;
            Data_tree->GetEntry(entry);
            R2off = Xoff*Xoff+Yoff*Yoff;
            ra_sky = mean_tele_point_ra+ra_sky;
            dec_sky = mean_tele_point_dec+dec_sky;
            int energy = Hist_ErecS.FindBin(ErecS*1000.)-1;
            int energy_fine = Hist_ErecS_fine.FindBin(ErecS*1000.)-1;
            if (energy<0) continue;
            if (energy>=N_energy_bins) continue;
            int e = energy;
            if (!SelectNImages(3,4)) continue;
            if (FoV() && GammaFoV()) raw_gamma_count[energy_fine] += 1.;
        }

    }

    double photon_weight = 1.0;
    for (int run=0;run<PhotonMC_runlist.size();run++)
    {
        char run_number[50];
        char Data_observation[50];
        sprintf(run_number, "%i", int(PhotonMC_runlist[run].second));
        sprintf(Data_observation, "%s", PhotonMC_runlist[run].first.c_str());
        string filename;
        filename = TString("$VERITAS_USER_DATA_DIR/"+TString(Data_observation)+"_V6_Moderate-TMVA-BDT.RB."+TString(run_number)+".root");


        TFile*  input_file = TFile::Open(filename.c_str());
	TH1* i_hEffAreaP = ( TH1* )getEffAreaHistogram(input_file,PhotonMC_runlist[run].second);
        TString root_file = "run_"+TString(run_number)+"/stereo/data_on";
        TTree* Data_tree = (TTree*) input_file->Get(root_file);
        Data_tree->SetBranchAddress("Xoff",&Xoff);
        Data_tree->SetBranchAddress("Yoff",&Yoff);
        Data_tree->SetBranchAddress("theta2",&theta2);
        Data_tree->SetBranchAddress("ra",&ra_sky);
        Data_tree->SetBranchAddress("dec",&dec_sky);
        Data_tree->SetBranchAddress("ErecS",&ErecS);
        Data_tree->SetBranchAddress("EChi2S",&EChi2S);
        Data_tree->SetBranchAddress("MSCW",&MSCW);
        Data_tree->SetBranchAddress("MSCL",&MSCL);
        Data_tree->SetBranchAddress("NImages",&NImages);
        Data_tree->SetBranchAddress("Xcore",&Xcore);
        Data_tree->SetBranchAddress("Ycore",&Ycore);
        Data_tree->SetBranchAddress("SizeSecondMax",&SizeSecondMax);
        Data_tree->SetBranchAddress("Time",&Time);
        Data_tree->SetBranchAddress("Shower_Ze",&Shower_Ze);
        Data_tree->SetBranchAddress("Shower_Az",&Shower_Az);

        // Get effective area and livetime and determine the cosmic electron counts for this run.
        Data_tree->GetEntry(0);
        double time_0 = Time;
        Data_tree->GetEntry(Data_tree->GetEntries()-1);
        double time_1 = Time;

        for (int entry=0;entry<Data_tree->GetEntries();entry++) 
        {
            ErecS = 0;
            EChi2S = 0;
            NImages = 0;
            Xcore = 0;
            Ycore = 0;
            SizeSecondMax = 0;
            MSCW = 0;
            MSCL = 0;
            R2off = 0;
            Data_tree->GetEntry(entry);
            R2off = Xoff*Xoff+Yoff*Yoff;
            Phioff = atan2(Yoff,Xoff)+M_PI;
            ra_sky = mean_tele_point_ra+ra_sky;
            dec_sky = mean_tele_point_dec+dec_sky;
            if (TString(target).Contains("Everything"))
            {
                ra_sky = Xoff;
                dec_sky = Yoff;
            }
            pair<double,double> evt_l_b = ConvertRaDecToGalactic(ra_sky,dec_sky);
            int energy = Hist_ErecS.FindBin(ErecS*1000.)-1;
            int energy_fine = Hist_ErecS_fine.FindBin(ErecS*1000.)-1;
            photon_weight = gamma_count[energy_fine]/raw_gamma_count[energy_fine];
            if (energy<0) continue;
            if (energy>=N_energy_bins) continue;
            int e = energy;
            if (!SelectNImages(3,4)) continue;
            if (theta2<0.3) Hist_GammaMC_MSCLW.at(e).Fill(MSCL,MSCW,1.);
            Hist_Data_ShowerDirection.Fill(Shower_Az,Shower_Ze,photon_weight);
            if (FoV() && GammaFoV())
            {
                Hist_Data_MSCLW_incl.Fill(MSCL,MSCW,photon_weight);
                Hist_Data_MSCLW.at(e).Fill(MSCL,MSCW,photon_weight);
            }
            if (SignalSelectionTheta2())
            {
                if (FoV() && GammaFoV())
                {
                    Hist_Data_SR_Energy.at(e).Fill(ErecS*1000.,photon_weight);
                    if (RoIFoV()) Hist_Data_SR_RoI_Energy.at(e).Fill(ErecS*1000.,photon_weight);
                    Hist_Data_SR_Skymap_Theta2.at(e).Fill(theta2,photon_weight);
                    Hist_Data_SR_Skymap.at(e).Fill(ra_sky,dec_sky,photon_weight);
                    Hist_Data_SR_Skymap_Galactic.at(e).Fill(evt_l_b.first,evt_l_b.second,photon_weight);
                }
            }
            if (ControlSelectionTheta2())
            {
                if (FoV() && GammaFoV())
                {
                    int binx = Hist_Dark_SR_CameraFoV.at(e).GetXaxis()->FindBin(R2off);
                    int biny = Hist_Dark_SR_CameraFoV.at(e).GetYaxis()->FindBin(Phioff);
                    double dark_cr_content = Hist_Dark_CR_CameraFoV.at(e).GetBinContent(binx,biny);
                    double dark_sr_content = Hist_Dark_SR_CameraFoV.at(e).GetBinContent(binx,biny);
                    double weight = 0.;
                    if (dark_cr_content>0.) weight = dark_sr_content/dark_cr_content;
                    int bin_e = Hist_Dark_SR_Energy.at(e).FindBin(ErecS*1000.);
                    double dark_cr_content_e = Hist_Dark_CR_Energy.at(e).GetBinContent(bin_e);
                    double dark_sr_content_e = Hist_Dark_SR_Energy.at(e).GetBinContent(bin_e);
                    double weight_e = 0.;
                    if (dark_cr_content_e>0.) weight_e = dark_sr_content_e/dark_cr_content_e;
                    Hist_Data_CR_Energy.at(e).Fill(ErecS*1000.,weight_e*photon_weight);
                    if (RoIFoV()) Hist_Data_CR_RoI_Energy.at(e).Fill(ErecS*1000.,weight_e*photon_weight);
                    Hist_Data_CR_Skymap_Theta2.at(e).Fill(theta2,weight*photon_weight);
                    Hist_Data_CR_Skymap_Theta2_Raw.at(e).Fill(theta2,photon_weight);
                    Hist_Data_CR_CameraFoV_Theta2.at(e).Fill(R2off,weight*photon_weight);
                    Hist_Data_CR_CameraFoV_Theta2_Raw.at(e).Fill(R2off,photon_weight);
                    Hist_Data_CR_Skymap.at(e).Fill(ra_sky,dec_sky,weight*photon_weight);
                    Hist_Data_CR_Skymap_Raw.at(e).Fill(ra_sky,dec_sky,photon_weight);
                    Hist_Data_CR_Skymap_Galactic.at(e).Fill(evt_l_b.first,evt_l_b.second,weight*photon_weight);
                }
            }
        }
        input_file->Close();
    }

    // Get Gamma ray data template
    for (int run=0;run<PhotonData_runlist.size();run++)
    {
        char run_number[50];
        char Data_observation[50];
        sprintf(run_number, "%i", int(PhotonData_runlist[run].second));
        sprintf(Data_observation, "%s", PhotonData_runlist[run].first.c_str());
        string filename;
        filename = TString("$VERITAS_USER_DATA_DIR/"+TString(Data_observation)+"_V6_Moderate-TMVA-BDT.RB."+TString(run_number)+".root");


        TFile*  input_file = TFile::Open(filename.c_str());
	TH1* i_hEffAreaP = ( TH1* )getEffAreaHistogram(input_file,PhotonData_runlist[run].second);
        TString root_file = "run_"+TString(run_number)+"/stereo/data_on";
        TTree* Data_tree = (TTree*) input_file->Get(root_file);
        Data_tree->SetBranchAddress("Xoff",&Xoff);
        Data_tree->SetBranchAddress("Yoff",&Yoff);
        Data_tree->SetBranchAddress("theta2",&theta2);
        Data_tree->SetBranchAddress("ra",&ra_sky);
        Data_tree->SetBranchAddress("dec",&dec_sky);
        Data_tree->SetBranchAddress("ErecS",&ErecS);
        Data_tree->SetBranchAddress("EChi2S",&EChi2S);
        Data_tree->SetBranchAddress("MSCW",&MSCW);
        Data_tree->SetBranchAddress("MSCL",&MSCL);
        Data_tree->SetBranchAddress("NImages",&NImages);
        Data_tree->SetBranchAddress("Xcore",&Xcore);
        Data_tree->SetBranchAddress("Ycore",&Ycore);
        Data_tree->SetBranchAddress("SizeSecondMax",&SizeSecondMax);
        Data_tree->SetBranchAddress("Time",&Time);
        Data_tree->SetBranchAddress("Shower_Ze",&Shower_Ze);
        Data_tree->SetBranchAddress("Shower_Az",&Shower_Az);

        // Get effective area and livetime and determine the cosmic electron counts for this run.
        Data_tree->GetEntry(0);
        double time_0 = Time;
        Data_tree->GetEntry(Data_tree->GetEntries()-1);
        double time_1 = Time;

        for (int entry=0;entry<Data_tree->GetEntries();entry++) 
        {
            ErecS = 0;
            EChi2S = 0;
            NImages = 0;
            Xcore = 0;
            Ycore = 0;
            SizeSecondMax = 0;
            MSCW = 0;
            MSCL = 0;
            R2off = 0;
            Data_tree->GetEntry(entry);
            R2off = Xoff*Xoff+Yoff*Yoff;
            int energy = Hist_ErecS.FindBin(ErecS*1000.)-1;
            int energy_fine = Hist_ErecS_fine.FindBin(ErecS*1000.)-1;
            if (energy<0) continue;
            if (energy>=N_energy_bins) continue;
            int e = energy;
            if (!SelectNImages(3,4)) continue;
            if (SizeSecondMax<600.) continue;
            if (pow(Xcore*Xcore+Ycore*Ycore,0.5)>350) continue;
            Hist_Data_ShowerDirection.Fill(Shower_Az,Shower_Ze,photon_weight);
            if (theta2<0.3)
            {
                Hist_GammaDataON_MSCLW.at(e).Fill(MSCL,MSCW);
            }
            else if (theta2>0.3 && theta2<0.5)
            {
                Hist_GammaDataOFF_MSCLW.at(e).Fill(MSCL,MSCW);
            }
        }
        input_file->Close();
    }
    for (int e=0;e<N_energy_bins;e++) 
    {
        int binx_blind = Hist_GammaDataON_MSCLW.at(e).GetXaxis()->FindBin(1.);
        int binx_lower = Hist_GammaDataON_MSCLW.at(e).GetXaxis()->FindBin(MSCL_plot_lower);
        int biny_blind = Hist_GammaDataON_MSCLW.at(e).GetYaxis()->FindBin(1.);
        int biny_lower = Hist_GammaDataON_MSCLW.at(e).GetYaxis()->FindBin(MSCW_plot_lower);
        double GammaDataON_SR_Integral = Hist_GammaDataON_MSCLW.at(e).Integral(binx_lower,binx_blind,biny_lower,biny_blind);
        double GammaDataOFF_SR_Integral = Hist_GammaDataOFF_MSCLW.at(e).Integral(binx_lower,binx_blind,biny_lower,biny_blind);
        double GammaDataON_all_Integral = Hist_GammaDataON_MSCLW.at(e).Integral();
        double GammaDataOFF_all_Integral = Hist_GammaDataOFF_MSCLW.at(e).Integral();
        double GammaDataON_CR_Integral = GammaDataON_all_Integral-GammaDataON_SR_Integral;
        double GammaDataOFF_CR_Integral = GammaDataOFF_all_Integral-GammaDataOFF_SR_Integral;
        double scale = GammaDataON_CR_Integral/GammaDataOFF_CR_Integral;
        Hist_GammaDataOFF_MSCLW.at(e).Scale(scale);
        Hist_GammaData_MSCLW.at(e).Add(&Hist_GammaDataON_MSCLW.at(e));
        Hist_GammaData_MSCLW.at(e).Add(&Hist_GammaDataOFF_MSCLW.at(e),-1.);
        for (int binx=0;binx<Hist_GammaData_MSCLW.at(e).GetNbinsX();binx++)
        {
            for (int biny=0;biny<Hist_GammaData_MSCLW.at(e).GetNbinsY();biny++)
            {
                double old_content = Hist_GammaData_MSCLW.at(e).GetBinContent(binx+1,biny+1);
                double old_error = Hist_GammaData_MSCLW.at(e).GetBinError(binx+1,biny+1);
                if (old_content<0)
                {
                    Hist_GammaData_MSCLW.at(e).SetBinContent(binx+1,biny+1,0);
                    Hist_GammaData_MSCLW.at(e).SetBinError(binx+1,biny+1,0);
                }
            }
        }
    }


    for (int e=0;e<N_energy_bins;e++) 
    {
        int binx_lower = Hist_Dark_MSCLW.at(e).GetXaxis()->FindBin(MSCL_cut_lower);
        int binx_blind = Hist_Dark_MSCLW.at(e).GetXaxis()->FindBin(MSCL_cut_blind)-1;
        int binx_upper = Hist_Dark_MSCLW.at(e).GetXaxis()->FindBin(1.)-1;
        int biny_lower = Hist_Dark_MSCLW.at(e).GetYaxis()->FindBin(MSCW_cut_lower);
        int biny_blind = Hist_Dark_MSCLW.at(e).GetYaxis()->FindBin(MSCW_cut_blind)-1;
        int biny_upper = Hist_Dark_MSCLW.at(e).GetYaxis()->FindBin(1.)-1;
        double Dark_SR_Integral = Hist_Dark_MSCLW.at(e).Integral(binx_lower,binx_blind,biny_lower,biny_blind);
        double Data_SR_Integral = Hist_Data_MSCLW.at(e).Integral(binx_lower,binx_blind,biny_lower,biny_blind);
        //double Dark_SR_Integral = Hist_Dark_MSCLW.at(e).Integral(binx_lower,binx_upper,biny_lower,biny_upper);
        //double Data_SR_Integral = Hist_Data_MSCLW.at(e).Integral(binx_lower,binx_upper,biny_lower,biny_upper);
        double Dark_Integral = Hist_Dark_MSCLW.at(e).Integral();
        double Data_Integral = Hist_Data_MSCLW.at(e).Integral();
        double Dark_CR_Integral = Dark_Integral-Dark_SR_Integral;
        double Data_CR_Integral = Data_Integral-Data_SR_Integral;
        double Dark_CR_Error = pow(Dark_CR_Integral,0.5);
        double Data_CR_Error = pow(Data_CR_Integral,0.5);
        double scale = Data_CR_Integral/Dark_CR_Integral;
        double scale_err = scale*pow(pow(Dark_CR_Error/Dark_CR_Integral,2)+pow(Data_CR_Error/Data_CR_Integral,2),0.5);
        Hist_Dark_MSCLW.at(e).Scale(scale);
        double gamma_total = Hist_Data_MSCLW.at(e).Integral(binx_lower,binx_blind,biny_lower,biny_blind)-Hist_Dark_MSCLW.at(e).Integral(binx_lower,binx_blind,biny_lower,biny_blind);
        gamma_total = max(0.,gamma_total);
        if (PercentCrab>0.)
        {
            Hist_GammaDark_MSCLW.at(e).Reset();
            Hist_GammaDark_MSCLW.at(e).Add(&Hist_GammaMC_MSCLW.at(e));
            double scale_gamma = double(gamma_total)/double(Hist_GammaDark_MSCLW.at(e).Integral(binx_lower,binx_blind,biny_lower,biny_blind));
            Hist_GammaDark_MSCLW.at(e).Scale(scale_gamma);
        }
        else
        {
            Hist_GammaDark_MSCLW.at(e).Reset();
            Hist_GammaDark_MSCLW.at(e).Add(&Hist_GammaData_MSCLW.at(e));
            double scale_gamma = double(gamma_total)/double(Hist_GammaDark_MSCLW.at(e).Integral(binx_lower,binx_blind,biny_lower,biny_blind));
            Hist_GammaDark_MSCLW.at(e).Scale(scale_gamma);
        }
    }

    NSB_avg = NSB_avg/exposure_hours;
    NSB_avg_dark = NSB_avg_dark/exposure_hours_dark;
    for (int e=0;e<N_energy_fine_bins;e++) 
    {
        Hist_EffArea.SetBinContent(e+1,Hist_EffArea.GetBinContent(e+1)/(3600.*exposure_hours));
    }

    TFile OutputFile("../Netflix_"+TString(target)+"_Crab"+std::to_string(int(PercentCrab))+"_TelElev"+std::to_string(int(TelElev_lower))+"to"+std::to_string(int(TelElev_upper))+"_"+file_tag+".root","recreate");
    TTree InfoTree("InfoTree","info tree");
    InfoTree.Branch("exposure_hours",&exposure_hours,"exposure_hours/D");
    InfoTree.Branch("exposure_hours_dark",&exposure_hours_dark,"exposure_hours_dark/D");
    InfoTree.Branch("MSCW_cut_blind",&MSCW_cut_blind,"MSCW_cut_blind/D");
    InfoTree.Branch("MSCL_cut_blind",&MSCL_cut_blind,"MSCL_cut_blind/D");
    InfoTree.Branch("Theta2_cut_lower",&Theta2_cut_lower,"Theta2_cut_lower/D");
    InfoTree.Branch("Theta2_cut_upper",&Theta2_cut_upper,"Theta2_cut_upper/D");
    InfoTree.Branch("NSB",&NSB_avg,"NSB/D");
    InfoTree.Branch("NSB_dark",&NSB_avg_dark,"NSB_dark/D");
    InfoTree.Branch("mean_tele_point_ra",&mean_tele_point_ra,"mean_tele_point_ra/D");
    InfoTree.Branch("mean_tele_point_dec",&mean_tele_point_dec,"mean_tele_point_dec/D");
    InfoTree.Branch("mean_tele_point_l",&mean_tele_point_l,"mean_tele_point_l/D");
    InfoTree.Branch("mean_tele_point_b",&mean_tele_point_b,"mean_tele_point_b/D");
    InfoTree.Fill();
    InfoTree.Write();
    Hist_Dark_NSB.Write();
    Hist_Data_NSB.Write();
    Hist_EffArea.Write();
    Hist_Dark_ShowerDirection.Write();
    Hist_Data_ShowerDirection.Write();
    Hist_Data_MSCLW_incl.Write();
    for (int e=0;e<N_energy_bins;e++)
    {
        Hist_GammaDark_MSCLW.at(e).Write();
        Hist_GammaMC_MSCLW.at(e).Write();
        Hist_GammaData_MSCLW.at(e).Write();
        Hist_Data_MSCLW.at(e).Write();
        Hist_TrueBkgd_MSCLW.at(e).Write();
        Hist_Dark_MSCLW.at(e).Write();
        Hist_TrueBkgd_SR_Energy.at(e).Write();
        Hist_Data_SR_Energy.at(e).Write();
        Hist_Data_CR_Energy.at(e).Write();
        Hist_Data_SR_RoI_Energy.at(e).Write();
        Hist_Data_CR_RoI_Energy.at(e).Write();
        Hist_TrueBkgd_SR_Skymap_Theta2.at(e).Write();
        Hist_Data_SR_Skymap_Theta2.at(e).Write();
        Hist_Data_CR_Skymap_Theta2.at(e).Write();
        Hist_Data_CR_Skymap_Theta2_Raw.at(e).Write();
        Hist_Data_CR_CameraFoV_Theta2.at(e).Write();
        Hist_Data_CR_CameraFoV_Theta2_Raw.at(e).Write();
        Hist_Data_SR_Skymap.at(e).Write();
        Hist_Data_CR_Skymap.at(e).Write();
        Hist_Data_CR_Skymap_Raw.at(e).Write();
        Hist_Data_SR_Skymap_Galactic.at(e).Write();
        Hist_Data_CR_Skymap_Galactic.at(e).Write();
        Hist_Data_SR_CameraFoV.at(e).Write();
        Hist_Data_CR_CameraFoV.at(e).Write();
        Hist_Data_CR_CameraFoV_Raw.at(e).Write();
    }
    OutputFile.Close();

}
