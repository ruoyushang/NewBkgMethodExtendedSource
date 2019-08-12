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

const int N_energy_bins = 12;
double energy_bins[N_energy_bins+1] = {200,237,282,335,398,473,562,794,1122,1585,2239,4467,8913};
bool use_this_energy_bin[N_energy_bins] = {false,false,false,false,true,false,false,false,false,false,false,false};

int N_bins_for_deconv = 40;
double MSCW_plot_lower = -1.;
double MSCW_plot_upper = 3.;
double MSCL_plot_lower = -1.;
double MSCL_plot_upper = 3.;
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
double R2off = 0;
double theta2 = 0;
double ra_sky = 0;
double dec_sky = 0;
double exposure_hours = 0.;

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

vector<pair<string,int>> SelectONRunList(vector<pair<string,int>> Data_runlist, double Elev_cut_lower, double Elev_cut_upper, double Azim_cut_lower, double Azim_cut_upper)
{
    vector<pair<string,int>> new_list;
    for (int run=0;run<Data_runlist.size();run++)
    {
        char run_number[50];
        char Data_observation[50];
        sprintf(run_number, "%i", int(Data_runlist[run].second));
        sprintf(Data_observation, "%s", Data_runlist[run].first.c_str());
        string filename;
        filename = TString("$VERITAS_USER_DATA_DIR/"+TString(Data_observation)+"_V6_Moderate-TMVA-BDT.RB."+TString(run_number)+".root");

        if (!PointingSelection(filename,int(Data_runlist[run].second),TelElev_lower,TelElev_upper,0,360)) continue;
        new_list.push_back(std::make_pair(Data_runlist[run].first,Data_runlist[run].second));
    }
    return new_list;
}

vector<pair<string,int>> SelectOFFRunList(vector<pair<string,int>> ON_runlist, vector<pair<string,int>> OFF_runlist)
{

    vector<pair<double,double>> ON_pointing;
    for (int on_run=0;on_run<ON_runlist.size();on_run++)
    {
        char ON_runnumber[50];
        char ON_observation[50];
        sprintf(ON_runnumber, "%i", int(ON_runlist[on_run].second));
        sprintf(ON_observation, "%s", ON_runlist[on_run].first.c_str());
        string ON_filename;
        ON_filename = TString("$VERITAS_USER_DATA_DIR/"+TString(ON_observation)+"_V6_Moderate-TMVA-BDT.RB."+TString(ON_runnumber)+".root");
        ON_pointing.push_back(GetRunElevAzim(ON_filename,int(ON_runlist[on_run].second)));
    }

    vector<pair<double,double>> OFF_pointing;
    for (int off_run=0;off_run<OFF_runlist.size();off_run++)
    {
        char OFF_runnumber[50];
        char OFF_observation[50];
        sprintf(OFF_runnumber, "%i", int(OFF_runlist[off_run].second));
        sprintf(OFF_observation, "%s", OFF_runlist[off_run].first.c_str());
        string OFF_filename;
        OFF_filename = TString("$VERITAS_USER_DATA_DIR/"+TString(OFF_observation)+"_V6_Moderate-TMVA-BDT.RB."+TString(OFF_runnumber)+".root");
        OFF_pointing.push_back(GetRunElevAzim(OFF_filename,int(OFF_runlist[off_run].second)));
    }

    vector<pair<string,int>> new_list;
    for (int n=0;n<1;n++)
    {
        for (int on_run=0;on_run<ON_runlist.size();on_run++)
        {

            pair<string,int> best_match;
            pair<double,double> best_pointing;
            double best_chi2 = 10000.;
            for (int off_run=0;off_run<OFF_runlist.size();off_run++)
            {
                if (ON_runlist[on_run].first.compare(OFF_runlist[off_run].first) == 0) continue;
                bool already_used_run = false;
                for (int new_run=0;new_run<new_list.size();new_run++)
                {
                    if (int(new_list[new_run].second)==int(OFF_runlist[off_run].second)) already_used_run = true;
                }
                if (already_used_run) continue;
                double chi2 = 100.*pow(ON_pointing[on_run].first-OFF_pointing[off_run].first,2);
                chi2 += pow(ON_pointing[on_run].second-OFF_pointing[off_run].second,2);
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
                std::cout << "add run:" << std::endl;
                std::cout << best_match.first << " " << best_match.second << std::endl;
                std::cout << best_pointing.first << " " << best_pointing.second << std::endl;
            }
        }
    }
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
    return std::make_pair(func->Eval(energy),func->Eval(energy));
    //
    //TGraph *func = new TGraph(&Hist_Flux);
    //return std::make_pair(1.*func->Eval(energy),1.*func->Eval(energy));
}
bool DarkFoV() {
    //if (R2off<Theta2_cut_lower) return false;
    //if (R2off>Theta2_cut_upper) return false;
    if (theta2<0.1) return false;
    if (theta2>Theta2_cut_upper) return false;
    return true;
}
bool FoV() {
    //if (R2off<Theta2_cut_lower) return false;
    //if (R2off>Theta2_cut_upper) return false;
    if (theta2<Theta2_cut_lower) return false;
    if (theta2>Theta2_cut_upper) return false;
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
    if (MSCL>MSCL_cut_blind*1.0) return false;
    if (MSCW>MSCW_cut_blind*3.0) return false;
    return true;
}
void NetflixMethodGetShowerImage(string target_data, double tel_elev_lower_input, double tel_elev_upper_input, double theta2_cut_lower_input, double theta2_cut_upper_input)
{

    TH1::SetDefaultSumw2();
    sprintf(target, "%s", target_data.c_str());
    Theta2_cut_lower = theta2_cut_lower_input;
    Theta2_cut_upper = theta2_cut_upper_input;
    TelElev_lower = tel_elev_lower_input;
    TelElev_upper = tel_elev_upper_input;

    TH1D Hist_ErecS = TH1D("Hist_ErecS","",N_energy_bins,energy_bins);
    TH2D Hist_Dark_ShowerDirection = TH2D("Hist_Dark_ShowerDirection","",180,0,360,90,0,90);
    TH2D Hist_Data_ShowerDirection = TH2D("Hist_Data_ShowerDirection","",180,0,360,90,0,90);
    vector<TH2D> Hist_Data_MSCLW;
    vector<TH2D> Hist_Ring_MSCLW;
    vector<TH2D> Hist_Ring_Syst_MSCLW;
    vector<TH2D> Hist_Dark_MSCLW;
    vector<TH2D> Hist_Dark_Syst_MSCLW;
    vector<TH1D> Hist_Data_SR_SelectFoV_Theta2;
    vector<TH1D> Hist_Data_CR_SelectFoV_Theta2;
    vector<TH2D> Hist_Data_SR_Skymap;
    vector<TH2D> Hist_Data_CR_Skymap;
    vector<TH2D> Hist_Dark_SR_Skymap;
    vector<TH2D> Hist_Dark_CR_Skymap;
    for (int e=0;e<N_energy_bins;e++) 
    {
        char e_low[50];
        sprintf(e_low, "%i", int(energy_bins[e]));
        char e_up[50];
        sprintf(e_up, "%i", int(energy_bins[e+1]));
        Hist_Data_MSCLW.push_back(TH2D("Hist_Data_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_Ring_MSCLW.push_back(TH2D("Hist_Ring_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_Ring_Syst_MSCLW.push_back(TH2D("Hist_Ring_Syst_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_Dark_MSCLW.push_back(TH2D("Hist_Dark_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_Dark_Syst_MSCLW.push_back(TH2D("Hist_Dark_Syst_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_Data_SR_SelectFoV_Theta2.push_back(TH1D("Hist_Data_SR_SelectFoV_Theta2_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",1024,0,10));
        Hist_Data_CR_SelectFoV_Theta2.push_back(TH1D("Hist_Data_CR_SelectFoV_Theta2_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",1024,0,10));
        Hist_Data_SR_Skymap.push_back(TH2D("Hist_Data_SR_Skymap_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",150,-3,3,150,-3,3));
        Hist_Data_CR_Skymap.push_back(TH2D("Hist_Data_CR_Skymap_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",150,-3,3,150,-3,3));
        Hist_Dark_SR_Skymap.push_back(TH2D("Hist_Dark_SR_Skymap_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",30,-3,3,30,-3,3));
        Hist_Dark_CR_Skymap.push_back(TH2D("Hist_Dark_CR_Skymap_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",30,-3,3,30,-3,3));
    }

    // Get a list of target observation runs
    vector<pair<string,int>> Data_runlist_init = GetRunList(target);
    vector<pair<string,int>> Data_runlist;
    Data_runlist = SelectONRunList(Data_runlist_init,TelElev_lower,TelElev_upper,0,360);
    std::cout << "Data_runlist size = " << Data_runlist.size() << std::endl;
    if (Data_runlist.size()==0) return;

    // Get a list of dark observation runs
    vector<pair<string,int>> Dark_runlist_init = GetRunList("Everything");
    vector<pair<string,int>> Dark_runlist;
    std::cout << "initial Dark_runlist size = " << Dark_runlist_init.size() << std::endl;
    Dark_runlist = SelectOFFRunList(Data_runlist, Dark_runlist_init);
    std::cout << "final Dark_runlist size = " << Dark_runlist.size() << std::endl;

    for (int run=0;run<Dark_runlist.size();run++)
    {
        char run_number[50];
        char Dark_observation[50];
        sprintf(run_number, "%i", int(Dark_runlist[run].second));
        sprintf(Dark_observation, "%s", Dark_runlist[run].first.c_str());
        string filename;
        filename = TString("$VERITAS_USER_DATA_DIR/"+TString(Dark_observation)+"_V6_Moderate-TMVA-BDT.RB."+TString(run_number)+".root");


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
            int energy = Hist_ErecS.FindBin(ErecS*1000.)-1;
            if (energy<0) continue;
            if (energy>=N_energy_bins) continue;
            int e = energy;
            if (!SelectNImages(2,4)) continue;
            if (SizeSecondMax<600.) continue;
            if (pow(Xcore*Xcore+Ycore*Ycore,0.5)>350) continue;
            if (R2off>4.) continue;
            Hist_Dark_ShowerDirection.Fill(Shower_Az,Shower_Ze);
            if (DarkFoV())
            {
                Hist_Dark_MSCLW.at(e).Fill(MSCL,MSCW);
                Hist_Dark_Syst_MSCLW.at(e).Fill(MSCL,MSCW);
                if (SignalSelectionTheta2())
                {
                    Hist_Dark_SR_Skymap.at(e).Fill(Xoff,Yoff);
                }
                if (ControlSelectionTheta2())
                {
                    Hist_Dark_CR_Skymap.at(e).Fill(Xoff,Yoff);
                }
            }
        }
        input_file->Close();
    }



    for (int run=0;run<Data_runlist.size();run++)
    {
        char run_number[50];
        char Data_observation[50];
        sprintf(run_number, "%i", int(Data_runlist[run].second));
        sprintf(Data_observation, "%s", Data_runlist[run].first.c_str());
        string filename;
        filename = TString("$VERITAS_USER_DATA_DIR/"+TString(Data_observation)+"_V6_Moderate-TMVA-BDT.RB."+TString(run_number)+".root");


        TFile*  input_file = TFile::Open(filename.c_str());
	TH1* i_hEffAreaP = ( TH1* )getEffAreaHistogram(input_file,Data_runlist[run].second);
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
        exposure_hours += (time_1-time_0)/3600.;

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
            if (energy<0) continue;
            if (energy>=N_energy_bins) continue;
            int e = energy;
            if (!SelectNImages(2,4)) continue;
            if (SizeSecondMax<600.) continue;
            if (pow(Xcore*Xcore+Ycore*Ycore,0.5)>350) continue;
            if (R2off>4.) continue;
            Hist_Data_ShowerDirection.Fill(Shower_Az,Shower_Ze);
            if (FoV())
            {
                Hist_Data_MSCLW.at(e).Fill(MSCL,MSCW);
            }
            if (RingFoV())
            {
                Hist_Ring_MSCLW.at(e).Fill(MSCL,MSCW);
                Hist_Ring_Syst_MSCLW.at(e).Fill(MSCL,MSCW);
            }
            if (SignalSelectionTheta2())
            {
                if (FoV())
                {
                    Hist_Data_SR_SelectFoV_Theta2.at(e).Fill(theta2);
                    Hist_Data_SR_Skymap.at(e).Fill(-1.*ra_sky,dec_sky);
                }
            }
            if (ControlSelectionTheta2())
            {
                if (FoV())
                {
                    int binx = Hist_Dark_SR_Skymap.at(e).GetXaxis()->FindBin(Xoff);
                    int biny = Hist_Dark_SR_Skymap.at(e).GetYaxis()->FindBin(Yoff);
                    double dark_cr_content = Hist_Dark_CR_Skymap.at(e).GetBinContent(binx,biny);
                    double dark_sr_content = Hist_Dark_SR_Skymap.at(e).GetBinContent(binx,biny);
                    double weight = 0.;
                    if (dark_cr_content>0.) weight = dark_sr_content/dark_cr_content;
                    Hist_Data_CR_SelectFoV_Theta2.at(e).Fill(theta2,weight);
                    Hist_Data_CR_Skymap.at(e).Fill(-1.*ra_sky,dec_sky,weight);
                }
            }
        }
        input_file->Close();
    }


    for (int e=0;e<N_energy_bins;e++) 
    {
        int binx_blind = Hist_Ring_MSCLW.at(e).GetXaxis()->FindBin(MSCL_cut_lower);
        int binx_upper = Hist_Ring_MSCLW.at(e).GetXaxis()->FindBin(MSCL_cut_blind)-1;
        int biny_blind = Hist_Ring_MSCLW.at(e).GetYaxis()->FindBin(MSCW_cut_blind);
        int biny_upper = Hist_Ring_MSCLW.at(e).GetYaxis()->FindBin(MSCW_cut_blind*3)-1;
        double Ring_CR_Integral = Hist_Ring_MSCLW.at(e).Integral(binx_blind,binx_upper,biny_blind,biny_upper);
        double Data_CR_Integral = Hist_Data_MSCLW.at(e).Integral(binx_blind,binx_upper,biny_blind,biny_upper);
        double Ring_CR_Error = pow(Ring_CR_Integral,0.5);
        double Data_CR_Error = pow(Data_CR_Integral,0.5);
        double scale = Data_CR_Integral/Ring_CR_Integral;
        double scale_err = scale*pow(pow(Ring_CR_Error/Ring_CR_Integral,2)+pow(Data_CR_Error/Data_CR_Integral,2),0.5);
        Hist_Ring_MSCLW.at(e).Scale(scale);
        Hist_Ring_Syst_MSCLW.at(e).Scale(scale_err);
    }
    for (int e=0;e<N_energy_bins;e++) 
    {
        int binx_lower = Hist_Dark_MSCLW.at(e).GetXaxis()->FindBin(MSCL_cut_lower);
        int binx_blind = Hist_Dark_MSCLW.at(e).GetXaxis()->FindBin(MSCL_cut_blind)-1;
        int binx_upper = Hist_Dark_MSCLW.at(e).GetXaxis()->FindBin(MSCL_cut_blind*3)-1;
        int biny_lower = Hist_Dark_MSCLW.at(e).GetYaxis()->FindBin(MSCW_cut_lower);
        int biny_blind = Hist_Dark_MSCLW.at(e).GetYaxis()->FindBin(MSCW_cut_blind)-1;
        int biny_upper = Hist_Dark_MSCLW.at(e).GetYaxis()->FindBin(MSCW_cut_blind*3)-1;
        double Dark_SR_Integral = Hist_Dark_MSCLW.at(e).Integral(binx_lower,binx_blind,biny_lower,biny_blind);
        double Data_SR_Integral = Hist_Data_MSCLW.at(e).Integral(binx_lower,binx_blind,biny_lower,biny_blind);
        double Dark_Integral = Hist_Dark_MSCLW.at(e).Integral(binx_lower,binx_upper,biny_lower,biny_upper);
        double Data_Integral = Hist_Data_MSCLW.at(e).Integral(binx_lower,binx_upper,biny_lower,biny_upper);
        double Dark_CR_Integral = Dark_Integral-Dark_SR_Integral;
        double Data_CR_Integral = Data_Integral-Data_SR_Integral;
        double Dark_CR_Error = pow(Dark_CR_Integral,0.5);
        double Data_CR_Error = pow(Data_CR_Integral,0.5);
        double scale = Data_CR_Integral/Dark_CR_Integral;
        double scale_err = scale*pow(pow(Dark_CR_Error/Dark_CR_Integral,2)+pow(Data_CR_Error/Data_CR_Integral,2),0.5);
        Hist_Dark_MSCLW.at(e).Scale(scale);
        Hist_Dark_Syst_MSCLW.at(e).Scale(scale_err);
    }


    TFile OutputFile("output_Jul16/Netflix_"+TString(target)+"_TelElev"+std::to_string(int(TelElev_lower))+"to"+std::to_string(int(TelElev_upper))+"_Theta2"+std::to_string(int(10.*Theta2_cut_lower))+"to"+std::to_string(int(10.*Theta2_cut_upper))+".root","recreate");
    TTree InfoTree("InfoTree","info tree");
    InfoTree.Branch("exposure_hours",&exposure_hours,"exposure_hours/D");
    InfoTree.Fill();
    InfoTree.Write();
    Hist_Dark_ShowerDirection.Write();
    Hist_Data_ShowerDirection.Write();
    for (int e=0;e<N_energy_bins;e++)
    {
        Hist_Data_MSCLW.at(e).Write();
        Hist_Ring_MSCLW.at(e).Write();
        Hist_Ring_Syst_MSCLW.at(e).Write();
        Hist_Dark_MSCLW.at(e).Write();
        Hist_Dark_Syst_MSCLW.at(e).Write();
        Hist_Data_SR_SelectFoV_Theta2.at(e).Write();
        Hist_Data_CR_SelectFoV_Theta2.at(e).Write();
        Hist_Data_SR_Skymap.at(e).Write();
        Hist_Data_CR_Skymap.at(e).Write();
    }
    OutputFile.Close();

}
