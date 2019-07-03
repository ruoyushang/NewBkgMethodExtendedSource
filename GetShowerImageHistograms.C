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

double MSCW_cut_lower = -1.0;
double MSCW_cut_blind = 1.0;
double MSCW_cut_upper = 1.0;
const int Number_of_SR = 6;
double MSCL_signal_cut_lower[Number_of_SR] = {0.75,0.50,0.25,0.00,-0.25,-0.50};
double MSCL_signal_cut_upper[Number_of_SR] = {1.00,0.75,0.50,0.25, 0.00,-0.25};
const int Number_of_CR = 6;
double MSCL_control_cut_lower[Number_of_CR] = {2.25,2.00,1.75,1.50,1.25,1.00};
double MSCL_control_cut_upper[Number_of_CR] = {2.50,2.25,2.00,1.75,1.50,1.25};
//const int Number_of_SR = 3;
//double MSCL_signal_cut_lower[Number_of_SR] = {0.50,0.00,-0.50};
//double MSCL_signal_cut_upper[Number_of_SR] = {1.00,0.50,0.00};
//const int Number_of_CR = 3;
//double MSCL_control_cut_lower[Number_of_CR] = {2.00,1.50,1.00};
//double MSCL_control_cut_upper[Number_of_CR] = {2.50,2.00,1.50};

const int N_energy_bins = 18;
double energy_bins[N_energy_bins+1] = {200,237,282,335,398,473,562,667,794,943,1122,1332,1585,1882,2239,3162,4467,6310,8913};
double electron_flux[N_energy_bins] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
double electron_flux_err[N_energy_bins] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
double electron_count[N_energy_bins] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
double electron_count_err[N_energy_bins] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
bool use_this_energy_bin[N_energy_bins] = {false,false,false,false,true,false,false,false,false,false,false,false,false,false,false,false,false,false};
//const int N_energy_bins = 1;
//double energy_bins[N_energy_bins+1] = {282,335};
//double electron_flux[N_energy_bins] = {0};
//double electron_flux_err[N_energy_bins] = {0};
//double electron_count[N_energy_bins] = {0};
//double electron_count_err[N_energy_bins] = {0};

int N_bins_for_deconv = 480;
double MSCW_plot_lower = -30.;
double MSCW_plot_upper = 30.;
double Theta2_cut_lower = 0;
double Theta2_cut_upper = 0;
double Theta2_upper_limit = 10;
char target[50] = "";

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
double exposure_hours = 0.;

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
}
bool FoV() {
    //if (R2off<Theta2_cut_lower) return false;
    //if (R2off>Theta2_cut_upper) return false;
    if (theta2<Theta2_cut_lower) return false;
    if (theta2>Theta2_cut_upper) return false;
    return true;
}
bool RingFoV() {
    //if (R2off<Theta2_cut_lower) return false;
    //if (R2off>Theta2_cut_upper) return false;
    if (theta2<Theta2_cut_upper) return false;
    return true;
}
bool SelectNImages(int Nmin, int Nmax)
{
    if (NImages<Nmin) return false;
    if (NImages>Nmax) return false;
    return true;
}
bool ControlSelectionTheta2()
{
    if (MSCL>MSCL_signal_cut_upper[0]) return false;
    if (MSCW<MSCW_cut_blind*1.0) return false;
    if (MSCW>MSCW_cut_blind*3.0) return false;
    return true;
}
bool SignalSelectionTheta2()
{
    if (MSCL>MSCL_signal_cut_upper[0]) return false;
    if (MSCW>MSCW_cut_blind) return false;
    if (MSCW<MSCW_cut_lower) return false;
    return true;
}
void GetShowerImageHistograms(string target_data, double theta2_cut_lower_input, double theta2_cut_upper_input)
{

    sprintf(target, "%s", target_data.c_str());
    Theta2_cut_lower = theta2_cut_lower_input;
    Theta2_cut_upper = theta2_cut_upper_input;

    TH1D Hist_ErecS = TH1D("Hist_ErecS","",N_energy_bins,energy_bins);
    TH1D Hist_EffAreaTime("Hist_EffAreaTime","",N_energy_bins,energy_bins);
    vector<TH2D> Hist_Data_MSCLW;
    vector<TH2D> Hist_Ring_MSCLW;
    vector<TH2D> Hist_Dark_MSCLW;
    vector<TH2D> Hist_CrabON_MSCLW;
    vector<TH2D> Hist_CrabOFF_MSCLW;
    vector<TH2D> Hist_Elec_MSCLW;
    vector<TH1D> Hist_Data_SR_FullFoV_Theta2;
    vector<TH1D> Hist_Data_SR_SelectFoV_Theta2;
    vector<TH1D> Hist_Data_CR_SelectFoV_Theta2;
    vector<TH2D> Hist_Data_SR_Skymap;
    vector<TH2D> Hist_Data_CR_Skymap;
    vector<TH1D> Hist_Dark_SR_FullFoV_Theta2;
    vector<TH1D> Hist_Dark_SR_SelectFoV_Theta2;
    vector<TH1D> Hist_Dark_CR_FullFoV_Theta2;
    for (int e=0;e<N_energy_bins;e++) 
    {
        char e_low[50];
        sprintf(e_low, "%i", int(energy_bins[e]));
        char e_up[50];
        sprintf(e_up, "%i", int(energy_bins[e+1]));
        if (energy_bins[e]>=200.) N_bins_for_deconv = 960;
        if (energy_bins[e]>=1000.) N_bins_for_deconv = 480;
        if (energy_bins[e]>=3200.) N_bins_for_deconv = 240;
        Hist_Data_MSCLW.push_back(TH2D("Hist_Data_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",Number_of_SR+Number_of_CR,MSCL_signal_cut_lower[Number_of_SR-1],MSCL_control_cut_upper[0],N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_Ring_MSCLW.push_back(TH2D("Hist_Ring_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",Number_of_SR+Number_of_CR,MSCL_signal_cut_lower[Number_of_SR-1],MSCL_control_cut_upper[0],N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_Dark_MSCLW.push_back(TH2D("Hist_Dark_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",Number_of_SR+Number_of_CR,MSCL_signal_cut_lower[Number_of_SR-1],MSCL_control_cut_upper[0],N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_CrabON_MSCLW.push_back(TH2D("Hist_CrabON_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",Number_of_SR+Number_of_CR,MSCL_signal_cut_lower[Number_of_SR-1],MSCL_control_cut_upper[0],N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_CrabOFF_MSCLW.push_back(TH2D("Hist_CrabOFF_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",Number_of_SR+Number_of_CR,MSCL_signal_cut_lower[Number_of_SR-1],MSCL_control_cut_upper[0],N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_Elec_MSCLW.push_back(TH2D("Hist_Elec_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",Number_of_SR+Number_of_CR,MSCL_signal_cut_lower[Number_of_SR-1],MSCL_control_cut_upper[0],N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_Data_SR_FullFoV_Theta2.push_back(TH1D("Hist_Data_SR_FullFoV_Theta2_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",1024,0,10));
        Hist_Data_SR_SelectFoV_Theta2.push_back(TH1D("Hist_Data_SR_SelectFoV_Theta2_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",1024,0,10));
        Hist_Data_CR_SelectFoV_Theta2.push_back(TH1D("Hist_Data_CR_SelectFoV_Theta2_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",1024,0,10));
        Hist_Data_SR_Skymap.push_back(TH2D("Hist_Data_SR_Skymap_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",100,-2,2,100,-2,2));
        Hist_Data_CR_Skymap.push_back(TH2D("Hist_Data_CR_Skymap_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",100,-2,2,100,-2,2));
        Hist_Dark_SR_FullFoV_Theta2.push_back(TH1D("Hist_Dark_SR_FullFoV_Theta2_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",1024,0,10));
        Hist_Dark_SR_SelectFoV_Theta2.push_back(TH1D("Hist_Dark_SR_SelectFoV_Theta2_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",1024,0,10));
        Hist_Dark_CR_FullFoV_Theta2.push_back(TH1D("Hist_Dark_CR_FullFoV_Theta2_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",1024,0,10));
    }

    // Get a list of dark observation runs for radial acceptance
    vector<pair<string,int>> Dark_runlist = GetRunList("Segue1V6");
    for (int run=0;run<Dark_runlist.size();run++)
    {
        char run_number[50];
        char Dark_observation[50];
        sprintf(run_number, "%i", int(Dark_runlist[run].second));
        sprintf(Dark_observation, "%s", Dark_runlist[run].first.c_str());
        string filename;
        filename = TString("$VERITAS_USER_DATA_DIR/"+TString(Dark_observation)+"_V6_Moderate-TMVA-BDT.RB."+TString(run_number)+".root");

        if (!PointingSelection(filename,int(Dark_runlist[run].second),60,90,0,360)) continue;

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
        Dark_tree->SetBranchAddress("Time",&Time);

        for (int entry=0;entry<Dark_tree->GetEntries();entry++) 
        {
            ErecS = 0;
            EChi2S = 0;
            NImages = 0;
            MSCW = 0;
            MSCL = 0;
            R2off = 0;
            Dark_tree->GetEntry(entry);
            R2off = Xoff*Xoff+Yoff*Yoff;
            int energy = Hist_ErecS.FindBin(ErecS*1000.)-1;
            if (energy<0) continue;
            if (energy>=N_energy_bins) continue;
            int e = energy;
            if (!SelectNImages(3,4)) continue;
            if (FoV())
            {
                Hist_Dark_MSCLW.at(e).Fill(MSCL,MSCW);
            }
            if (SignalSelectionTheta2())
            {
                Hist_Dark_SR_FullFoV_Theta2.at(e).Fill(R2off);
                if (FoV())
                {
                    Hist_Dark_SR_SelectFoV_Theta2.at(e).Fill(R2off);
                }
            }
            if (ControlSelectionTheta2())
            {
                Hist_Dark_CR_FullFoV_Theta2.at(e).Fill(R2off);
            }
        }
    }

    // Get a list of target observation runs
    vector<pair<string,int>> Data_runlist = GetRunList(target);
    for (int run=0;run<Data_runlist.size();run++)
    {
        char run_number[50];
        char Data_observation[50];
        sprintf(run_number, "%i", int(Data_runlist[run].second));
        sprintf(Data_observation, "%s", Data_runlist[run].first.c_str());
        string filename;
        filename = TString("$VERITAS_USER_DATA_DIR/"+TString(Data_observation)+"_V6_Moderate-TMVA-BDT.RB."+TString(run_number)+".root");

        if (!PointingSelection(filename,int(Data_runlist[run].second),60,90,0,360)) continue;

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
        Data_tree->SetBranchAddress("Time",&Time);

        // Get effective area and livetime and determine the cosmic electron counts for this run.
        Data_tree->GetEntry(0);
        double time_0 = Time;
        Data_tree->GetEntry(Data_tree->GetEntries()-1);
        double time_1 = Time;
        exposure_hours += (time_1-time_0)/3600.;
        for (int e=0;e<N_energy_bins;e++) 
        {
            double eff_area = i_hEffAreaP->GetBinContent( i_hEffAreaP->FindBin( log10(0.5*(energy_bins[e]+energy_bins[e+1])/1000.)));
            std::pair <double,double> mcgillflux = GetMcGillElectronFlux((energy_bins[e+1]+energy_bins[e])/2.);
            electron_flux[e] = mcgillflux.first;
            electron_flux_err[e] = mcgillflux.second;
            std::cout << "electron_flux[e] = " << electron_flux[e] << std::endl;
            double expected_electrons = 1e-12*10000.*electron_flux[e]*eff_area*(time_1-time_0)*(energy_bins[e+1]-energy_bins[e])/1000.;
            double expected_electrons_err = 1e-12*10000.*electron_flux_err[e]*eff_area*(time_1-time_0)*(energy_bins[e+1]-energy_bins[e])/1000.;
            std::cout << "expected_electrons = " << expected_electrons << std::endl;
            electron_count[e] += expected_electrons; // this is used to normalize MC electron template.
            electron_count_err[e] += expected_electrons_err; // correlated error, add in linear.
            Hist_EffAreaTime.Fill(energy_bins[e],eff_area*(time_1-time_0));
        }

        for (int entry=0;entry<Data_tree->GetEntries();entry++) 
        {
            ErecS = 0;
            EChi2S = 0;
            NImages = 0;
            MSCW = 0;
            MSCL = 0;
            R2off = 0;
            Data_tree->GetEntry(entry);
            R2off = Xoff*Xoff+Yoff*Yoff;
            int energy = Hist_ErecS.FindBin(ErecS*1000.)-1;
            if (energy<0) continue;
            if (energy>=N_energy_bins) continue;
            int e = energy;
            if (!SelectNImages(3,4)) continue;
            if (FoV())
            {
                Hist_Data_MSCLW.at(e).Fill(MSCL,MSCW);
            }
            if (RingFoV())
            {
                Hist_Ring_MSCLW.at(e).Fill(MSCL,MSCW);
            }
            if (SignalSelectionTheta2())
            {
                Hist_Data_SR_FullFoV_Theta2.at(e).Fill(theta2);
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
                    int bin = Hist_Dark_CR_FullFoV_Theta2.at(e).FindBin(R2off);
                    double dark_cr_content = Hist_Dark_CR_FullFoV_Theta2.at(e).GetBinContent(bin);
                    double dark_sr_content = Hist_Dark_SR_FullFoV_Theta2.at(e).GetBinContent(bin);
                    double weight = 0.;
                    if (dark_cr_content>0.) weight = dark_sr_content/dark_cr_content;
                    Hist_Data_CR_SelectFoV_Theta2.at(e).Fill(theta2,weight);
                    Hist_Data_CR_Skymap.at(e).Fill(-1.*ra_sky,dec_sky,weight);
                }
            }
        }
    }
    for (int e=0;e<N_energy_bins;e++) 
    {
        int binx_blind = Hist_Ring_MSCLW.at(e).GetXaxis()->FindBin(MSCL_signal_cut_lower[Number_of_SR-1]);
        int binx_upper = Hist_Ring_MSCLW.at(e).GetXaxis()->FindBin(MSCL_signal_cut_upper[0])-1;
        int biny_blind = Hist_Ring_MSCLW.at(e).GetYaxis()->FindBin(MSCW_cut_blind);
        int biny_upper = Hist_Ring_MSCLW.at(e).GetYaxis()->FindBin(MSCW_plot_upper)-1;
        double Ring_CR_Integral = Hist_Ring_MSCLW.at(e).Integral(binx_blind,binx_upper,biny_blind,biny_upper);
        double Data_CR_Integral = Hist_Data_MSCLW.at(e).Integral(binx_blind,binx_upper,biny_blind,biny_upper);
        double scale = Data_CR_Integral/Ring_CR_Integral;
        Hist_Ring_MSCLW.at(e).Scale(scale);
    }
    for (int e=0;e<N_energy_bins;e++) 
    {
        int binx_blind = Hist_Dark_MSCLW.at(e).GetXaxis()->FindBin(MSCL_signal_cut_lower[Number_of_SR-1]);
        int binx_upper = Hist_Dark_MSCLW.at(e).GetXaxis()->FindBin(MSCL_signal_cut_upper[0])-1;
        int biny_blind = Hist_Dark_MSCLW.at(e).GetYaxis()->FindBin(MSCW_cut_blind);
        int biny_upper = Hist_Dark_MSCLW.at(e).GetYaxis()->FindBin(MSCW_plot_upper)-1;
        double Dark_CR_Integral = Hist_Dark_MSCLW.at(e).Integral(binx_blind,binx_upper,biny_blind,biny_upper);
        double Data_CR_Integral = Hist_Data_MSCLW.at(e).Integral(binx_blind,binx_upper,biny_blind,biny_upper);
        double scale = Data_CR_Integral/Dark_CR_Integral;
        Hist_Dark_MSCLW.at(e).Scale(scale);
    }

    // Get e/gamma template from Crab
    vector<pair<string,int>> Crab_runlist = GetRunList("Crab");
    for (int run=0;run<Crab_runlist.size();run++)
    {
        char run_number[50];
        char Crab_observation[50];
        sprintf(run_number, "%i", int(Crab_runlist[run].second));
        sprintf(Crab_observation, "%s", Crab_runlist[run].first.c_str());
        string filename;
        filename = TString("$VERITAS_USER_DATA_DIR/"+TString(Crab_observation)+"_V6_Moderate-TMVA-BDT.RB."+TString(run_number)+".root");

        if (!PointingSelection(filename,int(Crab_runlist[run].second),60,90,0,360)) continue;

        TFile*  input_file = TFile::Open(filename.c_str());
        TH1* i_hEffAreaP = ( TH1* )getEffAreaHistogram(input_file,Crab_runlist[run].second);
        TString root_file = "run_"+TString(run_number)+"/stereo/data_on";
        TTree* Crab_tree = (TTree*) input_file->Get(root_file);
        Crab_tree->SetBranchAddress("Xoff",&Xoff);
        Crab_tree->SetBranchAddress("Yoff",&Yoff);
        Crab_tree->SetBranchAddress("theta2",&theta2);
        Crab_tree->SetBranchAddress("ra",&ra_sky);
        Crab_tree->SetBranchAddress("dec",&dec_sky);
        Crab_tree->SetBranchAddress("ErecS",&ErecS);
        Crab_tree->SetBranchAddress("EChi2S",&EChi2S);
        Crab_tree->SetBranchAddress("MSCW",&MSCW);
        Crab_tree->SetBranchAddress("MSCL",&MSCL);
        Crab_tree->SetBranchAddress("NImages",&NImages);
        Crab_tree->SetBranchAddress("Time",&Time);

        for (int entry=0;entry<Crab_tree->GetEntries();entry++) 
        {
            ErecS = 0;
            EChi2S = 0;
            NImages = 0;
            MSCW = 0;
            MSCL = 0;
            R2off = 0;
            Crab_tree->GetEntry(entry);
            R2off = Xoff*Xoff+Yoff*Yoff;
            int energy = Hist_ErecS.FindBin(ErecS*1000.)-1;
            if (energy<0) continue;
            if (energy>=N_energy_bins) continue;
            int e = energy;
            if (!SelectNImages(3,4)) continue;
            if (theta2<0.2)
            {
                Hist_CrabON_MSCLW.at(e).Fill(MSCL,MSCW);
            }
            else if (theta2>0.2 && theta2<0.4)
            {
                Hist_CrabOFF_MSCLW.at(e).Fill(MSCL,MSCW);
            }
        }
    }
    for (int e=0;e<N_energy_bins;e++) 
    {
        int binx_blind = Hist_CrabON_MSCLW.at(e).GetXaxis()->FindBin(MSCL_signal_cut_lower[Number_of_SR-1]);
        int binx_upper = Hist_CrabON_MSCLW.at(e).GetXaxis()->FindBin(MSCL_signal_cut_upper[0])-1;
        int biny_blind = Hist_CrabON_MSCLW.at(e).GetYaxis()->FindBin(MSCW_cut_blind);
        int biny_upper = Hist_CrabON_MSCLW.at(e).GetYaxis()->FindBin(MSCW_plot_upper)-1;
        double CrabON_CR_Integral = Hist_CrabON_MSCLW.at(e).Integral(binx_blind,binx_upper,biny_blind,biny_upper);
        double CrabOFF_CR_Integral = Hist_CrabOFF_MSCLW.at(e).Integral(binx_blind,binx_upper,biny_blind,biny_upper);
        double scale = CrabON_CR_Integral/CrabOFF_CR_Integral;
        Hist_CrabOFF_MSCLW.at(e).Scale(scale);
        Hist_Elec_MSCLW.at(e).Add(&Hist_CrabON_MSCLW.at(e));
        Hist_Elec_MSCLW.at(e).Add(&Hist_CrabOFF_MSCLW.at(e),-1.);
        for (int binx=0;binx<Hist_Elec_MSCLW.at(e).GetNbinsX();binx++)
        {
            for (int biny=0;biny<Hist_Elec_MSCLW.at(e).GetNbinsY();biny++)
            {
                double old_content = Hist_Elec_MSCLW.at(e).GetBinContent(binx+1,biny+1);
                double old_error = Hist_Elec_MSCLW.at(e).GetBinError(binx+1,biny+1);
                if (old_content<0)
                {
                    Hist_Elec_MSCLW.at(e).SetBinContent(binx+1,biny+1,0);
                    Hist_Elec_MSCLW.at(e).SetBinError(binx+1,biny+1,0);
                }
            }
        }
    }

    // Get a list of electron MC runs
    //vector<pair<string,int>> MC_runlist = GetRunList("Photon");
    //for (int run=0;run<MC_runlist.size();run++)
    //{
    //    char run_number[50];
    //    char MC_observation[50];
    //    sprintf(run_number, "%i", int(MC_runlist[run].second));
    //    sprintf(MC_observation, "%s", MC_runlist[run].first.c_str());
    //    string filename;
    //    filename = TString("$VERITAS_USER_DATA_DIR/"+TString(MC_observation)+"_V6_Moderate-TMVA-BDT.RB."+TString(run_number)+".root");
    //    TFile*  input_file = TFile::Open(filename.c_str());
    //    TString root_file = "run_"+TString(run_number)+"/stereo/data_on";
    //    TTree* MC_tree = (TTree*) input_file->Get(root_file);
    //    MC_tree->SetBranchAddress("Xoff",&Xoff);
    //    MC_tree->SetBranchAddress("Yoff",&Yoff);
    //    MC_tree->SetBranchAddress("theta2",&theta2);
    //    MC_tree->SetBranchAddress("ra",&ra_sky);
    //    MC_tree->SetBranchAddress("dec",&dec_sky);
    //    MC_tree->SetBranchAddress("ErecS",&ErecS);
    //    MC_tree->SetBranchAddress("EChi2S",&EChi2S);
    //    MC_tree->SetBranchAddress("MSCW",&MSCW);
    //    MC_tree->SetBranchAddress("MSCL",&MSCL);
    //    MC_tree->SetBranchAddress("NImages",&NImages);

    //    for (int entry=0;entry<MC_tree->GetEntries();entry++) 
    //    {
    //        ErecS = 0;
    //        EChi2S = 0;
    //        NImages = 0;
    //        MSCW = 0;
    //        MSCL = 0;
    //        R2off = 0;
    //        MC_tree->GetEntry(entry);
    //        R2off = Xoff*Xoff+Yoff*Yoff;
    //        //MSCW += 0.2;
    //        //MSCL += 0.2;
    //        int energy = Hist_ErecS.FindBin(ErecS*1000.)-1;
    //        if (energy<0) continue;
    //        if (energy>=N_energy_bins) continue;
    //        int e = energy;
    //        if (!SelectNImages(3,4)) continue;
    //        if (FoV())
    //        {
    //            Hist_Elec_MSCLW.at(e).Fill(MSCL,MSCW);
    //        }
    //    }
    //}
    
     
    // Normalize MC electron histograms with cosmic electron flux
    for (int e=0;e<N_energy_bins;e++) 
    {
        double old_integral = Hist_Elec_MSCLW.at(e).Integral();
        double scale = electron_count[e]/old_integral;
        //double scale_err = electron_count_err[e]/old_integral;
        double scale_err = 0.;
        double select_fov_integral = Hist_Dark_SR_SelectFoV_Theta2.at(e).Integral();
        double full_fov_integral = Hist_Dark_SR_FullFoV_Theta2.at(e).Integral();
        double scale_fov = select_fov_integral/full_fov_integral;
        for (int binx=0;binx<Hist_Elec_MSCLW.at(e).GetNbinsX();binx++)
        {
            for (int biny=0;biny<Hist_Elec_MSCLW.at(e).GetNbinsY();biny++)
            {
                double old_content = Hist_Elec_MSCLW.at(e).GetBinContent(binx+1,biny+1);
                double old_error = Hist_Elec_MSCLW.at(e).GetBinError(binx+1,biny+1);
                double new_content = old_content*scale*scale_fov;
                double new_error = old_error*old_error*scale*scale+(old_content*scale_err)*(old_content*scale_err);
                if (new_error>0) new_error = pow(new_error,0.5);
                Hist_Elec_MSCLW.at(e).SetBinContent(binx+1,biny+1,new_content);
                Hist_Elec_MSCLW.at(e).SetBinError(binx+1,biny+1,new_error);
            }
        }
    }


    TFile OutputFile("output_Jun15/Deconvolution_"+TString(target)+"_Theta2"+std::to_string(int(10.*Theta2_cut_lower))+"to"+std::to_string(int(10.*Theta2_cut_upper))+".root","recreate");
    TTree InfoTree("InfoTree","info tree");
    InfoTree.Branch("exposure_hours",&exposure_hours,"exposure_hours/D");
    InfoTree.Fill();
    InfoTree.Write();
    for (int e=0;e<N_energy_bins;e++)
    {
        Hist_Data_MSCLW.at(e).Write();
        Hist_Ring_MSCLW.at(e).Write();
        Hist_Dark_MSCLW.at(e).Write();
        Hist_Elec_MSCLW.at(e).Write();
        Hist_Data_SR_SelectFoV_Theta2.at(e).Write();
        Hist_Data_CR_SelectFoV_Theta2.at(e).Write();
        Hist_Data_SR_Skymap.at(e).Write();
        Hist_Data_CR_Skymap.at(e).Write();
    }
    OutputFile.Close();

}
