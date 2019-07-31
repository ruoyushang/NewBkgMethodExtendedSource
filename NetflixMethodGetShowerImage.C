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

int N_bins_for_deconv = 10;
double MSCW_plot_lower = -1.;
double MSCW_plot_upper = 9.;
double MSCL_plot_lower = -1.;
double MSCL_plot_upper = 9.;
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
    //
    //TGraph *func = new TGraph(&Hist_Flux);
    //return std::make_pair(1.*func->Eval(energy),1.*func->Eval(energy));
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
    vector<TH2D> Hist_Data_MSCLW;
    vector<TH2D> Hist_Ring_MSCLW;
    vector<TH2D> Hist_Ring_Syst_MSCLW;
    vector<TH2D> Hist_Dark_MSCLW;
    vector<TH2D> Hist_Dark_Syst_MSCLW;
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
    }

    // Get a list of dark observation runs
    vector<pair<string,int>> Dark_runlist = GetRunList("Segue1V6");

    for (int run=0;run<Dark_runlist.size();run++)
    {
        char run_number[50];
        char Dark_observation[50];
        sprintf(run_number, "%i", int(Dark_runlist[run].second));
        sprintf(Dark_observation, "%s", Dark_runlist[run].first.c_str());
        string filename;
        filename = TString("$VERITAS_USER_DATA_DIR/"+TString(Dark_observation)+"_V6_Moderate-TMVA-BDT.RB."+TString(run_number)+".root");

        //if (!PointingSelection(filename,int(Dark_runlist[run].second),TelElev_lower,TelElev_upper,0,360)) continue;
        if (!PointingSelection(filename,int(Dark_runlist[run].second),70,85,0,360)) continue;

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
            //if (FoV())
            //{
                Hist_Dark_MSCLW.at(e).Fill(MSCL,MSCW);
                Hist_Dark_Syst_MSCLW.at(e).Fill(MSCL,MSCW);
            //}
        }
        input_file->Close();
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

        if (!PointingSelection(filename,int(Data_runlist[run].second),TelElev_lower,TelElev_upper,0,360)) continue;
        //if (!PointingSelection(filename,int(Data_runlist[run].second),55,70,0,360)) continue;

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
                Hist_Ring_Syst_MSCLW.at(e).Fill(MSCL,MSCW);
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
    for (int e=0;e<N_energy_bins;e++)
    {
        Hist_Data_MSCLW.at(e).Write();
        Hist_Ring_MSCLW.at(e).Write();
        Hist_Ring_Syst_MSCLW.at(e).Write();
        Hist_Dark_MSCLW.at(e).Write();
        Hist_Dark_Syst_MSCLW.at(e).Write();
    }
    OutputFile.Close();

}
