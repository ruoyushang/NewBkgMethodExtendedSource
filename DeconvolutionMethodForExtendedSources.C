
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
//char target[50] = "PKS1424";
//char target[50] = "3C264";
//char target[50] = "IC443";
//char target[50] = "Ton599";
//char target[50] = "2ndCrab";
//char target[50] = "BrandonValidation";

char Region[50] = "SR";
//char Region[50] = "VR";

double Elev_cut_lower = 55;
double Elev_cut_upper = 85;

double Target_Elev_cut_lower = 75;
double Target_Elev_cut_upper = 85;
//double Target_Elev_cut_lower = 55;
//double Target_Elev_cut_upper = 85;
//double Target_Elev_cut_lower = 35;
//double Target_Elev_cut_upper = 55;

double MSCW_cut_lower = -1.;
double MSCW_cut_upper = 1.;
double MSCL_cut_lower = -1.;
double MSCL_cut_upper = 1.;

double Depth_cut_lower = 6;
double Depth_cut_upper = 14;
//double Depth_cut_width = 3;
double Depth_cut_width = 4;

double Control_MSCL_cut = 0.;

double Norm_Lower = 1.0;
double Norm_Upper = 10.0;

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

//const int N_energy_bins = 1;
//double energy_bins[N_energy_bins+1] = {4000,1e6};
const int N_energy_bins = 11;
double energy_bins[N_energy_bins+1] = {150,200,250,300,400,600,800,1000,2000,4000,6000,10000};
//const int N_energy_bins = 23;
//double energy_bins[N_energy_bins+1] = {pow(10,2.1),pow(10,2.2),pow(10,2.3),pow(10,2.4),pow(10,2.5),pow(10,2.6),pow(10,2.7),pow(10,2.8),pow(10,2.9),pow(10,3.0),pow(10,3.1),pow(10,3.2),pow(10,3.3),pow(10,3.4),pow(10,3.5),pow(10,3.6),pow(10,3.7),pow(10,3.8),pow(10,3.9),pow(10,4.0),pow(10,4.1),pow(10,4.2),pow(10,4.3),pow(10,4.4)};
const int N_energy_bins_log = 23;
double energy_bins_log[N_energy_bins_log+1] = {pow(10,2.1),pow(10,2.2),pow(10,2.3),pow(10,2.4),pow(10,2.5),pow(10,2.6),pow(10,2.7),pow(10,2.8),pow(10,2.9),pow(10,3.0),pow(10,3.1),pow(10,3.2),pow(10,3.3),pow(10,3.4),pow(10,3.5),pow(10,3.6),pow(10,3.7),pow(10,3.8),pow(10,3.9),pow(10,4.0),pow(10,4.1),pow(10,4.2),pow(10,4.3),pow(10,4.4)};

int N_bins_for_deconv = 1200;

vector<int> GetRunList(string source) {
        vector<int> list;
        if (source=="BrandonValidation") {
                list.push_back(38030);
                list.push_back(38194);
                list.push_back(38418);
                list.push_back(38421);
                list.push_back(38538);
                list.push_back(38559);
                list.push_back(38622);
                list.push_back(38919);
                list.push_back(39032);
                list.push_back(39111);
                list.push_back(39165);
                list.push_back(39197);
                list.push_back(39200);
                list.push_back(39230);
                list.push_back(39258);
                list.push_back(39359);
                list.push_back(39588);
                list.push_back(39614);
                list.push_back(39642);
                list.push_back(39645);
                list.push_back(39702);
                list.push_back(39710);
                list.push_back(39796);
                list.push_back(39857);
                list.push_back(39877);
                list.push_back(39917);
                list.push_back(39921);
                list.push_back(39947);
                list.push_back(39950);
                list.push_back(40151);
                list.push_back(40225);
                list.push_back(40281);
                list.push_back(40314);
                list.push_back(40361);
                list.push_back(40391);
                list.push_back(40424);
                list.push_back(40555);
                list.push_back(40558);
                list.push_back(40593);
                list.push_back(40736);
                list.push_back(40760);
                list.push_back(40785);
                list.push_back(40844);
                list.push_back(40847);
                list.push_back(40941);
                list.push_back(40958);
                list.push_back(40976);
                list.push_back(40989);
                list.push_back(41011);
                list.push_back(41122);
                list.push_back(41125);
                list.push_back(43901);
                list.push_back(44049);
                list.push_back(44090);
                list.push_back(44494);
                list.push_back(44603);
                list.push_back(44692);
                list.push_back(45319);
                list.push_back(45322);
                list.push_back(45364);
        }
        if (source=="Ton599") {
                list.push_back(88356);
                list.push_back(88357);
                list.push_back(88378);
                list.push_back(88379);
                list.push_back(88380);
                list.push_back(88397);
                list.push_back(88420);
                list.push_back(88421);
                list.push_back(88422);
                list.push_back(88423);
                list.push_back(88447);
                list.push_back(88448);
                list.push_back(88449);
                list.push_back(88475);
        }
        if (source=="3C264") {
                list.push_back(89177);
                list.push_back(89178);
                list.push_back(89181);
                list.push_back(89200);
                list.push_back(89201);
                list.push_back(89202);
                list.push_back(89312);
                list.push_back(89531);
                list.push_back(89532);
                list.push_back(89534);
                list.push_back(89535);
                list.push_back(89537);
                list.push_back(89538);
                list.push_back(89554);
                list.push_back(89555);
                list.push_back(89557);
                list.push_back(89558);
                list.push_back(89559);
                list.push_back(89577);
                list.push_back(89590);
                list.push_back(89591);
                list.push_back(89592);
                list.push_back(89599);
                list.push_back(89600);
                list.push_back(89601);
                list.push_back(89648);
                list.push_back(89649);
                list.push_back(89650);
                list.push_back(89662);
                list.push_back(89684);
                list.push_back(89685);
                list.push_back(89692);
                list.push_back(89700);
                list.push_back(89701);
                list.push_back(89707);
                list.push_back(89708);
                list.push_back(89755);
                list.push_back(89756);
                list.push_back(89761);
                list.push_back(89762);
                list.push_back(89779);
                list.push_back(89780);
                list.push_back(89781);
                list.push_back(89782);
                list.push_back(89816);
                list.push_back(89817);
                list.push_back(89818);
                list.push_back(89819);
                list.push_back(89827);
                list.push_back(89828);
                list.push_back(89840);
                list.push_back(89841);
                list.push_back(89846);
                list.push_back(89847);
                list.push_back(89848);
                list.push_back(89878);
                list.push_back(89879);
                list.push_back(89880);
                list.push_back(89915);
                list.push_back(89916);
                list.push_back(89935);
                list.push_back(89936);
        }
        if (source=="PKS1424") {
                list.push_back(89252);
                list.push_back(89291);
                list.push_back(89342);
                list.push_back(89358);
                list.push_back(89359);
                list.push_back(89360);
                list.push_back(89560);
                list.push_back(89561);
                list.push_back(89562);
                list.push_back(89563);
                list.push_back(89764);
                list.push_back(89765);
                list.push_back(89945);
                list.push_back(89946);
                list.push_back(90082);
                list.push_back(90083);
                list.push_back(90228);
                list.push_back(90229);
                list.push_back(90230);
                list.push_back(90231);
                list.push_back(90232);
                list.push_back(90357);
                list.push_back(90432);
                list.push_back(90433);
        }
        if (source=="2ndCrab") {
                list.push_back(83006);
                list.push_back(83031);
                list.push_back(83032);
                list.push_back(83103);
                list.push_back(83157);
                list.push_back(83158);
                list.push_back(83168);
                list.push_back(83276);
                list.push_back(83277);
                list.push_back(83308);
                list.push_back(83310);
                list.push_back(83311);
                list.push_back(83359);
                list.push_back(83525);
                list.push_back(83626);
                list.push_back(83627);
                list.push_back(83628);
                list.push_back(83630);
                list.push_back(83653);
                list.push_back(83654);
                list.push_back(83655);
                list.push_back(83656);
                list.push_back(83678);
                list.push_back(83679);
                list.push_back(83974);
                list.push_back(83988);
                list.push_back(84087);
                list.push_back(84108);
                list.push_back(84109);
                list.push_back(84119);
                list.push_back(84121);
                list.push_back(84135);
                list.push_back(84174);
                list.push_back(84175);
                list.push_back(84287);
                list.push_back(84360);
                list.push_back(84377);
                list.push_back(84486);
                list.push_back(84587);
                list.push_back(84588);
                list.push_back(84594);
                list.push_back(84595);
                list.push_back(84620);
                list.push_back(84644);
                list.push_back(84685);
                list.push_back(84686);
                list.push_back(84824);
                list.push_back(84868);
                list.push_back(84896);
                list.push_back(84898);
        }
        if (source=="IC443") {
                list.push_back(70442);
                list.push_back(70443);
                list.push_back(70464);
                list.push_back(70491);
                list.push_back(70492);
                list.push_back(70536);
                list.push_back(70537);
        }
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
bool FoV() {
    if (theta2>0.2) return true;
    //if (theta2<1.0) return true;
    //if (theta2<0.5) return true;
    //if (theta2<0.2) return true;  // Crab signal icontamination with this cut is too strong for the deconvolution method.
    return false;
}
bool QualitySelection() {
    //if (EmissionHeightChi2/EmissionHeight<0.2) return true;
    return true;
}
bool SignalSelectionDepth() {
    if (MSCL<MSCL_cut_lower) return false;
    if (MSCL>MSCL_cut_upper) return false;
    if (MSCW<MSCW_cut_lower) return false;
    if (MSCW>MSCW_cut_upper) return false;
    return true;
}
bool SignalSelectionMSCW() {
    if (MSCL<MSCL_cut_lower) return false;
    if (MSCL>MSCL_cut_upper) return false;
    double cut_mean = 5.+log2(pow(ErecS*1000./0.08,0.4));
    Depth_cut_lower = cut_mean-Depth_cut_width;
    Depth_cut_upper = cut_mean+Depth_cut_width;
    if (SlantDepth*100./37.<Depth_cut_lower) return false;
    if (SlantDepth*100./37.>Depth_cut_upper) return false;
    return true;
}
bool AuxSignalSelectionMSCW() {
    if (MSCL<MSCL_cut_lower) return false;
    if (MSCL>MSCL_cut_upper) return false;
    double cut_mean = 5.+log2(pow(ErecS*1000./0.08,0.4));
    Depth_cut_lower = cut_mean-Depth_cut_width;
    Depth_cut_upper = cut_mean+Depth_cut_width;
    if (SlantDepth*100./37.>Depth_cut_lower && SlantDepth*100./37.<Depth_cut_upper) return false;
    if (SlantDepth*100./37.>Depth_cut_upper+1.) return false;
    if (SlantDepth*100./37.<Depth_cut_lower-3.) return false;
    return true;
}
bool ControlSelectionMSCW() {
    if (MSCL<MSCL_cut_upper) return false;
    double cut_mean = 5.+log2(pow(ErecS*1000./0.08,0.4));
    Depth_cut_lower = cut_mean-Depth_cut_width;
    Depth_cut_upper = cut_mean+Depth_cut_width;
    if (SlantDepth*100./37.<Depth_cut_lower) return false;
    if (SlantDepth*100./37.>Depth_cut_upper) return false;
    return true;
}
bool AuxControlSelectionMSCW() {
    if (MSCL<MSCL_cut_upper) return false;
    double cut_mean = 5.+log2(pow(ErecS*1000./0.08,0.4));
    Depth_cut_lower = cut_mean-Depth_cut_width;
    Depth_cut_upper = cut_mean+Depth_cut_width;
    if (SlantDepth*100./37.>Depth_cut_lower && SlantDepth*100./37.<Depth_cut_upper) return false;
    if (SlantDepth*100./37.>Depth_cut_upper+1.) return false;
    if (SlantDepth*100./37.<Depth_cut_lower-3.) return false;
    return true;
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
        //sp.Deconvolution(source,response,N_bins,25,1,10000); // new best option
        n_iteration = 10;
        //if (Hist_response->GetRMS()<=2.0) n_iteration = 100;
        //if (Hist_response->GetRMS()<=1.3) n_iteration = 1000;
        sp.Deconvolution(source,response,N_bins,n_iteration,1,100000); // new best option
        for (int i=0;i<N_bins;i++) {
            Hist_Deconv->SetBinContent(i+1,max(source[i],0.));
        }
}

void Convolution(TH1* Hist_source, TH1* Hist_response, TH1* Hist_Conv) {
        for (int i=0;i<N_bins_for_deconv;i++) {
            double conv = 0;
            for (int j=0;j<N_bins_for_deconv;j++) {
                if (i>j)
                    conv += Hist_source->GetBinContent(j+1)*Hist_response->GetBinContent(i-j+1);
                else
                    conv += Hist_source->GetBinContent(j+1)*Hist_response->GetBinContent(i-j+1+N_bins_for_deconv);
            }
            Hist_Conv->SetBinContent(i+1,max(conv,0.));
        }
}

void DeconvolutionMethodForExtendedSources(string target_data, double elev_lower, double elev_upper) {

        TH1::SetDefaultSumw2();
        sprintf(target, "%s", target_data.c_str());
        Target_Elev_cut_lower = elev_lower;
        Target_Elev_cut_upper = elev_upper;
        Elev_cut_lower = elev_lower;
        Elev_cut_upper = elev_upper;

        if (TString(Region)=="VR") {
                MSCW_cut_lower = 1.5;
                MSCW_cut_upper = 2.5;
                MSCL_cut_lower = 1.5;
                MSCL_cut_upper = 2.5;
        }

        TRandom rnd;
        TH2D Hist_Dark_TelElevAzim("Hist_Dark_TelElevAzim","",90,0,90,360,0,360);
        TH2D Hist_Target_TelElevAzim("Hist_Target_TelElevAzim","",90,0,90,360,0,360);
        TH2D Hist_Target_ON_ErecSDepth("Hist_Target_ON_ErecSDepth","",N_energy_bins_log,energy_bins_log,60,0,30);
        TH2D Hist_Target_OFF_ErecSDepth("Hist_Target_OFF_ErecSDepth","",N_energy_bins_log,energy_bins_log,60,0,30);
        TH1D Hist_Target_ON_MSCW_Alpha("Hist_Target_ON_MSCW_Alpha","",100,0,10);
        TH1D Hist_Target_OFF_MSCW_Alpha("Hist_Target_OFF_MSCW_Alpha","",100,0,10);
        vector<TH1D> Hist_Dark_SR_ErecS;
        vector<TH1D> Hist_Dark_SR_MSCW;
        vector<TH1D> Hist_Dark_CR_MSCW;
        vector<TH1D> Hist_Dark_ASR_MSCW;
        vector<TH1D> Hist_Dark_ACR_MSCW;
        vector<TH1D> Hist_Dark_BkgTemp_MSCW;
        vector<TH1D> Hist_Dark_Bkg_MSCW;
        vector<TH1D> Hist_Dark_ABkgTemp_MSCW;
        vector<TH1D> Hist_Dark_ABkg_MSCW;
        vector<TH1D> Hist_Dark_Deconv_MSCW;
        vector<TH1D> Hist_Dark_TrueDeconv_MSCW;
        vector<TH1D> Hist_Dark_Elec_MSCW;
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
        vector<TH1D> Hist_Target_TrueDeconv_MSCW;
        vector<TH1D> Hist_Target_Elec_MSCW;
        vector<TH1D> Hist_Target_Ring_MSCW;
        for (int e=0;e<N_energy_bins;e++) {
            char e_low[50];
            sprintf(e_low, "%i", int(energy_bins[e]));
            char e_up[50];
            sprintf(e_up, "%i", int(energy_bins[e+1]));
            N_bins_for_deconv = 1200;
            //if (energy_bins[e]>=400) N_bins_for_deconv = 1200;
            Hist_Dark_SR_ErecS.push_back(TH1D("Hist_Dark_SR_ErecS_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_energy_bins,energy_bins));
            Hist_Dark_SR_MSCW.push_back(TH1D("Hist_Dark_SR_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-200,100));
            Hist_Dark_CR_MSCW.push_back(TH1D("Hist_Dark_CR_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-200,100));
            Hist_Dark_ASR_MSCW.push_back(TH1D("Hist_Dark_ASR_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-200,100));
            Hist_Dark_ACR_MSCW.push_back(TH1D("Hist_Dark_ACR_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-200,100));
            Hist_Dark_BkgTemp_MSCW.push_back(TH1D("Hist_Dark_BkgTemp_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-200,100));
            Hist_Dark_Bkg_MSCW.push_back(TH1D("Hist_Dark_Bkg_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-200,100));
            Hist_Dark_ABkgTemp_MSCW.push_back(TH1D("Hist_Dark_ABkgTemp_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-200,100));
            Hist_Dark_ABkg_MSCW.push_back(TH1D("Hist_Dark_ABkg_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-200,100));
            Hist_Dark_Deconv_MSCW.push_back(TH1D("Hist_Dark_Deconv_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-200,100));
            Hist_Dark_TrueDeconv_MSCW.push_back(TH1D("Hist_Dark_TrueDeconv_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-200,100));
            Hist_Dark_Elec_MSCW.push_back(TH1D("Hist_Dark_Elec_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-200,100));
            Hist_Target_SR_ErecS.push_back(TH1D("Hist_Target_SR_ErecS_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_energy_bins,energy_bins));
            Hist_Target_SR_MSCW.push_back(TH1D("Hist_Target_SR_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-200,100));
            Hist_Target_CR_MSCW.push_back(TH1D("Hist_Target_CR_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-200,100));
            Hist_Target_ASR_MSCW.push_back(TH1D("Hist_Target_ASR_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-200,100));
            Hist_Target_ACR_MSCW.push_back(TH1D("Hist_Target_ACR_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-200,100));
            Hist_Target_BkgTemp_MSCW.push_back(TH1D("Hist_Target_BkgTemp_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-200,100));
            Hist_Target_Bkg_MSCW.push_back(TH1D("Hist_Target_Bkg_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-200,100));
            Hist_Target_ABkgTemp_MSCW.push_back(TH1D("Hist_Target_ABkgTemp_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-200,100));
            Hist_Target_ABkg_MSCW.push_back(TH1D("Hist_Target_ABkg_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-200,100));
            Hist_Target_Deconv_MSCW.push_back(TH1D("Hist_Target_Deconv_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-200,100));
            Hist_Target_TrueDeconv_MSCW.push_back(TH1D("Hist_Target_TrueDeconv_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-200,100));
            Hist_Target_Elec_MSCW.push_back(TH1D("Hist_Target_Elec_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-200,100));
            Hist_Target_Ring_MSCW.push_back(TH1D("Hist_Target_Ring_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-200,100));
        }

        vector<int> Dark_runlist = GetRunList("Crab");
        for (int run=0;run<Dark_runlist.size();run++) {
                char run_number[50];
                sprintf(run_number, "%i", int(Dark_runlist[run]));
                filename = TString("$VERITAS_USER_DATA_DIR/Crab_V6_Moderate-TMVA-BDT.RB."+TString(run_number)+".root");
                TFile*  input_file = TFile::Open(filename.c_str());
                //std::cout << "Open " << filename << std::endl;
                TTree* pointing_tree = (TTree*) input_file->Get("run_"+TString(run_number)+"/stereo/pointingDataReduced");
                pointing_tree->SetBranchAddress("TelElevation",&TelElevation);
                pointing_tree->SetBranchAddress("TelAzimuth",&TelAzimuth);
                pointing_tree->GetEntry(0);
                //std::cout << TelElevation << std::endl;
                if (TelElevation<Elev_cut_lower) continue;
                if (TelElevation>Elev_cut_upper) continue;
                for (int entry=0;entry<pointing_tree->GetEntries();entry++) {
                        pointing_tree->GetEntry(entry);
                        Hist_Dark_TelElevAzim.Fill(TelElevation,TelAzimuth);
                }

                TTree* Dark_tree = (TTree*) input_file->Get("run_"+TString(run_number)+"/stereo/data_off");
                Dark_tree->SetBranchAddress("ErecS",&ErecS);
                Dark_tree->SetBranchAddress("EChi2S",&EChi2S);
                Dark_tree->SetBranchAddress("MSCW",&MSCW);
                Dark_tree->SetBranchAddress("MSCL",&MSCL);
                Dark_tree->SetBranchAddress("SlantDepth",&SlantDepth);
                Dark_tree->SetBranchAddress("EmissionHeight",&EmissionHeight);
                Dark_tree->SetBranchAddress("EmissionHeightChi2",&EmissionHeightChi2);
                Dark_tree->SetBranchAddress("theta2",&theta2);
                for (int entry=0;entry<Dark_tree->GetEntries();entry++) {
                        Dark_tree->GetEntry(entry);
                        int e = Hist_Dark_SR_ErecS.at(0).FindBin(ErecS*1000.)-1;
                        if (e>=N_energy_bins) continue;
                        if (e<0) continue;
                        if (!FoV()) continue;
                        if (!QualitySelection()) continue;
                        if (SignalSelectionMSCW()) Hist_Dark_SR_MSCW.at(e).Fill(MSCW);
                        if (ControlSelectionMSCW()) Hist_Dark_CR_MSCW.at(e).Fill(MSCW);
                        if (AuxSignalSelectionMSCW()) Hist_Dark_ASR_MSCW.at(e).Fill(MSCW);
                        if (AuxControlSelectionMSCW()) Hist_Dark_ACR_MSCW.at(e).Fill(MSCW);
                }
                input_file->Close();

        }

        vector<int> Target_runlist = GetRunList(target);
        char observation[50];
        sprintf(observation, "%s", target);
        if (TString(target)=="2ndCrab") sprintf(observation, "%s", "Crab");
        for (int run=0;run<Target_runlist.size();run++) {
                char run_number[50];
                sprintf(run_number, "%i", int(Target_runlist[run]));
                filename = TString("$VERITAS_USER_DATA_DIR/"+TString(observation)+"_V6_Moderate-TMVA-BDT.RB."+TString(run_number)+".root");
                TFile*  input_file = TFile::Open(filename.c_str());
                //std::cout << "Open " << filename << std::endl;
                TTree* pointing_tree = (TTree*) input_file->Get("run_"+TString(run_number)+"/stereo/pointingDataReduced");
                pointing_tree->SetBranchAddress("TelElevation",&TelElevation);
                pointing_tree->SetBranchAddress("TelAzimuth",&TelAzimuth);
                pointing_tree->GetEntry(0);
                //std::cout << TelElevation << std::endl;
                if (TelElevation<Target_Elev_cut_lower) continue;
                if (TelElevation>Target_Elev_cut_upper) continue;
                for (int entry=0;entry<pointing_tree->GetEntries();entry++) {
                        pointing_tree->GetEntry(entry);
                        Hist_Target_TelElevAzim.Fill(TelElevation,TelAzimuth);
                }

                TString root_file = "run_"+TString(run_number)+"/stereo/data_on";
                if (TString(observation)=="BrandonValidation") {
                    root_file = "run_"+TString(run_number)+"/stereo/data_on";
                }
                TTree* Target_tree = (TTree*) input_file->Get(root_file);
                Target_tree->SetBranchAddress("ErecS",&ErecS);
                Target_tree->SetBranchAddress("EChi2S",&EChi2S);
                Target_tree->SetBranchAddress("MSCW",&MSCW);
                Target_tree->SetBranchAddress("MSCL",&MSCL);
                Target_tree->SetBranchAddress("SlantDepth",&SlantDepth);
                Target_tree->SetBranchAddress("EmissionHeight",&EmissionHeight);
                Target_tree->SetBranchAddress("EmissionHeightChi2",&EmissionHeightChi2);
                Target_tree->SetBranchAddress("theta2",&theta2);
                for (int entry=0;entry<Target_tree->GetEntries();entry++) {
                        Target_tree->GetEntry(entry);
                        int e = Hist_Target_SR_ErecS.at(0).FindBin(ErecS*1000.)-1;
                        if (e>=N_energy_bins) continue;
                        if (e<0) continue;
                        if (!QualitySelection()) continue;
                        if (FoV()) {
                                if (SignalSelectionMSCW()) Hist_Target_SR_MSCW.at(e).Fill(MSCW);
                                if (ControlSelectionMSCW()) Hist_Target_CR_MSCW.at(e).Fill(MSCW);
                                if (AuxSignalSelectionMSCW()) Hist_Target_ASR_MSCW.at(e).Fill(MSCW);
                                if (AuxControlSelectionMSCW()) Hist_Target_ACR_MSCW.at(e).Fill(MSCW);
                                if (SignalSelectionDepth()&&theta2<0.2) 
                                    Hist_Target_ON_ErecSDepth.Fill(ErecS*1000.,SlantDepth*100./37.);
                                if (SignalSelectionMSCW()&&theta2<0.2) 
                                    Hist_Target_ON_MSCW_Alpha.Fill(MSCW);
                                if (SignalSelectionDepth()&&theta2>0.2) 
                                    Hist_Target_OFF_ErecSDepth.Fill(ErecS*1000.,SlantDepth*100./37.);
                                if (SignalSelectionMSCW()&&theta2>0.2) 
                                    Hist_Target_OFF_MSCW_Alpha.Fill(MSCW);
                        }
                        else {
                                if (SignalSelectionMSCW()) Hist_Target_Ring_MSCW.at(e).Fill(MSCW);
                        }
                }
                input_file->Close();

        }

        vector<int> N_iter;
        for (int e=0;e<N_energy_bins;e++) {
                int n_iteration = 20;
                //if (energy_bins[e]>=400) n_iteration = 10;
                //if (energy_bins[e]>=1000) n_iteration = 15;
                //if (energy_bins[e]>=4000) n_iteration = 15;
                N_iter.push_back(n_iteration);
        }
        // Target deconvolution method
        for (int e=0;e<N_energy_bins;e++) {
                int n_iteration = N_iter.at(e);
                Deconvolution(&Hist_Target_ACR_MSCW.at(e),&Hist_Target_ASR_MSCW.at(e),&Hist_Target_Deconv_MSCW.at(e),n_iteration);
                Deconvolution(&Hist_Target_ACR_MSCW.at(e),&Hist_Target_Deconv_MSCW.at(e),&Hist_Target_ABkgTemp_MSCW.at(e),n_iteration);
                Deconvolution(&Hist_Target_CR_MSCW.at(e),&Hist_Target_Deconv_MSCW.at(e),&Hist_Target_BkgTemp_MSCW.at(e),n_iteration);
                Deconvolution(&Hist_Target_CR_MSCW.at(e),&Hist_Target_SR_MSCW.at(e),&Hist_Target_TrueDeconv_MSCW.at(e),n_iteration);
        }

        // Dark run method

        for (int e=0;e<N_energy_bins;e++) {
                int n_iteration = N_iter.at(e);
                Deconvolution(&Hist_Dark_ACR_MSCW.at(e), &Hist_Dark_ASR_MSCW.at(e), &Hist_Dark_Deconv_MSCW.at(e),n_iteration);
                Deconvolution(&Hist_Dark_ACR_MSCW.at(e), &Hist_Dark_Deconv_MSCW.at(e), &Hist_Dark_ABkgTemp_MSCW.at(e),n_iteration);
                Deconvolution(&Hist_Dark_CR_MSCW.at(e), &Hist_Dark_Deconv_MSCW.at(e), &Hist_Dark_BkgTemp_MSCW.at(e),n_iteration);
                Deconvolution(&Hist_Dark_CR_MSCW.at(e), &Hist_Dark_SR_MSCW.at(e), &Hist_Dark_TrueDeconv_MSCW.at(e),n_iteration);
        }

        for (int e=0;e<N_energy_bins;e++) {
                double offset_begin = Hist_Dark_ASR_MSCW.at(e).GetMean()-Hist_Dark_ABkgTemp_MSCW.at(e).GetMean();
                double offset_fit = 0;
                double scale_fit = 0;
                double chi2 = 0.;
                for (int fit=0;fit<20;fit++) {
                        double offset = offset_begin-1.+fit*0.1;
                        for (int i=0;i<Hist_Dark_SR_MSCW.at(e).GetNbinsX();i++) {
                                int b = Hist_Dark_Bkg_MSCW.at(e).FindBin(Hist_Dark_Bkg_MSCW.at(e).GetBinCenter(i+1)-offset);
                                Hist_Dark_ABkg_MSCW.at(e).SetBinContent(i+1,Hist_Dark_ABkgTemp_MSCW.at(e).GetBinContent(b));
                                Hist_Dark_ABkg_MSCW.at(e).SetBinError(i+1,Hist_Dark_ABkgTemp_MSCW.at(e).GetBinError(b));
                                Hist_Dark_Bkg_MSCW.at(e).SetBinContent(i+1,Hist_Dark_BkgTemp_MSCW.at(e).GetBinContent(b));
                                Hist_Dark_Bkg_MSCW.at(e).SetBinError(i+1,Hist_Dark_BkgTemp_MSCW.at(e).GetBinError(b));
                        }
                        int norm_bin_low_Dark = Hist_Dark_Bkg_MSCW.at(e).FindBin(Norm_Lower);
                        int norm_bin_up_Dark = Hist_Dark_Bkg_MSCW.at(e).FindBin(Norm_Upper);
                        double scale_Dark = Hist_Dark_SR_MSCW.at(e).Integral(norm_bin_low_Dark,norm_bin_up_Dark)/Hist_Dark_Bkg_MSCW.at(e).Integral(norm_bin_low_Dark,norm_bin_up_Dark);
                        Hist_Dark_Bkg_MSCW.at(e).Scale(scale_Dark);
                        double chi2_temp = 0;
                        for (int i=0;i<Hist_Dark_SR_MSCW.at(e).GetNbinsX();i++) {
                            if (Hist_Dark_Bkg_MSCW.at(e).GetBinCenter(i+1)<Norm_Lower) continue;
                            chi2_temp += pow(Hist_Dark_Bkg_MSCW.at(e).GetBinContent(i+1)-Hist_Dark_SR_MSCW.at(e).GetBinContent(i+1),2);
                        }
                        chi2_temp = 1./chi2_temp;
                        if (chi2<chi2_temp) {
                            chi2 = chi2_temp;
                            offset_fit = offset;
                            scale_fit = scale_Dark;
                        } 
                }
                for (int i=0;i<Hist_Dark_SR_MSCW.at(e).GetNbinsX();i++) {
                        int b = Hist_Dark_Bkg_MSCW.at(e).FindBin(Hist_Dark_Bkg_MSCW.at(e).GetBinCenter(i+1)-offset_fit);
                        Hist_Dark_ABkg_MSCW.at(e).SetBinContent(i+1,Hist_Dark_ABkgTemp_MSCW.at(e).GetBinContent(b));
                        Hist_Dark_ABkg_MSCW.at(e).SetBinError(i+1,Hist_Dark_ABkgTemp_MSCW.at(e).GetBinError(b));
                        Hist_Dark_Bkg_MSCW.at(e).SetBinContent(i+1,Hist_Dark_BkgTemp_MSCW.at(e).GetBinContent(b));
                        Hist_Dark_Bkg_MSCW.at(e).SetBinError(i+1,Hist_Dark_BkgTemp_MSCW.at(e).GetBinError(b));
                }
                Hist_Dark_Bkg_MSCW.at(e).Scale(scale_fit);
                Hist_Dark_ABkg_MSCW.at(e).Scale(scale_fit);
                Hist_Dark_CR_MSCW.at(e).Scale(scale_fit);
        }

        for (int e=0;e<N_energy_bins;e++) {
                double offset_begin = Hist_Target_ASR_MSCW.at(e).GetMean()-Hist_Target_ABkgTemp_MSCW.at(e).GetMean();
                double offset_fit = 0;
                double scale_fit = 0;
                double chi2 = 0.;
                for (int fit=0;fit<20;fit++) {
                        double offset = offset_begin-1.+fit*0.1;
                        for (int i=0;i<Hist_Target_SR_MSCW.at(e).GetNbinsX();i++) {
                                int b = Hist_Target_Bkg_MSCW.at(e).FindBin(Hist_Target_Bkg_MSCW.at(e).GetBinCenter(i+1)-offset);
                                Hist_Target_ABkg_MSCW.at(e).SetBinContent(i+1,Hist_Target_ABkgTemp_MSCW.at(e).GetBinContent(b));
                                Hist_Target_ABkg_MSCW.at(e).SetBinError(i+1,Hist_Target_ABkgTemp_MSCW.at(e).GetBinError(b));
                                Hist_Target_Bkg_MSCW.at(e).SetBinContent(i+1,Hist_Target_BkgTemp_MSCW.at(e).GetBinContent(b));
                                Hist_Target_Bkg_MSCW.at(e).SetBinError(i+1,Hist_Target_BkgTemp_MSCW.at(e).GetBinError(b));
                        }
                        int norm_bin_low_target = Hist_Target_Bkg_MSCW.at(e).FindBin(Norm_Lower);
                        int norm_bin_up_target = Hist_Target_Bkg_MSCW.at(e).FindBin(Norm_Upper);
                        double scale_target = Hist_Target_SR_MSCW.at(e).Integral(norm_bin_low_target,norm_bin_up_target)/Hist_Target_Bkg_MSCW.at(e).Integral(norm_bin_low_target,norm_bin_up_target);
                        Hist_Target_Bkg_MSCW.at(e).Scale(scale_target);
                        double chi2_temp = 0;
                        for (int i=0;i<Hist_Target_SR_MSCW.at(e).GetNbinsX();i++) {
                            if (Hist_Target_Bkg_MSCW.at(e).GetBinCenter(i+1)<Norm_Lower) continue;
                            chi2_temp += pow(Hist_Target_Bkg_MSCW.at(e).GetBinContent(i+1)-Hist_Target_SR_MSCW.at(e).GetBinContent(i+1),2);
                        }
                        chi2_temp = 1./chi2_temp;
                        if (chi2<chi2_temp) {
                            chi2 = chi2_temp;
                            offset_fit = offset;
                            scale_fit = scale_target;
                        } 
                }
                for (int i=0;i<Hist_Target_SR_MSCW.at(e).GetNbinsX();i++) {
                        int b = Hist_Target_Bkg_MSCW.at(e).FindBin(Hist_Target_Bkg_MSCW.at(e).GetBinCenter(i+1)-offset_fit);
                        Hist_Target_ABkg_MSCW.at(e).SetBinContent(i+1,Hist_Target_ABkgTemp_MSCW.at(e).GetBinContent(b));
                        Hist_Target_ABkg_MSCW.at(e).SetBinError(i+1,Hist_Target_ABkgTemp_MSCW.at(e).GetBinError(b));
                        Hist_Target_Bkg_MSCW.at(e).SetBinContent(i+1,Hist_Target_BkgTemp_MSCW.at(e).GetBinContent(b));
                        Hist_Target_Bkg_MSCW.at(e).SetBinError(i+1,Hist_Target_BkgTemp_MSCW.at(e).GetBinError(b));
                }
                Hist_Target_Bkg_MSCW.at(e).Scale(scale_fit);
                Hist_Target_ABkg_MSCW.at(e).Scale(scale_fit);
                Hist_Target_CR_MSCW.at(e).Scale(scale_fit);
        }

        // simple dark template method
        for (int e=0;e<N_energy_bins;e++) {
                int norm_bin_low_target = Hist_Target_Bkg_MSCW.at(e).FindBin(Norm_Lower);
                int norm_bin_up_target = Hist_Target_Bkg_MSCW.at(e).FindBin(Norm_Upper);
                double scale_target = Hist_Target_SR_MSCW.at(e).Integral(norm_bin_low_target,norm_bin_up_target)/Hist_Dark_SR_MSCW.at(e).Integral(norm_bin_low_target,norm_bin_up_target);
                Hist_Dark_SR_MSCW.at(e).Scale(scale_target);
                Hist_Dark_Bkg_MSCW.at(e).Scale(scale_target);
        }


        // Get electron bkg
        for (int e=0;e<N_energy_bins;e++) {

                Hist_Dark_Elec_MSCW.at(e).Reset();
                Hist_Dark_Elec_MSCW.at(e).Add(&Hist_Dark_SR_MSCW.at(e));
                Hist_Dark_Elec_MSCW.at(e).Add(&Hist_Dark_Bkg_MSCW.at(e),-1.);
                //for (int i=0;i<Hist_Target_SR_MSCW.at(e).GetNbinsX();i++) {
                //        if (Hist_Dark_SR_MSCW.at(e).GetBinContent(i+1)<Hist_Dark_Bkg_MSCW.at(e).GetBinContent(i+1)) 
                //            Hist_Dark_Elec_MSCW.at(e).SetBinContent(i+1,0.);
                //}
                Hist_Dark_Elec_MSCW.at(e).Divide(&Hist_Dark_Bkg_MSCW.at(e));
                Hist_Target_Elec_MSCW.at(e).Reset();
                Hist_Target_Elec_MSCW.at(e).Add(&Hist_Target_Bkg_MSCW.at(e));
                Hist_Target_Elec_MSCW.at(e).Multiply(&Hist_Dark_Elec_MSCW.at(e));
                Hist_Dark_Elec_MSCW.at(e).Multiply(&Hist_Dark_Bkg_MSCW.at(e));

        }

        // Get ring bkg
        for (int e=0;e<N_energy_bins;e++) {
                int norm_bin_low_ring = Hist_Target_Ring_MSCW.at(e).FindBin(Norm_Lower);
                int norm_bin_up_ring = Hist_Target_Ring_MSCW.at(e).FindBin(Norm_Upper);
                double scale_ring = Hist_Target_SR_MSCW.at(e).Integral(norm_bin_low_ring,norm_bin_up_ring)/Hist_Target_Ring_MSCW.at(e).Integral(norm_bin_low_ring,norm_bin_up_ring);
                if (!(scale_ring>0)) scale_ring = 0;
                if (!(Hist_Target_Ring_MSCW.at(e).Integral(norm_bin_low_ring,norm_bin_up_ring)>0)) scale_ring = 0;
                Hist_Target_Ring_MSCW.at(e).Scale(scale_ring);
        }

        int norm_bin_low_ring = Hist_Target_ON_MSCW_Alpha.FindBin(Norm_Lower);
        int norm_bin_up_ring = Hist_Target_ON_MSCW_Alpha.FindBin(Norm_Upper);
        double scale_ring = Hist_Target_ON_MSCW_Alpha.Integral(norm_bin_low_ring,norm_bin_up_ring)/Hist_Target_OFF_MSCW_Alpha.Integral(norm_bin_low_ring,norm_bin_up_ring);
        Hist_Target_OFF_ErecSDepth.Scale(scale_ring);
        Hist_Target_ON_ErecSDepth.Add(&Hist_Target_OFF_ErecSDepth,-1.);

        TFile OutputFile("output/Deconvolution_"+TString(target)+"_Elev"+std::to_string(int(Target_Elev_cut_lower))+"to"+std::to_string(int(Target_Elev_cut_upper))+"_"+TString(Region)+".root","recreate"); 
        TTree InfoTree("InfoTree","info tree");
        InfoTree.Branch("MSCW_cut_lower",&MSCW_cut_lower,"MSCW_cut_lower/D");
        InfoTree.Branch("MSCW_cut_upper",&MSCW_cut_upper,"MSCW_cut_upper/D");
        InfoTree.Branch("MSCL_cut_lower",&MSCL_cut_lower,"MSCW_cut_lower/D");
        InfoTree.Branch("MSCL_cut_upper",&MSCL_cut_upper,"MSCW_cut_upper/D");
        InfoTree.Branch("Depth_cut_width",&Depth_cut_width,"Depth_cut_width/D");
        InfoTree.Branch("Elev_cut_lower",&Elev_cut_lower,"Elev_cut_lower/D");
        InfoTree.Branch("Elev_cut_upper",&Elev_cut_upper,"Elev_cut_upper/D");
        InfoTree.Branch("Target_Elev_cut_lower",&Target_Elev_cut_lower,"Target_Elev_cut_lower/D");
        InfoTree.Branch("Target_Elev_cut_upper",&Target_Elev_cut_upper,"Target_Elev_cut_upper/D");
        InfoTree.Fill();
        InfoTree.Write();
        Hist_Target_TelElevAzim.Write();
        Hist_Target_ON_ErecSDepth.Write();
        Hist_Target_OFF_ErecSDepth.Write();
        Hist_Dark_TelElevAzim.Write();
        for (int e=0;e<N_energy_bins;e++) {
                Hist_Dark_ASR_MSCW.at(e).Write();
                Hist_Dark_ACR_MSCW.at(e).Write();
                Hist_Dark_ABkg_MSCW.at(e).Write();
                Hist_Dark_SR_MSCW.at(e).Write();
                Hist_Dark_CR_MSCW.at(e).Write();
                Hist_Dark_Bkg_MSCW.at(e).Write();
                Hist_Dark_Deconv_MSCW.at(e).Write();
                Hist_Dark_TrueDeconv_MSCW.at(e).Write();
                Hist_Dark_Elec_MSCW.at(e).Write();
                Hist_Target_ASR_MSCW.at(e).Write();
                Hist_Target_ACR_MSCW.at(e).Write();
                Hist_Target_ABkg_MSCW.at(e).Write();
                Hist_Target_SR_MSCW.at(e).Write();
                Hist_Target_CR_MSCW.at(e).Write();
                Hist_Target_Bkg_MSCW.at(e).Write();
                Hist_Target_Deconv_MSCW.at(e).Write();
                Hist_Target_TrueDeconv_MSCW.at(e).Write();
                Hist_Target_Elec_MSCW.at(e).Write();
                Hist_Target_Ring_MSCW.at(e).Write();
        }
        OutputFile.Close();

}
