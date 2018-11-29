
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
char Region[50] = "SR";

double Elev_cut_lower = 75;
double Elev_cut_upper = 85;

double Target_Elev_cut_lower = 75;
double Target_Elev_cut_upper = 85;
//double Target_Elev_cut_lower = 55;
//double Target_Elev_cut_upper = 85;
//double Target_Elev_cut_lower = 35;
//double Target_Elev_cut_upper = 55;

double MSCW_cut_lower = -0.5;
double MSCW_cut_upper = 0.7;
double MSCL_cut_lower = -0.5;
double MSCL_cut_upper = 0.7;

double Depth_cut_lower = 6;
double Depth_cut_upper = 14;
double Depth_cut_width = 4.;

string  filename;
double TelElevation = 0;
double TelAzimuth = 0;
double ErecS = 0;
double EChi2S = 0;
double MSCW = 0;
double MSCL = 0;
double SlantDepth = 0;
double theta2 = 0;

//const int N_energy_bins = 1;
//double energy_bins[N_energy_bins+1] = {4000,1e10};
//const int N_energy_bins = 10;
//double energy_bins[N_energy_bins+1] = {pow(10,2.2),pow(10,2.4),pow(10,2.6),pow(10,2.8),pow(10,3.0),pow(10,3.2),pow(10,3.4),pow(10,3.6),pow(10,3.8),pow(10,4.0),pow(10,4.2)};
const int N_energy_bins = 23;
double energy_bins[N_energy_bins+1] = {pow(10,2.1),pow(10,2.2),pow(10,2.3),pow(10,2.4),pow(10,2.5),pow(10,2.6),pow(10,2.7),pow(10,2.8),pow(10,2.9),pow(10,3.0),pow(10,3.1),pow(10,3.2),pow(10,3.3),pow(10,3.4),pow(10,3.5),pow(10,3.6),pow(10,3.7),pow(10,3.8),pow(10,3.9),pow(10,4.0),pow(10,4.1),pow(10,4.2),pow(10,4.3),pow(10,4.4)};

//const int N_bins_for_deconv = int(pow(2,7));
//const int N_bins_for_deconv = 600;
const int N_bins_for_deconv = 2400;
//const int N_bins_for_deconv = 4800;

vector<int> GetRunList(string source) {
        vector<int> list;
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
    //if (theta2<0.4) return true;
    return false;
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
    return true;
}
bool ControlSelectionMSCW() {
    if (MSCL<MSCL_cut_upper*1.5) return false;
    if (MSCL>MSCL_cut_upper*15.) return false;
    double cut_mean = 5.+log2(pow(ErecS*1000./0.08,0.4));
    Depth_cut_lower = cut_mean-Depth_cut_width;
    Depth_cut_upper = cut_mean+Depth_cut_width;
    if (SlantDepth*100./37.<Depth_cut_lower) return false;
    if (SlantDepth*100./37.>Depth_cut_upper) return false;
    return true;
}
bool AuxControlSelectionMSCW() {
    if (MSCL<MSCL_cut_upper*1.5) return false;
    if (MSCL>MSCL_cut_upper*15.) return false;
    double cut_mean = 5.+log2(pow(ErecS*1000./0.08,0.4));
    Depth_cut_lower = cut_mean-Depth_cut_width;
    Depth_cut_upper = cut_mean+Depth_cut_width;
    if (SlantDepth*100./37.>Depth_cut_lower && SlantDepth*100./37.<Depth_cut_upper) return false;
    return true;
    //if (ControlSelectionMSCW() || AuxSignalSelectionMSCW()) return true;
    return false;
}
bool SignalSelectionMSCL() {
    if (MSCW<MSCW_cut_lower) return false;
    if (MSCW>MSCW_cut_upper) return false;
    double cut_mean = 5.+log2(pow(ErecS*1000./0.08,0.4));
    Depth_cut_lower = cut_mean-Depth_cut_width;
    Depth_cut_upper = cut_mean+Depth_cut_width;
    if (SlantDepth*100./37.<Depth_cut_lower) return false;
    if (SlantDepth*100./37.>Depth_cut_upper) return false;
    return true;
}
bool AuxSignalSelectionMSCL() {
    if (MSCW<MSCW_cut_lower) return false;
    if (MSCW>MSCW_cut_upper) return false;
    double cut_mean = 5.+log2(pow(ErecS*1000./0.08,0.4));
    Depth_cut_lower = cut_mean-Depth_cut_width;
    Depth_cut_upper = cut_mean+Depth_cut_width;
    if (SlantDepth*100./37.>Depth_cut_lower && SlantDepth*100./37.<Depth_cut_upper) return false;
    return true;
}
bool ControlSelectionMSCL() {
    if (MSCW<MSCW_cut_upper*1.5) return false;
    if (MSCW>MSCW_cut_upper*15.) return false;
    double cut_mean = 5.+log2(pow(ErecS*1000./0.08,0.4));
    Depth_cut_lower = cut_mean-Depth_cut_width;
    Depth_cut_upper = cut_mean+Depth_cut_width;
    if (SlantDepth*100./37.<Depth_cut_lower) return false;
    if (SlantDepth*100./37.>Depth_cut_upper) return false;
    return true;
}
bool AuxControlSelectionMSCL() {
    if (MSCW<MSCW_cut_upper*1.5) return false;
    if (MSCW>MSCW_cut_upper*15.) return false;
    double cut_mean = 5.+log2(pow(ErecS*1000./0.08,0.4));
    Depth_cut_lower = cut_mean-Depth_cut_width;
    Depth_cut_upper = cut_mean+Depth_cut_width;
    if (SlantDepth*100./37.>Depth_cut_lower && SlantDepth*100./37.<Depth_cut_upper) return false;
    return true;
    //if (ControlSelectionMSCL() || AuxSignalSelectionMSCL()) return true;
    return false;
}
bool SignalSelectionDepth() {
    if (MSCW<MSCW_cut_lower) return false;
    if (MSCW>MSCW_cut_upper) return false;
    if (MSCL<MSCL_cut_lower) return false;
    if (MSCL>MSCL_cut_upper) return false;
    return true;
}
bool AuxSignalSelectionDepth() {
    if (MSCW<MSCW_cut_upper*1.5) return false;
    if (MSCW>MSCW_cut_upper*2.5) return false;
    return true;
}
bool ControlSelectionDepth() {
    if (MSCL<MSCL_cut_upper*1.5) return false;
    if (MSCL>MSCL_cut_upper*2.5) return false;
    return true;
}
bool AuxControlSelectionDepth() {
    //if (MSCL<MSCL_cut_upper*1.5) return false;
    //if (MSCL>MSCL_cut_upper*2.5) return false;
    //if (MSCW<MSCW_cut_upper*1.5) return false;
    //if (MSCW>MSCW_cut_upper*2.5) return false;
    //return true;
    if (ControlSelectionDepth() || AuxSignalSelectionDepth()) return true;
    return false;
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
        //sp.Deconvolution(source,response,N_bins_for_deconv,100,1,100);  // best option
        //sp.Deconvolution(source,response,N_bins_for_deconv,100,1,1000);
        sp.Deconvolution(source,response,N_bins_for_deconv,100,1,10000); // new best option
        for (int i=0;i<N_bins_for_deconv;i++) {
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

void RatioMethodForExtendedSources() {

        TH1::SetDefaultSumw2();

        if (TString(Region)=="VR") {
                MSCW_cut_lower = 1.4;
                MSCW_cut_upper = 1.6;
                MSCL_cut_lower = 1.4;
                MSCL_cut_upper = 1.6;
        }

        TRandom rnd;
        TH2D Hist_Dark_TelElevAzim("Hist_Dark_TelElevAzim","",90,0,90,360,0,360);
        TH2D Hist_Target_TelElevAzim("Hist_Target_TelElevAzim","",90,0,90,360,0,360);
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
        vector<TH1D> Hist_Dark_SR_MSCL;
        vector<TH1D> Hist_Dark_CR_MSCL;
        vector<TH1D> Hist_Dark_ASR_MSCL;
        vector<TH1D> Hist_Dark_ACR_MSCL;
        vector<TH1D> Hist_Dark_BkgTemp_MSCL;
        vector<TH1D> Hist_Dark_Bkg_MSCL;
        vector<TH1D> Hist_Dark_ABkgTemp_MSCL;
        vector<TH1D> Hist_Dark_ABkg_MSCL;
        vector<TH1D> Hist_Dark_Deconv_MSCL;
        vector<TH1D> Hist_Dark_SR_Depth;
        vector<TH1D> Hist_Dark_CR_Depth;
        vector<TH1D> Hist_Dark_ASR_Depth;
        vector<TH1D> Hist_Dark_ACR_Depth;
        vector<TH1D> Hist_Dark_BkgTemp_Depth;
        vector<TH1D> Hist_Dark_Bkg_Depth;
        vector<TH1D> Hist_Dark_ABkgTemp_Depth;
        vector<TH1D> Hist_Dark_ABkg_Depth;
        vector<TH1D> Hist_Dark_Deconv_Depth;
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
        vector<TH1D> Hist_Target_SR_MSCL;
        vector<TH1D> Hist_Target_CR_MSCL;
        vector<TH1D> Hist_Target_ASR_MSCL;
        vector<TH1D> Hist_Target_ACR_MSCL;
        vector<TH1D> Hist_Target_BkgTemp_MSCL;
        vector<TH1D> Hist_Target_Bkg_MSCL;
        vector<TH1D> Hist_Target_ABkgTemp_MSCL;
        vector<TH1D> Hist_Target_ABkg_MSCL;
        vector<TH1D> Hist_Target_Deconv_MSCL;
        vector<TH1D> Hist_Target_SR_Depth;
        vector<TH1D> Hist_Target_CR_Depth;
        vector<TH1D> Hist_Target_ASR_Depth;
        vector<TH1D> Hist_Target_ACR_Depth;
        vector<TH1D> Hist_Target_BkgTemp_Depth;
        vector<TH1D> Hist_Target_Bkg_Depth;
        vector<TH1D> Hist_Target_ABkgTemp_Depth;
        vector<TH1D> Hist_Target_ABkg_Depth;
        vector<TH1D> Hist_Target_Deconv_Depth;
        for (int e=0;e<N_energy_bins;e++) {
            char e_low[50];
            sprintf(e_low, "%i", int(energy_bins[e]));
            char e_up[50];
            sprintf(e_up, "%i", int(energy_bins[e+1]));
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
            Hist_Dark_SR_MSCL.push_back(TH1D("Hist_Dark_SR_MSCL_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-200,100));
            Hist_Dark_CR_MSCL.push_back(TH1D("Hist_Dark_CR_MSCL_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-200,100));
            Hist_Dark_ASR_MSCL.push_back(TH1D("Hist_Dark_ASR_MSCL_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-200,100));
            Hist_Dark_ACR_MSCL.push_back(TH1D("Hist_Dark_ACR_MSCL_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-200,100));
            Hist_Dark_BkgTemp_MSCL.push_back(TH1D("Hist_Dark_BkgTemp_MSCL_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-200,100));
            Hist_Dark_Bkg_MSCL.push_back(TH1D("Hist_Dark_Bkg_MSCL_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-200,100));
            Hist_Dark_ABkgTemp_MSCL.push_back(TH1D("Hist_Dark_ABkgTemp_MSCL_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-200,100));
            Hist_Dark_ABkg_MSCL.push_back(TH1D("Hist_Dark_ABkg_MSCL_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-200,100));
            Hist_Dark_Deconv_MSCL.push_back(TH1D("Hist_Dark_Deconv_MSCL_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-200,100));
            Hist_Dark_SR_Depth.push_back(TH1D("Hist_Dark_SR_Depth_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-200,100));
            Hist_Dark_CR_Depth.push_back(TH1D("Hist_Dark_CR_Depth_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-200,100));
            Hist_Dark_ASR_Depth.push_back(TH1D("Hist_Dark_ASR_Depth_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-200,100));
            Hist_Dark_ACR_Depth.push_back(TH1D("Hist_Dark_ACR_Depth_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-200,100));
            Hist_Dark_BkgTemp_Depth.push_back(TH1D("Hist_Dark_BkgTemp_Depth_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-200,100));
            Hist_Dark_Bkg_Depth.push_back(TH1D("Hist_Dark_Bkg_Depth_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-200,100));
            Hist_Dark_ABkgTemp_Depth.push_back(TH1D("Hist_Dark_ABkgTemp_Depth_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-200,100));
            Hist_Dark_ABkg_Depth.push_back(TH1D("Hist_Dark_ABkg_Depth_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-200,100));
            Hist_Dark_Deconv_Depth.push_back(TH1D("Hist_Dark_Deconv_Depth_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-200,100));
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
            Hist_Target_SR_MSCL.push_back(TH1D("Hist_Target_SR_MSCL_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-200,100));
            Hist_Target_CR_MSCL.push_back(TH1D("Hist_Target_CR_MSCL_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-200,100));
            Hist_Target_ASR_MSCL.push_back(TH1D("Hist_Target_ASR_MSCL_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-200,100));
            Hist_Target_ACR_MSCL.push_back(TH1D("Hist_Target_ACR_MSCL_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-200,100));
            Hist_Target_BkgTemp_MSCL.push_back(TH1D("Hist_Target_BkgTemp_MSCL_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-200,100));
            Hist_Target_Bkg_MSCL.push_back(TH1D("Hist_Target_Bkg_MSCL_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-200,100));
            Hist_Target_ABkgTemp_MSCL.push_back(TH1D("Hist_Target_ABkgTemp_MSCL_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-200,100));
            Hist_Target_ABkg_MSCL.push_back(TH1D("Hist_Target_ABkg_MSCL_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-200,100));
            Hist_Target_Deconv_MSCL.push_back(TH1D("Hist_Target_Deconv_MSCL_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-200,100));
            Hist_Target_SR_Depth.push_back(TH1D("Hist_Target_SR_Depth_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-200,100));
            Hist_Target_CR_Depth.push_back(TH1D("Hist_Target_CR_Depth_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-200,100));
            Hist_Target_ASR_Depth.push_back(TH1D("Hist_Target_ASR_Depth_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-200,100));
            Hist_Target_ACR_Depth.push_back(TH1D("Hist_Target_ACR_Depth_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-200,100));
            Hist_Target_BkgTemp_Depth.push_back(TH1D("Hist_Target_BkgTemp_Depth_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-200,100));
            Hist_Target_Bkg_Depth.push_back(TH1D("Hist_Target_Bkg_Depth_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-200,100));
            Hist_Target_ABkgTemp_Depth.push_back(TH1D("Hist_Target_ABkgTemp_Depth_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-200,100));
            Hist_Target_ABkg_Depth.push_back(TH1D("Hist_Target_ABkg_Depth_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-200,100));
            Hist_Target_Deconv_Depth.push_back(TH1D("Hist_Target_Deconv_Depth_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-200,100));
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
                Dark_tree->SetBranchAddress("theta2",&theta2);
                for (int entry=0;entry<Dark_tree->GetEntries();entry++) {
                        Dark_tree->GetEntry(entry);
                        int e = Hist_Dark_SR_ErecS.at(0).FindBin(ErecS*1000.)-1;
                        if (e>=N_energy_bins) continue;
                        if (e<0) continue;
                        if (!FoV()) continue;
                        if (SignalSelectionMSCW()) Hist_Dark_SR_MSCW.at(e).Fill(MSCW);
                        if (ControlSelectionMSCW()) Hist_Dark_CR_MSCW.at(e).Fill(MSCW);
                        if (AuxSignalSelectionMSCW()) Hist_Dark_ASR_MSCW.at(e).Fill(MSCW);
                        if (AuxControlSelectionMSCW()) Hist_Dark_ACR_MSCW.at(e).Fill(MSCW);
                        if (SignalSelectionMSCL()) Hist_Dark_SR_MSCL.at(e).Fill(MSCL);
                        if (ControlSelectionMSCL()) Hist_Dark_CR_MSCL.at(e).Fill(MSCL);
                        if (AuxSignalSelectionMSCL()) Hist_Dark_ASR_MSCL.at(e).Fill(MSCL);
                        if (AuxControlSelectionMSCL()) Hist_Dark_ACR_MSCL.at(e).Fill(MSCL);
                        if (SignalSelectionDepth()) Hist_Dark_SR_Depth.at(e).Fill(SlantDepth*100./37.);
                        if (ControlSelectionDepth()) Hist_Dark_CR_Depth.at(e).Fill(SlantDepth*100./37.);
                        if (AuxSignalSelectionDepth()) Hist_Dark_ASR_Depth.at(e).Fill(SlantDepth*100./37.);
                        if (AuxControlSelectionDepth()) Hist_Dark_ACR_Depth.at(e).Fill(SlantDepth*100./37.);
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

                TTree* Target_tree = (TTree*) input_file->Get("run_"+TString(run_number)+"/stereo/data_on");
                Target_tree->SetBranchAddress("ErecS",&ErecS);
                Target_tree->SetBranchAddress("EChi2S",&EChi2S);
                Target_tree->SetBranchAddress("MSCW",&MSCW);
                Target_tree->SetBranchAddress("MSCL",&MSCL);
                Target_tree->SetBranchAddress("SlantDepth",&SlantDepth);
                Target_tree->SetBranchAddress("theta2",&theta2);
                for (int entry=0;entry<Target_tree->GetEntries();entry++) {
                        Target_tree->GetEntry(entry);
                        int e = Hist_Target_SR_ErecS.at(0).FindBin(ErecS*1000.)-1;
                        if (e>=N_energy_bins) continue;
                        if (e<0) continue;
                        if (!FoV()) continue;
                        if (SignalSelectionMSCW()) Hist_Target_SR_MSCW.at(e).Fill(MSCW);
                        if (ControlSelectionMSCW()) Hist_Target_CR_MSCW.at(e).Fill(MSCW);
                        if (AuxSignalSelectionMSCW()) Hist_Target_ASR_MSCW.at(e).Fill(MSCW);
                        if (AuxControlSelectionMSCW()) Hist_Target_ACR_MSCW.at(e).Fill(MSCW);
                        if (SignalSelectionMSCL()) Hist_Target_SR_MSCL.at(e).Fill(MSCL);
                        if (ControlSelectionMSCL()) Hist_Target_CR_MSCL.at(e).Fill(MSCL);
                        if (AuxSignalSelectionMSCL()) Hist_Target_ASR_MSCL.at(e).Fill(MSCL);
                        if (AuxControlSelectionMSCL()) Hist_Target_ACR_MSCL.at(e).Fill(MSCL);
                        if (SignalSelectionDepth()) Hist_Target_SR_Depth.at(e).Fill(SlantDepth*100./37.);
                        if (ControlSelectionDepth()) Hist_Target_CR_Depth.at(e).Fill(SlantDepth*100./37.);
                        if (AuxSignalSelectionDepth()) Hist_Target_ASR_Depth.at(e).Fill(SlantDepth*100./37.);
                        if (AuxControlSelectionDepth()) Hist_Target_ACR_Depth.at(e).Fill(SlantDepth*100./37.);
                }
                input_file->Close();

        }
        
        // Target deconvolution method
        for (int e=0;e<N_energy_bins;e++) {
                Deconvolution(&Hist_Target_ACR_MSCW.at(e), &Hist_Target_ASR_MSCW.at(e), &Hist_Target_Deconv_MSCW.at(e));
                Deconvolution(&Hist_Target_ACR_MSCW.at(e), &Hist_Target_Deconv_MSCW.at(e), &Hist_Target_ABkgTemp_MSCW.at(e));
                Deconvolution(&Hist_Target_CR_MSCW.at(e), &Hist_Target_Deconv_MSCW.at(e), &Hist_Target_BkgTemp_MSCW.at(e));
                Deconvolution(&Hist_Target_ACR_MSCL.at(e), &Hist_Target_ASR_MSCL.at(e), &Hist_Target_Deconv_MSCL.at(e));
                Deconvolution(&Hist_Target_ACR_MSCL.at(e), &Hist_Target_Deconv_MSCL.at(e), &Hist_Target_ABkgTemp_MSCL.at(e));
                Deconvolution(&Hist_Target_CR_MSCL.at(e), &Hist_Target_Deconv_MSCL.at(e), &Hist_Target_BkgTemp_MSCL.at(e));
                Deconvolution(&Hist_Target_ACR_Depth.at(e), &Hist_Target_CR_Depth.at(e), &Hist_Target_Deconv_Depth.at(e));
                Deconvolution(&Hist_Target_ACR_Depth.at(e), &Hist_Target_Deconv_Depth.at(e), &Hist_Target_ABkgTemp_Depth.at(e));
                Deconvolution(&Hist_Target_ASR_Depth.at(e), &Hist_Target_Deconv_Depth.at(e), &Hist_Target_BkgTemp_Depth.at(e));
        }

        for (int e=0;e<N_energy_bins;e++) {
                double offset = Hist_Target_ASR_MSCW.at(e).GetMean()-Hist_Target_ABkgTemp_MSCW.at(e).GetMean();
                for (int i=0;i<N_bins_for_deconv;i++) {
                        int b = Hist_Target_Bkg_MSCW.at(e).FindBin(Hist_Target_Bkg_MSCW.at(e).GetBinCenter(i+1)-offset);
                        Hist_Target_ABkg_MSCW.at(e).SetBinContent(i+1,Hist_Target_ABkgTemp_MSCW.at(e).GetBinContent(b));
                        Hist_Target_ABkg_MSCW.at(e).SetBinError(i+1,Hist_Target_ABkgTemp_MSCW.at(e).GetBinError(b));
                        Hist_Target_Bkg_MSCW.at(e).SetBinContent(i+1,Hist_Target_BkgTemp_MSCW.at(e).GetBinContent(b));
                        Hist_Target_Bkg_MSCW.at(e).SetBinError(i+1,Hist_Target_BkgTemp_MSCW.at(e).GetBinError(b));
                }
                int norm_bin_low_target = Hist_Target_Bkg_MSCW.at(e).FindBin(2*MSCW_cut_upper);
                int norm_bin_up_target = Hist_Target_Bkg_MSCW.at(e).FindBin(6*MSCW_cut_upper);
                double scale_target = Hist_Target_SR_MSCW.at(e).Integral(norm_bin_low_target,norm_bin_up_target)/Hist_Target_Bkg_MSCW.at(e).Integral(norm_bin_low_target,norm_bin_up_target);
                Hist_Target_Bkg_MSCW.at(e).Scale(scale_target);
                scale_target = Hist_Target_ASR_MSCW.at(e).Integral(norm_bin_low_target,norm_bin_up_target)/Hist_Target_ABkg_MSCW.at(e).Integral(norm_bin_low_target,norm_bin_up_target);
                Hist_Target_ABkg_MSCW.at(e).Scale(scale_target);
                scale_target = Hist_Target_SR_MSCW.at(e).Integral(norm_bin_low_target,norm_bin_up_target)/Hist_Target_CR_MSCW.at(e).Integral(norm_bin_low_target,norm_bin_up_target);
                if (!(scale_target>0)) scale_target = 0;
                if (!(Hist_Target_CR_MSCW.at(e).Integral(norm_bin_low_target,norm_bin_up_target)>0)) scale_target = 0;
                Hist_Target_CR_MSCW.at(e).Scale(scale_target);
        }

        for (int e=0;e<N_energy_bins;e++) {
                double offset = Hist_Target_ASR_MSCL.at(e).GetMean()-Hist_Target_ABkgTemp_MSCL.at(e).GetMean();
                for (int i=0;i<N_bins_for_deconv;i++) {
                        int b = Hist_Target_Bkg_MSCL.at(e).FindBin(Hist_Target_Bkg_MSCL.at(e).GetBinCenter(i+1)-offset);
                        Hist_Target_ABkg_MSCL.at(e).SetBinContent(i+1,Hist_Target_ABkgTemp_MSCL.at(e).GetBinContent(b));
                        Hist_Target_ABkg_MSCL.at(e).SetBinError(i+1,Hist_Target_ABkgTemp_MSCL.at(e).GetBinError(b));
                        Hist_Target_Bkg_MSCL.at(e).SetBinContent(i+1,Hist_Target_BkgTemp_MSCL.at(e).GetBinContent(b));
                        Hist_Target_Bkg_MSCL.at(e).SetBinError(i+1,Hist_Target_BkgTemp_MSCL.at(e).GetBinError(b));
                }
                std::cout << e << ", Deconv " << Hist_Target_Deconv_MSCL.at(e).Integral() << std::endl;
                std::cout << e << ", SR " << Hist_Target_SR_MSCL.at(e).Integral() << std::endl;
                std::cout << e << ", Bkg " << Hist_Target_Bkg_MSCL.at(e).Integral() << std::endl;
                int norm_bin_low_target = Hist_Target_Bkg_MSCL.at(e).FindBin(2*MSCL_cut_upper);
                int norm_bin_up_target = Hist_Target_Bkg_MSCL.at(e).FindBin(6*MSCL_cut_upper);
                double scale_target = 0;
                scale_target = Hist_Target_SR_MSCL.at(e).Integral(norm_bin_low_target,norm_bin_up_target)/Hist_Target_Bkg_MSCL.at(e).Integral(norm_bin_low_target,norm_bin_up_target);
                if (!(Hist_Target_SR_MSCL.at(e).Integral()>0)) scale_target = 0;
                if (!(Hist_Target_Bkg_MSCL.at(e).Integral()>0)) scale_target = 0;
                Hist_Target_Bkg_MSCL.at(e).Scale(scale_target);
                scale_target = Hist_Target_ASR_MSCL.at(e).Integral(norm_bin_low_target,norm_bin_up_target)/Hist_Target_ABkg_MSCL.at(e).Integral(norm_bin_low_target,norm_bin_up_target);
                if (!(Hist_Target_ASR_MSCL.at(e).Integral()>0)) scale_target = 0;
                if (!(Hist_Target_ABkg_MSCL.at(e).Integral()>0)) scale_target = 0;
                Hist_Target_ABkg_MSCL.at(e).Scale(scale_target);
                std::cout << e << ", Bkg (scaled) " << Hist_Target_Bkg_MSCL.at(e).Integral() << std::endl;
                scale_target = Hist_Target_SR_MSCL.at(e).Integral(norm_bin_low_target,norm_bin_up_target)/Hist_Target_CR_MSCL.at(e).Integral(norm_bin_low_target,norm_bin_up_target);
                if (!(scale_target>0)) scale_target = 0;
                if (!(Hist_Target_CR_MSCL.at(e).Integral(norm_bin_low_target,norm_bin_up_target)>0)) scale_target = 0;
                Hist_Target_CR_MSCL.at(e).Scale(scale_target);
        }

        for (int e=0;e<N_energy_bins;e++) {
                double offset = Hist_Target_CR_Depth.at(e).GetMean()-Hist_Target_ABkgTemp_Depth.at(e).GetMean();
                double scale_target = 0;
                for (int i=0;i<N_bins_for_deconv;i++) {
                        int b = Hist_Target_Bkg_Depth.at(e).FindBin(Hist_Target_Bkg_Depth.at(e).GetBinCenter(i+1)-offset);
                        Hist_Target_ABkg_Depth.at(e).SetBinContent(i+1,Hist_Target_ABkgTemp_Depth.at(e).GetBinContent(b));
                        Hist_Target_ABkg_Depth.at(e).SetBinError(i+1,Hist_Target_ABkgTemp_Depth.at(e).GetBinError(b));
                        Hist_Target_Bkg_Depth.at(e).SetBinContent(i+1,Hist_Target_BkgTemp_Depth.at(e).GetBinContent(b));
                        Hist_Target_Bkg_Depth.at(e).SetBinError(i+1,Hist_Target_BkgTemp_Depth.at(e).GetBinError(b));
                }
                std::cout << e << ", Deconv " << Hist_Target_Deconv_Depth.at(e).Integral() << std::endl;
                std::cout << e << ", SR " << Hist_Target_SR_Depth.at(e).Integral() << std::endl;
                std::cout << e << ", Bkg " << Hist_Target_Bkg_Depth.at(e).Integral() << std::endl;
                double cut_mean = 5.+log2(pow(ErecS*1000./0.08,0.4));
                Depth_cut_lower = cut_mean-Depth_cut_width;
                Depth_cut_upper = cut_mean+Depth_cut_width;
                int norm_bin_low_target = Hist_Target_Bkg_Depth.at(e).FindBin(Depth_cut_lower);
                int norm_bin_up_target = Hist_Target_Bkg_Depth.at(e).FindBin(Depth_cut_upper);
                double numerator = Hist_Target_SR_Depth.at(e).Integral(1,norm_bin_low_target)+Hist_Target_SR_Depth.at(e).Integral(norm_bin_up_target,Hist_Dark_Bkg_Depth.at(e).GetNbinsX());
                double denominator = Hist_Target_Bkg_Depth.at(e).Integral(1,norm_bin_low_target)+Hist_Target_Bkg_Depth.at(e).Integral(norm_bin_up_target,Hist_Dark_Bkg_Depth.at(e).GetNbinsX());
                scale_target = numerator/denominator;
                if (!(scale_target>0)) scale_target = 0;
                if (!(denominator>0)) scale_target = 0;
                Hist_Target_Bkg_Depth.at(e).Scale(scale_target);
                numerator = Hist_Target_ASR_Depth.at(e).Integral(1,norm_bin_low_target)+Hist_Target_ASR_Depth.at(e).Integral(norm_bin_up_target,Hist_Dark_Bkg_Depth.at(e).GetNbinsX());
                denominator = Hist_Target_ABkg_Depth.at(e).Integral(1,norm_bin_low_target)+Hist_Target_ABkg_Depth.at(e).Integral(norm_bin_up_target,Hist_Dark_Bkg_Depth.at(e).GetNbinsX());
                scale_target = numerator/denominator;
                if (!(scale_target>0)) scale_target = 0;
                if (!(denominator>0)) scale_target = 0;
                Hist_Target_ABkg_Depth.at(e).Scale(scale_target);
                std::cout << e << ", Bkg (scaled) " << Hist_Target_Bkg_Depth.at(e).Integral() << std::endl;
                numerator = Hist_Target_SR_Depth.at(e).Integral(1,norm_bin_low_target)+Hist_Target_SR_Depth.at(e).Integral(norm_bin_up_target,Hist_Dark_Bkg_Depth.at(e).GetNbinsX());
                denominator = Hist_Target_CR_Depth.at(e).Integral(1,norm_bin_low_target)+Hist_Target_CR_Depth.at(e).Integral(norm_bin_up_target,Hist_Dark_Bkg_Depth.at(e).GetNbinsX());
                scale_target = numerator/denominator;
                if (!(scale_target>0)) scale_target = 0;
                if (!(denominator>0)) scale_target = 0;
                Hist_Target_CR_Depth.at(e).Scale(scale_target);
        }

        // Dark run method
/*
        for (int e=0;e<N_energy_bins;e++) {
                Deconvolution(&Hist_Dark_ASR_MSCW.at(e), &Hist_Target_ASR_MSCW.at(e), &Hist_Dark_Deconv_MSCW.at(e));
                Deconvolution(&Hist_Dark_ASR_MSCW.at(e), &Hist_Dark_Deconv_MSCW.at(e), &Hist_Dark_ABkgTemp_MSCW.at(e));
                Deconvolution(&Hist_Dark_SR_MSCW.at(e), &Hist_Dark_Deconv_MSCW.at(e), &Hist_Dark_BkgTemp_MSCW.at(e));
                Deconvolution(&Hist_Dark_ASR_MSCL.at(e), &Hist_Target_ASR_MSCL.at(e), &Hist_Dark_Deconv_MSCL.at(e));
                Deconvolution(&Hist_Dark_ASR_MSCL.at(e), &Hist_Dark_Deconv_MSCL.at(e), &Hist_Dark_ABkgTemp_MSCL.at(e));
                Deconvolution(&Hist_Dark_SR_MSCL.at(e), &Hist_Dark_Deconv_MSCL.at(e), &Hist_Dark_BkgTemp_MSCL.at(e));
                Deconvolution(&Hist_Dark_ASR_Depth.at(e), &Hist_Target_ASR_Depth.at(e), &Hist_Dark_Deconv_Depth.at(e));
                Deconvolution(&Hist_Dark_ASR_Depth.at(e), &Hist_Dark_Deconv_Depth.at(e), &Hist_Dark_ABkgTemp_Depth.at(e));
                Deconvolution(&Hist_Dark_SR_Depth.at(e), &Hist_Dark_Deconv_Depth.at(e), &Hist_Dark_BkgTemp_Depth.at(e));
        }

        for (int e=0;e<N_energy_bins;e++) {
                double offset = Hist_Target_ASR_MSCW.at(e).GetMean()-Hist_Dark_ABkgTemp_MSCW.at(e).GetMean();
                for (int i=0;i<N_bins_for_deconv;i++) {
                        int b = Hist_Dark_Bkg_MSCW.at(e).FindBin(Hist_Dark_Bkg_MSCW.at(e).GetBinCenter(i+1)-offset);
                        Hist_Dark_ABkg_MSCW.at(e).SetBinContent(i+1,Hist_Dark_ABkgTemp_MSCW.at(e).GetBinContent(b));
                        Hist_Dark_ABkg_MSCW.at(e).SetBinError(i+1,Hist_Dark_ABkgTemp_MSCW.at(e).GetBinError(b));
                        Hist_Dark_Bkg_MSCW.at(e).SetBinContent(i+1,Hist_Dark_BkgTemp_MSCW.at(e).GetBinContent(b));
                        Hist_Dark_Bkg_MSCW.at(e).SetBinError(i+1,Hist_Dark_BkgTemp_MSCW.at(e).GetBinError(b));
                }
                int norm_bin_low_target = Hist_Dark_Bkg_MSCW.at(e).FindBin(2*MSCW_cut_upper);
                int norm_bin_up_target = Hist_Dark_Bkg_MSCW.at(e).FindBin(3*MSCW_cut_upper);
                double scale_target = Hist_Target_SR_MSCW.at(e).Integral(norm_bin_low_target,norm_bin_up_target)/Hist_Dark_Bkg_MSCW.at(e).Integral(norm_bin_low_target,norm_bin_up_target);
                if (!(scale_target>0)) scale_target = 0;
                if (!(Hist_Dark_Bkg_MSCW.at(e).Integral(norm_bin_low_target,norm_bin_up_target)>0)) scale_target = 0;
                Hist_Dark_Bkg_MSCW.at(e).Scale(scale_target);
                scale_target = Hist_Target_ASR_MSCW.at(e).Integral(norm_bin_low_target,norm_bin_up_target)/Hist_Dark_ABkg_MSCW.at(e).Integral(norm_bin_low_target,norm_bin_up_target);
                if (!(scale_target>0)) scale_target = 0;
                if (!(Hist_Dark_ABkg_MSCW.at(e).Integral(norm_bin_low_target,norm_bin_up_target)>0)) scale_target = 0;
                Hist_Dark_ABkg_MSCW.at(e).Scale(scale_target);
                scale_target = Hist_Target_SR_MSCW.at(e).Integral(norm_bin_low_target,norm_bin_up_target)/Hist_Dark_SR_MSCW.at(e).Integral(norm_bin_low_target,norm_bin_up_target);
                if (!(scale_target>0)) scale_target = 0;
                if (!(Hist_Dark_SR_MSCW.at(e).Integral(norm_bin_low_target,norm_bin_up_target)>0)) scale_target = 0;
                Hist_Dark_SR_MSCW.at(e).Scale(scale_target);
        }

        for (int e=0;e<N_energy_bins;e++) {
                double offset = Hist_Target_ASR_MSCL.at(e).GetMean()-Hist_Dark_ABkgTemp_MSCL.at(e).GetMean();
                for (int i=0;i<N_bins_for_deconv;i++) {
                        int b = Hist_Dark_Bkg_MSCL.at(e).FindBin(Hist_Dark_Bkg_MSCL.at(e).GetBinCenter(i+1)-offset);
                        Hist_Dark_ABkg_MSCL.at(e).SetBinContent(i+1,Hist_Dark_ABkgTemp_MSCL.at(e).GetBinContent(b));
                        Hist_Dark_ABkg_MSCL.at(e).SetBinError(i+1,Hist_Dark_ABkgTemp_MSCL.at(e).GetBinError(b));
                        Hist_Dark_Bkg_MSCL.at(e).SetBinContent(i+1,Hist_Dark_BkgTemp_MSCL.at(e).GetBinContent(b));
                        Hist_Dark_Bkg_MSCL.at(e).SetBinError(i+1,Hist_Dark_BkgTemp_MSCL.at(e).GetBinError(b));
                }
                int norm_bin_low_target = Hist_Dark_Bkg_MSCL.at(e).FindBin(2*MSCL_cut_upper);
                int norm_bin_up_target = Hist_Dark_Bkg_MSCL.at(e).FindBin(3*MSCL_cut_upper);
                double scale_target = Hist_Target_SR_MSCL.at(e).Integral(norm_bin_low_target,norm_bin_up_target)/Hist_Dark_Bkg_MSCL.at(e).Integral(norm_bin_low_target,norm_bin_up_target);
                if (!(scale_target>0)) scale_target = 0;
                if (!(Hist_Dark_Bkg_MSCL.at(e).Integral(norm_bin_low_target,norm_bin_up_target)>0)) scale_target = 0;
                Hist_Dark_Bkg_MSCL.at(e).Scale(scale_target);
                scale_target = Hist_Target_ASR_MSCL.at(e).Integral(norm_bin_low_target,norm_bin_up_target)/Hist_Dark_ABkg_MSCL.at(e).Integral(norm_bin_low_target,norm_bin_up_target);
                if (!(scale_target>0)) scale_target = 0;
                if (!(Hist_Dark_ABkg_MSCL.at(e).Integral(norm_bin_low_target,norm_bin_up_target)>0)) scale_target = 0;
                Hist_Dark_ABkg_MSCL.at(e).Scale(scale_target);
                scale_target = Hist_Target_SR_MSCL.at(e).Integral(norm_bin_low_target,norm_bin_up_target)/Hist_Dark_SR_MSCL.at(e).Integral(norm_bin_low_target,norm_bin_up_target);
                if (!(scale_target>0)) scale_target = 0;
                if (!(Hist_Dark_SR_MSCL.at(e).Integral(norm_bin_low_target,norm_bin_up_target)>0)) scale_target = 0;
                Hist_Dark_SR_MSCL.at(e).Scale(scale_target);
        }

        for (int e=0;e<N_energy_bins;e++) {
                double offset = Hist_Target_ASR_Depth.at(e).GetMean()-Hist_Dark_ABkgTemp_Depth.at(e).GetMean();
                for (int i=0;i<N_bins_for_deconv;i++) {
                        int b = Hist_Dark_Bkg_Depth.at(e).FindBin(Hist_Dark_Bkg_Depth.at(e).GetBinCenter(i+1)-offset);
                        Hist_Dark_ABkg_Depth.at(e).SetBinContent(i+1,Hist_Dark_ABkgTemp_Depth.at(e).GetBinContent(b));
                        Hist_Dark_ABkg_Depth.at(e).SetBinError(i+1,Hist_Dark_ABkgTemp_Depth.at(e).GetBinError(b));
                        Hist_Dark_Bkg_Depth.at(e).SetBinContent(i+1,Hist_Dark_BkgTemp_Depth.at(e).GetBinContent(b));
                        Hist_Dark_Bkg_Depth.at(e).SetBinError(i+1,Hist_Dark_BkgTemp_Depth.at(e).GetBinError(b));
                }
                int norm_bin_low_target = Hist_Dark_Bkg_Depth.at(e).FindBin(14);
                int norm_bin_up_target = Hist_Dark_Bkg_Depth.at(e).FindBin(100);
                double scale_target = Hist_Target_SR_Depth.at(e).Integral(norm_bin_low_target,norm_bin_up_target)/Hist_Dark_Bkg_Depth.at(e).Integral(norm_bin_low_target,norm_bin_up_target);
                if (!(scale_target>0)) scale_target = 0;
                if (!(Hist_Dark_Bkg_Depth.at(e).Integral(norm_bin_low_target,norm_bin_up_target)>0)) scale_target = 0;
                Hist_Dark_Bkg_Depth.at(e).Scale(scale_target);
                scale_target = Hist_Target_ASR_Depth.at(e).Integral(norm_bin_low_target,norm_bin_up_target)/Hist_Dark_ABkg_Depth.at(e).Integral(norm_bin_low_target,norm_bin_up_target);
                if (!(scale_target>0)) scale_target = 0;
                if (!(Hist_Dark_ABkg_Depth.at(e).Integral(norm_bin_low_target,norm_bin_up_target)>0)) scale_target = 0;
                Hist_Dark_ABkg_Depth.at(e).Scale(scale_target);
                scale_target = Hist_Target_SR_Depth.at(e).Integral(norm_bin_low_target,norm_bin_up_target)/Hist_Dark_SR_Depth.at(e).Integral(norm_bin_low_target,norm_bin_up_target);
                if (!(scale_target>0)) scale_target = 0;
                if (!(Hist_Dark_SR_Depth.at(e).Integral(norm_bin_low_target,norm_bin_up_target)>0)) scale_target = 0;
                Hist_Dark_SR_Depth.at(e).Scale(scale_target);
        }
*/

        for (int e=0;e<N_energy_bins;e++) {
                Hist_Dark_Deconv_MSCW.at(e).Reset();
                Hist_Dark_Deconv_MSCW.at(e).Add(&Hist_Target_ASR_MSCW.at(e));
                Hist_Dark_Deconv_MSCW.at(e).Divide(&Hist_Dark_ASR_MSCW.at(e));
                Hist_Dark_Bkg_MSCW.at(e).Reset();
                Hist_Dark_Bkg_MSCW.at(e).Add(&Hist_Dark_SR_MSCW.at(e));
                Hist_Dark_Bkg_MSCW.at(e).Multiply(&Hist_Dark_Deconv_MSCW.at(e));
                Hist_Dark_Deconv_MSCL.at(e).Reset();
                Hist_Dark_Deconv_MSCL.at(e).Add(&Hist_Target_ASR_MSCL.at(e));
                Hist_Dark_Deconv_MSCL.at(e).Divide(&Hist_Dark_ASR_MSCL.at(e));
                Hist_Dark_Bkg_MSCL.at(e).Reset();
                Hist_Dark_Bkg_MSCL.at(e).Add(&Hist_Dark_SR_MSCL.at(e));
                Hist_Dark_Bkg_MSCL.at(e).Multiply(&Hist_Dark_Deconv_MSCL.at(e));
                Hist_Dark_Deconv_Depth.at(e).Reset();
                Hist_Dark_Deconv_Depth.at(e).Add(&Hist_Target_ACR_Depth.at(e));
                Hist_Dark_Deconv_Depth.at(e).Divide(&Hist_Dark_ACR_Depth.at(e));
                Hist_Dark_Bkg_Depth.at(e).Reset();
                Hist_Dark_Bkg_Depth.at(e).Add(&Hist_Dark_SR_Depth.at(e));
                Hist_Dark_Bkg_Depth.at(e).Multiply(&Hist_Dark_Deconv_Depth.at(e));
        }
        for (int e=0;e<N_energy_bins;e++) {
                int norm_bin_low_target = Hist_Dark_Bkg_MSCW.at(e).FindBin(2*MSCW_cut_upper);
                int norm_bin_up_target = Hist_Dark_Bkg_MSCW.at(e).FindBin(3*MSCW_cut_upper);
                double scale_target = Hist_Target_SR_MSCW.at(e).Integral(norm_bin_low_target,norm_bin_up_target)/Hist_Dark_Bkg_MSCW.at(e).Integral(norm_bin_low_target,norm_bin_up_target);
                if (!(scale_target>0)) scale_target = 0;
                if (!(Hist_Dark_Bkg_MSCW.at(e).Integral(norm_bin_low_target,norm_bin_up_target)>0)) scale_target = 0;
                Hist_Dark_Bkg_MSCW.at(e).Scale(scale_target);
                scale_target = Hist_Target_SR_MSCW.at(e).Integral(norm_bin_low_target,norm_bin_up_target)/Hist_Dark_SR_MSCW.at(e).Integral(norm_bin_low_target,norm_bin_up_target);
                if (!(scale_target>0)) scale_target = 0;
                if (!(Hist_Dark_SR_MSCW.at(e).Integral(norm_bin_low_target,norm_bin_up_target)>0)) scale_target = 0;
                Hist_Dark_SR_MSCW.at(e).Scale(scale_target);

                norm_bin_low_target = Hist_Dark_Bkg_MSCL.at(e).FindBin(2*MSCL_cut_upper);
                norm_bin_up_target = Hist_Dark_Bkg_MSCL.at(e).FindBin(3*MSCL_cut_upper);
                scale_target = Hist_Target_SR_MSCL.at(e).Integral(norm_bin_low_target,norm_bin_up_target)/Hist_Dark_Bkg_MSCL.at(e).Integral(norm_bin_low_target,norm_bin_up_target);
                if (!(scale_target>0)) scale_target = 0;
                if (!(Hist_Dark_Bkg_MSCL.at(e).Integral(norm_bin_low_target,norm_bin_up_target)>0)) scale_target = 0;
                Hist_Dark_Bkg_MSCL.at(e).Scale(scale_target);
                scale_target = Hist_Target_SR_MSCL.at(e).Integral(norm_bin_low_target,norm_bin_up_target)/Hist_Dark_SR_MSCL.at(e).Integral(norm_bin_low_target,norm_bin_up_target);
                if (!(scale_target>0)) scale_target = 0;
                if (!(Hist_Dark_SR_MSCL.at(e).Integral(norm_bin_low_target,norm_bin_up_target)>0)) scale_target = 0;
                Hist_Dark_SR_MSCL.at(e).Scale(scale_target);

                double cut_mean = 5.+log2(pow(ErecS*1000./0.08,0.4));
                Depth_cut_lower = cut_mean-Depth_cut_width;
                Depth_cut_upper = cut_mean+Depth_cut_width;
                norm_bin_low_target = Hist_Dark_Bkg_Depth.at(e).FindBin(Depth_cut_lower);
                norm_bin_up_target = Hist_Dark_Bkg_Depth.at(e).FindBin(Depth_cut_upper);
                double numerator = Hist_Target_SR_Depth.at(e).Integral(1,norm_bin_low_target)+Hist_Target_SR_Depth.at(e).Integral(norm_bin_up_target,Hist_Dark_Bkg_Depth.at(e).GetNbinsX());
                double denominator = Hist_Dark_Bkg_Depth.at(e).Integral(1,norm_bin_low_target)+Hist_Dark_Bkg_Depth.at(e).Integral(norm_bin_up_target,Hist_Dark_Bkg_Depth.at(e).GetNbinsX());
                scale_target = numerator/denominator;
                if (!(scale_target>0)) scale_target = 0;
                if (!(denominator>0)) scale_target = 0;
                Hist_Dark_Bkg_Depth.at(e).Scale(scale_target);
                numerator = Hist_Target_SR_Depth.at(e).Integral(1,norm_bin_low_target)+Hist_Target_SR_Depth.at(e).Integral(norm_bin_up_target,Hist_Dark_Bkg_Depth.at(e).GetNbinsX());
                denominator = Hist_Dark_SR_Depth.at(e).Integral(1,norm_bin_low_target)+Hist_Dark_SR_Depth.at(e).Integral(norm_bin_up_target,Hist_Dark_Bkg_Depth.at(e).GetNbinsX());
                scale_target = numerator/denominator;
                if (!(scale_target>0)) scale_target = 0;
                if (!(denominator>0)) scale_target = 0;
                Hist_Dark_SR_Depth.at(e).Scale(scale_target);
        }


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
        Hist_Dark_TelElevAzim.Write();
        for (int e=0;e<N_energy_bins;e++) {
                Hist_Dark_ASR_MSCW.at(e).Write();
                Hist_Dark_ACR_MSCW.at(e).Write();
                Hist_Dark_ABkg_MSCW.at(e).Write();
                Hist_Dark_SR_MSCW.at(e).Write();
                Hist_Dark_CR_MSCW.at(e).Write();
                Hist_Dark_Bkg_MSCW.at(e).Write();
                Hist_Dark_Deconv_MSCW.at(e).Write();
                Hist_Dark_ASR_MSCL.at(e).Write();
                Hist_Dark_ACR_MSCL.at(e).Write();
                Hist_Dark_ABkg_MSCL.at(e).Write();
                Hist_Dark_SR_MSCL.at(e).Write();
                Hist_Dark_CR_MSCL.at(e).Write();
                Hist_Dark_Bkg_MSCL.at(e).Write();
                Hist_Dark_Deconv_MSCL.at(e).Write();
                Hist_Dark_ASR_Depth.at(e).Write();
                Hist_Dark_ACR_Depth.at(e).Write();
                Hist_Dark_ABkg_Depth.at(e).Write();
                Hist_Dark_SR_Depth.at(e).Write();
                Hist_Dark_CR_Depth.at(e).Write();
                Hist_Dark_Bkg_Depth.at(e).Write();
                Hist_Dark_Deconv_Depth.at(e).Write();
                Hist_Target_ASR_MSCW.at(e).Write();
                Hist_Target_ACR_MSCW.at(e).Write();
                Hist_Target_ABkg_MSCW.at(e).Write();
                Hist_Target_SR_MSCW.at(e).Write();
                Hist_Target_CR_MSCW.at(e).Write();
                Hist_Target_Bkg_MSCW.at(e).Write();
                Hist_Target_Deconv_MSCW.at(e).Write();
                Hist_Target_ASR_MSCL.at(e).Write();
                Hist_Target_ACR_MSCL.at(e).Write();
                Hist_Target_ABkg_MSCL.at(e).Write();
                Hist_Target_SR_MSCL.at(e).Write();
                Hist_Target_CR_MSCL.at(e).Write();
                Hist_Target_Bkg_MSCL.at(e).Write();
                Hist_Target_Deconv_MSCL.at(e).Write();
                Hist_Target_ASR_Depth.at(e).Write();
                Hist_Target_ACR_Depth.at(e).Write();
                Hist_Target_ABkg_Depth.at(e).Write();
                Hist_Target_SR_Depth.at(e).Write();
                Hist_Target_CR_Depth.at(e).Write();
                Hist_Target_Bkg_Depth.at(e).Write();
                Hist_Target_Deconv_Depth.at(e).Write();
        }
        OutputFile.Close();

}
