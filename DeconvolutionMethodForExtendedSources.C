
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

double MSCW_cut_lower = -0.5;
double MSCW_cut_upper = 0.5;
double MSCL_cut_lower = -0.5;
double MSCL_cut_upper = 0.5;

double Depth_cut_lower = 6;
double Depth_cut_upper = 14;
//double Depth_cut_width = 3;
double Depth_cut_width = 4;
//double Depth_cut_width = 0;

double Control_MSCL_cut = 0.;

//double Norm_Lower = 0.0;
double Norm_Lower = 1.0;
//double Norm_Upper = 8.5;
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
const int N_energy_bins = 7;
double energy_bins[N_energy_bins+1] = {300,400,600,1000,1500,2000,4000,10000};
//const int N_energy_bins = 4;
//double energy_bins[N_energy_bins+1] = {1000,1500,2000,4000,10000};
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
    //if (theta2<4.0) return true;
    //if (theta2>2.0) return true;
    //if (theta2<0.5) return true;
    if (theta2<1.0) return true;
    //if (theta2>0.2 && theta2<4.0) return true;
    return false;
}
bool RingFoV() {
    //if (theta2<0.5) return false;
    if (theta2<1.0) return false;
    //if (theta2>4.0) return false;
    return true;
}
bool QualitySelection() {
    if (EmissionHeightChi2/EmissionHeight>0.2) return false;
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
    //if (SlantDepth*100./37.<Depth_cut_lower) return false;
    //if (SlantDepth*100./37.>Depth_cut_upper) return false;
    return true;
}
bool Aux1SignalSelectionMSCW() {
    if (MSCL<MSCL_cut_lower) return false;
    if (MSCL>MSCL_cut_upper) return false;
    double cut_mean = 5.+log2(pow(ErecS*1000./0.08,0.4));
    Depth_cut_lower = cut_mean-Depth_cut_width;
    Depth_cut_upper = cut_mean+Depth_cut_width;
    if (SlantDepth*100./37.>Depth_cut_lower) return false;
    //if (SlantDepth*100./37.<Depth_cut_lower-4.) return false;
    return true;
}
bool Aux2SignalSelectionMSCW() {
    if (MSCL<MSCL_cut_lower) return false;
    if (MSCL>MSCL_cut_upper) return false;
    double cut_mean = 5.+log2(pow(ErecS*1000./0.08,0.4));
    Depth_cut_lower = cut_mean-Depth_cut_width;
    Depth_cut_upper = cut_mean+Depth_cut_width;
    if (SlantDepth*100./37.<Depth_cut_upper) return false;
    //if (SlantDepth*100./37.>Depth_cut_upper+4.) return false;
    return true;
}
bool ControlSelectionMSCW() {
    if (MSCL<1.0) return false;
    if (MSCL>1.2) return false;
    double cut_mean = 5.+log2(pow(ErecS*1000./0.08,0.4));
    Depth_cut_lower = cut_mean-Depth_cut_width;
    Depth_cut_upper = cut_mean+Depth_cut_width;
    //if (SlantDepth*100./37.<Depth_cut_lower) return false;
    //if (SlantDepth*100./37.>Depth_cut_upper) return false;
    return true;
}
bool Aux1ControlSelectionMSCW() {
    if (MSCL<1.0) return false;
    if (MSCL>1.2) return false;
    double cut_mean = 5.+log2(pow(ErecS*1000./0.08,0.4));
    Depth_cut_lower = cut_mean-Depth_cut_width;
    Depth_cut_upper = cut_mean+Depth_cut_width;
    if (SlantDepth*100./37.>Depth_cut_lower) return false;
    //if (SlantDepth*100./37.<Depth_cut_lower-4.) return false;
    return true;
}
bool Aux2ControlSelectionMSCW() {
    if (MSCL<1.0) return false;
    if (MSCL>1.2) return false;
    double cut_mean = 5.+log2(pow(ErecS*1000./0.08,0.4));
    Depth_cut_lower = cut_mean-Depth_cut_width;
    Depth_cut_upper = cut_mean+Depth_cut_width;
    if (SlantDepth*100./37.<Depth_cut_upper) return false;
    //if (SlantDepth*100./37.>Depth_cut_upper+4.) return false;
    return true;
}
double background(Double_t *x, Double_t *par) {
    return (1./par[0])*(x[0]+par[1])*exp(-pow((x[0]+par[1])/par[0],2));
}
std::pair<double, double> RatioSRCR(TH1D* Hist_Bkg, double SR_cut_low, double SR_cut_up, double CR_cut) {
    int b_low = Hist_Bkg->FindBin(SR_cut_low);
    int b_up = Hist_Bkg->FindBin(SR_cut_up);
    int b_cr = Hist_Bkg->FindBin(CR_cut);
    double NSR_err;
    double NSR = Hist_Bkg->IntegralAndError(b_low,b_up-1,NSR_err);
    double NCR_err;
    double NCR = Hist_Bkg->IntegralAndError(b_cr,Hist_Bkg->GetNbinsX(),NCR_err);
    //std::cout << "NSR_err" << NSR_err << std::endl;
    if (NCR==0) return std::make_pair(0,0);
    return std::make_pair(NSR/NCR,NSR/NCR*pow(pow(NSR_err/NSR,2)+pow(NCR_err/NCR,2),0.5));
}
double GetChi2(TH1* Hist_SR, TH1* Hist_Bkg, double norm_low, double norm_up, bool includeSR) {
    double chi2_temp = 0;
    double mean = Hist_SR->GetMean();
    double rms = Hist_SR->GetRMS();
    for (int i=0;i<Hist_SR->GetNbinsX();i++) {
        if (Hist_Bkg->GetBinCenter(i+1)<-2.) continue;
        if (Hist_Bkg->GetBinCenter(i+1)>norm_up) continue;
        //if (Hist_Bkg->GetBinCenter(i+1)>mean+8.*rms) continue;
        if (Hist_Bkg->GetBinCenter(i+1)>mean+4.*rms) continue;
        double bkg = Hist_Bkg->GetBinContent(i+1);
        double data = Hist_SR->GetBinContent(i+1);
        double bkg_err = Hist_Bkg->GetBinError(i+1);
        double data_err = Hist_SR->GetBinError(i+1);
        if (data_err==0) data_err = 0.1;
        if (!includeSR && Hist_Bkg->GetBinCenter(i+1)>-0.7 && Hist_Bkg->GetBinCenter(i+1)<Norm_Lower) {
            continue;
        }
        //chi2_temp += pow(bkg-data,2)/Hist_SR->Integral();
        //chi2_temp += pow(pow(bkg-data,2),0.5)/Hist_SR->Integral();
        //chi2_temp += pow(bkg-data,2)/(bkg+data);
        //chi2_temp += pow(pow(bkg-data,2),0.5)/(data_err);
        chi2_temp += pow(bkg-data,2)/(data_err*data_err+bkg_err*bkg_err);
    }
    chi2_temp = 1./chi2_temp;
    return chi2_temp;
}
double ShiftAndNormalize(TH1* Hist_SR, TH1* Hist_BkgTemp, TH1* Hist_Bkg, double shift_begin, double norm_low, double norm_up, bool doShift) {
    double shift_fit = 0;
    double scale_fit = 0;
    double chi2 = 0.;
    if (doShift) {
        //for (int fit=0;fit<40;fit++) {
                //double shift = shift_begin-1.0+double(fit)*0.05;
        for (int fit=0;fit<80;fit++) {
                double shift = shift_begin-2.0+double(fit)*0.05;
                for (int i=0;i<Hist_SR->GetNbinsX();i++) {
                        int b = Hist_SR->FindBin(Hist_SR->GetBinCenter(i+1)-shift);
                        Hist_Bkg->SetBinContent(i+1,Hist_BkgTemp->GetBinContent(b));
                        Hist_Bkg->SetBinError(i+1,Hist_BkgTemp->GetBinError(b));
                }
                int norm_bin_low = Hist_SR->FindBin(norm_low);
                int norm_bin_up = Hist_SR->FindBin(norm_up);
                double scale = Hist_SR->Integral(norm_bin_low,norm_bin_up)/Hist_Bkg->Integral(norm_bin_low,norm_bin_up);
                Hist_Bkg->Scale(scale);
                double chi2_temp = GetChi2(Hist_SR, Hist_Bkg, norm_low, norm_up,false);
                if (chi2<chi2_temp) {
                    chi2 = chi2_temp;
                    shift_fit = shift;
                    scale_fit = scale;
                } 
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
    int norm_bin_low = Hist_SR->FindBin(norm_low);
    int norm_bin_up = Hist_SR->FindBin(norm_up);
    double scale = Hist_SR->Integral(norm_bin_low,norm_bin_up)/Hist_Bkg->Integral(norm_bin_low,norm_bin_up);
    Hist_Bkg->Scale(scale);
    return shift_fit;
}

void BuildSRB(TH1D* Hist_SR, TH1D* Hist_ASR1, TH1D* Hist_ASR2, TH1D* Hist_SRB) {
    Hist_SRB->Reset();
    double chi2_best = 0;
    double scale_best = 0;
    //for (int i=0;i<20;i++) {
    //    double scale = double(i)*0.1;
    //    Hist_SRB->Reset();
    //    Hist_SRB->Add(Hist_ASR1);
    //    Hist_SRB->Add(Hist_ASR2,scale);
    //    int norm_bin_low = Hist_SRB->FindBin(Norm_Lower);
    //    int norm_bin_up = Hist_SRB->FindBin(Norm_Upper);
    //    double norm = Hist_SR->Integral(norm_bin_low,norm_bin_up)/Hist_SRB->Integral(norm_bin_low,norm_bin_up);
    //    Hist_SRB->Scale(norm);
    //    double chi2 = GetChi2(Hist_SR, Hist_SRB, Norm_Lower, Norm_Upper, false);
    //    if (chi2_best<chi2) {
    //        chi2_best = chi2;
    //        scale_best = scale;
    //    }
    //}
    //Hist_SRB->Reset();
    //Hist_SRB->Add(Hist_ASR1);
    //Hist_SRB->Add(Hist_ASR2,scale_best);
    //int norm_bin_low = Hist_SRB->FindBin(Norm_Lower);
    //int norm_bin_up = Hist_SRB->FindBin(Norm_Upper);
    //double norm = Hist_SR->Integral(norm_bin_low,norm_bin_up)/Hist_SRB->Integral(norm_bin_low,norm_bin_up);
    //Hist_SRB->Scale(norm);
    for (int b=0;b<Hist_SRB->GetNbinsX();b++) {
        if (Hist_SRB->GetBinCenter(b+1)>Norm_Lower)
        {
          Hist_SRB->SetBinContent(b+1,Hist_SR->GetBinContent(b+1));
          Hist_SRB->SetBinError(b+1,Hist_SR->GetBinError(b+1));
        }
        else {
          Hist_SRB->SetBinContent(b+1,0);
          Hist_SRB->SetBinError(b+1,1.);
        }
    }
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
        //n_iteration = 10;
        //n_iteration = 20; // for 35-55 elevation
        //if (Hist_response->GetRMS()<=2.0) n_iteration = 100;
        //if (Hist_response->GetRMS()<=1.3) n_iteration = 1000;
        sp.Deconvolution(source,response,N_bins,n_iteration,1,100000); // new best option
        //sp.Deconvolution(source,response,N_bins,n_iteration,2,1.2);
        for (int i=0;i<N_bins;i++) {
            Hist_Deconv->SetBinContent(i+1,max(source[i],0.));
        }
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
                MSCL_cut_lower = 1.5;
                MSCL_cut_upper = 2.5;
        }

        TRandom rnd;
        TH2D Hist_Dark_TelElevAzim("Hist_Dark_TelElevAzim","",90,0,90,360,0,360);
        TH2D Hist_Target_TelElevAzim("Hist_Target_TelElevAzim","",90,0,90,360,0,360);
        TH2D Hist_Target_ON_ErecSDepth("Hist_Target_ON_ErecSDepth","",N_energy_bins_log,energy_bins_log,60,0,30);
        TH2D Hist_Target_OFF_ErecSDepth("Hist_Target_OFF_ErecSDepth","",N_energy_bins_log,energy_bins_log,60,0,30);
        TH2D Hist_Target_ON_ElevDepth("Hist_Target_ON_ElevDepth","",12,30,90,40,-10,10);
        TH2D Hist_Target_OFF_ElevDepth("Hist_Target_OFF_ElevDepth","",12,30,90,40,-10,10);
        TH1D Hist_Target_ON_MSCW_Alpha("Hist_Target_ON_MSCW_Alpha","",100,0,10);
        TH1D Hist_Target_OFF_MSCW_Alpha("Hist_Target_OFF_MSCW_Alpha","",100,0,10);
        TH1D Hist_Target_SR_Energy("Hist_Target_SR_Energy","",N_energy_bins_log,energy_bins_log);
        TH1D Hist_Target_CR_Energy("Hist_Target_CR_Energy","",N_energy_bins_log,energy_bins_log);
        TH1D Hist_Target_Bkg_Energy("Hist_Target_Bkg_Energy","",N_energy_bins_log,energy_bins_log);
        vector<TH1D> Hist_Dark_SR_ErecS;
        vector<TH1D> Hist_Dark_SR_MSCW;
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
            //N_bins_for_deconv = 7680;
            //N_bins_for_deconv = 3840;
            //N_bins_for_deconv = 1920;
            //N_bins_for_deconv = 960;
            N_bins_for_deconv = 480;
            Hist_Dark_SR_ErecS.push_back(TH1D("Hist_Dark_SR_ErecS_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_energy_bins,energy_bins));
            Hist_Dark_SR_MSCW.push_back(TH1D("Hist_Dark_SR_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-15,15));
            Hist_Dark_SRB_MSCW.push_back(TH1D("Hist_Dark_SRB_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-15,15));
            Hist_Dark_CR_MSCW.push_back(TH1D("Hist_Dark_CR_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-15,15));
            Hist_Dark_ASR1_MSCW.push_back(TH1D("Hist_Dark_ASR1_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-15,15));
            Hist_Dark_ACR1_MSCW.push_back(TH1D("Hist_Dark_ACR1_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-15,15));
            Hist_Dark_ASR2_MSCW.push_back(TH1D("Hist_Dark_ASR2_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-15,15));
            Hist_Dark_ACR2_MSCW.push_back(TH1D("Hist_Dark_ACR2_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-15,15));
            Hist_Dark_BkgTemp_MSCW.push_back(TH1D("Hist_Dark_BkgTemp_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-15,15));
            Hist_Dark_Bkg_MSCW.push_back(TH1D("Hist_Dark_Bkg_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-15,15));
            Hist_Dark_ABkg1Temp_MSCW.push_back(TH1D("Hist_Dark_ABkg1Temp_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-15,15));
            Hist_Dark_ABkg1_MSCW.push_back(TH1D("Hist_Dark_ABkg1_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-15,15));
            Hist_Dark_ABkg2Temp_MSCW.push_back(TH1D("Hist_Dark_ABkg2Temp_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-15,15));
            Hist_Dark_ABkg2_MSCW.push_back(TH1D("Hist_Dark_ABkg2_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-15,15));
            Hist_Dark_Deconv_MSCW.push_back(TH1D("Hist_Dark_Deconv_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-15,15));
            Hist_Dark_TrueDeconv_MSCW.push_back(TH1D("Hist_Dark_TrueDeconv_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-15,15));
            Hist_Dark_Elec_MSCW.push_back(TH1D("Hist_Dark_Elec_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-15,15));
            Hist_Target_SR_ErecS.push_back(TH1D("Hist_Target_SR_ErecS_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_energy_bins,energy_bins));
            Hist_Target_SR_MSCW.push_back(TH1D("Hist_Target_SR_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-15,15));
            Hist_Target_SRB_MSCW.push_back(TH1D("Hist_Target_SRB_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-15,15));
            Hist_Target_CR_MSCW.push_back(TH1D("Hist_Target_CR_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-15,15));
            Hist_Target_ASR1_MSCW.push_back(TH1D("Hist_Target_ASR1_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-15,15));
            Hist_Target_ACR1_MSCW.push_back(TH1D("Hist_Target_ACR1_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-15,15));
            Hist_Target_ASR2_MSCW.push_back(TH1D("Hist_Target_ASR2_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-15,15));
            Hist_Target_ACR2_MSCW.push_back(TH1D("Hist_Target_ACR2_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-15,15));
            Hist_Target_BkgTemp_MSCW.push_back(TH1D("Hist_Target_BkgTemp_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-15,15));
            Hist_Target_Bkg_MSCW.push_back(TH1D("Hist_Target_Bkg_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-15,15));
            Hist_Target_ABkg1Temp_MSCW.push_back(TH1D("Hist_Target_ABkg1Temp_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-15,15));
            Hist_Target_ABkg1_MSCW.push_back(TH1D("Hist_Target_ABkg1_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-15,15));
            Hist_Target_ABkg2Temp_MSCW.push_back(TH1D("Hist_Target_ABkg2Temp_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-15,15));
            Hist_Target_ABkg2_MSCW.push_back(TH1D("Hist_Target_ABkg2_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-15,15));
            Hist_Target_Deconv_MSCW.push_back(TH1D("Hist_Target_Deconv_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-15,15));
            Hist_Target_TrueDeconv_MSCW.push_back(TH1D("Hist_Target_TrueDeconv_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-15,15));
            Hist_Target_Elec_MSCW.push_back(TH1D("Hist_Target_Elec_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-15,15));
            Hist_Target_Ring_MSCW.push_back(TH1D("Hist_Target_Ring_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-15,15));
            Hist_Target_TotalBkg_MSCW.push_back(TH1D("Hist_Target_TotalBkg_MSCW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,-15,15));
            Hist_Target_SR_SkyMap.push_back(TH2D("Hist_Target_SR_SkyMap_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",100,-2,2,100,-2,2));
            Hist_Target_CR_SkyMap.push_back(TH2D("Hist_Target_CR_SkyMap_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",100,-2,2,100,-2,2));
            Hist_Target_Bkg_SkyMap.push_back(TH2D("Hist_Target_Bkg_SkyMap_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",100,-2,2,100,-2,2));
            Hist_Target_SR_Theta2.push_back(TH1D("Hist_Target_SR_Theta2_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",40,0,4));
            Hist_Target_CR_Theta2.push_back(TH1D("Hist_Target_CR_Theta2_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",40,0,4));
            Hist_Target_Bkg_Theta2.push_back(TH1D("Hist_Target_Bkg_Theta2_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",40,0,4));
        }

        std::cout << "Getting dark runs... " << std::endl;
        vector<int> Dark_runlist = GetRunList("Crab");
        //vector<int> Dark_runlist = GetRunList("PKS1424");
        for (int run=0;run<Dark_runlist.size();run++) {
                char run_number[50];
                sprintf(run_number, "%i", int(Dark_runlist[run]));
                filename = TString("$VERITAS_USER_DATA_DIR/Crab_V6_Moderate-TMVA-BDT.RB."+TString(run_number)+".root");
                //filename = TString("$VERITAS_USER_DATA_DIR/PKS1424_V6_Moderate-TMVA-BDT.RB."+TString(run_number)+".root");
                TFile*  input_file = TFile::Open(filename.c_str());
                //std::cout << "Open " << filename << std::endl;
                TTree* pointing_tree = (TTree*) input_file->Get("run_"+TString(run_number)+"/stereo/pointingDataReduced");
                pointing_tree->SetBranchAddress("TelElevation",&TelElevation);
                pointing_tree->SetBranchAddress("TelAzimuth",&TelAzimuth);
                pointing_tree->GetEntry(0);
                //std::cout << TelElevation << std::endl;
                //if (TelElevation<Elev_cut_lower) continue;
                //if (TelElevation>Elev_cut_upper) continue;
                if (TelElevation<55.) continue;
                if (TelElevation>85.) continue;
                for (int entry=0;entry<pointing_tree->GetEntries();entry++) {
                        pointing_tree->GetEntry(entry);
                        Hist_Dark_TelElevAzim.Fill(TelElevation,TelAzimuth);
                }

                TTree* Dark_tree = (TTree*) input_file->Get("run_"+TString(run_number)+"/stereo/data_off");
                //TTree* Dark_tree = (TTree*) input_file->Get("run_"+TString(run_number)+"/stereo/data_on");
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
                        if (!QualitySelection()) continue;
                        if (FoV()) {
                            if (SignalSelectionMSCW()) Hist_Dark_SR_MSCW.at(e).Fill(MSCW);
                            if (ControlSelectionMSCW()) Hist_Dark_CR_MSCW.at(e).Fill(MSCW);
                            if (Aux1SignalSelectionMSCW()) Hist_Dark_ASR1_MSCW.at(e).Fill(MSCW);
                            if (Aux1ControlSelectionMSCW()) Hist_Dark_ACR1_MSCW.at(e).Fill(MSCW);
                            if (Aux2SignalSelectionMSCW()) Hist_Dark_ASR2_MSCW.at(e).Fill(MSCW);
                            if (Aux2ControlSelectionMSCW()) Hist_Dark_ACR2_MSCW.at(e).Fill(MSCW);
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
                Target_tree->SetBranchAddress("ra",&ra_sky);
                Target_tree->SetBranchAddress("dec",&dec_sky);
                for (int entry=0;entry<Target_tree->GetEntries();entry++) {
                        Target_tree->GetEntry(entry);
                        int e = Hist_Target_SR_ErecS.at(0).FindBin(ErecS*1000.)-1;
                        if (e>=N_energy_bins) continue;
                        if (e<0) continue;
                        if (!QualitySelection()) continue;
                        if (MSCL<1) {
                            double depth_mean = 5.+log2(pow(ErecS*1000./0.08,0.4));
                            if (theta2<0.1) {
                                if (MSCW<1) {
                                        Hist_Target_ON_ElevDepth.Fill(TelElevation,SlantDepth*100./37.-depth_mean);
                                        Hist_Target_ON_ErecSDepth.Fill(ErecS*1000.,SlantDepth*100./37.);
                                }
                                Hist_Target_ON_MSCW_Alpha.Fill(MSCW);
                            }
                            else {
                                if (MSCW<1) {
                                        Hist_Target_OFF_ElevDepth.Fill(TelElevation,SlantDepth*100./37.-depth_mean);
                                        Hist_Target_OFF_ErecSDepth.Fill(ErecS*1000.,SlantDepth*100./37.);
                                }
                                Hist_Target_OFF_MSCW_Alpha.Fill(MSCW);
                            }
                        }
                        if (FoV()) {
                                if (SignalSelectionMSCW()) Hist_Target_SR_MSCW.at(e).Fill(MSCW);
                                if (ControlSelectionMSCW()) Hist_Target_CR_MSCW.at(e).Fill(MSCW);
                                if (Aux1SignalSelectionMSCW()) Hist_Target_ASR1_MSCW.at(e).Fill(MSCW);
                                if (Aux1ControlSelectionMSCW()) Hist_Target_ACR1_MSCW.at(e).Fill(MSCW);
                                if (Aux2SignalSelectionMSCW()) Hist_Target_ASR2_MSCW.at(e).Fill(MSCW);
                                if (Aux2ControlSelectionMSCW()) Hist_Target_ACR2_MSCW.at(e).Fill(MSCW);
                                if (SignalSelectionMSCW()) {
                                    if (MSCW>MSCW_cut_lower && MSCW<MSCW_cut_upper) {
                                        Hist_Target_SR_Energy.Fill(ErecS*1000.);
                                    }
                                    if (MSCW>Norm_Lower) {
                                        Hist_Target_CR_Energy.Fill(ErecS*1000.);
                                    }
                                }
                        }
                        else if (RingFoV()) {
                                if (SignalSelectionMSCW()) Hist_Target_Ring_MSCW.at(e).Fill(MSCW);
                        }
                        //if (FoV()) {
                        if (SignalSelectionMSCW()) {
                            if (MSCW>MSCW_cut_lower && MSCW<MSCW_cut_upper) {
                                Hist_Target_SR_SkyMap.at(e).Fill(ra_sky,dec_sky);
                                Hist_Target_SR_Theta2.at(e).Fill(theta2);
                            }
                            if (MSCW>Norm_Lower) {
                                Hist_Target_CR_SkyMap.at(e).Fill(ra_sky,dec_sky);
                                Hist_Target_CR_Theta2.at(e).Fill(theta2);
                            }
                        }
                        //}
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
                if (hist_error/hist_integral>0.05) {
                //if (bin_width/rms<0.2 || hist_error/hist_integral>0.2) {
                        Hist_Dark_SR_MSCW.at(e).Rebin(2);
                        Hist_Dark_SRB_MSCW.at(e).Rebin(2);
                        Hist_Dark_CR_MSCW.at(e).Rebin(2);
                        Hist_Dark_ASR1_MSCW.at(e).Rebin(2);
                        Hist_Dark_ACR1_MSCW.at(e).Rebin(2);
                        Hist_Dark_ASR2_MSCW.at(e).Rebin(2);
                        Hist_Dark_ACR2_MSCW.at(e).Rebin(2);
                        Hist_Dark_Bkg_MSCW.at(e).Rebin(2);
                        Hist_Dark_BkgTemp_MSCW.at(e).Rebin(2);
                        Hist_Dark_ABkg1_MSCW.at(e).Rebin(2);
                        Hist_Dark_ABkg1Temp_MSCW.at(e).Rebin(2);
                        Hist_Dark_ABkg2_MSCW.at(e).Rebin(2);
                        Hist_Dark_ABkg2Temp_MSCW.at(e).Rebin(2);
                        Hist_Dark_Deconv_MSCW.at(e).Rebin(2);
                        Hist_Dark_TrueDeconv_MSCW.at(e).Rebin(2);
                        Hist_Dark_Elec_MSCW.at(e).Rebin(2);
                        Hist_Target_SR_MSCW.at(e).Rebin(2);
                        Hist_Target_SRB_MSCW.at(e).Rebin(2);
                        Hist_Target_CR_MSCW.at(e).Rebin(2);
                        Hist_Target_ASR1_MSCW.at(e).Rebin(2);
                        Hist_Target_ACR1_MSCW.at(e).Rebin(2);
                        Hist_Target_ASR2_MSCW.at(e).Rebin(2);
                        Hist_Target_ACR2_MSCW.at(e).Rebin(2);
                        Hist_Target_Bkg_MSCW.at(e).Rebin(2);
                        Hist_Target_BkgTemp_MSCW.at(e).Rebin(2);
                        Hist_Target_ABkg1_MSCW.at(e).Rebin(2);
                        Hist_Target_ABkg1Temp_MSCW.at(e).Rebin(2);
                        Hist_Target_ABkg2_MSCW.at(e).Rebin(2);
                        Hist_Target_ABkg2Temp_MSCW.at(e).Rebin(2);
                        Hist_Target_Deconv_MSCW.at(e).Rebin(2);
                        Hist_Target_TrueDeconv_MSCW.at(e).Rebin(2);
                        Hist_Target_Elec_MSCW.at(e).Rebin(2);
                        Hist_Target_Ring_MSCW.at(e).Rebin(2);
                        Hist_Target_TotalBkg_MSCW.at(e).Rebin(2);
                }
                else break;
            }
        }

        vector<int> N_iter;
        vector<double> N_rms;
        vector<double> N_mean;
        vector<double> N_shift;
        for (int e=0;e<N_energy_bins;e++) {
                N_iter.push_back(50);
                //if (energy_bins[e]>=400) N_iter.at(e) = 100;
                //if (energy_bins[e]>=600) N_iter.at(e) = 100;
                //if (energy_bins[e]>=800) N_iter.at(e) = 100;
                //if (energy_bins[e]>=1000) N_iter.at(e) = 100;
                //if (energy_bins[e]>=1200) N_iter.at(e) = 100;
                //if (energy_bins[e]>=1600) N_iter.at(e) = 100;
                //if (energy_bins[e]>=2000) N_iter.at(e) = 50;
                //if (energy_bins[e]>=3000) N_iter.at(e) = 50;
                //if (energy_bins[e]>=4000) N_iter.at(e) = 50;
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

            double n_iter_final = N_iter.at(e);
            double chi2_best = 0.;
            for (int delta_n_iter = 10;delta_n_iter<=100;delta_n_iter++) {
                  int n_iter = delta_n_iter;
                  //int n_iter = N_iter.at(e);
                  double offset_begin = 0;
                  double chi2 = 0;
                  Deconvolution(&Hist_Dark_CR_MSCW.at(e),&Hist_Dark_SR_MSCW.at(e),&Hist_Dark_Deconv_MSCW.at(e),n_iter);
                  double rms_begin = Hist_Dark_Deconv_MSCW.at(e).GetRMS();
                  double mean_begin = Hist_Dark_Deconv_MSCW.at(e).GetMean();
                  double rms = rms_begin;
                  double mean = mean_begin;
                  myfunc->SetParameters(10.,mean,rms);
                  Hist_Dark_Deconv_MSCW.at(e).Reset();
                  Hist_Dark_Deconv_MSCW.at(e).FillRandom("myfunc",Hist_Dark_SR_MSCW.at(e).Integral()*100);
                  Deconvolution(&Hist_Dark_CR_MSCW.at(e),&Hist_Dark_Deconv_MSCW.at(e),&Hist_Dark_BkgTemp_MSCW.at(e),n_iter);
                  offset_begin = Hist_Dark_SR_MSCW.at(e).GetMean()-Hist_Dark_BkgTemp_MSCW.at(e).GetMean();
                  offset_begin = ShiftAndNormalize(&Hist_Dark_SR_MSCW.at(e),&Hist_Dark_BkgTemp_MSCW.at(e),&Hist_Dark_Bkg_MSCW.at(e),offset_begin,-1.,Norm_Upper,true);
                  chi2 += GetChi2(&Hist_Dark_SR_MSCW.at(e), &Hist_Dark_Bkg_MSCW.at(e),-1.,Norm_Upper,true);
                  if (chi2_best<chi2) {
                      chi2_best = chi2;
                      N_rms.at(e) = rms;
                      N_mean.at(e) = mean;
                      n_iter_final = n_iter;
                      N_iter.at(e) = n_iter;
                  } 
            }
            std::cout << "Dark, e " << energy_bins[e] << ", initial mean = " << N_mean.at(e) << std::endl;
            std::cout << "Dark, e " << energy_bins[e] << ", initial rms = " << N_rms.at(e) << std::endl;
            std::cout << "Dark, e " << energy_bins[e] << ", n_iter_final = " << n_iter_final << std::endl;
            for (int n_rms = 0; n_rms<=20;n_rms++) {
                double offset_begin = 0;
                double chi2 = 0;
                double rms = N_rms.at(e)-0.5*N_rms.at(e)+double(n_rms)*1.0*N_rms.at(e)/20.;
                double mean = N_mean.at(e);
                myfunc->SetParameters(10.,mean,rms);
                Hist_Dark_Deconv_MSCW.at(e).Reset();
                Hist_Dark_Deconv_MSCW.at(e).FillRandom("myfunc",Hist_Dark_SR_MSCW.at(e).Integral()*100);
                Deconvolution(&Hist_Dark_CR_MSCW.at(e),&Hist_Dark_Deconv_MSCW.at(e),&Hist_Dark_BkgTemp_MSCW.at(e),n_iter_final);
                offset_begin = Hist_Dark_SR_MSCW.at(e).GetMean()-Hist_Dark_BkgTemp_MSCW.at(e).GetMean();
                offset_begin = ShiftAndNormalize(&Hist_Dark_SR_MSCW.at(e),&Hist_Dark_BkgTemp_MSCW.at(e),&Hist_Dark_Bkg_MSCW.at(e),offset_begin,-1.,Norm_Upper,true);
                chi2 += GetChi2(&Hist_Dark_SR_MSCW.at(e), &Hist_Dark_Bkg_MSCW.at(e),-1.,Norm_Upper,true);
                if (chi2_best<chi2) {
                    chi2_best = chi2;
                    N_rms.at(e) = rms;
                    N_mean.at(e) = mean;
                } 
            }
            double rms = N_rms.at(e);
            double mean = N_mean.at(e);
            myfunc->SetParameters(10.,mean,rms);
            Hist_Dark_Deconv_MSCW.at(e).Reset();
            Hist_Dark_Deconv_MSCW.at(e).FillRandom("myfunc",Hist_Dark_SR_MSCW.at(e).Integral()*100);
            Deconvolution(&Hist_Dark_ACR1_MSCW.at(e),&Hist_Dark_Deconv_MSCW.at(e),&Hist_Dark_ABkg1Temp_MSCW.at(e),n_iter_final);
            Deconvolution(&Hist_Dark_ACR2_MSCW.at(e),&Hist_Dark_Deconv_MSCW.at(e),&Hist_Dark_ABkg2Temp_MSCW.at(e),n_iter_final);
            Deconvolution(&Hist_Dark_CR_MSCW.at(e),&Hist_Dark_Deconv_MSCW.at(e),&Hist_Dark_BkgTemp_MSCW.at(e),n_iter_final);
            Deconvolution(&Hist_Dark_CR_MSCW.at(e),&Hist_Dark_SR_MSCW.at(e),&Hist_Dark_TrueDeconv_MSCW.at(e),n_iter_final);
            std::cout << "Dark, e " << energy_bins[e] << ", final mean = " << Hist_Dark_Deconv_MSCW.at(e).GetMean() << std::endl;
            std::cout << "Dark, e " << energy_bins[e] << ", final rms = " << Hist_Dark_Deconv_MSCW.at(e).GetRMS() << std::endl;
            std::cout << "Dark, e " << energy_bins[e] << ", true Mean = " << Hist_Dark_TrueDeconv_MSCW.at(e).GetMean() << std::endl;
            std::cout << "Dark, e " << energy_bins[e] << ", true RMS = " << Hist_Dark_TrueDeconv_MSCW.at(e).GetRMS() << std::endl;
            double offset_begin = 0;
            offset_begin = Hist_Dark_SR_MSCW.at(e).GetMean()-Hist_Dark_BkgTemp_MSCW.at(e).GetMean();
            std::cout << "offset_begin = " << offset_begin << std::endl;
            offset_begin = ShiftAndNormalize(&Hist_Dark_ASR1_MSCW.at(e),&Hist_Dark_ABkg1Temp_MSCW.at(e),&Hist_Dark_ABkg1_MSCW.at(e),offset_begin,-1.,Norm_Upper,true);
            offset_begin = ShiftAndNormalize(&Hist_Dark_ASR2_MSCW.at(e),&Hist_Dark_ABkg2Temp_MSCW.at(e),&Hist_Dark_ABkg2_MSCW.at(e),offset_begin,-1.,Norm_Upper,true);
            offset_begin = Hist_Dark_SR_MSCW.at(e).GetMean()-Hist_Dark_BkgTemp_MSCW.at(e).GetMean();
            std::cout << "Dark, e " << energy_bins[e] << ", offset_begin = " << offset_begin << std::endl;
            offset_begin = ShiftAndNormalize(&Hist_Dark_SR_MSCW.at(e),&Hist_Dark_BkgTemp_MSCW.at(e),&Hist_Dark_Bkg_MSCW.at(e),offset_begin,-1.,Norm_Upper,true);
            N_shift.at(e) = offset_begin;
            std::cout << "offset_begin = " << offset_begin << std::endl;
        }

        std::cout << "Target run deconvolution... " << std::endl;
        // Target deconvolution method
        for (int e=0;e<N_energy_bins;e++) {
            std::cout << "=================================================================" << std::endl;
            std::cout << "Target, e " << energy_bins[e] << ", N_iter.at(e) = " << N_iter.at(e) << std::endl;

            double n_iter_final = N_iter.at(e);
            double chi2_best = 0.;
            for (int delta_n_iter = 10;delta_n_iter<=100;delta_n_iter++) {
                  int n_iter = delta_n_iter;
                  //int n_iter = N_iter.at(e);
                  double offset_begin = N_shift.at(e);
                  double chi2 = 0;
                  double rms_begin = N_rms.at(e);
                  double mean_begin = N_mean.at(e);
                  double rms = rms_begin;
                  double mean = mean_begin;
                  myfunc->SetParameters(10.,mean,rms);
                  Hist_Target_Deconv_MSCW.at(e).Reset();
                  Hist_Target_Deconv_MSCW.at(e).FillRandom("myfunc",Hist_Target_SR_MSCW.at(e).Integral()*100);
                  Deconvolution(&Hist_Target_CR_MSCW.at(e),&Hist_Target_Deconv_MSCW.at(e),&Hist_Target_BkgTemp_MSCW.at(e),n_iter);
                  //offset_begin = Hist_Target_SR_MSCW.at(e).GetMean()-Hist_Target_BkgTemp_MSCW.at(e).GetMean();
                  offset_begin = ShiftAndNormalize(&Hist_Target_SR_MSCW.at(e),&Hist_Target_BkgTemp_MSCW.at(e),&Hist_Target_Bkg_MSCW.at(e),offset_begin,Norm_Lower,Norm_Upper,true);
                  chi2 += GetChi2(&Hist_Target_SR_MSCW.at(e), &Hist_Target_Bkg_MSCW.at(e),-1.,Norm_Upper,false);
                  if (chi2_best<chi2) {
                      chi2_best = chi2;
                      N_rms.at(e) = rms;
                      N_mean.at(e) = mean;
                      n_iter_final = n_iter;
                  } 
            }
            std::cout << "Target, e " << energy_bins[e] << ", initial mean = " << N_mean.at(e) << std::endl;
            std::cout << "Target, e " << energy_bins[e] << ", initial rms = " << N_rms.at(e) << std::endl;
            std::cout << "Target, e " << energy_bins[e] << ", n_iter_final = " << n_iter_final << std::endl;
            for (int n_rms = 0; n_rms<=100;n_rms++) {
                double offset_begin = N_shift.at(e);
                double chi2 = 0;
                double rms = N_rms.at(e)-0.5*N_rms.at(e)+double(n_rms)*1.0*N_rms.at(e)/100.;
                double mean = N_mean.at(e);
                myfunc->SetParameters(10.,mean,rms);
                Hist_Target_Deconv_MSCW.at(e).Reset();
                Hist_Target_Deconv_MSCW.at(e).FillRandom("myfunc",Hist_Target_SR_MSCW.at(e).Integral()*100);
                Deconvolution(&Hist_Target_CR_MSCW.at(e),&Hist_Target_Deconv_MSCW.at(e),&Hist_Target_BkgTemp_MSCW.at(e),n_iter_final);
                //offset_begin = Hist_Target_SR_MSCW.at(e).GetMean()-Hist_Target_BkgTemp_MSCW.at(e).GetMean();
                offset_begin = ShiftAndNormalize(&Hist_Target_SR_MSCW.at(e),&Hist_Target_BkgTemp_MSCW.at(e),&Hist_Target_Bkg_MSCW.at(e),offset_begin,Norm_Lower,Norm_Upper,true);
                chi2 += GetChi2(&Hist_Target_SR_MSCW.at(e), &Hist_Target_Bkg_MSCW.at(e),-1.,Norm_Upper,false);
                if (chi2_best<chi2) {
                    chi2_best = chi2;
                    N_rms.at(e) = rms;
                    N_mean.at(e) = mean;
                } 
            }
            double rms = N_rms.at(e);
            double mean = N_mean.at(e);
            myfunc->SetParameters(10.,mean,rms);
            Hist_Target_Deconv_MSCW.at(e).Reset();
            Hist_Target_Deconv_MSCW.at(e).FillRandom("myfunc",Hist_Target_SR_MSCW.at(e).Integral()*100);
            Deconvolution(&Hist_Target_ACR1_MSCW.at(e),&Hist_Target_Deconv_MSCW.at(e),&Hist_Target_ABkg1Temp_MSCW.at(e),n_iter_final);
            Deconvolution(&Hist_Target_ACR2_MSCW.at(e),&Hist_Target_Deconv_MSCW.at(e),&Hist_Target_ABkg2Temp_MSCW.at(e),n_iter_final);
            Deconvolution(&Hist_Target_CR_MSCW.at(e),&Hist_Target_Deconv_MSCW.at(e),&Hist_Target_BkgTemp_MSCW.at(e),n_iter_final);
            Deconvolution(&Hist_Target_CR_MSCW.at(e),&Hist_Target_SR_MSCW.at(e),&Hist_Target_TrueDeconv_MSCW.at(e),n_iter_final);
            std::cout << "Target, e " << energy_bins[e] << ", final mean = " << Hist_Target_Deconv_MSCW.at(e).GetMean() << std::endl;
            std::cout << "Target, e " << energy_bins[e] << ", final rms = " << Hist_Target_Deconv_MSCW.at(e).GetRMS() << std::endl;
            std::cout << "Target, e " << energy_bins[e] << ", true Mean = " << Hist_Target_TrueDeconv_MSCW.at(e).GetMean() << std::endl;
            std::cout << "Target, e " << energy_bins[e] << ", true RMS = " << Hist_Target_TrueDeconv_MSCW.at(e).GetRMS() << std::endl;
            double offset_begin = N_shift.at(e);
            //offset_begin = Hist_Target_SR_MSCW.at(e).GetMean()-Hist_Target_BkgTemp_MSCW.at(e).GetMean();
            std::cout << "offset_begin = " << offset_begin << std::endl;
            offset_begin = ShiftAndNormalize(&Hist_Target_ASR1_MSCW.at(e),&Hist_Target_ABkg1Temp_MSCW.at(e),&Hist_Target_ABkg1_MSCW.at(e),offset_begin,Norm_Lower,Norm_Upper,true);
            offset_begin = ShiftAndNormalize(&Hist_Target_ASR2_MSCW.at(e),&Hist_Target_ABkg2Temp_MSCW.at(e),&Hist_Target_ABkg2_MSCW.at(e),offset_begin,Norm_Lower,Norm_Upper,true);
            //offset_begin = Hist_Target_SR_MSCW.at(e).GetMean()-Hist_Target_BkgTemp_MSCW.at(e).GetMean();
            offset_begin = N_shift.at(e);
            std::cout << "Target, e " << energy_bins[e] << ", offset_begin = " << offset_begin << std::endl;
            offset_begin = ShiftAndNormalize(&Hist_Target_SR_MSCW.at(e),&Hist_Target_BkgTemp_MSCW.at(e),&Hist_Target_Bkg_MSCW.at(e),offset_begin,Norm_Lower,Norm_Upper,true);
            std::cout << "offset_begin = " << offset_begin << std::endl;
        }


        // simple dark template method
        std::cout << "Simple dark template... " << std::endl;
        for (int e=0;e<N_energy_bins;e++) {
                int norm_bin_low_target = Hist_Target_Bkg_MSCW.at(e).FindBin(Norm_Lower);
                int norm_bin_up_target = Hist_Target_Bkg_MSCW.at(e).FindBin(Norm_Upper);
                double scale_target = Hist_Target_SR_MSCW.at(e).Integral(norm_bin_low_target,norm_bin_up_target)/Hist_Dark_SR_MSCW.at(e).Integral(norm_bin_low_target,norm_bin_up_target);
                Hist_Dark_SR_MSCW.at(e).Scale(scale_target);
                Hist_Dark_SRB_MSCW.at(e).Scale(scale_target);
                Hist_Dark_Bkg_MSCW.at(e).Scale(scale_target);
        }


        // Get electron bkg
        std::cout << "Getting electrons... " << std::endl;
        for (int e=0;e<N_energy_bins;e++) {

                Hist_Dark_Elec_MSCW.at(e).Reset();
                Hist_Dark_Elec_MSCW.at(e).Add(&Hist_Dark_SR_MSCW.at(e));
                Hist_Dark_Elec_MSCW.at(e).Add(&Hist_Dark_Bkg_MSCW.at(e),-1.);
                for (int i=0;i<Hist_Target_SR_MSCW.at(e).GetNbinsX();i++) {
                        if (Hist_Dark_SR_MSCW.at(e).GetBinCenter(i+1)>Norm_Lower) 
                            Hist_Dark_Elec_MSCW.at(e).SetBinContent(i+1,0.);
                            Hist_Dark_Elec_MSCW.at(e).SetBinError(i+1,0.);
                }
                //Hist_Dark_Elec_MSCW.at(e).Divide(&Hist_Dark_Bkg_MSCW.at(e));
                //Hist_Target_Elec_MSCW.at(e).Reset();
                //Hist_Target_Elec_MSCW.at(e).Add(&Hist_Target_Bkg_MSCW.at(e));
                //Hist_Target_Elec_MSCW.at(e).Multiply(&Hist_Dark_Elec_MSCW.at(e));
                //Hist_Dark_Elec_MSCW.at(e).Multiply(&Hist_Dark_Bkg_MSCW.at(e));
                //Hist_Target_TotalBkg_MSCW.at(e).Reset();
                //Hist_Target_TotalBkg_MSCW.at(e).Add(&Hist_Target_Bkg_MSCW.at(e));
                //Hist_Target_TotalBkg_MSCW.at(e).Add(&Hist_Target_Elec_MSCW.at(e));
                double e2p_ratio = Hist_Dark_Elec_MSCW.at(e).Integral()/Hist_Dark_Bkg_MSCW.at(e).Integral();
                double target_elec_total = e2p_ratio*Hist_Target_Bkg_MSCW.at(e).Integral();
                Hist_Target_Elec_MSCW.at(e).Reset();
                Hist_Target_Elec_MSCW.at(e).Add(&Hist_Dark_Elec_MSCW.at(e));
                double scale = target_elec_total/Hist_Dark_Elec_MSCW.at(e).Integral();
                Hist_Target_Elec_MSCW.at(e).Scale(scale);
                Hist_Target_TotalBkg_MSCW.at(e).Reset();
                Hist_Target_TotalBkg_MSCW.at(e).Add(&Hist_Target_Bkg_MSCW.at(e));
                Hist_Target_TotalBkg_MSCW.at(e).Add(&Hist_Target_Elec_MSCW.at(e));
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

        // Get sky map bkg
        for (int e=0;e<N_energy_bins;e++) {
            double ratio, ratio_err;
            std::tie(ratio,ratio_err) = RatioSRCR(&Hist_Target_TotalBkg_MSCW.at(e), MSCW_cut_lower, MSCW_cut_upper, Norm_Lower);
            Hist_Target_Bkg_Theta2.at(e).Reset();
            for (int b=0;b<Hist_Target_Bkg_Theta2.at(e).GetNbinsX();b++) {
                double NCR = Hist_Target_CR_Theta2.at(e).GetBinContent(b+1);
                double NCR_err = Hist_Target_CR_Theta2.at(e).GetBinError(b+1);
                double NSR = NCR*ratio;
                double NSR_err = NCR*ratio*pow(pow(NCR_err/NCR,2)+pow(ratio_err/ratio,2),0.5);
                Hist_Target_Bkg_Theta2.at(e).SetBinContent(b+1,NSR);
                Hist_Target_Bkg_Theta2.at(e).SetBinError(b+1,NSR_err);
            }
            Hist_Target_Bkg_SkyMap.at(e).Reset();
            for (int bx=0;bx<Hist_Target_CR_SkyMap.at(e).GetNbinsX();bx++) {
                for (int by=0;by<Hist_Target_CR_SkyMap.at(e).GetNbinsY();by++) {
                    double NCR = Hist_Target_CR_SkyMap.at(e).GetBinContent(bx+1,by+1);
                    double NCR_err = Hist_Target_CR_SkyMap.at(e).GetBinError(bx+1,by+1);
                    double NSR = NCR*ratio;
                    double NSR_err = NCR*ratio*pow(pow(NCR_err/NCR,2)+pow(ratio_err/ratio,2),0.5);
                    Hist_Target_Bkg_SkyMap.at(e).SetBinContent(bx+1,by+1,NSR);
                    Hist_Target_Bkg_SkyMap.at(e).SetBinError(bx+1,by+1,NSR_err);
                }
            }
        }
        Hist_Target_Bkg_Energy.Reset();
        for (int b=0;b<Hist_Target_SR_Energy.GetNbinsX();b++) {
            int e = Hist_Target_SR_ErecS.at(0).FindBin(Hist_Target_SR_Energy.GetBinCenter(b+1))-1;
            if (e<0) continue;
            if (e>=Hist_Target_SR_ErecS.at(0).GetNbinsX()) continue;
            std::cout << "e = " << e << std::endl;
            std::cout << "Hist_Target_SR_Energy.GetBinCenter(b+1) = " << Hist_Target_SR_Energy.GetBinCenter(b+1) << std::endl;
            double ratio, ratio_err;
            std::tie(ratio,ratio_err) = RatioSRCR(&Hist_Target_TotalBkg_MSCW.at(e), MSCW_cut_lower, MSCW_cut_upper, Norm_Lower);
            double NCR = Hist_Target_CR_Energy.GetBinContent(b+1);
            double NCR_err = Hist_Target_CR_Energy.GetBinError(b+1);
            double NSR = NCR*ratio;
            double NSR_err = NCR*ratio*pow(pow(NCR_err/NCR,2)+pow(ratio_err/ratio,2),0.5);
            Hist_Target_Bkg_Energy.SetBinContent(b+1,NSR);
            Hist_Target_Bkg_Energy.SetBinError(b+1,NSR_err);
        }

        int norm_bin_low_ring = Hist_Target_ON_MSCW_Alpha.FindBin(Norm_Lower);
        int norm_bin_up_ring = Hist_Target_ON_MSCW_Alpha.FindBin(Norm_Upper);
        double scale_ring = Hist_Target_ON_MSCW_Alpha.Integral(norm_bin_low_ring,norm_bin_up_ring)/Hist_Target_OFF_MSCW_Alpha.Integral(norm_bin_low_ring,norm_bin_up_ring);
        Hist_Target_OFF_ErecSDepth.Scale(scale_ring);
        Hist_Target_OFF_ElevDepth.Scale(scale_ring);
        Hist_Target_ON_ErecSDepth.Add(&Hist_Target_OFF_ErecSDepth,-1.);
        Hist_Target_ON_ElevDepth.Add(&Hist_Target_OFF_ElevDepth,-1.);

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
        Hist_Target_ON_ElevDepth.Write();
        Hist_Target_OFF_ElevDepth.Write();
        Hist_Dark_TelElevAzim.Write();
        Hist_Target_SR_Energy.Write();
        Hist_Target_CR_Energy.Write();
        Hist_Target_Bkg_Energy.Write();
        for (int e=0;e<N_energy_bins;e++) {
                Hist_Dark_ASR1_MSCW.at(e).Write();
                Hist_Dark_ACR1_MSCW.at(e).Write();
                Hist_Dark_ABkg1_MSCW.at(e).Write();
                Hist_Dark_ASR2_MSCW.at(e).Write();
                Hist_Dark_ACR2_MSCW.at(e).Write();
                Hist_Dark_ABkg2_MSCW.at(e).Write();
                Hist_Dark_SR_MSCW.at(e).Write();
                Hist_Dark_SRB_MSCW.at(e).Write();
                Hist_Dark_CR_MSCW.at(e).Write();
                Hist_Dark_Bkg_MSCW.at(e).Write();
                Hist_Dark_Deconv_MSCW.at(e).Write();
                Hist_Dark_TrueDeconv_MSCW.at(e).Write();
                Hist_Dark_Elec_MSCW.at(e).Write();
                Hist_Target_ASR1_MSCW.at(e).Write();
                Hist_Target_ACR1_MSCW.at(e).Write();
                Hist_Target_ABkg1_MSCW.at(e).Write();
                Hist_Target_ASR2_MSCW.at(e).Write();
                Hist_Target_ACR2_MSCW.at(e).Write();
                Hist_Target_ABkg2_MSCW.at(e).Write();
                Hist_Target_SR_MSCW.at(e).Write();
                Hist_Target_SRB_MSCW.at(e).Write();
                Hist_Target_CR_MSCW.at(e).Write();
                Hist_Target_Bkg_MSCW.at(e).Write();
                Hist_Target_Deconv_MSCW.at(e).Write();
                Hist_Target_TrueDeconv_MSCW.at(e).Write();
                Hist_Target_Elec_MSCW.at(e).Write();
                Hist_Target_Ring_MSCW.at(e).Write();
                Hist_Target_TotalBkg_MSCW.at(e).Write();
                Hist_Target_SR_SkyMap.at(e).Write();
                Hist_Target_CR_SkyMap.at(e).Write();
                Hist_Target_Bkg_SkyMap.at(e).Write();
                Hist_Target_SR_Theta2.at(e).Write();
                Hist_Target_CR_Theta2.at(e).Write();
                Hist_Target_Bkg_Theta2.at(e).Write();
        }
        OutputFile.Close();

}
