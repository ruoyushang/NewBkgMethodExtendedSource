
import sys
import runlist
from runlist import *
import EmissionHeightMethodConfig
from EmissionHeightMethodConfig import *
#import MSCWMethodConfig
#from MSCWMethodConfig import *

ROOT.gStyle.SetOptStat(0)
ROOT.TH1.SetDefaultSumw2()
ROOT.gStyle.SetPaintTextFormat("0.3f")

NoSignalAssumption = False

target = sys.argv[1]

Elev_cut_lower = Elev_Bin[0]
Elev_cut_upper = Elev_Bin[len(Elev_Bin)-1]
MSCL_cut_lower = -2
MSCL_cut_upper = 0.7
MSCW_cut_lower = -2
MSCW_cut_upper = 0.7
#MSCL_cut_lower = 2.
#MSCL_cut_upper = 3.
#MSCW_cut_lower = 2.
#MSCW_cut_upper = 3.
Erec_cut_lower = 0.1
Erec_cut_upper = 1000

Hist_Gamma_Attenu = []
Hist_CR_Attenu = []
Hist_CR_Efficiency_list = []
Hist_Inflation = []
Hist_EmissionHeight_Source = TH1D("Hist_EmissionHeight_Source","",40,0,40)
Hist_EmissionHeight_Target = TH1D("Hist_EmissionHeight_Target","",40,0,40)
Hist_MSCW_Source = TH1D("Hist_MSCW_Source","",60,-5,10)
Hist_MSCW_Target = TH1D("Hist_MSCW_Target","",60,-5,10)

signal_width = 1000.0
control_width = 4.0

observation_time = 0
time_cut = 1e10

def AzimuthConverter(azimuth):

    new_azimuth = abs(azimuth-180.)
    return float(new_azimuth)

def ControlRegionHeightSelection(tree,elev):
    
    cut_upper = 10.+control_width
    cut_lower = 10.-control_width

    if (tree.EmissionHeight<cut_upper and tree.EmissionHeight>cut_lower): return False
    #if (tree.EmissionHeight>(cut_upper+1.) or tree.EmissionHeight<(cut_lower-1.)): return False
    return True

def SignalRegionHeightSelection(tree,elev):
    
    cut_upper = 10.+signal_width
    cut_lower = 10.-signal_width

    if (tree.EmissionHeight>cut_upper or tree.EmissionHeight<cut_lower): return False
    return True

def CommonSelection(tree):
    if (tree.MSCL<MSCL_cut_lower): return False
    if (tree.MSCL>MSCL_cut_upper): return False
    if (tree.Erec<Erec_cut_lower): return False
    if (tree.Erec>Erec_cut_upper): return False
    return True

def MakeATag():
    global tag
    tag = ''
    tag += '_Elev%sto%s'%(Elev_cut_lower,Elev_cut_upper)
    tag += '_MSCL%sto%s'%(MSCL_cut_lower,MSCL_cut_upper)
    tag += '_Erec%0.1fto%0.1f'%(Erec_cut_lower,Erec_cut_upper)
    #tag += '_Time%s'%(time_cut/(60.*60.))
    tag += '_'+Height_cut
    tag += '_'+RW_method
    tag += '_'+Fake_Removal
    tag += '_'+target+'_'+target_field
    tag += '_'+tag_method


def GetAnttenuationHistogram(source,region,run,data_type):

    if data_type=='anasum':
        if source == 'Crab':
            anasum_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/Crab_V6_Moderate-TMVA-BDT.RB.%s.root"%(run))
        if source == '2ndCrab':
            anasum_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/Crab_V6_Moderate-TMVA-BDT.RB.%s.root"%(run))
        if source == 'PKS1424':
            anasum_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/PKS1424_V6_Moderate-TMVA-BDT.RB.%s.root"%(run))
        if source == '3C264':
            anasum_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/3C264_V6_Moderate-TMVA-BDT.RB.%s.root"%(run))
        if source == 'H1426':
            anasum_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/H1426_V6_Moderate-TMVA-BDT.RB.%s.root"%(run))
        if region == 'on':
            anasum_tree = anasum_file.Get("run_%s/stereo/data_on"%(run))
            ring_tree = anasum_file.Get("run_%s/stereo/data_off"%(run))
        else:
            anasum_tree = anasum_file.Get("run_%s/stereo/data_off"%(run))
        pointing_tree = anasum_file.Get("run_%s/stereo/pointingDataReduced"%(run))
    else:
        anasum_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/analysis/Results/%s.mscw.root"%(run))
        anasum_tree = anasum_file.Get("data")
        pointing_tree = anasum_file.Get("data")

    pointing_tree.GetEntry(0)
    time_begin = pointing_tree.Time
    pointing_tree.GetEntry(pointing_tree.GetEntries()-1)
    time_end = pointing_tree.Time
    hist_Ele_vs_time = ROOT.TProfile("hist_Ele_vs_time","",100,time_begin,time_end+0.01*(time_end-time_begin),0,100)
    hist_Azi_vs_time = ROOT.TProfile("hist_Azi_vs_time","",100,time_begin,time_end+0.01*(time_end-time_begin),-400,400)
    hist_RA_vs_time = ROOT.TProfile("hist_RA_vs_time","",100,time_begin,time_end+0.01*(time_end-time_begin),0,100)
    hist_Dec_vs_time = ROOT.TProfile("hist_Dec_vs_time","",100,time_begin,time_end+0.01*(time_end-time_begin),0,100)
    for i in range(0,pointing_tree.GetEntries()):
        pointing_tree.GetEntry(i)
        hist_Ele_vs_time.Fill(pointing_tree.Time,pointing_tree.TelElevation)
        hist_Azi_vs_time.Fill(pointing_tree.Time,AzimuthConverter(pointing_tree.TelAzimuth))
        hist_RA_vs_time.Fill(pointing_tree.Time,pointing_tree.TelRAJ2000)
        hist_Dec_vs_time.Fill(pointing_tree.Time,pointing_tree.TelDecJ2000)

    for i in range(0,anasum_tree.GetEntries()):
        anasum_tree.GetEntry(i)
        #if not CommonSelection(anasum_tree): continue
        if (anasum_tree.MSCL<MSCL_cut_lower): continue
        if (anasum_tree.MSCL>MSCL_cut_upper): continue
        if (anasum_tree.MSCW<MSCW_cut_lower): continue
        if (anasum_tree.MSCW>MSCW_cut_upper and anasum_tree.MSCW<3.): continue # MSCW>3 is for normalization of off region
        if (anasum_tree.Erec<Erec_cut_lower): continue
        if (anasum_tree.Erec>Erec_cut_upper): continue
        elevation = hist_Ele_vs_time.GetBinContent(hist_Ele_vs_time.FindBin(anasum_tree.Time))
        azimuth = hist_Azi_vs_time.GetBinContent(hist_Azi_vs_time.FindBin(anasum_tree.Time))
        RA = hist_RA_vs_time.GetBinContent(hist_RA_vs_time.FindBin(anasum_tree.Time))
        Dec = hist_Dec_vs_time.GetBinContent(hist_Dec_vs_time.FindBin(anasum_tree.Time))
        if elevation<Elev_cut_lower or elevation>Elev_cut_upper: continue
        scale_elev = Hist_TelElev_Counts_scale.GetBinContent(Hist_TelElev_Counts_scale.FindBin(elevation))
        scale_azim = Hist_TelAzim_Counts_scale.GetBinContent(Hist_TelAzim_Counts_scale.FindBin(AzimuthConverter(azimuth)))
        scale = scale_elev
        if RW_method=='NoRW': 
            scale = 1
        if region == 'off':
            Hist_EmissionHeight_Source.Fill(anasum_tree.EmissionHeight,scale)
            Hist_TelElev_Counts_source_Sum.Fill(elevation,scale)
            Hist_TelAzim_Counts_source_Sum.Fill(AzimuthConverter(azimuth),scale)
        if ControlRegionHeightSelection(anasum_tree,elevation):
            if region == 'on':
                Hist_OnData_Control.Fill(anasum_tree.MSCW,scale)
            else:
                Hist_OffData_Control.Fill(anasum_tree.MSCW,scale)
        if SignalRegionHeightSelection(anasum_tree,elevation):
            if region == 'on':
                Hist_OnData_Signal.Fill(anasum_tree.MSCW,scale)
            else:
                Hist_OffData_Signal.Fill(anasum_tree.MSCW,scale)

def GetSourceElevation(run):

    if source == 'Crab':
        anasum_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/Crab_V6_Moderate-TMVA-BDT.RB.%s.root"%(run))
    pointing_tree = anasum_file.Get("run_%s/stereo/pointingDataReduced"%(run))

    for i in range(0,pointing_tree.GetEntries()):
        pointing_tree.GetEntry(i)
        Hist_TelElev_Counts_source.Fill(pointing_tree.TelElevation)
        Hist_TelAzim_Counts_source.Fill(AzimuthConverter(pointing_tree.TelAzimuth))

def GetTargetElevation(run):

    if target == '2ndCrab':
        target_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/Crab_V6_Moderate-TMVA-BDT.RB.%s.root"%(run))
    if target == 'PKS1424':
        target_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/PKS1424_V6_Moderate-TMVA-BDT.RB.%s.root"%(run))
    if target == '3C264':
        target_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/3C264_V6_Moderate-TMVA-BDT.RB.%s.root"%(run))
    if target == 'H1426':
        target_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/H1426_V6_Moderate-TMVA-BDT.RB.%s.root"%(run))
    if target == 'Ton599':
        target_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/Ton599_V6_Moderate-TMVA-BDT.RB.%s.root"%(run))
    if target == 'IC443':
        target_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/IC443_V6_Moderate-TMVA-BDT.RB.%s.root"%(run))
    target_pointing_tree = target_file.Get("run_%s/stereo/pointingDataReduced"%(run))

    for i in range(0,target_pointing_tree.GetEntries()):
        target_pointing_tree.GetEntry(i)
        Hist_TelElev_Counts_target.Fill(target_pointing_tree.TelElevation)
        Hist_TelAzim_Counts_target.Fill(AzimuthConverter(target_pointing_tree.TelAzimuth))


def MeasureAttenuation():

    Rate_CR_Attenuation = 0
    Err_Rate_CR_Attenuation = 0
    Rate_Gamma_Attenuation = 0
    Err_Rate_Gamma_Attenuation = 0
    
    if source == 'Crab':
        for run in Crab_runlist:
            GetSourceElevation(run)
    if target == '2ndCrab':
        for run in SecondCrab_runlist:
            GetTargetElevation(run)
    if target == 'PKS1424':
        for run in PKS1424_runlist:
            GetTargetElevation(run)
    if target == '3C264':
        for run in runlist_3C264:
            GetTargetElevation(run)
    if target == 'H1426':
        for run in runlist_H1426:
            GetTargetElevation(run)
    if target == 'Ton599':
        for run in runlist_Ton599:
            GetTargetElevation(run)
    if target == 'IC443':
        for run in runlist_IC443:
            GetTargetElevation(run)
    #print Hist_TelElev_Counts_target.Print('All')
    for b in range(0,Hist_TelElev_Counts_scale.GetNbinsX()):
        Hist_TelElev_Counts_scale.SetBinContent(b+1,Hist_TelElev_Counts_target.GetBinContent(b+1))
    Hist_TelElev_Counts_scale.Divide(Hist_TelElev_Counts_source)
    for b in range(0,Hist_TelAzim_Counts_scale.GetNbinsX()):
        Hist_TelAzim_Counts_scale.SetBinContent(b+1,Hist_TelAzim_Counts_target.GetBinContent(b+1))
    Hist_TelAzim_Counts_scale.Divide(Hist_TelAzim_Counts_source)
    #print Hist_TelElev_Counts_scale.Print('All')
    if source == 'Crab':
        for run in Crab_runlist:
            GetAnttenuationHistogram('Crab','on',run,'anasum')
            GetAnttenuationHistogram('Crab','off',run,'anasum')
    if source == '2ndCrab':
        for run in SecondCrab_runlist:
            GetAnttenuationHistogram('Crab','on',run,'anasum')
            GetAnttenuationHistogram('Crab','off',run,'anasum')
    Hist_TelElev_Counts_scale.Reset()
    Hist_TelElev_Counts_source.Reset()
    Hist_TelElev_Counts_target.Reset()
    Hist_TelAzim_Counts_scale.Reset()
    Hist_TelAzim_Counts_source.Reset()
    Hist_TelAzim_Counts_target.Reset()
    
    print 'Hist_OnData_Signal.GetEntries() = %s'%(Hist_OnData_Signal.GetEntries())
    print 'Hist_OnData_Control.GetEntries() = %s'%(Hist_OnData_Control.GetEntries())
    print 'Hist_OffData_Signal.GetEntries() = %s'%(Hist_OffData_Signal.GetEntries())
    print 'Hist_OffData_Control.GetEntries() = %s'%(Hist_OffData_Control.GetEntries())
    print 'Hist_OnData_Signal.Integral() = %s'%(Hist_OnData_Signal.GetBinContent(Hist_OnData_Signal.FindBin(0)))
    print 'Hist_OnData_Control.Integral() = %s'%(Hist_OnData_Control.GetBinContent(Hist_OnData_Signal.FindBin(0)))
    print 'Hist_OffData_Signal.Integral() = %s'%(Hist_OffData_Signal.GetBinContent(Hist_OnData_Signal.FindBin(0)))
    print 'Hist_OffData_Control.Integral() = %s'%(Hist_OffData_Control.GetBinContent(Hist_OnData_Signal.FindBin(0)))

    #if Hist_OffData_Signal.GetBinContent(len(MSCW_Bin)-1)<1:
    #    return 0,0,0,0
    #if Hist_OffData_Control.GetBinContent(len(MSCW_Bin)-1)<1:
    #    return 0,0,0,0
    Normalization_Factor_Signal = Hist_OnData_Signal.GetBinContent(len(MSCW_Bin)-1)
    Normalization_Factor_Signal = Normalization_Factor_Signal/Hist_OffData_Signal.GetBinContent(len(MSCW_Bin)-1)
    Hist_OffData_Signal.Scale(Normalization_Factor_Signal)
    Hist_OffData_Control.Scale(Normalization_Factor_Signal)
    print 'Normalization_Factor_Signal = %s'%(Normalization_Factor_Signal)
    #Normalization_Factor_Control = Hist_OnData_Control.GetBinContent(len(MSCW_Bin)-1)
    #Normalization_Factor_Control = Normalization_Factor_Control/Hist_OffData_Control.GetBinContent(len(MSCW_Bin)-1)
    #Hist_OffData_Signal.Scale(Normalization_Factor_Control)
    #Hist_OffData_Control.Scale(Normalization_Factor_Control)
    
    for b in range(0,Hist_CR_Attenuation.GetNbinsX()):
        Hist_CR_Attenuation.SetBinContent(b+1,Hist_OffData_Control.GetBinContent(b+1))
        Hist_CR_Attenuation.SetBinError(b+1,Hist_OffData_Control.GetBinError(b+1))
    Hist_CR_Attenuation.Divide(Hist_OffData_Signal)
    Rate_CR_Attenuation = Hist_CR_Attenuation.GetBinContent(Hist_CR_Attenuation.FindBin(0))
    Err_Rate_CR_Attenuation = Hist_CR_Attenuation.GetBinError(Hist_CR_Attenuation.FindBin(0))
    print 'Rate_CR_Attenuation = %s +/- %s'%(Rate_CR_Attenuation,Err_Rate_CR_Attenuation)
    
    for b in range(0,Hist_Gamma_Signal.GetNbinsX()):
        Hist_Gamma_Signal.SetBinContent(b+1,Hist_OnData_Signal.GetBinContent(b+1))
        Hist_Gamma_Signal.SetBinError(b+1,Hist_OnData_Signal.GetBinError(b+1))
    Hist_Gamma_Signal.Add(Hist_OffData_Signal,-1)
    for b in range(0,Hist_Gamma_Control.GetNbinsX()):
        Hist_Gamma_Control.SetBinContent(b+1,Hist_OnData_Control.GetBinContent(b+1))
        Hist_Gamma_Control.SetBinError(b+1,Hist_OnData_Control.GetBinError(b+1))
    Hist_Gamma_Control.Add(Hist_OffData_Control,-1)
    for b in range(0,Hist_Gamma_Attenuation.GetNbinsX()):
        Hist_Gamma_Attenuation.SetBinContent(b+1,Hist_Gamma_Control.GetBinContent(b+1))
        Hist_Gamma_Attenuation.SetBinError(b+1,Hist_Gamma_Control.GetBinError(b+1))
    Hist_Gamma_Attenuation.Divide(Hist_Gamma_Signal)
    Rate_Gamma_Attenuation = Hist_Gamma_Attenuation.GetBinContent(Hist_Gamma_Attenuation.FindBin(0))
    Err_Rate_Gamma_Attenuation = Hist_Gamma_Attenuation.GetBinError(Hist_Gamma_Attenuation.FindBin(0))
    print 'Rate_Gamma_Attenuation = %s +/- %s'%(Rate_Gamma_Attenuation,Err_Rate_Gamma_Attenuation)

    Hist_OnData_Signal.Reset()
    Hist_OnData_Control.Reset()
    Hist_OffData_Signal.Reset()
    Hist_OffData_Control.Reset()
    Hist_Gamma_Signal.Reset()
    Hist_Gamma_Control.Reset()
    Hist_CR_Attenuation.Reset()
    Hist_Gamma_Attenuation.Reset()

    return Rate_CR_Attenuation,Err_Rate_CR_Attenuation,Rate_Gamma_Attenuation,Err_Rate_Gamma_Attenuation

def GenerateExtendedSignalMSCW(source,region,run,data_type,scale):

    if data_type=='anasum':
        if source == 'Crab':
            anasum_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/Crab_V6_Moderate-TMVA-BDT.RB.%s.root"%(run))
        if source == '2ndCrab':
            anasum_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/Crab_V6_Moderate-TMVA-BDT.RB.%s.root"%(run))
        anasum_tree = anasum_file.Get("run_%s/stereo/data_on"%(run))
        pointing_tree = anasum_file.Get("run_%s/stereo/pointingDataReduced"%(run))
    else:
        anasum_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/analysis/Results/%s.mscw.root"%(run))
        anasum_tree = anasum_file.Get("data")
        pointing_tree = anasum_file.Get("pointingDataReduced")

    pointing_tree.GetEntry(0)
    time_begin = pointing_tree.Time
    pointing_tree.GetEntry(pointing_tree.GetEntries()-1)
    time_end = pointing_tree.Time
    hist_Ele_vs_time = ROOT.TProfile("hist_Ele_vs_time","",100,time_begin,time_end+0.01*(time_end-time_begin),0,100)
    hist_Azi_vs_time = ROOT.TProfile("hist_Azi_vs_time","",100,time_begin,time_end+0.01*(time_end-time_begin),-400,400)
    hist_RA_vs_time = ROOT.TProfile("hist_RA_vs_time","",100,time_begin,time_end+0.01*(time_end-time_begin),0,100)
    hist_Dec_vs_time = ROOT.TProfile("hist_Dec_vs_time","",100,time_begin,time_end+0.01*(time_end-time_begin),0,100)
    for i in range(0,pointing_tree.GetEntries()):
        pointing_tree.GetEntry(i)
        hist_Ele_vs_time.Fill(pointing_tree.Time,pointing_tree.TelElevation)
        hist_Azi_vs_time.Fill(pointing_tree.Time,AzimuthConverter(pointing_tree.TelAzimuth))
        hist_RA_vs_time.Fill(pointing_tree.Time,pointing_tree.TelRAJ2000)
        hist_Dec_vs_time.Fill(pointing_tree.Time,pointing_tree.TelDecJ2000)

    for i in range(0,anasum_tree.GetEntries()):
        anasum_tree.GetEntry(i)
        if not CommonSelection(anasum_tree): continue
        if region == 'on':
            if (anasum_tree.theta2>0.01): continue
        if region == 'off':
            if (anasum_tree.theta2<0.1): continue
            if (anasum_tree.theta2>1.0): continue
        elevation = hist_Ele_vs_time.GetBinContent(hist_Ele_vs_time.FindBin(anasum_tree.Time))
        azimuth = hist_Azi_vs_time.GetBinContent(hist_Azi_vs_time.FindBin(anasum_tree.Time))
        RA = hist_RA_vs_time.GetBinContent(hist_RA_vs_time.FindBin(anasum_tree.Time))
        Dec = hist_Dec_vs_time.GetBinContent(hist_Dec_vs_time.FindBin(anasum_tree.Time))
        if elevation<Elev_cut_lower or elevation>Elev_cut_upper: continue
        elev_weight = Hist_TelElev_Counts_scale.GetBinContent(Hist_TelElev_Counts_scale.FindBin(elevation))
        azim_weight = Hist_TelAzim_Counts_scale.GetBinContent(Hist_TelAzim_Counts_scale.FindBin(AzimuthConverter(azimuth)))
        pointing_weight = elev_weight
        if RW_method=='NoRW': 
            pointing_weight = 1
        signal_ra = MyRandom.Gaus(0.,0.2)
        signal_dec = MyRandom.Gaus(0.,0.2)
        signal_theta2 = signal_ra*signal_ra+signal_dec*signal_dec
        if ControlRegionMSCWSelection(anasum_tree,elevation):
            if (signal_theta2<FOV*FOV):
                Hist_Erec_CR_Raw_AddSignal.Fill(min(energy_bins[len(energy_bins)-2],anasum_tree.Erec),pointing_weight*scale)
            Hist2D_Xoff_vs_Yoff_CR_Raw_AddSignal.Fill(signal_ra,signal_dec,pointing_weight*scale)
            Hist_theta2_CR_Raw_AddSignal.Fill(signal_theta2,pointing_weight*scale)
            Hist_theta2_zoomin_CR_Raw_AddSignal.Fill(signal_theta2,pointing_weight*scale)
        if SignalRegionMSCWSelection(anasum_tree,elevation):
            if (signal_theta2<FOV*FOV):
                Hist_Erec_Data_AddSignal.Fill(min(energy_bins[len(energy_bins)-2],anasum_tree.Erec),pointing_weight*scale)
            Hist2D_Xoff_vs_Yoff_Data_AddSignal.Fill(signal_ra,signal_dec,pointing_weight*scale)
            Hist_theta2_Data_AddSignal.Fill(signal_theta2,pointing_weight*scale)
            Hist_theta2_zoomin_Data_AddSignal.Fill(signal_theta2,pointing_weight*scale)

def GenerateExtendedSignal(source,region,run,data_type,scale):

    if data_type=='anasum':
        if source == 'Crab':
            anasum_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/Crab_V6_Moderate-TMVA-BDT.RB.%s.root"%(run))
        if source == '2ndCrab':
            anasum_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/Crab_V6_Moderate-TMVA-BDT.RB.%s.root"%(run))
        anasum_tree = anasum_file.Get("run_%s/stereo/data_on"%(run))
        pointing_tree = anasum_file.Get("run_%s/stereo/pointingDataReduced"%(run))
    else:
        anasum_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/analysis/Results/%s.mscw.root"%(run))
        anasum_tree = anasum_file.Get("data")
        pointing_tree = anasum_file.Get("pointingDataReduced")

    pointing_tree.GetEntry(0)
    time_begin = pointing_tree.Time
    pointing_tree.GetEntry(pointing_tree.GetEntries()-1)
    time_end = pointing_tree.Time
    hist_Ele_vs_time = ROOT.TProfile("hist_Ele_vs_time","",100,time_begin,time_end+0.01*(time_end-time_begin),0,100)
    hist_Azi_vs_time = ROOT.TProfile("hist_Azi_vs_time","",100,time_begin,time_end+0.01*(time_end-time_begin),-400,400)
    hist_RA_vs_time = ROOT.TProfile("hist_RA_vs_time","",100,time_begin,time_end+0.01*(time_end-time_begin),0,100)
    hist_Dec_vs_time = ROOT.TProfile("hist_Dec_vs_time","",100,time_begin,time_end+0.01*(time_end-time_begin),0,100)
    for i in range(0,pointing_tree.GetEntries()):
        pointing_tree.GetEntry(i)
        hist_Ele_vs_time.Fill(pointing_tree.Time,pointing_tree.TelElevation)
        hist_Azi_vs_time.Fill(pointing_tree.Time,AzimuthConverter(pointing_tree.TelAzimuth))
        hist_RA_vs_time.Fill(pointing_tree.Time,pointing_tree.TelRAJ2000)
        hist_Dec_vs_time.Fill(pointing_tree.Time,pointing_tree.TelDecJ2000)

    for i in range(0,anasum_tree.GetEntries()):
        anasum_tree.GetEntry(i)
        if not CommonSelection(anasum_tree): continue
        if region == 'on':
            if (anasum_tree.theta2>0.01): continue
        if region == 'off':
            if (anasum_tree.theta2<0.1): continue
            if (anasum_tree.theta2>1.0): continue
        elevation = hist_Ele_vs_time.GetBinContent(hist_Ele_vs_time.FindBin(anasum_tree.Time))
        azimuth = hist_Azi_vs_time.GetBinContent(hist_Azi_vs_time.FindBin(anasum_tree.Time))
        RA = hist_RA_vs_time.GetBinContent(hist_RA_vs_time.FindBin(anasum_tree.Time))
        Dec = hist_Dec_vs_time.GetBinContent(hist_Dec_vs_time.FindBin(anasum_tree.Time))
        if elevation<Elev_cut_lower or elevation>Elev_cut_upper: continue
        elev_weight = Hist_TelElev_Counts_scale.GetBinContent(Hist_TelElev_Counts_scale.FindBin(elevation))
        azim_weight = Hist_TelAzim_Counts_scale.GetBinContent(Hist_TelAzim_Counts_scale.FindBin(AzimuthConverter(azimuth)))
        pointing_weight = elev_weight
        if RW_method=='NoRW': 
            pointing_weight = 1
        signal_ra = MyRandom.Gaus(0.,0.2)
        signal_dec = MyRandom.Gaus(0.,0.2)
        signal_theta2 = signal_ra*signal_ra+signal_dec*signal_dec
        if ControlRegionHeightSelection(anasum_tree,elevation):
            if (anasum_tree.MSCW>MSCW_cut_lower and anasum_tree.MSCW<MSCW_cut_upper):
                if (signal_theta2<FOV*FOV):
                    Hist_Erec_CR_Raw_AddSignal.Fill(min(energy_bins[len(energy_bins)-2],anasum_tree.Erec),pointing_weight*scale)
                Hist2D_Xoff_vs_Yoff_CR_Raw_AddSignal.Fill(signal_ra,signal_dec,pointing_weight*scale)
                Hist_theta2_CR_Raw_AddSignal.Fill(signal_theta2,pointing_weight*scale)
                Hist_theta2_zoomin_CR_Raw_AddSignal.Fill(signal_theta2,pointing_weight*scale)
        if SignalRegionHeightSelection(anasum_tree,elevation):
            if (anasum_tree.MSCW>MSCW_cut_lower and anasum_tree.MSCW<MSCW_cut_upper):
                if (signal_theta2<FOV*FOV):
                    Hist_Erec_Data_AddSignal.Fill(min(energy_bins[len(energy_bins)-2],anasum_tree.Erec),pointing_weight*scale)
                Hist2D_Xoff_vs_Yoff_Data_AddSignal.Fill(signal_ra,signal_dec,pointing_weight*scale)
                Hist_theta2_Data_AddSignal.Fill(signal_theta2,pointing_weight*scale)
                Hist_theta2_zoomin_Data_AddSignal.Fill(signal_theta2,pointing_weight*scale)

def SingleRunAnalysisMSCW(source,region,run,data_type):

    global observation_time
    global N_GammaMC_Fake
    global N_GammaData_Fake
    global N_GammaMC_ES
    global N_GammaData_ES

    if data_type=='anasum':
        if source == 'Crab':
            anasum_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/Crab_V6_Moderate-TMVA-BDT.RB.%s.root"%(run))
        if source == '2ndCrab':
            anasum_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/Crab_V6_Moderate-TMVA-BDT.RB.%s.root"%(run))
        if source == 'PKS1424':
            anasum_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/PKS1424_V6_Moderate-TMVA-BDT.RB.%s.root"%(run))
        if source == '3C264':
            anasum_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/3C264_V6_Moderate-TMVA-BDT.RB.%s.root"%(run))
        if source == 'H1426':
            anasum_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/H1426_V6_Moderate-TMVA-BDT.RB.%s.root"%(run))
        if source == 'Ton599':
            anasum_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/Ton599_V6_Moderate-TMVA-BDT.RB.%s.root"%(run))
        if source == 'IC443':
            anasum_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/IC443_V6_Moderate-TMVA-BDT.RB.%s.root"%(run))
        anasum_tree = anasum_file.Get("run_%s/stereo/data_on"%(run))
        ring_tree = anasum_file.Get("run_%s/stereo/data_off"%(run))
        pointing_tree = anasum_file.Get("run_%s/stereo/pointingDataReduced"%(run))
    else:
        anasum_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/analysis/Results/%s.mscw.root"%(run))
        anasum_tree = anasum_file.Get("data")
        pointing_tree = anasum_file.Get("pointingDataReduced")

    if not source=="MC":
        pointing_tree.GetEntry(0)
        time_begin = pointing_tree.Time
        pointing_tree.GetEntry(pointing_tree.GetEntries()-1)
        time_end = pointing_tree.Time
        hist_Ele_vs_time = ROOT.TProfile("hist_Ele_vs_time","",100,time_begin,time_end+0.01*(time_end-time_begin),0,100)
        hist_Azi_vs_time = ROOT.TProfile("hist_Azi_vs_time","",100,time_begin,time_end+0.01*(time_end-time_begin),-400,400)
        hist_RA_vs_time = ROOT.TProfile("hist_RA_vs_time","",100,time_begin,time_end+0.01*(time_end-time_begin),0,100)
        hist_Dec_vs_time = ROOT.TProfile("hist_Dec_vs_time","",100,time_begin,time_end+0.01*(time_end-time_begin),0,100)
        for i in range(0,pointing_tree.GetEntries()):
            pointing_tree.GetEntry(i)
            hist_Ele_vs_time.Fill(pointing_tree.Time,pointing_tree.TelElevation)
            hist_Azi_vs_time.Fill(pointing_tree.Time,AzimuthConverter(pointing_tree.TelAzimuth))
            hist_RA_vs_time.Fill(pointing_tree.Time,pointing_tree.TelRAJ2000)
            hist_Dec_vs_time.Fill(pointing_tree.Time,pointing_tree.TelDecJ2000)

        observation_time_this_run = 0
        #print 'time cut = %s s'%(time_cut)
        for i in range(0,anasum_tree.GetEntries()):
            anasum_tree.GetEntry(i)
            observation_time_this_run = anasum_tree.Time-time_begin
            if time_cut<(observation_time_this_run+observation_time): continue
            elevation = hist_Ele_vs_time.GetBinContent(hist_Ele_vs_time.FindBin(anasum_tree.Time))
            azimuth = hist_Azi_vs_time.GetBinContent(hist_Azi_vs_time.FindBin(anasum_tree.Time))
            RA = hist_RA_vs_time.GetBinContent(hist_RA_vs_time.FindBin(anasum_tree.Time))
            Dec = hist_Dec_vs_time.GetBinContent(hist_Dec_vs_time.FindBin(anasum_tree.Time))
            if not CommonSelection(anasum_tree): continue
            if region=="off":
                if anasum_tree.theta2<0.1: continue
            if (anasum_tree.MSCW>MSCW_cut_lower and anasum_tree.MSCW<MSCW_cut_upper):
                Hist2D_TelElev_vs_EmissionHeight_Data_All.Fill(elevation,anasum_tree.EmissionHeight)
            if elevation<Elev_cut_lower or elevation>Elev_cut_upper: continue
            if (anasum_tree.MSCW>MSCW_cut_lower and anasum_tree.MSCW<MSCW_cut_upper):
                Hist2D_Erec_vs_EmissionHeight_Data_All.Fill(anasum_tree.Erec,anasum_tree.EmissionHeight)
            Hist_MSCW_Target.Fill(anasum_tree.MSCW)
            Hist_TelElev_Counts_target_Sum.Fill(elevation)
            Hist_TelAzim_Counts_target_Sum.Fill(AzimuthConverter(azimuth))
            scale = 1
            Xoff_derot = anasum_tree.Xoff_derot
            Yoff_derot = anasum_tree.Yoff_derot
            Xoff = anasum_tree.Xoff
            Yoff = anasum_tree.Yoff
            Xcore = anasum_tree.Xcore
            Ycore = anasum_tree.Ycore
            if ControlRegionMSCWSelection(anasum_tree,elevation):
                if (anasum_tree.theta2<FOV*FOV):
                    Hist_Erec_CR_Raw.Fill(min(energy_bins[len(energy_bins)-2],anasum_tree.Erec),1.)
		Hist2D_Xoff_vs_Yoff_CR_Raw.Fill(anasum_tree.ra,anasum_tree.dec)
		Hist_theta2_CR_Raw.Fill(anasum_tree.theta2)
		Hist_theta2_zoomin_CR_Raw.Fill(anasum_tree.theta2)
            if SignalRegionMSCWSelection(anasum_tree,elevation):
                if (anasum_tree.theta2<FOV*FOV):
                    Hist_Erec_Data.Fill(min(energy_bins[len(energy_bins)-2],anasum_tree.Erec))
                elif (anasum_tree.theta2>Ring_Inner*Ring_Inner and anasum_tree.theta2<Ring_Outer*Ring_Outer):
                    Hist_Erec_Ring.Fill(min(energy_bins[len(energy_bins)-2],anasum_tree.Erec))
		Hist2D_Xoff_vs_Yoff_Data.Fill(anasum_tree.ra,anasum_tree.dec)
		Hist_theta2_Data.Fill(anasum_tree.theta2)
		Hist_theta2_zoomin_Data.Fill(anasum_tree.theta2)
            if (anasum_tree.MSCW>3.):
                if (anasum_tree.theta2<FOV*FOV):
                    Hist_Norm_Data.Fill(0)
                elif (anasum_tree.theta2>Ring_Inner*Ring_Inner and anasum_tree.theta2<Ring_Outer*Ring_Outer):
                    Hist_Norm_Ring.Fill(0)

    observation_time += observation_time_this_run
    #print 'current observation time = %s min.'%(observation_time/60.)

def SingleRunAnalysis(source,region,run,data_type):

    global observation_time
    global N_GammaMC_Fake
    global N_GammaData_Fake
    global N_GammaMC_ES
    global N_GammaData_ES

    if data_type=='anasum':
        if source == 'Crab':
            anasum_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/Crab_V6_Moderate-TMVA-BDT.RB.%s.root"%(run))
        if source == '2ndCrab':
            anasum_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/Crab_V6_Moderate-TMVA-BDT.RB.%s.root"%(run))
        if source == 'PKS1424':
            anasum_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/PKS1424_V6_Moderate-TMVA-BDT.RB.%s.root"%(run))
        if source == '3C264':
            anasum_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/3C264_V6_Moderate-TMVA-BDT.RB.%s.root"%(run))
        if source == 'H1426':
            anasum_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/H1426_V6_Moderate-TMVA-BDT.RB.%s.root"%(run))
        if source == 'Ton599':
            anasum_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/Ton599_V6_Moderate-TMVA-BDT.RB.%s.root"%(run))
        if source == 'IC443':
            anasum_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/IC443_V6_Moderate-TMVA-BDT.RB.%s.root"%(run))
        anasum_tree = anasum_file.Get("run_%s/stereo/data_on"%(run))
        ring_tree = anasum_file.Get("run_%s/stereo/data_off"%(run))
        pointing_tree = anasum_file.Get("run_%s/stereo/pointingDataReduced"%(run))
    else:
        anasum_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/analysis/Results/%s.mscw.root"%(run))
        anasum_tree = anasum_file.Get("data")
        pointing_tree = anasum_file.Get("pointingDataReduced")

    if not source=="MC":
        pointing_tree.GetEntry(0)
        time_begin = pointing_tree.Time
        pointing_tree.GetEntry(pointing_tree.GetEntries()-1)
        time_end = pointing_tree.Time
        hist_Ele_vs_time = ROOT.TProfile("hist_Ele_vs_time","",100,time_begin,time_end+0.01*(time_end-time_begin),0,100)
        hist_Azi_vs_time = ROOT.TProfile("hist_Azi_vs_time","",100,time_begin,time_end+0.01*(time_end-time_begin),-400,400)
        hist_RA_vs_time = ROOT.TProfile("hist_RA_vs_time","",100,time_begin,time_end+0.01*(time_end-time_begin),0,100)
        hist_Dec_vs_time = ROOT.TProfile("hist_Dec_vs_time","",100,time_begin,time_end+0.01*(time_end-time_begin),0,100)
        for i in range(0,pointing_tree.GetEntries()):
            pointing_tree.GetEntry(i)
            hist_Ele_vs_time.Fill(pointing_tree.Time,pointing_tree.TelElevation)
            hist_Azi_vs_time.Fill(pointing_tree.Time,AzimuthConverter(pointing_tree.TelAzimuth))
            hist_RA_vs_time.Fill(pointing_tree.Time,pointing_tree.TelRAJ2000)
            hist_Dec_vs_time.Fill(pointing_tree.Time,pointing_tree.TelDecJ2000)

        observation_time_this_run = 0
        for i in range(0,anasum_tree.GetEntries()):
            anasum_tree.GetEntry(i)
            observation_time_this_run = anasum_tree.Time-time_begin
            if time_cut<(observation_time_this_run+observation_time): continue
            elevation = hist_Ele_vs_time.GetBinContent(hist_Ele_vs_time.FindBin(anasum_tree.Time))
            azimuth = hist_Azi_vs_time.GetBinContent(hist_Azi_vs_time.FindBin(anasum_tree.Time))
            RA = hist_RA_vs_time.GetBinContent(hist_RA_vs_time.FindBin(anasum_tree.Time))
            Dec = hist_Dec_vs_time.GetBinContent(hist_Dec_vs_time.FindBin(anasum_tree.Time))
            if not CommonSelection(anasum_tree): continue
            if region=="off":
                if anasum_tree.theta2<0.1: continue
            if (anasum_tree.MSCW>MSCW_cut_lower and anasum_tree.MSCW<MSCW_cut_upper):
                Hist2D_TelElev_vs_EmissionHeight_Data_All.Fill(elevation,anasum_tree.EmissionHeight)
            if elevation<Elev_cut_lower or elevation>Elev_cut_upper: continue
            if (anasum_tree.MSCW>MSCW_cut_lower and anasum_tree.MSCW<MSCW_cut_upper):
                Hist2D_Erec_vs_EmissionHeight_Data_All.Fill(anasum_tree.Erec,anasum_tree.EmissionHeight)
                Hist_EmissionHeight_Target.Fill(anasum_tree.EmissionHeight)
                Hist_TelElev_Counts_target_Sum.Fill(elevation)
                Hist_TelAzim_Counts_target_Sum.Fill(AzimuthConverter(azimuth))
            scale = 1
            Xoff_derot = anasum_tree.Xoff_derot
            Yoff_derot = anasum_tree.Yoff_derot
            Xcore = anasum_tree.Xcore
            Ycore = anasum_tree.Ycore
            if ControlRegionHeightSelection(anasum_tree,elevation):
                if (anasum_tree.MSCW>MSCW_cut_lower and anasum_tree.MSCW<MSCW_cut_upper):
                    if (anasum_tree.theta2<FOV*FOV):
                        Hist_Erec_CR_Raw.Fill(min(energy_bins[len(energy_bins)-2],anasum_tree.Erec),1.)
		    Hist2D_Xoff_vs_Yoff_CR_Raw.Fill(anasum_tree.ra,anasum_tree.dec)
		    Hist_theta2_CR_Raw.Fill(anasum_tree.theta2)
		    Hist_theta2_zoomin_CR_Raw.Fill(anasum_tree.theta2)
            if SignalRegionHeightSelection(anasum_tree,elevation):
                if (anasum_tree.MSCW>MSCW_cut_lower and anasum_tree.MSCW<MSCW_cut_upper):
                    if (anasum_tree.theta2<FOV*FOV):
                        Hist_Erec_Data.Fill(min(energy_bins[len(energy_bins)-2],anasum_tree.Erec))
                    elif (anasum_tree.theta2>Ring_Inner*Ring_Inner and anasum_tree.theta2<Ring_Outer*Ring_Outer):
                        Hist_Erec_Ring.Fill(min(energy_bins[len(energy_bins)-2],anasum_tree.Erec))
		    Hist2D_Xoff_vs_Yoff_Data.Fill(anasum_tree.ra,anasum_tree.dec)
		    Hist_theta2_Data.Fill(anasum_tree.theta2)
		    Hist_theta2_zoomin_Data.Fill(anasum_tree.theta2)
                if (anasum_tree.MSCW>3.):
                    if (anasum_tree.theta2<FOV*FOV):
                        Hist_Norm_Data.Fill(0)
                    elif (anasum_tree.theta2>Ring_Inner*Ring_Inner and anasum_tree.theta2<Ring_Outer*Ring_Outer):
                        Hist_Norm_Ring.Fill(0)

    observation_time += observation_time_this_run
    #print 'current observation time = %s min.'%(observation_time/60.)

    #if region == 'on' and data_type=='anasum':
    #    for i in range(0,ring_tree.GetEntries()):
    #        ring_tree.GetEntry(i)
    #        if not CommonSelection(ring_tree): continue
    #        elevation = hist_Ele_vs_time.GetBinContent(hist_Ele_vs_time.FindBin(ring_tree.Time))
    #        azimuth = hist_Azi_vs_time.GetBinContent(hist_Azi_vs_time.FindBin(ring_tree.Time))
    #        if elevation<Elev_cut_lower or elevation>Elev_cut_upper: continue
    #        if SignalRegionHeightSelection(ring_tree,elevation):
    #            if (ring_tree.MSCW<1.):
    #                Hist_Erec_Ring.Fill(min(energy_bins[len(energy_bins)-2],ring_tree.Erec))
    #            if (ring_tree.MSCW>3.):
    #                Hist_Norm_Ring.Fill(0)

def CalculateBkgMSCWMethod(Hist_CR,Hist_Data,Hist_CR_Raw,CR_attenu,err_CR_attenu):

    Hist_CR.Reset()
    Hist_CR.Add(Hist_CR_Raw,1./CR_attenu)
    for b in range(0,Hist_CR.GetNbinsX()):
        old_err = Hist_CR.GetBinError(b+1)
        add_err = Hist_CR.GetBinContent(b+1)*(err_CR_attenu/CR_attenu)
        Hist_CR.SetBinError(b+1,pow(old_err*old_err+add_err*add_err,0.5))
    for b in range(0,Hist_CR.GetNbinsX()):
        Hist_CR.SetBinContent(b+1,max(0,Hist_CR.GetBinContent(b+1)))

def CalculateBkgHeightMethod(Hist_CR,Hist_Data,Hist_CR_Raw,CR_attenu,err_CR_attenu,gamma_attenu,err_gamma_attenu):

    #N_signal (number of events in signal region) = N_egamma + N_CR
    #N_control (number of events in control region) = a*N_egamma + b*N_CR, a = gamma attenuation, b = CR attenuation
    #N_control - b*N_signal = (a-b)*N_egamma
    #N_egamma = (N_control-b*N_signal)/(a-b)
    #N_CR = N_signal-N_egamma

    if NoSignalAssumption: gamma_attenu = 0

    Hist_CR.Reset()
    Hist_CR.Add(Hist_Data)
    Hist_CR.Add(Hist_CR_Raw,-1./CR_attenu)
    for b in range(0,Hist_CR.GetNbinsX()):
        old_err = Hist_CR.GetBinError(b+1)
        add_err = Hist_CR_Raw.GetBinContent(b+1)*err_CR_attenu
        Hist_CR.SetBinError(b+1,pow(old_err*old_err+add_err*add_err,0.5))
    Hist_CR.Scale(1./(1.-gamma_attenu/CR_attenu))
    if not gamma_attenu==0:
        for b in range(0,Hist_CR.GetNbinsX()):
            old_err = Hist_CR.GetBinError(b+1)
            numerator = pow(err_gamma_attenu/gamma_attenu,2)+pow(err_CR_attenu/CR_attenu,2)
            numerator = pow(numerator,0.5)*gamma_attenu/CR_attenu
            add_err = Hist_CR_Raw.GetBinContent(b+1)*numerator
            Hist_CR.SetBinError(b+1,pow(old_err*old_err+add_err*add_err,0.5))
    Hist_CR.Scale(-1.)
    Hist_CR.Add(Hist_Data)
    for b in range(0,Hist_CR.GetNbinsX()):
        Hist_CR.SetBinContent(b+1,max(0,Hist_CR.GetBinContent(b+1)))

def RunExtendedSourceAnalysis(data_type,CR_attenu,err_CR_attenu,gamma_attenu,err_gamma_attenu):

    print 'Erec_cut_lower = %s, Erec_cut_upper = %s'%(Erec_cut_lower,Erec_cut_upper)
    print 'CR_attenu = %s +/- %s, gamma_attenu = %s +/- %s'%(CR_attenu,err_CR_attenu,gamma_attenu,err_gamma_attenu)
    if gamma_attenu<0: gamma_attenu = 0

    if target == 'Crab':
        for run in Crab_runlist:
            SingleRunAnalysis('Crab',target_field,run,data_type)
    if target == '2ndCrab':
        for run in SecondCrab_runlist:
            SingleRunAnalysis('Crab',target_field,run,data_type)
    if target == 'PKS1424':
        for run in PKS1424_runlist:
            SingleRunAnalysis('PKS1424',target_field,run,data_type)
    if target == '3C264':
        for run in runlist_3C264:
            SingleRunAnalysis('3C264',target_field,run,data_type)
    if target == 'H1426':
        for run in runlist_H1426:
            SingleRunAnalysis('H1426',target_field,run,data_type)
    if target == 'Ton599':
        for run in runlist_Ton599:
            SingleRunAnalysis('Ton599',target_field,run,data_type)
    if target == 'IC443':
        for run in runlist_IC443:
            SingleRunAnalysis('IC443',target_field,run,data_type)

    # Here we add signals using Crab data
    if add_signal:
        for run in SecondCrab_runlist:
            #GenerateExtendedSignal('Crab','on',run,data_type,1.0)
            GenerateExtendedSignal('Crab','off',run,data_type,1.0)
            #GenerateExtendedSignal('Crab','off',run,data_type,scale*Normalization_Factor_Signal)
    

    if data_type=='anasum':
        Normalization_Factor_Ring = 0
        if not Hist_Norm_Ring.Integral()==0:
            Normalization_Factor_Ring = Hist_Norm_Data.Integral()
            Normalization_Factor_Ring = Normalization_Factor_Ring/Hist_Norm_Ring.Integral()
        Hist_Erec_Ring.Scale(Normalization_Factor_Ring)

    print 'Hist_Erec_CR_Raw.Integral() = %s'%(Hist_Erec_CR_Raw.Integral())
    print 'Hist_Erec_Data.Integral() = %s'%(Hist_Erec_Data.Integral())
    print 'Hist_Erec_CR_Raw_AddSignal.Integral() = %s'%(Hist_Erec_CR_Raw_AddSignal.Integral())
    print 'Hist_Erec_Data_AddSignal.Integral() = %s'%(Hist_Erec_Data_AddSignal.Integral())
    if not Fake_Removal=='NoFakeRemoved':
        if not add_signal:
            CalculateBkgHeightMethod(Hist_Erec_CR,Hist_Erec_Data,Hist_Erec_CR_Raw,CR_attenu,err_CR_attenu,gamma_attenu,err_gamma_attenu)
        else:
            Hist_Erec_Data_AddSignal.Add(Hist_Erec_Data)
            Hist_Erec_CR_Raw_AddSignal.Add(Hist_Erec_CR_Raw)
            CalculateBkgHeightMethod(Hist_Erec_CR,Hist_Erec_Data_AddSignal,Hist_Erec_CR_Raw_AddSignal,CR_attenu,err_CR_attenu,gamma_attenu,err_gamma_attenu)
        if not add_signal:
            CalculateBkgHeightMethod(Hist2D_Xoff_vs_Yoff_CR,Hist2D_Xoff_vs_Yoff_Data,Hist2D_Xoff_vs_Yoff_CR_Raw,CR_attenu,err_CR_attenu,gamma_attenu,err_gamma_attenu)
        else:
            Hist2D_Xoff_vs_Yoff_Data_AddSignal.Add(Hist2D_Xoff_vs_Yoff_Data)
            Hist2D_Xoff_vs_Yoff_CR_Raw_AddSignal.Add(Hist2D_Xoff_vs_Yoff_CR_Raw)
            CalculateBkgHeightMethod(Hist2D_Xoff_vs_Yoff_CR,Hist2D_Xoff_vs_Yoff_Data_AddSignal,Hist2D_Xoff_vs_Yoff_CR_Raw_AddSignal,CR_attenu,err_CR_attenu,gamma_attenu,err_gamma_attenu)
        if not add_signal:
            CalculateBkgHeightMethod(Hist_theta2_CR,Hist_theta2_Data,Hist_theta2_CR_Raw,CR_attenu,err_CR_attenu,gamma_attenu,err_gamma_attenu)
        else:
            Hist_theta2_Data_AddSignal.Add(Hist_theta2_Data)
            Hist_theta2_CR_Raw_AddSignal.Add(Hist_theta2_CR_Raw)
            CalculateBkgHeightMethod(Hist_theta2_CR,Hist_theta2_Data_AddSignal,Hist_theta2_CR_Raw_AddSignal,CR_attenu,err_CR_attenu,gamma_attenu,err_gamma_attenu)
        if not add_signal:
            CalculateBkgHeightMethod(Hist_theta2_zoomin_CR,Hist_theta2_zoomin_Data,Hist_theta2_zoomin_CR_Raw,CR_attenu,err_CR_attenu,gamma_attenu,err_gamma_attenu)
        else:
            Hist_theta2_zoomin_Data_AddSignal.Add(Hist_theta2_zoomin_Data)
            Hist_theta2_zoomin_CR_Raw_AddSignal.Add(Hist_theta2_zoomin_CR_Raw)
            CalculateBkgHeightMethod(Hist_theta2_zoomin_CR,Hist_theta2_zoomin_Data_AddSignal,Hist_theta2_zoomin_CR_Raw_AddSignal,CR_attenu,err_CR_attenu,gamma_attenu,err_gamma_attenu)
    Err_N_CR = ROOT.Double(0.)
    N_CR = Hist_Erec_CR.IntegralAndError(1,Hist_Erec_CR.GetNbinsX(),Err_N_CR)
    print 'Hist_Erec_CR.Integral() = (CR_raw - gamma_att*data)/(CR_att-gamma_att) = %s +/- %s'%(N_CR,Err_N_CR)
    print 'Hist_Erec_Ring.Integral() = %s'%(Hist_Erec_Ring.Integral())
    #print Hist_theta2_CR_Raw.Print("All")
    #print Hist_theta2_CR.Print("All")
    
    
    Hist_Erec_Data_Sum.Add(Hist_Erec_Data)
    if add_signal: 
        Hist_Erec_Data_Sum_AddSignal.Add(Hist_Erec_Data_AddSignal)
    Hist_Erec_Ring_Sum.Add(Hist_Erec_Ring)
    Hist_Erec_CR_Sum.Add(Hist_Erec_CR)

    Hist2D_Xoff_vs_Yoff_Data_Sum.Add(Hist2D_Xoff_vs_Yoff_Data)
    if add_signal: 
        Hist2D_Xoff_vs_Yoff_Data_Sum_AddSignal.Add(Hist2D_Xoff_vs_Yoff_Data_AddSignal)
    Hist2D_Xoff_vs_Yoff_CR_Sum.Add(Hist2D_Xoff_vs_Yoff_CR)

    Hist_theta2_Data_Sum.Add(Hist_theta2_Data)
    if add_signal: 
        Hist_theta2_Data_Sum_AddSignal.Add(Hist_theta2_Data_AddSignal)
    Hist_theta2_CR_Sum.Add(Hist_theta2_CR)

    Hist_theta2_zoomin_Data_Sum.Add(Hist_theta2_zoomin_Data)
    if add_signal: 
        Hist_theta2_zoomin_Data_Sum_AddSignal.Add(Hist_theta2_zoomin_Data_AddSignal)
    Hist_theta2_zoomin_CR_Sum.Add(Hist_theta2_zoomin_CR)

    MakeStackPlot(Hist_theta2_Data,Hist_theta2_CR,Hist_theta2_Ring_Sum,'theta2')
    #MakeStackPlot(Hist_theta2_zoomin_Data,Hist_theta2_zoomin_CR,Hist_theta2_zoomin_Ring_Sum,'theta2_zoomin')

    Hist_Norm_Data.Reset()
    Hist_Norm_Ring.Reset()
    Hist_Erec_CR.Reset()
    Hist_Erec_CR_Raw.Reset()
    Hist_Erec_CR_Raw_AddSignal.Reset()
    Hist_Erec_Data.Reset()
    Hist_Erec_Data_AddSignal.Reset()
    Hist_Erec_Ring.Reset()
    Hist2D_Xoff_vs_Yoff_Data.Reset()
    Hist2D_Xoff_vs_Yoff_Data_AddSignal.Reset()
    Hist2D_Xoff_vs_Yoff_CR.Reset()
    Hist2D_Xoff_vs_Yoff_CR_Raw.Reset()
    Hist2D_Xoff_vs_Yoff_CR_Raw_AddSignal.Reset()
    Hist_theta2_Data.Reset()
    Hist_theta2_Data_AddSignal.Reset()
    Hist_theta2_CR.Reset()
    Hist_theta2_CR_Raw.Reset()
    Hist_theta2_CR_Raw_AddSignal.Reset()
    Hist_theta2_zoomin_Data.Reset()
    Hist_theta2_zoomin_Data_AddSignal.Reset()
    Hist_theta2_zoomin_CR.Reset()
    Hist_theta2_zoomin_CR_Raw.Reset()
    Hist_theta2_zoomin_CR_Raw_AddSignal.Reset()

def AttenuationRateAtDifferentHeight():

    global Erec_cut_lower
    global Erec_cut_upper
    global MSCL_cut_lower
    global MSCL_cut_upper
    global control_width
    global Hist_Gamma_Attenu
    global Hist_CR_Attenu
    global Hist_Inflation

    for energy in range(0,len(energy_bins)-1):
        Hist_Gamma_Attenu += [Hist_ControlWidth_Bins.Clone()]
        Hist_CR_Attenu += [Hist_ControlWidth_Bins.Clone()]
        Hist_Inflation += [Hist_ControlWidth_Bins.Clone()]
        for width in range(1,Hist_ControlWidth_Bins.GetNbinsX()+1):
            control_width = Hist_ControlWidth_Bins.GetBinLowEdge(width)
            Erec_cut_lower = energy_bins[energy]
            Erec_cut_upper = energy_bins[energy+1]
            #if energy==len(energy_bins)-2: Erec_cut_upper = 1e10
            MakeATag()
            print '++++++++++++++++++++++++++++++++++++++++++++++'
            print 'selection = %s, control region = %s'%(tag, control_width)
            CR_attenu,err_CR_attenu,gamma_attenu,err_gamma_attenu = MeasureAttenuation()
            if CR_attenu<0: CR_attenu = 0
            if gamma_attenu<0: gamma_attenu = 0
            Hist_Gamma_Attenu[energy].SetBinContent(Hist_ControlWidth_Bins.FindBin(control_width),gamma_attenu)
            Hist_Gamma_Attenu[energy].SetBinError(Hist_ControlWidth_Bins.FindBin(control_width),err_gamma_attenu)
            Hist_CR_Attenu[energy].SetBinContent(Hist_ControlWidth_Bins.FindBin(control_width),CR_attenu)
            Hist_CR_Attenu[energy].SetBinError(Hist_ControlWidth_Bins.FindBin(control_width),err_CR_attenu)
            inflation = pow(1+gamma_attenu/CR_attenu,0.5)/(1-gamma_attenu/CR_attenu)
            inflation = inflation/pow(CR_attenu,0.5)
            print 'inflation = %s'%(inflation)
            Hist_Inflation[energy].SetBinContent(Hist_ControlWidth_Bins.FindBin(control_width),inflation)

    canvas2 = ROOT.TCanvas("canvas2","canvas2", 200, 10, 600, 600)
    pad2 = ROOT.TPad("pad2","pad2",0.1,0.8,1,1)
    pad2.SetBottomMargin(0.0)
    pad2.SetTopMargin(0.03)
    pad2.SetBorderMode(1)
    pad1 = ROOT.TPad("pad1","pad1",0.1,0,1,0.8)
    pad1.SetBottomMargin(0.2)
    pad1.SetTopMargin(0.0)
    pad1.SetBorderMode(0)
    pad1.Draw()
    pad2.Draw()
    legend = ROOT.TLegend(0.2,0.1,0.7,0.9)
    legend.SetTextFont(42)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.15)
    legend.SetFillColor(0)
    legend.SetFillStyle(0)
    legend.SetLineColor(0)

    pad1.cd()
    Hist_Gamma_Attenu[0].GetXaxis().SetTitle('Emmission Height of Control Region')
    Hist_Gamma_Attenu[0].GetYaxis().SetTitle('#gamma attenuation')
    maxContent = 0.
    for energy in range(0,len(energy_bins)-1):
        maxContent = max(Hist_Gamma_Attenu[energy].GetMaximum(),maxContent)
    Hist_Gamma_Attenu[0].SetMaximum(1.2*maxContent)
    Hist_Gamma_Attenu[0].SetMinimum(0)
    Hist_Gamma_Attenu[0].Draw("E")
    for energy in range(0,len(energy_bins)-1):
        Hist_Gamma_Attenu[energy].SetLineColor(energy+2)
        Hist_Gamma_Attenu[energy].Draw("E same")
    pad2.cd()
    legend.Clear()
    for energy in range(0,len(energy_bins)-1):
        legend.AddEntry(Hist_Gamma_Attenu[energy],'E %s-%s TeV'%(energy_bins[energy],energy_bins[energy+1]),"pl")
    legend.Draw("SAME")
    canvas2.SaveAs('output/Gamma_Attenu_%s.pdf'%(tag))

    pad1.cd()
    Hist_CR_Attenu[0].GetXaxis().SetTitle('Emmission Height of Control Region')
    Hist_CR_Attenu[0].GetYaxis().SetTitle('Bkg attenuation')
    maxContent = 0.
    for energy in range(0,len(energy_bins)-1):
        maxContent = max(Hist_CR_Attenu[energy].GetMaximum(),maxContent)
    Hist_CR_Attenu[0].SetMaximum(1.2*maxContent)
    Hist_CR_Attenu[0].SetMinimum(0)
    Hist_CR_Attenu[0].Draw("E")
    for energy in range(0,len(energy_bins)-1):
        Hist_CR_Attenu[energy].SetLineColor(energy+2)
        Hist_CR_Attenu[energy].Draw("E same")
    pad2.cd()
    legend.Clear()
    for energy in range(0,len(energy_bins)-1):
        legend.AddEntry(Hist_CR_Attenu[energy],'E %s-%s TeV'%(energy_bins[energy],energy_bins[energy+1]),"pl")
    legend.Draw("SAME")
    canvas2.SaveAs('output/CR_Attenu_%s.pdf'%(tag))

    pad1.cd()
    Hist_Inflation[0].GetXaxis().SetTitle('Emmission Height of Control Region')
    Hist_Inflation[0].GetYaxis().SetTitle('Unc. inflation factor')
    maxContent = 0.
    for energy in range(0,len(energy_bins)-1):
        maxContent = max(Hist_Inflation[energy].GetMaximum(),maxContent)
    Hist_Inflation[0].SetMaximum(1.2*maxContent)
    Hist_Inflation[0].SetMinimum(0)
    Hist_Inflation[0].Draw("E")
    for energy in range(0,len(energy_bins)-1):
        Hist_Inflation[energy].SetLineColor(energy+2)
        Hist_Inflation[energy].Draw("E same")
    pad2.cd()
    legend.Clear()
    for energy in range(0,len(energy_bins)-1):
        legend.AddEntry(Hist_Inflation[energy],'E %s-%s TeV'%(energy_bins[energy],energy_bins[energy+1]),"pl")
    legend.Draw("SAME")
    canvas2.SaveAs('output/Inflation_%s.pdf'%(tag))

def set_histStyle( hist , color):
    hist.SetFillColor(color)
    #hist.SetLineColor(1)
    hist.SetLineColor(color)
    hist.SetMarkerColor(color)
    hist.SetFillStyle(1001)
    pass

def MakeStackPlot(Hist_Data,Hist_CR,Hist_Ring,title):
    
    c_both = ROOT.TCanvas("c_both","c both", 200, 10, 600, 600)
    pad3 = ROOT.TPad("pad3","pad3",0,0.8,1,1)
    pad3.SetBottomMargin(0.0)
    pad3.SetTopMargin(0.03)
    pad3.SetBorderMode(1)
    pad1 = ROOT.TPad("pad1","pad1",0,0.3,1,0.8)
    pad1.SetBottomMargin(0.0)
    pad1.SetTopMargin(0.0)
    pad1.SetBorderMode(0)
    pad2 = ROOT.TPad("pad2","pad2",0,0,1,0.3)
    pad2.SetBottomMargin(0.39)
    pad2.SetTopMargin(0.0)
    pad2.SetBorderMode(0)
    pad2.Draw()
    pad1.Draw()
    pad3.Draw()

    pad1.cd()
    Hist_Data.SetLineColor(1)
    Hist_Data.SetLineWidth(3)
    set_histStyle( Hist_CR , 30)
    stack = ROOT.THStack("stack", "")
    stack.Add( Hist_CR )
    if Hist_Data.GetMaximum()>stack.GetMaximum():
        if Hist_Data.GetMinimum()>stack.GetMinimum():
            Hist_Data.SetMinimum(max(1.,stack.GetMinimum()/2.))
        Hist_Data.Draw("E")
    else:
        if Hist_Data.GetMinimum()<stack.GetMinimum():
            stack.SetMinimum(max(1.,Hist_Data.GetMinimum()/2.))
        stack.Draw("hist")
    stack.Draw("hist same")
    Hist_Data.Draw("E same")
    if not Hist_Ring.Integral()==0:
        if add_signal:
            Hist_Ring.SetLineColor(2)
            Hist_Ring.SetLineWidth(2)
            Hist_Ring.Draw("E same")
    Hist_All = Hist_CR.Clone()
    Hist_All.SetFillColor(1)
    Hist_All.SetFillStyle(3004)
    Hist_All.SetMarkerSize(0)
    Hist_All.Draw("e2 same")
    pad3.cd()
    legend = ROOT.TLegend(0.55,0.1,0.94,0.9)
    legend.SetTextFont(42)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.2)
    legend.SetFillColor(0)
    legend.SetFillStyle(0)
    legend.SetLineColor(0)
    legend.Clear()
    legend.AddEntry(Hist_Data,'Data',"pl")
    legend.AddEntry(Hist_CR,'Background',"f")
    if not Hist_Ring.Integral()==0:
        if add_signal:
            legend.AddEntry(Hist_Ring,'True Bkg',"pl")
        #else:
        #    legend.AddEntry(Hist_Ring,'Ring model',"pl")
    legend.Draw("SAME")
    lumilab1 = ROOT.TLatex(0.15,0.75,'   ' )
    lumilab1.SetNDC()
    lumilab1.SetTextSize(0.15)
    lumilab1.Draw()
    lumilab2 = ROOT.TLatex(0.15,0.6,'%0.2f < Tel elev. < %0.2f'%(Elev_cut_lower,Elev_cut_upper) )
    lumilab2.SetNDC()
    lumilab2.SetTextSize(0.15)
    lumilab2.Draw()
    lumilab3 = ROOT.TLatex(0.15,0.45,'%0.2f < MSCL < %0.2f'%(MSCL_cut_lower,MSCL_cut_upper) )
    lumilab3.SetNDC()
    lumilab3.SetTextSize(0.15)
    lumilab3.Draw()
    lumilab4 = ROOT.TLatex(0.15,0.3,'%0.2f < MSCW < %0.2f'%(MSCL_cut_lower,MSCL_cut_upper) )
    lumilab4.SetNDC()
    lumilab4.SetTextSize(0.15)
    lumilab4.Draw()
    pad2.cd()
    Hist_Band = Hist_All.Clone()
    #Hist_Band.Divide(Hist_All)
    for b in range(0,Hist_Band.GetNbinsX()):
        Hist_Band.SetBinContent(b+1,0)
        Hist_Band.SetBinError(b+1,3)
    Hist_Band.SetFillColor(1)
    Hist_Band.SetFillStyle(3004)
    Hist_Band.SetMarkerSize(0)
    Hist_Band.GetXaxis().SetTitle(title)
    Hist_Band.GetXaxis().SetTitleOffset(1.1)
    Hist_Band.GetXaxis().SetTitleSize(0.13)
    Hist_Band.GetXaxis().SetLabelSize(0.1)
    Hist_Band.GetYaxis().SetLabelSize(0.1)
    Hist_Band.GetYaxis().SetTitleOffset(0.3)
    Hist_Band.GetYaxis().SetTitle("#sigma(data vs. bg)")
    #if not Hist_Ring.Integral()==0:
    #    if add_signal:
    #        #Hist_Band.GetYaxis().SetTitle("True/Predicted")
    #        Hist_Band.GetYaxis().SetTitle("#sigma(true bg vs. bg)")
    #    else:
    #        #Hist_Band.GetYaxis().SetTitle("Ring/Bkg")
    #        Hist_Band.GetYaxis().SetTitle("#sigma(ring bg vs. bg)")
    #else:
    #    #Hist_Band.GetYaxis().SetTitle("Data/Bkg")
    #    Hist_Band.GetYaxis().SetTitle("#sigma(data vs. bg)")
    Hist_Band.GetYaxis().SetTitleSize(0.13)
    Hist_Band.GetYaxis().SetNdivisions(505)
    Hist_Band.SetMaximum(5)
    Hist_Band.SetMinimum(-5)
    Hist_Band.Draw("e2")
    Hist_Ratio = Hist_Data.Clone()
    #if not Hist_Ring.Integral()==0:
    #    Hist_Ratio = Hist_Ring.Clone()
    #else:
    #    Hist_Ratio = Hist_Data.Clone()
    #Hist_Ratio.Divide(Hist_All)
    Hist_Ratio.Add(Hist_All,-1.)
    for b in range(0,Hist_Ratio.GetNbinsX()):
        unc = Hist_All.GetBinError(b+1)*Hist_All.GetBinError(b+1) + Hist_Data.GetBinError(b+1)*Hist_Data.GetBinError(b+1)
        #if not Hist_Ring.Integral()==0:
        #    unc = Hist_All.GetBinError(b+1)*Hist_All.GetBinError(b+1) + Hist_Ring.GetBinError(b+1)*Hist_Ring.GetBinError(b+1)
        #else:
        #    unc = Hist_All.GetBinError(b+1)*Hist_All.GetBinError(b+1) + Hist_Data.GetBinError(b+1)*Hist_Data.GetBinError(b+1)
        unc = pow(unc,0.5)
        if not unc==0:
            Hist_Ratio.SetBinContent(b+1,Hist_Ratio.GetBinContent(b+1)/unc)
        else:
            Hist_Ratio.SetBinContent(b+1,0)
        Hist_Ratio.SetBinError(b+1,0)
    Hist_Ratio.SetFillColor(1)
    Hist_Ratio.Draw("B same")
    if title=='Erec':
        pad1.SetLogy()
        pad1.SetLogx()
        pad2.SetLogx()
    c_both.SaveAs('output/%s_%s.pdf'%(title,tag))
    pad1.SetLogy(0)
    pad1.SetLogx(0)
    pad2.SetLogx(0)

def ControlRegionMSCWSelection(tree,elev):
    
    if (tree.MSCW<control_width): return False
    return True

def SignalRegionMSCWSelection(tree,elev):
    
    if (tree.MSCW<MSCW_cut_lower): return False
    if (tree.MSCW>MSCW_cut_upper): return False
    return True

def GetEfficiencyHistogram(source,region,run,data_type):

    if data_type=='anasum':
        if source == 'Crab':
            anasum_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/Crab_V6_Moderate-TMVA-BDT.RB.%s.root"%(run))
        if source == '2ndCrab':
            anasum_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/Crab_V6_Moderate-TMVA-BDT.RB.%s.root"%(run))
        if source == 'PKS1424':
            anasum_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/PKS1424_V6_Moderate-TMVA-BDT.RB.%s.root"%(run))
        if source == '3C264':
            anasum_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/3C264_V6_Moderate-TMVA-BDT.RB.%s.root"%(run))
        if source == 'H1426':
            anasum_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/H1426_V6_Moderate-TMVA-BDT.RB.%s.root"%(run))
        if source == 'Ton599':
            anasum_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/Ton599_V6_Moderate-TMVA-BDT.RB.%s.root"%(run))
        if source == 'IC443':
            anasum_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/IC443_V6_Moderate-TMVA-BDT.RB.%s.root"%(run))
        if region == 'on':
            anasum_tree = anasum_file.Get("run_%s/stereo/data_on"%(run))
            ring_tree = anasum_file.Get("run_%s/stereo/data_off"%(run))
        else:
            anasum_tree = anasum_file.Get("run_%s/stereo/data_off"%(run))
        pointing_tree = anasum_file.Get("run_%s/stereo/pointingDataReduced"%(run))
    else:
        anasum_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/analysis/Results/%s.mscw.root"%(run))
        anasum_tree = anasum_file.Get("data")
        pointing_tree = anasum_file.Get("data")

    pointing_tree.GetEntry(0)
    time_begin = pointing_tree.Time
    pointing_tree.GetEntry(pointing_tree.GetEntries()-1)
    time_end = pointing_tree.Time
    hist_Ele_vs_time = ROOT.TProfile("hist_Ele_vs_time","",100,time_begin,time_end+0.01*(time_end-time_begin),0,100)
    hist_Azi_vs_time = ROOT.TProfile("hist_Azi_vs_time","",100,time_begin,time_end+0.01*(time_end-time_begin),-400,400)
    hist_RA_vs_time = ROOT.TProfile("hist_RA_vs_time","",100,time_begin,time_end+0.01*(time_end-time_begin),0,100)
    hist_Dec_vs_time = ROOT.TProfile("hist_Dec_vs_time","",100,time_begin,time_end+0.01*(time_end-time_begin),0,100)
    for i in range(0,pointing_tree.GetEntries()):
        pointing_tree.GetEntry(i)
        hist_Ele_vs_time.Fill(pointing_tree.Time,pointing_tree.TelElevation)
        hist_Azi_vs_time.Fill(pointing_tree.Time,AzimuthConverter(pointing_tree.TelAzimuth))
        hist_RA_vs_time.Fill(pointing_tree.Time,pointing_tree.TelRAJ2000)
        hist_Dec_vs_time.Fill(pointing_tree.Time,pointing_tree.TelDecJ2000)

    for i in range(0,anasum_tree.GetEntries()):
        anasum_tree.GetEntry(i)
        if not CommonSelection(anasum_tree): continue
        elevation = hist_Ele_vs_time.GetBinContent(hist_Ele_vs_time.FindBin(anasum_tree.Time))
        azimuth = hist_Azi_vs_time.GetBinContent(hist_Azi_vs_time.FindBin(anasum_tree.Time))
        RA = hist_RA_vs_time.GetBinContent(hist_RA_vs_time.FindBin(anasum_tree.Time))
        Dec = hist_Dec_vs_time.GetBinContent(hist_Dec_vs_time.FindBin(anasum_tree.Time))
        if elevation<Elev_cut_lower or elevation>Elev_cut_upper: continue
        elev_scale = Hist_TelElev_Counts_scale.GetBinContent(Hist_TelElev_Counts_scale.FindBin(elevation))
        azim_scale = Hist_TelAzim_Counts_scale.GetBinContent(Hist_TelAzim_Counts_scale.FindBin(AzimuthConverter(azimuth)))
        scale = elev_scale
        if RW_method=='NoRW': 
            scale = 1
        if region == 'off':
            Hist_MSCW_Source.Fill(anasum_tree.MSCW,scale)
            Hist_TelElev_Counts_source_Sum.Fill(elevation,scale)
            Hist_TelAzim_Counts_source_Sum.Fill(AzimuthConverter(azimuth),scale)
        if ControlRegionMSCWSelection(anasum_tree,elevation):
            Hist_OffData_Control.Fill(0,scale)
        if SignalRegionMSCWSelection(anasum_tree,elevation):
            Hist_OffData_Signal.Fill(0,scale)

def MeasureEfficiency():

    Rate_CR_Efficiency = 0
    Err_Rate_CR_Efficiency = 0
    
    if source == 'Crab':
        for run in Crab_runlist:
            GetSourceElevation(run)
    if target == 'Crab':
        for run in Crab_runlist:
            GetTargetElevation(run)
    if target == '2ndCrab':
        for run in SecondCrab_runlist:
            GetTargetElevation(run)
    if target == 'PKS1424':
        for run in PKS1424_runlist:
            GetTargetElevation(run)
    if target == '3C264':
        for run in runlist_3C264:
            GetTargetElevation(run)
    if target == 'H1426':
        for run in runlist_H1426:
            GetTargetElevation(run)
    if target == 'Ton599':
        for run in runlist_Ton599:
            GetTargetElevation(run)
    if target == 'IC443':
        for run in runlist_IC443:
            GetTargetElevation(run)
    #print Hist_TelElev_Counts_target.Print('All')
    for b in range(0,Hist_TelElev_Counts_scale.GetNbinsX()):
        Hist_TelElev_Counts_scale.SetBinContent(b+1,Hist_TelElev_Counts_target.GetBinContent(b+1))
    Hist_TelElev_Counts_scale.Divide(Hist_TelElev_Counts_source)
    for b in range(0,Hist_TelAzim_Counts_scale.GetNbinsX()):
        Hist_TelAzim_Counts_scale.SetBinContent(b+1,Hist_TelAzim_Counts_target.GetBinContent(b+1))
    Hist_TelAzim_Counts_scale.Divide(Hist_TelAzim_Counts_source)
    #print Hist_TelElev_Counts_scale.Print('All')
    if source == 'Crab':
        for run in Crab_runlist:
            GetEfficiencyHistogram('Crab','off',run,'anasum')
    if source == '2ndCrab':
        for run in SecondCrab_runlist:
            GetEfficiencyHistogram('Crab','off',run,'anasum')
    Hist_TelElev_Counts_scale.Reset()
    Hist_TelElev_Counts_source.Reset()
    Hist_TelElev_Counts_target.Reset()
    Hist_TelAzim_Counts_scale.Reset()
    Hist_TelAzim_Counts_source.Reset()
    Hist_TelAzim_Counts_target.Reset()
    
    print 'Hist_OffData_Signal.GetEntries() = %s'%(Hist_OffData_Signal.GetEntries())
    print 'Hist_OffData_Control.GetEntries() = %s'%(Hist_OffData_Control.GetEntries())
    print 'Hist_OffData_Signal.Integral() = %s'%(Hist_OffData_Signal.GetBinContent(Hist_OffData_Signal.FindBin(0)))

    for b in range(0,Hist_CR_Efficiency.GetNbinsX()):
        Hist_CR_Efficiency.SetBinContent(b+1,Hist_OffData_Control.GetBinContent(b+1))
        Hist_CR_Efficiency.SetBinError(b+1,Hist_OffData_Control.GetBinError(b+1))
    Hist_CR_Efficiency.Divide(Hist_OffData_Signal)
    Rate_CR_Efficiency = Hist_CR_Efficiency.GetBinContent(Hist_CR_Efficiency.FindBin(0))
    Err_Rate_CR_Efficiency = Hist_CR_Efficiency.GetBinError(Hist_CR_Efficiency.FindBin(0))
    print 'Rate_CR_Efficiency = %s +/- %s'%(Rate_CR_Efficiency,Err_Rate_CR_Efficiency)
    
    Hist_OffData_Signal.Reset()
    Hist_OffData_Control.Reset()
    Hist_CR_Efficiency.Reset()

    return Rate_CR_Efficiency,Err_Rate_CR_Efficiency

def RunExtendedSourceAnalysisMSCWMethod(data_type,CR_Efficiency,err_CR_Efficiency):

    print 'Erec_cut_lower = %s, Erec_cut_upper = %s'%(Erec_cut_lower,Erec_cut_upper)
    print 'CR_Efficiency = %s +/- %s'%(CR_Efficiency,err_CR_Efficiency)

    if target == 'Crab':
        for run in Crab_runlist:
            SingleRunAnalysisMSCW('Crab',target_field,run,data_type)
    if target == '2ndCrab':
        for run in SecondCrab_runlist:
            SingleRunAnalysisMSCW('Crab',target_field,run,data_type)
    if target == 'PKS1424':
        for run in PKS1424_runlist:
            SingleRunAnalysisMSCW('PKS1424',target_field,run,data_type)
    if target == '3C264':
        for run in runlist_3C264:
            SingleRunAnalysisMSCW('3C264',target_field,run,data_type)
    if target == 'H1426':
        for run in runlist_H1426:
            SingleRunAnalysisMSCW('H1426',target_field,run,data_type)
    if target == 'Ton599':
        for run in runlist_Ton599:
            SingleRunAnalysisMSCW('Ton599',target_field,run,data_type)
    if target == 'IC443':
        for run in runlist_IC443:
            SingleRunAnalysisMSCW('IC443',target_field,run,data_type)
    
    # Here we add signals using Crab data
    if add_signal:
        for run in SecondCrab_runlist:
            #GenerateExtendedSignalMSCW('Crab','on',run,data_type,1.0)
            GenerateExtendedSignalMSCW('Crab','off',run,data_type,1.0)
            #GenerateExtendedSignalMSCW('Crab','off',run,data_type,scale*Normalization_Factor_Signal)

    if data_type=='anasum':
        Normalization_Factor_Ring = 0
        if not Hist_Norm_Ring.Integral()==0:
            Normalization_Factor_Ring = Hist_Norm_Data.Integral()
            Normalization_Factor_Ring = Normalization_Factor_Ring/Hist_Norm_Ring.Integral()
        Hist_Erec_Ring.Scale(Normalization_Factor_Ring)

    #N_signal (number of events in signal region) = N_egamma + N_CR
    #N_control (number of events in control region) = 0*N_egamma + b*N_CR, b = CR efficiency
    #N_CR = N_control/b

    print 'Hist_Erec_CR_Raw.Integral() = %s'%(Hist_Erec_CR_Raw.Integral())
    print 'Hist_Erec_Data.Integral() = %s'%(Hist_Erec_Data.Integral())
    print 'Hist_Erec_CR_Raw_AddSignal.Integral() = %s'%(Hist_Erec_CR_Raw_AddSignal.Integral())
    print 'Hist_Erec_Data_AddSignal.Integral() = %s'%(Hist_Erec_Data_AddSignal.Integral())
    if not Fake_Removal=='NoFakeRemoved':
        if not add_signal:
            CalculateBkgMSCWMethod(Hist_Erec_CR,Hist_Erec_Data,Hist_Erec_CR_Raw,CR_Efficiency,err_CR_Efficiency)
        else:
            Hist_Erec_Data_AddSignal.Add(Hist_Erec_Data)
            Hist_Erec_CR_Raw_AddSignal.Add(Hist_Erec_CR_Raw)
            CalculateBkgMSCWMethod(Hist_Erec_CR,Hist_Erec_Data_AddSignal,Hist_Erec_CR_Raw_AddSignal,CR_Efficiency,err_CR_Efficiency)
        if not add_signal:
            CalculateBkgMSCWMethod(Hist2D_Xoff_vs_Yoff_CR,Hist2D_Xoff_vs_Yoff_Data,Hist2D_Xoff_vs_Yoff_CR_Raw,CR_Efficiency,err_CR_Efficiency)
        else:
            Hist2D_Xoff_vs_Yoff_Data_AddSignal.Add(Hist2D_Xoff_vs_Yoff_Data)
            Hist2D_Xoff_vs_Yoff_CR_Raw_AddSignal.Add(Hist2D_Xoff_vs_Yoff_CR_Raw)
            CalculateBkgMSCWMethod(Hist2D_Xoff_vs_Yoff_CR,Hist2D_Xoff_vs_Yoff_Data_AddSignal,Hist2D_Xoff_vs_Yoff_CR_Raw_AddSignal,CR_Efficiency,err_CR_Efficiency)
        if not add_signal:
            CalculateBkgMSCWMethod(Hist_theta2_CR,Hist_theta2_Data,Hist_theta2_CR_Raw,CR_Efficiency,err_CR_Efficiency)
        else:
            Hist_theta2_Data_AddSignal.Add(Hist_theta2_Data)
            Hist_theta2_CR_Raw_AddSignal.Add(Hist_theta2_CR_Raw)
            CalculateBkgMSCWMethod(Hist_theta2_CR,Hist_theta2_Data_AddSignal,Hist_theta2_CR_Raw_AddSignal,CR_Efficiency,err_CR_Efficiency)
        if not add_signal:
            CalculateBkgMSCWMethod(Hist_theta2_zoomin_CR,Hist_theta2_zoomin_Data,Hist_theta2_zoomin_CR_Raw,CR_Efficiency,err_CR_Efficiency)
        else:
            Hist_theta2_zoomin_Data_AddSignal.Add(Hist_theta2_zoomin_Data)
            Hist_theta2_zoomin_CR_Raw_AddSignal.Add(Hist_theta2_zoomin_CR_Raw)
            CalculateBkgMSCWMethod(Hist_theta2_zoomin_CR,Hist_theta2_zoomin_Data_AddSignal,Hist_theta2_zoomin_CR_Raw_AddSignal,CR_Efficiency,err_CR_Efficiency)
    Err_N_CR = ROOT.Double(0.)
    N_CR = Hist_Erec_CR.IntegralAndError(1,Hist_Erec_CR.GetNbinsX(),Err_N_CR)
    print 'Hist_Erec_CR.Integral() = %s +/- %s'%(N_CR,Err_N_CR)
    print 'Hist_Erec_Ring.Integral() = %s'%(Hist_Erec_Ring.Integral())
    #print Hist_theta2_CR_Raw.Print("All")
    #print Hist_theta2_CR.Print("All")
    
    
    Hist_Erec_Data_Sum.Add(Hist_Erec_Data)
    if add_signal: 
        Hist_Erec_Data_Sum_AddSignal.Add(Hist_Erec_Data)
        Hist_Erec_Data_Sum_AddSignal.Add(Hist_Erec_Data_AddSignal)
    Hist_Erec_Ring_Sum.Add(Hist_Erec_Ring)
    Hist_Erec_CR_Sum.Add(Hist_Erec_CR)

    Hist2D_Xoff_vs_Yoff_Data_Sum.Add(Hist2D_Xoff_vs_Yoff_Data)
    if add_signal: 
        Hist2D_Xoff_vs_Yoff_Data_Sum_AddSignal.Add(Hist2D_Xoff_vs_Yoff_Data)
        Hist2D_Xoff_vs_Yoff_Data_Sum_AddSignal.Add(Hist2D_Xoff_vs_Yoff_Data_AddSignal)
    Hist2D_Xoff_vs_Yoff_CR_Sum.Add(Hist2D_Xoff_vs_Yoff_CR)

    Hist_theta2_Data_Sum.Add(Hist_theta2_Data)
    if add_signal: 
        Hist_theta2_Data_Sum_AddSignal.Add(Hist_theta2_Data)
        Hist_theta2_Data_Sum_AddSignal.Add(Hist_theta2_Data_AddSignal)
    Hist_theta2_CR_Sum.Add(Hist_theta2_CR)

    Hist_theta2_zoomin_Data_Sum.Add(Hist_theta2_zoomin_Data)
    if add_signal: 
        Hist_theta2_zoomin_Data_Sum_AddSignal.Add(Hist_theta2_zoomin_Data)
        Hist_theta2_zoomin_Data_Sum_AddSignal.Add(Hist_theta2_zoomin_Data_AddSignal)
    Hist_theta2_zoomin_CR_Sum.Add(Hist_theta2_zoomin_CR)

    MakeStackPlot(Hist_theta2_Data,Hist_theta2_CR,Hist_theta2_Ring_Sum,'theta2')
    #MakeStackPlot(Hist_theta2_zoomin_Data,Hist_theta2_zoomin_CR,Hist_theta2_zoomin_Ring_Sum,'theta2_zoomin')

    Hist_Norm_Data.Reset()
    Hist_Norm_Ring.Reset()
    Hist_Erec_CR.Reset()
    Hist_Erec_CR_Raw.Reset()
    Hist_Erec_CR_Raw_AddSignal.Reset()
    Hist_Erec_Data.Reset()
    Hist_Erec_Data_AddSignal.Reset()
    Hist_Erec_Ring.Reset()
    Hist2D_Xoff_vs_Yoff_Data.Reset()
    Hist2D_Xoff_vs_Yoff_Data_AddSignal.Reset()
    Hist2D_Xoff_vs_Yoff_CR.Reset()
    Hist2D_Xoff_vs_Yoff_CR_Raw.Reset()
    Hist2D_Xoff_vs_Yoff_CR_Raw_AddSignal.Reset()
    Hist_theta2_Data.Reset()
    Hist_theta2_Data_AddSignal.Reset()
    Hist_theta2_CR.Reset()
    Hist_theta2_CR_Raw.Reset()
    Hist_theta2_CR_Raw_AddSignal.Reset()
    Hist_theta2_zoomin_Data.Reset()
    Hist_theta2_zoomin_Data_AddSignal.Reset()
    Hist_theta2_zoomin_CR.Reset()
    Hist_theta2_zoomin_CR_Raw.Reset()
    Hist_theta2_zoomin_CR_Raw_AddSignal.Reset()

def EfficiencyRateAtDifferentMSCW():

    global Erec_cut_lower
    global Erec_cut_upper
    global MSCL_cut_lower
    global MSCL_cut_upper
    global control_width
    global Hist_CR_Efficiency_list
    global Hist_Inflation

    for energy in range(0,len(energy_bins)-1):
        Hist_CR_Efficiency_list += [Hist_ControlWidth_Bins.Clone()]
        Hist_Inflation += [Hist_ControlWidth_Bins.Clone()]
        for width in range(1,Hist_ControlWidth_Bins.GetNbinsX()+1):
            control_width = Hist_ControlWidth_Bins.GetBinLowEdge(width)
            Erec_cut_lower = energy_bins[energy]
            Erec_cut_upper = energy_bins[energy+1]
            #if energy==len(energy_bins)-2: Erec_cut_upper = 1e10
            MakeATag()
            print '++++++++++++++++++++++++++++++++++++++++++++++'
            print 'selection = %s, control region = %s'%(tag, control_width)
            CR_Efficiency,err_CR_Efficiency = MeasureEfficiency()
            Hist_CR_Efficiency_list[energy].SetBinContent(Hist_ControlWidth_Bins.FindBin(control_width),CR_Efficiency)
            Hist_CR_Efficiency_list[energy].SetBinError(Hist_ControlWidth_Bins.FindBin(control_width),err_CR_Efficiency)
            inflation = 1./(CR_Efficiency)
            print 'inflation = %s'%(inflation)
            Hist_Inflation[energy].SetBinContent(Hist_ControlWidth_Bins.FindBin(control_width),inflation)

    canvas2 = ROOT.TCanvas("canvas2","canvas2", 200, 10, 600, 600)
    pad2 = ROOT.TPad("pad2","pad2",0.1,0.8,1,1)
    pad2.SetBottomMargin(0.0)
    pad2.SetTopMargin(0.03)
    pad2.SetBorderMode(1)
    pad1 = ROOT.TPad("pad1","pad1",0.1,0,1,0.8)
    pad1.SetBottomMargin(0.2)
    pad1.SetTopMargin(0.0)
    pad1.SetBorderMode(0)
    pad1.Draw()
    pad2.Draw()
    legend = ROOT.TLegend(0.2,0.1,0.7,0.9)
    legend.SetTextFont(42)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.15)
    legend.SetFillColor(0)
    legend.SetFillStyle(0)
    legend.SetLineColor(0)

    pad1.cd()
    Hist_CR_Efficiency_list[0].GetXaxis().SetTitle('Emmission Height of Control Region')
    Hist_CR_Efficiency_list[0].GetYaxis().SetTitle('Bkg Efficiency')
    maxContent = 0.
    for energy in range(0,len(energy_bins)-1):
        maxContent = max(Hist_CR_Efficiency_list[energy].GetMaximum(),maxContent)
    Hist_CR_Efficiency_list[0].SetMaximum(1.2*maxContent)
    Hist_CR_Efficiency_list[0].SetMinimum(0)
    Hist_CR_Efficiency_list[0].Draw("E")
    for energy in range(0,len(energy_bins)-1):
        Hist_CR_Efficiency_list[energy].SetLineColor(energy+2)
        Hist_CR_Efficiency_list[energy].Draw("E same")
    pad2.cd()
    legend.Clear()
    for energy in range(0,len(energy_bins)-1):
        legend.AddEntry(Hist_CR_Efficiency_list[energy],'E %s-%s TeV'%(energy_bins[energy],energy_bins[energy+1]),"pl")
    legend.Draw("SAME")
    canvas2.SaveAs('output/CR_Efficiency_%s.pdf'%(tag))

    pad1.cd()
    Hist_Inflation[0].GetXaxis().SetTitle('Emmission Height of Control Region')
    Hist_Inflation[0].GetYaxis().SetTitle('Unc. inflation factor')
    maxContent = 0.
    for energy in range(0,len(energy_bins)-1):
        maxContent = max(Hist_Inflation[energy].GetMaximum(),maxContent)
    Hist_Inflation[0].SetMaximum(1.2*maxContent)
    Hist_Inflation[0].SetMinimum(0)
    Hist_Inflation[0].Draw("E")
    for energy in range(0,len(energy_bins)-1):
        Hist_Inflation[energy].SetLineColor(energy+2)
        Hist_Inflation[energy].Draw("E same")
    pad2.cd()
    legend.Clear()
    for energy in range(0,len(energy_bins)-1):
        legend.AddEntry(Hist_Inflation[energy],'E %s-%s TeV'%(energy_bins[energy],energy_bins[energy+1]),"pl")
    legend.Draw("SAME")
    canvas2.SaveAs('output/Inflation_%s.pdf'%(tag))

def EnergySpectrum(method):

    global Erec_cut_lower
    global Erec_cut_upper
    global control_width
    global observation_time

    for energy in range(0,len(energy_bins)-1):
        observation_time = 0
        Erec_cut_lower = energy_bins[energy]
        Erec_cut_upper = energy_bins[energy+1]
        #if energy==len(energy_bins)-2: Erec_cut_upper = 1000
        MakeATag()
        cr_eff = 1
        control_width = ControlWidthAtThisEnergy[energy]
        print '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
        print 'control region = %s'%(control_width)
        if method == "Height":
            cr_att = 1
            gamma_att = 0
            control_width = ControlWidthAtThisEnergy[energy]
            cr_att = Hist_CR_Attenu[energy].GetBinContent(Hist_CR_Attenu[energy].FindBin(control_width))
            err_cr_att = Hist_CR_Attenu[energy].GetBinError(Hist_CR_Attenu[energy].FindBin(control_width))
            gamma_att = Hist_Gamma_Attenu[energy].GetBinContent(Hist_Gamma_Attenu[energy].FindBin(control_width))
            err_gamma_att = Hist_Gamma_Attenu[energy].GetBinError(Hist_Gamma_Attenu[energy].FindBin(control_width))
            inflation = pow(1+gamma_att/cr_att,0.5)/(1-gamma_att/cr_att)
            inflation = inflation/pow(cr_att,0.5)
            print 'inflation = %s'%(inflation)
            RunExtendedSourceAnalysis('anasum',cr_att,err_cr_att,gamma_att,err_gamma_att)
        if method == "MSCW":
            cr_eff = Hist_CR_Efficiency_list[energy].GetBinContent(Hist_CR_Efficiency_list[energy].FindBin(control_width))
            err_cr_eff = Hist_CR_Efficiency_list[energy].GetBinError(Hist_CR_Efficiency_list[energy].FindBin(control_width))
            inflation = 1./cr_eff
            print 'inflation = %s'%(inflation)
            RunExtendedSourceAnalysisMSCWMethod('anasum',cr_eff,err_cr_eff)
        N_CR = Hist_Erec_CR_Sum.GetBinContent(energy+1)
        Err_N_CR = Hist_Erec_CR_Sum.GetBinError(energy+1)
        N_Data = Hist_Erec_Data_Sum.GetBinContent(energy+1)
        Err_N_Data = Hist_Erec_Data_Sum.GetBinError(energy+1)
        Hist_Sensitivity_vs_Erec.SetBinContent(energy+1,5.*Err_N_CR/N_CR)
        Hist_Sensitivity_vs_Erec_Other.SetBinContent(energy+1,5.*pow(2,0.5)*Err_N_Data/N_Data)
        canvas = ROOT.TCanvas("canvas","canvas", 200, 10, 600, 600)
        pad1 = ROOT.TPad("pad1","pad1",0,0,1,1)
        pad1.SetBottomMargin(0.15)
        pad1.SetRightMargin(0.15)
        pad1.SetLeftMargin(0.15)
        pad1.SetTopMargin(0.15)
        pad1.SetBorderMode(0)
        pad1.Draw()
        pad1.cd()
        pad1.SetLogz(0)
        pad1.SetLogx(0)
        pad1.SetLogy(0)
        for bx in range(0,Hist2D_Xoff_vs_Yoff_Sig.GetNbinsX()):
            for by in range(0,Hist2D_Xoff_vs_Yoff_Sig.GetNbinsY()):
                this_data = Hist2D_Xoff_vs_Yoff_Data_Sum.GetBinContent(bx,by)
                this_data_err = Hist2D_Xoff_vs_Yoff_Data_Sum.GetBinError(bx,by)
                if add_signal: 
                    this_data = Hist2D_Xoff_vs_Yoff_Data_Sum_AddSignal.GetBinContent(bx,by)
                    this_data_err = Hist2D_Xoff_vs_Yoff_Data_Sum_AddSignal.GetBinError(bx,by)
                this_bkg = Hist2D_Xoff_vs_Yoff_CR_Sum.GetBinContent(bx,by)
                this_bkg_err = Hist2D_Xoff_vs_Yoff_CR_Sum.GetBinError(bx,by)
                value = 0
                if not this_bkg_err==0:
                    value = (this_data-this_bkg)/pow(this_bkg_err*this_bkg_err+this_data_err*this_data_err,0.5)
                Hist2D_Xoff_vs_Yoff_Sig.SetBinContent(bx,by,value)
        Hist2D_Xoff_vs_Yoff_Sig.GetYaxis().SetTitle('Yoff')
        Hist2D_Xoff_vs_Yoff_Sig.GetXaxis().SetTitle('Xoff')
        Hist2D_Xoff_vs_Yoff_Sig.GetZaxis().SetTitle('Significance')
        Hist2D_Xoff_vs_Yoff_Sig.SetMaximum(5)
        Hist2D_Xoff_vs_Yoff_Sig.SetMinimum(-3)
        Hist2D_Xoff_vs_Yoff_Sig.Draw("COL4Z")
        canvas.SaveAs('output/Xoff_vs_Yoff_Sig_%s.pdf'%(tag))
    Erec_cut_lower = energy_bins[0]
    MakeATag()
    MakePlots()

def SensitivityVsTime(method):

    global time_cut
    global observation_time

    for time in range(0,len(time_bins)-1):
        Hist_Erec_CR_Sum.Reset()
        Hist_Erec_Data_Sum.Reset()
        Hist_Erec_Data_Sum_AddSignal.Reset()
        Hist_Erec_Ring_Sum.Reset()
        Hist2D_Xoff_vs_Yoff_Data_Sum.Reset()
        Hist2D_Xoff_vs_Yoff_Data_Sum_AddSignal.Reset()
        Hist2D_Xoff_vs_Yoff_CR_Sum.Reset()
        Hist_theta2_Data_Sum.Reset()
        Hist_theta2_Data_Sum_AddSignal.Reset()
        Hist_theta2_CR_Sum.Reset()
        Hist_theta2_zoomin_Data_Sum.Reset()
        Hist_theta2_zoomin_Data_Sum_AddSignal.Reset()
        Hist_theta2_zoomin_CR_Sum.Reset()
        Hist_EmissionHeight_Source.Reset()
        Hist_MSCW_Source.Reset()
        Hist_EmissionHeight_Target.Reset()
        Hist_MSCW_Target.Reset()
        observation_time = 0
        time_cut = time_bins[time]*60.*60.
        #print 'time cut = %s hrs'%(time_bins[time])
        EnergySpectrum(method)
        Err_N_CR = ROOT.Double(0.)
        N_CR = Hist_Erec_CR_Sum.IntegralAndError(1,Hist_Erec_CR_Sum.GetNbinsX(),Err_N_CR)
        Err_N_Data = ROOT.Double(0.)
        N_Data = Hist_Erec_Data_Sum.IntegralAndError(1,Hist_Erec_Data_Sum.GetNbinsX(),Err_N_Data)
        Hist_NCR_vs_Time.SetBinContent(time_bins[time],N_CR)
        Hist_NCR_vs_Time.SetBinError(time_bins[time],Err_N_CR)
        Hist_NData_vs_Time.SetBinContent(time_bins[time],N_Data)
        Hist_NData_vs_Time.SetBinError(time_bins[time],Err_N_Data)
        Hist_Sensitivity_vs_Time.SetBinContent(time_bins[time],5.*Err_N_CR/N_CR)
        Hist_Sensitivity_vs_Time_Other.SetBinContent(time_bins[time],5.*pow(2,0.5)*Err_N_Data/N_Data)
    MakeStackPlot(Hist_NData_vs_Time,Hist_NCR_vs_Time,Hist_NRing_vs_Time,'Time')
    canvas = ROOT.TCanvas("canvas","canvas", 200, 10, 600, 600)
    Hist_Sensitivity_vs_Time.GetYaxis().SetTitle('Signal flux excluded (Bg unit)')
    Hist_Sensitivity_vs_Time.GetXaxis().SetTitle('Observation hours')
    Hist_Sensitivity_vs_Time.Draw()
    Hist_Sensitivity_vs_Time_Other.SetLineColor(2)
    Hist_Sensitivity_vs_Time_Other.Draw("same")
    canvas.SaveAs('output/Sensitivity_vs_ObsTime_%s.pdf'%(tag))

def MakePlots():
    MakeStackPlot(Hist_Erec_Data_Sum,Hist_Erec_CR_Sum,Hist_Erec_Ring_Sum,'Erec')
    MakeStackPlot(Hist_theta2_Data_Sum,Hist_theta2_CR_Sum,Hist_theta2_Ring_Sum,'theta2')
    MakeStackPlot(Hist_theta2_zoomin_Data_Sum,Hist_theta2_zoomin_CR_Sum,Hist_theta2_zoomin_Ring_Sum,'theta2_zoomin')
    if add_signal: 
        MakeStackPlot(Hist_Erec_Data_Sum_AddSignal,Hist_Erec_CR_Sum,Hist_Erec_Data_Sum,'Erec_addSignal')
        MakeStackPlot(Hist_theta2_Data_Sum_AddSignal,Hist_theta2_CR_Sum,Hist_theta2_Data_Sum,'theta2_addSignal')
        MakeStackPlot(Hist_theta2_zoomin_Data_Sum_AddSignal,Hist_theta2_zoomin_CR_Sum,Hist_theta2_zoomin_Data_Sum,'theta2_zoomin_addSignal')
    #print Hist_theta2_CR_Sum.Print("All")

    canvas = ROOT.TCanvas("canvas","canvas", 200, 10, 600, 600)
    pad1 = ROOT.TPad("pad1","pad1",0,0,1,1)
    pad1.SetBottomMargin(0.15)
    pad1.SetRightMargin(0.15)
    pad1.SetLeftMargin(0.15)
    pad1.SetTopMargin(0.15)
    pad1.SetBorderMode(0)
    pad1.Draw()
    pad1.cd()
    Hist2D_TelElev_vs_EmissionHeight_Data_All.GetYaxis().SetTitle('EmissionHeight')
    Hist2D_TelElev_vs_EmissionHeight_Data_All.GetXaxis().SetTitle('Tel. elevation')
    Hist2D_TelElev_vs_EmissionHeight_Data_All.Draw("COL4Z")
    for b in range(0,Hist1D_TelElev_vs_EmissionHeight_Data_All.GetNbinsX()):
        hist_temp = Hist2D_TelElev_vs_EmissionHeight_Data_All.ProjectionY("hist_temp",b+1,b+1)
        Hist1D_TelElev_vs_EmissionHeight_Data_All.SetBinContent(b+1,hist_temp.GetMean())
        Hist1D_TelElev_vs_EmissionHeight_Data_All.SetBinError(b+1,hist_temp.GetRMS())
    Hist1D_TelElev_vs_EmissionHeight_Data_All.SetLineColor(2)
    Hist1D_TelElev_vs_EmissionHeight_Data_All.Draw("E same")
    #f_cut_mean = ROOT.TF1("f_cut_mean","10.+10.*pow((x-90)/50.,4)",0,90)
    #f_cut_mean.SetLineColor(2)
    #f_cut_mean.SetLineStyle(2)
    #f_cut_mean.Draw("same")
    #f_cut_upper = ROOT.TF1("f_cut_upper","3.0+10.+10.*pow((x-90)/47.,4)",0,90)
    #f_cut_upper.SetLineColor(2)
    #f_cut_upper.Draw("same")
    #f_cut_lower = ROOT.TF1("f_cut_lower","-3.0+10.+10.*pow((x-90)/53.,4)",0,90)
    #f_cut_lower.SetLineColor(2)
    #f_cut_lower.Draw("same")
    pad1.SetLogz()
    canvas.SaveAs('output/TelElev_vs_EmissionHeight_All_%s.pdf'%(tag))

    Hist2D_Erec_vs_EmissionHeight_Data_All.GetYaxis().SetTitle('EmissionHeight')
    Hist2D_Erec_vs_EmissionHeight_Data_All.GetXaxis().SetTitle('Erec')
    Hist2D_Erec_vs_EmissionHeight_Data_All.Draw("COL4Z")
    for b in range(0,Hist1D_Erec_vs_EmissionHeight_Data_All.GetNbinsX()):
        hist_temp = Hist2D_Erec_vs_EmissionHeight_Data_All.ProjectionY("hist_temp",b+1,b+1)
        Hist1D_Erec_vs_EmissionHeight_Data_All.SetBinContent(b+1,hist_temp.GetMean())
        Hist1D_Erec_vs_EmissionHeight_Data_All.SetBinError(b+1,hist_temp.GetRMS())
    Hist1D_Erec_vs_EmissionHeight_Data_All.SetLineColor(2)
    Hist1D_Erec_vs_EmissionHeight_Data_All.Draw("E same")
    pad1.SetLogz()
    pad1.SetLogx()
    canvas.SaveAs('output/Erec_vs_EmissionHeight_All_%s.pdf'%(tag))

    pad1.SetLogz(0)
    pad1.SetLogx(0)
    pad1.SetLogy(0)
    for bx in range(0,Hist2D_Xoff_vs_Yoff_Sig.GetNbinsX()):
        for by in range(0,Hist2D_Xoff_vs_Yoff_Sig.GetNbinsY()):
            this_data = Hist2D_Xoff_vs_Yoff_Data_Sum.GetBinContent(bx,by)
            this_data_err = Hist2D_Xoff_vs_Yoff_Data_Sum.GetBinError(bx,by)
            if add_signal: 
                this_data = Hist2D_Xoff_vs_Yoff_Data_Sum_AddSignal.GetBinContent(bx,by)
                this_data_err = Hist2D_Xoff_vs_Yoff_Data_Sum_AddSignal.GetBinError(bx,by)
            this_bkg = Hist2D_Xoff_vs_Yoff_CR_Sum.GetBinContent(bx,by)
            this_bkg_err = Hist2D_Xoff_vs_Yoff_CR_Sum.GetBinError(bx,by)
            value = 0
            if not this_bkg_err==0:
                value = (this_data-this_bkg)/pow(this_bkg_err*this_bkg_err+this_data_err*this_data_err,0.5)
            Hist2D_Xoff_vs_Yoff_Sig.SetBinContent(bx,by,value)
    Hist2D_Xoff_vs_Yoff_Sig.GetYaxis().SetTitle('Yoff')
    Hist2D_Xoff_vs_Yoff_Sig.GetXaxis().SetTitle('Xoff')
    Hist2D_Xoff_vs_Yoff_Sig.GetZaxis().SetTitle('Significance')
    Hist2D_Xoff_vs_Yoff_Sig.SetMaximum(5)
    Hist2D_Xoff_vs_Yoff_Sig.SetMinimum(-3)
    Hist2D_Xoff_vs_Yoff_Sig.Draw("COL4Z")
    canvas.SaveAs('output/Xoff_vs_Yoff_Sig_%s.pdf'%(tag))
    Hist2D_Xoff_vs_Yoff_Data_Sum.GetYaxis().SetTitle('Yoff')
    Hist2D_Xoff_vs_Yoff_Data_Sum.GetXaxis().SetTitle('Xoff')
    Hist2D_Xoff_vs_Yoff_Data_Sum.Draw("COL4Z")
    canvas.SaveAs('output/Xoff_vs_Yoff_Data_%s.pdf'%(tag))
    Hist2D_Xoff_vs_Yoff_CR_Sum.GetYaxis().SetTitle('Yoff')
    Hist2D_Xoff_vs_Yoff_CR_Sum.GetXaxis().SetTitle('Xoff')
    Hist2D_Xoff_vs_Yoff_CR_Sum.Draw("COL4Z")
    canvas.SaveAs('output/Xoff_vs_Yoff_CR_%s.pdf'%(tag))

    Hist_Sensitivity_vs_Erec.GetYaxis().SetTitle('Signal flux excluded (Bg unit.)')
    Hist_Sensitivity_vs_Erec.GetXaxis().SetTitle('Erec')
    Hist_Sensitivity_vs_Erec.Draw()
    Hist_Sensitivity_vs_Erec_Other.SetLineColor(2)
    Hist_Sensitivity_vs_Erec_Other.Draw("same")
    pad1.SetLogx()
    canvas.SaveAs('output/Sensitivity_vs_Erec_%s.pdf'%(tag))

    pad1.SetLogx(0)
    Hist_TelElev_Counts_source_Sum.GetXaxis().SetTitle('Tel. elevation')
    Hist_TelElev_Counts_source_Sum.SetLineColor(4)
    Hist_TelElev_Counts_target_Sum.SetLineColor(2)
    numerator = Hist_TelElev_Counts_target_Sum.Integral()
    denominator = Hist_TelElev_Counts_source_Sum.Integral()
    scale = numerator/denominator
    Hist_TelElev_Counts_source_Sum.Scale(scale)
    if Hist_TelElev_Counts_target_Sum.GetMaximum()>Hist_TelElev_Counts_source_Sum.GetMaximum():
        Hist_TelElev_Counts_target_Sum.Draw("L")
        Hist_TelElev_Counts_source_Sum.Draw("L same")
    else:
        Hist_TelElev_Counts_source_Sum.Draw("L")
        Hist_TelElev_Counts_target_Sum.Draw("L same")
    canvas.SaveAs('output/TelElev_Counts_%s.pdf'%(tag))

    Hist_TelAzim_Counts_source_Sum.GetXaxis().SetTitle('Tel. azimuth')
    Hist_TelAzim_Counts_source_Sum.SetLineColor(4)
    Hist_TelAzim_Counts_target_Sum.SetLineColor(2)
    numerator = Hist_TelAzim_Counts_target_Sum.Integral()
    denominator = Hist_TelAzim_Counts_source_Sum.Integral()
    scale = numerator/denominator
    Hist_TelAzim_Counts_source_Sum.Scale(scale)
    if Hist_TelAzim_Counts_target_Sum.GetMaximum()>Hist_TelAzim_Counts_source_Sum.GetMaximum():
        Hist_TelAzim_Counts_target_Sum.Draw("L")
        Hist_TelAzim_Counts_source_Sum.Draw("L same")
    else:
        Hist_TelAzim_Counts_source_Sum.Draw("L")
        Hist_TelAzim_Counts_target_Sum.Draw("L same")
    canvas.SaveAs('output/TelAzim_Counts_%s.pdf'%(tag))

    if Hist_EmissionHeight_Source.Integral()!=0:
        numerator = Hist_EmissionHeight_Target.Integral(Hist_EmissionHeight_Target.FindBin(15),Hist_EmissionHeight_Target.FindBin(20))
        denominator = Hist_EmissionHeight_Source.Integral(Hist_EmissionHeight_Target.FindBin(15),Hist_EmissionHeight_Target.FindBin(20))
        scale = numerator/denominator
        Hist_EmissionHeight_Source.Scale(scale)
        Hist_EmissionHeight_Target.SetLineColor(2)
        Hist_EmissionHeight_Source.SetLineColor(4)
        if Hist_EmissionHeight_Target.GetMaximum()>Hist_EmissionHeight_Source.GetMaximum():
            Hist_EmissionHeight_Target.Draw("L")
            Hist_EmissionHeight_Source.Draw("L same")
        else:
            Hist_EmissionHeight_Source.Draw("L")
            Hist_EmissionHeight_Target.Draw("L same")
        pad1.SetLogx(0)
        canvas.SaveAs('output/Source_vs_Target_EmissionHeight_%s.pdf'%(tag))

    if Hist_MSCW_Source.Integral()!=0:
        numerator = Hist_MSCW_Target.Integral(Hist_MSCW_Target.FindBin(3),Hist_MSCW_Target.FindBin(9))
        denominator = Hist_MSCW_Source.Integral(Hist_MSCW_Target.FindBin(3),Hist_MSCW_Target.FindBin(9))
        scale = numerator/denominator
        Hist_MSCW_Source.Scale(scale)
        Hist_MSCW_Target.SetLineColor(2)
        Hist_MSCW_Source.SetLineColor(4)
        if Hist_MSCW_Target.GetMaximum()>Hist_MSCW_Source.GetMaximum():
            Hist_MSCW_Target.Draw("L")
            Hist_MSCW_Source.Draw("L same")
        else:
            Hist_MSCW_Source.Draw("L")
            Hist_MSCW_Target.Draw("L same")
        pad1.SetLogx(0)
        canvas.SaveAs('output/Source_vs_Target_MSCW_%s.pdf'%(tag))
