
import ROOT
import sys,os
import array
from math import *
from ROOT import *
from array import *

ROOT.gStyle.SetOptStat(0)
ROOT.TH1.SetDefaultSumw2()
ROOT.gStyle.SetPaintTextFormat("0.3f")

Hist_TelElev_Counts_target = TH1D("Hist_TelElev_Counts_target","",18,0,90)
Hist_TelElev_Counts_source = TH1D("Hist_TelElev_Counts_source","",18,0,90)
Hist_TelElev_Counts_scale = TH1D("Hist_TelElev_Counts_scale","",18,0,90)
Hist_TelAzim_Counts_target = TH1D("Hist_TelAzim_Counts_target","",18,0,180)
Hist_TelAzim_Counts_source = TH1D("Hist_TelAzim_Counts_source","",18,0,180)
Hist_TelAzim_Counts_scale = TH1D("Hist_TelAzim_Counts_scale","",18,0,180)
Hist_TelElevAzim_Counts_target = TH2D("Hist_TelElevAzim_Counts_target","",18,0,90,18,0,180)
Hist_TelElevAzim_Counts_source = TH2D("Hist_TelElevAzim_Counts_source","",18,0,90,18,0,180)
Hist_TelElevAzim_Counts_scale = TH2D("Hist_TelElevAzim_Counts_scale","",18,0,90,18,0,180)
Hist_TelElev_Counts_target_Sum = TH1D("Hist_TelElev_Counts_target_Sum","",18,0,90)
Hist_TelElev_Counts_source_Sum = TH1D("Hist_TelElev_Counts_source_Sum","",18,0,90)
Hist_TelAzim_Counts_target_Sum = TH1D("Hist_TelAzim_Counts_target_Sum","",18,0,180)
Hist_TelAzim_Counts_source_Sum = TH1D("Hist_TelAzim_Counts_source_Sum","",18,0,180)

tag = ''
tag_method = 'EmissionHeight'

FOV = 1.0
Ring_Inner = FOV+0.5
Ring_Outer = FOV+1.0

MyRandom = TRandom()

Hist_Norm_Data = TH1D("Hist_Norm_Data","",1,0,1)
Hist_Norm_Ring = TH1D("Hist_Norm_Ring","",1,0,1)

Elev_Bin = [50,60,70,80,90]
#Elev_Bin = [70,80]
Hist_Elev_Bins = TH1D("Hist_Elev_Bins","",len(Elev_Bin)-1,array('d',Elev_Bin))
Azim_Bin = [0,20,40,60,80]
#Azim_Bin = [60,80]
Hist_Azim_Bins = TH1D("Hist_Azim_Bins","",len(Azim_Bin)-1,array('d',Azim_Bin))

time_bins = [1,2,3,4,5,6,7,8,9,10,11]
Hist_NData_vs_Time = TH1D("Hist_NData_vs_Time","",len(time_bins)-1,array('d',time_bins))
Hist_NCR_vs_Time = TH1D("Hist_NCR_vs_Time","",len(time_bins)-1,array('d',time_bins))
Hist_NRing_vs_Time = TH1D("Hist_NRing_vs_Time","",len(time_bins)-1,array('d',time_bins))
Hist_Sensitivity_vs_Time = TH1D("Hist_Sensitivity_vs_Time","",len(time_bins)-1,array('d',time_bins))
Hist_Sensitivity_vs_Time_Other = TH1D("Hist_Sensitivity_vs_Time_Other","",len(time_bins)-1,array('d',time_bins))

#ControlWidth_Bin = [3,4,5,6]
ControlWidth_Bin = [4,5]
#ControlWidth_Bin = [1,2,3,4,5,6,7,8,9,10]
Hist_ControlWidth_Bins = TH1D("Hist_ControlWidth_Bins","",len(ControlWidth_Bin)-1,array('d',ControlWidth_Bin))

#energy_bins = [pow(10,-0.9),pow(10,-0.8),pow(10,-0.7),pow(10,-0.6),pow(10,-0.5)]
##ControlWidthAtThisEnergy = [3,3,3,4,4]
#ControlWidthAtThisEnergy = [4,4,4,4,4]
#energy_bins += [pow(10,-0.4),pow(10,-0.3),pow(10,-0.2),pow(10,-0.1),pow(10,0),pow(10,0.1)]
##ControlWidthAtThisEnergy += [4,4,5,5,5,5]
#ControlWidthAtThisEnergy += [4,4,4,4,4,4]

energy_bins = [0.16,0.2,0.25,0.32]
ControlWidthAtThisEnergy = [4,4,4,4]

#MSCW_Bin = [-100,3,100]  # the last bin has to be the CR normalization region!!
MSCW_Bin = [-100,100]  # the last bin has to be the CR normalization region!!
Hist_OnData_Signal = TH1D("Hist_OnData_Signal","",len(MSCW_Bin)-1,array('d',MSCW_Bin))
Hist_OnData_Control = TH1D("Hist_OnData_Control","",len(MSCW_Bin)-1,array('d',MSCW_Bin))
Hist_OffData_Signal = TH1D("Hist_OffData_Signal","",len(MSCW_Bin)-1,array('d',MSCW_Bin))
Hist_OffData_Control = TH1D("Hist_OffData_Control","",len(MSCW_Bin)-1,array('d',MSCW_Bin))
Hist_Gamma_Signal = TH1D("Hist_Gamma_Signal","",len(MSCW_Bin)-1,array('d',MSCW_Bin))
Hist_Gamma_Control = TH1D("Hist_Gamma_Control","",len(MSCW_Bin)-1,array('d',MSCW_Bin))
Hist_CR_Attenuation = TH1D("Hist_CR_Attenuation","",len(MSCW_Bin)-1,array('d',MSCW_Bin))
Hist_Gamma_Attenuation = TH1D("Hist_Gamma_Attenuation","",len(MSCW_Bin)-1,array('d',MSCW_Bin))

Hist_Erec_CR = TH1D("Hist_Erec_CR","",len(energy_bins)-1,array('d',energy_bins))
Hist_Erec_CR_Raw = TH1D("Hist_Erec_CR_Raw","",len(energy_bins)-1,array('d',energy_bins))
Hist_Erec_CR_Raw_AddSignal = TH1D("Hist_Erec_CR_Raw_AddSignal","",len(energy_bins)-1,array('d',energy_bins))
Hist_Erec_Data = TH1D("Hist_Erec_Data","",len(energy_bins)-1,array('d',energy_bins))
Hist_Erec_Data_AddSignal = TH1D("Hist_Erec_Data_AddSignal","",len(energy_bins)-1,array('d',energy_bins))
Hist_Erec_Ring = TH1D("Hist_Erec_Ring","",len(energy_bins)-1,array('d',energy_bins))
Hist_Erec_Data_Sum = TH1D("Hist_Erec_Data_Sum","",len(energy_bins)-1,array('d',energy_bins))
Hist_Erec_Data_Sum_AddSignal = TH1D("Hist_Erec_Data_Sum_AddSignal","",len(energy_bins)-1,array('d',energy_bins))
Hist_Erec_Ring_Sum = TH1D("Hist_Erec_Ring_Sum","",len(energy_bins)-1,array('d',energy_bins))
Hist_Erec_CR_Sum = TH1D("Hist_Erec_CR_Sum","",len(energy_bins)-1,array('d',energy_bins))
Hist_Sensitivity_vs_Erec = TH1D("Hist_Sensitivity_vs_Erec","",len(energy_bins)-1,array('d',energy_bins))
Hist_Sensitivity_vs_Erec_Other = TH1D("Hist_Sensitivity_vs_Erec_Other","",len(energy_bins)-1,array('d',energy_bins))

Xcenter = 0
Ycenter = 0
Window = 2.
Nbins = 20
Hist2D_Xoff_vs_Yoff_Data = TH2D("Hist2D_Xoff_vs_Yoff_Data","",Nbins,Xcenter-Window,Xcenter+Window,Nbins,Ycenter-Window,Ycenter+Window)
Hist2D_Xoff_vs_Yoff_Data_AddSignal = TH2D("Hist2D_Xoff_vs_Yoff_Data_AddSignal","",Nbins,Xcenter-Window,Xcenter+Window,Nbins,Ycenter-Window,Ycenter+Window)
Hist2D_Xoff_vs_Yoff_CR = TH2D("Hist2D_Xoff_vs_Yoff_CR","",Nbins,Xcenter-Window,Xcenter+Window,Nbins,Ycenter-Window,Ycenter+Window)
Hist2D_Xoff_vs_Yoff_CR_Raw = TH2D("Hist2D_Xoff_vs_Yoff_CR_Raw","",Nbins,Xcenter-Window,Xcenter+Window,Nbins,Ycenter-Window,Ycenter+Window)
Hist2D_Xoff_vs_Yoff_CR_Raw_AddSignal = TH2D("Hist2D_Xoff_vs_Yoff_CR_Raw_AddSignal","",Nbins,Xcenter-Window,Xcenter+Window,Nbins,Ycenter-Window,Ycenter+Window)
Hist2D_Xoff_vs_Yoff_Data_Sum = TH2D("Hist2D_Xoff_vs_Yoff_Data_Sum","",Nbins,Xcenter-Window,Xcenter+Window,Nbins,Ycenter-Window,Ycenter+Window)
Hist2D_Xoff_vs_Yoff_Data_Sum_AddSignal = TH2D("Hist2D_Xoff_vs_Yoff_Data_Sum_AddSignal","",Nbins,Xcenter-Window,Xcenter+Window,Nbins,Ycenter-Window,Ycenter+Window)
Hist2D_Xoff_vs_Yoff_CR_Sum = TH2D("Hist2D_Xoff_vs_Yoff_CR_Sum","",Nbins,Xcenter-Window,Xcenter+Window,Nbins,Ycenter-Window,Ycenter+Window)
Hist2D_Xoff_vs_Yoff_Sig = TH2D("Hist2D_Xoff_vs_Yoff_Sig","",Nbins,Xcenter-Window,Xcenter+Window,Nbins,Ycenter-Window,Ycenter+Window)

Hist_theta2_Data = TH1D("Hist_theta2_Data","",20,0,4.0)
Hist_theta2_Data_AddSignal = TH1D("Hist_theta2_Data_AddSignal","",20,0,4.0)
Hist_theta2_CR = TH1D("Hist_theta2_CR","",20,0,4.0)
Hist_theta2_CR_Raw = TH1D("Hist_theta2_CR_Raw","",20,0,4.0)
Hist_theta2_CR_Raw_AddSignal = TH1D("Hist_theta2_CR_Raw_AddSignal","",20,0,4.0)
Hist_theta2_Data_Sum = TH1D("Hist_theta2_Data_Sum","",20,0,4.0)
Hist_theta2_Data_Sum_AddSignal = TH1D("Hist_theta2_Data_Sum_AddSignal","",20,0,4.0)
Hist_theta2_CR_Sum = TH1D("Hist_theta2_CR_Sum","",20,0,4.0)
Hist_theta2_Ring_Sum = TH1D("Hist_theta2_Ring_Sum","",20,0,4.0)

Hist_theta2_zoomin_Data = TH1D("Hist_theta2_zoomin_Data","",20,0,0.4)
Hist_theta2_zoomin_Data_AddSignal = TH1D("Hist_theta2_zoomin_Data_AddSignal","",20,0,0.4)
Hist_theta2_zoomin_CR = TH1D("Hist_theta2_zoomin_CR","",20,0,0.4)
Hist_theta2_zoomin_CR_Raw = TH1D("Hist_theta2_zoomin_CR_Raw","",20,0,0.4)
Hist_theta2_zoomin_CR_Raw_AddSignal = TH1D("Hist_theta2_zoomin_CR_Raw_AddSignal","",20,0,0.4)
Hist_theta2_zoomin_Data_Sum = TH1D("Hist_theta2_zoomin_Data_Sum","",20,0,0.4)
Hist_theta2_zoomin_Data_Sum_AddSignal = TH1D("Hist_theta2_zoomin_Data_Sum_AddSignal","",20,0,0.4)
Hist_theta2_zoomin_CR_Sum = TH1D("Hist_theta2_zoomin_CR_Sum","",20,0,0.4)
Hist_theta2_zoomin_Ring_Sum = TH1D("Hist_theta2_zoomin_Ring_Sum","",20,0,0.4)

variable_bins = []
b = -1
while (b<2):
    variable_bins += [pow(10,b)]
    b += 0.2
Hist2D_Erec_vs_EmissionHeight_Data_All = TH2D("Hist2D_Erec_vs_EmissionHeight_Data_All","",len(variable_bins)-1,array('d',variable_bins),40,0,40)
Hist1D_Erec_vs_EmissionHeight_Data_All = TH1D("Hist1D_Erec_vs_EmissionHeight_Data_All","",len(variable_bins)-1,array('d',variable_bins))
Hist2D_TelElev_vs_EmissionHeight_Data_All = TH2D("Hist2D_TelElev_vs_EmissionHeight_Data_All","",18,0,90,40,0,40)
Hist1D_TelElev_vs_EmissionHeight_Data_All = TH1D("Hist1D_TelElev_vs_EmissionHeight_Data_All","",18,0,90)
Hist2D_TelAzim_vs_EmissionHeight_Data_All = TH2D("Hist2D_TelAzim_vs_EmissionHeight_Data_All","",18,0,180,40,0,40)
Hist1D_TelAzim_vs_EmissionHeight_Data_All = TH1D("Hist1D_TelAzim_vs_EmissionHeight_Data_All","",18,0,180)

source = 'Crab'
field = 'on'

#add_signal = True
add_signal = False

#target = 'Crab'
#target = '2ndCrab'
#target = '3C264'
target = 'PKS1424'
#target = 'H1426'
#target = 'Ton599'
#target = 'IC443'
#target_field = 'on'
target_field = 'off'

Height_cut = ''
RW_method = ''
Fake_Removal = ''
Height_cut = ''
RW_method = 'NoRW'
#Fake_Removal = 'NoFakeRemoved'

