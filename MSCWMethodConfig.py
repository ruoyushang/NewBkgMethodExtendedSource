
import ROOT
import sys,os
import array
from math import *
from ROOT import *
from array import *

ROOT.gStyle.SetOptStat(0)
ROOT.TH1.SetDefaultSumw2()
ROOT.gStyle.SetPaintTextFormat("0.3f")

Hist_TelElev_Counts = TH1D("Hist_TelElev_Counts","",45,0,90)
Hist_TelElev_Counts_target = TH1D("Hist_TelElev_Counts_target","",18,0,90)
Hist_TelElev_Counts_source = TH1D("Hist_TelElev_Counts_source","",18,0,90)
Hist_TelElev_Counts_scale = TH1D("Hist_TelElev_Counts_scale","",18,0,90)

tag = ''

FOV = 1.0
Ring_Inner = FOV+0.5
Ring_Outer = FOV+1.0

MyRandom = TRandom()

Hist_Norm_Data = TH1D("Hist_Norm_Data","",1,0,1)
Hist_Norm_Ring = TH1D("Hist_Norm_Ring","",1,0,1)

Elev_Bin = [55,85]
#Elev_Bin = [30,40,50,60,70,80,90]
Hist_Elev_Bins = TH1D("Hist_Elev_Bins","",len(Elev_Bin)-1,array('d',Elev_Bin))

time_bins = [1,2,3,4,5]
Hist_NData_vs_Time = TH1D("Hist_NData_vs_Time","",len(time_bins)-1,array('d',time_bins))
Hist_NCR_vs_Time = TH1D("Hist_NCR_vs_Time","",len(time_bins)-1,array('d',time_bins))
Hist_NRing_vs_Time = TH1D("Hist_NRing_vs_Time","",len(time_bins)-1,array('d',time_bins))
Hist_Sensitivity_vs_Time = TH1D("Hist_Sensitivity_vs_Time","",len(time_bins)-1,array('d',time_bins))

ControlWidth_Bin = [3,4]
Hist_ControlWidth_Bins = TH1D("Hist_ControlWidth_Bins","",len(ControlWidth_Bin)-1,array('d',ControlWidth_Bin))

#energy_bins = [pow(10,0),pow(10,0.2),pow(10,0.4),pow(10,0.6),pow(10,0.8),pow(10,1.0)]
#ControlWidthAtThisEnergy = [3,3,3,3,3,3]
energy_bins = [pow(10,0),pow(10,0.2)]
ControlWidthAtThisEnergy = [3,3]

MSCW_Bin = [0,1]
Hist_OffData_Signal = TH1D("Hist_OffData_Signal","",len(MSCW_Bin)-1,array('d',MSCW_Bin))
Hist_OffData_Control = TH1D("Hist_OffData_Control","",len(MSCW_Bin)-1,array('d',MSCW_Bin))
Hist_CR_Efficiency = TH1D("Hist_CR_Efficiency","",len(MSCW_Bin)-1,array('d',MSCW_Bin))

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

Xcenter = 0
Ycenter = 0
Window = 2.
Hist2D_Xoff_vs_Yoff_Data = TH2D("Hist2D_Xoff_vs_Yoff_Data","",40,Xcenter-Window,Xcenter+Window,40,Ycenter-Window,Ycenter+Window)
Hist2D_Xoff_vs_Yoff_Data_AddSignal = TH2D("Hist2D_Xoff_vs_Yoff_Data_AddSignal","",40,Xcenter-Window,Xcenter+Window,40,Ycenter-Window,Ycenter+Window)
Hist2D_Xoff_vs_Yoff_CR = TH2D("Hist2D_Xoff_vs_Yoff_CR","",40,Xcenter-Window,Xcenter+Window,40,Ycenter-Window,Ycenter+Window)
Hist2D_Xoff_vs_Yoff_CR_Raw = TH2D("Hist2D_Xoff_vs_Yoff_CR_Raw","",40,Xcenter-Window,Xcenter+Window,40,Ycenter-Window,Ycenter+Window)
Hist2D_Xoff_vs_Yoff_CR_Raw_AddSignal = TH2D("Hist2D_Xoff_vs_Yoff_CR_Raw_AddSignal","",40,Xcenter-Window,Xcenter+Window,40,Ycenter-Window,Ycenter+Window)
Hist2D_Xoff_vs_Yoff_Data_Sum = TH2D("Hist2D_Xoff_vs_Yoff_Data_Sum","",40,Xcenter-Window,Xcenter+Window,40,Ycenter-Window,Ycenter+Window)
Hist2D_Xoff_vs_Yoff_Data_Sum_AddSignal = TH2D("Hist2D_Xoff_vs_Yoff_Data_Sum_AddSignal","",40,Xcenter-Window,Xcenter+Window,40,Ycenter-Window,Ycenter+Window)
Hist2D_Xoff_vs_Yoff_CR_Sum = TH2D("Hist2D_Xoff_vs_Yoff_CR_Sum","",40,Xcenter-Window,Xcenter+Window,40,Ycenter-Window,Ycenter+Window)
Hist2D_Xoff_vs_Yoff_Sig = TH2D("Hist2D_Xoff_vs_Yoff_Sig","",40,Xcenter-Window,Xcenter+Window,40,Ycenter-Window,Ycenter+Window)
Hist_theta2_Data = TH1D("Hist_theta2_Data","",20,0,0.5)
Hist_theta2_Data_AddSignal = TH1D("Hist_theta2_Data_AddSignal","",20,0,0.5)
Hist_theta2_CR = TH1D("Hist_theta2_CR","",20,0,0.5)
Hist_theta2_CR_Raw = TH1D("Hist_theta2_CR_Raw","",20,0,0.5)
Hist_theta2_CR_Raw_AddSignal = TH1D("Hist_theta2_CR_Raw_AddSignal","",20,0,0.5)
Hist_theta2_Data_Sum = TH1D("Hist_theta2_Data_Sum","",20,0,0.5)
Hist_theta2_Data_Sum_AddSignal = TH1D("Hist_theta2_Data_Sum_AddSignal","",20,0,0.5)
Hist_theta2_CR_Sum = TH1D("Hist_theta2_CR_Sum","",20,0,0.5)
Hist_theta2_Ring_Sum = TH1D("Hist_theta2_Ring_Sum","",20,0,0.5)

variable_bins = []
b = -1
while (b<2):
    variable_bins += [pow(10,b)]
    b += 0.2
Hist2D_Erec_vs_EmissionHeight_Data_All = TH2D("Hist2D_Erec_vs_EmissionHeight_Data_All","",len(variable_bins)-1,array('d',variable_bins),40,0,40)
Hist1D_Erec_vs_EmissionHeight_Data_All = TH1D("Hist1D_Erec_vs_EmissionHeight_Data_All","",len(variable_bins)-1,array('d',variable_bins))
Hist2D_TelElev_vs_EmissionHeight_Data_All = TH2D("Hist2D_TelElev_vs_EmissionHeight_Data_All","",18,0,90,40,0,40)
Hist1D_TelElev_vs_EmissionHeight_Data_All = TH1D("Hist1D_TelElev_vs_EmissionHeight_Data_All","",18,0,90)

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
target_field = 'on'
#target_field = 'off'

Height_cut = ''
RW_method = ''
Fake_Removal = ''
Height_cut = ''
#RW_method = 'NoRW'
#Fake_Removal = 'NoFakeRemoved'


