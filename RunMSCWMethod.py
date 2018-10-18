
import ExtendedSourceBkgTemplate_Oct10
from ExtendedSourceBkgTemplate_Oct10 import *
import MSCWMethodConfig
from MSCWMethodConfig import *

MakeATag()
EfficiencyRateAtDifferentMSCW()


#Do_analysis = False
Do_analysis = True
if Do_analysis:
    EnergySpectrum("MSCW")
    #SensitivityVsTime("MSCW")

