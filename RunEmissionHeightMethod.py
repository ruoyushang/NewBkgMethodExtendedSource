
import ExtendedSourceBkgTemplate_Oct10
from ExtendedSourceBkgTemplate_Oct10 import *
import EmissionHeightMethodConfig
from EmissionHeightMethodConfig import *

MakeATag()
AttenuationRateAtDifferentHeight()

#Do_analysis = False
Do_analysis = True
if Do_analysis:
    EnergySpectrum("Height")
    #SensitivityVsTime("Height")
