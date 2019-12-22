
energy_list = []
energy_list += [200]
energy_list += [237]
energy_list += [282]
energy_list += [335]
energy_list += [398]
energy_list += [473]
energy_list += [562]
energy_list += [667]
energy_list += [794]
energy_list += [943]
energy_list += [1122]
energy_list += [1332]
energy_list += [1585]
energy_list += [1882]
energy_list += [2239]
energy_list += [3162]
energy_list += [4467]
energy_list += [6310]
energy_list += [8913]

source = []
for_syst = []

rank = 4
MSCW_cut = 0.4
MSCL_cut = 0.8
NSB_diff = 0.
Elev_diff = 0.

Theta2_cut = 1.0

folder = 'output_unblind_4x4_nominal_tight'
rank = 4
MSCW_cut = 0.3
MSCL_cut = 0.5
NSB_diff = 0.
Elev_diff = 0.
Theta2_cut = 10.
#folder = 'output_unblind_4x4_nominal_medium'
#rank = 4
#MSCW_cut = 0.5
#MSCL_cut = 0.7
#NSB_diff = 0.
#Elev_diff = 0.
#Theta2_cut = 10.
#folder = 'output_unblind_4x4_nominal_loose'
#rank = 4
#MSCW_cut = 1.0
#MSCL_cut = 1.0
#NSB_diff = 0.
#Elev_diff = 0.
#Theta2_cut = 10.
#folder = 'output_unblind_4x4_dNSBp3_loose'
#rank = 4
#MSCW_cut = 1.0
#MSCL_cut = 1.0
#NSB_diff = 3.
#Elev_diff = 0.
#Theta2_cut = 10.
#folder = 'output_unblind_4x4_dNSBm3_loose'
#rank = 4
#MSCW_cut = 1.0
#MSCL_cut = 1.0
#NSB_diff = -3.
#Elev_diff = 0.
#Theta2_cut = 10.
#folder = 'output_unblind_4x4_dElevm20_loose'
#rank = 4
#MSCW_cut = 1.0
#MSCL_cut = 1.0
#NSB_diff = 0.
#Elev_diff = -20.
#Theta2_cut = 10.
#folder = 'output_unblind_4x4_dElevm10_loose'
#rank = 4
#MSCW_cut = 1.0
#MSCL_cut = 1.0
#NSB_diff = 0.
#Elev_diff = -10.
#Theta2_cut = 10.
#folder = 'output_unblind_4x4_smallfov_tight'
#rank = 4
#MSCW_cut = 0.3
#MSCL_cut = 0.5
#NSB_diff = 0.
#Elev_diff = 0.
#Theta2_cut = 1.
#folder = 'output_test'
#rank = 4
#MSCW_cut = 0.3
#MSCL_cut = 0.5
#NSB_diff = 0.
#Elev_diff = 0.
#Theta2_cut = 10.

source += ['Segue1V6']
for_syst += [True]
source += ['IC443HotSpot']
for_syst += [True]
source += ['Mrk421']
for_syst += [True]
source += ['1ES0229']
for_syst += [True]
source += ['PKS1424']
for_syst += [True]
source += ['H1426']
for_syst += [True]
source += ['3C264']
for_syst += [True]
source += ['Crab']
for_syst += [True]
source += ['G079']
for_syst += [True]
source += ['M82']
for_syst += [True]
source += ['CasA']
for_syst += [True]
source += ['RGBJ0710']
for_syst += [True]
source += ['OJ287V6']
for_syst += [True]
source += ['1ES1011V6']
for_syst += [True]
source += ['WComaeV6']
for_syst += [True]
source += ['1ES1218V6']
for_syst += [False]
source += ['NGC1275V6']
for_syst += [True]
source += ['1ES0647V6']
for_syst += [True]
source += ['1ES1440V6']
for_syst += [True]
source += ['1ES1741V6']
for_syst += [True]
source += ['RBS0413V6']
for_syst += [True]
source += ['PG1553V6']
for_syst += [True]
source += ['S3_1227_V6']
for_syst += [True]
source += ['MS1221V6']
for_syst += [True]
source += ['PKS1441V6']
for_syst += [True]
source += ['MGRO_J1908_V6']
for_syst += [False]
#source += ['Segue1V5']
#for_syst += [False]
source += ['MGRO_J1908_V5']
for_syst += [False]
source += ['IC443HotSpotV5']
for_syst += [False]
#source += ['VA_Segue1']
#for_syst += [False]
source += ['Proton_NSB200']
for_syst += [False]
source += ['Proton_NSB750']
for_syst += [False]
source += ['GemingaV6']
for_syst += [False]
source += ['GemingaV5']
for_syst += [False]
source += ['Everything']
for_syst += [True]

elev = []
elev += [85,75,65,55,45,35,25]

#gamma = []
#gamma += [0,10,20,50,100,200]
gamma = []
gamma += [0]

#if folder=='output_test':
#    elev = []
#    elev += [85,75,65]
if not 'nominal' in folder:
    gamma = []
    gamma += [0]

for s in range(0,len(source)):
    file = open("run/run_RDBM_%s_LargeOFF.sh"%(source[s]),"w") 
    file.write('cd /home/rshang/EventDisplay/NewBkgMethodExtendedSource\n')
    file.write("root -b -l -q 'GetShowerImageHistograms.C+(\"%s\",70,85,0.2,10.0)'\n"%(source[s])) 
    file.write("root -b -l -q 'RecurrentDeconvolution.C+(\"%s\",70,85,0.2,10.0)'\n"%(source[s])) 
    file.write("root -b -l -q 'GetShowerImageHistograms.C+(\"%s\",55,70,0.2,10.0)'\n"%(source[s])) 
    file.write("root -b -l -q 'RecurrentDeconvolution.C+(\"%s\",55,70,0.2,10.0)'\n"%(source[s])) 
    file.write("root -b -l -q 'GetShowerImageHistograms.C+(\"%s\",40,55,0.2,10.0)'\n"%(source[s])) 
    file.write("root -b -l -q 'RecurrentDeconvolution.C+(\"%s\",40,55,0.2,10.0)'\n"%(source[s])) 
    file.close() 
    file = open("run/run_RDBM_%s_LargeON.sh"%(source[s]),"w") 
    file.write('cd /home/rshang/EventDisplay/NewBkgMethodExtendedSource\n')
    file.write("root -b -l -q 'GetShowerImageHistograms.C+(\"%s\",70,85,0.0,10.0)'\n"%(source[s])) 
    file.write("root -b -l -q 'RecurrentDeconvolution.C+(\"%s\",70,85,0.0,10.0)'\n"%(source[s])) 
    file.write("root -b -l -q 'GetShowerImageHistograms.C+(\"%s\",55,70,0.0,10.0)'\n"%(source[s])) 
    file.write("root -b -l -q 'RecurrentDeconvolution.C+(\"%s\",55,70,0.0,10.0)'\n"%(source[s])) 
    file.write("root -b -l -q 'GetShowerImageHistograms.C+(\"%s\",40,55,0.0,10.0)'\n"%(source[s])) 
    file.write("root -b -l -q 'RecurrentDeconvolution.C+(\"%s\",40,55,0.0,10.0)'\n"%(source[s])) 
    file.close() 
    file = open("run/run_RDBM_%s_SmallON.sh"%(source[s]),"w") 
    file.write('cd /home/rshang/EventDisplay/NewBkgMethodExtendedSource\n')
    file.write("root -b -l -q 'GetShowerImageHistograms.C+(\"%s\",70,85,0.0,0.2)'\n"%(source[s])) 
    file.write("root -b -l -q 'RecurrentDeconvolution.C+(\"%s\",70,85,0.0,0.2)'\n"%(source[s])) 
    file.write("root -b -l -q 'GetShowerImageHistograms.C+(\"%s\",55,70,0.0,0.2)'\n"%(source[s])) 
    file.write("root -b -l -q 'RecurrentDeconvolution.C+(\"%s\",55,70,0.0,0.2)'\n"%(source[s])) 
    file.write("root -b -l -q 'GetShowerImageHistograms.C+(\"%s\",40,55,0.0,0.2)'\n"%(source[s])) 
    file.write("root -b -l -q 'RecurrentDeconvolution.C+(\"%s\",40,55,0.0,0.2)'\n"%(source[s])) 
    file.close() 

for s in range(0,len(source)):
    file = open("run/run_Netflix_%s.sh"%(source[s]),"w") 
    file.write('cd /home/rshang/EventDisplay/NewBkgMethodExtendedSource\n')
    file.write('rm -r %s/%s_LOFF\n'%(folder,source[s]))
    file.write('mkdir %s/%s_LOFF\n'%(folder,source[s]))
    file.write('cp GetRunList.h %s/%s_LOFF\n'%(folder,source[s]))
    file.write('cp NetflixMethodGetShowerImage.C %s/%s_LOFF\n'%(folder,source[s]))
    file.write('cp NetflixMethodPrediction.C %s/%s_LOFF\n'%(folder,source[s]))
    file.write('rm -r %s/%s_LON\n'%(folder,source[s]))
    file.write('mkdir %s/%s_LON\n'%(folder,source[s]))
    file.write('cp GetRunList.h %s/%s_LON\n'%(folder,source[s]))
    file.write('cp NetflixMethodGetShowerImage.C %s/%s_LON\n'%(folder,source[s]))
    file.write('cp NetflixMethodPrediction.C %s/%s_LON\n'%(folder,source[s]))
    file.write('cd /home/rshang/EventDisplay/NewBkgMethodExtendedSource\n')
    file.write('rm -r %s/%s_LPHO\n'%(folder,source[s]))
    file.write('mkdir %s/%s_LPHO\n'%(folder,source[s]))
    file.write('cp GetRunList.h %s/%s_LPHO\n'%(folder,source[s]))
    file.write('cp NetflixMethodGetShowerImage.C %s/%s_LPHO\n'%(folder,source[s]))
    file.write('cp NetflixMethodPrediction.C %s/%s_LPHO\n'%(folder,source[s]))
    file.write('cd /home/rshang/EventDisplay/NewBkgMethodExtendedSource\n')
    file.write('cd %s/%s_LOFF\n'%(folder,source[s]))
    file.write('rm *_C*\n')
    for g in range(0,len(gamma)):
        for e in range(0,len(elev)-1):
            if 'Proton' in source[s]: continue
            if for_syst[s]==False: continue
            file.write("root -b -l -q 'NetflixMethodGetShowerImage.C+(\"%s\",%s,%s,%s,%s,%s,false,%s,%s,%s)'\n"%(source[s],gamma[g],elev[e+1],elev[e],Elev_diff,NSB_diff,MSCW_cut,MSCL_cut,Theta2_cut)) 
            file.write("root -b -l -q 'NetflixMethodPrediction.C+(\"%s\",%s,%s,%s,false,%s,%s,%s)'\n"%(source[s],gamma[g],elev[e+1],elev[e],MSCW_cut,MSCL_cut,rank))
    file.write('cd /home/rshang/EventDisplay/NewBkgMethodExtendedSource\n')
    file.write('cd %s/%s_LON\n'%(folder,source[s]))
    file.write('rm *_C*\n')
    for e in range(0,len(elev)-1):
        file.write("root -b -l -q 'NetflixMethodGetShowerImage.C+(\"%s\",0,%s,%s,%s,%s,true,%s,%s,%s)'\n"%(source[s],elev[e+1],elev[e],Elev_diff,NSB_diff,MSCW_cut,MSCL_cut,Theta2_cut)) 
        file.write("root -b -l -q 'NetflixMethodPrediction.C+(\"%s\",0,%s,%s,true,%s,%s,%s)'\n"%(source[s],elev[e+1],elev[e],MSCW_cut,MSCL_cut,rank))
    file.close() 

qfile = open("run/qsub_RDBM.sh","w") 
for s in range(0,len(source)):
    qfile.write('qsub -V -N job_%s_LOFF run_RDBM_%s_LargeOFF.sh\n'%(source[s],source[s]))
    qfile.write('qsub -V -N job_%s_SON run_RDBM_%s_SmallON.sh\n'%(source[s],source[s]))
    qfile.write('qsub -V -N job_%s_LON run_RDBM_%s_LargeON.sh\n'%(source[s],source[s]))
qfile.close() 

qfile = open("run/qsub_Netflix.sh","w") 
for s in range(0,len(source)):
    qfile.write('qsub -V -N job_%s run_Netflix_%s.sh\n'%(source[s],source[s]))
qfile.close() 

