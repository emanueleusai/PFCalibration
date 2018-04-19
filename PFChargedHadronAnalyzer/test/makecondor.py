#infiles=['800','1600']
#infiles=['DY']
#infiles=['800_1', '800_2', '800_3', '800_4', '800_5', '800_6', '800_7', '800_8', '800_9', '800_10']
infiles=['1600_1', '1600_2', '1600_3', '1600_4', '1600_5', '1600_6', '1600_7', '1600_8', '1600_9', '1600_10']
#infiles=['DY_1', 'DY_2', 'DY_3', 'DY_4', 'DY_5', 'DY_6', 'DY_7', 'DY_8', 'DY_9', 'DY_10', 'DY_11', 'DY_12', 'DY_13', 'DY_14', 'DY_15', 'DY_16', 'DY_17', 'DY_18', 'DY_19', 'DY_20', 'DY_21', 'DY_22', 'DY_23', 'DY_24', 'DY_25', 'DY_26', 'DY_27', 'DY_28', 'DY_29', 'DY_30', 'DY_31', 'DY_32', 'DY_33', 'DY_34', 'DY_35', 'DY_36', 'DY_37', 'DY_38']
cards=[

#['CMS_PhaseII_0PU_v03_HG_1p33_new.tcl','_PU0_G1p33X_new'],
# ['CMS_PhaseII_200PU_v03_HG_1p33_new.tcl','_PU200_G1p33X_new'],
# ['CMS_PhaseII_200PU_v03_HG_2x_new.tcl','_PU200_G2X_new'],
# ['CMS_PhaseII_200PU_v03_HGE_1p33_new.tcl','_PU200_EG1p33X_new'],
# ['CMS_PhaseII_200PU_v03_HGE_2x_new.tcl','_PU200_EG2X_new'],
# ['CMS_PhaseII_200PU_v03_new.tcl','_PU200_G1X_new'],

# ['CMS_PhaseII_200PU_v03_HG_4x_new.tcl','_PU200_G4X_new'],
# ['CMS_PhaseII_200PU_v03_HG_6x_new.tcl','_PU200_G6X_new'],

# ['CMS_PhaseII_200PU_v03_HGE_4x_new.tcl','_PU200_EG4X_new'],
# ['CMS_PhaseII_200PU_v03_HGE_6x_new.tcl','_PU200_EG6X_new'],
['CMS_PhaseII_0PU_v03_nofilter_new.tcl','_PU0_nofilter_new'],

]
commands=[]
for card in cards:
    for infile in infiles:
        filename='WZ_M'+infile+card[1]
        file=open(filename+'.jdl','w')
        file.write(
'universe = vanilla\n\
Executable = run.sh\n\
Requirements = OpSys == "LINUX" && (Arch != "DUMMY" )\n\
request_disk = 10000000\n\
request_memory = 2100\n\
Should_Transfer_Files = YES\n\
WhenToTransferOutput = ON_EXIT\n\
Transfer_Input_Files = run.sh,delphes.tar.gz,'+card[0]+'\n\
Output = $(Cluster)_$(Process).stdout\n\
Error = $(Cluster)_$(Process).stderr\n\
Log = $(Cluster)_$(Process).condor\n\
notification = Never\n\
x509userproxy = $ENV(X509_USER_PROXY)\n\
Arguments = '+infile+'.root '+filename+'.root root://cmseos.fnal.gov//store/user/eusai/ '+card[0]+'\n\
Queue 1\n')
        file.close()
        commands.append('condor_submit '+filename+'.jdl')
for i in commands:
    print i