import subprocess

# folder='/store/user/eusai/CRAB_PrivateMC/step1_pigun01/180417_002722/0000/'
# outname='step1_pigun01_list'
# folder='/store/user/eusai/CRAB_PrivateMC/step1_pigun02/180417_010257/0000/'
# outname='step1_pigun02_list'
# folder='/store/user/eusai/CRAB_PrivateMC/step1_pigun03/180417_010544/0000/'
# outname='step1_pigun03_list'
# folder='/store/user/eusai/CRAB_PrivateMC/step1_pigun04/180417_010927/0000/'
# outname='step1_pigun04_list'
# folder='/store/user/eusai/CRAB_PrivateMC/step1_pigun05/180417_012449/0000/'
# outname='step1_pigun05_list'
# folder='/store/user/eusai/CRAB_UserFiles/step3_pigun01/180421_014640/0000/'
# outname='step3_pigun01_list'

# folder='/store/user/eusai/CRAB_PrivateMC/step1_dipigun01_0p5_0p2/180609_011839/0000/'
# outname='step1_dipigun01_0p5_0p2_list'
# folder='/store/user/eusai/CRAB_PrivateMC/step1_dipigun01_0p2_0p1/180609_013428/0000/'
# outname='step1_dipigun01_0p2_0p1_list'
# folder='/store/user/eusai/CRAB_PrivateMC/step1_dipigun01_0p1_0p05/180609_015005/0000/'
# outname='step1_dipigun01_0p1_0p05_list'
# folder='/store/user/eusai/CRAB_PrivateMC/step1_dipigun01_0p05_0p02/180609_015654/0000/'
# outname='step1_dipigun01_0p05_0p02_list'
# folder='/store/user/eusai/CRAB_PrivateMC/step1_dipigun01_0p02_0p01/180609_015810/0000/'
# outname='step1_dipigun01_0p02_0p01_list'
# folder='/store/user/eusai/CRAB_PrivateMC/step1_dipigun01_0p01_0p00/180609_015908/0000/'
# outname='step1_dipigun01_0p01_0p00_list'

folder='/store/user/eusai/CRAB_UserFiles/step2_dipigun01/180612_020934/0001/'
outname='step2_dipigun01_list2'

xrootd='root://cmsxrootd-site.fnal.gov/'
eosls=subprocess.check_output(['eos','root://cmseos.fnal.gov','ls',folder])
names=eosls.splitlines()
if 'log' in names:
	names.remove('log')
if 'failed' in names:
	names.remove('failed')
outfile=open(outname,'w')
for i in names:
	outfile.write(xrootd+folder+i+'\n')
outfile.close()