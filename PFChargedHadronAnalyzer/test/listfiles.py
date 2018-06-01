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

folder='/store/user/eusai/CRAB_UserFiles/step4_pigun01/180424_025713/0000/'
outname='step4_pigun01_list'

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