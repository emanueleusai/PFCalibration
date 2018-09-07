import subprocess

folder='/store/user/eusai/CRAB_UserFiles/step2_dipigun01/180612_020934/0000/'
select='step2_4'

eosls=subprocess.check_output(['eos','root://cmseos.fnal.gov','ls',folder])
names=eosls.splitlines()
for name in names:
	if select in name:
		command=['eos','root://cmseos.fnal.gov','rm',folder+name]
		print command
		print subprocess.check_output(command)
