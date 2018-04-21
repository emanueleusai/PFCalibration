import subprocess

folder='/store/user/eusai/'
select='crab3checkwrite'

eosls=subprocess.check_output(['eos','root://cmseos.fnal.gov','ls',folder])
names=eosls.splitlines()
for name in names:
	if select in name:
		command=['eos','root://cmseos.fnal.gov','rm',folder+name]
		print command
		print subprocess.check_output(command)
