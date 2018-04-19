l800=[]
tipo=['EG1p33X','EG2X','G1p33X','G2X','G1X','G4X','G6X','EG4X','EG6X','nofilter']
pu='0'
index_tipo=9
hadd800='hadd /uscms_data/d3/eusai/WZ_M800_PU'+pu+'_'+tipo[index_tipo]+'_new_hadd.root '
for i in range(10):
	l800.append("800_"+str(i+1))
	hadd800=hadd800+'root://cmseos.fnal.gov//store/user/eusai/WZ_M800_'+str(i+1)+'_PU'+pu+'_'+tipo[index_tipo]+'_new.root '
	print "xrdcp -f /uscms_data/d3/eusai/800_"+str(i+1)+".root root://cmseos.fnal.gov//store/user/eusai/800_"+str(i+1)+".root"
l1600=[]
hadd1600='hadd /uscms_data/d3/eusai/WZ_M1600_PU'+pu+'_'+tipo[index_tipo]+'_new_hadd.root '
for i in range(10):
	l1600.append("1600_"+str(i+1))
	hadd1600=hadd1600+'root://cmseos.fnal.gov//store/user/eusai/WZ_M1600_'+str(i+1)+'_PU'+pu+'_'+tipo[index_tipo]+'_new.root '
	print "xrdcp -f /uscms_data/d3/eusai/1600_"+str(i+1)+".root root://cmseos.fnal.gov//store/user/eusai/1600_"+str(i+1)+".root"
lDY=[]
haddDY='hadd /uscms_data/d3/eusai/WZ_MDY_PU'+pu+'_'+tipo[index_tipo]+'_new_hadd.root '
for i in range(38):
	lDY.append("DY_"+str(i+1))
	haddDY=haddDY+'root://cmseos.fnal.gov//store/user/eusai/WZ_MDY_'+str(i+1)+'_PU'+pu+'_'+tipo[index_tipo]+'_new.root '
	print "xrdcp -f /uscms_data/d3/eusai/DY_"+str(i+1)+".root root://cmseos.fnal.gov//store/user/eusai/DY_"+str(i+1)+".root"

print l800
print '\n'
print l1600
print '\n'
print lDY
print '\n'
print hadd800
print '\n'
print hadd1600
print '\n'
print haddDY