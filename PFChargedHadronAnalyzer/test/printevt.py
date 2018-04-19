from DataFormats.FWLite import Handle, Events
import ROOT

basename='DY'
events = Events('/uscms_data/d3/eusai/'+basename+'.root')
aux = Handle("edm::EventAuxiliary")
files=[]
step=3000
ranges=range(0,events.size(),step)
for i in range(len(ranges)):
	files.append(open(basename+'_'+str(i+1)+'.txt', 'w'))
i=0
for event in events:
	i=i+1
	for j in range(len(ranges)):
		if i>ranges[j] and i<=ranges[j]+step:
			files[j].write(str(event.eventAuxiliary().run())+':'+str(event.eventAuxiliary().event())+'\n')
#':'+str(event.eventAuxiliary().luminosityBlock())+
for i in range(len(ranges)):
	files[i].close()
	print 'edmCopyPickMerge outputFile=/uscms_data/d3/eusai/'+basename+'_'+str(i+1)+'.root eventsToProcess_load='+basename+'_'+str(i+1)+'.txt inputFiles=file:///uscms_data/d3/eusai/'+basename+'.root'
    

