import ROOT
import array
ROOT.gROOT.SetBatch()
ROOT.gStyle.SetOptStat(0)
file = ROOT.TFile('file:step4_hadd.root')
#tree = file.Get("s")
outfile=ROOT.TFile('file:out.root','RECREATE')

ROOT.ROOT.EnableImplicitMT()
d=ROOT.Experimental.TDataFrame("s", file)

def plot1D(name,variable,xaxis,yaxis,nbins,xlow,xhigh):
	h_p=d.Define(name, variable).Histo1D(ROOT.Experimental.TDF.TH1DModel(name+"_p",";"+xaxis+";"+yaxis,nbins,xlow,xhigh),name)
	h_p.Write()
	h_c=ROOT.TCanvas(name+"_c", "", 10, 10, 700, 500)
	h_p.Draw()
	h_c.SaveAs(name+'.pdf')

plot1D("hcal1","(hcal_energy-hcal_energy_raw)/hcal_energy","(hcal_energy-hcal_energy_raw)/hcal_energy","N",100,0,5)

# eta = d.Histo1D("eta")
# eta.Write()

# hcal1_p=d.Define("hcal1", "(hcal_energy-hcal_energy_raw)/hcal_energy").Histo1D(ROOT.Experimental.TH1DModel("hcal1_p"),"hcal1")
# hcal1_p.Write()
# hcal1_c=ROOT.TCanvas("hcal1_c", "", 10, 10, 700, 500)
# hcal1_p.Draw()


outfile.Close()
# for event in tree:
# 	print event.p
