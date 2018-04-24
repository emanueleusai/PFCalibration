import ROOT

file = ROOT.TFile('file:step4.root')
#tree = file.Get("s")

ROOT.ROOT.EnableImplicitMT()
d=ROOT.Experimental.TDataFrame("s", file)
myHisto = d.Histo1D("eta")
myHisto.Draw()
# for event in tree:
# 	print event.p
