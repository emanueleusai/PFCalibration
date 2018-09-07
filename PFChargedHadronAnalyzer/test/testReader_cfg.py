# #!/usr/bin/env cmsRun
# import FWCore.ParameterSet.Config as cms
# import sys

# #name=sys.argv[2]

# process = cms.Process("LHE")

# process.source = cms.Source("LHESource",
# 	fileNames = cms.untracked.vstring('file:/uscms/homes/e/eusai/cmsgrid_final.lhe')
# )

# process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))


# process.LHE = cms.OutputModule("PoolOutputModule",
# 	dataset = cms.untracked.PSet(dataTier = cms.untracked.string('LHE')),
# 	fileName = cms.untracked.string('lhe.root')
# )


#!/usr/bin/env cmsRun
import FWCore.ParameterSet.Config as cms

process = cms.Process("LHE")

process.source = cms.Source("LHESource",
	fileNames = cms.untracked.vstring('file:/uscms_data/d3/eusai/27tev/lhe/2.lhe'),
        processCode = cms.int32(-11361)                       

)

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.LHE = cms.OutputModule("PoolOutputModule",
	dataset = cms.untracked.PSet(dataTier = cms.untracked.string('LHE')),
	fileName = cms.untracked.string('Zeventslhe.root')
)

process.outpath = cms.EndPath(process.LHE)

