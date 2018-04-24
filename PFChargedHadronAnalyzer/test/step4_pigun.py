# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step3 --conditions auto:phase1_2017_realistic -n 10 --era Run2_2017 --eventcontent RECOSIM,MINIAODSIM,DQM --runUnscheduled -s RAW2DIGI,L1Reco,RECO,EI,PAT,VALIDATION:@standardValidation+@miniAODValidation,DQM:@standardDQM+@miniAODDQM --datatier GEN-SIM-RECO,MINIAODSIM,DQMIO --geometry DB:Extended --conditions 90X_upgrade2017_realistic_v15 --filein file:step2.root --fileout file:step3.root
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras
#from step3_pigun01_list import step3_files

process = cms.Process('RECO2',eras.Run2_2017)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
#process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.Geometry.GeometryExtended2023D17Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('CommonTools.ParticleFlow.EITopPAG_cff')
process.load('PhysicsTools.PatAlgos.slimming.metFilterPaths_cff')
process.load('Configuration.StandardSequences.Validation_cff')
process.load('DQMOffline.Configuration.DQMOfflineMC_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("PoolSource",
    #fileNames = cms.untracked.vstring('file:step2.root'),
    #fileNames = cms.untracked.vstring('file:rereco.root'),
    fileNames = cms.untracked.vstring('file:step3.root'),

    #fileNames = cms.untracked.vstring(*step3_files),

    #fileNames = cms.untracked.vstring('file:E20ADE15-A5AE-E711-9434-0023AEEEB55F.root'),
    #fileNames = cms.untracked.vstring('root://se01.indiacms.res.in//store/user/spandey/step2/PGun_step2_DIGI_902_Apr_14_FULL/CRAB_UserFiles/crab_PGun_step2_DIGI_902_Apr_14_FULL/170414_125708/0000/step2_2.root'),
    #fileNames = cms.untracked.vstring(),
    secondaryFileNames = cms.untracked.vstring(),
     duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
)

process.options = cms.untracked.PSet(
    allowUnscheduled = cms.untracked.bool(True)
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step3 nevts:10'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)


#Setup FWK for multithreaded
# process.options.numberOfThreads=cms.untracked.uint32(8)
# process.options.numberOfStreams=cms.untracked.uint32(0)

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '93X_upgrade2023_realistic_v2', '')



process.pfChargedHadronAnalyzer = cms.EDAnalyzer(
    "PFHgcalAnalyzer",
    PFCandidates = cms.InputTag("particleFlow"),
    # SimTracks = cms.InputTag("g4SimHits"),
    GenParticles = cms.InputTag("genParticles"),
    ptMin = cms.double(1.),                     # Minimum pt                                                                         
    pMin = cms.double(1.),                      # Minimum p                                                                          
    nPixMin = cms.int32(2),                     # Nb of pixel hits                                                                   
    nHitMin = cms.vint32(14,17,20,17,10),       # Nb of track hits                                                                   
    nEtaMin = cms.vdouble(1.4,1.6,2.0,2.4,2.6), # in these eta ranges                                                                
    calMin = cms.double(0.5),                   # Minimum hcal energy                                                               
    ecalMax = cms.double(1E9),                  # Maximum ecal energy                                                                                                                                                                    
    rootOutputFile = cms.string("step4.root"),# the root tree                                                                                                                                                
)


# process.load("RecoParticleFlow.PFProducer.particleFlowSimParticle_cfi")
#process.load("RecoParticleFlow.Configuration.HepMCCopy_cfi")                                                                        


# process.genReReco = cms.Sequence(process.particleFlowSimParticle)


process.bla = cms.EndPath(process.pfChargedHadronAnalyzer)
# process.blo = cms.EndPath(process.genReReco)

process.schedule = cms.Schedule(process.bla)#,process.RECOSIMoutput_step)process.blo,


