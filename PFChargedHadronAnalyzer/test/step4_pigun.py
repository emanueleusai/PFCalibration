# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step3 --conditions auto:phase1_2017_realistic -n 10 --era Run2_2017 --eventcontent RECOSIM,MINIAODSIM,DQM --runUnscheduled -s RAW2DIGI,L1Reco,RECO,EI,PAT,VALIDATION:@standardValidation+@miniAODValidation,DQM:@standardDQM+@miniAODDQM --datatier GEN-SIM-RECO,MINIAODSIM,DQMIO --geometry DB:Extended --conditions 90X_upgrade2017_realistic_v15 --filein file:step2.root --fileout file:step3.root
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

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
    fileNames = cms.untracked.vstring(
    'file:step3.root',

    ),

    #fileNames = cms.untracked.vstring('file:E20ADE15-A5AE-E711-9434-0023AEEEB55F.root'),
    #fileNames = cms.untracked.vstring('root://se01.indiacms.res.in//store/user/spandey/step2/PGun_step2_DIGI_902_Apr_14_FULL/CRAB_UserFiles/crab_PGun_step2_DIGI_902_Apr_14_FULL/170414_125708/0000/step2_2.root'),
    #fileNames = cms.untracked.vstring(),
    secondaryFileNames = cms.untracked.vstring()
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
process.options.numberOfThreads=cms.untracked.uint32(8)
process.options.numberOfStreams=cms.untracked.uint32(0)

# Output definition

# process.RECOSIMoutput = cms.OutputModule("PoolOutputModule",
#     dataset = cms.untracked.PSet(
#         dataTier = cms.untracked.string('GEN-SIM-RECO'),
#         filterName = cms.untracked.string('')
#     ),
#     eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
#     fileName = cms.untracked.string('file:step3cms.root'),
#     outputCommands = cms.untracked.vstring('keep *'),
#     splitLevel = cms.untracked.int32(0)
# )

# process.MINIAODSIMoutput = cms.OutputModule("PoolOutputModule",
#     compressionAlgorithm = cms.untracked.string('LZMA'),
#     compressionLevel = cms.untracked.int32(4),
#     dataset = cms.untracked.PSet(
#         dataTier = cms.untracked.string('MINIAODSIM'),
#         filterName = cms.untracked.string('')
#     ),
#     dropMetaData = cms.untracked.string('ALL'),
#     eventAutoFlushCompressedSize = cms.untracked.int32(15728640),
#     fastCloning = cms.untracked.bool(False),
#     fileName = cms.untracked.string('file:step3_inMINIAODSIM.root'),
#     outputCommands = process.MINIAODSIMEventContent.outputCommands,
#     overrideInputFileSplitLevels = cms.untracked.bool(True)
# )

# process.DQMoutput = cms.OutputModule("DQMRootOutputModule",
#     dataset = cms.untracked.PSet(
#         dataTier = cms.untracked.string('DQMIO'),
#         filterName = cms.untracked.string('')
#     ),
#     fileName = cms.untracked.string('file:step3_inDQM.root'),
#     outputCommands = process.DQMEventContent.outputCommands,
#     splitLevel = cms.untracked.int32(0)
# )

# Additional output definition

# # Other statements
# process.mix.playback = True
# process.mix.digitizers = cms.PSet()
# for a in process.aliases: delattr(process, a)
# process.RandomNumberGeneratorService.restoreStateLabel=cms.untracked.string("randomEngineStateProducer")
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2017_realistic', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, '92X_upgrade2017_realistic_forJetMET_HCALScaleStudies', '')
process.GlobalTag = GlobalTag(process.GlobalTag, '93X_upgrade2023_realistic_v2', '')



process.pfChargedHadronAnalyzer = cms.EDAnalyzer(
    "PFChargedHadronAnalyzer",
    #"PFHgcalAnalyzer",
    PFCandidates = cms.InputTag("particleFlow"),
    PFSimParticles = cms.InputTag("particleFlowSimParticle"),
    EcalPFClusters = cms.InputTag("particleFlowClusterECAL"),
    HcalPFClusters = cms.InputTag("particleFlowClusterHCAL"),
    ptMin = cms.double(1.),                     # Minimum pt                                                                         
    pMin = cms.double(1.),                      # Minimum p                                                                          
    nPixMin = cms.int32(2),                     # Nb of pixel hits                                                                   
    nHitMin = cms.vint32(14,17,20,17,10),       # Nb of track hits                                                                   
    nEtaMin = cms.vdouble(1.4,1.6,2.0,2.4,2.6), # in these eta ranges                                                                
    hcalMin = cms.double(0.5),                   # Minimum hcal energy                                                               
    ecalMax = cms.double(1E9),                  # Maximum ecal energy                                                                
    verbose = cms.untracked.bool(True),         # not used.                                                                          
    #rootOutputFile = cms.string("PGun__2_200GeV__81X_upgrade2017_realistic_v22.root"),# the root tree                               
    rootOutputFile = cms.string("step4.root"),# the root tree                                                       
#    IsMinBias = cms.untracked.bool(False)                                                                                           
)




process.load("RecoParticleFlow.PFProducer.particleFlowSimParticle_cfi")
#process.load("RecoParticleFlow.Configuration.HepMCCopy_cfi")                                                                        

process.particleFlowSimParticle.ParticleFilter = cms.PSet(
        # Allow *ALL* protons with energy > protonEMin  
        verbose = cms.untracked.bool(True),                                                                    
        protonEMin = cms.double(5000.0),
        # Particles must have abs(eta) < etaMax (if close enough to 0,0,0)                                                           
        etaMax = cms.double(5.3),
        # Charged particles with pT < pTMin (GeV/c) are not simulated                                                                
        chargedPtMin = cms.double(0.0),
        # Particles must have energy greater than EMin [GeV]                                                                         
        EMin = cms.double(0.0),
        ecalRecHitsEB = cms.InputTag('ecalRecHit','EcalRecHitsEB'),
        ecalRecHitsEE = cms.InputTag('ecalRecHit','EcalRecHitsEE'),
        )



# process.load("RecoParticleFlow.PFProducer.simPFProducer_cfi")
# process.load("RecoParticleFlow.PFTracking.hgcalTrackCollection_cfi")
# process.load("RecoParticleFlow.PFTracking.pfTrack_cfi")
# process.load("RecoMuon.MuonIdentification.muons1stStep_cfi")

process.genReReco = cms.Sequence(#process.generator+                                                                                 
                                 #process.genParticles+                                                                              
                                 #process.genJetParticles+                                                                           
                                 #process.recoGenJets+                                                                               
                                 #process.genMETParticles+                                                                           
                                 #process.recoGenMET+
                                 # process.muons1stStep+
                                 # process.pfTrack+
                                 # process.hgcalTrackCollection+
                                 # process.simPFProducer+                                                                                
                                 process.particleFlowSimParticle)








# Path and EndPath definitions
# process.raw2digi_step = cms.Path(process.RawToDigi) 
# process.L1Reco_step = cms.Path(process.L1Reco)
# process.reconstruction_step = cms.Path(process.reconstruction)
# process.eventinterpretaion_step = cms.Path(process.EIsequence)
#process.endjob_step = cms.EndPath(process.endOfProcess)
process.bla = cms.EndPath(process.pfChargedHadronAnalyzer)
process.blo = cms.EndPath(process.genReReco)


# process.Flag_trackingFailureFilter = cms.Path(process.goodVertices+process.trackingFailureFilter)
# process.Flag_goodVertices = cms.Path(process.primaryVertexFilter)
# process.Flag_CSCTightHaloFilter = cms.Path(process.CSCTightHaloFilter)
# process.Flag_trkPOGFilters = cms.Path(process.trkPOGFilters)
# process.Flag_trkPOG_logErrorTooManyClusters = cms.Path(~process.logErrorTooManyClusters)
# process.Flag_EcalDeadCellTriggerPrimitiveFilter = cms.Path(process.EcalDeadCellTriggerPrimitiveFilter)
# process.Flag_ecalLaserCorrFilter = cms.Path(process.ecalLaserCorrFilter)
# process.Flag_globalSuperTightHalo2016Filter = cms.Path(process.globalSuperTightHalo2016Filter)
# process.Flag_eeBadScFilter = cms.Path(process.eeBadScFilter)
# process.Flag_METFilters = cms.Path(process.metFilters)
# process.Flag_chargedHadronTrackResolutionFilter = cms.Path(process.chargedHadronTrackResolutionFilter)
# process.Flag_globalTightHalo2016Filter = cms.Path(process.globalTightHalo2016Filter)
# process.Flag_CSCTightHaloTrkMuUnvetoFilter = cms.Path(process.CSCTightHaloTrkMuUnvetoFilter)
# process.Flag_HBHENoiseIsoFilter = cms.Path(process.HBHENoiseFilterResultProducer+process.HBHENoiseIsoFilter)
# process.Flag_hcalLaserEventFilter = cms.Path(process.hcalLaserEventFilter)
# process.Flag_HBHENoiseFilter = cms.Path(process.HBHENoiseFilterResultProducer+process.HBHENoiseFilter)
# process.Flag_trkPOG_toomanystripclus53X = cms.Path(~process.toomanystripclus53X)
# process.Flag_EcalDeadCellBoundaryEnergyFilter = cms.Path(process.EcalDeadCellBoundaryEnergyFilter)
# process.Flag_trkPOG_manystripclus53X = cms.Path(~process.manystripclus53X)
# process.Flag_HcalStripHaloFilter = cms.Path(process.HcalStripHaloFilter)
# process.Flag_muonBadTrackFilter = cms.Path(process.muonBadTrackFilter)
# process.Flag_CSCTightHalo2015Filter = cms.Path(process.CSCTightHalo2015Filter)
# process.prevalidation_step = cms.Path(process.prevalidation)
# process.prevalidation_step1 = cms.Path(process.prevalidationMiniAOD)
# process.validation_step = cms.EndPath(process.validation)
# process.validation_step1 = cms.EndPath(process.validationMiniAOD)
# process.dqmoffline_step = cms.EndPath(process.DQMOffline)
# process.dqmoffline_1_step = cms.EndPath(process.DQMOfflineMiniAOD)
# process.dqmofflineOnPAT_step = cms.EndPath(process.PostDQMOffline)
# process.dqmofflineOnPAT_1_step = cms.EndPath(process.PostDQMOfflineMiniAOD)

# process.RECOSIMoutput_step = cms.EndPath(process.RECOSIMoutput)

# process.MINIAODSIMoutput_step = cms.EndPath(process.MINIAODSIMoutput)
# process.DQMoutput_step = cms.EndPath(process.DQMoutput)

# Schedule definition
#process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.eventinterpretaion_step,process.Flag_HBHENoiseFilter,process.Flag_HBHENoiseIsoFilter,process.Flag_CSCTightHaloFilter,process.Flag_CSCTightHaloTrkMuUnvetoFilter,process.Flag_CSCTightHalo2015Filter,process.Flag_globalTightHalo2016Filter,process.Flag_globalSuperTightHalo2016Filter,process.Flag_HcalStripHaloFilter,process.Flag_hcalLaserEventFilter,process.Flag_EcalDeadCellTriggerPrimitiveFilter,process.Flag_EcalDeadCellBoundaryEnergyFilter,process.Flag_goodVertices,process.Flag_eeBadScFilter,process.Flag_ecalLaserCorrFilter,process.Flag_trkPOGFilters,process.Flag_chargedHadronTrackResolutionFilter,process.Flag_muonBadTrackFilter,process.Flag_trkPOG_manystripclus53X,process.Flag_trkPOG_toomanystripclus53X,process.Flag_trkPOG_logErrorTooManyClusters,process.Flag_METFilters,process.prevalidation_step,process.prevalidation_step1,process.validation_step,process.validation_step1,process.dqmoffline_step,process.dqmoffline_1_step,process.dqmofflineOnPAT_step,process.dqmofflineOnPAT_1_step,process.RECOSIMoutput_step,process.MINIAODSIMoutput_step,process.DQMoutput_step)



#process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.eventinterpretaion_step,process.endjob_step,process.blo,process.bla)

#process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.eventinterpretaion_step,process.blo,process.bla)
process.schedule = cms.Schedule(process.blo,process.bla)#,process.RECOSIMoutput_step)

# customisation of the process.

# Automatic addition of the customisation function from SimGeneral.MixingModule.fullMixCustomize_cff
#from SimGeneral.MixingModule.fullMixCustomize_cff import setCrossingFrameOn 

#call to customisation function setCrossingFrameOn imported from SimGeneral.MixingModule.fullMixCustomize_cff
#process = setCrossingFrameOn(process)

# End of customisation functions
#do not add changes to your config after this point (unless you know what you are doing)
# from FWCore.ParameterSet.Utilities import convertToUnscheduled
# process=convertToUnscheduled(process)
# process.load('Configuration.StandardSequences.PATMC_cff')
# from FWCore.ParameterSet.Utilities import cleanUnscheduled
# process=cleanUnscheduled(process)

# customisation of the process.

# # Automatic addition of the customisation function from PhysicsTools.PatAlgos.slimming.miniAOD_tools
# from PhysicsTools.PatAlgos.slimming.miniAOD_tools import miniAOD_customizeAllMC 

# #call to customisation function miniAOD_customizeAllMC imported from PhysicsTools.PatAlgos.slimming.miniAOD_tools
# process = miniAOD_customizeAllMC(process)

# End of customisation functions

# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
# from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
# process = customiseEarlyDelete(process)
# End adding early deletion
