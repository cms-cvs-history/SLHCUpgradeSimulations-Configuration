import FWCore.ParameterSet.Config as cms

process = cms.Process("DigiReco")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Services_cff")

#process.load('Configuration.StandardSequences.MixingNoPileUp_cff')
process.load('SLHCUpgradeSimulations.Configuration.mixLowLumPU_HCal_P1Trk_cff')
#process.load('Configuration.StandardSequences.GeometryExtended_cff')
process.load('SLHCUpgradeSimulations.Geometry.PhaseI_cmsSimIdealGeometryXML_R39F16_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load("Configuration.StandardSequences.VtxSmearedBetafuncNominalCollision_cff")

process.load("Configuration.StandardSequences.SimulationRandomNumberGeneratorSeeds_cff")

process.load("Configuration.StandardSequences.Sim_cff")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
#process.load("Configuration.StandardSequences.Digi_cff")
process.load('SLHCUpgradeSimulations.Geometry.Digi_Phase1_cff')

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.load("Configuration.StandardSequences.FakeConditions_cff")
process.GlobalTag.globaltag = 'DESIGN_36_V10::All'

# use hardcoded values
process.es_hardcode.toGet.extend(['Gains', 'Pedestals', 'PedestalWidths', 'QIEData', 'ElectronicsMap','ChannelQuality','RespCorrs','ZSThresholds','LutMetadata','L1TriggerObjects','TimeCorrs','PFCorrs','LUTCorrs'])
process.es_hardcode.H2Mode = cms.untracked.bool(False)
process.es_hardcode.SLHCMode = cms.untracked.bool(True)
process.es_prefer_hcalHardcode = cms.ESPrefer("HcalHardcodeCalibrations", "es_hardcode")

# Event input
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
	'/store/relval/CMSSW_3_6_3_SLHC3_patch1/RelValFourMuPt_1_50/GEN-SIM/DESIGN_36_V10_Gauss_special-v1/0013/7A48047C-9351-E011-933F-00261894393F.root'
	)
)   

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)

# Event output
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load("Configuration.EventContent.EventContent_cff")
#myOutputCommands = cms.untracked.vstring('keep *')
myOutputCommands = cms.untracked.vstring()
myOutputCommands.extend(process.FEVTSIMEventContent.outputCommands)
myOutputCommands.extend([
    'keep *_hcalDigis_*_*', 'keep *_simHcalUnsuppressedDigis_*_*',
    'keep *_towerMakerWithHO_*_*'
    ])

process.output = cms.OutputModule("PoolOutputModule",
    outputCommands = myOutputCommands,
    fileName = cms.untracked.string('test_PU.root'),
#    SelectEvents = cms.untracked.PSet(
#        SelectEvents = cms.vstring('p0')
#    )
)

doRelabeling = True

#turn on test numbering
process.g4SimHits.Physics.type = 'SimG4Core/Physics/QGSP_BERT_EMV'
process.g4SimHits.HCalSD.TestNumberingScheme = doRelabeling

#turn off zero suppression hopefully ?
process.simHcalDigis.HBlevel = -1000
process.simHcalDigis.HElevel = -1000
process.simHcalDigis.HOlevel = -1000

#turn on SiPMs in HO
process.hcalSimParameters.ho.siPMCode = 1
process.hcalSimParameters.ho.pixels = cms.int32(2500)
process.hcalSimParameters.ho.photoelectronsToAnalog = 3.0

#turn on SiPMs in HB/HE
process.hcalSimParameters.hb.siPMCells = [1]
process.hcalSimParameters.hb.pixels = cms.int32(4500*4)
process.hcalSimParameters.hb.photoelectronsToAnalog = 10.0
process.hcalSimParameters.he.pixels = cms.int32(4500*4)
process.hcalSimParameters.he.photoelectronsToAnalog = 10.0

#turn on SLHC topology
process.HcalTopologyIdealEP.SLHCMode = cms.untracked.bool(True)

#turn on hit relabeling and set depth segmentation
process.simHcalUnsuppressedDigis.RelabelHits = cms.untracked.bool(doRelabeling)
process.simHcalUnsuppressedDigis.RelabelRules = cms.untracked.PSet(
    # Eta1 = cms.untracked.vint32(1,1,2,2,2,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4),
    Eta1 = cms.untracked.vint32(1,2,2,2,2,3,3,3,3,4,4,4,4,4,4,4,4,5,5),
    #Eta17 = cms.untracked.vint32(1,1,1,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4)
    Eta17 = cms.untracked.vint32(1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,5,5,5,5,5)
    )

process.pgen.remove(process.genJetMET)
process.pdigi.remove(process.simHcalTriggerPrimitiveDigis)

process.mix.input.nbPileupEvents = cms.PSet(
  averageNumber = cms.double(25.0)
)
process.Timing =  cms.Service("Timing")
process.simSiPixelDigis.AddPixelInefficiency = 20

# ############################################################### -->
# Keep the RecoValidation steps for the Phase1 tracker together for now...
process.load('Configuration.StandardSequences.Reconstruction_cff')

process.load("SLHCUpgradeSimulations.Geometry.fakeConditions_Phase1_cff")
process.load("SLHCUpgradeSimulations.Geometry.fakeConditions_Phase1_R39F16_cff")
process.load("SLHCUpgradeSimulations.Geometry.recoFromSimDigis_cff")
process.load("SLHCUpgradeSimulations.Geometry.upgradeTracking_phase1_cff")

process.ctfWithMaterialTracks.TTRHBuilder = 'WithTrackAngle'
process.PixelCPEGenericESProducer.UseErrorsFromTemplates = cms.bool(False)
process.PixelCPEGenericESProducer.TruncatePixelCharge = cms.bool(False)
process.PixelCPEGenericESProducer.LoadTemplatesFromDB = cms.bool(False)
process.PixelCPEGenericESProducer.Upgrade = cms.bool(True)
process.PixelCPEGenericESProducer.SmallPitch = False
process.PixelCPEGenericESProducer.IrradiationBiasCorrection = False
process.PixelCPEGenericESProducer.DoCosmics = False

## CPE for other steps
process.siPixelRecHits.CPE = cms.string('PixelCPEGeneric')
process.newPixelRecHits.CPE = cms.string('PixelCPEGeneric')
process.secPixelRecHits.CPE = cms.string('PixelCPEGeneric')
process.thPixelRecHits.CPE = cms.string('PixelCPEGeneric')
process.preFilterZeroStepTracks.TTRHBuilder = cms.string('WithTrackAngle')
process.preFilterStepOneTracks.TTRHBuilder = cms.string('WithTrackAngle')
process.secWithMaterialTracks.TTRHBuilder = cms.string('WithTrackAngle')
process.thWithMaterialTracks.TTRHBuilder = cms.string('WithTrackAngle')
process.fourthWithMaterialTracks.TTRHBuilder = cms.string('WithTrackAngle')
process.fifthWithMaterialTracks.TTRHBuilder = cms.string('WithTrackAngle')

# Need these lines to stop some errors about missing siStripDigis collections.
# should add them to fakeConditions_Phase1_cff
process.MeasurementTracker.UseStripStripQualityDB      = cms.bool(False)
process.MeasurementTracker.UsePixelModuleQualityDB     = cms.bool(False)
process.MeasurementTracker.UsePixelROCQualityDB        = cms.bool(False)
process.newMeasurementTracker.inactiveStripDetectorLabels = cms.VInputTag()
process.newMeasurementTracker.UseStripModuleQualityDB     = cms.bool(False)
process.newMeasurementTracker.UseStripAPVFiberQualityDB   = cms.bool(False)
process.newMeasurementTracker.UseStripStripQualityDB      = cms.bool(False)
process.newMeasurementTracker.UsePixelModuleQualityDB     = cms.bool(False)
process.newMeasurementTracker.UsePixelROCQualityDB        = cms.bool(False)
process.secMeasurementTracker.inactiveStripDetectorLabels = cms.VInputTag()
process.secMeasurementTracker.UseStripModuleQualityDB     = cms.bool(False)
process.secMeasurementTracker.UseStripAPVFiberQualityDB   = cms.bool(False)
process.secMeasurementTracker.UseStripStripQualityDB      = cms.bool(False)
process.secMeasurementTracker.UsePixelModuleQualityDB     = cms.bool(False)
process.secMeasurementTracker.UsePixelROCQualityDB        = cms.bool(False)
process.thMeasurementTracker.inactiveStripDetectorLabels = cms.VInputTag()
process.thMeasurementTracker.UseStripModuleQualityDB     = cms.bool(False)
process.thMeasurementTracker.UseStripAPVFiberQualityDB   = cms.bool(False)
process.thMeasurementTracker.UseStripStripQualityDB      = cms.bool(False)
process.thMeasurementTracker.UsePixelModuleQualityDB     = cms.bool(False)
process.thMeasurementTracker.UsePixelROCQualityDB        = cms.bool(False)
process.fourthMeasurementTracker.inactiveStripDetectorLabels = cms.VInputTag()
process.fifthMeasurementTracker.inactiveStripDetectorLabels = cms.VInputTag()

# Now Validation and other user functions #######################
process.load("Validation.RecoTrack.cutsTPEffic_cfi")
process.load("Validation.RecoTrack.cutsTPFake_cfi")

process.load("SimTracker.TrackAssociation.TrackAssociatorByChi2_cfi")
process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi")
process.load('SimTracker.TrackAssociation.quickTrackAssociatorByHits_cfi')

process.load('Configuration.StandardSequences.Validation_cff')
### look look at OOTB generalTracks and high purity collections
### for high purity also look at 6 and 8 hit requirements
### some definitions in Validation/RecoTrack/python/TrackValidation_cff.py

import PhysicsTools.RecoAlgos.recoTrackSelector_cfi

process.cutsRecoTracksHpw6hits = PhysicsTools.RecoAlgos.recoTrackSelector_cfi.recoTrackSelector.clone()
process.cutsRecoTracksHpw6hits.quality=cms.vstring("highPurity")
process.cutsRecoTracksHpw6hits.minHit=cms.int32(6)

process.cutsRecoTracksHpw8hits = PhysicsTools.RecoAlgos.recoTrackSelector_cfi.recoTrackSelector.clone()
process.cutsRecoTracksHpw8hits.quality=cms.vstring("highPurity")
process.cutsRecoTracksHpw8hits.minHit=cms.int32(8)

process.trackValidator.label=cms.VInputTag(cms.InputTag("generalTracks"),
                                           cms.InputTag("cutsRecoTracksHp")
#                                           cms.InputTag("cutsRecoTracksHpw6hits"),
#                                           cms.InputTag("cutsRecoTracksHpw8hits"),
#                                           cms.InputTag("cutsRecoTracksZeroHp"),
#                                           cms.InputTag("cutsRecoTracksFirstHp")
                                           )
#process.trackValidator.associators = ['TrackAssociatorByHits']
process.trackValidator.associators = cms.vstring('quickTrackAssociatorByHits')
process.trackValidator.UseAssociators = True
process.trackValidator.nint = cms.int32(20)
process.trackValidator.nintpT = cms.int32(100)
process.trackValidator.maxpT = cms.double(200.0)
process.trackValidator.outputFile = "validfullP1.root"
process.trackValidator.skipHistoFit = False

process.slhcTracksValidation = cms.Sequence(process.cutsRecoTracksHp*
                                 process.cutsRecoTracksHpw6hits*
                                 process.cutsRecoTracksHpw8hits*
                                 process.cutsRecoTracksZeroHp*
                                 process.cutsRecoTracksFirstHp*
#                                 process.cutsRecoTracksSecondHp*
#                                 process.cutsRecoTracksThirdHp*
                                 process.trackValidator)

process.ReadLocalMeasurement = cms.EDAnalyzer("StdHitNtuplizer",
   src = cms.InputTag("siPixelRecHits"),
   stereoRecHits = cms.InputTag("siStripMatchedRecHits","stereoRecHit"),
   rphiRecHits = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
   matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
   ### if using simple (non-iterative) or old (as in 1_8_4) tracking
   trackProducer = cms.InputTag("generalTracks"),
   OutputFile = cms.string("stdgrechitfullph1g_ntuple.root"),
   ### for using track hit association
   associatePixel = cms.bool(True),
   associateStrip = cms.bool(False),
   associateRecoTracks = cms.bool(False),
   ROUList = cms.vstring('g4SimHitsTrackerHitsPixelBarrelLowTof',
                         'g4SimHitsTrackerHitsPixelBarrelHighTof',
                         'g4SimHitsTrackerHitsPixelEndcapLowTof',
                         'g4SimHitsTrackerHitsPixelEndcapHighTof')
)
process.anal = cms.EDAnalyzer("EventContentAnalyzer")

## need this at the end as the validation config redefines random seed with just mix
process.load("IOMC.RandomEngine.IOMC_cff")
# ############################################################### <--

process.digitisation_step       = cms.Path(process.pdigi)
process.reconstruction_step     = cms.Path(process.trackerlocalreco*
                                                process.offlineBeamSpot+
                                                process.recopixelvertexing*process.ckftracks_wodEdXandSteps2345)
#process.reconstruction_step     = cms.Path(process.reconstruction)
process.debug_step              = cms.Path(process.anal)
process.validation_step         = cms.Path(process.cutsTPEffic*
                                                process.cutsTPFake*
                                                process.slhcTracksValidation)
process.user_step               = cms.Path(process.ReadLocalMeasurement)
process.endjob_step             = cms.Path(process.endOfProcess)
process.outpath = cms.EndPath(process.output)
process.schedule = cms.Schedule(process.digitisation_step,
				process.reconstruction_step,
				#process.debug_step,
				process.validation_step,
				process.user_step,
				process.endjob_step,
				process.outpath)
