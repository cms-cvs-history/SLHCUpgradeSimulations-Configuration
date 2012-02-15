# Generate a Gen-Sim Digi file with cmsDriver.py with the following options
# cmsDriver.py SLHCUpgradeSimulations/Configuration/python/FourMuPt_1_200_cfi -n 5 -s GEN,SIM,DIGI --conditions DESIGN42_V17::All --eventcontent FEVTDEBUG --beamspot Gauss --slhc Phase1_R30F12_HCal --datatier GEN-SIM-DIGI --no_exec
# Then run Reco over it with this
# Passed down a couple times with changes, but it originally came from /afs/cern.ch/user/r/rpw/public/SLHC
import FWCore.ParameterSet.Config as cms

process = cms.Process("Reco")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Services_cff")

process.load('SimGeneral.MixingModule.mixNoPU_cfi')
#process.load('SLHCUpgradeSimulations.Geometry.mixLowLumPU_Phase1_R30F12_HCal_cff')
process.load('SLHCUpgradeSimulations.Geometry.Phase1_R30F12_HCal_cmsSimIdealGeometryXML_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load('SLHCUpgradeSimulations.Geometry.Digi_Phase1_R30F12_HCal_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'DESIGN42_V17::All'

# use hardcoded values
process.es_hardcode.toGet.extend(['Gains', 'Pedestals', 'PedestalWidths', 'QIEData', 'ElectronicsMap','ChannelQuality','RespCorrs','ZSThresholds','L1TriggerObjects','TimeCorrs','PFCorrs','LUTCorrs', "RecoParams"])
process.es_hardcode.H2Mode = cms.untracked.bool(False)
process.es_hardcode.SLHCMode = cms.untracked.bool(True)
# Use text file for the LutMetadata
process.es_ascii = cms.ESSource("HcalTextCalibrations",
                                input = cms.VPSet (
        cms.PSet (
            object = cms.string ('LutMetadata'),
            file = cms.FileInPath('CondFormats/HcalObjects/data/HcalLutMetadataLSB50.txt')
            )
        )
)

# ESPrefers
process.es_prefer_hcalAscii    = cms.ESPrefer("HcalTextCalibrations"    , "es_ascii")
process.es_prefer_hcalHardcode = cms.ESPrefer("HcalHardcodeCalibrations", "es_hardcode")

process.source = cms.Source("PoolSource",
    firstEvent = cms.untracked.uint32(1),
    fileNames = cms.untracked.vstring(
      'file:FourMuPt_1_200_cfi_GEN_SIM_DIGI.root'
  )
)

# Event output
process.load("Configuration.EventContent.EventContent_cff")
process.output = cms.OutputModule("PoolOutputModule",
    #outputCommands = process.RECOEventContent.outputCommands,
    outputCommands = cms.untracked.vstring('keep *'),
    fileName = cms.untracked.string('file:FourMuPt_1_200_cfi_GEN_SIM_DIGI_RECO.root'),
)


doRelabeling = False


#turn on SLHC topology

process.HcalTopologyIdealEP.SLHCMode = cms.untracked.bool(True)
process.hcalRecHitDump = cms.EDAnalyzer("HcalRecHitDump")
process.caloTowerDump = cms.EDAnalyzer("CaloTowersDump")
process.eventContent = cms.EDAnalyzer("EventContentAnalyzer")
process.load("SimCalorimetry.HcalTrigPrimProducers.hcalupgradetpdigi_cff")

#turn on SLHC topology
process.HcalTopologyIdealEP.SLHCMode = cms.untracked.bool(True)
#process.simHcalUnsuppressedDigis.RelabelRules = cms.untracked.PSet(
#    Eta1 = cms.untracked.vint32(1,2,2,2,2,3,3,3,3,4,4,4,4,4,4,4,4,5,5),
#    Eta17 = cms.untracked.vint32(1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,5,5,5,5,5)
#    )


process.load("RecoLocalCalo.HcalRecProducers.HcalUpgradeReconstructor_cff")
process.load("RecoJets/Configuration/CaloTowersRec_cff")
process.load("RecoLocalCalo.HcalRecAlgos.hcalRecAlgoESProd_cfi")
process.load("RecoLocalCalo.Configuration.RecoLocalCalo_cff")
process.ecalGlobalUncalibRecHit.EBdigiCollection = cms.InputTag("simEcalDigis","ebDigis")
process.ecalGlobalUncalibRecHit.EEdigiCollection = cms.InputTag("simEcalDigis","eeDigis")

#process.load("RecoLocalCalo.CaloTowersCreator.calotowermaker_cfi")
#process.load("Geometry.CaloEventSetup.CaloTowerConstituents_cfi")
#process.mix.input.nbPileupEvents.averageNumber = 0.
process.hcalupgradereco.digiLabel = "simHcalUnsuppressedDigis"

process.p0 = cms.Path(process.ecalLocalRecoSequence_nopreshower+process.hcalupgradereco+process.hcalRecHitDump+process.calotowermaker+process.caloTowerDump)
process.outpath = cms.EndPath(process.output)
process.schedule = cms.Schedule(process.p0,process.outpath)

#
#----- Summary -----
#-------------------
process.options = cms.untracked.PSet(
           wantSummary = cms.untracked.bool(True)
           )
