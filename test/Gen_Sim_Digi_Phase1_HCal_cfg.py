import FWCore.ParameterSet.Config as cms

process = cms.Process("GenSimDigi")

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

# Event output
process.load("Configuration.EventContent.EventContent_cff")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(3)
)

ceniEta = 1
ceniPhi = 1
offset = 0

cenEta = (ceniEta - 0.5)*0.087
if (ceniPhi < 0):
    cenEta = (ceniEta + 0.5)*0.087
cenPhi = (ceniPhi - 0.5 + offset)*0.087
energy = 100.

process.source = cms.Source("EmptySource")

process.generator = cms.EDProducer(
    "FlatRandomEGunProducer",
    PGunParameters = cms.PSet(
        PartID = cms.vint32(13),
        MaxEta = cms.double(cenEta),
        MaxPhi = cms.double(cenPhi),
        MinEta = cms.double(cenEta),
        MinPhi = cms.double(cenPhi),
        MinE = cms.double(energy),
        MaxE = cms.double(energy)
        ),
    Verbosity = cms.untracked.int32(0),
    psethack = cms.string('single particle E ' + str(energy)),
    AddAntiParticle = cms.bool(False),
    firstRun = cms.untracked.uint32(1)
    )

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
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('p0')
    )
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

process.p0 = cms.Path(process.pgen)
process.p1 = cms.Path(process.psim)
#process.p2 = cms.Path(process.pdigi)
process.p2 = cms.Path(process.randomEngineStateProducer+process.mix+process.simHcalUnsuppressedDigis+process.simHcalDigis)
process.outpath = cms.EndPath(process.output)
process.schedule = cms.Schedule(process.p0,
                                process.p1,
                                process.p2,
                                process.outpath)

# special treatment in case of production filter sequence  
for path in process.paths: 
    getattr(process,path)._seq = process.generator*getattr(process,path)._seq
