import FWCore.ParameterSet.Config as cms

process = cms.Process("GenSimDigi")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Services_cff")

process.load('SimGeneral.MixingModule.mixNoPU_cfi')
#process.load('Configuration.StandardSequences.GeometryExtended_cff')
process.load('SLHCUpgradeSimulations.Geometry.Phase1_R39F16_cmsSimIdealGeometryXML_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedGauss_cfi')

process.load("Configuration.StandardSequences.SimulationRandomNumberGeneratorSeeds_cff")

process.load("Configuration.StandardSequences.Sim_cff")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.load("Configuration.StandardSequences.FakeConditions_cff")
process.GlobalTag.globaltag = 'DESIGN42_V11::All'
process.load("SLHCUpgradeSimulations.Geometry.fakeConditions_Phase1_R39F16_cff")

# use hardcoded values
process.es_hardcode.toGet.extend(['Gains', 'Pedestals', 'PedestalWidths', 'QIEData', 'ElectronicsMap','ChannelQuality','RespCorrs','ZSThresholds','LutMetadata','L1TriggerObjects','TimeCorrs','PFCorrs','LUTCorrs'])
process.es_hardcode.H2Mode = cms.untracked.bool(False)
process.es_hardcode.SLHCMode = cms.untracked.bool(True)
process.es_prefer_hcalHardcode = cms.ESPrefer("HcalHardcodeCalibrations", "es_hardcode")

# Event output
process.load("Configuration.EventContent.EventContent_cff")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(25)
)

process.source = cms.Source("EmptySource")

process.generator = cms.EDFilter("Pythia6GeneratorFilter",
    pythiaPylistVerbosity = cms.untracked.int32(0),
    filterEfficiency = cms.untracked.double(1.0),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    comEnergy = cms.double(14000.0),
    maxEventsToPrint = cms.untracked.int32(0),
    PythiaParameters = cms.PSet(
        pythiaUESettings = cms.vstring('MSTJ(11)=3     ! Choice of the fragmentation function',
            'MSTJ(22)=2     ! Decay those unstable particles', 
            'PARJ(71)=10 .  ! for which ctau  10 mm', 
            'MSTP(2)=1      ! which order running alphaS',
            'MSTP(33)=0     ! no K factors in hard cross sections',
            'MSTP(51)=10042 ! structure function chosen (external PDF CTEQ6L1)',
            'MSTP(52)=2     ! work with LHAPDF', 
            'MSTP(81)=1     ! multiple parton interactions 1 is Pythia default',
            'MSTP(82)=4     ! Defines the multi-parton model', 
            'MSTU(21)=1     ! Check on possible errors during program execution',
            'PARP(82)=1.8387   ! pt cutoff for multiparton interactions', 
            'PARP(89)=1960. ! sqrts for which PARP82 is set', 
            'PARP(83)=0.5   ! Multiple interactions: matter distrbn parameter',
            'PARP(84)=0.4   ! Multiple interactions: matter distribution parameter',
            'PARP(90)=0.16  ! Multiple interactions: rescaling power', 
            'PARP(67)=2.5    ! amount of initial-state radiation', 
            'PARP(85)=1.0  ! gluon prod. mechanism in MI', 
            'PARP(86)=1.0  ! gluon prod. mechanism in MI',
            'PARP(62)=1.25   ! ', 
            'PARP(64)=0.2    ! ',
            'MSTP(91)=1      !', 
            'PARP(91)=2.1   ! kt distribution',
            'PARP(93)=15.0  ! '),
        processParameters = cms.vstring('MSEL=0         ! User defined processes',
            'MSUB(11)=1     ! Min bias process',
            'MSUB(12)=1     ! Min bias process',
            'MSUB(13)=1     ! Min bias process',
            'MSUB(28)=1     ! Min bias process',
            'MSUB(53)=1     ! Min bias process',
            'MSUB(68)=1     ! Min bias process',
            'MSUB(92)=1     ! Min bias process, single diffractive',
            'MSUB(93)=1     ! Min bias process, single diffractive',
            'MSUB(94)=1     ! Min bias process, double diffractive',
            'MSUB(95)=1     ! Min bias process'),
        parameterSets = cms.vstring('pythiaUESettings',
            'processParameters')
    )
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
    fileName = cms.untracked.string('MinBias_Phase1_R39F16_HCal.root'),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('p0')
    )
)

doRelabeling = True

#turn on test numbering
process.g4SimHits.Physics.type = 'SimG4Core/Physics/QGSP_BERT_EMV'
process.g4SimHits.HCalSD.TestNumberingScheme = doRelabeling

#turn on SLHC topology
process.HcalTopologyIdealEP.SLHCMode = cms.untracked.bool(True)

process.pgen.remove(process.genJetMET)

process.p0 = cms.Path(process.pgen)
process.p1 = cms.Path(process.psim)
process.outpath = cms.EndPath(process.output)
process.schedule = cms.Schedule(process.p0,
                                process.outpath)

# special treatment in case of production filter sequence  
for path in process.paths: 
    getattr(process,path)._seq = process.generator*getattr(process,path)._seq
