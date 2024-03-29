import FWCore.ParameterSet.Config as cms
from SLHCUpgradeSimulations.Configuration.postLS1Customs import customisePostLS1
from SLHCUpgradeSimulations.Configuration.phase2TkCustomsBE import customise as customiseBE
from SLHCUpgradeSimulations.Configuration.phase2TkCustoms_LB_6PS import customise as customiseLB6PS
from SLHCUpgradeSimulations.Configuration.phase2TkCustoms_LB_4LPS_2L2S import customise as customiseLB4LPS_2L2S
from SLHCUpgradeSimulations.Configuration.phase2TkCustomsBE import l1EventContent as customise_ev_BE
from SLHCUpgradeSimulations.Configuration.phase2TkCustoms_LB_6PS import l1EventContent as customise_ev_LB6PS
from SLHCUpgradeSimulations.Configuration.phase2TkCustoms_LB_4LPS_2L2S import l1EventContent as customise_ev_LB4LPS_2L2S
from SLHCUpgradeSimulations.Configuration.customise_mixing import customise_NoCrossing
from SLHCUpgradeSimulations.Configuration.phase1TkCustoms import customise as customisePhase1Tk
from SLHCUpgradeSimulations.Configuration.HCalCustoms import customise_HcalPhase1, customise_HcalPhase0

import SLHCUpgradeSimulations.Configuration.aging as aging

def cust_phase2_BE(process):
    process=customisePostLS1(process)
    process=customiseBE(process)
    process=customise_ev_BE(process)
    return process

def cust_phase2_LB6PS(process):
    process=customisePostLS1(process)
    process=customiseLB6PS(process)
    process=customise_ev_LB6PS(process)
    return process

def cust_phase2_LB4LPS_2L2S(process):
    process=customisePostLS1(process)
    process=customiseLB4LPS_2L2S(process)
    process=customise_ev_LB4LPS_2L2S(process)
    return process

def cust_2017(process):
    process=customisePostLS1(process)
    process=customisePhase1Tk(process)
    process=customise_HcalPhase0(process)
    process=fixRPCConditions(process)
    return process

def cust_2019(process):
    process=customisePostLS1(process)
    process=customisePhase1Tk(process)
    process=customise_HcalPhase1(process)
    process=fixRPCConditions(process)
    return process

def noCrossing(process):
    process=customise_NoCrossing(process)
    return process


def fixRPCConditions(process):
    if not hasattr(process.GlobalTag,'toGet'):
        process.GlobalTag.toGet=cms.VPSet()
    process.GlobalTag.toGet.extend( cms.VPSet(
        cms.PSet(record = cms.string("RPCStripNoisesRcd"),
                 tag = cms.string("RPCStripNoise_upscope_mc_v2"),
                 connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_31X_RPC")
                 )
        )
                                    )
    return process

##### clone aging.py here 
def agePixel(process,lumi):
    process=process.agePixel(process,lumi)
    return process

def ageHcal(process,lumi):
    process=aging.ageHcal(process,lumi)
    return process

def ageEcal(process,lumi):
    process=aging.ageEcal(process,lumi)
    return process

def customise_aging_100(process):
    process=aging.customise_aging_100(process)
    return process

def customise_aging_200(process):
    process=aging.customise_aging_200(process)
    return process

def customise_aging_300(process):
    process=aging.customise_aging_300(process)
    return process

def customise_aging_500(process):
    process=aging.customise_aging_500(process)
    return process

def customise_aging_1000(process):
    process=aging.customise_aging_1000(process)
    return process

def customise_aging_3000(process):
    process=aging.customise_aging_3000(process)
    return process

def customise_aging_ecalonly_300(process):
    process=aging.customise_aging_ecalonly_300(process)
    return process

def customise_aging_ecalonly_1000(process):
    process=aging.customise_aging_ecalonly_1000(process)
    return process

def customise_aging_ecalonly_3000(process):
    process=aging.customise_aging_ecalonly_3000(process)
    return process

def customise_aging_newpixel_1000(process):
    process=aging.customise_aging_newpixel_1000(process)
    return process

def customise_aging_newpixel_3000(process):
    process=aging.customise_aging_newpixel_3000(process)
    return process

def ecal_complete_aging(process):
    process=aging.ecal_complete_aging(process)
    return process

def turn_off_HE_aging(process):
    process=aging.turn_off_HE_aging(process)
    return process

def turn_off_HF_aging(process):
    process=aging.turn_off_HF_aging(process)
    return process

def turn_off_Pixel_aging(process):
    process=aging.turn_off_Pixel_aging(process)
    return process

def turn_on_Pixel_aging_1000(process):
    process=aging.turn_on_Pixel_aging_1000(process)
    return process

def hf_complete_aging(process):
    process=aging.hf_complete_aging(process)
    return process
    
def ecal_complete_aging_300(process):
    process=aging.ecal_complete_aging_300(process)
    return process

def ecal_complete_aging_1000(process):
    process=aging.ecal_complete_aging_1000(process)
    return process

def ecal_complete_aging_3000(process):
    process=aging.ecal_complete_aging_3000(process)
    return process

def fixEcalConditions_150(process):
    if not hasattr(process.GlobalTag,'toGet'):
        process.GlobalTag.toGet=cms.VPSet()
    process.GlobalTag.toGet.extend( cms.VPSet(
        cms.PSet(record = cms.string("EcalSRSettingsRcd"),
                 tag = cms.string("EcalSRSettings_TL150_mc"),
                 connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_34X_ECAL")
                 )
        )
                                    )
    process.GlobalTag.toGet.extend( cms.VPSet(
        cms.PSet(record = cms.string("EcalTPGLutIdMapRcd"),
                 tag = cms.string("EcalTPGLutIdMap_beamv5_upgrade_mc"),
                 connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_34X_ECAL")
                 )
        )
                                    )
    return process

def fixEcalConditions_100(process):
    if not hasattr(process.GlobalTag,'toGet'):
        process.GlobalTag.toGet=cms.VPSet()
    process.GlobalTag.toGet.extend( cms.VPSet(
        cms.PSet(record = cms.string("EcalSRSettingsRcd"),
                 tag = cms.string("EcalSRSettings_TL100_mc"),
                 connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_34X_ECAL")
                 )
        )
                                    )
    process.GlobalTag.toGet.extend( cms.VPSet(
        cms.PSet(record = cms.string("EcalTPGLutIdMapRcd"),
                 tag = cms.string("EcalTPGLutIdMap_beamv5_upgrade_mc"),
                 connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_34X_ECAL")
                 )
        )
                                    )
    return process

def fixEcalConditions_200(process):
    if not hasattr(process.GlobalTag,'toGet'):
        process.GlobalTag.toGet=cms.VPSet()
    process.GlobalTag.toGet.extend( cms.VPSet(
        cms.PSet(record = cms.string("EcalSRSettingsRcd"),
                 tag = cms.string("EcalSRSettings_TL200_mc"),
                 connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_34X_ECAL")
                 )
        )
                                    )
    process.GlobalTag.toGet.extend( cms.VPSet(
        cms.PSet(record = cms.string("EcalTPGLutIdMapRcd"),
                 tag = cms.string("EcalTPGLutIdMap_beamv5_upgrade_mc"),
                 connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_34X_ECAL")
                 )
        )
                                    )
    return process

def fixEcalConditions_300(process):
    if not hasattr(process.GlobalTag,'toGet'):
        process.GlobalTag.toGet=cms.VPSet()
    process.GlobalTag.toGet.extend( cms.VPSet(
        cms.PSet(record = cms.string("EcalSRSettingsRcd"),
                 tag = cms.string("EcalSRSettings_TL300_mc"),
                 connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_34X_ECAL")
                 )
        )
                                    )
    process.GlobalTag.toGet.extend( cms.VPSet(
        cms.PSet(record = cms.string("EcalTPGLutIdMapRcd"),
                 tag = cms.string("EcalTPGLutIdMap_beamv5_upgrade_mc"),
                 connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_34X_ECAL")
                 )
        )
                                    )
    return process

def fixEcalConditions_500(process):
    if not hasattr(process.GlobalTag,'toGet'):
        process.GlobalTag.toGet=cms.VPSet()
    process.GlobalTag.toGet.extend( cms.VPSet(
        cms.PSet(record = cms.string("EcalSRSettingsRcd"),
                 tag = cms.string("EcalSRSettings_TL500_mc"),
                 connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_34X_ECAL")
                 )
        )
                                    )
    process.GlobalTag.toGet.extend( cms.VPSet(
        cms.PSet(record = cms.string("EcalTPGLutIdMapRcd"),
                 tag = cms.string("EcalTPGLutIdMap_beamv5_upgrade_mc"),
                 connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_34X_ECAL")
                 )
        )
                                    )
    return process

def fixEcalConditions_1000(process):
    if not hasattr(process.GlobalTag,'toGet'):
        process.GlobalTag.toGet=cms.VPSet()
    process.GlobalTag.toGet.extend( cms.VPSet(
        cms.PSet(record = cms.string("EcalSRSettingsRcd"),
                 tag = cms.string("EcalSRSettings_TL1000_mc"),
                 connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_34X_ECAL")
                 )
        )
                                    )
    process.GlobalTag.toGet.extend( cms.VPSet(
        cms.PSet(record = cms.string("EcalTPGLutIdMapRcd"),
                 tag = cms.string("EcalTPGLutIdMap_beamv5_3GeV_upgrade_mc"),
                 connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_34X_ECAL")
                 )
        )
                                    )
    return process

def fixEcalConditions_3000(process):
    if not hasattr(process.GlobalTag,'toGet'):
        process.GlobalTag.toGet=cms.VPSet()
    process.GlobalTag.toGet.extend( cms.VPSet(
        cms.PSet(record = cms.string("EcalSRSettingsRcd"),
                 tag = cms.string("EcalSRSettings_TL3000_mc"),
                 connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_34X_ECAL")
                 )
        )
                                    )
    process.GlobalTag.toGet.extend( cms.VPSet(
        cms.PSet(record = cms.string("EcalTPGLutIdMapRcd"),
                 tag = cms.string("EcalTPGLutIdMap_beamv5_4GeV_upgrade_mc"),
                 connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_34X_ECAL")
                 )
        )
                                    )
    return process
