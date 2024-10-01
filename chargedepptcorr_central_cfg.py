import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Era_Run3_2023_UPC_cff import Run3_2023_UPC
process = cms.Process("NtrkDistribution", Run3_2023_UPC)

process.load('Configuration.Geometry.GeometryDB_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')

# __________________ General _________________

# Configure the logger
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

# Configure the number of maximum event the analyser run on in interactive mode
# -1 == ALL
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))


#__________________ I/O files _________________

#----- Testing one One file -------------------
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
                            '/store/hidata/HIRun2023A/HIForward0/AOD/16Jan2024-v1/40000/7efdd2e1-1597-424f-b77f-22935545155d.root'
                            ),
)
#----------------------------------------------

# Define output file name
import os
process.TFileService = cms.Service("TFileService",
            fileName = cms.string('Output.root')
)
#___________________Global Tag________________________

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '132X_dataRun3_Prompt_HI_LowPtPhotonReg_v2', '')

#____________________HFCUT___________________________________________

process.towersAboveThreshold = cms.EDFilter("GenericPFCandidateSelector",
                src = cms.InputTag("particleFlow"),
                cut = cms.string("(pdgId() == 1 || pdgId == 2) && abs(eta) >= 3.0 && energy() >= 10")
)
# select HF+ towers above threshold                                                                      
process.hfPosTowers = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("towersAboveThreshold"),
    cut = cms.string("eta() >= 3.0 && eta() <= 6.0")
)
# select HF- towers above threshold                                                                      
process.hfNegTowers = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("towersAboveThreshold"),
    cut = cms.string("eta() >= -6.0 && eta() <= -3.0")
)
# require at least one HF+ tower above threshold                                                         
process.hfPosFilter = cms.EDFilter("CandCountFilter",
    src = cms.InputTag("hfPosTowers"),
    minNumber = cms.uint32(1)
)
# require at least one HF- tower above threshold                                                         
process.hfNegFilter = cms.EDFilter("CandCountFilter",
    src = cms.InputTag("hfNegTowers"),
    minNumber = cms.uint32(1)
)
process.hfFilterN_seq = cms.Sequence(
    process.towersAboveThreshold *
    process.hfPosTowers *
    process.hfNegTowers 
    ~process.hfPosFilter *
    ~process.hfNegFilter
)

#__________________________JSON___________________________ 
import FWCore.PythonUtilities.LumiList as LumiList

#Json file is not affecting, why?
process.source.lumisToProcess = LumiList.LumiList(filename = '/eos/user/c/cmsdqm/www/CAF/certification/Collisions23HI/Cert_Collisions2023HI_374288_375823_Muon.json').getVLuminosityBlockRange()

# trigger selection
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltSelect = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltSelect.HLTPaths = [
    #"HLT_HIL1NotBptxOR_v*" #Empty BX
    #"HLT_HIZeroBias_HighRateRAW_v*" #Zero Bias dataset
    "HLT_HIUPC_ZeroBias_SinglePixelTrackLowPt_MaxPixelCluster400_v*"
    #"HLT_HIUPC_ZDC1nOR_SinglePixelTrackLowPt_MaxPixelCluster400_v*"
    #"HLT_HIUPC_SingleMuOpen_*"
]
process.hltSelect.andOr = cms.bool(True)
process.hltSelect.throw = cms.bool(False)                    

### centrality ###                                                                       
process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
process.centralityBin.Centrality = cms.InputTag("hiCentrality")
process.centralityBin.centralityVariable = cms.string("HFtowers")

# ZDC RecHit Producer
process.load('HeavyIonsAnalysis.ZDCAnalysis.QWZDC2018Producer_cfi')
process.load('HeavyIonsAnalysis.ZDCAnalysis.QWZDC2018RecHit_cfi')
process.load('HeavyIonsAnalysis.ZDCAnalysis.zdcanalyzer_cfi')

process.zdcdigi.SOI = cms.untracked.int32(2)
process.zdcanalyzer.doZDCRecHit = False
process.zdcanalyzer.doZDCDigi = True
process.zdcanalyzer.zdcRecHitSrc = cms.InputTag("QWzdcreco")
process.zdcanalyzer.zdcDigiSrc = cms.InputTag("hcalDigis", "ZDC")
process.zdcanalyzer.calZDCDigi = False
process.zdcanalyzer.verbose = False
process.zdcanalyzer.nZdcTs = cms.int32(6)

from CondCore.CondDB.CondDB_cfi import *
process.es_pool = cms.ESSource("PoolDBESSource",
                    timetype = cms.string('runnumber'),
                    toGet = cms.VPSet(
                    cms.PSet(
                    record = cms.string("HcalElectronicsMapRcd"),
                    tag = cms.string("HcalElectronicsMap_2021_v2.0_data")
                )
                ),
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
        authenticationMethod = cms.untracked.uint32(1)
)

process.es_prefer = cms.ESPrefer('HcalTextCalibrations', 'es_ascii')
process.es_ascii = cms.ESSource(
'HcalTextCalibrations',
input = cms.VPSet(
    cms.PSet(
        object = cms.string('ElectronicsMap'),
        file = cms.FileInPath("emap_2023_newZDC_v3.txt")
    )
)
)

#CCF and PVF:
# __________________ Event selection _________________

process.load('HeavyIonsAnalysis.EventAnalysis.collisionEventSelection_cff')
process.load('HeavyIonsAnalysis.TrackAnalysis.TrackAnalyzer_cfi')
process.eventFilters = cms.Sequence(
    process.primaryVertexFilter *
    process.clusterCompatibilityFilter
)



# __________________ Analyse Sequence _________________

# Load you analyzer with initial configuration
process.load("Analyzers.NtrkDistribution.ntrkdist_cfi")
process.defaultAnalysis = process.defaultNtrkDist.clone()

process.p = cms.Path(
                     process.eventFilters*
                     process.hltSelect*
                     process.hfFilterN_seq*
                     process.defaultAnalysis)

