import FWCore.ParameterSet.Config as cms

process = cms.Process("NtrkDistribution")
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
#process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')                                                               
process.load('Configuration.StandardSequences.ReconstructionHeavyIons_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load("RecoVertex.PrimaryVertexProducer.OfflinePrimaryVerticesRecovery_cfi")
#For using the merged pixel and general track
process.load('MergingProducer.generalAndHiPixelTracks.MergingPixAndGenProducer_cfi')

# __________________ General _________________

# Configure the logger
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10    

# Configure the number of maximum event the analyser run on in interactive mode
# -1 == ALL
process.maxEvents = cms.untracked.PSet( 
    input = cms.untracked.int32(-1) 
    #input = cms.untracked.int32(5000) 
    #input = cms.untracked.int32(500) 
)

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))


#__________________ I/O files _________________

#----- Testing one One file -------------------
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
#                                '/store/hidata/HIRun2018A/HIForward/AOD/04Apr2019-v1/00000/00DB29C1-E3F8-6546-BB32-120BFD2C7FEA.root'
'/store/group/phys_heavyions/prverma/Starlight/LheGenSim/LHE_GEN_SIM20230711/AOD_SIM_20230712/230712_073718/0000/aodsimstep_1.root'
                            ),
                        )
#----------------------------------------------

# Define output file name
import os
process.TFileService = cms.Service("TFileService",
                                   #fileName = cms.string(os.getenv('CMSSW_BASE') + '/src/Analyzers/ChargeDepAndPtCorr/test/chargeptdepcorr.root')
                                   fileName = cms.string('rhomc_genreco_july28_2.root')
                                   #fileName = cms .string('cmssw/PhysicsTools/Utilities/scripts/edmPickEvents.py')      
                               )

# __________________ Event selection _________________
process.load('HeavyIonsAnalysis.Configuration.collisionEventSelection_cff')

#process.anaStep = cms.Path(
   # process.offlinePrimaryVerticesRecovery)
#+ process.HiForest )

process.eventSelections = cms.Sequence(
    process.primaryVertexFilter +
    process.clusterCompatibilityFilter 
)

# trigger selection                                                                                                                    

import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltSelect = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltSelect.HLTPaths = [
    #"HLT_HIUPC_ZeroBias_SinglePixelTrack_v*",
    "HLT_HIZeroBias_v*"
    #"HLT_HIUPC_ZeroBias_SinglePixelTrack_MaxPixelTrack_v2"   
    # "HLT_HIUPC_ZeroBias_SinglePixelTrack_MaxPixelCluster160_v1"
]
process.hltSelect.andOr = cms.bool(True)
process.hltSelect.throw = cms.bool(False)                    

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.Geometry.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')

from Configuration.AlCa.GlobalTag import GlobalTag

process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')

process.GlobalTag = GlobalTag(process.GlobalTag, '103X_dataRun2_Prompt_v2', '')
#process.HiForest.GlobalTagLabel = process.GlobalTag.globaltag                                          

print('\n\033[31m~*~ USING CENTRALITY TABLE FOR PbPb 2018 ~*~\033[0m\n')
process.GlobalTag.snapshotTime = cms.string("9999-12-31 23:59:59.000")
process.GlobalTag.toGet.extend([
    cms.PSet(record = cms.string("HeavyIonRcd"),
             tag = cms.string("CentralityTable_HFtowers200_DataPbPb_periHYDJETshape_run2v1031x02_offline"),
             connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
             label = cms.untracked.string("HFtowers")
         ),
])

#==================================================================


# __________________ Analyse Sequence _________________

# Load you analyzer with initial configuration
process.load("Analyzers.NtrkDistribution.ntrkdist_cfi")
#process.defaultAnalysis_010   = process.CPDC010.clone()
#process.defaultAnalysis_510   = process.CPDC510.clone()
process.defaultAnalysis = process.defaultNtrkDist.clone()

# Load centrality producer for centrality calculation                                                                                                 
process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")

#here * means multiply
process.p = cms.Path(#process.defaultTrigSel *            # Select MB events                             
    #process.collisionEventSelectionPA * # PA event selection                           
    #process.olvFilter_pPb8TeV_dz1p0*    # PU filter                                    
    #process.pACentrality *
    #process.anaStep*
    process.hltSelect* 
    process.offlinePrimaryVerticesRecovery*
    #process.centralityBin*               # Centrality  
    process.eventSelections*
    #process.generalAndHiPixelTracks*    #generalAndHiPixelTracks
    process.defaultAnalysis)            # Run the analyzer  


from HLTrigger.Configuration.CustomConfigs import MassReplaceInputTag
process = MassReplaceInputTag(process,"offlinePrimaryVertices","offlinePrimaryVerticesRecovery")
process.offlinePrimaryVerticesRecovery.oldVertexLabel = "offlinePrimaryVertices"


