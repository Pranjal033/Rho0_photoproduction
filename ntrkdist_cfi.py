import FWCore.ParameterSet.Config as cms

defaultNtrkDist = cms.EDAnalyzer('NtrkDistribution', #Analyzer named: Correspond to the class name in 'plugin' folder
                             #Track collection
                                 tracks    = cms.InputTag('generalTracks'),
                                 #tracks     = cms.InputTag('hiConformalPixelTracks'),
                                 #tracks    = cms.InputTag('generalAndHiPixelTracks'),
                                 tracksgen = cms.InputTag("genParticles"),
                                 #Vertex collection
                                 vertex    = cms.InputTag('offlinePrimaryVerticesRecovery'),
                                 #Calorimeter tower collection
                                 caloTower = cms.InputTag('towerMaker'),
                                 #Centrality
                                 centralitySrc    = cms.InputTag("hiCentrality"),
                                 centralityBinSrc = cms.InputTag("centralityBin","HFtowers"),
                                 #Vertex selection
                                 minvz         = cms.untracked.double(-15.0), 
                                 maxvz         = cms.untracked.double(15.0),
                                 maxrho        = cms.untracked.double(0.2),
                                 isBVselByMult = cms.untracked.bool(False),
                                 #Multiplicity selection
                                 noffmin       = cms.untracked.int32(0),
                                 noffmax       = cms.untracked.int32(10000),
                                 ptnoffmin     = cms.untracked.double(0.0),
                                 ptnoffmax     = cms.untracked.double(5),
                                 dzdzerrornoff = cms.untracked.double(3.0),
                                 d0d0errornoff = cms.untracked.double(3.0),
                                 pterrorptnoff = cms.untracked.double(0.1),
                                 #Track selection
                                 etamin    = cms.untracked.double(-2.4),
                                 etamax    = cms.untracked.double(2.4),
                                 ptmin     = cms.untracked.double(0.5),
                                 ptmax     = cms.untracked.double(5.0),
                                 dzdzerror = cms.untracked.double(3.0),
                                 d0d0error = cms.untracked.double(3.0),
                                 pterrorpt = cms.untracked.double(0.1),
                                 #Acc X Eff
                                 fname = cms.untracked.InputTag('Hijing_8TeV_dataBS.root'),
                                 effmultbin = cms.untracked.vint32(0,10000)
                             )
