from CRABClient.UserUtilities import config
config = config()

#Do read the below link, for CRAB parameters:
#https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3ConfigurationFile#CRAB_configuration_parameters
#Don't chenge this folder name while submitting jobs
config.General.workArea = 'CrabProjects'
config.General.transferOutputs = True
config.General.transferLogs = False
config.JobType.allowUndistributedCMSSW = True
config.JobType.pluginName = 'Analysis'
config.Data.inputDBS = 'global'
config.Data.unitsPerJob = 8
config.Data.totalUnits = -1
config.Data.splitting = 'LumiBased'
#config.Data.splitting = 'FileBased'
config.Data.ignoreLocality = False
config.Data.publication = False
config.Site.storageSite = 'T3_CH_CERNBOX'
#config.JobType.maxJobRuntimeMin = 500
#config.JobType.maxMemoryMB = 4000
config.General.requestName = 'Rho_delete'
config.JobType.psetName = '../cfg/chargedepptcorr_central_cfg.py'
#config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PromptReco/Cert_318939-319488_13TeV_PromptReco_SpecialCollisions18_JSON_LOWPU.txt'
#config.Data.runRange = '326381-327564'
config.Data.runRange = '326381-326400'
config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/HI/PromptReco/Cert_326381-327564_HI_PromptReco_Collisions18_JSON.txt'
config.Data.inputDataset = '/HIForward/HIRun2018A-04Apr2019-v1/AOD'
config.Data.outputDatasetTag = 'rho'
config.Data.outLFNDirBase = '/store/user/prverma'
#config.Data.outLFNDirBase = '/eos/user/p/prverma'
