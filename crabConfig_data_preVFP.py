from CRABAPI.RawCommand import crabCommand
from CRABClient.ClientExceptions import ClientException
#from CRABClient.UserUtilities #import config, getUsernameFromSiteDB
from httplib import HTTPException
from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
#config.General.requestName = 'SingleMuonB'
config.General.workArea = '/afs/cern.ch/user/j/jgomespi/private/workspace/master_thesis/FullLep_Studies/2016_analysis/Ntuple_build/tutorial/Ntuple_building/CMSSW_10_6_20/src/MakeNTuple/MakeNTuple/crab_projects_Data'
config.General.transferOutputs = True
config.General.transferLogs = False

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'ConfFile_data_preVFP_cfg.py'
config.JobType.inputFiles = ['MyDataPileupHistogram.root', 'PileupMC.root','RoccoR2016aUL.txt','RoccoR2016bUL.txt' ]
config.JobType.outputFiles = ['out.root']
config.JobType.allowUndistributedCMSSW = True

config.section_("Data")
#config.Data.inputDataset = '/SingleMuon/Run2016B-17Jul2018_ver2-v1/MINIAOD'
#config.Data.useParent = True
config.Data.ignoreLocality = True
config.Data.inputDBS = 'global'
#config.Data.splitting = 'Automatic'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 50
config.Data.outLFNDirBase = '/store/user/jgomespi' 
config.Data.lumiMask = '/afs/cern.ch/user/j/jgomespi/private/workspace/master_thesis/FullLep_Studies/2016_analysis/Ntuple_build/tutorial/Ntuple_building/CMSSW_10_6_20/src/MakeNTuple/MakeNTuple/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON_PPSruns.txt'
config.Data.publication = False 
config.Data.outputDatasetTag = 'data'

config.section_("Site")
config.Site.whitelist = ['T2_CH_*','T2_DE_*','T2_IT_*','T2_US_*']
config.Site.storageSite = 'T2_BR_UERJ'

def submit(config):
    try:
        crabCommand('submit', config = config)
    except HTTPException as hte:
        print "Failed submitting task: %s" % (hte.headers)
    except ClientException as cle:
        print "Failed submitting task: %s" % (cle)

#############################################################################################
## From now on that's what users should modify: this is the a-la-CRAB2 configuration part. ##
#############################################################################################

config.General.requestName = 'SingleMuon_B'
config.Data.inputDataset = '/SingleMuon/Run2016B-21Feb2020_ver1_UL2016_HIPM-v1/MINIAOD'
submit(config)

config.General.requestName = 'SingleMuon_B_2'
config.Data.inputDataset = '/SingleMuon/Run2016B-21Feb2020_ver2_UL2016_HIPM-v1/MINIAOD'
submit(config)

config.General.requestName = 'SingleMuon_C'
config.Data.inputDataset = '/SingleMuon/Run2016C-21Feb2020_UL2016_HIPM-v1/MINIAOD'
submit(config)

config.General.requestName = 'SingleElectron_B'
config.Data.inputDataset = '/SingleElectron/Run2016B-21Feb2020_ver1_UL2016_HIPM-v1/MINIAOD'
submit(config)

config.General.requestName = 'SingleElectron_B_2'
config.Data.inputDataset = '/SingleElectron/Run2016B-21Feb2020_ver2_UL2016_HIPM-v1/MINIAOD'
submit(config)

config.General.requestName = 'SingleElectron_C'
config.Data.inputDataset = '/SingleElectron/Run2016C-21Feb2020_UL2016_HIPM-v1/MINIAOD'
submit(config)

