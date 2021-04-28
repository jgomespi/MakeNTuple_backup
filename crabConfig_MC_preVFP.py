#if __name__ == '__main__':

from CRABAPI.RawCommand import crabCommand
from CRABClient.ClientExceptions import ClientException
from CRABClient.UserUtilities import config 
from httplib import HTTPException
#from WMCore.Configuration import Configuration
config = config()

config.section_("General")
config.General.workArea = '/afs/cern.ch/user/j/jgomespi/private/workspace/master_thesis/FullLep_Studies/2016_analysis/Ntuple_build/tutorial/Ntuple_building/CMSSW_10_6_20/src/MakeNTuple/MakeNTuple/crab_projects_MC/'
config.General.transferOutputs = True
config.General.transferLogs = False

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'ConfFile_MC_preVFP_cfg.py'
config.JobType.inputFiles = ['MyDataPileupHistogram.root', 'PileupMC.root','direct_simu_reco_cff.py','RoccoR2016aUL.txt','RoccoR2016bUL.txt' ]
config.JobType.outputFiles = ['out.root']
config.JobType.allowUndistributedCMSSW = True	
config.JobType.maxMemoryMB = 2500
#config.JobType.numCores = 8

config.section_("Data")
config.Data.inputDBS = 'global'
#config.Data.splitting = 'Automatic'
config.Data.splitting = 'FileBased'  #'Automatic'
config.Data.unitsPerJob = 1
config.Data.outLFNDirBase = '/store/user/jgomespi' #%s/' % (getUsernameFromSiteDB())
config.Data.publication = False 
config.Data.outputDatasetTag = 'MC'
	
config.section_("Site")
config.Site.storageSite = 'T2_BR_UERJ'

def submit(config):
	try:
        	crabCommand('submit', config = config)
	except HTTPException as hte:
		print "Failed submitting task: %s" % (hte.headers)
	except ClientException as cle:
		print "Failed submitting task: %s" % (cle)

config.General.requestName = 'DYJetsToLL_M-50_preVFP'
config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/MINIAODSIM'
submit(config)

#config.General.requestName = 'DYJetsToLL_M-50_postVFP'
#config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL16MiniAOD-106X_mcRun2_asymptotic_v13-v2/MINIAODSIM'
#submit(config)

config.General.requestName = 'WW_preVFP'
config.Data.inputDataset = '/WW_TuneCP5_13TeV-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/MINIAODSIM'
submit(config)

#config.General.requestName = 'WW_postVFP'
#config.Data.inputDataset = '/WW_TuneCP5_13TeV-pythia8/RunIISummer20UL16MiniAOD-106X_mcRun2_asymptotic_v13-v2/MINIAODSIM'
#submit(config)
