import FWCore.ParameterSet.Config as cms
import copy
from Configuration.Eras.Era_Run2_2016_cff import *
from Validation.CTPPS.simu_config.base_cff import *

# Define collections
muon             = 'slimmedMuons'
electron        = 'slimmedElectrons'

MC=cms.bool(True)
MC_Signal=cms.bool(False)
DATA=cms.bool(False)
preVFP=cms.bool(False)
postVFP=cms.bool(True)
from Configuration.StandardSequences.Eras import eras
if MC:
	import direct_simu_reco_cff as profile
	process = cms.Process('CTPPSTestAcceptance', profile.era)
	profile.LoadConfig(process)
	process.load("CalibPPS.ESProducers.ctppsBeamParametersESSource_cfi")
	process.load("Validation.CTPPS.simu_config.year_2016_cff")
	process.load('direct_simu_reco_cff')
else:
	process = cms.Process("CTPPSTestProtonReconstruction", eras.ctpps_2016)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Geometry.CaloEventSetup.CaloTowerConstituents_cfi")
process.load('Configuration.StandardSequences.GeometryDB_cff')
# Added in order to get the proton collections (https://twiki.cern.ch/twiki/bin/view/CMS/TaggedProtonsDirectSimulation#Full_step_by_step_recipe_to_setu)
process.load('Configuration.EventContent.EventContent_cff')

process.maxEvents = cms.untracked.PSet( 
			input = cms.untracked.int32(-1) 
)

process.options=cms.untracked.PSet(
	wantSummary=cms.untracked.bool(True)
	, SkipEvent = cms.untracked.vstring('ProductNotFound')
#	, numberOfThreads = cms.untracked.uint32( 8 )
)
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring("file:MC_preVFP.root"),#
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
)

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
#Joao Pedro changed:
if MC and preVFP:
        process.GlobalTag = GlobalTag(process.GlobalTag, '106X_mcRun2_asymptotic_preVFP_v9','') # MC: https://twiki.cern.ch/twiki/bin/view/CMS/PdmVLegacy2016preVFPAnalysis
elif MC and postVFP:
        process.GlobalTag = GlobalTag(process.GlobalTag, '106X_mcRun2_asymptotic_v15', '') # MC: https://twiki.cern.ch/twiki/bin/view/CMS/PdmVLegacy2016postVFPAnalysis
else:
        process.GlobalTag = GlobalTag(process.GlobalTag, '106X_dataRun2_v32') # data

# local RP reconstruction chain with standard settings
process.load("RecoCTPPS.Configuration.recoCTPPS_cff")


##################      P R O T O N    ####################
if MC:
	# override LHCInfo source         
	UseCrossingAngle(185., process)

	# update settings of beam-smearing module
	process.beamDivergenceVtxGenerator.src = cms.InputTag("")
	process.beamDivergenceVtxGenerator.srcGenParticle = cms.VInputTag(
    		cms.InputTag("genPUProtons","genPUProtons")
  		, cms.InputTag("prunedGenParticles")	
	)

	# do not apply vertex smearing again                                                                                                                   
	process.ctppsBeamParametersESSource.vtxStddevX = 0
	process.ctppsBeamParametersESSource.vtxStddevY = 0
	process.ctppsBeamParametersESSource.vtxStddevZ = 0

	# undo CMS vertex shift                                                                                                                                
	process.ctppsBeamParametersESSource.vtxOffsetX45 = -0.1077 
	process.ctppsBeamParametersESSource.vtxOffsetY45 = -0.418
	process.ctppsBeamParametersESSource.vtxOffsetZ45 = +1.576

###################    T R I G G E R    ###################
from HLTrigger.HLTfilters.hltHighLevel_cfi import *
process.hltFilter = copy.deepcopy(hltHighLevel)
process.hltFilter.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
process.hltFilter.HLTPaths = ['HLT_Ele27_WPTight_Gsf_v*','HLT_IsoMu24_v*']
process.hltFilter.throw = cms.bool(False)
process.hltFilter.andOr = cms.bool(True) # how to deal with multiple triggers: True (OR) accept if ANY is true, False (AND) accept if ALL are true

from HLTrigger.HLTfilters.hltHighLevel_cfi import *
process.METFilter = copy.deepcopy(hltHighLevel) #HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
if MC:
	process.METFilter.TriggerResultsTag =  cms.InputTag("TriggerResults","","PAT") # cms.InputTag("TriggerResults","","RECO")
else:
	process.METFilter.TriggerResultsTag =  cms.InputTag("TriggerResults","","RECO")
process.METFilter.throw = cms.bool(False) # throw exception on unknown path names
process.METFilter.andOr = cms.bool(False) # how to deal with multiple triggers: True (OR) accept if ANY is true, False (AND) accept if ALL are true
# Filters from here: https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2
# MET filter recommendations will be updated for UltraLegacy datasets. The filters are currently under review and it is likely that the recommendations for ECAL related filters will change. Updates recommendations to be released as soon as they are available.
if MC:
	process.METFilter.HLTPaths = ['Flag_globalSuperTightHalo2016Filter','Flag_HBHENoiseFilter','Flag_HBHENoiseIsoFilter','Flag_globalTightHalo2016Filter','Flag_EcalDeadCellTriggerPrimitiveFilter','Flag_goodVertices'] #MC
else:
        process.METFilter.HLTPaths = ['Flag_globalSuperTightHalo2016Filter','Flag_HBHENoiseFilter','Flag_HBHENoiseIsoFilter','Flag_globalTightHalo2016Filter','Flag_EcalDeadCellTriggerPrimitiveFilter','Flag_goodVertices','Flag_eeBadScFilter'] #data

###########################################################

# COMMENTED BY JOAO PEDRO 18/02/21

#################################
  ###  JET TOOLBOX FOR CHS ###
#################################
# AK R=0.8 jets from CHS inputs with basic grooming, W tagging, and top tagging                                                            
from JMEAnalysis.JetToolbox.jetToolbox_cff import *
jetToolbox( process, 'ak8', 'ak8JetSubs', 'noOutput',
#jetToolbox( process, 'ak8', 'jetSequence', 'noOutput',
                PUMethod='CHS',runOnMC=MC,
                addPruning=True, addSoftDrop=False ,           # add basic grooming                                                            
                addTrimming=False, addFiltering=False,
                addSoftDropSubjets=False,
                addNsub=True, maxTau=4,                       # add Nsubjettiness tau1, tau2, tau3, tau4                                      
                Cut='pt > 100.0',
                bTagDiscriminators=['pfCombinedInclusiveSecondaryVertexV2BJetTags'],
                #bTagDiscriminators=['pfCombinedSecondaryVertexV2BJetTags'],
                # added L1FastJet on top of the example config file
                JETCorrPayload = 'AK8PFchs', JETCorrLevels = ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']
                )
jetToolbox( process, 'ak4', 'ak4JetSubs', 'noOutput',
#jetToolbox( process, 'ak4', 'jetSequence', 'noOutput',
                PUMethod='CHS',runOnMC=MC,
                addPruning=True, addSoftDrop=False ,           # add basic grooming                                                            
                addTrimming=False, addFiltering=False,
                addSoftDropSubjets=False,
                addNsub=True, maxTau=4,                       # add Nsubjettiness tau1, tau2, tau3, tau4                                      
                Cut='pt > 10.0',
                bTagDiscriminators=['pfCombinedInclusiveSecondaryVertexV2BJetTags'],
                #bTagDiscriminators=['pfCombinedSecondaryVertexV2BJetTags'],
                # added L1FastJet on top of the example config file
                JETCorrPayload = 'AK4PFchs', JETCorrLevels = ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']
                )

if MC:
    #################################
    ###  JER SMEARING AK8###
    #################################
    from RecoMET.METProducers.METSigParams_cfi import *
    process.slimmedAK8JetsSmeared = cms.EDProducer('SmearedPATJetProducer',
        src = cms.InputTag("selectedPatJetsAK8PFCHS"),
        enabled = cms.bool(True),
        rho = cms.InputTag("fixedGridRhoFastjetAll"),
        algo = cms.string('AK8PFchs'),
        algopt = cms.string('AK8PFchs_pt'),
        genJets = cms.InputTag('slimmedGenJets'),
        dRMax = cms.double(0.2),
        dPtMaxFactor = cms.double(3),
        seed = cms.uint32(37428479),
        debug = cms.untracked.bool(False),
    # Systematic variation
    # 0: Nominal
    # -1: -1 sigma (down variation)
    # 1: +1 sigma (up variation)
     variation = cms.int32(0)  # If not specified, default to 0
       )
    #################################
    ###  JER SMEARING AK4###
    #################################
    from RecoMET.METProducers.METSigParams_cfi import *
    process.slimmedAK4JetsSmeared = cms.EDProducer('SmearedPATJetProducer',
        src = cms.InputTag("selectedPatJetsAK4PFCHS"),
        enabled = cms.bool(True),
        rho = cms.InputTag("fixedGridRhoFastjetAll"),
        algo = cms.string('AK4PFchs'),
        algopt = cms.string('AK4PFchs_pt'),
        genJets = cms.InputTag('slimmedGenJets'),
        dRMax = cms.double(0.2),
        dPtMaxFactor = cms.double(3),
        seed = cms.uint32(37424479),
        debug = cms.untracked.bool(False),
    # Systematic variation
    # 0: Nominal
    # -1: -1 sigma (down variation)
    # 1: +1 sigma (up variation)
    variation = cms.int32(0)  # If not specified, default to 0
)

if MC:
	JetAK4Collection="slimmedAK4JetsSmeared"
	JetAK8Collection="slimmedAK8JetsSmeared"
else:
        JetAK4Collection="selectedPatJetsAK4PFCHS"
        JetAK8Collection="selectedPatJetsAK8PFCHS"

###################    C L E A N E R    ###################
process.cleanPatMuons = cms.EDProducer("PATMuonCleaner",
    src = cms.InputTag(muon),
    # preselection (any string-based cut for pat::Muon)
    preselection = cms.string("passed('CutBasedIdTight') && passed('PFIsoTight') && pt() > 20. && abs(eta()) < 2.4"),
    # overlap checking configurables
    checkOverlaps = cms.PSet(),
    # finalCut (any string-based cut for pat::Muon)
    finalCut = cms.string(''),
)

process.cleanPatElectrons = cms.EDProducer("PATElectronCleaner",
    ## pat electron input source
    src = cms.InputTag(electron),
    # preselection (any string-based cut for pat::Electron)
    preselection = cms.string("electronID('cutBasedElectronID-Fall17-94X-V1-tight') && pt() > 20. && abs(eta()) < 2.4"), #electronID('egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-tight
    # overlap checking configurables
    checkOverlaps = cms.PSet(),
    # finalCut (any string-based cut for pat::Electron)
    finalCut = cms.string(''),
)

process.cleanJets = cms.EDProducer("PATJetCleaner",
    src = cms.InputTag(JetAK4Collection),

    # preselection (any string-based cut on pat::Jet)
    preselection = cms.string(''),

    # overlap checking configurables
    checkOverlaps = cms.PSet(
        muons = cms.PSet(
           src       = cms.InputTag("cleanPatMuons"),
           algorithm = cms.string("byDeltaR"),
           preselection        = cms.string(""),
           deltaR              = cms.double(0.3),
           checkRecoComponents = cms.bool(False), # don't check if they share some AOD object ref
           pairCut             = cms.string(""),
           requireNoOverlaps   = cms.bool(True), # overlaps don't cause the jet to be discared
        ),
        electrons = cms.PSet(
           src       = cms.InputTag("cleanPatElectrons"),
           algorithm = cms.string("byDeltaR"),
           preselection        = cms.string(""),
           deltaR              = cms.double(0.3),
           checkRecoComponents = cms.bool(False), # don't check if they share some AOD object ref
           pairCut             = cms.string(""),
           requireNoOverlaps   = cms.bool(True), # overlaps don't cause the jet to be discared
        ),
    ),
    # finalCut (any string-based cut on pat::Jet)
    finalCut = cms.string(''),
)
process.cleanJetsAK8 = cms.EDProducer("PATJetCleaner",
    src = cms.InputTag(JetAK8Collection), # selectedPatJetsAK8PFCHS
    # preselection (any string-based cut on pat::Jet)
    preselection = cms.string(''),

    # overlap checking configurables
    checkOverlaps = cms.PSet(
        muons = cms.PSet(
           src       = cms.InputTag("cleanPatMuons"),
           algorithm = cms.string("byDeltaR"),
           preselection        = cms.string(""),
           deltaR              = cms.double(1.0),
           checkRecoComponents = cms.bool(False), # don't check if they share some AOD object ref
           pairCut             = cms.string(""),
           requireNoOverlaps   = cms.bool(True), # overlaps don't cause the jet to be discared
        ),
        electrons = cms.PSet(
           src       = cms.InputTag("cleanPatElectrons"),
           algorithm = cms.string("byDeltaR"),
           preselection        = cms.string(""),
           deltaR              = cms.double(1.0),
           checkRecoComponents = cms.bool(False), # don't check if they share some AOD object ref
           pairCut             = cms.string(""),
           requireNoOverlaps   = cms.bool(True), # overlaps don't cause the jet to be discared
        ),
    ),
    # finalCut (any string-based cut on pat::Jet)
    finalCut = cms.string(''),
)
###########################################################

#################### P R E F I R I N G ####################
from PhysicsTools.PatUtils.l1ECALPrefiringWeightProducer_cfi import l1ECALPrefiringWeightProducer
process.prefiringweight = l1ECALPrefiringWeightProducer.clone(
    TheJets = cms.InputTag(JetAK8Collection) #this should be the slimmedJets collection with up to date JECs !
    , DataEra = cms.string("2016BtoH") #Use 2016BtoH for 2016
    , UseJetEMPt = cms.bool(False)
    , PrefiringRateSystematicUncty = cms.double(0.2)
    , SkipWarnings = False)

############# E / G A M M A   P O S T   R E C O ###########
from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(
		process
		,era='2016postVFP-UL' # change era
) 

###################   M E T   C O R R     #################

# MET
if MC:
	process.genMet = cms.EDProducer("GenMETExtractor",
        	metSource = cms.InputTag("slimmedMETs", "", "@skipCurrentProcess")
        	)

# Raw MET
process.uncorrectedMet = cms.EDProducer("RecoMETExtractor",
        correctionLevel = cms.string('raw'),
        metSource = cms.InputTag("slimmedMETs", "", "@skipCurrentProcess")
        )

# Raw PAT MET
from PhysicsTools.PatAlgos.tools.metTools import addMETCollection
addMETCollection(process, labelName="uncorrectedPatMet", metSource="uncorrectedMet")
if MC:
	process.uncorrectedPatMet.genMETSource = cms.InputTag('genMet')
else: 
	process.uncorrectedPatMet.addGenMET = False

#Type-1 correction
process.Type1CorrForNewJEC = cms.EDProducer("PATPFJetMETcorrInputProducer",
        src = cms.InputTag("selectedPatJetsAK4PFCHS"),
        jetCorrLabel = cms.InputTag("L3Absolute"),
        jetCorrLabelRes = cms.InputTag("L2L3Residual"),
        offsetCorrLabel = cms.InputTag("L1FastJet"),
        skipEM = cms.bool(True),
        skipEMfractionThreshold = cms.double(0.9),
        skipMuonSelection = cms.string('isGlobalMuon | isStandAloneMuon'),
        skipMuons = cms.bool(True),
        type1JetPtThreshold = cms.double(15.0)
        )
if MC:
	process.slimmedMETsNewJEC = cms.EDProducer('CorrectedPATMETProducer',
       	 	src = cms.InputTag('uncorrectedPatMet'),
        	srcCorrections = cms.VInputTag(cms.InputTag('Type1CorrForNewJEC', 'type1'))
        	)
else:
        process.slimmedMETsNewJEC = cms.EDProducer('CorrectedPATMETProducer',
 		src = cms.InputTag('uncorrectedPatMet'),
                srcCorrections = cms.VInputTag(cms.InputTag('Type1CorrForNewJEC', 'type1')),
		applyType2Corrections = cms.bool(False)
                )

if MC:
	process.shiftedMETCorrModuleForSmearedJets = cms.EDProducer('ShiftedParticleMETcorrInputProducer',
        	srcOriginal = cms.InputTag("selectedPatJetsAK4PFCHS"),
         	srcShifted = cms.InputTag(JetAK4Collection) 
        	)
	process.slimmedMETsSmeared = cms.EDProducer('CorrectedPATMETProducer',
       	 	src = cms.InputTag('slimmedMETsNewJEC'),
       	 	srcCorrections = cms.VInputTag(cms.InputTag('shiftedMETCorrModuleForSmearedJets'))
       	 	)

if MC:
	METCollection="slimmedMETsSmeared"
else:
	METCollection="slimmedMETsNewJEC"

# Create a EDFilter to choose for 2 oposite charge leptons with |eta|<2.4 and pt>38GeV and loose ID, primary vertex at |z| < 15cm.
process.demo2 = cms.EDFilter('MyUserSelector'
        , muons                    = cms.InputTag(muon)
        , electrons                = cms.InputTag(electron)
        , vertices                 = cms.InputTag('offlineSlimmedPrimaryVertices')
	, Z_vtx			   = cms.double(15.)
	, pT_min		   = cms.double(38.)
	, eta_max		   = cms.double(2.4)	
        , filter = cms.bool(True)
)

# demo process, which processes MakeNTuple:
process.demo = cms.EDAnalyzer('MakeNTuple'
        , muons                    = cms.InputTag(muon)
        , electrons                = cms.InputTag(electron)
        , MET                      = cms.InputTag(METCollection)
        , PFCand                   = cms.InputTag('packedPFCandidates')
        , PileupSumInfoInputTag = cms.InputTag('slimmedAddPileupInfo')
        , vertices                      = cms.InputTag('offlineSlimmedPrimaryVertices')
	, MC			   = MC
        , MC_Signal                = MC_Signal
        , DATA                     = DATA
        , preVFP_                  = preVFP
        , postVFP_                 = postVFP
        , ppsRecoProtonSingleRPTag = cms.InputTag("ctppsProtons", "singleRP")
        , ppsRecoProtonMultiRPTag  = cms.InputTag("ctppsProtons", "multiRP")
	, TriggerResults = cms.InputTag('TriggerResults', '', 'HLT')
        , genTag  = cms.InputTag('prunedGenParticles')
	, MCEvent = cms.untracked.InputTag('source','generator')
)

# Output
process.TFileService = cms.Service("TFileService",
                                       fileName = cms.string("out.root")
)

process.ctppsLocalTrackLiteProducer.includePixels = cms.bool(False)

# After the modification by JOAO PEDRO 18/02/21:
if (MC_Signal==False and MC==True):
        process.p = cms.Path(
                 process.hltFilter
                *process.METFilter
                *process.slimmedAK8JetsSmeared
                *process.slimmedAK4JetsSmeared
                *process.cleanPatElectrons
                *process.cleanPatMuons
                *process.cleanJets
                *process.cleanJetsAK8
                *process.genMet
                *process.uncorrectedMet
                *process.uncorrectedPatMet
                *process.Type1CorrForNewJEC
                *process.slimmedMETsNewJEC
                *process.shiftedMETCorrModuleForSmearedJets
                *process.slimmedMETsSmeared
                *process.egammaPostRecoSeq
		*process.prefiringweight
                *process.demo2
                *process.demo)
elif MC_Signal==True:
        process.p = cms.Path(
                 process.beamDivergenceVtxGenerator
                *process.ctppsDirectProtonSimulation
                *process.reco_local
                *process.ctppsProtons
                *process.hltFilter
                *process.METFilter
                *process.slimmedAK8JetsSmeared
                *process.slimmedAK4JetsSmeared
                *process.cleanPatElectrons
                *process.cleanPatMuons
                *process.cleanJets
                *process.cleanJetsAK8
                *process.genMet
                *process.uncorrectedMet
                *process.uncorrectedPatMet
                *process.Type1CorrForNewJEC
                *process.slimmedMETsNewJEC
                *process.shiftedMETCorrModuleForSmearedJets
                *process.slimmedMETsSmeared
                *process.egammaPostRecoSeq
		*process.prefiringweight
                *process.demo2
                *process.demo)
else:
	process.p = cms.Path(
		process.hltFilter
		*process.METFilter
		*process.cleanPatElectrons
		*process.cleanPatMuons
		*process.cleanJets
		*process.cleanJetsAK8
		*process.uncorrectedMet
		*process.Type1CorrForNewJEC
		*process.slimmedMETsNewJEC
		*process.egammaPostRecoSeq
		*process.prefiringweight
		*process.demo2
		*process.demo)

