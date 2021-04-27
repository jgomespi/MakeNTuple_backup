import FWCore.ParameterSet.Config as cms
import copy

# Define collections
muon             = 'slimmedMuons'
electron        = 'slimmedElectrons'

MC=cms.bool(True)
MC_Signal=cms.bool(False)
DATA=cms.bool(False)


process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(3000) )

process.options=cms.untracked.PSet(wantSummary=cms.untracked.bool(True))

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring('file:output_numEvent10000.root'),
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
    secondaryFileNames = cms.untracked.vstring()
)


process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
if MC:
	process.GlobalTag = GlobalTag(process.GlobalTag, '94X_mcRun2_asymptotic_v3', '') # MC
else:
	process.GlobalTag = GlobalTag(process.GlobalTag, '94X_dataRun2_v10', '') # data




###################    T R I G G E R    ###################
from HLTrigger.HLTfilters.hltHighLevel_cfi import *
process.hltFilter = copy.deepcopy(hltHighLevel)
process.hltFilter.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
process.hltFilter.HLTPaths = ['HLT_Mu50_v*','HLT_Ele27_WPTight_Gsf_v*','HLT_IsoMu24_v*']
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
#process.METFilter.HLTPaths = ['Flag_goodVertices','Flag_globalSuperTightHalo2016Filter','Flag_HBHENoiseFilter','Flag_HBHENoiseIsoFilter','Flag_EcalDeadCellTriggerPrimitiveFilter','Flag_BadPFMuonFilter','Flag_BadChargedCandidateFilter','Flag_eeBadScFilteri']
if MC:
	process.METFilter.HLTPaths = ['Flag_globalSuperTightHalo2016Filter','Flag_HBHENoiseFilter','Flag_HBHENoiseIsoFilter','Flag_globalTightHalo2016Filter','Flag_EcalDeadCellTriggerPrimitiveFilter','Flag_goodVertices'] #MC
else:
        process.METFilter.HLTPaths = ['Flag_globalSuperTightHalo2016Filter','Flag_HBHENoiseFilter','Flag_HBHENoiseIsoFilter','Flag_globalTightHalo2016Filter','Flag_EcalDeadCellTriggerPrimitiveFilter','Flag_goodVertices','Flag_eeBadScFilter'] #data
#process.METFilter.HLTPaths = ['Flag_HBHENoiseFilter','Flag_HBHENoiseIsoFilter','Flag_globalTightHalo2016Filter','Flag_EcalDeadCellTriggerPrimitiveFilter','Flag_goodVertices'] #FastSim

process.load('RecoMET.METFilters.ecalBadCalibFilter_cfi')

baddetEcallist = cms.vuint32(
    [872439604,872422825,872420274,872423218,
     872423215,872416066,872435036,872439336,
     872420273,872436907,872420147,872439731,
     872436657,872420397,872439732,872439339,
     872439603,872422436,872439861,872437051,
     872437052,872420649,872422436,872421950,
     872437185,872422564,872421566,872421695,
     872421955,872421567,872437184,872421951,
     872421694,872437056,872437057,872437313])


process.ecalBadCalibReducedMINIAODFilter = cms.EDFilter(
    "EcalBadCalibFilter",
    EcalRecHitSource = cms.InputTag("reducedEgamma:reducedEERecHits"),
    ecalMinEt        = cms.double(50.),
    baddetEcal    = baddetEcallist, 
    taggingMode = cms.bool(True),
    debug = cms.bool(False)
    )
###########################################################




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
    preselection = cms.string("passed('CutBasedIdTight') && passed('PFIsoTight') && pt() > 50 && abs(eta()) < 2.4"),
#    preselection = cms.string("muonID('GlobalMuonPromptTight') && numberOfMatchedStations() > 1 && dB() < 0.2 && pt() > 50 && abs(eta()) < 2.4 && ((pfIsolationR04().sumChargedHadronPt + max(0., pfIsolationR04().sumNeutralHadronEt + pfIsolationR04().sumPhotonEt - 0.5*pfIsolationR04().sumPUPt))/pt() < 0.1)"),  # "&& ((isolationR03().sumPt)/(pt()) < 0.1)",
    # overlap checking configurables
    checkOverlaps = cms.PSet(),
    # finalCut (any string-based cut for pat::Muon)
    finalCut = cms.string(''),
)

process.cleanPatElectrons = cms.EDProducer("PATElectronCleaner",
    ## pat electron input source
    src = cms.InputTag(electron),
    # preselection (any string-based cut for pat::Electron)
    preselection = cms.string("electronID('cutBasedElectronID-Summer16-80X-V1-tight') && pt() > 50 && abs(eta()) < 2.4"), #electronID('egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-tight
        #(pfIsolationVariables().sumChargedHadronPt()+max(0.0,pfIsolationVariables().sumNeutralHadronEt()+pfIsolationVariables().sumPhotonEt()-rho*effArea))/pt())
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

# Type-1 correction
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
#        METCollection="Type1CorrForNewJEC"

###########################################################



###################  J E T  F I L T E R   #################
process.filterJets = cms.EDFilter( 'FilterAK8Jet'
        , jetsAK8                  = cms.InputTag(JetAK8Collection)
)
###########################################################


process.demo = cms.EDAnalyzer('MakeNTuple'
        , jetsAK8                  = cms.InputTag(JetAK8Collection)
        , jetsAK4                  = cms.InputTag(JetAK4Collection)
	, genJets		   = cms.InputTag('slimmedGenJets')
        , genJetsAK8               = cms.InputTag('slimmedGenJetsAK8')
        , muons                    = cms.InputTag(muon)
        , electrons                = cms.InputTag(electron)
        , MET                      = cms.InputTag(METCollection)
        , PFCand                   = cms.InputTag('packedPFCandidates')
        , PileupSumInfoInputTag = cms.InputTag('slimmedAddPileupInfo')
        , vertices                      = cms.InputTag('offlineSlimmedPrimaryVertices')
	, MC			   = MC
        , MC_Signal                = MC_Signal
        , DATA                     = DATA
	, totemRPLocalTrackLabel        = cms.InputTag("ctppsLocalTrackLiteProducer")
	, TriggerResults = cms.InputTag('TriggerResults', '', 'HLT')
        , MCEvent = cms.untracked.InputTag("LHCTransport")
        , psimHitTag = cms.InputTag('CTPPSSimHits','CTPPSHits')
        , recHitTag = cms.InputTag("CTPPSFastRecHits","CTPPSFastRecHits")
        , tracksPPSTag = cms.InputTag("CTPPSFastTracks","CTPPSFastTrack")
        , tagLocalTrack = cms.InputTag("totemRPLocalTrackFitter")
)


process.TFileService = cms.Service("TFileService",
                                       fileName = cms.string("out.root")
                                                                        )
if MC:
	process.p = cms.Path(process.hltFilter*process.METFilter*process.ecalBadCalibReducedMINIAODFilter*process.slimmedAK8JetsSmeared*process.slimmedAK4JetsSmeared*process.cleanPatElectrons*process.cleanPatMuons*process.cleanJets*process.cleanJetsAK8*process.genMet*process.uncorrectedMet*process.uncorrectedPatMet*process.Type1CorrForNewJEC*process.slimmedMETsNewJEC*process.shiftedMETCorrModuleForSmearedJets*process.slimmedMETsSmeared*process.filterJets*process.demo)
else:
        process.p = cms.Path(process.hltFilter*process.METFilter*process.ecalBadCalibReducedMINIAODFilter*process.cleanPatElectrons*process.cleanPatMuons*process.cleanJets*process.cleanJetsAK8*process.uncorrectedMet*process.Type1CorrForNewJEC*process.slimmedMETsNewJEC*process.filterJets*process.demo)






