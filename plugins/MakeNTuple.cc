// -*- C++ -*-
//
// Package:    MakeNTuple/MakeNTuple
// Class:      MakeNTuple
// 
/**\class MakeNTuple MakeNTuple.cc MakeNTuple/MakeNTuple/plugins/MakeNTuple.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Mauricio Thiel
//         Created:  Tue, 09 Apr 2019 16:31:41 GMT
//
// Modified by JOAO PEDRO (since 17/02/21)
// The backup of this analyzer is at /afs/cern.ch/user/j/jgomespi/private/workspace/master_thesis/SemiLep_Studies/backup


// system include files
#include <memory>
#include <iostream> 
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "TRandom.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "TLorentzVector.h" 
#include <vector>
#include <TMath.h>

//HLT
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"


//PAT
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

//ROOT
#include "TTree.h" 
#include <TGraph.h>

// P P S
/*
#include "DataFormats/CTPPSReco/interface/CTPPSLocalTrackLite.h"
#include "FastSimDataFormats/CTPPSFastSim/interface/CTPPSFastRecHit.h"
#include "FastSimDataFormats/CTPPSFastSim/interface/CTPPSFastRecHitContainer.h"
#include "FastSimDataFormats/CTPPSFastSim/interface/CTPPSFastTrack.h"
#include "FastSimDataFormats/CTPPSFastSim/interface/CTPPSFastTrackContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "DataFormats/CTPPSReco/interface/TotemRPLocalTrack.h"
*/
#include "DataFormats/CTPPSDetId/interface/CTPPSDetId.h"
#include "DataFormats/CTPPSReco/interface/CTPPSLocalTrackLite.h"
#include "DataFormats/ProtonReco/interface/ForwardProton.h"

#include "CondFormats/RunInfo/interface/LHCInfo.h"
#include "CondFormats/DataRecord/interface/LHCInfoRcd.h"

#include "RoccoR.h"
#include "TRandom3.h"

using namespace std; 
using namespace reco;
using namespace edm;

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class MakeNTuple : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
	public:
		explicit MakeNTuple(const edm::ParameterSet&);
		~MakeNTuple();

		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


	private:
		virtual void beginJob() override;
		virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
		virtual void endJob() override;

		// ----------member data ---------------------------
		edm::Handle<pat::MuonCollection> muons;
		edm::EDGetTokenT<pat::MuonCollection> muonsToken;
		edm::Handle<pat::ElectronCollection> electrons;
		edm::EDGetTokenT<pat::ElectronCollection> electronsToken;
		edm::Handle<pat::METCollection> MET;
		edm::EDGetTokenT<pat::METCollection> MetToken;
		edm::Handle<std::vector<pat::PackedCandidate>> PFCand;
		edm::EDGetTokenT<std::vector<pat::PackedCandidate>> PFCandToken;
		edm::Handle<std::vector<PileupSummaryInfo> > PupInfo;
		edm::EDGetTokenT<std::vector<PileupSummaryInfo> > PileupSumInfoInputTag;
		edm::Handle<reco::VertexCollection> vertices;
		edm::EDGetTokenT<reco::VertexCollection> verticesToken;
		edm::EDGetTokenT<std::vector<reco::ForwardProton> > recoProtonsSingleRPToken_;
		edm::EDGetTokenT<std::vector<reco::ForwardProton> > recoProtonsMultiRPToken_;

		edm::LumiReWeighting *LumiWeights_;

		bool MC, MC_Signal, DATA;
		bool preVFP, postVFP;

		edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken_;

		// FOR GEN INFO
		edm::EDGetTokenT<std::vector<reco::GenParticle> > genToken;
		edm::Handle< std::vector<reco::GenParticle> > genHandle;
		edm::EDGetTokenT< edm::HepMCProduct > mcEventToken;
		edm::Handle< edm::HepMCProduct > EvtHandle ;
		void Get_t_and_xi(const TLorentzVector* p,double& t, double& xi);
		const double ProtonMass = 0.93827;
		const double ProtonMassSQ = pow(ProtonMass,2);
		double fBeamEnergy;
		double fBeamMomentum;
		void set_BeamEnergy(double e) {fBeamEnergy=e;fBeamMomentum = sqrt(fBeamEnergy*fBeamEnergy - ProtonMassSQ);};

		//FOR PU REWGT

		TTree* EventBranchs;
		int BX, Run, LumiSection, EventNum, xangle;
		double PUWeight, nPU;

		//VERTEX INFO
		int nVtx;
		bool vtx_isValid, vtx_isFake;
		double vtx_z;

		HLTConfigProvider hltConfig_;
		HLTPrescaleProvider hltPrescaleProvider_;

		std::vector<double> *ArmF_xi_gen        = new std::vector<double> ();
		std::vector<double> *ArmB_xi_gen        = new std::vector<double> ();
		std::vector<double> *ArmF_t_gen        = new std::vector<double> ();
		std::vector<double> *ArmB_t_gen        = new std::vector<double> ();
		std::vector<double> *ArmF_thx_gen        = new std::vector<double> ();
		std::vector<double> *ArmB_thx_gen        = new std::vector<double> ();
		std::vector<double> *ArmF_thy_gen        = new std::vector<double> ();
		std::vector<double> *ArmB_thy_gen        = new std::vector<double> ();

		std::vector<double> *ProtCand_xi		= new std::vector<double> ();
		std::vector<double> *ProtCand_t			= new std::vector<double> ();
		std::vector<double> *ProtCand_ThX		= new std::vector<double> ();
		std::vector<double> *ProtCand_ThY		= new std::vector<double> ();
		std::vector<double> *ProtCand_rpid		= new std::vector<double> ();
		std::vector<double> *ProtCand_arm       	= new std::vector<double> ();
		std::vector<double> *ProtCand_ismultirp		= new std::vector<double> ();
                std::vector<double> *ProtCand_x                = new std::vector<double> ();
                std::vector<double> *ProtCand_y                = new std::vector<double> ();
                std::vector<double> *ProtCand_xn                = new std::vector<double> ();
                std::vector<double> *ProtCand_yn                = new std::vector<double> ();
                std::vector<double> *ProtCand_xf                = new std::vector<double> ();
                std::vector<double> *ProtCand_yf                = new std::vector<double> ();

		std::vector<double> *muon_px        = new std::vector<double> ();
		std::vector<double> *muon_py        = new std::vector<double> ();
		std::vector<double> *muon_pz        = new std::vector<double> ();
		std::vector<double> *muon_pt        = new std::vector<double> ();
		std::vector<double> *muon_E         = new std::vector<double> ();
		std::vector<double> *muon_vtxZ      = new std::vector<double> ();
		std::vector<double> *muon_phi       = new std::vector<double> ();
		std::vector<double> *muon_eta       = new std::vector<double> ();
		std::vector<double> *muon_charge     = new std::vector<double> ();
		std::vector<double> *muon_PFBasedIso     = new std::vector<double> ();
		std::vector<double> *muon_TrackBasedIso     = new std::vector<double> ();
		std::vector<bool> *muon_isTightMuon     = new std::vector<bool> ();
		std::vector<bool> *muon_isMediumMuon     = new std::vector<bool> ();
		std::vector<bool> *muon_isLooseMuon     = new std::vector<bool> ();
		std::vector<bool> *muon_isHighPtMuon     = new std::vector<bool> ();
		std::vector<double> *muon_Roccor		= new std::vector<double> ();

		std::vector<double> *electron_px        = new std::vector<double> ();
		std::vector<double> *electron_py        = new std::vector<double> ();
		std::vector<double> *electron_pz        = new std::vector<double> ();
		std::vector<double> *electron_pt        = new std::vector<double> ();
		std::vector<double> *electron_E         = new std::vector<double> ();
		std::vector<double> *electron_vtxZ      = new std::vector<double> ();
		std::vector<double> *electron_phi       = new std::vector<double> ();
		std::vector<double> *electron_eta       = new std::vector<double> ();
		std::vector<double> *electron_charge       = new std::vector<double> ();
		std::vector<bool> *electron_isTightElectron     = new std::vector<bool> ();
		std::vector<bool> *electron_isMediumElectron     = new std::vector<bool> ();
		std::vector<bool> *electron_isLooseElectron     = new std::vector<bool> ();
		std::vector<bool> *electron_isVetoElectron     = new std::vector<bool> ();

                std::vector<double> *electron_ecalTrkEnergyPreCorr       = new std::vector<double> ();
		std::vector<double> *electron_ecalTrkEnergyErrPreCorr       = new std::vector<double> ();
		std::vector<double> *electron_ecalTrkEnergyPostCorr       = new std::vector<double> ();
                std::vector<double> *electron_ecalTrkEnergyErrPostCorr       = new std::vector<double> ();

		double prefweight;
		double prefweightup;
		double prefweightdown;


		double METPx, METPy, METPt, METphi;

		std::vector<double> *pfphi      = new std::vector<double> ();
		std::vector<double> *pfeta      = new std::vector<double> ();

		std::vector<int> *pffromPV      = new std::vector<int> ();
		std::vector<double> *pfdz      = new std::vector<double> ();
		std::vector<double> *pfpt      = new std::vector<double> ();

		std::vector<std::string> *HLT_name        = new std::vector<std::string> ();
		std::vector<bool> *HLT_pass        = new std::vector<bool> ();
		std::vector<int> *HLT_prescale        = new std::vector<int> ();
		std::vector<std::string> HLT_list = {"HLT_IsoMu24_v", "HLT_Ele27_WPTight_Gsf_v"};

		std::vector<double> *muon_gen_px        = new std::vector<double> ();
		std::vector<double> *muon_gen_py       = new std::vector<double> ();
		std::vector<double> *muon_gen_pz        = new std::vector<double> ();
		std::vector<double> *muon_gen_pt        = new std::vector<double> ();
                std::vector<double> *muon_gen_charge        = new std::vector<double> ();
		std::vector<double> *muon_gen_E        = new std::vector<double> ();
		std::vector<double> *muon_gen_phi        = new std::vector<double> ();
		std::vector<double> *muon_gen_eta        = new std::vector<double> ();

		std::vector<double> *neut_gen_px        = new std::vector<double> ();
		std::vector<double> *neut_gen_py       = new std::vector<double> ();
		std::vector<double> *neut_gen_pz        = new std::vector<double> ();
		std::vector<double> *neut_gen_pt        = new std::vector<double> ();
		std::vector<double> *neut_gen_E        = new std::vector<double> ();
		std::vector<double> *neut_gen_phi        = new std::vector<double> ();
		std::vector<double> *neut_gen_eta        = new std::vector<double> ();

		std::vector<double> *ele_gen_px        = new std::vector<double> ();
		std::vector<double> *ele_gen_py       = new std::vector<double> ();
		std::vector<double> *ele_gen_pz        = new std::vector<double> ();
		std::vector<double> *ele_gen_pt        = new std::vector<double> ();
                std::vector<double> *ele_gen_charge        = new std::vector<double> ();
		std::vector<double> *ele_gen_E        = new std::vector<double> ();
		std::vector<double> *ele_gen_phi        = new std::vector<double> ();
		std::vector<double> *ele_gen_eta        = new std::vector<double> ();

		std::vector<double> *Wp_gen_px        = new std::vector<double> ();
		std::vector<double> *Wp_gen_py       = new std::vector<double> ();
		std::vector<double> *Wp_gen_pz        = new std::vector<double> ();
		std::vector<double> *Wp_gen_pt        = new std::vector<double> ();
		std::vector<double> *Wp_gen_E        = new std::vector<double> ();
		std::vector<double> *Wp_gen_phi        = new std::vector<double> ();
		std::vector<double> *Wp_gen_eta        = new std::vector<double> ();

		std::vector<double> *Wm_gen_px        = new std::vector<double> ();
		std::vector<double> *Wm_gen_py       = new std::vector<double> ();
		std::vector<double> *Wm_gen_pz        = new std::vector<double> ();
		std::vector<double> *Wm_gen_pt        = new std::vector<double> ();
		std::vector<double> *Wm_gen_E        = new std::vector<double> ();
		std::vector<double> *Wm_gen_phi        = new std::vector<double> ();
		std::vector<double> *Wm_gen_eta        = new std::vector<double> ();

                edm::EDGetTokenT< double > prefweight_token;
                edm::EDGetTokenT< double > prefweightup_token;
                edm::EDGetTokenT< double > prefweightdown_token;

};

void MakeNTuple::Get_t_and_xi(const TLorentzVector* proton,double& t,double& xi) {
	set_BeamEnergy(13000/2.);
	t = 0.;
	xi = -1.;
	if (!proton) return;
	double mom    = proton->P();
	if (mom>fBeamMomentum) mom=fBeamMomentum;
	double energy = proton->E();
	double theta  = (proton->Pz()>0)?proton->Theta():TMath::Pi()-proton->Theta();
	t      = -2.*(ProtonMassSQ-fBeamEnergy*energy+fBeamMomentum*mom*cos(theta));
	xi     = (1.0-energy/fBeamEnergy);
}

MakeNTuple::MakeNTuple(const edm::ParameterSet& iConfig):
	  muonsToken (consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons")))
	, electronsToken (consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons")))
	, MetToken (consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("MET")))
	, PFCandToken (consumes<std::vector<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("PFCand")))
	, PileupSumInfoInputTag (consumes<std::vector<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("PileupSumInfoInputTag")))
	, verticesToken (consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices")))
	, recoProtonsSingleRPToken_   ( consumes<std::vector<reco::ForwardProton> >      ( iConfig.getParameter<edm::InputTag>( "ppsRecoProtonSingleRPTag" ) ) )
	, recoProtonsMultiRPToken_   ( consumes<std::vector<reco::ForwardProton> >      ( iConfig.getParameter<edm::InputTag>( "ppsRecoProtonMultiRPTag" ) ) )
	, LumiWeights_(0)
	, MC(iConfig.getParameter<bool>("MC"))
	, MC_Signal(iConfig.getParameter<bool>("MC_Signal"))
	, DATA(iConfig.getParameter<bool>("DATA"))
        , preVFP(iConfig.getParameter<bool>("preVFP_"))
        , postVFP(iConfig.getParameter<bool>("postVFP_"))
	, triggerResultsToken_(consumes<edm::TriggerResults>              (iConfig.getParameter<edm::InputTag>("TriggerResults")))
	, genToken(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genTag")))
	, mcEventToken (consumes<edm::HepMCProduct>(iConfig.getUntrackedParameter<edm::InputTag>("MCEvent",std::string(""))))
	, hltPrescaleProvider_(iConfig, consumesCollector(), *this)
	, prefweight_token (consumes< double >(edm::InputTag("prefiringweight:nonPrefiringProb")))
	, prefweightup_token  (consumes< double >(edm::InputTag("prefiringweight:nonPrefiringProbUp")))
	, prefweightdown_token (consumes< double >(edm::InputTag("prefiringweight:nonPrefiringProbDown")))
{
	//now do what ever initialization is needed
	usesResource("TFileService");
	edm::Service<TFileService> fs;
	EventBranchs                            = fs->make<TTree>( "Events","Events" );
	EventBranchs->Branch("Run", &Run, "Run/I");
	EventBranchs->Branch("LumiSection", &LumiSection, "LumiSection/I");
	EventBranchs->Branch("BX", &BX, "BX/I");
	EventBranchs->Branch("xangle", &xangle, "xangle/I");

	EventBranchs->Branch("EventNum", &EventNum, "EventNum/I");
	EventBranchs->Branch("nVtx", &nVtx, "nVtx/I");

	EventBranchs->Branch("PUWeight",&PUWeight,"PUWeight/D");
	EventBranchs->Branch("nPU",&nPU,"nPU/D");

	EventBranchs->Branch("vtx_z",&vtx_z,"vtx_z/D");
	EventBranchs->Branch("vtx_isValid",&vtx_isValid,"vtx_isValid/B");
	EventBranchs->Branch("vtx_isFake",&vtx_isFake,"vtx_isFake/B");

	EventBranchs->Branch("PUWeight",&PUWeight,"PUWeight/D");
	EventBranchs->Branch("nPU",&nPU,"nPU/D");

	const edm::InputTag& PileupSumInfoInputTags = edm::InputTag( "slimmedAddPileupInfo" ) ;

	LumiWeights_ = new edm::LumiReWeighting("PileupMC.root", "MyDataPileupHistogram.root", "input_Event/N_TrueInteractions", "pileup", PileupSumInfoInputTags);

	EventBranchs->Branch("ArmF_xi_gen","std::vector<double>",&ArmF_xi_gen);
	EventBranchs->Branch("ArmB_xi_gen","std::vector<double>",&ArmB_xi_gen);
	EventBranchs->Branch("ArmF_t_gen","std::vector<double>",&ArmF_t_gen);
	EventBranchs->Branch("ArmB_t_gen","std::vector<double>",&ArmB_t_gen);
	EventBranchs->Branch("ArmF_thx_gen","std::vector<double>",&ArmF_thx_gen);
	EventBranchs->Branch("ArmB_thx_gen","std::vector<double>",&ArmB_thx_gen);
	EventBranchs->Branch("ArmF_thy_gen","std::vector<double>",&ArmF_thy_gen);
	EventBranchs->Branch("ArmB_thy_gen","std::vector<double>",&ArmB_thy_gen);

	EventBranchs->Branch("ProtCand_xi","std::vector<double>",&ProtCand_xi);
	EventBranchs->Branch("ProtCand_t","std::vector<double>",&ProtCand_t);
	EventBranchs->Branch("ProtCand_ThX","std::vector<double>",&ProtCand_ThX);
	EventBranchs->Branch("ProtCand_ThY","std::vector<double>",&ProtCand_ThY);
	EventBranchs->Branch("ProtCand_rpid","std::vector<double>",&ProtCand_rpid);
	EventBranchs->Branch("ProtCand_arm","std::vector<double>",&ProtCand_arm);
	EventBranchs->Branch("ProtCand_ismultirp","std::vector<double>",&ProtCand_ismultirp);
        EventBranchs->Branch("ProtCand_x","std::vector<double>",&ProtCand_x);
        EventBranchs->Branch("ProtCand_y","std::vector<double>",&ProtCand_y);
        EventBranchs->Branch("ProtCand_xn","std::vector<double>",&ProtCand_xn);
        EventBranchs->Branch("ProtCand_yn","std::vector<double>",&ProtCand_yn);
        EventBranchs->Branch("ProtCand_xf","std::vector<double>",&ProtCand_xf);
        EventBranchs->Branch("ProtCand_yf","std::vector<double>",&ProtCand_yf);

	EventBranchs->Branch("muon_px","std::vector<double>",&muon_px);
	EventBranchs->Branch("muon_py","std::vector<double>",&muon_py);
	EventBranchs->Branch("muon_pz","std::vector<double>",&muon_pz);
	EventBranchs->Branch("muon_pt","std::vector<double>",&muon_pt);
	EventBranchs->Branch("muon_E","std::vector<double>",&muon_E);
	EventBranchs->Branch("muon_vtxZ","std::vector<double>",&muon_vtxZ);
	EventBranchs->Branch("muon_phi","std::vector<double>",&muon_phi);
	EventBranchs->Branch("muon_eta","std::vector<double>",&muon_eta);
	EventBranchs->Branch("muon_charge","std::vector<double>",&muon_charge);
	EventBranchs->Branch("muon_PFBasedIso","std::vector<double>",&muon_PFBasedIso);
	EventBranchs->Branch("muon_TrackBasedIso","std::vector<double>",&muon_TrackBasedIso);
	EventBranchs->Branch("muon_isTightMuon","std::vector<bool>",&muon_isTightMuon);
	EventBranchs->Branch("muon_isMediumMuon","std::vector<bool>",&muon_isMediumMuon);
	EventBranchs->Branch("muon_isLooseMuon","std::vector<bool>",&muon_isLooseMuon);
	EventBranchs->Branch("muon_isHighPtMuon","std::vector<bool>",&muon_isHighPtMuon);
        EventBranchs->Branch("muon_Roccor","std::vector<double>",&muon_Roccor);

	EventBranchs->Branch("electron_px","std::vector<double>",&electron_px);
	EventBranchs->Branch("electron_py","std::vector<double>",&electron_py);
	EventBranchs->Branch("electron_pz","std::vector<double>",&electron_pz);
	EventBranchs->Branch("electron_pt","std::vector<double>",&electron_pt);
	EventBranchs->Branch("electron_E","std::vector<double>",&electron_E);
	EventBranchs->Branch("electron_vtxZ","std::vector<double>",&electron_vtxZ);
	EventBranchs->Branch("electron_phi","std::vector<double>",&electron_phi);
	EventBranchs->Branch("electron_eta","std::vector<double>",&electron_eta);
	EventBranchs->Branch("electron_charge","std::vector<double>",&electron_charge);
	EventBranchs->Branch("electron_isTightElectron","std::vector<bool>",&electron_isTightElectron);
	EventBranchs->Branch("electron_isMediumElectron","std::vector<bool>",&electron_isMediumElectron);
	EventBranchs->Branch("electron_isLooseElectron","std::vector<bool>",&electron_isLooseElectron);
	EventBranchs->Branch("electron_isVetoElectron","std::vector<bool>",&electron_isVetoElectron);

        EventBranchs->Branch("electron_ecalTrkEnergyPreCorr","std::vector<double>",&electron_ecalTrkEnergyPreCorr);
        EventBranchs->Branch("electron_ecalTrkEnergyErrPreCorr","std::vector<double>",&electron_ecalTrkEnergyErrPreCorr);
        EventBranchs->Branch("electron_ecalTrkEnergyPostCorr","std::vector<double>",&electron_ecalTrkEnergyPostCorr);
        EventBranchs->Branch("electron_ecalTrkEnergyErrPostCorr","std::vector<double>",&electron_ecalTrkEnergyErrPostCorr);

	EventBranchs->Branch("prefweight",&prefweight, "prefweight/D");
        EventBranchs->Branch("prefweightup",&prefweightup, "prefweightup/D");
        EventBranchs->Branch("prefweightdown",&prefweightdown, "prefweightdown/D");
	
	EventBranchs->Branch("pfphi","std::vector<double>",&pfphi);
	EventBranchs->Branch("pfeta","std::vector<double>",&pfeta);

	EventBranchs->Branch("pffromPV","std::vector<int>",&pffromPV);
	EventBranchs->Branch("pfdz","std::vector<double>",&pfdz);
	EventBranchs->Branch("pfpt","std::vector<double>",&pfpt);

	EventBranchs->Branch("METPx", &METPx, "METPx/D");
	EventBranchs->Branch("METPy", &METPy, "METPy/D");
	EventBranchs->Branch("METPt", &METPt, "METPt/D");
	EventBranchs->Branch("METphi", &METphi, "METphi/D");

	EventBranchs->Branch("HLT_name","std::vector<std::string>",&HLT_name);
	EventBranchs->Branch("HLT_pass","std::vector<bool>",&HLT_pass);
	EventBranchs->Branch("HLT_prescale","std::vector<int>",&HLT_prescale);

	EventBranchs->Branch("muon_gen_px","std::vector<double>",&muon_gen_px);
	EventBranchs->Branch("muon_gen_py","std::vector<double>",&muon_gen_py);
	EventBranchs->Branch("muon_gen_pz","std::vector<double>",&muon_gen_pz);
	EventBranchs->Branch("muon_gen_pt","std::vector<double>",&muon_gen_pt);
        EventBranchs->Branch("muon_gen_charge","std::vector<double>",&muon_gen_charge);
	EventBranchs->Branch("muon_gen_E","std::vector<double>",&muon_gen_E);
	EventBranchs->Branch("muon_gen_phi","std::vector<double>",&muon_gen_phi);
	EventBranchs->Branch("muon_gen_eta","std::vector<double>",&muon_gen_eta);

	EventBranchs->Branch("neut_gen_px","std::vector<double>",&neut_gen_px);
	EventBranchs->Branch("neut_gen_py","std::vector<double>",&neut_gen_py);
	EventBranchs->Branch("neut_gen_pz","std::vector<double>",&neut_gen_pz);
	EventBranchs->Branch("neut_gen_pt","std::vector<double>",&neut_gen_pt);
	EventBranchs->Branch("neut_gen_E","std::vector<double>",&neut_gen_E);
	EventBranchs->Branch("neut_gen_phi","std::vector<double>",&neut_gen_phi);
	EventBranchs->Branch("neut_gen_eta","std::vector<double>",&neut_gen_eta);

	EventBranchs->Branch("ele_gen_px","std::vector<double>",&ele_gen_px);
	EventBranchs->Branch("ele_gen_py","std::vector<double>",&ele_gen_py);
	EventBranchs->Branch("ele_gen_pz","std::vector<double>",&ele_gen_pz);
	EventBranchs->Branch("ele_gen_pt","std::vector<double>",&ele_gen_pt);
        EventBranchs->Branch("ele_gen_charge","std::vector<double>",&ele_gen_charge);
	EventBranchs->Branch("ele_gen_E","std::vector<double>",&ele_gen_E);
	EventBranchs->Branch("ele_gen_phi","std::vector<double>",&ele_gen_phi);
	EventBranchs->Branch("ele_gen_eta","std::vector<double>",&ele_gen_eta);

	EventBranchs->Branch("Wp_gen_px","std::vector<double>",&Wp_gen_px);
	EventBranchs->Branch("Wp_gen_py","std::vector<double>",&Wp_gen_py);
	EventBranchs->Branch("Wp_gen_pz","std::vector<double>",&Wp_gen_pz);
	EventBranchs->Branch("Wp_gen_pt","std::vector<double>",&Wp_gen_pt);
	EventBranchs->Branch("Wp_gen_E","std::vector<double>",&Wp_gen_E);
	EventBranchs->Branch("Wp_gen_phi","std::vector<double>",&Wp_gen_phi);
	EventBranchs->Branch("Wp_gen_eta","std::vector<double>",&Wp_gen_eta);

	EventBranchs->Branch("Wm_gen_px","std::vector<double>",&Wm_gen_px);
	EventBranchs->Branch("Wm_gen_py","std::vector<double>",&Wm_gen_py);
	EventBranchs->Branch("Wm_gen_pz","std::vector<double>",&Wm_gen_pz);
	EventBranchs->Branch("Wm_gen_pt","std::vector<double>",&Wm_gen_pt);
	EventBranchs->Branch("Wm_gen_E","std::vector<double>",&Wm_gen_E);
	EventBranchs->Branch("Wm_gen_phi","std::vector<double>",&Wm_gen_phi);
	EventBranchs->Branch("Wm_gen_eta","std::vector<double>",&Wm_gen_eta);

}


MakeNTuple::~MakeNTuple()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
	void
MakeNTuple::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	ArmF_xi_gen->clear();
	ArmB_xi_gen->clear();
	ArmF_t_gen->clear();
	ArmB_t_gen->clear();
	ArmF_thx_gen->clear();
	ArmB_thx_gen->clear();
	ArmF_thy_gen->clear();
	ArmB_thy_gen->clear();

	ProtCand_xi->clear();
	ProtCand_t->clear();
	ProtCand_ThX->clear();
	ProtCand_ThY->clear();
	ProtCand_rpid->clear();
	ProtCand_arm->clear();
	ProtCand_ismultirp->clear();
        ProtCand_x->clear();
        ProtCand_y->clear();
        ProtCand_xn->clear();
        ProtCand_yn->clear();
        ProtCand_xf->clear();
        ProtCand_yf->clear();
	
	muon_px->clear();
	muon_py->clear();
	muon_pz->clear();
	muon_pt->clear();
	muon_E->clear();
	muon_vtxZ->clear();
	muon_phi->clear();
	muon_eta->clear();
        muon_charge->clear();
	muon_PFBasedIso->clear();
	muon_TrackBasedIso->clear();
	muon_isTightMuon->clear();
	muon_isMediumMuon->clear();
	muon_isLooseMuon->clear();
	muon_isHighPtMuon->clear();
	muon_Roccor->clear();

	electron_px->clear();
	electron_py->clear();
	electron_pz->clear();
	electron_pt->clear();
	electron_E->clear();
	electron_vtxZ->clear();
	electron_phi->clear();
	electron_eta->clear();
        electron_charge->clear();
	electron_isTightElectron->clear();
	electron_isMediumElectron->clear();
	electron_isLooseElectron->clear();
	electron_isVetoElectron->clear();

	pfphi->clear();
	pfeta->clear();
	pffromPV->clear();
	pfdz->clear();
	pfpt->clear();
	HLT_name->clear();
	HLT_pass->clear();
	HLT_prescale->clear();

	muon_gen_px->clear();
	muon_gen_py->clear();
	muon_gen_pz->clear();
	muon_gen_pt->clear();
        muon_gen_charge->clear();
	muon_gen_E->clear();
	muon_gen_phi->clear();
	muon_gen_eta->clear();
	neut_gen_px->clear();
	neut_gen_py->clear();
	neut_gen_pz->clear();
	neut_gen_pt->clear();
	neut_gen_E->clear();
	neut_gen_phi->clear();
	neut_gen_eta->clear();
	ele_gen_px->clear();
	ele_gen_py->clear();
	ele_gen_pz->clear();
	ele_gen_pt->clear();
        ele_gen_charge->clear();
	ele_gen_E->clear();
	ele_gen_phi->clear();
	ele_gen_eta->clear();

	Wp_gen_px->clear();
	Wp_gen_py->clear();
	Wp_gen_pz->clear();
	Wp_gen_pt->clear();
	Wp_gen_E->clear();
	Wp_gen_phi->clear();
	Wp_gen_eta->clear();
	Wm_gen_px->clear();
	Wm_gen_py->clear();
	Wm_gen_pz->clear();
	Wm_gen_pt->clear();
	Wm_gen_E->clear();
	Wm_gen_phi->clear();
	Wm_gen_eta->clear();

	electron_ecalTrkEnergyPreCorr->clear();
        electron_ecalTrkEnergyErrPreCorr->clear();
        electron_ecalTrkEnergyPostCorr->clear();
        electron_ecalTrkEnergyErrPostCorr->clear();

	prefweight=1.; prefweightup=1.; prefweightdown=1.;
	using namespace edm;
	using namespace std;
	iEvent.getByToken(muonsToken, muons);
	iEvent.getByToken(electronsToken, electrons);
	iEvent.getByToken(MetToken, MET);
	iEvent.getByToken(PFCandToken, PFCand);
	iEvent.getByToken(verticesToken, vertices);

	edm::Handle< double > theprefweight;
	iEvent.getByToken(prefweight_token, theprefweight ) ;
	prefweight =(*theprefweight);

	edm::Handle< double > theprefweightup;
	iEvent.getByToken(prefweightup_token, theprefweightup ) ;
	prefweightup =(*theprefweightup);

	edm::Handle< double > theprefweightdown;
	iEvent.getByToken(prefweightdown_token, theprefweightdown ) ;
	prefweightdown =(*theprefweightdown);

	if(MC){
		iEvent.getByToken(PileupSumInfoInputTag, PupInfo);
	}

	//	if(MC_Signal || DATA){
	edm::Handle<vector<reco::ForwardProton>> recoMultiRPProtons;
	iEvent.getByToken(recoProtonsMultiRPToken_, recoMultiRPProtons);
	edm::Handle<vector<reco::ForwardProton>> recoSingleRPProtons;
	iEvent.getByToken(recoProtonsSingleRPToken_, recoSingleRPProtons);
	//	}

	// Event Info
	BX = iEvent.bunchCrossing();
	Run = iEvent.id().run();
	LumiSection = iEvent.luminosityBlock();
	EventNum = iEvent.id().event();

	nVtx = vertices->size();

	vtx_isValid = false;
	vtx_isFake = true;
	vtx_z = -999;

	if (nVtx>0){
		vtx_isValid = vertices->at(0).isValid();
		vtx_isFake = vertices->at(0).isFake();
		vtx_z = vertices->at(0).z();
	}

	// X A N G L E   I N F O
	/*	if(DATA){
		edm::ESHandle<LHCInfo> pSetup;
		const string label = "";
		iSetup.get<LHCInfoRcd>().get(label, pSetup);
		const LHCInfo* pInfo = pSetup.product();
		std::cout << pInfo->crossingAngle() << std::endl;
		xangle = pInfo->crossingAngle();
		}
		*/


	// F O R  P I L E U P  R E W E I G T I N G
	if(MC){
		const edm::EventBase* iEventB = dynamic_cast<const edm::EventBase*>(&iEvent);
		if (LumiWeights_) PUWeight = LumiWeights_->weight( (*iEventB) );
		std::vector<PileupSummaryInfo>::const_iterator PVI;
		int npv = -1;
		for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
			int BXs = PVI->getBunchCrossing();
			if(BXs == 0) {
				npv = PVI->getTrueNumInteractions();
				continue;
			}
		}
		nPU = npv;
	}

	// F O R  T R I G G E R  I N F O
	// The trigger information was included by Joao Pedro at signal samples
	//if(!MC_Signal){
		edm::Handle<edm::TriggerResults> hltResults;
		iEvent.getByToken(triggerResultsToken_, hltResults);
		const edm::TriggerNames& trigNames = iEvent.triggerNames(*hltResults);


		for (size_t i=0; i<trigNames.size(); i++) {
			for (size_t j =0; j<HLT_list.size(); j++) {
				if (trigNames.triggerNames().at(i).find(HLT_list.at(j)) != std::string::npos){
					HLT_name->push_back(trigNames.triggerNames().at(i));
//					std::cout << "trigNames.triggerNames().at(i): " << trigNames.triggerNames().at(i) << std::endl;
					HLT_pass->push_back(hltResults->accept(i));
					HLT_prescale->push_back(hltPrescaleProvider_.prescaleValue(iEvent, iSetup,trigNames.triggerName(i)));
				}
			}
		}
	//}



	// FOR GEN INFO	
	if (MC_Signal) {

           iEvent.getByToken( mcEventToken, EvtHandle );
           const HepMC::GenEvent* Evt = EvtHandle->GetEvent() ;

	   //iEvent.getByToken(genToken, genHandle);
	   std::vector<math::XYZTLorentzVector> protonCTPPS;
	   protonCTPPS.clear();

	   //for (auto&& genPart : *(genHandle.product())){
	   for(HepMC::GenEvent::particle_const_iterator i=Evt->particles_begin(); i != Evt->particles_end();++i) {
	      int myId = (*i)->pdg_id();
	      if(myId==2212){
	 	//const pat::GenParticle& mcProton = (*genPart);
		HepMC::FourVector momentum=(*i)->momentum();
	        const HepMC::FourVector p((*i)->momentum());
		protonCTPPS.push_back(math::XYZTLorentzVector(p.x(),p.y(),p.z(),p.t()));
		double t,xi;
		double px = momentum.x();
		double py = momentum.y();
		double pz = momentum.z();
		double thx = atan(px/pz);
		double thy = atan(py/pz);
		double e = sqrt(px*px+py*py+pz*pz+ProtonMassSQ);
		TLorentzVector* proton = new TLorentzVector(px,py,pz,e);
		
		if((*i)->status() == 1) {
		   Get_t_and_xi(proton,t,xi);
		   if (pz > 5000.0) {
		      ArmF_xi_gen->push_back(xi);
		      ArmF_t_gen->push_back(t);
		      ArmF_thx_gen->push_back(thx);
		      ArmF_thy_gen->push_back(thy);
		   }
		   if (pz < 5000.0) {
		      ArmB_xi_gen->push_back(xi);
		      ArmB_t_gen->push_back(t);
		      ArmB_thx_gen->push_back(thx);
		      ArmB_thy_gen->push_back(thy);
		   }
		 }
	       }
	    
	       int pId = abs((*i)->pdg_id());

               if ( pId == 24){
		  for ( HepMC::GenVertex::particle_iterator dau  =(*i)->end_vertex()->particles_begin(HepMC::children); dau != (*i)->end_vertex()->particles_end(HepMC::children); ++dau ) {
		     if ( abs((*dau)->pdg_id()) == 24 ){
			for ( HepMC::GenVertex::particle_iterator ddau  =(*dau)->end_vertex()->particles_begin(HepMC::children); ddau != (*dau)->end_vertex()->particles_end(HepMC::children); ++ddau ) {
			    if ((*ddau)->pdg_id()==-13 || (*ddau)->pdg_id()==-11){
				Wp_gen_px->push_back((*dau)->momentum().x());
				Wp_gen_py->push_back((*dau)->momentum().y());
				Wp_gen_pz->push_back((*dau)->momentum().z());
				Wp_gen_pt->push_back( pow((*dau)->momentum().x()*(*dau)->momentum().x()+(*dau)->momentum().y()*(*dau)->momentum().y(),0.5));
				Wp_gen_E->push_back( pow((*dau)->momentum().x()*(*dau)->momentum().x()+(*dau)->momentum().y()*(*dau)->momentum().y()+(*dau)->momentum().z()*(*dau)->momentum().z(),0.5));
				Wp_gen_phi->push_back((*dau)->momentum().phi());
				Wp_gen_eta->push_back((*dau)->momentum().eta());
			     }
			     if ((*ddau)->pdg_id()==13 || (*ddau)->pdg_id()==11){
				Wm_gen_px->push_back((*dau)->momentum().x());
				Wm_gen_py->push_back((*dau)->momentum().y());
				Wm_gen_pz->push_back((*dau)->momentum().z());
				Wm_gen_pt->push_back( pow((*dau)->momentum().x()*(*dau)->momentum().x()+(*dau)->momentum().y()*(*dau)->momentum().y(),0.5));
				Wm_gen_E->push_back( pow((*dau)->momentum().x()*(*dau)->momentum().x()+(*dau)->momentum().y()*(*dau)->momentum().y()+(*dau)->momentum().z()*(*dau)->momentum().z(),0.5));
				Wm_gen_phi->push_back((*dau)->momentum().phi());
				Wm_gen_eta->push_back((*dau)->momentum().eta());
   			     }
			     if (abs((*ddau)->pdg_id()) == 13 ){
				muon_gen_px->push_back((*ddau)->momentum().x());
				muon_gen_py->push_back((*ddau)->momentum().y());
				muon_gen_pz->push_back((*ddau)->momentum().z());
                                if ((*ddau)->pdg_id()>0){muon_gen_charge->push_back(-1);}
				if ((*ddau)->pdg_id()<0){muon_gen_charge->push_back(+1);}
				muon_gen_pt->push_back( pow((*ddau)->momentum().x()*(*ddau)->momentum().x()+(*ddau)->momentum().y()*(*ddau)->momentum().y(),0.5));
				muon_gen_E->push_back( pow((*ddau)->momentum().x()*(*ddau)->momentum().x()+(*ddau)->momentum().y()*(*ddau)->momentum().y()+(*ddau)->momentum().z()*(*ddau)->momentum().z(),0.5));
				muon_gen_phi->push_back((*ddau)->momentum().phi());
				muon_gen_eta->push_back((*ddau)->momentum().eta());
			      }
			      if (abs((*ddau)->pdg_id()) == 11 ){
				 ele_gen_px->push_back((*ddau)->momentum().x());
				 ele_gen_py->push_back((*ddau)->momentum().y());
				 ele_gen_pz->push_back((*ddau)->momentum().z());
				 if ((*ddau)->pdg_id()>0){ele_gen_charge->push_back(-1);}
                                 if ((*ddau)->pdg_id()<0){ele_gen_charge->push_back(+1);}
				 ele_gen_pt->push_back( pow((*ddau)->momentum().x()*(*ddau)->momentum().x()+(*ddau)->momentum().y()*(*ddau)->momentum().y(),0.5));
				 ele_gen_E->push_back( pow((*ddau)->momentum().x()*(*ddau)->momentum().x()+(*ddau)->momentum().y()*(*ddau)->momentum().y()+(*ddau)->momentum().z()*(*ddau)->momentum().z(),0.5));
				 ele_gen_phi->push_back((*ddau)->momentum().phi());
				 ele_gen_eta->push_back((*ddau)->momentum().eta());
			       }
			       if (abs((*ddau)->pdg_id()) == 12 || abs((*ddau)->pdg_id()) == 14){
			          neut_gen_px->push_back((*ddau)->momentum().x());
				  neut_gen_py->push_back((*ddau)->momentum().y());
				  neut_gen_pz->push_back((*ddau)->momentum().z());
				  neut_gen_pt->push_back( pow((*ddau)->momentum().x()*(*ddau)->momentum().x()+(*ddau)->momentum().y()*(*ddau)->momentum().y(),0.5));
				  neut_gen_E->push_back( pow((*ddau)->momentum().x()*(*ddau)->momentum().x()+(*ddau)->momentum().y()*(*ddau)->momentum().y()+(*ddau)->momentum().z()*(*ddau)->momentum().z(),0.5));
				  neut_gen_phi->push_back((*ddau)->momentum().phi());
				  neut_gen_eta->push_back((*ddau)->momentum().eta());
				}
		           }  
			}
		      }
	            }
	     }
	 }

        if (MC && !MC_Signal) {
	   edm::Handle<reco::GenParticleCollection> genColl; 
	   iEvent.getByToken(genToken, genColl);
	   for (auto&& genPart : *(genColl.product())) { 
 	       if (genPart.status()==1 && std::abs(genPart.pdgId()) == 13) { 
                  muon_gen_px->push_back(genPart.px());
                  muon_gen_py->push_back(genPart.py());
                  muon_gen_pz->push_back(genPart.pz());	
		  muon_gen_charge->push_back(genPart.charge());
                  muon_gen_pt->push_back(pow(genPart.px()*genPart.px()+genPart.py()*genPart.py(),0.5));
                  muon_gen_E->push_back( pow(genPart.px()*genPart.px()+genPart.py()*genPart.py()+genPart.pz()*genPart.pz(),0.5));
                  muon_gen_phi->push_back(genPart.phi());
                  muon_gen_eta->push_back(genPart.eta());
	       } 
               if (genPart.status()==1 && std::abs(genPart.pdgId()) == 11) {
                  ele_gen_px->push_back(genPart.px());
                  ele_gen_py->push_back(genPart.py());
                  ele_gen_pz->push_back(genPart.pz());
                  ele_gen_charge->push_back(genPart.charge());
                  ele_gen_pt->push_back(pow(genPart.px()*genPart.px()+genPart.py()*genPart.py(),0.5));
                  ele_gen_E->push_back( pow(genPart.px()*genPart.px()+genPart.py()*genPart.py()+genPart.pz()*genPart.pz(),0.5));
                  ele_gen_phi->push_back(genPart.phi());
                  ele_gen_eta->push_back(genPart.eta());
               }
               if (genPart.status()==1 && (std::abs(genPart.pdgId()) == 12 || std::abs(genPart.pdgId()) == 14)) {
                  neut_gen_px->push_back(genPart.px());
                  neut_gen_py->push_back(genPart.py());
                  neut_gen_pz->push_back(genPart.pz());
                  neut_gen_pt->push_back(pow(genPart.px()*genPart.px()+genPart.py()*genPart.py(),0.5));
                  neut_gen_E->push_back( pow(genPart.px()*genPart.px()+genPart.py()*genPart.py()+genPart.pz()*genPart.pz(),0.5));
                  neut_gen_phi->push_back(genPart.phi());
                  neut_gen_eta->push_back(genPart.eta());
               }
	   }    

	}



	//if (DATA || MC_Signal){
		for(size_t i = 0;i<recoMultiRPProtons->size();i++){
			if(recoMultiRPProtons->at(i).validFit()){
				CTPPSDetId rpId((*recoMultiRPProtons->at(i).contributingLocalTracks().begin())->getRPId());
				//			int decRPId = rpId.arm()*100 + rpId.station()*10 + rpId.rp();
				int armId = rpId.arm();
				ProtCand_xi->push_back(recoMultiRPProtons->at(i).xi());
				//			ProtCand_t->push_back(recoMultiRPProtons->at(i).t());
				ProtCand_ThX->push_back(recoMultiRPProtons->at(i).thetaX());
				ProtCand_ThY->push_back(recoMultiRPProtons->at(i).thetaY());
				ProtCand_rpid->push_back(-999);
				ProtCand_arm->push_back(armId);
				ProtCand_ismultirp->push_back(1);
                                ProtCand_x->push_back(-999.);
                                ProtCand_y->push_back(-999.);
                                int cont = 0;
                                for (const auto &tr : recoMultiRPProtons->at(i).contributingLocalTracks()) {
                                        if(cont==0){
                                                ProtCand_xn->push_back(tr->getX());
                                                ProtCand_yn->push_back(tr->getY());
                                        }
                                        if(cont==1){
                                                ProtCand_xf->push_back(tr->getX());
                                                ProtCand_yf->push_back(tr->getY());
                                        }
                                        cont++;
                                }	
			}
		}

		for(size_t i = 0;i<recoSingleRPProtons->size();i++){
			if(recoSingleRPProtons->at(i).validFit()){
				CTPPSDetId rpId((*recoSingleRPProtons->at(i).contributingLocalTracks().begin())->getRPId());
				int decRPId = rpId.arm()*100 + rpId.station()*10 + rpId.rp();
				//                      int armId = rpId.arm();
				ProtCand_xi->push_back(recoSingleRPProtons->at(i).xi());
				//                      ProtCand_t->push_back(recoSingleRPProtons->at(i).t());
				ProtCand_ThX->push_back(recoSingleRPProtons->at(i).thetaX());
				ProtCand_ThY->push_back(recoSingleRPProtons->at(i).thetaY());
				ProtCand_rpid->push_back(decRPId);
				ProtCand_arm->push_back(-999);
				ProtCand_ismultirp->push_back(0);
                                for (const auto &tr : recoSingleRPProtons->at(i).contributingLocalTracks()) {
                                        ProtCand_x->push_back(tr->getX());
                                        ProtCand_y->push_back(tr->getY());
                                }
                                ProtCand_xn->push_back(-999.);
                                ProtCand_yn->push_back(-999.);
                                ProtCand_xf->push_back(-999.);
                                ProtCand_yf->push_back(-999.);
			}
		}
	//}


	/*
	   if (MC && !MC_Signal){
	   double strips_xi_arm0_N;
	   double strips_xi_arm0_F;
	   double strips_xi_arm1_N;
	   double strips_xi_arm1_F;
	   float a, b, c, d;


	   Int_t ncols;
	   Int_t nlines = 0;
	   FILE *fp = fopen("RPs_xi.txt","r");

	   int rand_no = rand() % 557673;

	   while (nlines<rand_no+1) {

	   ncols = fscanf(fp,"%f %f %f %f",&a, &b, &c, &d);


	   if (ncols < 0) break;

	   if (nlines==rand_no){
	   strips_xi_arm0_N = a;
	   strips_xi_arm0_F = b;
	   strips_xi_arm1_N = c;
	   strips_xi_arm1_F = d;
	   }

	   nlines++;

	   }

	   fclose(fp);
	   if (strips_xi_arm0_N>0.005) {ArmF_N_xi->push_back(strips_xi_arm0_N);}
	   if (strips_xi_arm0_F>0.005) {ArmF_F_xi->push_back(strips_xi_arm0_F);}
	   if (strips_xi_arm1_N>0.005) {ArmB_N_xi->push_back(strips_xi_arm1_N);}
	   if (strips_xi_arm1_F>0.005) {ArmB_F_xi->push_back(strips_xi_arm1_F);}

	   }
	   */

	TRandom3 rnd(1234);
	RoccoR rcA, rcB;
	//RoccoR rcA("/afs/cern.ch/user/j/jgomespi/private/workspace/master_thesis/FullLep_Studies/2016_analysis/Ntuple_build/CMSSW_10_6_20/src/MakeNTuple/MakeNTuple/plugins/RoccoR2016aUL.txt"); 
	//RoccoR rcB("/afs/cern.ch/user/j/jgomespi/private/workspace/master_thesis/FullLep_Studies/2016_analysis/Ntuple_build/CMSSW_10_6_20/src/MakeNTuple/MakeNTuple/plugins/RoccoR2016bUL.txt"); 

	//rcA.init(edm::FileInPath("/afs/cern.ch/user/j/jgomespi/private/workspace/master_thesis/FullLep_Studies/2016_analysis/Ntuple_build/CMSSW_10_6_20/src/MakeNTuple/MakeNTuple/plugins/RoccoR2016aUL.txt").fullPath());
        //rcB.init(edm::FileInPath("/afs/cern.ch/user/j/jgomespi/private/workspace/master_thesis/FullLep_Studies/2016_analysis/Ntuple_build/CMSSW_10_6_20/src/MakeNTuple/MakeNTuple/RoccoR2016bUL.txt").fullPath());
	
	rcA.init("RoccoR2016aUL.txt");
	rcB.init("RoccoR2016bUL.txt");

	for(size_t i = 0;i<muons->size();i++){
		if(muons->at(i).pt()>20 && abs(muons->at(i).eta())<2.4){
		    muon_px->push_back(muons->at(i).px());
		    muon_py->push_back(muons->at(i).py());
		    muon_pz->push_back(muons->at(i).pz());
		    muon_pt->push_back(muons->at(i).pt());
		    muon_E->push_back(muons->at(i).energy());
		    muon_vtxZ->push_back(muons->at(i).vertex().z());
		    muon_phi->push_back(muons->at(i).phi());
		    muon_eta->push_back(muons->at(i).eta());
	            muon_charge->push_back(muons->at(i).charge());
		    muon_PFBasedIso->push_back((muons->at(i).pfIsolationR04().sumChargedHadronPt + max(0., muons->at(i).pfIsolationR04().sumNeutralHadronEt + muons->at(i).pfIsolationR04().sumPhotonEt - 0.5*muons->at(i).pfIsolationR04().sumPUPt))/muons->at(i).pt());
		    muon_TrackBasedIso->push_back((muons->at(i).isolationR03().sumPt)/(muons->at(i).pt()));
		    muon_isTightMuon->push_back(muon::isTightMuon(muons->at(i), vertices->at(0)));
		    muon_isMediumMuon->push_back(muon::isMediumMuon(muons->at(i)));
		    muon_isLooseMuon->push_back(muon::isLooseMuon(muons->at(i)));
		    muon_isHighPtMuon->push_back(muon::isHighPtMuon(muons->at(i), vertices->at(0)));

		    float bestdr(9999.);
                    float bestmatch_pt(-1);
		
		    if (MC && MC_Signal) {
                        iEvent.getByToken( mcEventToken, EvtHandle );
                        const HepMC::GenEvent* Evt = EvtHandle->GetEvent() ;
 		        for(HepMC::GenEvent::particle_const_iterator j=Evt->particles_begin(); j != Evt->particles_end();++j) {
			   if (std::abs((*j)->pdg_id()) != 13 || (*j)->status() != 1) continue;
			   HepMC::FourVector momentum=(*j)->momentum();
			   const HepMC::FourVector p((*j)->momentum());
			   double px = momentum.x(),  py = momentum.y(), pz = momentum.z(), e = sqrt(px*px+py*py+pz*pz);
			   TLorentzVector* GenPart = new TLorentzVector(px,py,pz,e);					     
			   float dr = deltaR2(muons->at(i).eta(),muons->at(i).phi(), GenPart->Eta(),GenPart->Phi());
      			   if (dr < 0.1 && dr < bestdr) {
         		      bestmatch_pt = GenPart->Pt();
         		      bestdr = dr;
      			   } 
  		        }
		    }

			
                    if (MC && !MC_Signal) {
           	       edm::Handle<reco::GenParticleCollection> genColl;
           	       iEvent.getByToken(genToken, genColl);
                       for (auto&& genPart : *(genColl.product())) {
               		   if (genPart.status()!=1 || std::abs(genPart.pdgId()) != 13) continue;
                           float dr = deltaR2(muons->at(i).eta(),muons->at(i).phi(), genPart.eta(),genPart.phi());
                           if (dr < 0.1 && dr < bestdr) {
                              bestmatch_pt = genPart.pt();
                              bestdr = dr;
                           }
                        }
		    }


		    double muSF(1.);

		    if (MC && preVFP){
		       if (bestmatch_pt>0){muSF = rcA.kSpreadMC(muons->at(i).charge(), muons->at(i).pt(), muons->at(i).eta(), muons->at(i).phi(),bestmatch_pt);}
		       else if (muons->at(i).innerTrack().isNonnull()){muSF = rcA.kSmearMC(muons->at(i).charge(), muons->at(i).pt(), muons->at(i).eta(), muons->at(i).phi(),muons->at(i).innerTrack()->hitPattern().trackerLayersWithMeasurement(), gRandom->Rndm());}
		    }
                    if (MC && postVFP){
                       if (bestmatch_pt>0){muSF = rcB.kSpreadMC(muons->at(i).charge(), muons->at(i).pt(), muons->at(i).eta(), muons->at(i).phi(),bestmatch_pt);}
                       else if (muons->at(i).innerTrack().isNonnull()){muSF = rcB.kSmearMC(muons->at(i).charge(), muons->at(i).pt(), muons->at(i).eta(), muons->at(i).phi(),muons->at(i).innerTrack()->hitPattern().trackerLayersWithMeasurement(), gRandom->Rndm());}
                    } 
                    if (DATA){
                       if (preVFP) {muSF = rcA.kScaleDT(muons->at(i).charge(), muons->at(i).pt(), muons->at(i).eta(), muons->at(i).phi());}
                       if (postVFP){muSF = rcB.kScaleDT(muons->at(i).charge(), muons->at(i).pt(), muons->at(i).eta(), muons->at(i).phi());}
                    }
 		    muon_Roccor->push_back(muSF);
		    cout << "muSF = " << muSF << endl;
		
		}
		if (muons->at(i).innerTrack().isNonnull()){cout << "trackerLayersWithMeasurement = " << muons->at(i).innerTrack()->hitPattern().trackerLayersWithMeasurement()<<endl;}
	}

	for(size_t i = 0;i<electrons->size();i++){
		if(electrons->at(i).pt()>20 && abs(electrons->at(i).eta())<2.4){
			electron_px->push_back(electrons->at(i).px());
			electron_py->push_back(electrons->at(i).py());
			electron_pz->push_back(electrons->at(i).pz());
			electron_pt->push_back(electrons->at(i).pt());
			electron_E->push_back(electrons->at(i).energy());
			electron_vtxZ->push_back(electrons->at(i).vertex().z());
			electron_phi->push_back(electrons->at(i).phi());
			electron_eta->push_back(electrons->at(i).eta());
			electron_charge->push_back(electrons->at(i).charge());	
			electron_isTightElectron->push_back(electrons->at(i).electronID("cutBasedElectronID-Fall17-94X-V1-tight"));
                        electron_isMediumElectron->push_back(electrons->at(i).electronID("cutBasedElectronID-Fall17-94X-V1-medium"));
                        electron_isLooseElectron->push_back(electrons->at(i).electronID("cutBasedElectronID-Fall17-94X-V1-loose"));
                        electron_isVetoElectron->push_back(electrons->at(i).electronID("cutBasedElectronID-Fall17-94X-V1-veto"));
			electron_ecalTrkEnergyPreCorr->push_back(electrons->at(i).userFloat("ecalTrkEnergyPreCorr"));
                        electron_ecalTrkEnergyErrPreCorr->push_back(electrons->at(i).userFloat("ecalTrkEnergyErrPreCorr"));
                        electron_ecalTrkEnergyPostCorr->push_back(electrons->at(i).userFloat("ecalTrkEnergyPostCorr"));
                        electron_ecalTrkEnergyErrPostCorr->push_back(electrons->at(i).userFloat("ecalTrkEnergyErrPostCorr"));
		}
	}

	METPx = (MET->front()).px();
	METPy = (MET->front()).py();
	METPt = (MET->front()).pt();
	METphi = (MET->front()).phi();

	int npfVtx = 0;
	for (size_t pf = 0; pf < PFCand->size(); pf++) {
		if (abs(PFCand->at(pf).pdgId()) == 211 || abs(PFCand->at(pf).pdgId()) == 11 || abs(PFCand->at(pf).pdgId()) == 13){
//			if (PFCand->at(pf).fromPV(0)>1) {
//				if (abs(PFCand->at(pf).dz())<0.15){
					npfVtx++;
					pfphi->push_back(PFCand->at(pf).phiAtVtx());
					pfeta->push_back(PFCand->at(pf).eta());
					pffromPV->push_back(PFCand->at(pf).fromPV(0));
					pfdz->push_back(PFCand->at(pf).dz());
					pfpt->push_back(PFCand->at(pf).pt());
//				}
//			}
		}
	}

	EventBranchs->Fill();

#ifdef THIS_IS_AN_EVENT_EXAMPLE
	Handle<ExampleData> pIn;
	iEvent.getByLabel("example",pIn);
#endif

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
	ESHandle<SetupData> pSetup;
	iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
	void 
MakeNTuple::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
	void 
MakeNTuple::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MakeNTuple::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MakeNTuple);
