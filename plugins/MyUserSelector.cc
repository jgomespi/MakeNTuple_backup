// -*- C++ -*-
//
// Package:    UserCode/MyUserSelector
// Class:      MyUserSelector
// 
/**\class MyUserSelector MyUserSelector.cc UserCode/MyUserSelector/plugins/MyUserSelector.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Joao Pedro Gomes Pinheiro
//         Created:  Fri, 09 Apr 2021 16:59:23 GMT
//
//


// system include files
#include <memory>
#include <iostream>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/VertexReco/interface/Vertex.h"

#include <vector>
#include <TMath.h>

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"

//PAT
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"

//ROOT
#include "TTree.h"
#include <TGraph.h>

#include "CondFormats/RunInfo/interface/LHCInfo.h"
#include "CondFormats/DataRecord/interface/LHCInfoRcd.h"

using namespace std;
//using namespace reco;
using namespace edm;

//
// class declaration
//

class MyUserSelector : public edm::stream::EDFilter<> {
   public:
      explicit MyUserSelector(const edm::ParameterSet&);
      ~MyUserSelector();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual bool filter(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------

         edm::Handle<pat::MuonCollection> muons;
         edm::EDGetTokenT<pat::MuonCollection> muonsToken;
         edm::Handle<pat::ElectronCollection> electrons;
         edm::EDGetTokenT<pat::ElectronCollection> electronsToken;
         edm::Handle<reco::VertexCollection> vertices;
         edm::EDGetTokenT<reco::VertexCollection> verticesToken;

	double Z_vtx, pT_min, eta_max;

	// ----------------------------------- SELECTION VARIABLES -------------------------
	
	bool vtx_quality; // |z|<15.cm and vxtIsVald
	bool lepton_pair; // 2 lepton of opposite charge, pT > 38GeV, |eta|<2.4 and loose ID.
	bool select_event; // Final condition to select the event (vtx_quality && lepton_pair)

	// ---------------------------------------------------------------------------------

         //VERTEX INFO
         int nVtx;
         bool vtx_isValid, vtx_isFake;
         double vtx_z;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
MyUserSelector::MyUserSelector(const edm::ParameterSet& iConfig):
          muonsToken (consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons")))
        , electronsToken (consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons")))
        , verticesToken (consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices")))
	, Z_vtx (iConfig.getParameter<double>("Z_vtx"))
	, pT_min(iConfig.getParameter<double>("pT_min"))
	, eta_max(iConfig.getParameter<double>("eta_max"))
{
   //now do what ever initialization is needed

}

MyUserSelector::~MyUserSelector()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
MyUserSelector::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

        iEvent.getByToken(muonsToken, muons);
        iEvent.getByToken(electronsToken, electrons);
        iEvent.getByToken(verticesToken, vertices);

        nVtx = vertices->size();

        vtx_isValid = false;
        vtx_isFake = true;
        vtx_z = -999;

        if (nVtx>0){
                vtx_isValid = vertices->at(0).isValid();
                vtx_isFake = vertices->at(0).isFake();
                vtx_z = vertices->at(0).z();
        }

	//cout << "Z_vtx: " << Z_vtx << "; pT_min = " << pT_min << "; eta_max = " << eta_max << endl;

	vtx_quality=false;
	if (vtx_isValid && abs(vtx_z)<Z_vtx) {vtx_quality=true;} // FIRST pre-SELECTION CRITERIA	
	if (vtx_quality==true){
	lepton_pair=false;
        for(size_t i = 0;i<muons->size();i++){
           if(muons->at(i).pt()>pT_min && 
	     abs(muons->at(i).eta())<eta_max && 
	     muon::isLooseMuon(muons->at(i))){
             for(size_t j = 0;j<muons->size();j++){
                if(muons->at(j).pt()>pT_min && 
		  abs(muons->at(j).eta())<eta_max && 
		  muon::isLooseMuon(muons->at(j))){
		  if (muons->at(i).charge()*muons->at(j).charge()<0){lepton_pair=true;break;}
	     	}
		if (lepton_pair == true){break;}
	     }
 	   }
	   if (lepton_pair == true){break;}
	}
	if (lepton_pair == false){
	for(size_t i = 0;i<electrons->size();i++){
           if(electrons->at(i).pt()>pT_min 
	     && abs(electrons->at(i).eta())<eta_max
	     && electrons->at(i).electronID("cutBasedElectronID-Fall17-94X-V1-loose")){
             for(size_t j = 0;j<electrons->size();j++){
                if(electrons->at(j).pt()>pT_min && 
		  abs(electrons->at(j).eta())<eta_max && 
		  electrons->at(j).electronID("cutBasedElectronID-Fall17-94X-V1-loose")){
       		  if (electrons->at(i).charge()*electrons->at(j).charge()<0){lepton_pair=true;break;}
		}
		if (lepton_pair == true){break;}
	      }
           }
	   if (lepton_pair == true){break;}
	}
	}
        if (lepton_pair == false){
        for(size_t i = 0;i<muons->size();i++){
           if(muons->at(i).pt()>pT_min &&
             abs(muons->at(i).eta())<eta_max &&
             muon::isLooseMuon(muons->at(i))){             
	      for(size_t j = 0;j<electrons->size();j++){
                if(electrons->at(j).pt()>pT_min &&
                  abs(electrons->at(j).eta())<eta_max &&
                  electrons->at(j).electronID("cutBasedElectronID-Fall17-94X-V1-loose")){
                  if (muons->at(i).charge()*electrons->at(j).charge()<0){lepton_pair=true;break;}
                }
                if (lepton_pair == true){break;}
              }
           }
           if (lepton_pair == true){break;}
        }
        }
	}
	select_event = false;
	select_event = vtx_quality && lepton_pair;

#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
   return select_event;
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
MyUserSelector::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
MyUserSelector::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
MyUserSelector::beginRun(edm::Run const&, edm::EventSetup const&)
{ 
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
MyUserSelector::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
MyUserSelector::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
MyUserSelector::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MyUserSelector::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(MyUserSelector);
