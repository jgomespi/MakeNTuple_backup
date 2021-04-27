#include <memory>
#include <iostream>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/PatCandidates/interface/Jet.h"


class FilterAK8Jet : public edm::stream::EDFilter<> {
	public:
		explicit FilterAK8Jet(const edm::ParameterSet&);


	private:
		virtual bool filter(edm::Event&, const edm::EventSetup&) override;

		edm::Handle<pat::JetCollection> jetsAK8;
		edm::EDGetTokenT<pat::JetCollection> jetsAK8Token;

};


FilterAK8Jet::FilterAK8Jet(const edm::ParameterSet& iConfig):
	jetsAK8Token (consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jetsAK8"))){}

bool FilterAK8Jet::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

using namespace edm;
using namespace std;

iEvent.getByToken(jetsAK8Token, jetsAK8);

   size_t j = 0;
   for(size_t i=0; i<jetsAK8->size(); i++) {
	if(jetsAK8->at(i).pt()>170 && abs(jetsAK8->at(i).eta())<2.4) j++;
	}

if (j<1) {return false;} else return true;

}
//define this as a plug-in
DEFINE_FWK_MODULE(FilterAK8Jet);
