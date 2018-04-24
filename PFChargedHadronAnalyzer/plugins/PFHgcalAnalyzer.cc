#include "PFCalibration/PFChargedHadronAnalyzer/plugins/PFHgcalAnalyzer.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementTrack.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementCluster.h"

#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFraction.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h" 
#include "DataFormats/EcalRecHit/interface/EcalUncalibratedRecHit.h" 
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h" 
#include "DataFormats/HcalRecHit/interface/HcalRecHitFwd.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "SimDataFormats/CaloHit/interface/PCaloHit.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"


#include <TROOT.h>
#include <TVector3.h>

using namespace std;
using namespace edm;
using namespace reco;

PFHgcalAnalyzer::PFHgcalAnalyzer(const edm::ParameterSet& iConfig) {

  inputTagPFCandidates_  = iConfig.getParameter<InputTag>("PFCandidates");
  tokenPFCandidates_ = consumes<reco::PFCandidateCollection>(inputTagPFCandidates_);

  // inputTagSimTracks_ = iConfig.getParameter<InputTag>("SimTracks");
  // tokenSimTracks_ = consumes<std::vector<SimTrack> >(inputTagSimTracks_);

  inputTagGenParticles_ = iConfig.getParameter<InputTag>("GenParticles");
  tokenGenParticles_ = consumes<std::vector<reco::GenParticle> >(inputTagGenParticles_);

  // Smallest track pt
  ptMin_ = iConfig.getParameter<double>("ptMin");

  // Smallest track p
  pMin_ = iConfig.getParameter<double>("pMin");

  // Smallest raw calo energy linked to the track
  calMin_ = iConfig.getParameter<double>("calMin");

  // Largest ECAL energy linked to the track to define a MIP
  ecalMax_ = iConfig.getParameter<double>("ecalMax");

  // Smallest number of pixel hits
  nPixMin_ = iConfig.getParameter<int>("nPixMin");

  // Smallest number of track hits in different eta ranges
  nHitMin_ = iConfig.getParameter< std::vector<int> > ("nHitMin");
  nEtaMin_ = iConfig.getParameter< std::vector<double> > ("nEtaMin");

  //Is minbias from simulation
  // isMBMC_ = iConfig.getUntrackedParameter<bool>("isMinBiasMC",false);

  verbose_ = 
    iConfig.getUntrackedParameter<bool>("verbose",false);

  LogDebug("PFHgcalAnalyzer")
    <<" input collection : "<<inputTagPFCandidates_ ;
   

  // The root tuple
  outputfile_ = iConfig.getParameter<std::string>("rootOutputFile"); 
  tf1 = new TFile(outputfile_.c_str(), "RECREATE");  
  s = new TTree("s"," PFCalibration");

  // s->Branch("true",&true_,"true/F");  
  // s->Branch("p",&p_,"p/F");  
  // s->Branch("ecal",&ecal_,"ecal/F");  
  // s->Branch("hcal",&hcal_,"hcal/F");  
  // s->Branch("ho",&ho_,"ho/F");  
  // s->Branch("eta",&eta_,"eta/F");  
  // s->Branch("phi",&phi_,"phi/F");
  // s->Branch("charge",&charge_,"charge/I");

  

 s->Branch("ecal_energy",&ecal_energy_,"ecal_energy/F");
 s->Branch("ecal_energy_raw",&ecal_energy_raw_,"ecal_energy_raw/F");
 s->Branch("hcal_energy",&hcal_energy_,"hcal_energy/F");
 s->Branch("hcal_energy_raw",&hcal_energy_raw_,"hcal_energy_raw/F");
 s->Branch("ho_energy",&ho_energy_,"ho_energy/F");
 s->Branch("ho_energy_raw",&ho_energy_raw_,"ho_energy_raw/F");
 s->Branch("p",&p_,"p/F");
 s->Branch("pt",&pt_,"pt/F");
 s->Branch("eta",&eta_,"eta/F");
 s->Branch("phi",&phi_,"phi/F");
 s->Branch("p_gen",&p_gen_,"p_gen/F");
 s->Branch("pt_gen",&pt_gen_,"pt_gen/F");
 s->Branch("eta_gen",&eta_gen_,"eta_gen/F");
 s->Branch("phi_gen",&phi_gen_,"phi_gen/F");
 s->Branch("e_gen",&e_gen_,"e_gen/F");




  s->Branch("run",&orun,"orun/l");
  s->Branch("evt",&oevt,"orun/l");
  s->Branch("lumiBlock",&olumiBlock,"orun/l");
  s->Branch("time",&otime,"orun/l");

   

}



PFHgcalAnalyzer::~PFHgcalAnalyzer() { 



  tf1->cd();
  s->Write();
  tf1->Write();
  tf1->Close();  


}



void 
PFHgcalAnalyzer::beginRun(const edm::Run& run, 
				  const edm::EventSetup & es) { }


void 
PFHgcalAnalyzer::analyze(const Event& iEvent, 
				 const EventSetup& iSetup) {


  orun  = (size_t)iEvent.id().run();
  oevt  = (size_t)iEvent.id().event();
  olumiBlock = (size_t)iEvent.id().luminosityBlock();
  otime = (size_t)((iEvent.time().value())>>32);


  
  // get PFCandidates
  Handle<PFCandidateCollection> pfCandidates;
  iEvent.getByToken(tokenPFCandidates_, pfCandidates);

  // Handle<std::vector<SimTrack> > simTracks;
  // iEvent.getByToken(tokenSimTracks_, simTracks);

  Handle<std::vector<reco::GenParticle> > genParticles;
  iEvent.getByToken(tokenGenParticles_, genParticles);



 ecal_energy_ = 0.0;
 ecal_energy_raw_ = 0.0;
 hcal_energy_ = 0.0;
 hcal_energy_raw_ = 0.0;
 ho_energy_ = 0.0;
 ho_energy_raw_ = 0.0;
 p_ = 0.0;
 pt_ = 0.0;
 eta_ = 0.0;
 phi_ = 0.0;
 p_gen_ = 0.0;
 pt_gen_ = 0.0;
 eta_gen_ = 0.0;
 phi_gen_ = 0.0;
 e_gen_ = 0.0;

  // count number of charged hadrons;
  std::vector<const reco::PFCandidate*> charged_hadrons;
  for (const auto& pfc : *pfCandidates)
  {
    if ( pfc.particleId() == reco::PFCandidate::h)
    {
      charged_hadrons.emplace_back(&pfc);
      // if ( fabs(pfc.eta())>1.5 && fabs(pfc.eta())<3.0)
      // {
      //   std::cout <<" raw Ecal:"<<pfc.rawEcalEnergy() <<" raw Hcal:"<<pfc.rawHcalEnergy() <<" Hcal:"<<pfc.hcalEnergy() <<" p:"<<pfc.p()<<" pt:"<<pfc.pt() <<" eta:"<<pfc.eta() <<std::endl;
      // }
    }
  }
  //select only events with exactly one charged hadron
  //change to 2 id the antiparticle is generated as well
  if (charged_hadrons.size()!=1) return; 

  for (const auto& pfc : charged_hadrons)
  {

    // Charged hadron minimum pt (the track pt, to an excellent approximation)
    if ( pfc->pt() < ptMin_ ) continue;


    // double ecalRaw = pfc->rawEcalEnergy();
    // double hcalRaw = pfc->rawHcalEnergy();
    // double hoRaw = pfc->rawHoEnergy();
    if ( pfc->rawEcalEnergy() + pfc->rawHcalEnergy() < calMin_ ) continue;


    // Find the corresponding PF block elements
    const auto& theElements = pfc->elementsInBlocks();
    if( theElements.empty() ) continue;


    const auto& track = pfc->bestTrack();


    // Characteristics of the track
    //const auto& et = dynamic_cast<const reco::PFBlockElementTrack *>( iTrack.at(0) );
    // const auto& p = pfc->p();  
    // const auto& pt = pfc->pt(); 
    // const auto& eta = pfc->eta();

    // A minimum p and pt selection
    if ( pfc->p() < pMin_ || pfc->pt() < ptMin_ ) continue;
    
    //track part

    // Count the number of valid hits (first three iteration only)
    //unsigned int nHits = et.trackRef()->found();
    unsigned int tobN = 0;
    unsigned int tecN = 0;
    unsigned int tibN = 0;
    unsigned int tidN = 0;
    unsigned int pxbN = 0;
    unsigned int pxdN = 0;
    const auto& hp = track->hitPattern();

    //selecting only certain types of tracks
    switch ( track->algo() )
    {
      case TrackBase::initialStep:
      case TrackBase::lowPtTripletStep:
      case TrackBase::pixelPairStep:
      case TrackBase::detachedTripletStep:
        tobN += hp.numberOfValidStripTOBHits();
        tecN += hp.numberOfValidStripTECHits();
        tibN += hp.numberOfValidStripTIBHits();
        tidN += hp.numberOfValidStripTIDHits();
        pxbN += hp.numberOfValidPixelBarrelHits(); 
        pxdN += hp.numberOfValidPixelEndcapHits(); 
        break;
      case TrackBase::mixedTripletStep:
      case TrackBase::pixelLessStep:
      case TrackBase::tobTecStep:
      case TrackBase::jetCoreRegionalStep:
      case TrackBase::muonSeededStepInOut:
      case TrackBase::muonSeededStepOutIn:
      default:
        break;
    }
    int inner = pxbN+pxdN;
    int outer = tibN+tobN+tidN+tecN;
    
    // selecting Number of pixel hits
    if ( inner < nPixMin_ ) continue;

    
    // Number of tracker hits (eta-dependent cut)
    bool trackerHitOK = false;
    double etaMin = 0.;
    for ( unsigned int ieta=0; ieta<nEtaMin_.size(); ++ieta )
    { 
      if ( fabs(pfc->eta()) < etaMin ) break;
      double etaMax = nEtaMin_[ieta];
      trackerHitOK = fabs(pfc->eta())>etaMin && fabs(pfc->eta())<etaMax && inner+outer>nHitMin_[ieta]; 
      if ( trackerHitOK ) break;
      etaMin = etaMax;
    }

    if ( !trackerHitOK ) continue;

    
    // Selects only ECAL MIPs
    //threshold for ecal energy in the charged case?
    if ( pfc->rawEcalEnergy() > ecalMax_ ) continue;


    // Fill the root-tuple
    // p_ = p;
    // ecal_ = pfc->rawEcalEnergy();
    // hcal_ = hcalRaw;
    // ho_ = hoRaw;

    // charge_ = pfc->charge();


    ecal_energy_ = pfc->ecalEnergy();
    ecal_energy_raw_ = pfc->rawEcalEnergy();
    hcal_energy_ = pfc->hcalEnergy();
    hcal_energy_raw_ = pfc->rawHcalEnergy();
    ho_energy_ = pfc->hoEnergy();
    ho_energy_raw_ = pfc->rawHoEnergy();
    p_ = pfc->p();
    pt_ = pfc->pt();
    eta_ = pfc->eta();
    phi_ = pfc->phi();
    p_gen_ = genParticles->at(0).p();
    pt_gen_ = genParticles->at(0).pt();
    eta_gen_ = genParticles->at(0).eta();
    phi_gen_ = genParticles->at(0).phi();
    e_gen_ = genParticles->at(0).energy();

    //Cluster characteristics
   
    //getting trajectory point at ecal entrance
    // std::cout<<"simtrack "<<simTracks->size()<<std::endl;
    // for (const auto& tr: *simTracks)
    // {

    //   std::cout<<tr.momentum().Pt()<<" "<<tr.type()<<std::endl;

    // }

    // for (const auto& tr: *genParticles)
    // {

    //   std::cout<<tr.pdgId()<<" "<<tr.pt()<<" "<<tr.numberOfMothers()<<std::endl;

    // }
    //reco::PFTrajectoryPoint::LayerType ecalEntrance = reco::PFTrajectoryPoint::ECALEntrance;
    // const reco::PFTrajectoryPoint& tpatecal = ((*trueParticles)[0]).extrapolatedPoint( ecalEntrance );
    // eta_ = tpatecal.positionREP().Eta();
    // phi_ = tpatecal.positionREP().Phi();
    // true_ = std::sqrt(tpatecal.momentum().Vect().Mag2());



    s->Fill();
        
    
  }
}


DEFINE_FWK_MODULE(PFHgcalAnalyzer);
