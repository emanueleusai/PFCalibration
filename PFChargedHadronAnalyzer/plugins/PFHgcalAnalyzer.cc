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

  inputTagPFCandidates_ 
    = iConfig.getParameter<InputTag>("PFCandidates");
  tokenPFCandidates_ = consumes<reco::PFCandidateCollection>(inputTagPFCandidates_);

  inputTagEcalPFClusters_ 
    = iConfig.getParameter<InputTag>("EcalPFClusters");
  tokenEcalPFClusters_ = consumes<reco::PFClusterCollection>(inputTagEcalPFClusters_);

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
  isMBMC_ = iConfig.getUntrackedParameter<bool>("isMinBiasMC",false);

  verbose_ = 
    iConfig.getUntrackedParameter<bool>("verbose",false);

  LogDebug("PFHgcalAnalyzer")
    <<" input collection : "<<inputTagPFCandidates_ ;
   

  // The root tuple
  outputfile_ = iConfig.getParameter<std::string>("rootOutputFile"); 
  tf1 = new TFile(outputfile_.c_str(), "RECREATE");  
  s = new TTree("s"," PFCalibration");

  s->Branch("true",&true_,"true/F");  
  s->Branch("p",&p_,"p/F");  
  s->Branch("ecal",&ecal_,"ecal/F");  
  s->Branch("hcal",&hcal_,"hcal/F");  
  s->Branch("ho",&ho_,"ho/F");  
  s->Branch("eta",&eta_,"eta/F");  
  s->Branch("phi",&phi_,"phi/F");
  s->Branch("charge",&charge_,"charge/I");

  s->Branch("dr",&dr_);  //spandey Apr_27 dR
  s->Branch("Eecal",&Eecal_);  //spandey Apr_27 dR
  s->Branch("Ehcal",&Ehcal_);  //spandey Apr_27 dR
  s->Branch("pfcID",&pfcID_);  //spandey Apr_27 dR

  s->Branch("pfcs",&pfcsID);

  


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
  
  LogDebug("PFHgcalAnalyzer")<<"START event: "<<iEvent.id().event()
			 <<" in run "<<iEvent.id().run()<<endl;
  

   edm::ESHandle<CaloGeometry> pCalo;
   iSetup.get<CaloGeometryRecord>().get( pCalo );
   theCaloGeom = pCalo.product();



  run  = iEvent.id().run();
  evt  = iEvent.id().event();
  lumiBlock = iEvent.id().luminosityBlock();
  time = iEvent.time();

  orun = (size_t)run;
  oevt = (size_t)evt;
  olumiBlock = (size_t)lumiBlock;
  otime = (size_t)((iEvent.time().value())>>32);

  
  // get PFCandidates
  Handle<PFCandidateCollection> pfCandidates;
  iEvent.getByToken(tokenPFCandidates_, pfCandidates);

  //get Ecal PFClusters
  Handle<reco::PFClusterCollection> pfClustersEcal;
  iEvent.getByToken(tokenEcalPFClusters_, pfClustersEcal);

  // Handle<PFSimParticleCollection> trueParticles;
  // //FIXME
  // // bool isMBMC=true;
  // bool isSimu = iEvent.getByToken(tokenPFSimParticles_, trueParticles);

  //simHits
  EcalSimHits.clear();
  ESSimHits.clear();
  HcalSimHits.clear();
  
  //recHits
  EcalRecHits.clear();
  ESRecHits.clear();
  HcalRecHits.clear();
  EcalRecHitsDr.clear();
  ESRecHitsDr.clear();
  HcalRecHitsDr.clear();

  pfcsID.clear();

  charge_=0;
  dr_.clear();
  Eecal_.clear();
  Ehcal_.clear();  
  pfcID_.clear();
  p_ = 0.;
  charge_=0;
  ecal_ = 0.;
  hcal_ = 0.;

  // count number of charged hadrons;
  // auto charged_hadrons = std::make_unique<reco::PFCandidateCollection>();
  std::vector<const reco::PFCandidate*> charged_hadrons;
  for (const auto& pfc : *pfCandidates)
  {
    if ( pfc.particleId() == reco::PFCandidate::h)
    {
      charged_hadrons.emplace_back(&pfc);
      if ( fabs(pfc.eta())>1.5 && fabs(pfc.eta())<3.0)
      {
        std::cout<<" "<< pfc.rawEcalEnergy() <<" "<<pfc.rawHcalEnergy() <<" "<<pfc.hcalEnergy() <<" "<<pfc.pt() <<" "<<pfc.eta() <<std::endl;
      }
    }
  }
  //select only events with exactly one charged hadron
  //change to 2 id the antiparticle is generated as well
  if (charged_hadrons.size()!=1) return; 

  for (const auto& pfc : charged_hadrons)
  {

    // Charged hadron minimum pt (the track pt, to an excellent approximation)
    if ( pfc->pt() < ptMin_ ) continue;


    if ( pfc->rawEcalEnergy() + pfc->rawHcalEnergy() < calMin_ ) continue;

    // Find the corresponding PF block elements
    const auto& theElements = pfc->elementsInBlocks();
    if( theElements.empty() ) continue;
    
    const auto&  blockRef = theElements[0].first;
    const auto&  linkData =  blockRef->linkData();
    const auto&  elements = blockRef->elements();

    // Check that there is only one track in the block.

    std::vector<const reco::PFBlockElement*> iTrack;
    std::vector<const reco::PFBlockElement*> iECAL;
    std::vector<const reco::PFBlockElement*> iHCAL;

    for(const auto& iEle: elements)
    {
      // Find the tracks in the block
      const auto& type = iEle.type();        
      switch( type )
      {
        case PFBlockElement::TRACK:
        	iTrack.emplace_back(&iEle);
	        break;
        case PFBlockElement::ECAL:
	        iECAL.emplace_back(&iEle);
	        break;
        case PFBlockElement::HCAL:
	        iHCAL.emplace_back(&iEle);
	        break;
        default:
	        continue;
      }

    }

    //bypass for neutrals
    //don't want to match to neutral PF clusters?
    if ( iTrack.size() != 1 ) continue;

    // Characteristics of the track
    const auto& et = dynamic_cast<const reco::PFBlockElementTrack *>( iTrack.at(0) );
    const auto& p = et->trackRef()->p();  
    const auto& pt = et->trackRef()->pt(); 
    const auto& eta = et->trackRef()->eta();
    const auto& phi = et->trackRef()->phi();
    
    //ECAL element
    //loop over ecal components of the pf particle
    //for(unsigned int ii=0;ii<nEcal;ii++)
    for(const auto& ii: iECAL)
    {
      const auto& eecal = dynamic_cast<const reco::PFBlockElementCluster *>( ii );
      const auto& E_ECAL = eecal->clusterRef()->energy();  
      const auto& eta_ECAL = eecal->clusterRef()->eta();
      const auto& phi_ECAL = eecal->clusterRef()->phi();

      cluEcalE.push_back( E_ECAL );
      cluEcalEta.push_back( eta_ECAL );
      cluEcalPhi.push_back( phi_ECAL );
      
      const auto& d = blockRef->dist(*iTrack.at(0), *ii, linkData);	
      distEcalTrk.push_back( d );
      vector<float> tmp;
      emHitF.push_back( tmp );
      emHitE.push_back( tmp );
      emHitX.push_back( tmp );
      emHitY.push_back( tmp );
      emHitZ.push_back( tmp );

	    const std::vector< reco::PFRecHitFraction > erh=eecal.clusterRef()->recHitFractions();
         //loop over rechit fractions
	       for(unsigned int ieh=0;ieh<erh.size();ieh++)
         {
	         emHitF[ii].push_back( erh[ieh].fraction() );
	         emHitE[ii].push_back(  erh[ieh].recHitRef()->energy() );
	         bool isEB= erh[ieh].recHitRef()->layer()==-1;
	         emHitX[ii].push_back( isEB?erh[ieh].recHitRef()->position().eta() :erh[ieh].recHitRef()->position().x() );
	         emHitY[ii].push_back( isEB?erh[ieh].recHitRef()->position().phi() :erh[ieh].recHitRef()->position().y() );
	         emHitZ[ii].push_back( isEB?0:erh[ieh].recHitRef()->position().z() );
	       }
    }//ecal element loop

    //HCAL element
    for(unsigned int ii=0;ii<nHcal;ii++)
    {
	    const reco::PFBlockElementCluster& ehcal =
	    dynamic_cast<const reco::PFBlockElementCluster &>( elements[iHCAL[ii] ] );
	    double E_HCAL = ehcal.clusterRef()->energy();  
	    double eta_HCAL = ehcal.clusterRef()->eta();
	    double phi_HCAL = ehcal.clusterRef()->phi();

	    cluHcalE.push_back( E_HCAL );
	    cluHcalEta.push_back( eta_HCAL );
	    cluHcalPhi.push_back( phi_HCAL );

	    double d = blockRef->dist(iTrack, iHCAL[ii], linkData);	
	    distHcalTrk.push_back( d );

	    //ECAL-HCAL distance
	    vector<float> tmp;
	    distHcalEcal.push_back(tmp);
      //loop over ecal elements and store distance
	    for(unsigned int ij=0;ij<nEcal;ij++)
      {
	       d = blockRef->dist(iECAL[ij], iHCAL[ii], linkData);	
	       distHcalEcal[ii].push_back( d );
	    }
	    hadHitF.push_back( tmp );
	    hadHitE.push_back( tmp );
	    hadHitX.push_back( tmp );
	    hadHitY.push_back( tmp );
	    hadHitZ.push_back( tmp );

	    if(isMBMC_ || isSimu)
      {
	       const std::vector< reco::PFRecHitFraction > erh=ehcal.clusterRef()->recHitFractions();
         //loop on fractions
	       for(unsigned int ieh=0;ieh<erh.size();ieh++)
         {
	         hadHitF[ii].push_back( erh[ieh].fraction() );
	         hadHitE[ii].push_back(  erh[ieh].recHitRef()->energy() );

	         bool isHB= erh[ieh].recHitRef()->layer()==1;
	         hadHitX[ii].push_back( isHB?erh[ieh].recHitRef()->position().eta() :erh[ieh].recHitRef()->position().x() );
	         hadHitY[ii].push_back( isHB?erh[ieh].recHitRef()->position().phi() :erh[ieh].recHitRef()->position().y() );
	         hadHitZ[ii].push_back( isHB?0:erh[ieh].recHitRef()->position().z() );
	       }
	     }

    }//hcal element loop

    
    // A minimum p and pt selection
      if ( p < pMin_ || pt < ptMin_ ) continue;
      nCh[5]++;
    
    //track part

    // Count the number of valid hits (first three iteration only)
    //unsigned int nHits = et.trackRef()->found();
    unsigned int tobN = 0;
    unsigned int tecN = 0;
    unsigned int tibN = 0;
    unsigned int tidN = 0;
    unsigned int pxbN = 0;
    unsigned int pxdN = 0;
    const reco::HitPattern& hp = et.trackRef()->hitPattern();

    //selecting only certain types of tracks
    switch ( et.trackRef()->algo() )
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
    nCh[6]++;
    
    // Number of tracker hits (eta-dependent cut)
    bool trackerHitOK = false;
    double etaMin = 0.;
    for ( unsigned int ieta=0; ieta<nEtaMin_.size(); ++ieta )
    { 
      if ( fabs(eta) < etaMin ) break;
      double etaMax = nEtaMin_[ieta];
      trackerHitOK = fabs(eta)>etaMin && fabs(eta)<etaMax && inner+outer>nHitMin_[ieta]; 
      if ( trackerHitOK ) break;
      etaMin = etaMax;
    }

    if ( !trackerHitOK ) continue;
    nCh[7]++;
    
    // Selects only ECAL MIPs
    //threshold for ecal energy in the charged case?
    if ( ecalRaw > ecalMax_ ) continue;
    nCh[8]++;

    
    //extrapolate track to ECAL --> impact position
    etaEcal_ = et.positionAtECALEntrance().Eta();
    phiEcal_ = et.positionAtECALEntrance().Phi();


    // Fill the root-tuple
    p_ = p;
    ecal_ = ecalRaw;
    hcal_ = hcalRaw;
    ho_ = hoRaw;

    charge_ = pfc.charge();

    //Cluster characteristics
   
    //getting trajectory point at ecal entrance
    if( isSimu )
    {
      reco::PFTrajectoryPoint::LayerType ecalEntrance = reco::PFTrajectoryPoint::ECALEntrance;
      const reco::PFTrajectoryPoint& tpatecal = ((*trueParticles)[0]).extrapolatedPoint( ecalEntrance );
      eta_ = tpatecal.positionREP().Eta();
      phi_ = tpatecal.positionREP().Phi();
      true_ = std::sqrt(tpatecal.momentum().Vect().Mag2());
    }
    else
    {
      eta_ = eta; 
      phi_ = phi; 
      true_ = p; 
    }


    s->Fill();

    addDr.clear();
    addPdgId.clear();
    addEmE.clear();
    addHadE.clear();
    addEta.clear();
    addPhi.clear();
    
    cluEcalE.clear();
    cluEcalEta.clear();
    cluEcalPhi.clear();

    distEcalTrk.clear();

    cluHcalE.clear();
    cluHcalEta.clear();
    cluHcalPhi.clear();

    distHcalTrk.clear();
    distHcalEcal.clear();

    genDr.clear();
    genPdgId.clear();
    genE.clear();
    genEta.clear();
    genPhi.clear();
  
    emHitF.clear();
    emHitE.clear();
    emHitX.clear();
    emHitY.clear();
    emHitZ.clear();
    hadHitF.clear();
    hadHitE.clear();
    hadHitX.clear();
    hadHitY.clear();
    hadHitZ.clear();
    
    bcEcalE.clear();
    bcEcalEta.clear();
    bcEcalPhi.clear();
    
    
  }
}


float PFHgcalAnalyzer::dR(float eta1, float eta2, float phi1, float phi2 ) {

  TVector3 v1(0,0,0),v2(0,0,0);
  
  v1.SetPtEtaPhi(1, eta1, phi1);
  v2.SetPtEtaPhi(1, eta2, phi2);

  return v1.DrEtaPhi( v2 );
  
}


DEFINE_FWK_MODULE(PFHgcalAnalyzer);
