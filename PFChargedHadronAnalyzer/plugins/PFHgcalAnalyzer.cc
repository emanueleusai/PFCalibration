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
  auto charged_hadrons = std::make_unique<reco::PFCandidateCollection>();
  for (const auto& pfc : *pfCandidates)
  {
    if ( pfc.particleId() == reco::PFCandidate::h)
    {
      charged_hadrons.emplace_back(pfc);
      if ( fabs(pfc.eta())>1.5 && fabs(pfc.eta())<3.0)
      {
        std::cout<<" "<< pfc.rawEcalEnergy() <<" "<<pfc.rawHcalEnergy() <<" "<<pfc.hcalEnergy() <<" "<<pfc.pt() <<" "<<pfc.eta() <<std::endl;
      }
    }
  }
  //select only events with exactly two charges hadrons
  if (charged_hadrons.size()!=2) return;


  for (const auto& pfc : charged_hadrons)
  {

    // Charged hadron minimum pt (the track pt, to an excellent approximation)
    if ( pfc.pt() < ptMin_ ) continue;


    if ( pfc.rawEcalEnergy() + pfc.rawHcalEnergy() < calMin_ ) continue;

    // Find the corresponding PF block elements
    const auto& theElements = pfc.elementsInBlocks();
    if( theElements.empty() ) continue;
    
    const reco::PFBlockRef blockRef = theElements[0].first;
    PFBlock::LinkData linkData =  blockRef->linkData();
    const edm::OwnVector<reco::PFBlockElement>& elements = blockRef->elements();

    // Check that there is only one track in the block.
    unsigned int nTracks = 0;
    unsigned int nEcal = 0;
    unsigned int nHcal = 0;
    unsigned iTrack = 999;
    vector<unsigned> iECAL;// =999;
    vector<unsigned> iHCAL;// =999;
    for(unsigned iEle=0; iEle<elements.size(); iEle++)
    {
      // Find the tracks in the block
      PFBlockElement::Type type = elements[iEle].type();        
      switch( type )
      {
        case PFBlockElement::TRACK:
        	iTrack = iEle;
	        nTracks++;
	        break;
        case PFBlockElement::ECAL:
	        iECAL.push_back( iEle );
	        nEcal++;
	        break;
        case PFBlockElement::HCAL:
	        iHCAL.push_back( iEle );
	        nHcal++;
	        break;
        default:
	        continue;
      }

    }

    //bypass for neutrals
    //don't want to match to neutral PF clusters?
    if ( nTracks != 1 ) continue;
    nCh[4]++;

    // Characteristics of the track
    const reco::PFBlockElementTrack& et = dynamic_cast<const reco::PFBlockElementTrack &>( elements[iTrack] );
    double p = et.trackRef()->p();  
    double pt = et.trackRef()->pt(); 
    double eta = et.trackRef()->eta();
    double phi = et.trackRef()->phi();
    
    //ECAL element
    //loop over ecal components of the pf particle
    for(unsigned int ii=0;ii<nEcal;ii++)
    {
      const reco::PFBlockElementCluster& eecal =
	    dynamic_cast<const reco::PFBlockElementCluster &>( elements[ iECAL[ii] ] );
      double E_ECAL = eecal.clusterRef()->energy();  
      double eta_ECAL = eecal.clusterRef()->eta();
      double phi_ECAL = eecal.clusterRef()->phi();

      cluEcalE.push_back( E_ECAL );
      cluEcalEta.push_back( eta_ECAL );
      cluEcalPhi.push_back( phi_ECAL );
      
      double d = blockRef->dist(iTrack, iECAL[ii], linkData);	
      distEcalTrk.push_back( d );
      vector<float> tmp;
      emHitF.push_back( tmp );
      emHitE.push_back( tmp );
      emHitX.push_back( tmp );
      emHitY.push_back( tmp );
      emHitZ.push_back( tmp );

      if(isMBMC_ || isSimu)
      {
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


void PFHgcalAnalyzer::SaveSimHit(const edm::Event& iEvent,  float eta_, float phi_) {

  //Access to simHits informations
  Handle<PCaloHitContainer> h_PCaloHitsEB;
  iEvent.getByLabel("g4SimHits","EcalHitsEB", h_PCaloHitsEB);

  Handle<PCaloHitContainer> h_PCaloHitsEE;
  iEvent.getByLabel("g4SimHits","EcalHitsEE", h_PCaloHitsEE);

  Handle<PCaloHitContainer> h_PCaloHitsES;
  iEvent.getByLabel("g4SimHits","EcalHitsES", h_PCaloHitsES);
  
  Handle<PCaloHitContainer> h_PCaloHitsH;
  iEvent.getByLabel("g4SimHits","HcalHits", h_PCaloHitsH);

  //iterator
  PCaloHitContainer::const_iterator genSH;

  //match hits... dR 0.2, should contains all simHits
  
  //ECAL
  if( fabs(eta_) <1.5 ) { //barrel
    
    for(genSH = h_PCaloHitsEB->begin(); genSH != h_PCaloHitsEB->end(); genSH++) {
      // float theta = genSH->thetaAtEntry();
      // float phi = genSH->phiAtEntry();
      // float eta = Eta( theta );
      // float dr = dR( eta, eta_, phi, phi_ );
      
      // if(dr > 0.2 ) continue;
      //cout<<" ecal hit : "<<genSH->energy()<<endl;
      EcalSimHits.push_back( genSH->energy() );
    }
  }
  else {
    
    for(genSH = h_PCaloHitsEE->begin(); genSH != h_PCaloHitsEE->end(); genSH++) {
      // float theta = genSH->thetaAtEntry();
      // float phi = genSH->phiAtEntry();
      // float eta = Eta( theta );
      // float dr = dR( eta, eta_, phi, phi_ );

      // if(dr > 0.2 ) continue;
      EcalSimHits.push_back( genSH->energy() );
    }    

    for(genSH = h_PCaloHitsES->begin(); genSH != h_PCaloHitsES->end(); genSH++) {
       // float theta = genSH->thetaAtEntry();
       // float phi = genSH->phiAtEntry();
       // float eta = Eta( theta );
       // float dr = dR( eta, eta_, phi, phi_ );

       // if(dr > 0.2 ) continue;
       ESSimHits.push_back( genSH->energy() );
    }    
  }
  
  //Hcal
  float sH=0; 
  for(genSH = h_PCaloHitsH->begin(); genSH != h_PCaloHitsH->end(); genSH++) {
    // float theta = genSH->thetaAtEntry();
    // float phi = genSH->phiAtEntry();
    // float eta = Eta( theta );
    // float dr = dR( eta, eta_, phi, phi_ );
       // if(dr > 0.2 ) continue;
    sH += genSH->energy();
    //cout<<" ecal hit : "<<genSH->energy()<<"    "<<genSH->energyEM()<<"   "<<genSH->energyHad()<<"   "<<sH<<endl;
       HcalSimHits.push_back( genSH->energy() );
    }

}


float PFHgcalAnalyzer::Eta( float theta_ ) {
  if( sin(theta_/2.)==0 ) return 10000.*cos(theta_/2.);
  return -log(tan(theta_/2.0));
}


void PFHgcalAnalyzer::SaveRecHits(const edm::Event& iEvent, float eta_, float phi_) {

  //get rechits
  edm::Handle< EcalRecHitCollection > ebRecHits_h;
  edm::Handle< EcalRecHitCollection > eeRecHits_h;
  edm::Handle< EcalRecHitCollection > esRecHits_h;
 // Barrel
  iEvent.getByLabel( "ecalRecHit","EcalRecHitsEB", ebRecHits_h );
  // Endcaps
  iEvent.getByLabel( "ecalRecHit","EcalRecHitsEE", eeRecHits_h );
  // Preshower
  iEvent.getByLabel( "ecalRecHit","EcalRecHitsES", esRecHits_h );
  // Hcal
  edm::Handle< HBHERecHitCollection > hbheRecHits_h;
  iEvent.getByLabel( "hbhereco","", hbheRecHits_h );
  

  for( size_t ii =0; ii < ebRecHits_h->size(); ++ii )
    {
      EcalRecHitRef recHitRef( ebRecHits_h, ii );
      EBDetId id = recHitRef->id();

      const  GlobalPoint & rhPos = theCaloGeom->getPosition( id );
      float eta = rhPos.eta();
      float phi = rhPos.phi();
      float dr = dR( eta, eta_, phi, phi_ );
      if(dr > 0.1 ) continue;
      //cout<<"EB : "<<dr<<"  "<<recHitRef->energy()<<endl;
      EcalRecHits.push_back( recHitRef->energy() );
      EcalRecHitsDr.push_back( dr );
    }

  for( size_t ii =0; ii < eeRecHits_h->size(); ++ii )
    {
      EcalRecHitRef recHitRef( eeRecHits_h, ii );
      EEDetId id = recHitRef->id();

      const  GlobalPoint & rhPos = theCaloGeom->getPosition( id );
      float eta = rhPos.eta();
      float phi = rhPos.phi();
      float dr = dR( eta, eta_, phi, phi_ );
      if(dr > 0.1 ) continue;
      EcalRecHits.push_back( recHitRef->energy() );
      EcalRecHitsDr.push_back( dr );
    }


  for( size_t ii =0; ii < hbheRecHits_h->size(); ++ii )
    {
      HBHERecHitRef recHitRef( hbheRecHits_h, ii );
      HcalDetId id = recHitRef->id();

      const  GlobalPoint & rhPos = theCaloGeom->getPosition( id );
      float eta = rhPos.eta();
      float phi = rhPos.phi();
      float dr = dR( eta, eta_, phi, phi_ );
      if(dr > 0.15 ) continue;
      //cout<<"Hcal : "<<dr<<"  "<<recHitRef->energy()<<endl;
      HcalRecHits.push_back( recHitRef->energy() );
      HcalRecHitsDr.push_back( dr );
    }

  

}

float 
PFHgcalAnalyzer::phi( float x, float y ) {
  float phi_ =atan2(y, x);
  return (phi_>=0) ?  phi_ : phi_ + 2*3.141592;
}

float 
PFHgcalAnalyzer::dPhi( float phi1, float phi2 )
{
  float phi1_= phi( cos(phi1), sin(phi1) );
  float phi2_= phi( cos(phi2), sin(phi2) );
  float dphi_= phi1_-phi2_;
  if( dphi_> 3.141592 ) dphi_-=2*3.141592;
  if( dphi_<-3.141592 ) dphi_+=2*3.141592;
  return dphi_;
}

DEFINE_FWK_MODULE(PFHgcalAnalyzer);
