#ifndef RecoParticleFlow_PFPatProducer_PFHgcalAnalyzer_
#define RecoParticleFlow_PFPatProducer_PFHgcalAnalyzer_

// system include files
#include <memory>
#include <string>
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

#include "DataFormats/ParticleFlowReco/interface/PFSimParticle.h"
#include "DataFormats/ParticleFlowReco/interface/PFSimParticleFwd.h"

#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"


#include <TFile.h>
#include <TTree.h>
#include <TVector3.h>

#include <math.h>

#include "SimDataFormats/CaloHit/interface/PCaloHit.h"
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"


#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalDetId/interface/ESDetId.h"
#include "DataFormats/HcalRecHit/interface/HBHERecHit.h"




class PFHgcalAnalyzer : public edm::EDAnalyzer {
 public:

  explicit PFHgcalAnalyzer(const edm::ParameterSet&);

  ~PFHgcalAnalyzer();
  
  virtual void analyze(const edm::Event&, const edm::EventSetup&);

  virtual void beginRun(const edm::Run & r, const edm::EventSetup & c);

 private:
  

  /// PFCandidates in which we'll look for pile up particles 
  edm::InputTag   inputTagPFCandidates_;
  edm::InputTag   inputTagPFSimParticles_;
  
  edm::InputTag   inputTagEcalPFClusters_;
  edm::EDGetTokenT<reco::PFCandidateCollection> tokenPFCandidates_;
  edm::EDGetTokenT<reco::PFSimParticleCollection>   tokenPFSimParticles_;
  edm::EDGetTokenT<reco::PFClusterCollection>   tokenEcalPFClusters_;

  /// Min pt for charged hadrons
  double ptMin_;
  
  /// Min p for charged hadrons
  double pMin_;

  /// Min hcal raw energy for charged hadrons
  double calMin_;
  
  /// Max ecal raw energy to define a MIP
  double ecalMax_;
  
  /// Min number of pixel hits for charged hadrons
  int nPixMin_;

  /// Min number of track hits for charged hadrons
  std::vector<int> nHitMin_;
  std::vector<double> nEtaMin_;
  
  // Number of tracks after cuts
  std::vector<unsigned int> nCh;
  std::vector<unsigned int> nEv;
  
  std::string outputfile_;
  TFile *tf1;
  TTree* s;
  
  float true_,p_,ecal_,hcal_,eta_,phi_,ho_;
  float etaEcal_,phiEcal_;
  int charge_;
  std::vector<float> dr_,Eecal_,Ehcal_,pfcID_;
  size_t orun,oevt,olumiBlock,otime;
  
  edm::RunNumber_t run;
  edm::EventNumber_t evt;
  edm::LuminosityBlockNumber_t lumiBlock;
  edm::Timestamp time;


  std::vector<float> cluEcalE;
  std::vector<float> cluEcalEta;
  std::vector<float> cluEcalPhi;

  std::vector<float> distEcalTrk;


  std::vector<float> cluHcalE;
  std::vector<float> cluHcalEta;
  std::vector<float> cluHcalPhi;

  std::vector<float> distHcalTrk;
  std::vector<std::vector<float> > distHcalEcal;

  std::vector<int> pfcsID;

  std::vector<float> cluHadE;



  

  const CaloGeometry*    theCaloGeom;


  bool   verbose_;
  
  float dR(float eta1, float eta2, float phi1, float phi2 );


};

#endif
