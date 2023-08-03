// -*- Header -*-
//
// Package:    Analyzers/NtrkDistribution
// Class:      NtrkDistribution
// 
/**\class NtrkDistribution NtrkDistribution.h Analyzers/NtrkDistribution/interface/NtrkDistribution.h
   
 Description: [one line class summary]
 
 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Maxime Guilbaud
//         Created:  Thu, 01 Jun 2017 16:56:11 GMT
//
//


// system include files
#include <memory>
//#include <iostream>
#include <vector>

// CMSSW include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "DataFormats/HeavyIonEvent/interface/Centrality.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/HeavyIonEvent/interface/CentralityBins.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
// user include files
#include "TH1F.h"
#include "TH2D.h"
#include "TTree.h"
#include "TLorentzVector.h"

//using namespace std;

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class NtrkDistribution : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
  explicit NtrkDistribution(const edm::ParameterSet&);
  ~NtrkDistribution();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  int getEffNoffIndex();
  
 private:
 
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  // ----------member data ---------------------------
  // ## tracks ##
  // used to select what tracks to read from configuration file
  edm::EDGetTokenT<reco::TrackCollection> trackTags_;
  
  // ## vertex ##
  // used to select what vertex to read from configuration file
  edm::EDGetTokenT<reco::VertexCollection> vtxTags_; 
  
  // ## calotower ##
  // used to select what calo tower to read from configuration file
  edm::EDGetTokenT<CaloTowerCollection> caloTowersTags_; 
  
  // ## centrality ##
  // used to select what centrality collection to read from configuration file
  edm::EDGetTokenT<reco::Centrality> centralityTags_;
  // used to access centrality bins 
  edm::EDGetTokenT<int> centralityBinTags_;
  int cent_;
  
  //Gen Particles
  edm::EDGetTokenT<reco::GenParticleCollection> trackTagsgen_;
  
  // ## multiplicity selection (Noff)
  int noffmin_;          //minimum multiplicity of an event to be considered
  int noffmax_;          //maximum multiplicity of an event to be considered
  double ptnoffmin_;     //minimum pt cut to compute Noff
  double ptnoffmax_;     //maximum pt cut to compute Noff
  double dzdzerrornoff_; //cut on dz/dzerror of the tracks to compute Noff 
  double d0d0errornoff_; //cut on d0/d0error of the tracks to compute Noff
  double pterrorptnoff_; //cut on pterror/pt of the tracks to compute Noff
  int noff_;             //ntrk offline value for a given event
  int noffcorr_;         //ntrk offline corrected value for a given event

  //Defining Histograms
  TH1F* hPtTrk_;
      
  
  //Defining Trees
  // std::vector<double> Xbestvtxtree;
  //std::vector<double> Ybestvtxtree;
  //std::vector<double> Zbestvtxtree;
  //std::vector<double> Rhobestvtxtree;
  std::vector<double>dzvtx_;
  std::vector<double>dxyvtx_;
  std::vector<double>dzerror_;
  std::vector<double>dxyerror_;
  std::vector<double>pterror_;
  std::vector<double> trackstree;
  std::vector<double> multtree;
  std::vector<double> centralitytree;
  std::vector<double> Pxtree;
  std::vector<double> Pytree;
  std::vector<double> Pztree;
  std::vector<double> Pttree;
  std::vector<double> energytree;
  std::vector<double> Etatree;
  std::vector<double> Phitree;
  std::vector<double> chargetree;
  std::vector<double> nhitstree;
  std::vector<double> PtNPtree;
  std::vector<double> EtaNPtree;
  std::vector<double> PhiNPtree;
  std::vector<double> MNPtree;
  std::vector<double> PtPPtree;
  std::vector<double> EtaPPtree;
  std::vector<double> PhiPPtree;
  std::vector<double> MPPtree;
  std::vector<double> Ptsumtree;
  std::vector<double> Etasumtree;
  std::vector<double> Phisumtree;
  std::vector<double> rhopttree;
  std::vector<double> rhoetatree;
  std::vector<double> rhophitree;
  std::vector<double> rhomtree;
  std::vector<double> Msumtree;
  std::vector<double> hfptree;
  std::vector<double> hfmtree;
  std::vector<double> ZDCPlustree;
  std::vector<double> ZDCMinustree;
  std::vector<double> ZDCenergytree;
  std::vector<double> Npixeltree;
  std::vector<double> Ntrkofftree;
  std::vector<double> ettower;
  std::vector<double> etatower;
  std::vector<double> phitower;
  std::vector<double> HFMinustree;
  std::vector<double> HFPlustree; 
  std::vector<double> genetatree;
  std::vector<double> genpttree;
  std::vector<double> genphitree;
  std::vector<double> genchargetree;
  std::vector<double> genpdgidtree;
  
  // ## track selection ##
  double etamin_;    //min eta of the tracks
  double etamax_;    //max eta of the tracks
  double ptmin_;     //min pt of the tracks
  double ptmax_;     //max pt of the tracks
  double dzdzerror_; //cut on dz/dzerror of the tracks
  double d0d0error_; //cut on d0/d0error of the tracks
  double pterrorpt_; //cut on pterror/pt of the tracks
  int mult_;         //multiplicity (Nref) in a given event
  int multcorr_;     //multiplicity corrected (Nch) in a given event
      
  double p1;
  double p2;
  double e1;
  double e2;
  
  // ## vertex selection ##
  double  minvz_;         //minimum z distance wrt (0,0,0) for the vertex       
  double  maxvz_;         //maximum z distance wrt (0,0,0) for the vertex
  double  maxrho_;        //cut on rho distance for the vertex position 
  bool    isBVselByMult_; //sel best vertex based on vertex multiplicity (true) or sum p_T^2 (false)
  int     nvtx_;          //number of reconstructed vertices in a given events
  double  xBestVtx_;      //x coordinate of the best vertex
  double  yBestVtx_;      //y coordinate of the best vertex
  double  rhoBestVtx_;    //r coordinate of the best vertex
  double  zBestVtx_;      //z coordinate of the best vertex
  double  xBestVtxError_; //x coordinate error of the best vertex
  double  yBestVtxError_; //y coordinate error of the best vertex
  double  zBestVtxError_; //z coordinate error of the best vertex
  int  run;
  int  event;
  // ## file acc & eff & fake ##
  edm::InputTag fname_;         //file name that contains acc X eff corrections
  std::vector<int> effmultbin_; //Multiplicity binning of the correction 
  TFile* feff_;                 //TFile that contains 2D histos (pt, eta) with eff/(1-fak) 
  std::vector<TH2D*> heff_;     //TH2D histograms used for correction
  TTree* ptree;
    
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

