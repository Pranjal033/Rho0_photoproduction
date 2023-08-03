// Package:    Analyzers/NtrkDistribution
// Class:      NtrkDistribution
 
/*class NtrkDistribution NtrkDistribution.cc Analyzers/NtrkDistribution/plugins/NtrkDistribution.cc
 Description: [one line class summary]
 Implementation:
     [Notes on implementation]
 Original Author:  Maxime Guilbaud
 Created:  Tue, 24 Oct 2017 13:23:04 GMT
*/

// system include files
#include <memory>

// CMSSW include files
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// user include files
#include "Analyzers/NtrkDistribution/interface/NtrkDistribution.h"

// constructors and destructor
NtrkDistribution::NtrkDistribution(const edm::ParameterSet& iConfig) :
  //tracks
  trackTags_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("tracks"))),
  //vertex
  vtxTags_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertex"))),
  //caloTower
  caloTowersTags_(consumes<CaloTowerCollection>(iConfig.getParameter<edm::InputTag>("caloTower"))),
  //centrality
  centralityTags_(consumes<reco::Centrality>(iConfig.getParameter<edm::InputTag>("centralitySrc"))),
  centralityBinTags_(consumes<int>(iConfig.getParameter<edm::InputTag>("centralityBinSrc"))),
  //Generatorinfo
  trackTagsgen_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("tracksgen"))),

  //track selection
  etamin_(iConfig.getUntrackedParameter<double>("etamin")),
  etamax_(iConfig.getUntrackedParameter<double>("etamax")),
  ptmin_(iConfig.getUntrackedParameter<double>("ptmin")),
  ptmax_(iConfig.getUntrackedParameter<double>("ptmax")),
  dzdzerror_(iConfig.getUntrackedParameter<double>("dzdzerror")),
  d0d0error_(iConfig.getUntrackedParameter<double>("d0d0error")),
  pterrorpt_(iConfig.getUntrackedParameter<double>("pterrorpt")),
  //vertex selection
  minvz_(iConfig.getUntrackedParameter<double>("minvz")),
  maxvz_(iConfig.getUntrackedParameter<double>("maxvz")),
  maxrho_(iConfig.getUntrackedParameter<double>("maxrho")),
  isBVselByMult_(iConfig.getUntrackedParameter<bool>("isBVselByMult")),
  nvtx_(iConfig.getUntrackedParameter<int>("nvtx",-1)),
  xBestVtx_(iConfig.getUntrackedParameter<double>("xVtx",-99999)),
  yBestVtx_(iConfig.getUntrackedParameter<double>("yVtx",-99999)),
  rhoBestVtx_(iConfig.getUntrackedParameter<double>("rhoVtx",-99999)),
  zBestVtx_(iConfig.getUntrackedParameter<double>("zVtx",-99999)),
  xBestVtxError_(iConfig.getUntrackedParameter<double>("xVtxError",-99999)),
  yBestVtxError_(iConfig.getUntrackedParameter<double>("yVtxError",-99999)),
  zBestVtxError_(iConfig.getUntrackedParameter<double>("zVtxError",-99999)),
//file acc & eff & fake
//change - removing efficiency file.
fname_(iConfig.getUntrackedParameter<edm::InputTag>("fname")),
effmultbin_(iConfig.getUntrackedParameter< std::vector<int> >("effmultbin"))
  
  //now do what ever initialization is needed
  /* TString filename(fname_.label().c_str());
     feff_ = 0x0;
     if(!filename.IsNull())
     {
     edm::FileInPath fip(Form("Analyzers/NtrkDistribution/data/EFF/%s",filename.Data()));
     feff_ = new TFile(fip.fullPath().c_str(),"READ");
     
      heff_.resize(feff_->GetNkeys());
      if(heff_.size() != effmultbin_.size() - 1)
      {
      edm::LogWarning ("Inconsitent binning") << " Inconsistent binning for the acc X eff correction..."
      << " You might have wrong setting here";
      }
      
      for(unsigned short ik = 0; ik < heff_.size(); ++ik)
      {
      heff_[ik] = (TH2D*) feff_->Get(feff_->GetListOfKeys()->At(ik)->GetName());
      }
      }*/
{
  //Output
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  //Histogram for check
  TFileDirectory fTrkHist  = fs->mkdir("Tracks");
  hPtTrk_  = fTrkHist.make<TH1F>("hPttrk",  "", 100,  0.,  10.);
  //Converting in tree format

  ptree = fs->make<TTree>("ptree","Pion's Tree");

  //Vertex Branches
  ptree->Branch("XVtx",&xBestVtx_,"xBestVtx_/D");
  ptree->Branch("YVtx",&yBestVtx_,"yBestVtx_/D");
  ptree->Branch("ZVtx",&zBestVtx_,"xBestVtx_/D");
  ptree->Branch("RhoVtx",&rhoBestVtx_,"rhoBestVtx_/D");  
  ptree->Branch("run_no",&run,"run/I");
  ptree->Branch("event_no",&event,"event/I");
  ptree->Branch("dzvtx",&dzvtx_);
  ptree->Branch("dxyvtx",&dxyvtx_);
  ptree->Branch("dzerror",&dzerror_);
  ptree->Branch("dxyerror",&dxyerror_);
  ptree->Branch("pterror",&pterror_);

  //Track Branches
  ptree->Branch("tracks",&trackstree);
  ptree->Branch("mult",&multtree);
  ptree->Branch("Px",&Pxtree);
  ptree->Branch("Py",&Pytree);
  ptree->Branch("Pz",&Pztree); 
  ptree->Branch("Pt",&Pttree);
  ptree->Branch("Energy",&energytree);
  ptree->Branch("Eta",&Etatree);
  ptree->Branch("Phi",&Phitree);
  ptree->Branch("Charge",&chargetree);
  ptree->Branch("nhits",&nhitstree);
  //Pi+ + Pi-
  ptree->Branch("sumPt",&Ptsumtree);
  ptree->Branch("sumEta",&Etasumtree);
  ptree->Branch("sumPhi",&Phisumtree);
  ptree->Branch("sumM",&Msumtree);
  //Rho
  //ptree->Branch("rhoPt",&rhopttree);
  //ptree->Branch("rhoEta",&rhoetatree);
  //ptree->Branch("rhoPhi",&rhophitree);
  //ptree->Branch("rhoM",&rhomtree);
  //Centrality
  ptree->Branch("Centrality",&centralitytree);
  ptree->Branch("HFplus",&hfptree);
  ptree->Branch("HFminus",&hfmtree);
  ptree->Branch("ZDCplus",&ZDCPlustree);
  ptree->Branch("ZDCminus",&ZDCMinustree);
  ptree->Branch("Npixel",&Npixeltree);
  ptree->Branch("Ntrkoffline",&Ntrkofftree);
  //Calotower
  ptree->Branch("HFMinustower",&HFMinustree);
  ptree->Branch("HFPlustower",&HFPlustree);
  ptree->Branch("Ettow",&ettower);
  ptree->Branch("Etatow",&etatower);
  ptree->Branch("Phitow",&phitower);
  //GenParticles
  ptree->Branch("gpt",&genpttree);
  ptree->Branch("geta",&genetatree);
  ptree->Branch("gphi",&genphitree);
  ptree->Branch("gcharge",&genchargetree);
  ptree->Branch("pdgid",&genpdgidtree);
}

NtrkDistribution::~NtrkDistribution()
{
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}

// ------------ method called for each event  ------------
void
NtrkDistribution::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  
  // ----- Vertex selection -----
  // Get vertex collection by token
  edm::Handle< reco::VertexCollection > vertices;
  iEvent.getByToken(vtxTags_, vertices);
  reco::VertexCollection verticesColl = *vertices;
  if( !vertices->size() ) 
    { 
      edm::LogWarning ("Missing Collection") <<"Invalid or empty vertex collection!";
      return; 
    }
  nvtx_          = 0;     //N valid vertex in collection
  xBestVtx_      = -999.;
  yBestVtx_      = -999.; 
  rhoBestVtx_    = -999.; 
  zBestVtx_      = -999.; //Best vtx coordinates
  xBestVtxError_ = -999.; 
  yBestVtxError_ = -999.; 
  zBestVtxError_ = -999.; //Best vtx error
  
  // Sort vertex collection if you want to select highest multiplicity vertex as best vertex
  /*if(isBVselByMult_)
    {
      std::sort(verticesColl.begin(), verticesColl.end(), [](const reco::Vertex &a, const reco::Vertex &b)
		{
		   if ( a.tracksSize() == b.tracksSize() ) 
		     return a.chi2() < b.chi2();
		    return a.tracksSize() > b.tracksSize();
		});
		}*/
 
  // Loop over vertices
  for( reco::VertexCollection::const_iterator itVtx = verticesColl.begin();
       itVtx != verticesColl.end();
       ++itVtx )
    {
      // Drop fake vertex and vertex with less than 2 tracks attached to it
      if( !itVtx->isFake() && itVtx->tracksSize() >= 2 )
        {
	  // x,y,z vertex position
	  double xVtx = itVtx->x();
	  double yVtx = itVtx->y();
	  double zVtx = itVtx->z();
	  // x,y,z vertex position error
	  double xVtxError = itVtx->xError();
	  double yVtxError = itVtx->yError();
	  double zVtxError = itVtx->zError();
	  // Radial vertex position                                                         
	  double rho = sqrt(xVtx*xVtx + yVtx*yVtx);
	  // Increase N valid vertex in the collection
	  ++nvtx_;
	  
	  //Get the first vertex as the best one (greatest sum p_{T}^{2}) 
	  if( itVtx == verticesColl.begin() )
	    {
	      xBestVtx_ = xVtx; 
	      yBestVtx_ = yVtx; 
	      zBestVtx_ = zVtx; 
	      xBestVtxError_ = xVtxError; 
	      yBestVtxError_ = yVtxError; 
	      zBestVtxError_ = zVtxError;
	      rhoBestVtx_ = rho; 
	    }
        }
    }
  
  //Select event using vertex properties
  if ( nvtx_ <= 0 )                               return;
  //if ( zBestVtx_ < minvz_ || zBestVtx_ > maxvz_ ) return;
  //if ( rhoBestVtx_ > maxrho_ )                    return;

  run = iEvent.id().run();
  event = iEvent.id().event();  


  // fill centrality information
  edm::Handle<reco::Centrality> cent;
  iEvent.getByToken(centralityTags_, cent);
  if (cent.isValid())
    {   
      double hfp = cent->EtHFtowerSumPlus();
      hfptree.push_back(hfp);
      double hfm = cent->EtHFtowerSumMinus();
      hfmtree.push_back(hfm);
      double zdcp =cent->zdcSumPlus();
      ZDCPlustree.push_back(zdcp);
      double zdcm =cent->zdcSumMinus();
      ZDCMinustree.push_back(zdcm);
      Npixeltree.push_back(cent->multiplicityPixel());
      Ntrkofftree.push_back(cent->Ntracks()); 
    }

  edm::Handle<int> centBin;
  iEvent.getByToken(centralityBinTags_, centBin);
  if (centBin.isValid())
    {
      centralitytree.push_back(*centBin);
    }

  //HF Energy --------------------------------------------

  double etHFtowerSumPlus=0;
  double etHFtowerSumMinus=0;
  double etHFtowerSum=0;
  double max_HFP=0;
  double max_HFM=0;
  Handle<CaloTowerCollection> towers;
  iEvent.getByLabel("towerMaker",towers);
  for( size_t i = 0; i<towers->size(); ++ i)
    {
    const CaloTower & tower = (*towers)[ i ];
    double aeta = abs(tower.eta());
    bool isHF = aeta>3 && aeta <6;
    if(isHF && tower.eta() > 0){
      etHFtowerSumPlus += tower.energy();
      HFPlustree.push_back(tower.energy());
      if(tower.energy()>max_HFP)  
	max_HFP = tower.energy();
    }
    if(isHF && tower.eta() < 0){
      etHFtowerSumMinus += tower.energy();
      HFMinustree.push_back(tower.energy());
      if(tower.energy()>max_HFM)
	max_HFM = tower.energy();
    }
    }
  etHFtowerSum=etHFtowerSumPlus + etHFtowerSumMinus;
  std::cout<<"calculated HF energy = "<<etHFtowerSum<<", HF energy from centrality = "<<cent->EtHFtowerSum()<<std::endl;
  // std::cout<<"max_HFP = " <<max_HFP << "Max_HFMM = " << max_HFM<<std::endl;
  //Energy cuts
  if(max_HFP>=7.6) return;
  if(max_HFM>=7.3) return;   

  // ----- Track selection -----
  // Get track collection by token
  edm::Handle< reco::TrackCollection > tracks;
  iEvent.getByToken(trackTags_, tracks);
  if( !tracks->size() )
    {
      edm::LogWarning ("Missing Collection") <<"Invalid or empty track collection!";
      // return;
    }
  
  // Pranjal -Ntrack offline event selection -DELETED 
  
  std::vector<double> Vpt,Veta,Vphi,Vcharge,Vp1,Vp2,Ve1,Ve2, Vpx1,Vpx2,Vpy1,Vpy2,Vpz1,Vpz2;
  for( reco::TrackCollection::const_iterator itTrk = tracks->begin();
       itTrk != tracks->end();
       ++itTrk )
    {
      // std::cout<< "My Trackssssssssssssssss pt = " << itTrk->pt() <<std::endl;
      // Select tracks based on proximity to best vertex
      math::XYZPoint bestvtx(xBestVtx_,yBestVtx_,zBestVtx_);
      double dzvtx    = itTrk->dz(bestvtx);
      double dxyvtx   = itTrk->dxy(bestvtx);
      double dzerror  = sqrt(itTrk->dzError()*itTrk->dzError() + zBestVtxError_*zBestVtxError_);
      double dxyerror = sqrt(itTrk->d0Error()*itTrk->d0Error() + xBestVtxError_*yBestVtxError_);
      double pterror  = itTrk->ptError();
      
      // Get eta, pt, phi and charge of the track      
      double eta    = itTrk->eta();
      double pt     = itTrk->pt();
      double phi    = itTrk->phi();
      double charge = itTrk->charge();
      double nhits  = itTrk->numberOfValidHits();
      double px = itTrk->px();
      double py= itTrk->py();
      double pz= itTrk->pz();
  
      // Select track based on quality
      if( !itTrk->quality(reco::TrackBase::highPurity) ) continue;
      //if( fabs(pterror) / pt      >= pterrorpt_ ) continue;
      //if( fabs(dzvtx / dzerror)   >= dzdzerror_ ) continue;
      //if( fabs(dxyvtx / dxyerror) >= d0d0error_ ) continue;
      //if( pt <= 0.4 ) continue;
      //if( charge == 0 ) continue;
      //if(nhits < 6) continue;
      //if( pt < 0.0001 ) continue;                                                                    
      
      // Track selection for analysis(No other dca cuts for merged general and pixel track condition.

      if(abs(eta) > 2.5) continue;
      //if(pt < ptmin_ || pt > ptmax_) continue;

      //Filling histograms
      hPtTrk_ ->Fill(pt);
      //Filling Trees
      Pxtree.push_back(px);
      Pytree.push_back(py);
      Pztree.push_back(pz);
      Pttree.push_back(pt);
      energytree.push_back(sqrt(px*px+py*py+pz*pz));
      Etatree.push_back(eta);
      Phitree.push_back(phi);
      chargetree.push_back(charge);
      nhitstree.push_back(nhits);
      multtree.push_back(mult_);
      dzvtx_.push_back(dzvtx);
      dxyvtx_.push_back(dxyvtx);
      dzerror_.push_back(dzerror);
      dxyerror_.push_back(dxyerror);
      pterror_.push_back(pterror);

      //Filling Vectors
      Vpt.push_back(pt);
      Veta.push_back(eta);
      Vphi.push_back(phi);
      Vcharge.push_back(charge);


        }
  //Particle loop ends here .............................................                              
  trackstree.push_back(Vpt.size());

  //if(Vpt.size()!=2) return;
  if(Vpt.size()==2)
    {
      for(unsigned int i=0 ; i<Vpt.size() ; i++)
	{      
	  TLorentzVector TL1;
	  TL1.SetPtEtaPhiM(Vpt.at(i),Veta.at(i),Vphi.at(i),0.13957039);
	  for(unsigned int j=i+1 ; j<Vpt.size() ; j++)
	    {
	      if(i==j) continue;
	      if(Vcharge.at(i)!=Vcharge.at(j))
		{
		  TLorentzVector TL2;
		  TL2.SetPtEtaPhiM(Vpt.at(j),Veta.at(j),Vphi.at(j),0.13957039);
		  auto sum = TL1+TL2;
		  
		  Ptsumtree.push_back(sum.Pt());
		  Etasumtree.push_back(sum.Eta());
		  Phisumtree.push_back(sum.Phi());
		  Msumtree.push_back(sum.M());
		  
	      }
	    }
	}
    }


  // ----- Calotower selection -----                                                                                                                    
  // Get calo tower collection by token                                                                                                                 
  edm::Handle< CaloTowerCollection > calotowers;
  iEvent.getByToken(caloTowersTags_, calotowers);
  if( !calotowers->size() )
    {
      edm::LogWarning ("Missing Collection") <<"Invalid or empty caloTower collection!";
      //return;
    }
  // Loop over caloTowers                                                                                                                            
  for( CaloTowerCollection::const_iterator itCTow = calotowers->begin();
  itCTow != calotowers->end();
  ++itCTow )
  {
      // Get eta, pt and phi of the calo tower                                                                                                      
      // Select calo tower based on quality                                                                                                         
      if( itCTow->et() < 0.01 ) continue;
      ettower.push_back(itCTow->et());
      etatower.push_back(itCTow->eta());
      phitower.push_back(itCTow->phi());
   }
  ptree->Fill();
  
  //Clearing Vector
  dzvtx_.clear();
  dxyvtx_.clear();
  dzerror_.clear();
  dxyerror_.clear();
  pterror_.clear();
  Vpt.clear();
  Veta.clear();
  Vphi.clear();
  Vcharge.clear();
  
  Pxtree.clear();
  Pytree.clear();
  Pztree.clear();
  Pttree.clear();
  centralitytree.clear();
  energytree.clear();
  multtree.clear();
  Etatree.clear();
  Phitree.clear();
  chargetree.clear();
  nhitstree.clear();
 
  Ptsumtree.clear();
  Phisumtree.clear();
  Etasumtree.clear();
  Msumtree.clear();
  //rhopttree.clear();
  //rhoetatree.clear();
  //rhophitree.clear();
  //rhomtree.clear();
  
  HFMinustree.clear();
  HFPlustree.clear();
  hfptree.clear();
  hfmtree.clear();
  ZDCPlustree.clear();
  ZDCMinustree.clear();
  Npixeltree.clear();
  Ntrkofftree.clear();
  ettower.clear();
  etatower.clear();
  phitower.clear();
  trackstree.clear();
  //}

//Reco/Data Event loop ends here ...............
// ------------ method called once each job just before starting event loop  ------------
//void 
//NtrkDistribution::Gen(const edm::Event& iEvent, const edm::EventSetup& iSetup)
//{
   edm::Handle< reco::GenParticleCollection > tracksgen;
  iEvent.getByToken(trackTagsgen_, tracksgen);
  if( !tracksgen->size() )
    {

      edm::LogWarning ("Missing MC Gen Collection") <<"Invalid or empty MC-Gen track collection!";
      return;
    }
  for( reco::GenParticleCollection::const_iterator itTrk = tracksgen->begin();
       itTrk != tracksgen->end();
	 ++itTrk )
    {

// Get eta, pt, phi and charge of the track                                                                                                                                    
      if( itTrk->status() != 1 ) continue;
      if( itTrk->charge() == 0 ) continue;
      // Get eta, pt, phi and charge of the track                                                                                                                                    
      double geta      = itTrk->eta();
      double gpt       = itTrk->pt();
      double gphi      = itTrk->phi();
      int gcharge      = itTrk->charge();
      int pdgid        = itTrk->pdgId();
      //std::cout<< "My GEN Trackssssssssssssssss pt 000000000000000000 = " << itTrk->pt() <<std::endl;

      genetatree.push_back(geta);
      genpttree.push_back(gpt);
      genphitree.push_back(gphi);
      genchargetree.push_back(gcharge);
      genpdgidtree.push_back(pdgid);
    }
  ptree->Fill();
  genetatree.clear();
  genpttree.clear();
  genphitree.clear();
  genchargetree.clear();
  genpdgidtree.clear();
}



void
NtrkDistribution::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
NtrkDistribution::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
NtrkDistribution::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(NtrkDistribution);
