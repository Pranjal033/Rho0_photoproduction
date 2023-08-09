#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TH1D.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include <math.h>
#include <fstream>

using namespace std;
void myrhocode(){
  

// Defining histograms
  TH1D * PT = new TH1D("pt","pt",1000,0,2);
  TH1D * ETA = new TH1D("eta","eta",1000,-5,5);
  TH1D * PHI = new TH1D("phi","phi",1000,-4,4);
  TH1D * CHRG = new TH1D("charge","charge",100,-4,4);  
  TH1D * SUM_PT = new TH1D("sumpt","sumpt",100,0,5);
  TH1D * SUM_ETA = new TH1D("sumeta","sumeta",100,-5,5);
  TH1D * SUM_PHI = new TH1D("sumphi","sumphi",100,-4,4);
  TH1D * SUM_IM = new TH1D("sumIM","sumIM",200,0,3);
  TH1D * SUM_GPT = new TH1D("gsumpt","gsumpt",1000,0,0.5);
  TH1D * SUM_GETA = new TH1D("gsumeta","gsumeta",1000,-10,10);
  TH1D * SUM_GPHI = new TH1D("gsumphi","gsumphi",100,-4,4);
  TH1D * SUM_GIM = new TH1D("gsumIM","gsumIM",500,0,3);
  
// filelist: list of all my root files
  ifstream myfile("filelist.txt", ios::in);
  TString str;
  TString str2 = "";
  myfile >> str;
  
  while(str != "")
    {

      TFile * file = TFile::Open(str,"READ");
      TTree * tree = (TTree *)(file->Get("defaultAnalysis/ptree"));
   
      std::vector<double> *pt = 0;
      TBranch  *b_pt;
      std::vector<double> *eta = 0;
      TBranch  *b_eta;
      std::vector<double> *phi = 0;
      TBranch  *b_phi;
      std::vector<double> *charge = 0;
      TBranch  *b_charge;
      std::vector<double> *dzv = 0;
      TBranch  *b_dzv;
      std::vector<double> *dxyv = 0;
      TBranch  *b_dxyv;
      std::vector<double> *dze = 0;
      TBranch  *b_dze;
      std::vector<double> *dxye = 0;
      TBranch  *b_dxye;
      std::vector<double> *pte = 0;
      TBranch  *b_pte;
      std::vector<double> *nhit = 0;
      TBranch  *b_nhit;
      std::vector<double> *gpt = 0;
      TBranch  *b_gpt;
      std::vector<double> *geta = 0;
      TBranch  *b_geta;
      std::vector<double> *gphi = 0;
      TBranch  *b_gphi;
      std::vector<double> *gcharge = 0;
      TBranch  *b_gcharge;
      std::vector<double> *pdgid = 0;
      TBranch  *b_pdgid;
   


    //getting branches from the tree
    tree->SetBranchAddress("Pt",&pt,&b_pt);
    tree->SetBranchAddress("Eta",&eta,&b_eta);
    tree->SetBranchAddress("Phi",&phi,&b_phi);
    tree->SetBranchAddress("Charge",&charge,&b_charge);
    tree->SetBranchAddress("dzvtx",&dzv,&b_dzv);
    tree->SetBranchAddress("dxyvtx",&dxyv,&b_dxyv);
    tree->SetBranchAddress("dzerror",&dze,&b_dze);
    tree->SetBranchAddress("dxyerror",&dxye,&b_dxye);
    tree->SetBranchAddress("pterror",&pte,&b_pte);
    tree->SetBranchAddress("nhits",&nhit,&b_nhit);
    tree->SetBranchAddress("gpt",&gpt,&b_gpt);
    tree->SetBranchAddress("geta",&geta,&b_geta);
    tree->SetBranchAddress("gphi",&gphi,&b_gphi);
    tree->SetBranchAddress("gcharge",&gcharge,&b_gcharge);
    tree->SetBranchAddress("pdgid",&pdgid,&b_pdgid);

   
 int entries = tree->GetEntries();
			   
			   
			   for(int i = 0; i < entries; i++)
			     {
			       tree->GetEntry(i);
			       b_pt->GetEntry(i);
			       b_eta->GetEntry(i);
			       b_phi->GetEntry(i);
			       b_charge->GetEntry(i);
			       b_dzv->GetEntry(i);
			       b_dxyv->GetEntry(i);
			       b_dze->GetEntry(i);
			       b_dxye->GetEntry(i);
			       b_pte->GetEntry(i);
			       b_nhit->GetEntry(i);
			       b_geta->GetEntry(i);
			       b_gphi->GetEntry(i);
			       b_gpt->GetEntry(i);
			       b_gcharge->GetEntry(i);
			       b_pdgid->GetEntry(i);


			       //FILLING THE RECO INFORMATION:
			       if((*pt).size()==2)
				 {
				   
				   for(unsigned int i=0 ; i<(*pt).size() ; i++)
				     {
					// Track conditions
				       //if( fabs((*pte)[i]/(*pt)[i]) >= 0.1 ) continue;
				       //if( fabs((*dzv)[i] / (*dze)[i]) >= 20 ) continue;
				       //if( fabs((*dxyv)[i] / (*dxye)[i]) >= 7 ) continue;
				       //if((*charge)[i] == 0 ) return;
				       //if((*nhit)[i] >=6) continue;
		
				       TLorentzVector TL1;
				       TL1.SetPtEtaPhiM((*pt)[i],(*eta)[i],(*phi)[i],0.139);
				       for(unsigned int j=i+1 ; j<(*pt).size() ; j++)
					 {
					   if(i==j) continue;
					   if(((*charge)[i]>0 && (*charge)[j]<0) || ((*charge)[i]<0 && (*charge)[j]>0))
					     {
					       TLorentzVector sum,TL2;
					       TL2.SetPtEtaPhiM((*pt)[j],(*eta)[j],(*phi)[j],0.139);
					       sum = TL1+TL2;                          
					       
					       SUM_PT->Fill(sum.Pt());
					       SUM_ETA->Fill(sum.Eta());
					       SUM_PHI->Fill(sum.Phi());
					       SUM_IM->Fill(sum.M());
					     }
					 }
				     }
				 }
			     
			     
			       
			       for(int k=0;k<(*pt).size();k++)
				 {
				   float ptt =(*pt)[k];
				   PT->Fill(ptt);
				   float etaa =(*eta)[k];
				   ETA->Fill(etaa);
				   float phii =(*phi)[k];
				   PHI->Fill(phii);
				   float chargee =(*charge)[k];
				   CHRG->Fill(chargee);
				 }
			       


			//FILLING THE GENERATOR LEVEL INFORMATION
			//METHOD 1 (with PID)
			       if((*gpt).size()==2)
				{
			       for(int l=0;l<(*gpt).size();l++)
                                 {
				   if(abs((*pdgid)[l])!=211) continue;
				
				   TLorentzVector TL3;
				   TL3.SetPtEtaPhiM((*gpt)[l],(*geta)[l],(*gphi)[l],0.13957018);
				 
                                 
			       for(int m=l+1;m<(*gpt).size();m++)
                                 {
                                   if((*pdgid)[l] + (*pdgid)[m] != 0) continue;
                      
				   TLorentzVector TL4;
                                   TL4.SetPtEtaPhiM((*gpt)[m],(*geta)[m],(*gphi)[m],0.13957018);

				   TLorentzVector sum2;
				   sum2=TL3+TL4;
				   SUM_GPT->Fill(sum2.Pt());                                         
				   SUM_GETA->Fill(sum2.Eta());                                       
				   SUM_GPHI->Fill(sum2.Phi());                                       
				   SUM_GIM->Fill(sum2.M());
				 
				 }
				 }
				 }

			//METHOD 2 (with charge)
			       
			      /* if((*gpt).size()==2)
			       {
				   for(unsigned int i=0 ; i<(*gpt).size() ; i++)
				     {
				       TLorentzVector TL3;
				       
					 TL3.SetPtEtaPhiM((*gpt)[i],(*geta)[i],(*gphi)[i],0.13957018);
					 for(unsigned int j=i+1 ; j<(*gpt).size() ; j++)
					   {
					     if(((*gcharge)[i]>0 && (*gcharge)[j]<0) || ((*gcharge)[i]<0 && (*gcharge)[j]>0))
					       {
						 TLorentzVector TL4,sum2;
						 TL4.SetPtEtaPhiM((*gpt)[j],(*geta)[j],(*gphi)[j],0.13957018);
						 
						 sum2=TL3+TL4;
						 SUM_GPT->Fill(sum2.Pt());
						 SUM_GETA->Fill(sum2.Eta());
						 SUM_GPHI->Fill(sum2.Phi());
						 SUM_GIM->Fill(sum2.M());
					       }
					   }
				     }
			       }*/

		//NOTE: Method 1 and method 2 are giving same results
	   
			     }	   
				   
			   //delete tree;
				   myfile >> str;
				   file->Close();
			     }
  
  
  
    TFile * file1 = new TFile("Aug9_2.root","RECREATE");
    PT->Write();
    ETA->Write();
    PHI->Write();
    CHRG->Write();
    SUM_PT->Write();
    SUM_ETA->Write();
    SUM_PHI->Write();
    SUM_IM->Write();
    SUM_GPT->Write();
    SUM_GETA->Write();
    SUM_GPHI->Write();
    SUM_GIM->Write();
}
