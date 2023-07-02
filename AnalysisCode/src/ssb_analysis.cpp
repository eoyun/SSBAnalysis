#define ssb_analysis_cxx

#include <iostream>
#include <sstream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>

#include "TMath.h"
#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TMatrixD.h"
#include "TMatrixDEigen.h"

#include "./../interface/ssb_analysis.hpp"
#include "./../CommonTools.hpp"

using namespace std;

void ssb_analysis::TLVInitial()
{




   Top = new TLorentzVector();
   AnTop = new TLorentzVector();
   bJet = new TLorentzVector();
   AnbJet = new TLorentzVector();
   Nu = new TLorentzVector();
   AnNu = new TLorentzVector();
   W1 = new TLorentzVector();
   W2 = new TLorentzVector();
   Met = new TLorentzVector();
}

void ssb_analysis::Loop( char *logfile, string sampleName )
{
   

   //////////
   if (fChain == 0) return;
   //////////
   crosssection =1.;
   int isData=0;
   int isSingle=0;
   vector<TLorentzVector> TightEle;
   vector<TLorentzVector> TightMuon;
   if (sampleName.find("TTJets_Signal")!=string::npos) crosssection = 831.76E-12;
   if (sampleName.find("TTJets_others")!=string::npos) crosssection = 831.76E-12;
   if (sampleName.find("DYJetsToLL_M_10To50")!=string::npos) crosssection = 18810.0E-12;
   if (sampleName.find("DYJetsToLL_M_50")!=string::npos) crosssection = 5941.0E-12;
   if (sampleName.find("ST_tW_antitop")!=string::npos) crosssection = 35.6E-12;
   if (sampleName.find("ST_tW_top")!=string::npos) crosssection = 35.6E-12;
   if (sampleName.find("TTbar_WJetToLNu")!=string::npos) crosssection = 0.2043E-12;
   if (sampleName.find("TTbar_WQQ")!=string::npos) crosssection = 0.4062E-12;
   if (sampleName.find("TTbar_ZToLLNuNu")!=string::npos) crosssection = 0.2529E-12;
   if (sampleName.find("TTbar_ZQQ")!=string::npos) crosssection = 0.5297E-12;
   if (sampleName.find("WW")!=string::npos) crosssection = 118.7E-12;
   if (sampleName.find("WZ")!=string::npos) crosssection = 65.9E-12;
   if (sampleName.find("ZZ")!=string::npos) crosssection = 31.8E-12;
   if (sampleName.find("WJetsToLNu")!=string::npos) crosssection = 61526E-12;

   if (sampleName.find("Data")!=string::npos) isData=1;
   if (sampleName.find("Single")!=string::npos) isSingle=1;
   cout<<"dgb"<<endl;


   luminosity = 35.9;
   //////////
   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   float bdisccut =0.8484;
   float muonisobetacut = 0.15;
   float elecisocut = 0.0821;
   //////////

   ///My variables
   Long64_t __tot_evt = 0;

   /// Check Total Event
 
   ////////////////////////
   /// start event loop ///
   ////////////////////////
   int muon_eta_0 =0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) 
   {
      double w_pileup = puweight->weight(PileUp_Count_Intime)*Gen_EventWeight;
      if (isData) w_pileup=1.;
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0)
      {
         printf("ERROR: Could not load tree!!!\n");
         break;
      }
      TightEle.clear();
      TightMuon.clear();
  
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      int trigger_flag_mumu =0;
      int trigger_flag_elmu =0;
      int trigger_flag_elel =0;
      
      if (jentry % 10000 == 0) printf("Event %lld\n", jentry); //%lld supports Long64_t
      /*cout << "size of Trigger Name : "<<Trigger_Name->size()<<
      " size of Trigger Prescale : "<<Trigger_PreScale->size()<<
      " size of Trigger isError : "<<Trigger_isError->size()<<
      " size of Trigger isPass : "<<Trigger_isPass->size()<< 
      " size of Trigger isRun : "<<Trigger_isRun->size()<<endl;*/
      //cout << "evt num "<<jentry<<endl;
      FillHisto(h_cf_pu_intime, PileUp_Count_Intime,Gen_EventWeight);
      FillHisto(h_cf_pu_pv, PV_Count,Gen_EventWeight);
      FillHisto(h_cf_pu_interaction, PileUp_Count_Interaction,Gen_EventWeight);
      FillHisto(h_cf_pu_intime_w, PileUp_Count_Intime,w_pileup);
      FillHisto(h_cf_pu_pv_w, PV_Count,w_pileup);
      FillHisto(h_cf_pu_interaction_w, PileUp_Count_Interaction,w_pileup);
      //cout<<isData<<endl;
      if (sampleName.find("TTJets_Signal")!=string::npos && Channel_Idx != 22) continue;
      if (sampleName.find("TTJets_others")!=string::npos && Channel_Idx == 22) continue;
      if (isData==0){  
        for(int i=0;i < Trigger_Name->size();i++){
          if (Trigger_Name->at(i).find("HLT_IsoMu24_v")!=string::npos){ trigger_flag_elmu+=Trigger_isPass->at(i); trigger_flag_mumu+=Trigger_isPass->at(i);}
          if (Trigger_Name->at(i).find("HLT_IsoTKMu24_v4")!=string::npos){trigger_flag_elmu+=Trigger_isPass->at(i); trigger_flag_mumu+=Trigger_isPass->at(i);}
          if (Trigger_Name->at(i).find("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v")!=string::npos) trigger_flag_mumu+=Trigger_isPass->at(i);
          if (Trigger_Name->at(i).find("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v")!=string::npos) trigger_flag_mumu+=Trigger_isPass->at(i);
          if (Trigger_Name->at(i).find("HLT_Ele27_WPTight_Gsf_v")!=string::npos){trigger_flag_elmu+=Trigger_isPass->at(i); trigger_flag_elel+=Trigger_isPass->at(i);}
          if (Trigger_Name->at(i).find("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v")!=string::npos) trigger_flag_elel+=Trigger_isPass->at(i);
          if (Trigger_Name->at(i).find("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v")!=string::npos) trigger_flag_elmu+=Trigger_isPass->at(i);
          if (Trigger_Name->at(i).find("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v")!=string::npos) trigger_flag_elmu+=Trigger_isPass->at(i);
          /*cout<< "Name " << Trigger_Name->at(i)<<
          " | PreScale " << Trigger_PreScale->at(i)<<
          " | is Error " << Trigger_isError->at(i)<<
          " | is Pass " << Trigger_isPass->at(i)<<
          " | is Run " << Trigger_isRun->at(i)<< endl;*/
        }
      } else {
        for(int i=0;i < Trigger_Name->size();i++){
          if (Trigger_Name->at(i).find("HLT_IsoMu24_v")!=string::npos){ trigger_flag_elmu+=Trigger_isPass->at(i); trigger_flag_mumu+=Trigger_isPass->at(i);}
          if (Trigger_Name->at(i).find("HLT_IsoTKMu24_v4")!=string::npos){trigger_flag_elmu+=Trigger_isPass->at(i); trigger_flag_mumu+=Trigger_isPass->at(i);}
          if (Trigger_Name->at(i).find("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v")!=string::npos) trigger_flag_mumu+=Trigger_isPass->at(i);
          if (Trigger_Name->at(i).find("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v")!=string::npos) trigger_flag_mumu+=Trigger_isPass->at(i);
          if (Trigger_Name->at(i).find("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v")!=string::npos) trigger_flag_mumu+=Trigger_isPass->at(i);
          if (Trigger_Name->at(i).find("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v")!=string::npos) trigger_flag_mumu+=Trigger_isPass->at(i);
          if (Trigger_Name->at(i).find("HLT_Ele27_WPTight_Gsf_v")!=string::npos){trigger_flag_elmu+=Trigger_isPass->at(i); trigger_flag_elel+=Trigger_isPass->at(i);}
          if (Trigger_Name->at(i).find("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v")!=string::npos) trigger_flag_elel+=Trigger_isPass->at(i);
          if (Trigger_Name->at(i).find("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v")!=string::npos) trigger_flag_elmu+=Trigger_isPass->at(i);
          if (Trigger_Name->at(i).find("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v")!=string::npos) trigger_flag_elmu+=Trigger_isPass->at(i);
          if (Trigger_Name->at(i).find("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v")!=string::npos) trigger_flag_elmu+=Trigger_isPass->at(i);
          if (Trigger_Name->at(i).find("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v")!=string::npos) trigger_flag_elmu+=Trigger_isPass->at(i);
          /*cout<< "Name " << Trigger_Name->at(i)<<
          " | PreScale " << Trigger_PreScale->at(i)<<
          " | is Error " << Trigger_isError->at(i)<<
          " | is Pass " << Trigger_isPass->at(i)<<
          " | is Run " << Trigger_isRun->at(i)<< endl;*/
        }

      }
      int metfilter_flag = 0;
      //cout<<"diel : "<< trigger_flag_elel <<" | is single "<<isSingle <<endl;
      //cout<<"########MET filter line ##########"<<endl;
      for (int i=0;i<METFilter_Name->size();i++){
	if (METFilter_Name->at(i).find("Flag_HBHENoiseFilter")!=string::npos) {
	  metfilter_flag+=METFilter_isPass->at(i);
	  //cout<<METFilter_isPass->at(i)<<METFilter_Name->at(i).compare("Flag_HBHENoiseFilter")<<endl;
	}
	if (METFilter_Name->at(i).find("Flag_HBHENoiseIsoFilter")!=string::npos) metfilter_flag+=METFilter_isPass->at(i);
	if (METFilter_Name->at(i).find("Flag_muonBadTrackFilter")!=string::npos) metfilter_flag+=METFilter_isPass->at(i);
	if (METFilter_Name->at(i).find("Flag_CSCTightHaloFilter")!=string::npos) metfilter_flag+=METFilter_isPass->at(i);
	if (METFilter_Name->at(i).find("Flag_goodVertices")!=string::npos) metfilter_flag+=METFilter_isPass->at(i);
	if (METFilter_Name->at(i).find("Flag_chargedHadronTrackResolutionFilter")!=string::npos) metfilter_flag+=METFilter_isPass->at(i);
	if (METFilter_Name->at(i).find("Flag_EcalDeadCellTriggerPrimitiveFilter")!=string::npos) metfilter_flag+=METFilter_isPass->at(i);
	if (isData) {if (METFilter_Name->at(i).find("Flag_eeBadScFilter")!=string::npos) metfilter_flag+=METFilter_isPass->at(i);}
	/*cout<< "Name " << METFilter_Name->at(i)<<
	" | is Error " << METFilter_isError->at(i)<<
	" | is Pass " << METFilter_isPass->at(i)<<
	" | is Run " << METFilter_isRun->at(i)<< endl;*/
	
      }
      //cout<<metfilter_flag <<" : MET filter"<<std::endl;
      //cout<<trigger_flag_elel <<" : triggerfilter"<<std::endl;
      if(isData==0 && metfilter_flag!=7) continue;
      if(isData==1 && metfilter_flag!=8) continue;
      if (trigger_flag_elel==0) continue;
      if (isSingle&&trigger_flag_elel==2) continue;
      TLVInitial();
      //cout<<"hello"<<endl;
      //cout<<"########MET filter Add line ##########"<<endl;
      for (int i=0;i<METFilterAdd_Name->size();i++){
	/*cout<< "Name "<<METFilterAdd_Name->at(i)<<
	" | is Pass " <<METFilterAdd_isPass->at(i)<<endl;*/
      }
      //cout<<"dgb"<<endl;
      
      __tot_evt++;
      ////////////////////////////////////////
      /// start Main Loop Function and Cut ///
      ////////////////////////////////////////
      TLorentzVector* ele;
      TLorentzVector mother_ele;
      TLorentzVector* leading_ele;
      TLorentzVector* sub_ele;
      int leadpt_index_ele=-1;
      float leadpt=0.0;
      int subpt_index_ele=-1;
      float subpt=0.0;
      int number_of_tight_ele=0;
      for(int i = 0; i < Elec->GetEntries(); ++i ){
	
         ele = (TLorentzVector*)Elec->At(i);
	 if (ele->Pt()>20 && abs(ele->Eta())<2.4 && Elec_SCB_Tight->at(i) == true && Elec_PFIsoRho03->at(i) < elecisocut && !(TMath::Abs(Elec_Supercluster_Eta->at(i)) > 1.4442 && TMath::Abs(Elec_Supercluster_Eta->at(i)) < 1.5660)){
	   number_of_tight_ele+=1;
	   if (leadpt<ele->Pt()){
	     subpt_index_ele = leadpt_index_ele;
	     subpt = leadpt;
	     leadpt=ele->Pt();
	     leadpt_index_ele=i;
	   }
	   else if (subpt<ele->Pt()) {
	     subpt_index_ele = i;
	     subpt=ele->Pt();
	   }
	   TightEle.push_back(*ele);
	 } 
      }
	   
      int third_ele_veto=0;
      for(int i = 0; i < Elec->GetEntries(); ++i ){
         if(i==subpt_index_ele) continue;
         if(i==leadpt_index_ele) continue;
         ele =(TLorentzVector*)Elec->At(i);
	 if(Elec_SCB_Veto->at(i)==true&&ele->Pt()>20&&abs(ele->Eta())<2.4){
	   third_ele_veto=1;
	 }
      }
      for(int i = 0; i < Elec->GetEntries(); ++i )
      {
         ele = (TLorentzVector*)Elec->At(i);
         //cout << "Electron pT " << ele->Pt() << endl;
         FillHisto(h_cf_elept, ele->Pt());
         FillHisto(h_cf_eleeta, ele->Eta());
         FillHisto(h_cf_elephi, ele->Phi());
	 if (ele->Pt()>20 && abs(ele->Eta())<2.4 && Elec_SCB_Tight->at(i) == true && Elec_PFIsoRho03->at(i) < elecisocut){
           FillHisto(h_cf_elept_select, ele->Pt());
           FillHisto(h_cf_eleeta_select, ele->Eta());
           FillHisto(h_cf_elephi_select, ele->Phi());
	   
 	 }
      }
      //cout<<"dgbele"<<endl;
      TLorentzVector* met;
      for(int i = 0; i < MET->GetEntries(); ++i )
      {
         met = (TLorentzVector*)MET->At(i);
         //cout << "METtron pT " << met->Pt() << endl;
         FillHisto(h_cf_metpt, met->Pt());
         FillHisto(h_cf_metphi, met->Phi());
      }
      TLorentzVector* metcleancor;
      for(int i = 0; i < METMUCleanCor->GetEntries(); ++i )
      {
	 if (isData) metcleancor = (TLorentzVector*)METMUEGCleanCor->At(i);
         else metcleancor = (TLorentzVector*)METMUCleanCor->At(i);
         //cout << "METtron pT " << met->Pt() << endl;
         FillHisto(h_cf_metpt_select, metcleancor->Pt());
         FillHisto(h_cf_metphi_select, metcleancor->Phi());
      }
       
      if(isData){ if (METMUEGCleanCor->GetEntries()>0)metcleancor = (TLorentzVector*)METMUEGCleanCor->At(0);}
      else {if (METMUCleanCor->GetEntries()>0)metcleancor = (TLorentzVector*)METMUCleanCor->At(0);}
      //cout<<metcleancor->Pt()<<"the met "<<endl; 
      //cout<<"dgbmet"<<endl;
      Met = (TLorentzVector*)metcleancor;

      //cout<<Met->Pt()<<": MET | "<<metcleancor->Pt()<<endl;
      TLorentzVector* muon;
      TLorentzVector* leading_muon;
      TLorentzVector* sub_muon;
      int leadpt_index_muon=-1;
      float leadpt_muon=0.0;
      int subpt_index_muon=-1;
      float subpt_muon=0.0;
      int third_muon_veto=0;
      for(int i = 0; i < Muon->GetEntries(); ++i )
      {
         muon = (TLorentzVector*)Muon->At(i);
	 if (Muon_isLoose->at(i)==true&& Muon_PFIsodBeta04->at(i) < 0.25&& muon->Pt()>20 && abs(muon->Eta())<2.4 ) third_muon_veto=1;
	 if (leadpt_muon<muon->Pt()){
	   subpt_index_muon = leadpt_index_muon;
	   subpt_muon = leadpt_muon;
	   leadpt_muon=muon->Pt();
	   leadpt_index_muon=i;
	 }
	 else if (subpt_muon<muon->Pt()) {
	   subpt_index_muon = i;
	   subpt_muon=muon->Pt();
	 }
         if (muon->Eta()>=0&&muon->Eta()<0.01){
           FillHisto(h_cf_muonpt_zeroeta, muon->Pt());
           FillHisto(h_cf_muoneta_zeroeta,muon->Eta());
           FillHisto(h_cf_muonphi_zeroeta, muon->Phi());
	   //pt eta phi
	 }
         FillHisto(h_cf_muonpt, muon->Pt());
         FillHisto(h_cf_muoneta,muon->Eta());
         FillHisto(h_cf_muonphi, muon->Phi());
	 if (muon->Pt()<2){
           FillHisto(h_cf_muonpt_check, muon->Pt());
           FillHisto(h_cf_muoneta_check,muon->Eta());
           FillHisto(h_cf_muonphi_check, muon->Phi());
	 }
	 if (muon->Pt()>20 && abs(muon->Eta())<2.4 && Muon_isTight->at(i) == true && Muon_PFIsodBeta04->at(i) < muonisobetacut){
           FillHisto(h_cf_muonpt_select, muon->Pt());
           FillHisto(h_cf_muoneta_select,muon->Eta());
           FillHisto(h_cf_muonphi_select, muon->Phi());
	 }
	 if (muon->Pt()>20 && abs(muon->Eta())<2.4 && Muon_isLoose->at(i) == true && Muon_PFIsodBeta04->at(i) < muonisobetacut){
           FillHisto(h_cf_muonpt_loose, muon->Pt());
           FillHisto(h_cf_muoneta_loose,muon->Eta());
           FillHisto(h_cf_muonphi_loose, muon->Phi());
	 }
	 if (muon->Pt()>20 && abs(muon->Eta())<2.4 && Muon_isSoft->at(i) == true && Muon_PFIsodBeta04->at(i) < muonisobetacut){
           FillHisto(h_cf_muonpt_soft, muon->Pt());
           FillHisto(h_cf_muoneta_soft,muon->Eta());
           FillHisto(h_cf_muonphi_soft, muon->Phi());
	 }
	 if (muon->Pt()>20 && abs(muon->Eta())<2.4 && Muon_isMedium->at(i) == true && Muon_PFIsodBeta04->at(i) < muonisobetacut){
           FillHisto(h_cf_muonpt_medium, muon->Pt());
           FillHisto(h_cf_muoneta_medium,muon->Eta());
           FillHisto(h_cf_muonphi_medium, muon->Phi());
	 }
	 if (muon->Pt()>20 && abs(muon->Eta())<2.4 && Muon_isHighPt->at(i) == true && Muon_PFIsodBeta04->at(i) < muonisobetacut){
           FillHisto(h_cf_muonpt_highpt, muon->Pt());
           FillHisto(h_cf_muoneta_highpt,muon->Eta());
           FillHisto(h_cf_muonphi_highpt, muon->Phi());
	 }
      }
      if (Muon->GetEntries()>2){
	 leading_muon = (TLorentzVector*) Muon->At(leadpt_index_muon);	 
	 sub_muon = (TLorentzVector*) Muon->At(subpt_index_muon);
	 if (leading_muon->Pt()>20 && sub_muon->Pt()>20 && Muon_Charge->at(leadpt_index_muon)!=Muon_Charge->at(subpt_index_muon)){
           FillHisto(h_cf_muonpt_leading, leading_muon->Pt());
           FillHisto(h_cf_muoneta_leading, leading_muon->Eta());
           FillHisto(h_cf_muonphi_leading, leading_muon->Phi());
           FillHisto(h_cf_muonpt_sub, sub_muon->Pt());
           FillHisto(h_cf_muoneta_sub, sub_muon->Eta());
           FillHisto(h_cf_muonphi_sub, sub_muon->Phi());
	   //cout<< "leading muon pt : "<<leading_muon->Pt()<<" | subleading muon pt : "<<sub_muon->Pt()<<endl; 
	   TightMuon.push_back(*muon);
	 }
      }	 
      //cout<<"dgbmuon"<<endl;
      TLorentzVector* jet;
      TLorentzVector* leading_jet;
      TLorentzVector* sub_jet;
      int leadpt_index_jet=-1;
      float leadpt_jet=0.0;
      int subpt_index_jet=-1;
      float subpt_jet=0.0;
      int number_of_jet=0;
      int number_of_bjet=0;
      int jet_clearing_flag;
      v_jet_idx.clear();
      v_jet_TL.clear();
      v_bjet_idx.clear();
      for(int i = 0; i < Jet->GetEntries(); ++i )
      {
         jet = (TLorentzVector*)Jet->At(i);
         jet_clearing_flag=0;
	 //cout<<"???"<<jet->DeltaR()<<endl;
         //cout << "METtron pT " << met->Pt() << endl;
         for (int j =0;j<TightEle.size();j++){
	   if (jet->DeltaR(TightEle.at(j))<0.4){ jet_clearing_flag =1; break;}
	 }
         if(jet_clearing_flag) continue;
         for (int j =0;j<TightMuon.size();j++){
	   if (jet->DeltaR(TightMuon.at(j))<0.4){ jet_clearing_flag =1; break;}
	 }
         if(jet_clearing_flag) continue;
	 if (Jet_PFId->at(i) > 0 && jet->Pt()>30 && TMath::Abs(jet->Eta())<2.4){
	   v_jet_idx.push_back(i);
	   v_jet_TL.push_back(jet);
	   if (leadpt_jet<jet->Pt()){
	     subpt_index_jet = leadpt_index_jet;
	     subpt_jet = leadpt_jet;
	     leadpt_jet=jet->Pt();
	     leadpt_index_jet=i;
	   }
	   else if (subpt_jet<jet->Pt()) {
	     subpt_index_jet = i;
	     subpt_jet=jet->Pt();
	   }
	   number_of_jet+=1;
           FillHisto(h_cf_jetpt, jet->Pt());
           FillHisto(h_cf_jeteta,jet->Eta());
           FillHisto(h_cf_jetphi,jet->Phi());
	   if (Jet_bDisc->at(i)>0.8484){
	     v_bjet_idx.push_back(i);
	     number_of_bjet+=1;
             FillHisto(h_cf_bjetpt_select, jet->Pt());
             FillHisto(h_cf_bjeteta_select,jet->Eta());
             FillHisto(h_cf_bjetphi_select,jet->Phi());
	   }
	 }
	 else {
           FillHisto(h_cf_nonbjetpt, jet->Pt());
           FillHisto(h_cf_nonbjeteta,jet->Eta());
           FillHisto(h_cf_nonbjetphi,jet->Phi());
	 }
	 
      }

      //double w_pileup = puweight->weight(PileUp_Count_Intime);
      //cout<<"dbgstart"<<endl;
      if(isData==0 && METMUCleanCor->GetEntries()==0) continue;
      if(isData==1 && METMUEGCleanCor->GetEntries()==0) continue;
      //cout<<"dbgMET"<<endl;
      if(number_of_tight_ele<2)continue;
      //cout<<"dbgtightele"<<endl;
      if(third_ele_veto==1)continue;
      //cout<<"dbgthirdele"<<endl;
      if(third_muon_veto==1)continue;
      //cout<<"dbgthirdmuon"<<endl;
      if(number_of_jet<2)continue;
      //cout<<"dbgnjet"<<endl;
      if(number_of_bjet<1)continue;
      //cout<<"dbgnbjet"<<endl;
      leading_ele = (TLorentzVector*) Elec->At(leadpt_index_ele);	 
      sub_ele = (TLorentzVector*) Elec->At(subpt_index_ele);
      leading_jet = (TLorentzVector*) Jet->At(leadpt_index_jet);	 
      sub_jet = (TLorentzVector*) Jet->At(subpt_index_jet);
      if (Elec_Charge->at(leadpt_index_ele)!=Elec_Charge->at(subpt_index_ele)){
	if(leading_ele->Pt()>25&&sub_ele->Pt()>20){
	  mother_ele.SetVect(leading_ele->Vect()+sub_ele->Vect());	  
	  mother_ele.SetE(leading_ele->E()+sub_ele->E());
	  if(mother_ele.M()>20 &&( mother_ele.M()<76 || mother_ele.M()>106)){
	    if (metcleancor->Pt()>40){
	      if(Elec_Charge->at(leadpt_index_ele)>0){
                 Lep = (TLorentzVector*)Elec->At(leadpt_index_ele);
                 AnLep = (TLorentzVector*)Elec->At(subpt_index_ele);
              } else {
                 Lep = (TLorentzVector*)Elec->At(subpt_index_ele);
                 AnLep = (TLorentzVector*)Elec->At(leadpt_index_ele);
	      }	
 	      SetUpKINObs();
	      FillHisto(h_cf_elept_leading,leading_ele->Pt(),w_pileup);	      
	      FillHisto(h_cf_eleeta_leading,leading_ele->Eta(),w_pileup);	      
	      FillHisto(h_cf_elephi_leading,leading_ele->Phi(),w_pileup);
		      
	      FillHisto(h_cf_elept_sub,sub_ele->Pt(),w_pileup);	      
	      FillHisto(h_cf_eleeta_sub,sub_ele->Eta(),w_pileup);	      
	      FillHisto(h_cf_elephi_sub,sub_ele->Phi(),w_pileup);
	      FillHisto(h_cf_ele_invm,mother_ele.M(),w_pileup);	      
	      
	      FillHisto(h_cf_jetpt_leading,leading_jet->Pt(),w_pileup);	      
	      FillHisto(h_cf_jeteta_leading,leading_jet->Eta(),w_pileup);	      
	      FillHisto(h_cf_jetphi_leading,leading_jet->Phi(),w_pileup);
		      
	      FillHisto(h_cf_jetpt_sub,sub_jet->Pt(),w_pileup);	      
	      FillHisto(h_cf_jeteta_sub,sub_jet->Eta(),w_pileup);	      
	      FillHisto(h_cf_jetphi_sub,sub_jet->Phi(),w_pileup);

	      
              FillHisto(h_cf_metpt_evt, metcleancor->Pt(),w_pileup);
              FillHisto(h_cf_metphi_evt, metcleancor->Phi(),w_pileup);
	      FillHisto(h_cf_PVcount,PV_Count);
	      if(isKinSol){
		//cout<<Top->M()<<endl;
                FillHisto(h_cf_KIN_toppt, Top->Pt(), w_pileup);
                FillHisto(h_cf_KIN_topmass, Top->M(), w_pileup);
                FillHisto(h_cf_KIN_topeta, Top->Eta(), w_pileup);
                FillHisto(h_cf_KIN_topphi, Top->Phi(), w_pileup);
                FillHisto(h_cf_KIN_antitoppt, AnTop->Pt(), w_pileup);
                FillHisto(h_cf_KIN_antitopmass, AnTop->M(), w_pileup);
                FillHisto(h_cf_KIN_antitopeta, AnTop->Eta(), w_pileup);
                FillHisto(h_cf_KIN_antitopphi, AnTop->Phi(), w_pileup);

                FillHisto(h_cf_KIN_bjetpt, bJet->Pt(), w_pileup);
                FillHisto(h_cf_KIN_bjeteta, bJet->Eta(), w_pileup);
                FillHisto(h_cf_KIN_bjetphi, bJet->Phi(), w_pileup);
                FillHisto(h_cf_KIN_antibjetpt, AnbJet->Pt(), w_pileup);
                FillHisto(h_cf_KIN_antibjeteta, AnbJet->Eta(), w_pileup);
                FillHisto(h_cf_KIN_antibjetphi, AnbJet->Phi(), w_pileup);

                FillHisto(h_cf_KIN_Nupt, Nu->Pt(), w_pileup);
                FillHisto(h_cf_KIN_Nueta, Nu->Eta(), w_pileup);
                FillHisto(h_cf_KIN_Nuphi, Nu->Phi(), w_pileup);
                FillHisto(h_cf_KIN_antiNupt, AnNu->Pt(), w_pileup);
                FillHisto(h_cf_KIN_antiNueta, AnNu->Eta(), w_pileup);
                FillHisto(h_cf_KIN_antiNuphi, AnNu->Phi(), w_pileup);
                
                FillHisto(h_cf_KIN_leading_elept, leading_ele->Pt(), w_pileup);
                FillHisto(h_cf_KIN_leading_eleeta, leading_ele->Eta(), w_pileup);
                FillHisto(h_cf_KIN_leading_elephi, leading_ele->Phi(), w_pileup);
                FillHisto(h_cf_KIN_subleading_elept, sub_ele->Pt(), w_pileup);
                FillHisto(h_cf_KIN_subleading_eleeta, sub_ele->Eta(), w_pileup);
                FillHisto(h_cf_KIN_subleading_elephi, sub_ele->Phi(), w_pileup);

                FillHisto(h_cf_KIN_mu_invm_ll, mother_ele.M(), w_pileup);

                FillHisto(h_cf_KIN_metpt, met->Pt(), w_pileup);
                FillHisto(h_cf_KIN_metphi, met->Phi(), w_pileup);

                FillHisto(h_cf_KIN_leading_jetpt, leading_jet->Pt(), w_pileup);
                FillHisto(h_cf_KIN_leading_jeteta, leading_jet->Eta(), w_pileup);
                FillHisto(h_cf_KIN_leading_jetphi, leading_jet->Phi(), w_pileup);
                FillHisto(h_cf_KIN_subleading_jetpt, sub_jet->Pt(), w_pileup);
                FillHisto(h_cf_KIN_subleading_jeteta, sub_jet->Eta(), w_pileup);
                FillHisto(h_cf_KIN_subleading_jetphi, sub_jet->Phi(), w_pileup);
	        GetCPViolVar();
	        FillHisto(h_cf_KIN_O1,varO1,w_pileup);
	        FillHisto(h_cf_KIN_O3,varO3,w_pileup);
		
	      }		      
	    }	  
	  }	  
	}
      }
   }//event loop
   std::cout<< "the number of muon eta is 0 evt is "<<muon_eta_0<<std::endl; 
   printf("Total processed number of events: %lld\n", __tot_evt);


}//end Loop function

void ssb_analysis::Start( int genLoopon )
{
   if      ( genLoopon == 0 ){ fout = new TFile(Form("output/%s",outfile),"RECREATE");}
   else if      ( genLoopon == 1 ){ fout = new TFile(Form("output/%s",outfile),"UPDATE");}
   else {cout << "genLoopon error" << endl;}
   fout->cd("");

   TDirectory *dir = gDirectory;
   dir->cd();

   DeclareHistos();

}

void ssb_analysis::DeclareHistos()
{

   /// Test For Systematic All-in-One Code ///
   h_cf_elept          = new TH1D(Form("_h_cf_elept_"),Form("Electron pT"), 1000, 0, 1000); h_cf_elept->Sumw2();
   h_cf_eleeta         = new TH1D(Form("_h_cf_eleeta_"),Form("Electron eta"), 1000, -5, 5); h_cf_eleeta->Sumw2();
   h_cf_elephi         = new TH1D(Form("_h_cf_elephi_"),Form("Electron phi"), 1000, -5 ,5); h_cf_elephi->Sumw2();
   h_cf_elept_select   = new TH1D(Form("_h_cf_elept_select_"),Form("select Electron pT"), 1000, 0, 1000); h_cf_elept_select->Sumw2();
   h_cf_eleeta_select  = new TH1D(Form("_h_cf_eleeta_select_"),Form("select Electron eta"), 1000, -5, 5); h_cf_eleeta_select->Sumw2();
   h_cf_elephi_select  = new TH1D(Form("_h_cf_elephi_select_"),Form("select Electron phi"), 1000, -5 ,5); h_cf_elephi_select->Sumw2();
   h_cf_elept_leading  = new TH1D(Form("_h_cf_elept_leading_"),Form("leading Electron pT"), 1000, 0, 1000); h_cf_elept_leading->Sumw2();
   h_cf_eleeta_leading = new TH1D(Form("_h_cf_eleeta_leading_"),Form("leading Electron eta"), 1000, -5, 5); h_cf_eleeta_leading->Sumw2();
   h_cf_elephi_leading = new TH1D(Form("_h_cf_elephi_leading_"),Form("leading Electron phi"), 1000, -5 ,5); h_cf_elephi_leading->Sumw2();
   h_cf_elept_sub      = new TH1D(Form("_h_cf_elept_sub_leading_"),Form("sub leading Electron pT"), 1000, 0, 1000); h_cf_elept_sub->Sumw2();
   h_cf_eleeta_sub     = new TH1D(Form("_h_cf_eleeta_sub_leading_"),Form("sub leading Electron eta"), 1000, -5, 5); h_cf_eleeta_sub->Sumw2();
   h_cf_elephi_sub     = new TH1D(Form("_h_cf_elephi_sub_leading_"),Form("sub leading Electron phi"), 1000, -5 ,5); h_cf_elephi_sub->Sumw2();
   h_cf_metpt          = new TH1D(Form("_h_cf_metpt_"),Form("MET"), 1000, 0, 1000); h_cf_metpt->Sumw2();
   h_cf_metphi         = new TH1D(Form("_h_cf_metphi_"),Form("MET phi"), 1000, -5, 5); h_cf_metphi->Sumw2();
   h_cf_metpt_select   = new TH1D(Form("_h_cf_metpt_select_"),Form("MET"), 1000, 0, 1000); h_cf_metpt_select->Sumw2();
   h_cf_metphi_select  = new TH1D(Form("_h_cf_metphi_select_"),Form("MET phi"), 1000, -5, 5); h_cf_metphi_select->Sumw2();
   h_cf_metpt_evt      = new TH1D(Form("_h_cf_metpt_evt_"),Form("MET"), 1000, 0, 1000); h_cf_metpt_evt->Sumw2();
   h_cf_metphi_evt     = new TH1D(Form("_h_cf_metphi_evt_"),Form("MET phi"), 1000, -5, 5); h_cf_metphi_evt->Sumw2();
   h_cf_muonpt         = new TH1D(Form("_h_cf_muonpt_"),Form("Muon pT"), 1000, 0, 1000); h_cf_muonpt->Sumw2();
   h_cf_muoneta        = new TH1D(Form("_h_cf_muoneta_"),Form("Muon eta"), 1000, -5, 5); h_cf_muoneta->Sumw2();
   h_cf_muonphi        = new TH1D(Form("_h_cf_muonphi_"),Form("Muon phi"), 1000, -5, 5); h_cf_muonphi->Sumw2();
   h_cf_muonpt_zeroeta = new TH1D(Form("_h_cf_muonpt_zeroeta_"),Form("Muon pT (eta = 0)"), 1000, 0, 1000); h_cf_muonpt_zeroeta->Sumw2();
   h_cf_muoneta_zeroeta= new TH1D(Form("_h_cf_muoneta_zeroeta_"),Form("Muon eta (eta = 0)"), 1000, -5, 5); h_cf_muoneta_zeroeta->Sumw2();
   h_cf_muonphi_zeroeta= new TH1D(Form("_h_cf_muonphi_zeroeta_"),Form("Muon phi (eta = 0)"), 1000, -5, 5); h_cf_muonphi_zeroeta->Sumw2();
   h_cf_muonpt_check   = new TH1D(Form("_h_cf_muonpt_check_"),Form("Muon pT"), 1000, 0, 1000); h_cf_muonpt_check->Sumw2();
   h_cf_muoneta_check  = new TH1D(Form("_h_cf_muoneta_check_"),Form("Muon eta"), 1000, -5, 5); h_cf_muoneta_check->Sumw2();
   h_cf_muonphi_check  = new TH1D(Form("_h_cf_muonphi_check_"),Form("Muon phi"), 1000, -5, 5); h_cf_muonphi_check->Sumw2();
   h_cf_muonpt_select  = new TH1D(Form("_h_cf_muonpt_select_"),Form("select Muon pT"), 1000, 0, 1000); h_cf_muonpt_select->Sumw2();
   h_cf_muoneta_select = new TH1D(Form("_h_cf_muoneta_select_"),Form("select Muon eta"), 1000, -5, 5); h_cf_muoneta_select->Sumw2();
   h_cf_muonphi_select = new TH1D(Form("_h_cf_muonphi_select_"),Form("select Muon phi"), 1000, -5, 5); h_cf_muonphi_select->Sumw2();
   h_cf_muonpt_loose   = new TH1D(Form("_h_cf_muonpt_loose_"),Form("loose Muon pT"), 1000, 0, 1000); h_cf_muonpt_loose->Sumw2();
   h_cf_muoneta_loose  = new TH1D(Form("_h_cf_muoneta_loose_"),Form("loose Muon eta"), 1000, -5, 5); h_cf_muoneta_loose->Sumw2();
   h_cf_muonphi_loose  = new TH1D(Form("_h_cf_muonphi_loose_"),Form("loose Muon phi"), 1000, -5, 5); h_cf_muonphi_loose->Sumw2();
   h_cf_muonpt_medium  = new TH1D(Form("_h_cf_muonpt_medium_"),Form("medium Muon pT"), 1000, 0, 1000); h_cf_muonpt_medium->Sumw2();
   h_cf_muoneta_medium = new TH1D(Form("_h_cf_muoneta_medium_"),Form("medium Muon eta"), 1000, -5, 5); h_cf_muoneta_medium->Sumw2();
   h_cf_muonphi_medium = new TH1D(Form("_h_cf_muonphi_medium_"),Form("medium Muon phi"), 1000, -5, 5); h_cf_muonphi_medium->Sumw2();
   h_cf_muonpt_soft    = new TH1D(Form("_h_cf_muonpt_soft_"),Form("soft Muon pT"), 1000, 0, 1000); h_cf_muonpt_soft->Sumw2();
   h_cf_muoneta_soft   = new TH1D(Form("_h_cf_muoneta_soft_"),Form("soft Muon eta"), 1000, -5, 5); h_cf_muoneta_soft->Sumw2();
   h_cf_muonphi_soft   = new TH1D(Form("_h_cf_muonphi_soft_"),Form("soft Muon phi"), 1000, -5, 5); h_cf_muonphi_soft->Sumw2();
   h_cf_muonpt_highpt  = new TH1D(Form("_h_cf_muonpt_highpt_"),Form("highpt Muon pT"), 1000, 0, 1000); h_cf_muonpt_highpt->Sumw2();
   h_cf_muoneta_highpt = new TH1D(Form("_h_cf_muoneta_highpt_"),Form("highpt Muon eta"), 1000, -5, 5); h_cf_muoneta_highpt->Sumw2();
   h_cf_muonphi_highpt = new TH1D(Form("_h_cf_muonphi_highpt_"),Form("highpt Muon phi"), 1000, -5, 5); h_cf_muonphi_highpt->Sumw2();
   h_cf_muonpt_leading = new TH1D(Form("_h_cf_muonpt_leading_"),Form("leading Muon pT"), 1000, 0, 1000); h_cf_muonpt_leading->Sumw2();
   h_cf_muoneta_leading= new TH1D(Form("_h_cf_muoneta_leading_"),Form("leading Muon eta"), 1000, -5, 5); h_cf_muoneta_leading->Sumw2();
   h_cf_muonphi_leading= new TH1D(Form("_h_cf_muonphi_leading_"),Form("leading Muon phi"), 1000, -5, 5); h_cf_muonphi_leading->Sumw2();
   h_cf_muonpt_sub     = new TH1D(Form("_h_cf_muonpt_sub_leading_"),Form("sub leading Muon pT"), 1000, 0, 1000); h_cf_muonpt_sub->Sumw2();
   h_cf_muoneta_sub    = new TH1D(Form("_h_cf_muoneta_sub_leading_"),Form("sub leading Muon eta"), 1000, -5, 5); h_cf_muoneta_sub->Sumw2();
   h_cf_muonphi_sub    = new TH1D(Form("_h_cf_muonphi_sub_leading_"),Form("sub leading Muon phi"), 1000, -5, 5); h_cf_muonphi_sub->Sumw2();
   h_cf_jetpt          = new TH1D(Form("_h_cf_jetpt_"),Form("Jet pT"), 1000, 0, 1000); h_cf_jetpt->Sumw2();
   h_cf_jeteta         = new TH1D(Form("_h_cf_jeteta_"),Form("Jet eta"), 1000, -5, 5); h_cf_jeteta->Sumw2();
   h_cf_jetphi         = new TH1D(Form("_h_cf_jetphi_"),Form("Jet phi"), 1000, -5, 5); h_cf_jetphi->Sumw2();
   h_cf_bjetpt         = new TH1D(Form("_h_cf_bjetpt_"),Form("BJet pT"), 1000, 0, 1000); h_cf_bjetpt->Sumw2();
   h_cf_bjeteta        = new TH1D(Form("_h_cf_bjeteta_"),Form("BJet eta"), 1000, -5, 5); h_cf_bjeteta->Sumw2();
   h_cf_bjetphi        = new TH1D(Form("_h_cf_bjetphi_"),Form("BJet phi"), 1000, -5, 5); h_cf_bjetphi->Sumw2();
   h_cf_bjetpt_select  = new TH1D(Form("_h_cf_bjetpt_select_"),Form("select BJet pT"), 1000, 0, 1000); h_cf_bjetpt_select->Sumw2();
   h_cf_bjeteta_select = new TH1D(Form("_h_cf_bjeteta_select_"),Form("select BJet eta"), 1000, -5, 5); h_cf_bjeteta_select->Sumw2();
   h_cf_bjetphi_select = new TH1D(Form("_h_cf_bjetphi_select_"),Form("select BJet phi"), 1000, -5, 5); h_cf_bjetphi_select->Sumw2();
   h_cf_nonbjetpt      = new TH1D(Form("_h_cf_non_bjetpt_"),Form("non BJet pT"), 1000, 0, 1000); h_cf_nonbjetpt->Sumw2();
   h_cf_nonbjeteta     = new TH1D(Form("_h_cf_non_bjeteta_"),Form("non BJet eta"), 1000, -5, 5); h_cf_nonbjeteta->Sumw2();
   h_cf_nonbjetphi     = new TH1D(Form("_h_cf_non_bjetphi_"),Form("non BJet phi"), 1000, -5, 5); h_cf_nonbjetphi->Sumw2();
   h_cf_jetpt_leading  = new TH1D(Form("_h_cf_jetpt_leading_"),Form("leading Jet pT"), 1000, 0, 1000); h_cf_jetpt_leading->Sumw2();
   h_cf_jeteta_leading = new TH1D(Form("_h_cf_jeteta_leading_"),Form("leading Jet eta"), 1000, -5, 5); h_cf_jeteta_leading->Sumw2();
   h_cf_jetphi_leading = new TH1D(Form("_h_cf_jetphi_leading_"),Form("leading Jet phi"), 1000, -5 ,5); h_cf_jetphi_leading->Sumw2();
   h_cf_jetpt_sub      = new TH1D(Form("_h_cf_jetpt_sub_leading_"),Form("sub leading Jet pT"), 1000, 0, 1000); h_cf_jetpt_sub->Sumw2();
   h_cf_jeteta_sub     = new TH1D(Form("_h_cf_jeteta_sub_leading_"),Form("sub leading Jet eta"), 1000, -5, 5); h_cf_jeteta_sub->Sumw2();
   h_cf_jetphi_sub     = new TH1D(Form("_h_cf_jetphi_sub_leading_"),Form("sub leading Jet phi"), 1000, -5 ,5); h_cf_jetphi_sub->Sumw2();
   h_cf_ele_invm       = new TH1D(Form("_h_cf_ele_invm_"),"",1000,0,1000); h_cf_ele_invm->Sumw2();   
   h_cf_pu_intime      = new TH1D(Form("_h_cf_pu_intime"),"",100,0,100); h_cf_pu_intime->Sumw2();   
   h_cf_pu_pv          = new TH1D(Form("_h_cf_pu_pv"),"",100,0,100); h_cf_pu_pv->Sumw2();   
   h_cf_pu_interaction = new TH1D(Form("_h_cf_pu_interaction"),"",100,0,100); h_cf_pu_interaction->Sumw2();   
   h_cf_pu_intime_w    = new TH1D(Form("_h_cf_pu_intime_w"),"",100,0,100); h_cf_pu_intime_w->Sumw2();   
   h_cf_pu_pv_w        = new TH1D(Form("_h_cf_pu_pv_w"),"",100,0,100); h_cf_pu_pv_w->Sumw2();   
   h_cf_pu_interaction_w= new TH1D(Form("_h_cf_pu_interaction_w"),"",100,0,100); h_cf_pu_interaction_w->Sumw2();   
   h_cf_PVcount        = new TH1D(Form("_h_cf_pv_count"),"",100,0,100); h_cf_PVcount->Sumw2();   

   h_cf_KIN_toppt       = new TH1D(Form("_h_evt_KIN_toppt_"),Form("Top pT"), 1000, 0, 1000); h_cf_KIN_toppt->Sumw2();
   h_cf_KIN_topmass     = new TH1D(Form("_h_evt_KIN_topmass_"),Form("Top mass"), 1000, 0, 1000); h_cf_KIN_topmass->Sumw2();
   h_cf_KIN_topeta      = new TH1D(Form("_h_evt_KIN_topeta_"),Form("Top eta"), 1000, -5, 5); h_cf_KIN_topeta->Sumw2();
   h_cf_KIN_topphi      = new TH1D(Form("_h_evt_KIN_topphi_"),Form("Top phi"), 1000, -5, 5);
   h_cf_KIN_antitoppt   = new TH1D(Form("_h_evt_KIN_antitoppt_"),Form("AntiTop pT"), 1000, 0, 1000);
   h_cf_KIN_antitopeta  = new TH1D(Form("_h_evt_KIN_antitopeta_"),Form("AntiTop eta"), 1000, -5, 5);
   h_cf_KIN_antitopphi  = new TH1D(Form("_h_evt_KIN_antitopphi_"),Form("AntiTop phi"), 1000, -5, 5);
   h_cf_KIN_antitopmass     = new TH1D(Form("_h_evt_KIN_antitopmass_"),Form("AntiTop mass"), 1000, 0, 1000);

   h_cf_KIN_bjetpt        = new TH1D(Form("_h_evt_KIN_bjetpt_"),Form("Bjet pT "), 1000, 0, 1000);
   h_cf_KIN_bjeteta       = new TH1D(Form("_h_evt_KIN_bjeteta_"),Form("Bjet eta "), 1000, -5, 5);
   h_cf_KIN_bjetphi       = new TH1D(Form("_h_evt_KIN_bjetphi_"),Form("Bjet phi "), 1000, -5, 5);
   h_cf_KIN_antibjetpt    = new TH1D(Form("_h_evt_KIN_antibjetpt_"),Form("Bjet pT "), 1000, 0, 1000);
   h_cf_KIN_antibjeteta   = new TH1D(Form("_h_evt_KIN_antibjeteta_"),Form("Bjet eta "), 1000, -5, 5);
   h_cf_KIN_antibjetphi   = new TH1D(Form("_h_evt_KIN_antibjetphi_"),Form("Bjet phi "), 1000, -5, 5);

   h_cf_KIN_Nupt        = new TH1D(Form("_h_evt_KIN_Nupt_"),Form("Nu pT "), 1000, 0, 1000);
   h_cf_KIN_Nueta       = new TH1D(Form("_h_evt_KIN_Nueta_"),Form("Nu eta "), 1000, -5, 5);
   h_cf_KIN_Nuphi       = new TH1D(Form("_h_evt_KIN_Nuphi_"),Form("Nu phi "), 1000, -5, 5);
   h_cf_KIN_antiNupt    = new TH1D(Form("_h_evt_KIN_AntiNupt_"),Form("Nu pT "), 1000, 0, 1000);
   h_cf_KIN_antiNueta   = new TH1D(Form("_h_evt_KIN_AntiNueta_"),Form("Nu eta "), 1000, -5, 5);
   h_cf_KIN_antiNuphi   = new TH1D(Form("_h_evt_KIN_AntiNuphi_"),Form("Nu phi "), 1000, -5, 5);

   h_cf_KIN_leading_elept    = new TH1D(Form("_h_evt_KIN_leading_elept_"),Form("Elec pT (af cut)"), 500, 0, 500);
   h_cf_KIN_leading_eleeta   = new TH1D(Form("_h_evt_KIN_leading_eleeta_"),Form("Elec eta (af cut)"), 1000, -5, 5); 
   h_cf_KIN_leading_elephi   = new TH1D(Form("_h_evt_KIN_leading_elephi_"),Form("Elec phi (af cut)"), 1000, -5, 5);
   h_cf_KIN_subleading_elept    = new TH1D(Form("_h_evt_KIN_subleading_elept_"),Form("Elec pT (af cut)"), 500, 0, 500);
   h_cf_KIN_subleading_eleeta   = new TH1D(Form("_h_evt_KIN_subleading_eleeta_"),Form("Elec eta (af cut)"), 1000, -5, 5); 
   h_cf_KIN_subleading_elephi   = new TH1D(Form("_h_evt_KIN_subleading_elephi_"),Form("Elec phi (af cut)"), 1000, -5, 5);

   h_cf_KIN_mu_invm_ll      = new TH1D(Form("_h_evt_KIN_mu_invm_ll_"),Form("Invariant Mass of Dimuon"), 500, 0, 500);

   h_cf_KIN_metpt    = new TH1D(Form("_h_evt_KIN_metpt_"), Form("MET pT (evt)"), 1000, 0, 500);
   h_cf_KIN_metphi   = new TH1D(Form("_h_evt_KIN_metphi_"), Form("MET phi (evt)"), 1000, -5, 5);

   h_cf_KIN_leading_jetpt       = new TH1D(Form("_h_evt_KIN_leading_jetpt_"),Form("Leading jet pT"), 1000, 0, 1000);
   h_cf_KIN_leading_jeteta      = new TH1D(Form("_h_evt_KIN_leading_jeteta_"),Form("Leading jet eta"), 1000, -5, 5);
   h_cf_KIN_leading_jetphi      = new TH1D(Form("_h_evt_KIN_leading_jetphi_"),Form("Leading jet phi"), 1000, -5, 5);
   h_cf_KIN_subleading_jetpt    = new TH1D(Form("_h_evt_KIN_subleading_jetpt_"),Form("Subleading jet pT"), 1000, 0, 1000);
   h_cf_KIN_subleading_jeteta   = new TH1D(Form("_h_evt_KIN_subleading_jeteta_"),Form("Subleading jet eta"), 1000, -5, 5);
   h_cf_KIN_subleading_jetphi   = new TH1D(Form("_h_evt_KIN_subleading_jetphi_"),Form("Subleading jet phi"), 1000, -5, 5);

   h_cf_KIN_O1       = new TH1D(Form("_h_evt_KIN_O1_"),Form("O1"), 1000, -5, 5); h_cf_KIN_O3->Sumw2();
   h_cf_KIN_O3       = new TH1D(Form("_h_evt_KIN_O3_"),Form("O3"), 1000, -5, 5); h_cf_KIN_O1->Sumw2();

}


void ssb_analysis::End()
{
   fout->Write();
   fout->Close();
}

void ssb_analysis::SetOutputFileName(char *outname)
{   
   outfile = outname;
}


void ssb_analysis::SetUpKINObs()
{
   isKinSol=false;
   v_leptons_VLV.clear();
   v_jets_VLV.clear();
   v_bjets_VLV.clear();
   v_lepidx_KIN.clear();
   v_anlepidx_KIN.clear();
   v_jetidx_KIN.clear();
   v_bjetidx_KIN.clear();
   v_btagging_KIN.clear();
   /// lepton ///
   v_leptons_VLV.push_back(common::TLVtoLV(*Lep));
   v_lepidx_KIN.push_back(0);
   v_leptons_VLV.push_back(common::TLVtoLV(*AnLep));
   v_anlepidx_KIN.push_back(1);
   //cout << "jet info. for kinematic solver"  << endl;
   //cout << v_jet_idx.size() << " " << v_jet_TL.size() << endl;
   const KinematicReconstruction* kinematicReconstruction(0); 
   kinematicReconstruction = new KinematicReconstruction(1, true);
   
   const LV met_LV = common::TLVtoLV(*Met);
   
   for (int ijet = 0; ijet < v_jet_idx.size(); ++ijet)
   {
      int idx_jet = v_jet_idx[ijet];
      if(Jet_bDisc->at(idx_jet)  > 0.8484 ){
         v_bjets_VLV.push_back(common::TLVtoLV(*v_jet_TL[ijet]));
         //v_bjetidx_KIN.push_back(v_jet_idx[ijet]);
         v_bjetidx_KIN.push_back(ijet);
      }
      v_jets_VLV.push_back(common::TLVtoLV(*v_jet_TL[ijet]));
      v_jetidx_KIN.push_back(ijet);
      v_btagging_KIN.push_back(Jet_bDisc->at(idx_jet));
    }
    KinematicReconstructionSolutions kinematicReconstructionSolutions = kinematicReconstruction->solutions(v_lepidx_KIN,v_anlepidx_KIN, v_jetidx_KIN,  v_bjetidx_KIN,  v_leptons_VLV, v_jets_VLV, v_btagging_KIN, met_LV);
    //cout << "Num Sol : " << kinematicReconstructionSolutions.numberOfSolutions() << endl;
    //cout << "MET ? " << met_LV.pt() << endl;
    //cout << "MET Pt ? " << Met->Pt() << endl;
    if (kinematicReconstructionSolutions.numberOfSolutions())
    {
       isKinSol= true;
       LV top1 = kinematicReconstructionSolutions.solution().top();
       LV top2 = kinematicReconstructionSolutions.solution().antiTop();
       LV bjet1 = kinematicReconstructionSolutions.solution().bjet();
       LV bjet2 = kinematicReconstructionSolutions.solution().antiBjet();
       LV neutrino1 = kinematicReconstructionSolutions.solution().neutrino();
       LV neutrino2 = kinematicReconstructionSolutions.solution().antiNeutrino();
       //kinematicReconstructionSolutions.solution().print();
       //Top = new TLorentzVector(common::LVtoTLV(top1));   
       Top       = new TLorentzVector(common::LVtoTLV(top1));
       AnTop     = new TLorentzVector(common::LVtoTLV(top2));
       bJet      = new TLorentzVector(common::LVtoTLV(bjet1));
       AnbJet    = new TLorentzVector(common::LVtoTLV(bjet2));
       Nu        = new TLorentzVector(common::LVtoTLV(neutrino1));
       AnNu      = new TLorentzVector(common::LVtoTLV(neutrino2));
   
       (*W1)        = (*Lep)  + (*AnNu);
       (*W2)        = (*AnLep)  + (*Nu);
   
    }
    //delete 
    delete kinematicReconstruction;
}
void ssb_analysis::GetCPViolVar(){
    varO1=ssbcpviol->getO1Vari( Top, AnTop, AnLep, Lep);
    varO3=ssbcpviol->getO3Vari( Top, AnTop, bJet, AnbJet);
}
