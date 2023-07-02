#ifndef _ssb_analysis_

#define _ssb_analysis_
  
#include <set>
#include <string>
#include <fstream>
#include <cassert>
#include <sstream>
#include <map>
#include <vector>
#include <string>
#include <TH1.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TGraphErrors.h>
#include "TLorentzVector.h"
#include "TEnv.h"
  
#include "./../analysis/SSBTree.h"

// Textreader
#include "./../TextReader/TextReader.hpp"

#include "./../interface/ssb_eff.hpp"
#include "./../interface/ssb_cpviol.hpp"
#include "./../kinsol/TtFullLepKinSolver.hpp"
#include "./../interface/BTagCalibrationStandalone.h"
#include "./../interface/LumiReWeighting.h"

#include "./../roccor.2016.v3/RoccoR.h"

#include "./../KinSolv/analysisUtils.h"
#include "./../KinSolv/KinematicReconstruction.h"
#include "./../KinSolv/KinematicReconstructionSolution.h"

using namespace std;

class ssb_analysis : public SSBTree 
{
   public:
      //declare functions
      ssb_analysis(TTree *tree=0);
      virtual ~ssb_analysis();

      //basic frame
      virtual void Loop( char *logfile ,string sampleName);
      void Start( int genLoopon );
      void End();
      
      //user define functions
      void TLVInitial();
      void SetUpKINObs();
      void SetOutputFileName(char *outname);
      void DeclareHistos();
      vector<int> v_jet_idx;
      vector<int> v_bjet_idx;
      vector<TLorentzVector*> v_jet_TL;
      void GetCPViolVar();

      void SetUpKINObsSyst(vector<int> v_jetsys_idx, vector<TLorentzVector*> v_jetsys_TL, TLorentzVector* metsys);
      bool isKinSol;
      VLV v_leptons_VLV; 
      VLV v_jets_VLV; 
      VLV v_bjets_VLV; 
      vector<int> v_lepidx_KIN; 
      vector<int> v_anlepidx_KIN; 
      vector<int> v_jetidx_KIN;
      vector<int> v_bjetidx_KIN; 
      vector<double> v_btagging_KIN; 
      SSBCPViol* ssbcpviol;
      double varO1,varO3;
      

   private:
      //put variables that you want
      char *outfile;
      TFile *fout;
      edm::LumiReWeighting  *puweight;
      double luminosity;
      double crosssection;
      // vector for ChargeMisId
      TLorentzVector *Lep, *AnLep, *Met;
      TLorentzVector *W1, *W2;
      TLorentzVector *Top, *AnTop, *bJet, *AnbJet, *Nu, *AnNu;

   public:

      //declare histograms

      // Cut Flow
      TH1D *h_cf_elept; 
      TH1D *h_cf_eleeta; 
      TH1D *h_cf_elephi; 
      TH1D *h_cf_elept_select; 
      TH1D *h_cf_eleeta_select; 
      TH1D *h_cf_elephi_select; 
      TH1D *h_cf_elept_leading; 
      TH1D *h_cf_eleeta_leading; 
      TH1D *h_cf_elephi_leading; 
      TH1D *h_cf_elept_sub; 
      TH1D *h_cf_eleeta_sub;
      TH1D *h_cf_elephi_sub; 
      TH1D *h_cf_metpt; 
      TH1D *h_cf_metphi; 
      TH1D *h_cf_metpt_evt; 
      TH1D *h_cf_metphi_evt; 
      TH1D *h_cf_metpt_select; 
      TH1D *h_cf_metphi_select; 
      TH1D *h_cf_muonpt; 
      TH1D *h_cf_muoneta; 
      TH1D *h_cf_muonphi; 
      TH1D *h_cf_muonpt_zeroeta; 
      TH1D *h_cf_muoneta_zeroeta; 
      TH1D *h_cf_muonphi_zeroeta; 
      TH1D *h_cf_muonpt_check; 
      TH1D *h_cf_muoneta_check; 
      TH1D *h_cf_muonphi_check; 
      TH1D *h_cf_muonpt_select; 
      TH1D *h_cf_muoneta_select; 
      TH1D *h_cf_muonphi_select; 
      TH1D *h_cf_muonpt_medium; 
      TH1D *h_cf_muoneta_medium; 
      TH1D *h_cf_muonphi_medium; 
      TH1D *h_cf_muonpt_soft; 
      TH1D *h_cf_muoneta_soft; 
      TH1D *h_cf_muonphi_soft; 
      TH1D *h_cf_muonpt_loose; 
      TH1D *h_cf_muoneta_loose; 
      TH1D *h_cf_muonphi_loose; 
      TH1D *h_cf_muonpt_highpt; 
      TH1D *h_cf_muoneta_highpt; 
      TH1D *h_cf_muonphi_highpt; 
      TH1D *h_cf_muonpt_leading; 
      TH1D *h_cf_muoneta_leading; 
      TH1D *h_cf_muonphi_leading; 
      TH1D *h_cf_muonpt_sub; 
      TH1D *h_cf_muoneta_sub; 
      TH1D *h_cf_muonphi_sub; 
      TH1D *h_cf_jetpt; 
      TH1D *h_cf_jeteta; 
      TH1D *h_cf_jetphi; 
      TH1D *h_cf_bjetpt; 
      TH1D *h_cf_bjeteta; 
      TH1D *h_cf_bjetphi; 
      TH1D *h_cf_bjetpt_select; 
      TH1D *h_cf_bjeteta_select; 
      TH1D *h_cf_bjetphi_select; 
      TH1D *h_cf_nonbjetpt; 
      TH1D *h_cf_nonbjeteta; 
      TH1D *h_cf_nonbjetphi; 
      TH1D *h_cf_jetpt_leading; 
      TH1D *h_cf_jeteta_leading; 
      TH1D *h_cf_jetphi_leading; 
      TH1D *h_cf_jetpt_sub; 
      TH1D *h_cf_jeteta_sub;
      TH1D *h_cf_jetphi_sub; 
      TH1D *h_cf_ele_invm; 
      TH1D *h_cf_pu_intime;
      TH1D *h_cf_pu_pv;
      TH1D *h_cf_pu_interaction;
      TH1D *h_cf_pu_intime_w;
      TH1D *h_cf_pu_pv_w;
      TH1D *h_cf_pu_interaction_w;
      TH1D *h_cf_PVcount;
      
      TH1D *h_cf_KIN_toppt;
      TH1D *h_cf_KIN_topmass;
      TH1D *h_cf_KIN_topeta;
      TH1D *h_cf_KIN_topphi;
      TH1D *h_cf_KIN_antitoppt;
      TH1D *h_cf_KIN_antitopmass;
      TH1D *h_cf_KIN_antitopeta;
      TH1D *h_cf_KIN_antitopphi;

      TH1D *h_cf_KIN_W1pt;
      TH1D *h_cf_KIN_W1eta;
      TH1D *h_cf_KIN_W1phi;
      TH1D *h_cf_KIN_W2toppt;
      TH1D *h_cf_KIN_W2topeta;
      TH1D *h_cf_KIN_W2topphi;

      TH1D *h_cf_KIN_bjetpt;
      TH1D *h_cf_KIN_bjeteta;
      TH1D *h_cf_KIN_bjetphi;
      TH1D *h_cf_KIN_antibjetpt;
      TH1D *h_cf_KIN_antibjeteta;
      TH1D *h_cf_KIN_antibjetphi;

      TH1D *h_cf_KIN_Nupt;
      TH1D *h_cf_KIN_Nueta;
      TH1D *h_cf_KIN_Nuphi;
      TH1D *h_cf_KIN_antiNupt;
      TH1D *h_cf_KIN_antiNueta;
      TH1D *h_cf_KIN_antiNuphi;

      TH1D *h_cf_KIN_leading_elept;
      TH1D *h_cf_KIN_leading_eleeta;
      TH1D *h_cf_KIN_leading_elephi;
      TH1D *h_cf_KIN_subleading_elept;
      TH1D *h_cf_KIN_subleading_eleeta;
      TH1D *h_cf_KIN_subleading_elephi;

      TH1D *h_cf_KIN_mu_invm_ll;

      TH1D *h_cf_KIN_metpt;
      TH1D *h_cf_KIN_metphi;

      TH1D *h_cf_KIN_leading_jetpt;
      TH1D *h_cf_KIN_leading_jeteta;
      TH1D *h_cf_KIN_leading_jetphi;
      TH1D *h_cf_KIN_subleading_jetpt;
      TH1D *h_cf_KIN_subleading_jeteta;
      TH1D *h_cf_KIN_subleading_jetphi;
      
      TH1D *h_cf_KIN_O1;
      TH1D *h_cf_KIN_O3;
};
#endif

#ifdef ssb_analysis_cxx

ssb_analysis::ssb_analysis(TTree *tree)
{
   if (tree == 0)
   {
      printf("ERROR: Can't find any input tree.\n");
   }
   Init(tree);
   puweight = new edm::LumiReWeighting("./pileuInfo/MC_Moriond.root","./pileuInfo/PU_2016_69p2_36000_XSecCentral.root","pileup","pileup");
   
}

ssb_analysis::~ssb_analysis()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
   delete fout;
   delete puweight;
}

#endif
