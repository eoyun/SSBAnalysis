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
      void SetOutputFileName(char *outname);
      void DeclareHistos();

   private:
      //put variables that you want
      char *outfile;
      TFile *fout;
      edm::LumiReWeighting  *puweight;
      double luminosity;
      double crosssection;
      // vector for ChargeMisId
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
      TH1D *h_cf_PVcount;
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
