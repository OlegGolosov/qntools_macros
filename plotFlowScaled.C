bool verbose=false;

#include "plotFlow.h"

void plotFlowScaled(const char *inputFileName="scaled.corr.root", const char *outputFileName="scaled.plots.root", int pBeam=13, const char *ref="psd1_psd2_psd3", const char *comp="X")
{
  cout << "START!\n";
  auto mp=0.938, mNucl=0.931;
  auto eBeam=sqrt(pBeam*pBeam+mp*mp);
  auto yBeam=0.25*log((eBeam+pBeam)/(eBeam-pBeam));
  auto sqrt2MoverT=sqrt(2*mNucl/eBeam);
  auto ut0coefProton=sqrt2MoverT/0.938;
  auto ut0coefPion=sqrt2MoverT/0.140;
  printf("pBeam = %i \t yBeam = %f \n", pBeam, yBeam);
  cout << "ut0coefProton = " << ut0coefProton << "\tut0coefPion = " << ut0coefPion << endl;

  cout << "\nReading input...\n";
  TFile inputFile(inputFileName);
  for (const auto &k:*inputFile.GetListOfKeys()) 
    corrMap.emplace(k->GetName(), *(Qn::DataContainerStatCollect*)inputFile.Get(k->GetName()));  

  cout << "\nProcessing...\n";
  Qn::AxisD centAxOptimal ("centrality", {0,5,10,15,25,35,45,60,80,100});
  R1_3S("R1_psd1_3S_X", "psd1_RECENTERED_psd3_RECENTERED_XX_centrality", "psd1_RECENTERED_psd2_RECENTERED_XX_centrality", "psd2_RECENTERED_psd3_RECENTERED_XX_centrality", centAxOptimal); 
  R1_3S("R1_psd1_3S_Y", "psd1_RECENTERED_psd3_RECENTERED_YY_centrality", "psd1_RECENTERED_psd2_RECENTERED_YY_centrality", "psd2_RECENTERED_psd3_RECENTERED_YY_centrality", centAxOptimal); 
  R1_3S("R1_psd2_3S_X", "psd1_RECENTERED_psd2_RECENTERED_XX_centrality", "psd2_RECENTERED_psd3_RECENTERED_XX_centrality", "psd1_RECENTERED_psd3_RECENTERED_XX_centrality", centAxOptimal); 
  R1_3S("R1_psd2_3S_Y", "psd1_RECENTERED_psd2_RECENTERED_YY_centrality", "psd2_RECENTERED_psd3_RECENTERED_YY_centrality", "psd1_RECENTERED_psd3_RECENTERED_YY_centrality", centAxOptimal); 
  R1_3S("R1_psd3_3S_X", "psd2_RECENTERED_psd3_RECENTERED_XX_centrality", "psd1_RECENTERED_psd3_RECENTERED_XX_centrality", "psd1_RECENTERED_psd2_RECENTERED_XX_centrality", centAxOptimal); 
  R1_3S("R1_psd3_3S_Y", "psd2_RECENTERED_psd3_RECENTERED_YY_centrality", "psd1_RECENTERED_psd3_RECENTERED_YY_centrality", "psd1_RECENTERED_psd2_RECENTERED_YY_centrality", centAxOptimal); 

  string res;
  string refPart4S;
  float resScale4S;
  if (pBeam==13)
  {
    rebin("proton_RESCALED_psd._RECENTERED_(XX|YY|XY|YX)_centrality", {"trY0", {0.5, 0.7}}, "4S");
    rebin("proton_RESCALED_psd._RECENTERED_(XX|YY|XY|YX)_centrality_rebin_trY0_4S", {"trUt0", {0.2, 1}}, "4S");
    project("proton_RESCALED_psd._RECENTERED_(XX|YY|XY|YX)_centrality_rebin_trY0_4S_rebin_trUt0_4S", {"centrality"});
    res="4S";
    refPart4S="proton";
    resScale4S=1;
  }
  if (pBeam==30 || pBeam==41)
  {
    rebin("pionpos_RESCALED_psd._RECENTERED_(XX|YY|XY|YX)_centrality", {"trY0", {0.2, 0.5}}, "4S");
    rebin("pionpos_RESCALED_psd._RECENTERED_(XX|YY|XY|YX)_centrality_rebin_trY0_4S", {"trUt0", {1.2, 1.8}}, "4S");
    project("pionpos_RESCALED_psd._RECENTERED_(XX|YY|XY|YX)_centrality_rebin_trY0_4S_rebin_trUt0_4S", {"centrality"});
    res="3S";
    refPart4S="pionpos";
    resScale4S=-1;
  }
  R1_3S("R1_tpc_3S_X", refPart4S+"_RESCALED_psd1_RECENTERED_XX_centrality_rebin_trY0_4S_rebin_trUt0_4S_proj_centrality", refPart4S+"_RESCALED_psd3_RECENTERED_XX_centrality_rebin_trY0_4S_rebin_trUt0_4S_proj_centrality", "psd1_RECENTERED_psd3_RECENTERED_XX_centrality", centAxOptimal);
  R1_3S("R1_tpc_3S_Y", refPart4S+"_RESCALED_psd1_RECENTERED_YY_centrality_rebin_trY0_4S_rebin_trUt0_4S_proj_centrality", refPart4S+"_RESCALED_psd3_RECENTERED_YY_centrality_rebin_trY0_4S_rebin_trUt0_4S_proj_centrality", "psd1_RECENTERED_psd3_RECENTERED_YY_centrality", centAxOptimal);
  R1_3S("R1_psd1_4S_X", "psd1_RECENTERED_psd3_RECENTERED_XX_centrality", refPart4S+"_RESCALED_psd1_RECENTERED_XX_centrality_rebin_trY0_4S_rebin_trUt0_4S_proj_centrality", refPart4S+"_RESCALED_psd3_RECENTERED_XX_centrality_rebin_trY0_4S_rebin_trUt0_4S_proj_centrality", centAxOptimal); 
  R1_3S("R1_psd1_4S_Y", "psd1_RECENTERED_psd3_RECENTERED_YY_centrality", refPart4S+"_RESCALED_psd1_RECENTERED_YY_centrality_rebin_trY0_4S_rebin_trUt0_4S_proj_centrality", refPart4S+"_RESCALED_psd3_RECENTERED_YY_centrality_rebin_trY0_4S_rebin_trUt0_4S_proj_centrality", centAxOptimal); 
  R1_3S("R1_psd3_4S_X", "psd1_RECENTERED_psd3_RECENTERED_XX_centrality", refPart4S+"_RESCALED_psd3_RECENTERED_XX_centrality_rebin_trY0_4S_rebin_trUt0_4S_proj_centrality", refPart4S+"_RESCALED_psd1_RECENTERED_XX_centrality_rebin_trY0_4S_rebin_trUt0_4S_proj_centrality", centAxOptimal); 
  R1_3S("R1_psd3_4S_Y", "psd1_RECENTERED_psd3_RECENTERED_YY_centrality", refPart4S+"_RESCALED_psd3_RECENTERED_YY_centrality_rebin_trY0_4S_rebin_trUt0_4S_proj_centrality", refPart4S+"_RESCALED_psd1_RECENTERED_YY_centrality_rebin_trY0_4S_rebin_trUt0_4S_proj_centrality", centAxOptimal); 
  R1_4S("R1_psd2_4S_X", refPart4S+"_RESCALED_psd2_RECENTERED_XX_centrality_rebin_trY0_4S_rebin_trUt0_4S_proj_centrality", "R1_tpc_3S_X", resScale4S, centAxOptimal); 
  R1_4S("R1_psd2_4S_Y", refPart4S+"_RESCALED_psd2_RECENTERED_YY_centrality_rebin_trY0_4S_rebin_trUt0_4S_proj_centrality", "R1_tpc_3S_Y", resScale4S, centAxOptimal); 

  rebin("(pion|proton).*", centAxOptimal, "optimal");
  for (auto &part_:{"pionneg", "pionpos", "proton"})
    for (auto &ref_:{"psd1", "psd2", "psd3"})
      for (auto &comp_:{"X", "Y"})
        v1(Form("v1_%s_%s_%s", part_, ref_, comp_), 
           Form("%s_RESCALED_%s_RECENTERED_%s%s_centrality_rebin_centrality_optimal", part_, ref_, comp_, comp_), 
           Form("R1_%s_%s_%s", ref_, res.c_str(), comp_)); 

  merge("v1_proton_psd1_(X|Y)");
  merge("v1_proton_psd2_(X|Y)");
  merge("v1_proton_psd3_(X|Y)");
  merge("v1_proton_(psd1|psd2|psd3)_X");
  merge("v1_proton_(psd1|psd2|psd3)_Y");
  merge("v1_proton_(psd1|psd2|psd3)_(X|Y)");

  merge("v1_pionneg_psd1_(X|Y)");
  merge("v1_pionneg_psd2_(X|Y)");
  merge("v1_pionneg_psd3_(X|Y)");
  merge("v1_pionneg_(psd1|psd2|psd3)_X");
  merge("v1_pionneg_(psd1|psd2|psd3)_Y");
  merge("v1_pionneg_(psd1|psd2|psd3)_(X|Y)");

  merge("v1_pionpos_psd1_(X|Y)");
  merge("v1_pionpos_psd2_(X|Y)");
  merge("v1_pionpos_psd3_(X|Y)");
  merge("v1_pionpos_(psd1|psd2|psd3)_X");
  merge("v1_pionpos_(psd1|psd2|psd3)_Y");
  merge("v1_pionpos_(psd1|psd2|psd3)_(X|Y)");

  vector<double> 
    ptBins_proton={0.0, 0.4, 1.2},
    yBins_proton={0.4, 0.8},
    ptBins_pion={0.0, 1, 4.4},
    yBins_pion={0.0, 1.7};

  rebin("v1_proton_(.*)_(X|Y|X_Y)", {"trY0", yBins_proton},   "proton");
  rebin("v1_proton_(.*)_(X|Y|X_Y)", {"trUt0", ptBins_proton}, "proton");
  rebin("v1_pionpos_(.*)_(X|Y|X_Y)", {"trY0", yBins_pion}, "pionpos");
  rebin("v1_pionpos_(.*)_(X|Y|X_Y)", {"trUt0", ptBins_pion}, "pionpos");
  rebin("v1_pionneg_(.*)_(X|Y|X_Y)", {"trY0", yBins_pion}, "pionneg");
  rebin("v1_pionneg_(.*)_(X|Y|X_Y)", {"trUt0", ptBins_pion}, "pionneg");

  vector<double> centBinsWide={0,10,25,45,60};
  rebin("v1_.*_rebin_tr.*", {"centrality", centBinsWide}, "wide");
  unfold("v1.*rebin_tr.*", "centrality");

  unfold("v1_.*_rebin_trY0_.*_rebin_centrality_wide_unfold_centrality.*", "trY0", "%1.1f");
  unfold("v1_.*_rebin_trUt0_.*_rebin_centrality_wide_unfold_centrality.*", "trUt0", "%1.1f");
  unfold("v1_.*_rebin_trUt0_(pionneg|pionpos|proton)_unfold_centrality.*", "trUt0", "%1.1f");
  
  project("v1_.*rebin_centrality_wide.*_unfold_trY0.*", {"trUt0"});
  project("v1_.*_unfold_trUt0.*", {"trY0"});
  
  auto fpol3=new TF1("fpol3","[1]*x+[2]*x*x*x", 0, 0.6);
  fit("v1_.*_proj_trY0", "fpol3");
  for(int i=1; i<ptBins_pion.size();i++)
  {
    string s_pt=Form("%.1f-%.1f", ptBins_pion.at(i-1), ptBins_pion.at(i));
    const char *pt=s_pt.c_str();
    projectFitParameters(Form("v1_pionneg_%s_X_rebin_trUt0_pionneg_rebin_centrality_wide_unfold_(centrality)_(.*)-(.*)_unfold_trUt0_%s_proj_trY0_(fit_fpol3)", ref, pt));
    projectFitParameters(Form("v1_pionneg_%s_X_rebin_trUt0_pionneg_unfold_(centrality)_(.*)-(.*)_unfold_trUt0_%s_proj_trY0_(fit_fpol3)", ref, pt));
    projectFitParameters(Form("v1_pionpos_%s_X_rebin_trUt0_pionpos_rebin_centrality_wide_unfold_(centrality)_(.*)-(.*)_unfold_trUt0_%s_proj_trY0_(fit_fpol3)", ref, pt));
    projectFitParameters(Form("v1_pionpos_%s_X_rebin_trUt0_pionpos_unfold_(centrality)_(.*)-(.*)_unfold_trUt0_%s_proj_trY0_(fit_fpol3)", ref, pt));
    projectFitParameters(Form("v1_pionneg_%s_Y_rebin_trUt0_pionneg_rebin_centrality_wide_unfold_(centrality)_(.*)-(.*)_unfold_trUt0_%s_proj_trY0_(fit_fpol3)", ref, pt));
    projectFitParameters(Form("v1_pionneg_%s_Y_rebin_trUt0_pionneg_unfold_(centrality)_(.*)-(.*)_unfold_trUt0_%s_proj_trY0_(fit_fpol3)", ref, pt));
    projectFitParameters(Form("v1_pionpos_%s_Y_rebin_trUt0_pionpos_rebin_centrality_wide_unfold_(centrality)_(.*)-(.*)_unfold_trUt0_%s_proj_trY0_(fit_fpol3)", ref, pt));
    projectFitParameters(Form("v1_pionpos_%s_Y_rebin_trUt0_pionpos_unfold_(centrality)_(.*)-(.*)_unfold_trUt0_%s_proj_trY0_(fit_fpol3)", ref, pt));
  }
  for(int i=1; i<ptBins_proton.size();i++)
  {
    string s_pt=Form("%.1f-%.1f", ptBins_proton.at(i-1), ptBins_proton.at(i));
    const char *pt=s_pt.c_str();
    projectFitParameters(Form("v1_proton_%s_X_rebin_trUt0_proton_unfold_(centrality)_(.*)-(.*)_unfold_trUt0_%s_proj_trY0_(fit_fpol3)", ref, pt));
    projectFitParameters(Form("v1_proton_%s_X_rebin_trUt0_proton_rebin_centrality_wide_unfold_(centrality)_(.*)-(.*)_unfold_trUt0_%s_proj_trY0_(fit_fpol3)", ref, pt));
    projectFitParameters(Form("v1_proton_%s_Y_rebin_trUt0_proton_unfold_(centrality)_(.*)-(.*)_unfold_trUt0_%s_proj_trY0_(fit_fpol3)", ref, pt));
    projectFitParameters(Form("v1_proton_%s_Y_rebin_trUt0_proton_rebin_centrality_wide_unfold_(centrality)_(.*)-(.*)_unfold_trUt0_%s_proj_trY0_(fit_fpol3)", ref, pt));
  }

  TFile outputFile(outputFileName, "recreate");
  cout << "\nPlotting...\n";
//  writeMG({"psd1_RECENTERED_psd2_RECENTERED_(XX|YY|XY|YX)_centrality"});
//  writeMG({"psd2_RECENTERED_psd3_RECENTERED_(XX|YY|XY|YX)_centrality"});
//  writeMG({"psd1_RECENTERED_psd3_RECENTERED_(XX|YY|XY|YX)_centrality"});
//  writeMG({refPart4S+"_RESCALED_psd1_RECENTERED_(XX|YY|XY|YX)_centrality_rebin_trY0_4S_rebin_trUt0_4S_proj_centrality"});
//  writeMG({refPart4S+"_RESCALED_psd3_RECENTERED_(XX|YY|XY|YX)_centrality_rebin_trY0_4S_rebin_trUt0_4S_proj_centrality"});
//  writeMG({"R1_(psd.)_(3S)_(X)", "R1_(psd.)_(3S)_(Y)"});
//  writeMG({"R1_(psd.)_(4S)_(X)", "R1_(psd.)_(4S)_(Y)"});
//  writeMG({"R1_tpc_3S_(X)", "R1_tpc_3S_(Y)"});
//  writeMG({"R1_(psd.)_(4S)_(X)", "R1_(psd.)_(4S)_(Y)", "R1_(psd.)_(3S)_(X)", "R1_(psd.)_(3S)_(Y)"});

  for (int i=1; i<ptBins_pion.size(); i++)
  {
    string s_pt=Form("%.1f-%.1f", ptBins_pion.at(i-1), ptBins_pion.at(i));
    const char *pt=s_pt.c_str();
    writeMG({Form("v1_(pionneg)_%s_(X)_rebin_trUt0_pionneg_rebin_centrality_wide_unfold_centrality_(.*)_unfold_trUt0_%s_proj_trY0", ref, pt),
             Form("v1_(pionpos)_%s_(X)_rebin_trUt0_pionpos_rebin_centrality_wide_unfold_centrality_(.*)_unfold_trUt0_%s_proj_trY0", ref, pt),
             Form("v1_(pionneg)_%s_(Y)_rebin_trUt0_pionneg_rebin_centrality_wide_unfold_centrality_(.*)_unfold_trUt0_%s_proj_trY0", ref, pt),
             Form("v1_(pionpos)_%s_(Y)_rebin_trUt0_pionpos_rebin_centrality_wide_unfold_centrality_(.*)_unfold_trUt0_%s_proj_trY0", ref, pt),
	    });
    writeMG({Form("v1_(pionneg)_%s_X_rebin_trUt0_pionneg_rebin_centrality_wide_unfold_centrality_(.*)_unfold_trUt0_%s_proj_trY0_fit_fpol3", ref, pt),
             Form("v1_(pionpos)_%s_X_rebin_trUt0_pionpos_rebin_centrality_wide_unfold_centrality_(.*)_unfold_trUt0_%s_proj_trY0_fit_fpol3", ref, pt),
	    });
//    writeMG({Form("v1_(pionneg)_%s_Y_rebin_trUt0_pionneg_rebin_centrality_wide_unfold_centrality_(.*)_unfold_trUt0_%s_proj_trY0_fit_fpol3", ref, pt),
//             Form("v1_(pionpos)_%s_Y_rebin_trUt0_pionpos_rebin_centrality_wide_unfold_centrality_(.*)_unfold_trUt0_%s_proj_trY0_fit_fpol3", ref, pt),
//	    });
    writeMG({Form("v1_(pionneg)_%s_X_rebin_trUt0_pionneg_rebin_centrality_wide_unfold_centrality_unfold_trUt0_%s_proj_trY0_fit_fpol3_par1", ref, pt),
    	     Form("v1_(pionpos)_%s_X_rebin_trUt0_pionpos_rebin_centrality_wide_unfold_centrality_unfold_trUt0_%s_proj_trY0_fit_fpol3_par1", ref, pt),
	    });                                                                  
//    writeMG({Form("v1_(pionneg)_%s_Y_rebin_trUt0_pionneg_rebin_centrality_wide_unfold_centrality_unfold_trUt0_%s_proj_trY0_fit_fpol3_par1", ref, pt),
//    	       Form("v1_(pionpos)_%s_Y_rebin_trUt0_pionpos_rebin_centrality_wide_unfold_centrality_unfold_trUt0_%s_proj_trY0_fit_fpol3_par1", ref, pt),
//	    });
  }
  for (int i=1; i<ptBins_proton.size(); i++)
  {
    string s_pt=Form("%.1f-%.1f", ptBins_proton.at(i-1), ptBins_proton.at(i));
    const char *pt=s_pt.c_str();
    writeMG({Form("v1_proton_%s_X_rebin_trUt0_proton_rebin_centrality_wide_unfold_centrality_(.*)_unfold_trUt0_%s_proj_trY0", ref, pt),
	     "",
	     Form("v1_proton_%s_Y_rebin_trUt0_proton_rebin_centrality_wide_unfold_centrality_(.*)_unfold_trUt0_%s_proj_trY0", ref, pt),
	     });
    writeMG({Form("v1_proton_%s_X_rebin_trUt0_proton_rebin_centrality_wide_unfold_centrality_(.*)_unfold_trUt0_%s_proj_trY0_fit_fpol3", ref, pt)});
//    writeMG({Form("v1_proton_%s_Y_rebin_trUt0_proton_rebin_centrality_wide_unfold_centrality_(.*)_unfold_trUt0_%s_proj_trY0_fit_fpol3", ref, pt)});
    writeMG({Form("v1_proton_%s_X_rebin_trUt0_proton_rebin_centrality_wide_unfold_centrality_unfold_trUt0_%s_proj_trY0_fit_fpol3_par1", ref, pt)});
//    writeMG({Form("v1_proton_%s_Y_rebin_trUt0_proton_rebin_centrality_wide_unfold_centrality_unfold_trUt0_%s_proj_trY0_fit_fpol3_par1", ref, pt)});
  }
  for (int i=1; i<yBins_pion.size(); i++)
  {
    string s_y=Form("%.1f-%.1f", yBins_pion.at(i-1), yBins_pion.at(i));
    const char *y=s_y.c_str();
    writeMG({Form("v1_(pionpos)_%s_X_rebin_trY0_pionpos_rebin_centrality_wide_unfold_centrality_(.*)_unfold_trY0_%s_proj_trUt0", ref, y),
             Form("v1_(pionneg)_%s_X_rebin_trY0_pionneg_rebin_centrality_wide_unfold_centrality_(.*)_unfold_trY0_%s_proj_trUt0", ref, y)
	     });
//    writeMG({Form("v1_(pionpos)_%s_Y_rebin_trY0_pionpos_rebin_centrality_wide_unfold_centrality_(.*)_unfold_trY0_%s_proj_trUt0", ref, y),
//             Form("v1_(pionneg)_%s_Y_rebin_trY0_pionneg_rebin_centrality_wide_unfold_centrality_(.*)_unfold_trY0_%s_proj_trUt0", ref, y)
//	     });
  }
  for (int i=1; i<yBins_proton.size(); i++)
  {
    string s_y=Form("%.1f-%.1f", yBins_proton.at(i-1), yBins_proton.at(i));
    const char *y=s_y.c_str();
    writeMG({Form("v1_proton_%s_X_rebin_trY0_proton_rebin_centrality_wide_unfold_centrality_(.*)_unfold_trY0_%s_proj_trUt0", ref, y)});
//    writeMG({Form("v1_proton_%s_Y_rebin_trY0_proton_rebin_centrality_wide_unfold_centrality_(.*)_unfold_trY0_%s_proj_trUt0", ref, y)});
  }

  outputFile.Close();
}
