bool verbose=false;

#include "plotFlow.h"

void plotFlow(const char *inputFileName="pbpb13_16_011_goodEvent_16_011_mergedCorr.corr.root", const char *outputFileName="plots_mergedCorr.root")
{
  cout << "START!\n";
  TFile inputFile(inputFileName);
  for (const auto &k:*inputFile.GetListOfKeys()) 
    corrMap.emplace(k->GetName(), *(Qn::DataContainerStatCollect*)inputFile.Get(k->GetName()));  
  R1_3S("R1_psd1_3S_X", "psd1_RECENTERED_psd3_RECENTERED_XX_centrality", "psd1_RECENTERED_psd2_RECENTERED_XX_centrality", "psd2_RECENTERED_psd3_RECENTERED_XX_centrality"); 
  R1_3S("R1_psd1_3S_Y", "psd1_RECENTERED_psd3_RECENTERED_YY_centrality", "psd1_RECENTERED_psd2_RECENTERED_YY_centrality", "psd2_RECENTERED_psd3_RECENTERED_YY_centrality"); 
  R1_3S("R1_psd2_3S_X", "psd1_RECENTERED_psd2_RECENTERED_XX_centrality", "psd2_RECENTERED_psd3_RECENTERED_XX_centrality", "psd1_RECENTERED_psd3_RECENTERED_XX_centrality"); 
  R1_3S("R1_psd2_3S_Y", "psd1_RECENTERED_psd2_RECENTERED_YY_centrality", "psd2_RECENTERED_psd3_RECENTERED_YY_centrality", "psd1_RECENTERED_psd3_RECENTERED_YY_centrality"); 
  R1_3S("R1_psd3_3S_X", "psd2_RECENTERED_psd3_RECENTERED_XX_centrality", "psd1_RECENTERED_psd3_RECENTERED_XX_centrality", "psd1_RECENTERED_psd2_RECENTERED_XX_centrality"); 
  R1_3S("R1_psd3_3S_Y", "psd2_RECENTERED_psd3_RECENTERED_YY_centrality", "psd1_RECENTERED_psd3_RECENTERED_YY_centrality", "psd1_RECENTERED_psd2_RECENTERED_YY_centrality"); 
  rebin("proton_y_RESCALED_psd._RECENTERED_(XX|YY|XY|YX)_centrality", {"trY", {0.8, 1.2}}, "4S");
  project("proton_y_RESCALED_psd._RECENTERED_(XX|YY|XY|YX)_centrality_rebin_trY_4S", {"centrality"});
  R1_3S("R1_tpc_3S_X", "proton_y_RESCALED_psd1_RECENTERED_XX_centrality_rebin_trY_4S_proj_centrality", "proton_y_RESCALED_psd3_RECENTERED_XX_centrality_rebin_trY_4S_proj_centrality", "psd1_RECENTERED_psd3_RECENTERED_XX_centrality");
  R1_3S("R1_tpc_3S_Y", "proton_y_RESCALED_psd1_RECENTERED_YY_centrality_rebin_trY_4S_proj_centrality", "proton_y_RESCALED_psd3_RECENTERED_YY_centrality_rebin_trY_4S_proj_centrality", "psd1_RECENTERED_psd3_RECENTERED_YY_centrality");
  R1_3S("R1_psd1_3Sgap_X", "psd1_RECENTERED_psd3_RECENTERED_XX_centrality", "proton_y_RESCALED_psd1_RECENTERED_XX_centrality_rebin_trY_4S_proj_centrality", "proton_y_RESCALED_psd3_RECENTERED_XX_centrality_rebin_trY_4S_proj_centrality"); 
  R1_3S("R1_psd1_3Sgap_Y", "psd1_RECENTERED_psd3_RECENTERED_YY_centrality", "proton_y_RESCALED_psd1_RECENTERED_YY_centrality_rebin_trY_4S_proj_centrality", "proton_y_RESCALED_psd3_RECENTERED_YY_centrality_rebin_trY_4S_proj_centrality"); 
  R1_3S("R1_psd3_3Sgap_X", "psd1_RECENTERED_psd3_RECENTERED_XX_centrality", "proton_y_RESCALED_psd3_RECENTERED_XX_centrality_rebin_trY_4S_proj_centrality", "proton_y_RESCALED_psd1_RECENTERED_XX_centrality_rebin_trY_4S_proj_centrality"); 
  R1_3S("R1_psd3_3Sgap_Y", "psd1_RECENTERED_psd3_RECENTERED_YY_centrality", "proton_y_RESCALED_psd3_RECENTERED_YY_centrality_rebin_trY_4S_proj_centrality", "proton_y_RESCALED_psd1_RECENTERED_YY_centrality_rebin_trY_4S_proj_centrality"); 
  R1_4S("R1_psd1_4S_X", "psd1_RECENTERED_psd3_RECENTERED_XX_centrality", "R1_tpc_3S_X", "proton_y_RESCALED_psd3_RECENTERED_XX_centrality_rebin_trY_4S_proj_centrality"); 
  R1_4S("R1_psd1_4S_Y", "psd1_RECENTERED_psd3_RECENTERED_YY_centrality", "R1_tpc_3S_Y", "proton_y_RESCALED_psd3_RECENTERED_YY_centrality_rebin_trY_4S_proj_centrality"); 
  R1_4S("R1_psd3_4S_X", "psd1_RECENTERED_psd3_RECENTERED_XX_centrality", "R1_tpc_3S_X", "proton_y_RESCALED_psd1_RECENTERED_XX_centrality_rebin_trY_4S_proj_centrality"); 
  R1_4S("R1_psd3_4S_Y", "psd1_RECENTERED_psd3_RECENTERED_YY_centrality", "R1_tpc_3S_Y", "proton_y_RESCALED_psd1_RECENTERED_YY_centrality_rebin_trY_4S_proj_centrality"); 
  R1_4S("R1_psd2_4S_X", "proton_y_RESCALED_psd2_RECENTERED_XX_centrality_rebin_trY_4S_proj_centrality", "R1_tpc_3S_X"); 
  R1_4S("R1_psd2_4S_Y", "proton_y_RESCALED_psd2_RECENTERED_YY_centrality_rebin_trY_4S_proj_centrality", "R1_tpc_3S_Y"); 

  for (auto &part:{"proton", "pionneg", "pionpos"})
    for (auto &ax:{"y", "pt"})
      for (auto &ref:{"psd1", "psd2", "psd3"})
        for (auto &comp:{"X", "Y"})
          v1(Form("v1_%s_%s_%s_4S_%s", part, ax, ref, comp), 
             Form("%s_%s_RESCALED_%s_RECENTERED_%s%s_centrality", part, ax, ref, comp, comp), 
             Form("R1_%s_4S_%s", ref, comp)); 

  merge("v1_proton_y_psd1_4S_(X|Y)");
  merge("v1_proton_y_psd2_4S_(X|Y)");
  merge("v1_proton_y_psd3_4S_(X|Y)");
  merge("v1_proton_y_(psd1|psd2|psd3)_4S_X");
  merge("v1_proton_y_(psd1|psd2|psd3)_4S_Y");
  merge("v1_proton_y_(psd1|psd2|psd3)_4S_(X|Y)");

  merge("v1_pionneg_y_psd1_4S_(X|Y)");
  merge("v1_pionneg_y_psd2_4S_(X|Y)");
  merge("v1_pionneg_y_psd3_4S_(X|Y)");
  merge("v1_pionneg_y_(psd1|psd2|psd3)_4S_X");
  merge("v1_pionneg_y_(psd1|psd2|psd3)_4S_Y");
  merge("v1_pionneg_y_(psd1|psd2|psd3)_4S_(X|Y)");

  merge("v1_pionpos_y_psd1_4S_(X|Y)");
  merge("v1_pionpos_y_psd2_4S_(X|Y)");
  merge("v1_pionpos_y_psd3_4S_(X|Y)");
  merge("v1_pionpos_y_(psd1|psd2|psd3)_4S_X");
  merge("v1_pionpos_y_(psd1|psd2|psd3)_4S_Y");
  merge("v1_pionpos_y_(psd1|psd2|psd3)_4S_(X|Y)");

  merge("v1_proton_pt_psd1_4S_(X|Y)");
  merge("v1_proton_pt_psd2_4S_(X|Y)");
  merge("v1_proton_pt_psd3_4S_(X|Y)");
  merge("v1_proton_pt_(psd1|psd2|psd3)_4S_X");
  merge("v1_proton_pt_(psd1|psd2|psd3)_4S_Y");
  merge("v1_proton_pt_(psd1|psd2|psd3)_4S_(X|Y)");

  merge("v1_pionneg_pt_psd1_4S_(X|Y)");
  merge("v1_pionneg_pt_psd2_4S_(X|Y)");
  merge("v1_pionneg_pt_psd3_4S_(X|Y)");
  merge("v1_pionneg_pt_(psd1|psd2|psd3)_4S_X");
  merge("v1_pionneg_pt_(psd1|psd2|psd3)_4S_Y");
  merge("v1_pionneg_pt_(psd1|psd2|psd3)_4S_(X|Y)");

  merge("v1_pionpos_pt_psd1_4S_(X|Y)");
  merge("v1_pionpos_pt_psd2_4S_(X|Y)");
  merge("v1_pionpos_pt_psd3_4S_(X|Y)");
  merge("v1_pionpos_pt_(psd1|psd2|psd3)_4S_X");
  merge("v1_pionpos_pt_(psd1|psd2|psd3)_4S_Y");
  merge("v1_pionpos_pt_(psd1|psd2|psd3)_4S_(X|Y)");
  
  rebin("v1_.*", {"centrality",{0,10,25,45,60}}, "wide");
  rebin("v1_.*_pt_.*(X|Y)", {"centrality",{15,35}}, "wide2");
 
  unfold("v1.*", "centrality");
  project("v1_.*_y_.*", {"trY"});
  project("v1_.*_pt_.*", {"trPt"});
  fit("v1_.*_unfold_centrality_.*_proj_trY", "pol1", -0.4, 0.4);

  projectFitParameters("v1_proton_y_psd1_psd2_psd3_4S_X_unfold_(centrality)_(.*)_(.*)_proj_trY_(fit_pol1)_-0.4_0.4");
  projectFitParameters("v1_proton_y_psd1_psd2_psd3_4S_Y_unfold_(centrality)_(.*)_(.*)_proj_trY_(fit_pol1)_-0.4_0.4");
  projectFitParameters("v1_proton_y_psd1_psd2_psd3_4S_XY_unfold_(centrality)_(.*)_(.*)_proj_trY_(fit_pol1)_-0.4_0.4");

  projectFitParameters("v1_pionneg_y_psd1_psd2_psd3_4S_X_unfold_(centrality)_(.*)_(.*)_proj_trY_(fit_pol1)_-0.4_0.4");
  projectFitParameters("v1_pionneg_y_psd1_psd2_psd3_4S_Y_unfold_(centrality)_(.*)_(.*)_proj_trY_(fit_pol1)_-0.4_0.4");
  projectFitParameters("v1_pionneg_y_psd1_psd2_psd3_4S_XY_unfold_(centrality)_(.*)_(.*)_proj_trY_(fit_pol1)_-0.4_0.4");

  projectFitParameters("v1_pionpos_y_psd1_psd2_psd3_4S_X_unfold_(centrality)_(.*)_(.*)_proj_trY_(fit_pol1)_-0.4_0.4");
  projectFitParameters("v1_pionpos_y_psd1_psd2_psd3_4S_Y_unfold_(centrality)_(.*)_(.*)_proj_trY_(fit_pol1)_-0.4_0.4");
  projectFitParameters("v1_pionpos_y_psd1_psd2_psd3_4S_XY_unfold_(centrality)_(.*)_(.*)_proj_trY_(fit_pol1)_-0.4_0.4");

  for (auto &part:{"proton", "pionneg", "pionpos"})
    for (auto &ax:{"y", "pt"})
      for (auto &ref1:{"psd1", "psd2"})
        for (auto &ref2:{"psd2", "psd3"})
          for (auto &ucomp:{"X", "Y"})
            for (auto &q1comp:{"X", "Y"})
              for (auto &q2comp:{"X", "Y"})
                v2(Form("v2_%s_%s_%s%s_4S_%s2%s%s", part, ax, ref1, ref2, ucomp, q1comp, q2comp), 
                   Form("%s_%s_RESCALED_%s_RECENTERED_%s_RECENTERED_%s2%s%s_centrality", part, ax, ref1, ref2, ucomp, q1comp, q2comp), 
                   Form("R1_%s_4S_%s", ref1, q1comp),
                   Form("R1_%s_4S_%s", ref2, q2comp));

  scale("v2_.*X2YY", -1);

  for (auto &part:{"pionneg", "pionpos", "proton"})
  {
    merge(Form("v2_%s_pt_psd1psd2_4S_(Y2XY|Y2YX)", part));
    merge(Form("v2_%s_pt_psd2psd3_4S_(Y2XY|Y2YX)", part));
    merge(Form("v2_%s_pt_psd1psd3_4S_(Y2XY|Y2YX)", part));
    merge(Form("v2_%s_pt_psd1psd2_4S_(X2XX|X2YY)", part));
    merge(Form("v2_%s_pt_psd2psd3_4S_(X2XX|X2YY)", part));
    merge(Form("v2_%s_pt_psd1psd3_4S_(X2XX|X2YY)", part));
    merge(Form("v2_%s_pt_psd1psd2_4S_(X2XX|X2YY|Y2XY|Y2YX)", part));
    merge(Form("v2_%s_pt_psd2psd3_4S_(X2XX|X2YY|Y2XY|Y2YX)", part));
    merge(Form("v2_%s_pt_psd1psd3_4S_(X2XX|X2YY|Y2XY|Y2YX)", part));
    merge(Form("v2_%s_pt_(psd1psd2|psd2psd3|psd1psd3)_4S_X2XX", part));
    merge(Form("v2_%s_pt_(psd1psd2|psd2psd3|psd1psd3)_4S_X2YY", part));
    merge(Form("v2_%s_pt_(psd1psd2|psd2psd3|psd1psd3)_4S_Y2XY", part));
    merge(Form("v2_%s_pt_(psd1psd2|psd2psd3|psd1psd3)_4S_Y2YX", part));
    merge(Form("v2_%s_pt_(psd1psd2|psd2psd3|psd1psd3)_4S_(Y2XY|Y2YX)", part));
    merge(Form("v2_%s_pt_(psd1psd2|psd2psd3|psd1psd3)_4S_(X2XX|X2YY)", part));
    merge(Form("v2_%s_pt_(psd1psd2|psd2psd3|psd1psd3)_4S_(X2XX|X2YY|Y2XY|Y2YX)", part));
  }
 
  rebin("v2_.*", {"centrality",{0,10,25,45,60}}, "wide");
  rebin("v2_.*_pt_.*(X|Y)", {"centrality",{15,35}}, "wide2");

  unfold("v2.*", "centrality");
  project("v2_.*_y_.*", {"trY"});
  project("v2_.*_pt_.*", {"trPt"});

  TFile outputFile(outputFileName, "recreate");

  writeMG({"psd1_RECENTERED_psd2_RECENTERED_(XX|YY|XY|YX)_centrality"});
  writeMG({"psd2_RECENTERED_psd3_RECENTERED_(XX|YY|XY|YX)_centrality"});
  writeMG({"psd1_RECENTERED_psd3_RECENTERED_(XX|YY|XY|YX)_centrality"});
  writeMG({"proton_y_RESCALED_psd1_RECENTERED_(XX|YY|XY|YX)_centrality_rebin_trY_4S_proj_centrality"});
  writeMG({"proton_y_RESCALED_psd3_RECENTERED_(XX|YY|XY|YX)_centrality_rebin_trY_4S_proj_centrality"});
  writeMG({"R1_(psd.)_(3S)_(X)", "R1_(psd.)_(3S)_(Y)"});
  writeMG({"R1_(psd.)_(4S)_(X)", "R1_(psd.)_(4S)_(Y)"});
  writeMG({"R1_(psd1|psd3)_(4S)_(X)", "R1_(psd1|psd3)_(4S)_(Y)", "R1_(psd1|psd3)_(3Sgap)_(X)", "R1_(psd1|psd3)_(3Sgap)_(Y)"});

  writeMG({"v1_proton_pt_psd1_psd2_psd3_4S_X_rebin_centrality_wide_unfold_centrality_(.*)_proj_trPt"});
  writeMG({"v1_proton_y_psd1_psd2_psd3_4S_X_rebin_centrality_wide_unfold_centrality_(.*)_proj_trY"});
  writeMG({"v1_proton_pt_psd1_psd2_psd3_4S_(X|Y)_rebin_centrality_wide_unfold_centrality_10_25_proj_trPt"});
  writeMG({"v1_proton_pt_(psd.)_4S_X_Y_rebin_centrality_wide_unfold_centrality_10_25_proj_trPt"});
  writeMG({"v1_proton_y_psd1_psd2_psd3_4S_X_unfold_centrality_(.*)_proj_trY_fit_pol1_-0.4_0.4"});

  writeMG({"v1_pionneg_pt_psd1_psd2_psd3_4S_X_rebin_centrality_wide_unfold_centrality_(.*)_proj_trPt"});
  writeMG({"v1_pionneg_y_psd1_psd2_psd3_4S_X_rebin_centrality_wide_unfold_centrality_(.*)_proj_trY"});
  writeMG({"v1_pionneg_pt_psd1_psd2_psd3_4S_(X|Y)_rebin_centrality_wide_unfold_centrality_10_25_proj_trPt"});
  writeMG({"v1_pionneg_pt_(psd.)_4S_X_Y_rebin_centrality_wide_unfold_centrality_10_25_proj_trPt"});
  writeMG({"v1_pionneg_y_psd1_psd2_psd3_4S_X_unfold_centrality_(.*)_proj_trY_fit_pol1_-0.4_0.4"});
  
  writeMG({"v1_pionpos_pt_psd1_psd2_psd3_4S_X_rebin_centrality_wide_unfold_centrality_(.*)_proj_trPt"});
  writeMG({"v1_pionpos_y_psd1_psd2_psd3_4S_X_rebin_centrality_wide_unfold_centrality_(.*)_proj_trY"});
  writeMG({"v1_pionpos_pt_psd1_psd2_psd3_4S_(X|Y)_rebin_centrality_wide_unfold_centrality_10_25_proj_trPt"});
  writeMG({"v1_pionpos_pt_(psd.)_4S_X_Y_rebin_centrality_wide_unfold_centrality_10_25_proj_trPt"});
  writeMG({"v1_pionpos_y_psd1_psd2_psd3_4S_X_unfold_centrality_(.*)_proj_trY_fit_pol1_-0.4_0.4"});
  
  writeMG({"v1_(proton|pionneg|pionpos)_pt_psd1_psd2_psd3_4S_X_rebin_centrality_wide2_unfold_centrality_15_35_proj_trPt"});
  writeMG({"v1_(proton|pionneg|pionpos)_y_psd1_psd2_psd3_4S_X_unfold_centrality_proj_trY_fit_pol1_-0.4_0.4_par1"});

  writeMG({"v2_proton_pt_(psd1psd2|psd2psd3|psd1psd3)_4S_X2XX_X2YY_Y2XY_Y2YX_rebin_centrality_wide2_unfold_centrality_15_35_proj_trPt"});
  writeMG({"v2_proton_pt_(psd1psd2|psd2psd3|psd1psd3)_4S_X2XX_X2YY_rebin_centrality_wide2_unfold_centrality_15_35_proj_trPt"});
  writeMG({"v2_proton_pt_(psd1psd2|psd2psd3|psd1psd3)_4S_Y2XY_Y2YX_rebin_centrality_wide2_unfold_centrality_15_35_proj_trPt"});
  writeMG({"v2_proton_pt_psd1psd2_psd2psd3_psd1psd3_4S_(X2XX|X2YY|Y2XY|Y2YX)_rebin_centrality_wide2_unfold_centrality_15_35_proj_trPt"});
  writeMG({"v2_proton_pt_psd1psd2_psd2psd3_psd1psd3_4S_Y2XY_Y2YX_rebin_centrality_wide_unfold_centrality_(.*)_proj_trPt"});

  writeMG({"v2_pionneg_pt_(psd1psd2|psd2psd3|psd1psd3)_4S_X2XX_X2YY_Y2XY_Y2YX_rebin_centrality_wide2_unfold_centrality_15_35_proj_trPt"});
  writeMG({"v2_pionneg_pt_(psd1psd2|psd2psd3|psd1psd3)_4S_X2XX_X2YY_rebin_centrality_wide2_unfold_centrality_15_35_proj_trPt"});
  writeMG({"v2_pionneg_pt_(psd1psd2|psd2psd3|psd1psd3)_4S_Y2XY_Y2YX_rebin_centrality_wide2_unfold_centrality_15_35_proj_trPt"});
  writeMG({"v2_pionneg_pt_psd1psd2_psd2psd3_psd1psd3_4S_(X2XX|X2YY|Y2XY|Y2YX)_rebin_centrality_wide2_unfold_centrality_15_35_proj_trPt"});
  writeMG({"v2_pionneg_pt_psd1psd2_psd2psd3_psd1psd3_4S_Y2XY_Y2YX_rebin_centrality_wide_unfold_centrality_(.*)_proj_trPt"});

  writeMG({"v2_pionpos_pt_(psd1psd2|psd2psd3|psd1psd3)_4S_X2XX_X2YY_Y2XY_Y2YX_rebin_centrality_wide2_unfold_centrality_15_35_proj_trPt"});
  writeMG({"v2_pionpos_pt_(psd1psd2|psd2psd3|psd1psd3)_4S_X2XX_X2YY_rebin_centrality_wide2_unfold_centrality_15_35_proj_trPt"});
  writeMG({"v2_pionpos_pt_(psd1psd2|psd2psd3|psd1psd3)_4S_Y2XY_Y2YX_rebin_centrality_wide2_unfold_centrality_15_35_proj_trPt"});
  writeMG({"v2_pionpos_pt_psd1psd2_psd2psd3_psd1psd3_4S_(X2XX|X2YY|Y2XY|Y2YX)_rebin_centrality_wide2_unfold_centrality_15_35_proj_trPt"});
  writeMG({"v2_pionpos_pt_psd1psd2_psd2psd3_psd1psd3_4S_Y2XY_Y2YX_rebin_centrality_wide_unfold_centrality_(.*)_proj_trPt"});

  writeMG({"v2_(proton|pionneg|pionpos)_pt_psd1psd2_psd2psd3_psd1psd3_4S_Y2XY_Y2YX_rebin_centrality_wide2_unfold_centrality_15_35_proj_trPt"});

  outputFile.Close();
}
