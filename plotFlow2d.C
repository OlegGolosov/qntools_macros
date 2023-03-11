bool verbose=false;

#include "plotFlow.h"

void plotFlow2d(const char *inputFileName="pbpb13_goodEvent_16_011.corr.root", const char *outputFileName="plots2d.root")
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
  rebin("proton_RESCALED_psd._RECENTERED_(XX|YY|XY|YX)_centrality", {"trY", {0.8, 1.2}}, "4S");
  rebin("proton_RESCALED_psd._RECENTERED_(XX|YY|XY|YX)_centrality_rebin_trY_4S", {"trPt", {0., 2.2}}, "4S");
  project("proton_RESCALED_psd._RECENTERED_(XX|YY|XY|YX)_centrality_rebin_trY_4S_rebin_trPt_4S", {"centrality"});
  R1_3S("R1_tpc_3S_X", "proton_RESCALED_psd1_RECENTERED_XX_centrality_rebin_trY_4S_rebin_trPt_4S_proj_centrality", "proton_RESCALED_psd3_RECENTERED_XX_centrality_rebin_trY_4S_rebin_trPt_4S_proj_centrality", "psd1_RECENTERED_psd3_RECENTERED_XX_centrality");
  R1_3S("R1_tpc_3S_Y", "proton_RESCALED_psd1_RECENTERED_YY_centrality_rebin_trY_4S_rebin_trPt_4S_proj_centrality", "proton_RESCALED_psd3_RECENTERED_YY_centrality_rebin_trY_4S_rebin_trPt_4S_proj_centrality", "psd1_RECENTERED_psd3_RECENTERED_YY_centrality");
  R1_3S("R1_psd1_3Sgap_X", "psd1_RECENTERED_psd3_RECENTERED_XX_centrality", "proton_RESCALED_psd1_RECENTERED_XX_centrality_rebin_trY_4S_rebin_trPt_4S_proj_centrality", "proton_RESCALED_psd3_RECENTERED_XX_centrality_rebin_trY_4S_rebin_trPt_4S_proj_centrality"); 
  R1_3S("R1_psd1_3Sgap_Y", "psd1_RECENTERED_psd3_RECENTERED_YY_centrality", "proton_RESCALED_psd1_RECENTERED_YY_centrality_rebin_trY_4S_rebin_trPt_4S_proj_centrality", "proton_RESCALED_psd3_RECENTERED_YY_centrality_rebin_trY_4S_rebin_trPt_4S_proj_centrality"); 
  R1_3S("R1_psd3_3Sgap_X", "psd1_RECENTERED_psd3_RECENTERED_XX_centrality", "proton_RESCALED_psd3_RECENTERED_XX_centrality_rebin_trY_4S_rebin_trPt_4S_proj_centrality", "proton_RESCALED_psd1_RECENTERED_XX_centrality_rebin_trY_4S_rebin_trPt_4S_proj_centrality"); 
  R1_3S("R1_psd3_3Sgap_Y", "psd1_RECENTERED_psd3_RECENTERED_YY_centrality", "proton_RESCALED_psd3_RECENTERED_YY_centrality_rebin_trY_4S_rebin_trPt_4S_proj_centrality", "proton_RESCALED_psd1_RECENTERED_YY_centrality_rebin_trY_4S_rebin_trPt_4S_proj_centrality"); 
  R1_4S("R1_psd1_4S_X", "psd1_RECENTERED_psd3_RECENTERED_XX_centrality", "R1_tpc_3S_X", "proton_RESCALED_psd3_RECENTERED_XX_centrality_rebin_trY_4S_rebin_trPt_4S_proj_centrality"); 
  R1_4S("R1_psd1_4S_Y", "psd1_RECENTERED_psd3_RECENTERED_YY_centrality", "R1_tpc_3S_Y", "proton_RESCALED_psd3_RECENTERED_YY_centrality_rebin_trY_4S_rebin_trPt_4S_proj_centrality"); 
  R1_4S("R1_psd3_4S_X", "psd1_RECENTERED_psd3_RECENTERED_XX_centrality", "R1_tpc_3S_X", "proton_RESCALED_psd1_RECENTERED_XX_centrality_rebin_trY_4S_rebin_trPt_4S_proj_centrality"); 
  R1_4S("R1_psd3_4S_Y", "psd1_RECENTERED_psd3_RECENTERED_YY_centrality", "R1_tpc_3S_Y", "proton_RESCALED_psd1_RECENTERED_YY_centrality_rebin_trY_4S_rebin_trPt_4S_proj_centrality"); 
  R1_4S("R1_psd2_4S_X", "proton_RESCALED_psd2_RECENTERED_XX_centrality_rebin_trY_4S_rebin_trPt_4S_proj_centrality", "R1_tpc_3S_X"); 
  R1_4S("R1_psd2_4S_Y", "proton_RESCALED_psd2_RECENTERED_YY_centrality_rebin_trY_4S_rebin_trPt_4S_proj_centrality", "R1_tpc_3S_Y"); 

  for (auto &part:{"pionneg", "pionpos", "proton"})
    for (auto &ref:{"psd1", "psd2", "psd3"})
      for (auto &comp:{"X", "Y"})
        v1(Form("v1_%s_%s_4S_%s", part, ref, comp), 
           Form("%s_RESCALED_%s_RECENTERED_%s%s_centrality", part, ref, comp, comp), 
           Form("R1_%s_4S_%s", ref, comp)); 

  merge("v1_proton_psd1_4S_(X|Y)");
  merge("v1_proton_psd2_4S_(X|Y)");
  merge("v1_proton_psd3_4S_(X|Y)");
  merge("v1_proton_(psd1|psd2|psd3)_4S_X");
  merge("v1_proton_(psd1|psd2|psd3)_4S_Y");
  merge("v1_proton_(psd1|psd2|psd3)_4S_(X|Y)");

  merge("v1_pionneg_psd1_4S_(X|Y)");
  merge("v1_pionneg_psd2_4S_(X|Y)");
  merge("v1_pionneg_psd3_4S_(X|Y)");
  merge("v1_pionneg_(psd1|psd2|psd3)_4S_X");
  merge("v1_pionneg_(psd1|psd2|psd3)_4S_Y");
  merge("v1_pionneg_(psd1|psd2|psd3)_4S_(X|Y)");

  merge("v1_pionpos_psd1_4S_(X|Y)");
  merge("v1_pionpos_psd2_4S_(X|Y)");
  merge("v1_pionpos_psd3_4S_(X|Y)");
  merge("v1_pionpos_(psd1|psd2|psd3)_4S_X");
  merge("v1_pionpos_(psd1|psd2|psd3)_4S_Y");
  merge("v1_pionpos_(psd1|psd2|psd3)_4S_(X|Y)");

  rebin("v1_.*", {"centrality",{0,10,25,45,60}}, "wide");

  unfold("v1.*", "centrality");
  rebin("v1_(pionneg|pionpos|proton)_(.*)_4S_(X|Y|X_Y)_unfold_centrality_(.._..)", {"trPt",{0,2.2}}, "0_2.2");
  rebin("v1_(proton)_(.*)_4S_(X|Y|X_Y)_rebin_centrality_wide_unfold_centrality_(.._..)", {"trY",{0.0, 0.2, 0.4, 0.6, 1.0, 1.2, 1.4, 1.6, 2.2}}, "proton");
  rebin("v1_(proton)_(.*)_4S_(X|Y|X_Y)_rebin_centrality_wide_unfold_centrality_(.._..)", {"trPt",{0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.4, 1.8}}, "proton");
  rebin("v1_(pionneg|pionpos)_(.*)_4S_(X|Y|X_Y)_rebin_centrality_wide_unfold_centrality_(.._..)", {"trY",{0.0, 0.4, 0.8, 1.2, 2.0}}, "pion");
  rebin("v1_(pionneg|pionpos)_(.*)_4S_(X|Y|X_Y)_rebin_centrality_wide_unfold_centrality_(.._..)", {"trPt",{0.0, 0.2, 0.4, 0.6, 1.0, 1.8}}, "pion");
  unfold("v1_(pionneg|pionpos|proton)_(.*)_4S_(X|Y|X_Y)_rebin_centrality_wide_unfold_centrality_(.._..)_rebin_trPt_(pion|proton)", "trPt", "%1.1f");
  unfold("v1_(pionneg|pionpos|proton)_(.*)_4S_(X|Y|X_Y)_rebin_centrality_wide_unfold_centrality_(.._..)_rebin_trY_(pion|proton)", "trY", "%1.1f");
  project("v1_.*_rebin_trPt_(pion|proton)_unfold_trPt.*", {"trY"});
  project("v1_.*_rebin_trY_(pion|proton)_unfold_trY.*", {"trPt"});
  project("v1_.*_rebin_trPt_0_2.2", {"trY"});
  fit("v1_(pionneg|pionpos|proton)_psd1_psd2_psd3_4S_X_unfold_centrality_.*_rebin_trPt_0_2.2_proj_trY", "pol1", -0.4, 0.4);
  projectFitParameters("v1_pionneg_psd1_psd2_psd3_4S_X_unfold_(centrality)_(.*)_(.*)_rebin_trPt_0_2.2_proj_trY_(fit_pol1)_-0.4_0.4");
  projectFitParameters("v1_pionpos_psd1_psd2_psd3_4S_X_unfold_(centrality)_(.*)_(.*)_rebin_trPt_0_2.2_proj_trY_(fit_pol1)_-0.4_0.4");
  projectFitParameters("v1_proton_psd1_psd2_psd3_4S_X_unfold_(centrality)_(.*)_(.*)_rebin_trPt_0_2.2_proj_trY_(fit_pol1)_-0.4_0.4");

  TFile outputFile(outputFileName, "recreate");
  writeMG({"psd1_RECENTERED_psd2_RECENTERED_(XX|YY|XY|YX)_centrality"});
  writeMG({"psd2_RECENTERED_psd3_RECENTERED_(XX|YY|XY|YX)_centrality"});
  writeMG({"psd1_RECENTERED_psd3_RECENTERED_(XX|YY|XY|YX)_centrality"});
  writeMG({"proton_RESCALED_psd1_RECENTERED_(XX|YY|XY|YX)_centrality_rebin_trY_4S_rebin_trPt_4S_proj_centrality"});
  writeMG({"proton_RESCALED_psd3_RECENTERED_(XX|YY|XY|YX)_centrality_rebin_trY_4S_rebin_trPt_4S_proj_centrality"});
  writeMG({"R1_(psd.)_(3S)_(X)", "R1_(psd.)_(3S)_(Y)"});
  writeMG({"R1_(psd.)_(4S)_(X)", "R1_(psd.)_(4S)_(Y)"});
  writeMG({"R1_(psd1|psd3)_(4S)_(X)", "R1_(psd1|psd3)_(4S)_(Y)", "R1_(psd1|psd3)_(3Sgap)_(X)", "R1_(psd1|psd3)_(3Sgap)_(Y)"});

//  writeMG({"v1_proton_psd1_psd2_psd3_4S_X_rebin_centrality_wide_unfold_centrality_(.._..)_rebin_trY_0_1.2_proj_trPt"});
//  writeMG({"v1_proton_psd1_psd2_psd3_4S_X_rebin_centrality_wide_unfold_centrality_(.._..)_rebin_trPt_0_2.2_proj_trY"});
//  writeMG({"v1_proton_psd1_psd2_psd3_4S_(X|Y)_rebin_centrality_wide_unfold_centrality_10_25_rebin_trY_0_1.2_proj_trPt"});
//  writeMG({"v1_proton_(psd.)_4S_X_Y_rebin_centrality_wide_unfold_centrality_10_25_rebin_trY_0_1.2_proj_trPt"});
  writeMG({"v1_proton_psd1_psd2_psd3_4S_X_Y_rebin_centrality_wide_unfold_centrality_10_25_rebin_trPt_proton_unfold_trPt_(.*)_proj_trY"});
  writeMG({"v1_proton_psd1_psd2_psd3_4S_X_Y_rebin_centrality_wide_unfold_centrality_10_25_rebin_trY_proton_unfold_trY_(.*)_proj_trPt"});
  writeMG({"v1_proton_psd1_psd2_psd3_4S_X_unfold_centrality_(.._..)_rebin_trPt_0_2.2_proj_trY_fit_pol1_-0.4_0.4"});

//  writeMG({"v1_pionneg_psd1_psd2_psd3_4S_X_rebin_centrality_wide_unfold_centrality_(.._..)_rebin_trY_0_1.2_proj_trPt"});
//  writeMG({"v1_pionneg_psd1_psd2_psd3_4S_X_rebin_centrality_wide_unfold_centrality_(.._..)_rebin_trPt_0_2.2_proj_trY"});
//  writeMG({"v1_pionneg_psd1_psd2_psd3_4S_(X|Y)_rebin_centrality_wide_unfold_centrality_10_25_rebin_trY_0_1.2_proj_trPt"});
//  writeMG({"v1_pionneg_(psd.)_4S_X_Y_rebin_centrality_wide_unfold_centrality_10_25_rebin_trY_0_1.2_proj_trPt"});
  writeMG({"v1_pionneg_psd1_psd2_psd3_4S_X_Y_rebin_centrality_wide_unfold_centrality_10_25_rebin_trPt_pion_unfold_trPt_(.*)_proj_trY"});
  writeMG({"v1_pionneg_psd1_psd2_psd3_4S_X_Y_rebin_centrality_wide_unfold_centrality_10_25_rebin_trY_pion_unfold_trY_(.*)_proj_trPt"});
  writeMG({"v1_pionneg_psd1_psd2_psd3_4S_X_unfold_centrality_(.._..)_rebin_trPt_0_2.2_proj_trY_fit_pol1_-0.4_0.4"});

//  writeMG({"v1_pionpos_psd1_psd2_psd3_4S_X_rebin_centrality_wide_unfold_centrality_(.._..)_rebin_trY_0_1.2_proj_trPt"});
//  writeMG({"v1_pionpos_psd1_psd2_psd3_4S_X_rebin_centrality_wide_unfold_centrality_(.._..)_rebin_trPt_0_2.2_proj_trY"});
//  writeMG({"v1_pionpos_psd1_psd2_psd3_4S_(X|Y)_rebin_centrality_wide_unfold_centrality_10_25_rebin_trY_0_1.2_proj_trPt"});
//  writeMG({"v1_pionpos_(psd.)_4S_X_Y_rebin_centrality_wide_unfold_centrality_10_25_rebin_trY_0_1.2_proj_trPt"});
  writeMG({"v1_pionpos_psd1_psd2_psd3_4S_X_Y_rebin_centrality_wide_unfold_centrality_10_25_rebin_trPt_pion_unfold_trPt_(.*)_proj_trY"});
  writeMG({"v1_pionpos_psd1_psd2_psd3_4S_X_Y_rebin_centrality_wide_unfold_centrality_10_25_rebin_trY_pion_unfold_trY_(.*)_proj_trPt"});
  writeMG({"v1_pionpos_psd1_psd2_psd3_4S_X_unfold_centrality_(.._..)_rebin_trPt_0_2.2_proj_trY_fit_pol1_-0.4_0.4"});

  writeMG({"v1_(pionneg|pionpos|proton)_psd1_psd2_psd3_4S_X_unfold_centrality_rebin_trPt_0_2.2_proj_trY_fit_pol1_-0.4_0.4_par1"});

//  for (auto &part:{"proton", "pionneg", "pionpos"})
//    for (auto &ref1:{"psd1", "psd2"})
//      for (auto &ref2:{"psd2", "psd3"})
//        for (auto &ucomp:{"X", "Y"})
//          for (auto &q1comp:{"X", "Y"})
//            for (auto &q2comp:{"X", "Y"})
//              v2(Form("v2_%s_%s%s_4S_%s2%s%s", part, ref1, ref2, ucomp, q1comp, q2comp), 
//                 Form("%s_RESCALED_%s_RECENTERED_%s_RECENTERED_%s2%s%s_centrality", part, ref1, ref2, ucomp, q1comp, q2comp), 
//                 Form("R1_%s_4S_%s", ref1, q1comp),
//                 Form("R1_%s_4S_%s", ref2, q2comp));
//
//  scale("v2_.*X2YY", -1);
//
//  for (auto &part:{"pionneg", "pionpos", "proton"})
//  {
//    merge(Form("v2_%s_psd1psd2_4S_(Y2XY|Y2YX)", part));
//    merge(Form("v2_%s_psd2psd3_4S_(Y2XY|Y2YX)", part));
//    merge(Form("v2_%s_psd1psd3_4S_(Y2XY|Y2YX)", part));
//    merge(Form("v2_%s_psd1psd2_4S_(X2XX|X2YY)", part));
//    merge(Form("v2_%s_psd2psd3_4S_(X2XX|X2YY)", part));
//    merge(Form("v2_%s_psd1psd3_4S_(X2XX|X2YY)", part));
//    merge(Form("v2_%s_psd1psd2_4S_(X2XX|X2YY|Y2XY|Y2YX)", part));
//    merge(Form("v2_%s_psd2psd3_4S_(X2XX|X2YY|Y2XY|Y2YX)", part));
//    merge(Form("v2_%s_psd1psd3_4S_(X2XX|X2YY|Y2XY|Y2YX)", part));
//    merge(Form("v2_%s_(psd1psd2|psd2psd3|psd1psd3)_4S_X2XX", part));
//    merge(Form("v2_%s_(psd1psd2|psd2psd3|psd1psd3)_4S_X2YY", part));
//    merge(Form("v2_%s_(psd1psd2|psd2psd3|psd1psd3)_4S_Y2XY", part));
//    merge(Form("v2_%s_(psd1psd2|psd2psd3|psd1psd3)_4S_Y2YX", part));
//    merge(Form("v2_%s_(psd1psd2|psd2psd3|psd1psd3)_4S_(Y2XY|Y2YX)", part));
//    merge(Form("v2_%s_(psd1psd2|psd2psd3|psd1psd3)_4S_(X2XX|X2YY)", part));
//    merge(Form("v2_%s_(psd1psd2|psd2psd3|psd1psd3)_4S_(X2XX|X2YY|Y2XY|Y2YX)", part));
//  }
// 
//  rebin("v2_.*", {"centrality",{0,10,25,45,60}}, "wide");
//  rebin("v2_(pionneg|pionpos|proton)_psd1psd2_psd2psd3_psd1psd3_4S_Y2XY_Y2YX", {"centrality",{15,35}}, "wide2");
//
//  unfold("v2.*", "centrality");
//  rebin("v2_.*", {"trPt",{0,2.2}}, "0_2.2");
//  rebin("v2_.*", {"trY",{0,1.2}}, "0_1.2");
//  project("v2_.*_rebin_trPt_0_2.2", {"trY"});
//  project("v2_.*_rebin_trY_0_1.2", {"trPt"});


//  writeMG({"v2_proton_(psd1psd2|psd2psd3|psd1psd3)_4S_X2XX_X2YY_Y2XY_Y2YX_rebin_centrality_wide2_unfold_centrality_15_35_rebin_trY_0_1.2_proj_trPt"});
//  writeMG({"v2_proton_(psd1psd2|psd2psd3|psd1psd3)_4S_X2XX_X2YY_rebin_centrality_wide2_unfold_centrality_15_35_rebin_trY_0_1.2_proj_trPt"});
//  writeMG({"v2_proton_(psd1psd2|psd2psd3|psd1psd3)_4S_Y2XY_Y2YX_rebin_centrality_wide2_unfold_centrality_15_35_rebin_trY_0_1.2_proj_trPt"});
//  writeMG({"v2_proton_psd1psd2_psd2psd3_psd1psd3_4S_(X2XX|X2YY|Y2XY|Y2YX)_rebin_centrality_wide2_unfold_centrality_15_35_rebin_trY_0_1.2_proj_trPt"});
//  writeMG({"v2_proton_psd1psd2_psd2psd3_psd1psd3_4S_Y2XY_Y2YX_rebin_centrality_wide_unfold_centrality_(.._..)_rebin_trY_0_1.2_proj_trPt"});
//
//  writeMG({"v2_pionneg_(psd1psd2|psd2psd3|psd1psd3)_4S_X2XX_X2YY_Y2XY_Y2YX_rebin_centrality_wide2_unfold_centrality_15_35_rebin_trY_0_1.2_proj_trPt"});
//  writeMG({"v2_pionneg_(psd1psd2|psd2psd3|psd1psd3)_4S_X2XX_X2YY_rebin_centrality_wide2_unfold_centrality_15_35_rebin_trY_0_1.2_proj_trPt"});
//  writeMG({"v2_pionneg_(psd1psd2|psd2psd3|psd1psd3)_4S_Y2XY_Y2YX_rebin_centrality_wide2_unfold_centrality_15_35_rebin_trY_0_1.2_proj_trPt"});
//  writeMG({"v2_pionneg_psd1psd2_psd2psd3_psd1psd3_4S_(X2XX|X2YY|Y2XY|Y2YX)_rebin_centrality_wide2_unfold_centrality_15_35_rebin_trY_0_1.2_proj_trPt"});
//  writeMG({"v2_pionneg_psd1psd2_psd2psd3_psd1psd3_4S_Y2XY_Y2YX_rebin_centrality_wide_unfold_centrality_(.._..)_rebin_trY_0_1.2_proj_trPt"});
//
//  writeMG({"v2_pionpos_(psd1psd2|psd2psd3|psd1psd3)_4S_X2XX_X2YY_Y2XY_Y2YX_rebin_centrality_wide2_unfold_centrality_15_35_rebin_trY_0_1.2_proj_trPt"});
//  writeMG({"v2_pionpos_(psd1psd2|psd2psd3|psd1psd3)_4S_X2XX_X2YY_rebin_centrality_wide2_unfold_centrality_15_35_rebin_trY_0_1.2_proj_trPt"});
//  writeMG({"v2_pionpos_(psd1psd2|psd2psd3|psd1psd3)_4S_Y2XY_Y2YX_rebin_centrality_wide2_unfold_centrality_15_35_rebin_trY_0_1.2_proj_trPt"});
//  writeMG({"v2_pionpos_psd1psd2_psd2psd3_psd1psd3_4S_(X2XX|X2YY|Y2XY|Y2YX)_rebin_centrality_wide2_unfold_centrality_15_35_rebin_trY_0_1.2_proj_trPt"});
//  writeMG({"v2_pionpos_psd1psd2_psd2psd3_psd1psd3_4S_Y2XY_Y2YX_rebin_centrality_wide_unfold_centrality_(.._..)_rebin_trY_0_1.2_proj_trPt"});
//
//  writeMG({"v2_(proton|pionneg|pionpos)_psd1psd2_psd2psd3_psd1psd3_4S_Y2XY_Y2YX_rebin_centrality_wide2_unfold_centrality_15_35_rebin_trY_0_1.2_proj_trPt"});
  outputFile.Close();
}
