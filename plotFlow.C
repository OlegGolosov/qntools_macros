bool verbose=false;

#include "plotFlow.h"

void plotFlow(const char *inputFileName="pbpb13_goodEvent_16_011.corr.root", const char *outputFileName="plots2d.root", int pBeam=13, const char *ref="psd1_psd2_psd3", const char *comp="X")
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
  string refPart;
  float scale;
  string res="3S";
  if (pBeam==13)
  {
    rebin("proton_RESCALED_psd._RECENTERED_(XX|YY|XY|YX)_centrality", {"trY", {0.8, 1.2}}, "4S");
    rebin("proton_RESCALED_psd._RECENTERED_(XX|YY|XY|YX)_centrality_rebin_trY_4S", {"trPt", {0.4, 2.2}}, "4S");
    project("proton_RESCALED_psd._RECENTERED_(XX|YY|XY|YX)_centrality_rebin_trY_4S_rebin_trPt_4S", {"centrality"});
    refPart="proton";
    scale=1;
    res="4S";
  }
  if (pBeam==30 || pBeam==41)
  {
    rebin("pionpos_RESCALED_psd._RECENTERED_(XX|YY|XY|YX)_centrality", {"trY", {0.4, 1.0}}, "4S");
    rebin("pionpos_RESCALED_psd._RECENTERED_(XX|YY|XY|YX)_centrality_rebin_trY_4S", {"trPt", {0.4, 1.0}}, "4S");
    project("pionpos_RESCALED_psd._RECENTERED_(XX|YY|XY|YX)_centrality_rebin_trY_4S_rebin_trPt_4S", {"centrality"});
    refPart="pionpos";
    scale=-1;
  }
  R1_3S("R1_tpc_3S_X", refPart+"_RESCALED_psd1_RECENTERED_XX_centrality_rebin_trY_4S_rebin_trPt_4S_proj_centrality", refPart+"_RESCALED_psd3_RECENTERED_XX_centrality_rebin_trY_4S_rebin_trPt_4S_proj_centrality", "psd1_RECENTERED_psd3_RECENTERED_XX_centrality");
  R1_3S("R1_tpc_3S_Y", refPart+"_RESCALED_psd1_RECENTERED_YY_centrality_rebin_trY_4S_rebin_trPt_4S_proj_centrality", refPart+"_RESCALED_psd3_RECENTERED_YY_centrality_rebin_trY_4S_rebin_trPt_4S_proj_centrality", "psd1_RECENTERED_psd3_RECENTERED_YY_centrality");
  R1_3S("R1_psd1_4S_X", "psd1_RECENTERED_psd3_RECENTERED_XX_centrality", refPart+"_RESCALED_psd1_RECENTERED_XX_centrality_rebin_trY_4S_rebin_trPt_4S_proj_centrality", refPart+"_RESCALED_psd3_RECENTERED_XX_centrality_rebin_trY_4S_rebin_trPt_4S_proj_centrality"); 
  R1_3S("R1_psd1_4S_Y", "psd1_RECENTERED_psd3_RECENTERED_YY_centrality", refPart+"_RESCALED_psd1_RECENTERED_YY_centrality_rebin_trY_4S_rebin_trPt_4S_proj_centrality", refPart+"_RESCALED_psd3_RECENTERED_YY_centrality_rebin_trY_4S_rebin_trPt_4S_proj_centrality"); 
  R1_3S("R1_psd3_4S_X", "psd1_RECENTERED_psd3_RECENTERED_XX_centrality", refPart+"_RESCALED_psd3_RECENTERED_XX_centrality_rebin_trY_4S_rebin_trPt_4S_proj_centrality", refPart+"_RESCALED_psd1_RECENTERED_XX_centrality_rebin_trY_4S_rebin_trPt_4S_proj_centrality"); 
  R1_3S("R1_psd3_4S_Y", "psd1_RECENTERED_psd3_RECENTERED_YY_centrality", refPart+"_RESCALED_psd3_RECENTERED_YY_centrality_rebin_trY_4S_rebin_trPt_4S_proj_centrality", refPart+"_RESCALED_psd1_RECENTERED_YY_centrality_rebin_trY_4S_rebin_trPt_4S_proj_centrality"); 
  R1_4S("R1_psd2_4S_X", refPart+"_RESCALED_psd2_RECENTERED_XX_centrality_rebin_trY_4S_rebin_trPt_4S_proj_centrality", "R1_tpc_3S_X", scale); 
  R1_4S("R1_psd2_4S_Y", refPart+"_RESCALED_psd2_RECENTERED_YY_centrality_rebin_trY_4S_rebin_trPt_4S_proj_centrality", "R1_tpc_3S_Y", scale); 

  for (auto &part_:{"pionneg", "pionpos", "proton"})
    for (auto &ref_:{"psd1", "psd2", "psd3"})
      for (auto &comp_:{"X", "Y"})
        v1(Form("v1_%s_%s_%s", part_, ref_, comp_), 
           Form("%s_RESCALED_%s_RECENTERED_%s%s_centrality", part_, ref_, comp_, comp_), 
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

  vector<double> centBinsWide={0,10,25,45,60};
  rebin("v1_.*", {"centrality", centBinsWide}, "wide");
  unfold("v1.*", "centrality");

  Qn::AxisD axPt_proton, axY_proton, axPt_pion, axY_pion;
  vector<double> ptBins_dv1dy={0., 0.6, 1.0, 1.8};
  if(pBeam==13)
  {
    axY_proton={"trY",{0.0, 0.2, 0.4, 0.6, 1.0, 1.2, 1.4, 1.6, 2.2}};
    axPt_proton={"trPt",{0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.4, 1.8}};
    axY_pion={"trY",{0.0, 0.4, 0.8, 1.2, 2.0}};
    axPt_pion={"trPt",{0.0, 0.2, 0.4, 0.6, 1.0, 1.8}};
  }
  if(pBeam==30)
  {
    axY_proton={"trY",{0.0, 0.2, 0.4, 0.6, 1.0, 1.2, 1.4, 1.6, 2.2}};
    axPt_proton={"trPt",{0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.4, 1.8}};
    axY_pion={"trY",{0.0, 0.4, 0.8, 1.2, 2.0}};
    axPt_pion={"trPt",{0.0, 0.2, 0.4, 0.6, 1.0, 1.8}};
  }
  if(pBeam==41)
  {
    ptBins_dv1dy={0., 0.4, 2.2};
    axY_proton={"trY",{0.0, 0.6, 1.0, 1.2, 1.4, 1.6, 2.2}};
    axPt_proton={"trPt",{0.0, 0.6, 0.8, 1.0, 1.4, 1.8}};
    axY_pion={"trY",{0.0, 0.4, 0.8, 1.2, 2.0}};
    axPt_pion={"trPt",{0.0, 0.2, 0.4, 1.0, 1.8}};
  }

  rebin("v1_(pionneg|pionpos|proton)_(.*)_(X|Y|X_Y)_unfold_centrality_(..-..)", {"trPt",{0.4,2.2}}, "star");
  rebin("v1_(pionneg|pionpos|proton)_(.*)_(X|Y|X_Y)_unfold_centrality_(..-..)", {"trPt", ptBins_dv1dy}, "dv1dy");
  rebin("v1_proton_(.*)_(X|Y|X_Y)_rebin_centrality_wide_unfold_centrality_(..-..)", axY_proton,  "proton");
  rebin("v1_proton_(.*)_(X|Y|X_Y)_rebin_centrality_wide_unfold_centrality_(..-..)", axPt_proton, "proton");
  rebin("v1_pionpos_(.*)_(X|Y|X_Y)_rebin_centrality_wide_unfold_centrality_(..-..)", axY_pion, "pionpos");
  rebin("v1_pionpos_(.*)_(X|Y|X_Y)_rebin_centrality_wide_unfold_centrality_(..-..)", axPt_pion, "pionpos");
  rebin("v1_pionneg_(.*)_(X|Y|X_Y)_rebin_centrality_wide_unfold_centrality_(..-..)", axY_pion, "pionneg");
  rebin("v1_pionneg_(.*)_(X|Y|X_Y)_rebin_centrality_wide_unfold_centrality_(..-..)", axPt_pion, "pionneg");
  unfold("v1_(pionneg|pionpos|proton)_(.*)_(X|Y|X_Y)_rebin_centrality_wide_unfold_centrality_(..-..)_rebin_trPt_(pionneg|pionpos|proton)", "trPt", "%1.1f");
  unfold("v1_(pionneg|pionpos|proton)_(.*)_(X|Y|X_Y)_rebin_centrality_wide_unfold_centrality_(..-..)_rebin_trY_(pionneg|pionpos|proton)", "trY", "%1.1f");
  unfold("v1_(pionneg|pionpos|proton)_(.*)_(X|Y|X_Y)_unfold_centrality_(..-..)_rebin_trPt_dv1dy", "trPt", "%1.1f");
  project("v1_.*_rebin_trPt_(pionneg|pionpos|proton|dv1dy)_unfold_trPt.*", {"trY"});
  project("v1_.*_rebin_trY_(pionneg|pionpos|proton)_unfold_trY.*", {"trPt"});
  project("v1_.*_rebin_trPt_star", {"trY"});
  auto fpol3=new TF1("fpol3","[0]+[1]*x+[2]*x*x*x", -0.5, 1.5);
  fit(Form("v1_(pionneg|pionpos|proton)_%s_%s_unfold_centrality_.*_rebin_trPt_star_proj_trY", ref, comp), "fpol3");
  projectFitParameters(Form("v1_pionneg_%s_%s_unfold_(centrality)_(.*)-(.*)_rebin_trPt_star_proj_trY_(fit_fpol3)", ref, comp));
  projectFitParameters(Form("v1_pionpos_%s_%s_unfold_(centrality)_(.*)-(.*)_rebin_trPt_star_proj_trY_(fit_fpol3)", ref, comp));
  projectFitParameters(Form("v1_proton_%s_%s_unfold_(centrality)_(.*)-(.*)_rebin_trPt_star_proj_trY_(fit_fpol3)", ref, comp));
  for (int i=1; i<ptBins_dv1dy.size();i++)
  {
    const char *pt=Form("%.1f-%.1f", ptBins_dv1dy.at(i-1), ptBins_dv1dy.at(i));
    fit(Form("v1_(pionneg|pionpos|proton)_%s_%s_unfold_centrality_.*_rebin_trPt_dv1dy_unfold_trPt_%s_proj_trY", ref, comp, pt), "fpol3");
    projectFitParameters(Form("v1_pionneg_%s_%s_unfold_(centrality)_(.*)-(.*)_rebin_trPt_dv1dy_unfold_trPt_%s_proj_trY_(fit_fpol3)", ref, comp, pt));
    projectFitParameters(Form("v1_pionpos_%s_%s_unfold_(centrality)_(.*)-(.*)_rebin_trPt_dv1dy_unfold_trPt_%s_proj_trY_(fit_fpol3)", ref, comp, pt));
    projectFitParameters(Form("v1_proton_%s_%s_unfold_(centrality)_(.*)-(.*)_rebin_trPt_dv1dy_unfold_trPt_%s_proj_trY_(fit_fpol3)", ref, comp, pt));
  }

  TFile outputFile(outputFileName, "recreate");
  writeMG({"psd1_RECENTERED_psd2_RECENTERED_(XX|YY|XY|YX)_centrality"});
  writeMG({"psd2_RECENTERED_psd3_RECENTERED_(XX|YY|XY|YX)_centrality"});
  writeMG({"psd1_RECENTERED_psd3_RECENTERED_(XX|YY|XY|YX)_centrality"});
  writeMG({refPart+"_RESCALED_psd1_RECENTERED_(XX|YY|XY|YX)_centrality_rebin_trY_4S_rebin_trPt_4S_proj_centrality"});
  writeMG({refPart+"_RESCALED_psd3_RECENTERED_(XX|YY|XY|YX)_centrality_rebin_trY_4S_rebin_trPt_4S_proj_centrality"});
  writeMG({"R1_(psd.)_(3S)_(X)", "R1_(psd.)_(3S)_(Y)"});
  writeMG({"R1_(psd.)_(4S)_(X)", "R1_(psd.)_(4S)_(Y)"});
  writeMG({"R1_tpc_3S_(X)", "R1_tpc_3S_(Y)"});
  writeMG({"R1_(psd.)_(4S)_(X)", "R1_(psd.)_(4S)_(Y)", "R1_(psd.)_(3S)_(X)", "R1_(psd.)_(3S)_(Y)"});

  for (auto &part:{"pionneg", "pionpos", "proton"})
  {
    for (int i=1; i<centBinsWide.size(); i++)
    {
      const char *cent=Form("%02.0f-%02.0f", centBinsWide.at(i-1), centBinsWide.at(i));
      writeMG({Form("v1_%s_%s_%s_rebin_centrality_wide_unfold_centrality_%s_rebin_trPt_%s_unfold_trPt_(.*)_proj_trY", part, ref, comp, cent, part)});
      writeMG({Form("v1_%s_%s_%s_rebin_centrality_wide_unfold_centrality_%s_rebin_trY_%s_unfold_trY_(.*)_proj_trPt", part, ref, comp, cent, part)});
    }

    for (int i=1; i<ptBins_dv1dy.size();i++)
    {
      const char *pt=Form("%.1f-%.1f", ptBins_dv1dy.at(i-1), ptBins_dv1dy.at(i));
      writeMG({Form("v1_%s_%s_%s_unfold_centrality_(..-..)_rebin_trPt_dv1dy_unfold_trPt_%s_proj_trY_fit_fpol3", part, ref, comp, pt)});
    }

    writeMG({Form("v1_%s_%s_%s_unfold_centrality_(..-..)_rebin_trPt_star_proj_trY_fit_fpol3", part, ref, comp)});
    writeMG({Form("v1_%s_%s_%s_unfold_centrality_rebin_trPt_dv1dy_unfold_trPt_(.*)_proj_trY_fit_fpol3_par1", part, ref, comp)});
  }
  writeMG({Form("v1_(pionneg|pionpos|proton)_%s_%s_unfold_centrality_rebin_trPt_star_proj_trY_fit_fpol3_par1", ref, comp)});


  outputFile.Close();
}
