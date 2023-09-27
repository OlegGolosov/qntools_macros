bool verbose=false;

#include "plotFlow.h"

void plotFlow(const char *inputFileName="corr.root", const char *outputFileName="plots.root", int pBeam=13, const char *ref="psd1_psd2_psd3", const char *comp="X")
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
    rebin("proton_RESCALED_psd._RECENTERED_(XX|YY|XY|YX)_centrality", {"trY", {0.8, 1.2}}, "4S");
    rebin("proton_RESCALED_psd._RECENTERED_(XX|YY|XY|YX)_centrality_rebin_trY_4S", {"trPt", {0.4, 2.2}}, "4S");
    project("proton_RESCALED_psd._RECENTERED_(XX|YY|XY|YX)_centrality_rebin_trY_4S_rebin_trPt_4S", {"centrality"});
    res="4S";
    refPart4S="proton";
    resScale4S=1;
  }
  if (pBeam==30 || pBeam==41)
  {
    rebin("pionpos_RESCALED_psd._RECENTERED_(XX|YY|XY|YX)_centrality", {"trY", {0.4, 1.0}}, "4S");
    rebin("pionpos_RESCALED_psd._RECENTERED_(XX|YY|XY|YX)_centrality_rebin_trY_4S", {"trPt", {0.4, 1.0}}, "4S");
    project("pionpos_RESCALED_psd._RECENTERED_(XX|YY|XY|YX)_centrality_rebin_trY_4S_rebin_trPt_4S", {"centrality"});
    res="3S";
    refPart4S="pionpos";
    resScale4S=-1;
  }
  R1_3S("R1_tpc_3S_X", refPart4S+"_RESCALED_psd1_RECENTERED_XX_centrality_rebin_trY_4S_rebin_trPt_4S_proj_centrality", refPart4S+"_RESCALED_psd3_RECENTERED_XX_centrality_rebin_trY_4S_rebin_trPt_4S_proj_centrality", "psd1_RECENTERED_psd3_RECENTERED_XX_centrality", centAxOptimal);
  R1_3S("R1_tpc_3S_Y", refPart4S+"_RESCALED_psd1_RECENTERED_YY_centrality_rebin_trY_4S_rebin_trPt_4S_proj_centrality", refPart4S+"_RESCALED_psd3_RECENTERED_YY_centrality_rebin_trY_4S_rebin_trPt_4S_proj_centrality", "psd1_RECENTERED_psd3_RECENTERED_YY_centrality", centAxOptimal);
  R1_3S("R1_psd1_4S_X", "psd1_RECENTERED_psd3_RECENTERED_XX_centrality", refPart4S+"_RESCALED_psd1_RECENTERED_XX_centrality_rebin_trY_4S_rebin_trPt_4S_proj_centrality", refPart4S+"_RESCALED_psd3_RECENTERED_XX_centrality_rebin_trY_4S_rebin_trPt_4S_proj_centrality", centAxOptimal); 
  R1_3S("R1_psd1_4S_Y", "psd1_RECENTERED_psd3_RECENTERED_YY_centrality", refPart4S+"_RESCALED_psd1_RECENTERED_YY_centrality_rebin_trY_4S_rebin_trPt_4S_proj_centrality", refPart4S+"_RESCALED_psd3_RECENTERED_YY_centrality_rebin_trY_4S_rebin_trPt_4S_proj_centrality", centAxOptimal); 
  R1_3S("R1_psd3_4S_X", "psd1_RECENTERED_psd3_RECENTERED_XX_centrality", refPart4S+"_RESCALED_psd3_RECENTERED_XX_centrality_rebin_trY_4S_rebin_trPt_4S_proj_centrality", refPart4S+"_RESCALED_psd1_RECENTERED_XX_centrality_rebin_trY_4S_rebin_trPt_4S_proj_centrality", centAxOptimal); 
  R1_3S("R1_psd3_4S_Y", "psd1_RECENTERED_psd3_RECENTERED_YY_centrality", refPart4S+"_RESCALED_psd3_RECENTERED_YY_centrality_rebin_trY_4S_rebin_trPt_4S_proj_centrality", refPart4S+"_RESCALED_psd1_RECENTERED_YY_centrality_rebin_trY_4S_rebin_trPt_4S_proj_centrality", centAxOptimal); 
  R1_4S("R1_psd2_4S_X", refPart4S+"_RESCALED_psd2_RECENTERED_XX_centrality_rebin_trY_4S_rebin_trPt_4S_proj_centrality", "R1_tpc_3S_X", resScale4S, centAxOptimal); 
  R1_4S("R1_psd2_4S_Y", refPart4S+"_RESCALED_psd2_RECENTERED_YY_centrality_rebin_trY_4S_rebin_trPt_4S_proj_centrality", "R1_tpc_3S_Y", resScale4S, centAxOptimal); 

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

  vector<double> centBinsWide={0,10,25,45,60};
  vector<double> centBinsStar={10,35};
  rebin("v1_.*", {"centrality", centBinsWide}, "wide");
  rebin("v1_.*_(X|Y)", {"centrality", centBinsStar}, "star");
  unfold("v1.*", "centrality");

  Qn::AxisD axPt_proton, axY_proton, axPt_pion, axY_pion;
  vector<double> ptBins_dv1dy;
  if(pBeam==13)
  {
    //axY_proton={"trY",{0.0, 0.2, 0.4, 0.6, 1.0, 1.2, 1.4, 1.6, 2.2}};
    //axPt_proton={"trPt",{0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.4, 1.8}};
    //axY_pion={"trY",{0.0, 0.4, 0.8, 1.2, 2.0}};
    ptBins_dv1dy={0., 0.4, 2.2};
    axY_proton={"trY",{0.0, 0.6, 1.0, 1.2, 1.4, 1.6, 2.2}};
    axPt_proton={"trPt",{0.0, 0.6, 1.0, 1.4, 1.8}};
    axY_pion={"trY",{0.0, 0.6, 1.0, 2.0}};
    axPt_pion={"trPt",{0.0, 0.2, 0.4, 0.6, 1.0, 1.8}};
  }
  if(pBeam==30 || pBeam==41)
  {
    ptBins_dv1dy={0., 0.4, 2.2};
    axY_proton={"trY",{0.0, 0.6, 1.0, 1.2, 1.4, 1.6, 2.2}};
    axPt_proton={"trPt",{0.0, 0.6, 1.0, 1.4, 1.8}};
    axY_pion={"trY",{0.0, 0.6, 1.0, 2.0}};
    axPt_pion={"trPt",{0.0, 0.2, 0.4, 0.6, 1.0, 1.8}};
  }

  rebin("v1_proton_(.*)_(X|Y|X_Y)_rebin_centrality_star_unfold_centrality_(..-..)", {"trPt",{0.4,2.2}}, "star");
  rebin("v1_(pionneg|pionpos)_(.*)_(X|Y|X_Y)_rebin_centrality_star_unfold_centrality_(..-..)", {"trPt",{0.2,2.2}}, "star");
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
//  changeBinEdges("v1_.*_proj_trY", "trY", [yBeam](double y) {return y/yBeam;}, "trY0");
//  changeBinEdges("v1_pion.*_proj_trPt", "trPt", [ut0coefPion](double pt) {return pt*ut0coefPion;}, "trUt0");
//  changeBinEdges("v1_proton.*_proj_trPt", "trPt", [ut0coefProton](double pt) {return pt*ut0coefProton;}, "trUt0");
  project("v1_.*_rebin_trPt_star", {"trY"});

  auto fpol3=new TF1("fpol3","[0]+[1]*x+[2]*x*x*x", 0, 1.5);
  auto fpol3star=new TF1("fpol3star","[1]*x+[2]*x*x*x", 0, 1.5);
  fit(Form("v1_(pionneg|pionpos|proton)_%s_%s_rebin_centrality_star_unfold_centrality_.*_rebin_trPt_star_proj_trY", ref, comp), "fpol3star");
  projectFitParameters(Form("v1_pionneg_%s_%s_rebin_centrality_star_unfold_(centrality)_(.*)-(.*)_rebin_trPt_star_proj_trY_(fit_fpol3star)", ref, comp));
  projectFitParameters(Form("v1_pionpos_%s_%s_rebin_centrality_star_unfold_(centrality)_(.*)-(.*)_rebin_trPt_star_proj_trY_(fit_fpol3star)", ref, comp));
  projectFitParameters(Form("v1_proton_%s_%s_rebin_centrality_star_unfold_(centrality)_(.*)-(.*)_rebin_trPt_star_proj_trY_(fit_fpol3star)", ref, comp));
  for (int i=1; i<ptBins_dv1dy.size();i++)
  {
    const char *pt=Form("%.1f-%.1f", ptBins_dv1dy.at(i-1), ptBins_dv1dy.at(i));
    fit(Form("v1_(pionneg|pionpos|proton)_%s_%s_unfold_centrality_.*_rebin_trPt_dv1dy_unfold_trPt_%s_proj_trY", ref, comp, pt), "fpol3");
    projectFitParameters(Form("v1_pionneg_%s_%s_unfold_(centrality)_(.*)-(.*)_rebin_trPt_dv1dy_unfold_trPt_%s_proj_trY_(fit_fpol3)", ref, comp, pt));
    projectFitParameters(Form("v1_pionpos_%s_%s_unfold_(centrality)_(.*)-(.*)_rebin_trPt_dv1dy_unfold_trPt_%s_proj_trY_(fit_fpol3)", ref, comp, pt));
    projectFitParameters(Form("v1_proton_%s_%s_unfold_(centrality)_(.*)-(.*)_rebin_trPt_dv1dy_unfold_trPt_%s_proj_trY_(fit_fpol3)", ref, comp, pt));
  }

  TFile outputFile(outputFileName, "recreate");
  cout << "\nPlotting...\n";
  writeMG({"psd1_RECENTERED_psd2_RECENTERED_(XX|YY|XY|YX)_centrality"});
  writeMG({"psd2_RECENTERED_psd3_RECENTERED_(XX|YY|XY|YX)_centrality"});
  writeMG({"psd1_RECENTERED_psd3_RECENTERED_(XX|YY|XY|YX)_centrality"});
  writeMG({refPart4S+"_RESCALED_psd1_RECENTERED_(XX|YY|XY|YX)_centrality_rebin_trY_4S_rebin_trPt_4S_proj_centrality"});
  writeMG({refPart4S+"_RESCALED_psd3_RECENTERED_(XX|YY|XY|YX)_centrality_rebin_trY_4S_rebin_trPt_4S_proj_centrality"});
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

    writeMG({Form("v1_%s_%s_%s_rebin_centrality_star_unfold_centrality_(..-..)_rebin_trPt_star_proj_trY_fit_fpol3star", part, ref, comp)});
    writeMG({Form("v1_%s_%s_%s_unfold_centrality_rebin_trPt_dv1dy_unfold_trPt_(.*)_proj_trY_fit_fpol3_par1", part, ref, comp)});
  }

  for (int i=1; i<centBinsWide.size(); i++)
  {
    const char *cent=Form("%02.0f-%02.0f", centBinsWide.at(i-1), centBinsWide.at(i));
    writeMG({Form("v1_(pionneg)_%s_%s_rebin_centrality_wide_unfold_centrality_%s_rebin_trPt_pionneg_unfold_trPt_(.*)_proj_trY", ref, comp, cent),
             Form("v1_(pionpos)_%s_%s_rebin_centrality_wide_unfold_centrality_%s_rebin_trPt_pionpos_unfold_trPt_(.*)_proj_trY", ref, comp, cent)});
    writeMG({Form("v1_(pionneg)_%s_%s_rebin_centrality_wide_unfold_centrality_%s_rebin_trY_pionneg_unfold_trY_(.*)_proj_trPt", ref, comp, cent),
             Form("v1_(pionpos)_%s_%s_rebin_centrality_wide_unfold_centrality_%s_rebin_trY_pionpos_unfold_trY_(.*)_proj_trPt", ref, comp, cent)});
  }

  for (int i=1; i<ptBins_dv1dy.size();i++)
  {
    const char *pt=Form("%.1f-%.1f", ptBins_dv1dy.at(i-1), ptBins_dv1dy.at(i));
    writeMG({Form("v1_(pionneg|pionpos|proton)_%s_%s_unfold_centrality_rebin_trPt_dv1dy_unfold_trPt_%s_proj_trY_fit_fpol3_par1", ref, comp, pt)});
  }

  writeMG({Form("v1_(pionneg|pionpos|proton)_%s_%s_rebin_centrality_star_unfold_centrality_rebin_trPt_star_proj_trY_fit_fpol3star_par1", ref, comp)});


  outputFile.Close();
}
