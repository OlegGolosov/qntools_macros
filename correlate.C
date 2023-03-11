#include "QnDataFrame.hpp"
#include "utils.h"

vector <vector<string>> Q1Q1=
{
  {"psd0_RECENTERED", "psd1_RECENTERED"},
  {"psd0_RECENTERED", "psd2_RECENTERED"},
  {"psd0_RECENTERED", "psd3_RECENTERED"},
  {"psd1_RECENTERED", "psd2_RECENTERED"},
  {"psd1_RECENTERED", "psd3_RECENTERED"},
  {"psd2_RECENTERED", "psd3_RECENTERED"},
};

vector <vector<string>> u1Q1=
{
  {"proton_pt_RESCALED", "psd0_RECENTERED"},
  {"proton_y_RESCALED", "psd0_RECENTERED"},
  {"proton_pt_RESCALED", "psd1_RECENTERED"},
  {"proton_y_RESCALED", "psd1_RECENTERED"},
  {"proton_pt_RESCALED", "psd2_RECENTERED"},
  {"proton_y_RESCALED", "psd2_RECENTERED"},
  {"proton_pt_RESCALED", "psd3_RECENTERED"},
  {"proton_y_RESCALED", "psd3_RECENTERED"},

  {"pionneg_pt_RESCALED", "psd0_RECENTERED"},
  {"pionneg_y_RESCALED", "psd0_RECENTERED"},
  {"pionneg_pt_RESCALED", "psd1_RECENTERED"},
  {"pionneg_y_RESCALED", "psd1_RECENTERED"},
  {"pionneg_pt_RESCALED", "psd2_RECENTERED"},
  {"pionneg_y_RESCALED", "psd2_RECENTERED"},
  {"pionneg_pt_RESCALED", "psd3_RECENTERED"},
  {"pionneg_y_RESCALED", "psd3_RECENTERED"},

  {"proton_RESCALED", "psd0_RECENTERED"},
  {"proton_RESCALED", "psd1_RECENTERED"},
  {"proton_RESCALED", "psd2_RECENTERED"},
  {"proton_RESCALED", "psd3_RECENTERED"},
  {"pionneg_RESCALED", "psd0_RECENTERED"},
  {"pionneg_RESCALED", "psd1_RECENTERED"},
  {"pionneg_RESCALED", "psd2_RECENTERED"},
  {"pionneg_RESCALED", "psd3_RECENTERED"},
};

void correlate(string inputFiles="qn.root", string outputFile="corr.root")
{
  int nSamples = 50;
  Qn::AxisD centAxis("centrality", {0, 5, 10, 15, 25, 35, 45, 60, 80, 100.});
  auto axes_correlation = Qn::MakeAxes(centAxis);
  TChain *c=makeChain(inputFiles, "tree");
  ROOT::RDataFrame d(*c);
  auto d_samples = Qn::Correlation::Resample(d, nSamples);

  namespace P2 = Qn::Correlation::TwoParticle;
  namespace P3 = Qn::Correlation::MixedHarmonics;
  auto wn = Qn::Correlation::UseWeights::No;
  auto wy = Qn::Correlation::UseWeights::Yes;
  auto wUnity = [](const Qn::QVector &a, const Qn::QVector &b) { return 1; };
  auto wSumWu = [](const Qn::QVector &a, const Qn::QVector &b) { return a.sumweights(); };

  auto corrBuilder =
    Qn::Correlation::CorrelationBuilder{&d_samples, nSamples, axes_correlation};
  for (auto &corr:Q1Q1)
  {
    std::array<std::string, 2> qn{corr.at(0), corr.at(1)};
    string corrName=corr.at(0)+"_"+corr.at(1); 
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_XX_", P2::xx(1, 1), wUnity, wn, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_YY_", P2::yy(1, 1), wUnity, wn, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_XY_", P2::xy(1, 1), wUnity, wn, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_YX_", P2::yx(1, 1), wUnity, wn, qn, qn);
  }
  for (auto &corr:u1Q1)
  {
    std::array<std::string, 2> qn{corr.at(0), corr.at(1)};
    string corrName=corr.at(0)+"_"+corr.at(1); 
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_XX_", P2::xx(1, 1), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_YY_", P2::yy(1, 1), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_XY_", P2::xy(1, 1), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_YX_", P2::yx(1, 1), wSumWu, wy, qn, qn);
  }

  // ---------------- //
  // saving to output //
  // ---------------- //
  auto corrFile = TFile::Open(outputFile.c_str(), "RECREATE");
  corrFile->cd();
  auto results = corrBuilder.GetResults();
  for (auto &res : results) {
    res->Write();
  }
  corrFile->Close();
}
