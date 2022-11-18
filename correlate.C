#include "QnDataFrame.hpp"
#include "utils.h"

vector <vector<string>> Q1Q1=
{
  {"fhcal1_RECENTERED", "fhcal2_RECENTERED"},
  {"fhcal2_RECENTERED", "fhcal3_RECENTERED"},
  {"fhcal3_RECENTERED", "fhcal1_RECENTERED"},
};

vector <vector<string>> u1Q1=
{
  {"tr_RESCALED", "fhcal1_RECENTERED"},
  {"tr_RESCALED", "fhcal2_RECENTERED"},
  {"tr_RESCALED", "fhcal3_RECENTERED"},
};

void correlate(string inputFiles="qn.root", string outputFile="corr.root")
{
  int nSamples = 50;
  Qn::AxisD centAxis({"b", 4,0,10});
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
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_XX", P2::xx(1, 1), wUnity, wn, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_YY", P2::yy(1, 1), wUnity, wn, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_XY", P2::xy(1, 1), wUnity, wn, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_YX", P2::yx(1, 1), wUnity, wn, qn, qn);
  }
  for (auto &corr:u1Q1)
  {
    std::array<std::string, 2> qn{corr.at(0), corr.at(1)};
    string corrName=corr.at(0)+"_"+corr.at(1); 
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_XX", P2::xx(1, 1), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_YY", P2::yy(1, 1), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_XY", P2::xy(1, 1), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_YX", P2::yx(1, 1), wSumWu, wy, qn, qn);
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
