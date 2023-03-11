#include "utils.h"
#include "makeQvectors.h"
#include "makeDF.h"

void setupQvectors();

void makeQvectors(string inputFiles="tree.root", string calibFilePath="qa.root", string outFilePath="qn.root")
{
  ROOT::RDataFrame d(*makeChain(inputFiles, "t"));
  auto dd=makeDataFrame(d)
         //.Range(1000)
         .Filter("goodEvent");
  varPatterns=
  {
    "centrality",                                             // kEvent
 //   "psdMod(X|Y|Phi|E|Id|Sub.)",        // kChannel 
    "psdMod(Phi|E|Sub0wo45|Sub1wo45|Sub2|Sub3*)",        // kChannel 
    "tr(Pt|Y|Phi|Good)|proton|pionpos|pionneg|trTrackingEffCorr|trPidEffCorr|trTrackingPidEffCorr", // kRecParticle  
    "sim(Pt|Eta|Phi|Pdg|MotherId)"                   // kSimParticle  
  };
  init(dd, outFilePath, calibFilePath);
  setupQvectors(); 
  run(dd);
  cout << "Done!\n";
}

void setupQvectors()
{
  Qn::AxisD axCent("centrality", {0, 5, 10, 15, 25, 35, 45, 60, 80, 100.});
  vector<Qn::AxisD> corrAxesEvent = {axCent};

  Qn::AxisD axPt("trPt", {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.8, 2.2, 2.6, 3.0});
  Qn::AxisD axY("trY", {-0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4});
  vector<Qn::AxisD> corrAxesParticle={axPt, axY};
//  string trWeight="Ones";
//  string trWeight="trTrackingEffCorr";
//  string trWeight="trPidEffCorr";
  string trWeight="trTrackingPidEffCorr";

  for (auto &axis:corrAxesEvent)
    man.AddCorrectionAxis(axis);

  Qn::Recentering recentering;
  recentering.SetApplyWidthEqualization(false);
  Qn::TwistAndRescale twistRescale;
  twistRescale.SetApplyRescale(true);
  twistRescale.SetTwistAndRescaleMethod(Qn::TwistAndRescale::Method::DOUBLE_HARMONIC);
  
  auto sumW=Qn::QVector::Normalization::M;
  auto track=Qn::DetectorType::TRACK;
  auto channel=Qn::DetectorType::CHANNEL;
  auto plain=Qn::QVector::CorrectionStep::PLAIN;
  auto recentered=Qn::QVector::CorrectionStep::RECENTERED;
  auto twisted=Qn::QVector::CorrectionStep::TWIST;
  auto rescaled=Qn::QVector::CorrectionStep::RESCALED;
  
  man.AddDetector("psd0", channel, "psdModPhi", "psdModE", {}, {1}, sumW);
  man.AddCorrectionOnQnVector("psd0", recentering);
  //man.AddCorrectionOnQnVector("psd0", twistRescale);
  man.SetOutputQVectors("psd0", {plain, recentered});
  man.AddCutOnDetector("psd0", {"psdModSub0wo45"}, equal(1), "psd0");
//  man.AddHisto2D("psd0", {{"psdModId", 45, 0., 45}, {"psdModE", 400, 0., 4000}}, "psdModSub0wo45");
//  man.AddHisto2D("psd0", {{"psdModX", 100, -100, 100}, {"psdModY", 100, -100, 100}}, "psdModSub0wo45");
  
  man.AddDetector("psd1", channel, "psdModPhi", "psdModE", {}, {1}, sumW);
  man.AddCorrectionOnQnVector("psd1", recentering);
  //man.AddCorrectionOnQnVector("psd1", twistRescale);
  man.SetOutputQVectors("psd1", {plain, recentered});
  man.AddCutOnDetector("psd1", {"psdModSub1wo45"}, equal(1), "psd1");
//  man.AddHisto2D("psd1", {{"psdModId", 45, 0., 45}, {"psdModE", 400, 0., 4000}}, "psdModSub1wo45");
//  man.AddHisto2D("psd1", {{"psdModX", 100, -100, 100}, {"psdModY", 100, -100, 100}}, "psdModSub1wo45");
  
  man.AddDetector("psd2", channel, "psdModPhi", "psdModE", {}, {1}, sumW);
  man.AddCorrectionOnQnVector("psd2", recentering);
  //man.AddCorrectionOnQnVector("psd2", twistRescale);
  man.SetOutputQVectors("psd2", {plain, recentered});
  man.AddCutOnDetector("psd2", {"psdModSub2"}, equal(1), "psd2");
//  man.AddHisto2D("psd2", {{"psdModId", 45, 0., 45}, {"psdModE", 400, 0., 4000}}, "psdModSub2");
//  man.AddHisto2D("psd2", {{"psdModX", 100, -100, 100}, {"psdModY", 100, -100, 100}}, "psdModSub2");
  
  man.AddDetector("psd3", channel, "psdModPhi", "psdModE", {}, {1}, sumW);
  man.AddCorrectionOnQnVector("psd3", recentering);
  //man.AddCorrectionOnQnVector("psd3", twistRescale);
  man.SetOutputQVectors("psd3", {plain, recentered});
  man.AddCutOnDetector("psd3", {"psdModSub3"}, equal(1), "psd3");
//  man.AddHisto2D("psd3", {{"psdModId", 45, 0., 45}, {"psdModE", 400, 0., 4000}}, "psdModSub3");
//  man.AddHisto2D("psd3", {{"psdModX", 100, -100, 100}, {"psdModY", 100, -100, 100}}, "psdModSub3");

  man.AddDetector("proton", track, "trPhi", trWeight, {axPt, axY}, {1,2}, sumW);
  man.AddCutOnDetector("proton", {"particleType"}, equal(kRecParticle), "recParticle");
  man.AddCutOnDetector("proton", {"proton"}, equal(1), "proton");
  man.AddCorrectionOnQnVector("proton", recentering);
  man.AddCorrectionOnQnVector("proton", twistRescale);
  man.SetOutputQVectors("proton", {plain, recentered, twisted, rescaled});
  man.AddHisto2D("proton", {{"trY", 300, -1., 3.}, {"trPt",  200, 0., 5.}}, "Ones");
  man.AddHisto2D("proton", {{"trPhi", 200, -3.15, 3.15}, {"trPt",  200, 0., 5.}}, "Ones");
  man.AddHisto2D("proton", {{"trPhi", 200, -3.15, 3.15}, {"trY", 300, -1., 3.}}, "Ones");
  
  man.AddDetector("pionneg", track, "trPhi", trWeight, {axPt, axY}, {1,2}, sumW);
  man.AddCutOnDetector("pionneg", {"particleType"}, equal(kRecParticle), "recParticle");
  man.AddCutOnDetector("pionneg", {"pionneg"}, equal(1), "pionneg");
  man.AddCorrectionOnQnVector("pionneg", recentering);
  man.AddCorrectionOnQnVector("pionneg", twistRescale);
  man.SetOutputQVectors("pionneg", {plain, recentered, twisted, rescaled});
  man.AddHisto2D("pionneg", {{"trY", 300, -1., 3.}, {"trPt",  200, 0., 5.}}, "Ones");
  man.AddHisto2D("pionneg", {{"trPhi", 200, -3.15, 3.15}, {"trPt",  200, 0., 5.}}, "Ones");
  man.AddHisto2D("pionneg", {{"trPhi", 200, -3.15, 3.15}, {"trY", 300, -1., 3.}}, "Ones");
  
  man.AddDetector("pionpos", track, "trPhi", trWeight, {axPt, axY}, {1,2}, sumW);
  man.AddCutOnDetector("pionpos", {"particleType"}, equal(kRecParticle), "recParticle");
  man.AddCutOnDetector("pionpos", {"pionpos"}, equal(1), "pionpos");
  man.AddCorrectionOnQnVector("pionpos", recentering);
  man.AddCorrectionOnQnVector("pionpos", twistRescale);
  man.SetOutputQVectors("pionpos", {plain, recentered, twisted, rescaled});
  man.AddHisto2D("pionpos", {{"trY", 300, -1., 3.}, {"trPt",  200, 0., 5.}}, "Ones");
  man.AddHisto2D("pionpos", {{"trPhi", 200, -3.15, 3.15}, {"trPt",  200, 0., 5.}}, "Ones");
  man.AddHisto2D("pionpos", {{"trPhi", 200, -3.15, 3.15}, {"trY", 300, -1., 3.}}, "Ones");
}
