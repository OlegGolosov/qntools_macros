#include "utils.h"
#include "makeQvectors.h"
#include "makeDF.h"

void setupQvectors();

void makeQvectors(string inputFiles="tree.root", string calibFilePath="qa.root", string outFilePath="qn.root")
{
  ROOT::RDataFrame d(*makeChain(inputFiles, "t"));
  auto dd=makeDataFrame(d)
         //.Range(1000)
         .Filter("goodEventOld");
  varPatterns=
  {
    "centrality",                                             // kEvent
 //   "psdMod(X|Y|Phi|E|Id|Sub.)",        // kChannel 
    "psdMod(Phi|E|Sub0wo45|Sub1wo45|Sub2|Sub3*)",        // kChannel 
    "tr(Pt|Y|Phi)|proton|pionneg", // kRecParticle  
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

  Qn::AxisD axPt("trPt", {0.0f, 0.2f, 0.4f, 0.6f, 0.8f, 1.0f, 1.2f, 1.4f, 1.8f, 2.2f, 2.6f, 3.0f});
  Qn::AxisD axY("trY", 15, -0.6, 2.4);
  vector<Qn::AxisD> corrAxesParticle={axPt, axY};

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

  man.AddDetector("proton_pt", track, "trPhi", "Ones", {axPt}, {1,2}, sumW);
  man.AddCutOnDetector("proton_pt", {"particleType"}, equal(kRecParticle), "recParticle");
  man.AddCutOnDetector("proton_pt", {"proton"}, equal(1), "proton");
  man.AddCutOnDetector("proton_pt", {"trY"}, range(0, 1.2), "yCut");
  man.AddCorrectionOnQnVector("proton_pt", recentering);
  man.AddCorrectionOnQnVector("proton_pt", twistRescale);
  man.SetOutputQVectors("proton_pt", {plain, recentered, twisted, rescaled});
  man.AddHisto2D("proton_pt", {{"trY", 300, -1., 3.}, {"trPt",  200, 0., 5.}}, "Ones");
  man.AddHisto2D("proton_pt", {{"trPhi", 200, -3.15, 3.15}, {"trPt",  200, 0., 5.}}, "Ones");
  man.AddHisto2D("proton_pt", {{"trPhi", 200, -3.15, 3.15}, {"trY", 300, -1., 3.}}, "Ones");
  
  man.AddDetector("proton_y", track, "trPhi", "Ones", {axY}, {1,2}, sumW);
  man.AddCutOnDetector("proton_y", {"particleType"}, equal(kRecParticle), "recParticle");
  man.AddCutOnDetector("proton_y", {"proton"}, equal(1), "proton");
  man.AddCutOnDetector("proton_y", {"trPt"}, range(0, 1.2), "ptCut");
  man.AddCorrectionOnQnVector("proton_y", recentering);
  man.AddCorrectionOnQnVector("proton_y", twistRescale);
  man.SetOutputQVectors("proton_y", {plain, recentered, twisted, rescaled});
  man.AddHisto2D("proton_y", {{"trY", 300, -1., 3.}, {"trPt",  200, 0., 5.}}, "Ones");
  man.AddHisto2D("proton_y", {{"trPhi", 200, -3.15, 3.15}, {"trPt",  200, 0., 5.}}, "Ones");
  man.AddHisto2D("proton_y", {{"trPhi", 200, -3.15, 3.15}, {"trY", 300, -1., 3.}}, "Ones");
  
  man.AddDetector("pionneg_pt", track, "trPhi", "Ones", {axPt}, {1,2}, sumW);
  man.AddCutOnDetector("pionneg_pt", {"particleType"}, equal(kRecParticle), "recParticle");
  man.AddCutOnDetector("pionneg_pt", {"pionneg"}, equal(1), "pionneg");
  man.AddCutOnDetector("pionneg_pt", {"trY"}, range(0, 1.2), "yCut");
  man.AddCorrectionOnQnVector("pionneg_pt", recentering);
  man.AddCorrectionOnQnVector("pionneg_pt", twistRescale);
  man.SetOutputQVectors("pionneg_pt", {plain, recentered, twisted, rescaled});
  man.AddHisto2D("pionneg_pt", {{"trY", 300, -1., 3.}, {"trPt",  200, 0., 5.}}, "Ones");
  man.AddHisto2D("pionneg_pt", {{"trPhi", 200, -3.15, 3.15}, {"trPt",  200, 0., 5.}}, "Ones");
  man.AddHisto2D("pionneg_pt", {{"trPhi", 200, -3.15, 3.15}, {"trY", 300, -1., 3.}}, "Ones");
  
  man.AddDetector("pionneg_y", track, "trPhi", "Ones", {axY}, {1,2}, sumW);
  man.AddCutOnDetector("pionneg_y", {"particleType"}, equal(kRecParticle), "recParticle");
  man.AddCutOnDetector("pionneg_y", {"pionneg"}, equal(1), "pionneg");
  man.AddCutOnDetector("pionneg_y", {"trPt"}, range(0, 1.2), "ptCut");
  man.AddCorrectionOnQnVector("pionneg_y", recentering);
  man.AddCorrectionOnQnVector("pionneg_y", twistRescale);
  man.SetOutputQVectors("pionneg_y", {plain, recentered, twisted, rescaled});
  man.AddHisto2D("pionneg_y", {{"trY", 300, -1., 3.}, {"trPt",  200, 0., 5.}}, "Ones");
  man.AddHisto2D("pionneg_y", {{"trPhi", 200, -3.15, 3.15}, {"trPt",  200, 0., 5.}}, "Ones");
  man.AddHisto2D("pionneg_y", {{"trPhi", 200, -3.15, 3.15}, {"trY", 300, -1., 3.}}, "Ones");
}
