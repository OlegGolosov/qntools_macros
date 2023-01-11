#include "utils.h"
#include "makeQvectors.h"

filteredDF defineVariables(definedDF &d);
void setupQvectors();

void makeQvectors(string inputFiles="/home/ogolosov/desktop/bman/data/run8/sim/dcm_4gev.root", string calibFilePath="qa.root", string outFilePath="qn.root")
{
  ROOT::RDataFrame d(*makeChain(inputFiles, "t"));
  auto dd=defineVariables(d);
  init(dd, outFilePath, calibFilePath);
  setupQvectors(); 
  run(dd);
  cout << "Done!\n";
}

filteredDF defineVariables(definedDF &d)
{
  auto dd=d
    .Define("fhcalModPhi","RVec<float> phi; for(auto& pos:fhcalModPos) phi.push_back(pos.phi()); return phi;")
    .Define("fhcalModX","RVec<float> x; for(auto& pos:fhcalModPos) x.push_back(pos.x()); return x;")
    .Define("fhcalModY","RVec<float> y; for(auto& pos:fhcalModPos) y.push_back(pos.y()); return y;")
    .Define("fhcalModInSub1","fhcalModId>=0&&fhcalModId<=16")
    .Define("fhcalModInSub2","fhcalModId>=17&&fhcalModId<=34")
    .Define("fhcalModInSub3","fhcalModId>=35&&fhcalModId<=54")
    .Define("scwallModPhi","RVec<float> phi; for(auto& pos:scwallModPos) phi.push_back(pos.phi()); return phi;")
    .Define("trPt","RVec<float> pt; for(auto& mom:trMom) pt.push_back(mom.pt()); return pt;")
    .Define("trEta","RVec<float> eta; for(auto& mom:trMom) eta.push_back(mom.eta()); return eta;")
    .Define("trPhi","RVec<float> phi;for(auto& mom:trMom) phi.push_back(mom.phi()); return phi;")
    .Define("simPt","RVec<float> pt;for(auto& mom:simMom) pt.push_back(mom.pt()); return pt;")
    .Define("simEta","RVec<float> eta;for(auto& mom:simMom) eta.push_back(mom.eta()); return eta;")
    .Define("simPhi","RVec<float> phi;for(auto& mom:simMom) phi.push_back(mom.phi()); return phi;")
    .Define("eventId", "evtId") // for progress indicator 
    .Filter("eventId>=0") // at least one filter is mandatory!!!
    .Filter("vtxChi2>0.0001")
  ;
  
  varPatterns=
  {
    "b",                                             // kEvent
    "(fhcal|scwall)Mod(X|Y|Phi|E|Id|InSub.)",        // kChannel 
    "tr(Pt|Eta|Phi|BetaTof400|BetaTof700|SimIndex)", // kRecParticle  
    "",//"sim(Pt|Eta|Phi|Pdg|MotherId)"                   // kSimParticle  
  };

  return dd; 
}

void setupQvectors()
{
  vector<Qn::AxisD> corrAxesEvent=
  {
    {"b", 4,0,10},
  };

  vector<Qn::AxisD> corrAxesParticle=
  {
    {"trPt",4,0,2},
    {"trEta",4,0,2},
  };

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
  
  man.AddDetector("fhcal1", channel, "fhcalModPhi", "fhcalModE", {}, {1}, sumW);
  man.AddCorrectionOnQnVector("fhcal1", recentering);
  man.AddCorrectionOnQnVector("fhcal1", twistRescale);
  man.SetOutputQVectors("fhcal1", {plain, recentered});
  man.AddCutOnDetector("fhcal1", {"fhcalModInSub1"}, equal(1), "fhcal1");
  man.AddHisto2D("fhcal1", {{"fhcalModId", 100, 0., 100}, {"fhcalModE", 100, 0., 10}}, "fhcalModInSub1");
  man.AddHisto2D("fhcal1", {{"fhcalModX", 100, -100, 100}, {"fhcalModY", 100, -100, 100}}, "fhcalModInSub1");
  
  man.AddDetector("fhcal2", channel, "fhcalModPhi", "fhcalModE", {}, {1}, sumW);
  man.AddCorrectionOnQnVector("fhcal2", recentering);
  man.AddCorrectionOnQnVector("fhcal2", twistRescale);
  man.SetOutputQVectors("fhcal2", {plain, recentered});
  man.AddCutOnDetector("fhcal2", {"fhcalModInSub2"}, equal(1), "fhcal2");
  man.AddHisto2D("fhcal2", {{"fhcalModId", 100, 0., 100}, {"fhcalModE", 100, 0., 10}}, "fhcalModInSub2");
  man.AddHisto2D("fhcal2", {{"fhcalModX", 100, -100, 100}, {"fhcalModY", 100, -100, 100}}, "fhcalModInSub2");
  
  man.AddDetector("fhcal3", channel, "fhcalModPhi", "fhcalModE", {}, {1}, sumW);
  man.AddCorrectionOnQnVector("fhcal3", recentering);
  man.AddCorrectionOnQnVector("fhcal3", twistRescale);
  man.SetOutputQVectors("fhcal3", {plain, recentered});
  man.AddCutOnDetector("fhcal3", {"fhcalModInSub3"}, equal(1), "fhcal3");
  man.AddHisto2D("fhcal3", {{"fhcalModId", 100, 0., 100}, {"fhcalModE", 100, 0., 10}}, "fhcalModInSub3");
  man.AddHisto2D("fhcal3", {{"fhcalModX", 100, -100, 100}, {"fhcalModY", 100, -100, 100}}, "fhcalModInSub3");

  man.AddDetector("tr", track, "trPhi", "Ones", corrAxesParticle, {1,2}, sumW);
  man.AddCutOnDetector("tr", {"particleType"}, equal(kRecParticle), "recParticle");
//  man.AddCutOnDetector("tr", {"simPdg"}, [](float pdg){return int(pdg)==2212;}, "proton");
  man.AddCorrectionOnQnVector("tr", recentering);
  man.AddCorrectionOnQnVector("tr", twistRescale);
  man.SetOutputQVectors("tr", {plain, recentered, twisted, rescaled});
  man.AddHisto1D("tr", {"trPhi", 100, -3.15, 3.15}, "Ones");
//  man.AddHisto2D("tr", {{"trEta", 100, 0., 6.}, {"trPt",  100, 0., 3.}}, "Ones");
}
