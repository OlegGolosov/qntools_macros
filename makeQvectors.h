using filteredDF=ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void>;
using definedDF=ROOT::RDF::RInterface<ROOT::Detail::RDF::RLoopManager, void>;
using ROOT::VecOps::RVec;

enum varType
{
  kEvent=0,
  kChannel,
  kRecParticle,
  kSimParticle,
  kNVarTypes
};

string qaFilePath;
Qn::CorrectionManager man; 
TFile *outFile;
TTree *outTree;
short particleTypePosition;

vector <string> varPatterns(kNVarTypes);
vector <vector<string>> varNames(kNVarTypes);
vector<short> vcInitialPositions(kNVarTypes, 0);

void InitVariables(filteredDF &d, vector <string> &patterns, vector<vector<string>> &varNames) 
{
  // Add all needed variables
  short ivar=0;
  for (int type=0;type<kNVarTypes;type++)
  {
    vector<short> varSizes;
    for (auto &name : d.GetColumnNames()) 
    {
      short size{1};
      auto columnType=d.GetColumnType(name);
      if (type==kEvent && columnType.find("RVec")<columnType.size()) 
        continue;
      if (!regex_match(name, regex(patterns.at(type))))
        continue;
      if (type==kChannel)
        size=*d.Range(0,0).Define("n", name+".size()").Mean("n");
      varNames.at(type).push_back(name);
      varSizes.push_back(size);
        
    }

    vcInitialPositions.at(type)=ivar;
    cout << "\nPattern: " << patterns.at(type) << endl;
    for (int i=0;i<varNames.at(type).size();i++)
    {
      auto name=varNames.at(type).at(i);
      auto size=varSizes.at(i);
      man.AddVariable(name, ivar, size);
      printf("%s: %i+%i\n", name.c_str(), ivar, size);
      ivar+=size;
      if (type==kEvent)
        man.AddEventVariable(name);
    }
  }
  particleTypePosition=ivar;
  printf("%s: %i+%i\n", "particleType", particleTypePosition, 1);
  man.AddVariable("particleType", particleTypePosition, 1);
}

void DefineVariableFilling(filteredDF &d, vector <vector<string>> &varNames)
{
  string vcEventFillExpr="vector<float> vcEvent; ";
  for (auto& var:varNames.at(kEvent))
    vcEventFillExpr+=Form("vcEvent.push_back(%s); ",var.c_str());
  vcEventFillExpr+="return vcEvent;";
  cout << endl << vcEventFillExpr << endl; 
  
  string vcModuleFillExpr="vector<float> vcModule; ";
  for (auto& var:varNames.at(kChannel))
    vcModuleFillExpr+=Form("for(auto& val:%s) vcModule.push_back(val); ",var.c_str());
  vcModuleFillExpr+="return vcModule;";
  cout << endl << vcModuleFillExpr << endl; 
  
  string nRecPart="0";
  if (varNames.at(kRecParticle).size()>0)
    nRecPart=varNames.at(kRecParticle).front() + ".size()";
  string vcRecPartFillExpr=Form("vector<vector<float>> vcRecParticle(%s); for(int i=0;i<vcRecParticle.size();i++) { ", nRecPart.c_str());
  for (auto& var:varNames.at(kRecParticle))
    vcRecPartFillExpr+=Form("vcRecParticle.at(i).push_back(%s.at(i)); ", var.c_str());
  vcRecPartFillExpr+=" } return vcRecParticle;";
  cout << endl << vcRecPartFillExpr << endl; 
  
  string nSimPart="0";
  if (varNames.at(kSimParticle).size()>0)
    nSimPart=varNames.at(kSimParticle).front() + ".size()";
  string vcSimPartFillExpr=Form("vector<vector<float>> vcSimParticle(%s); for(int i=0;i<vcSimParticle.size();i++) { ", nSimPart.c_str());
  for (auto& var:varNames.at(kSimParticle))
    vcSimPartFillExpr+=Form("vcSimParticle.at(i).push_back(%s.at(i)); ", var.c_str());
  vcSimPartFillExpr+=" } return vcSimParticle;";
  cout << endl << vcSimPartFillExpr << endl; 
  
  d=d.Define("vcEvent", vcEventFillExpr)
     .Define("vcModule", vcModuleFillExpr)
     .Define("vcRecParticle", vcRecPartFillExpr)
     .Define("vcSimParticle", vcSimPartFillExpr)
  ;
//  d.Display({/*"vcEvent","vcModule","vcRecParticle",*/"vcSimParticle"},2)->Print();
}

void processEvent(const ULong64_t eventId, const vector<float> vcEvent, const vector<float> vcModule, const vector<vector<float>> vcRecParticle, const vector<vector<float>> vcSimParticle)
{
  if (eventId%1000==0)
    cout << "\r" << eventId;
  man.Reset();
  double *vc = man.GetVariableContainer();
  for(int ivar=0; ivar<vcEvent.size();ivar++)
    vc[vcInitialPositions.at(kEvent)+ivar]=vcEvent.at(ivar); 
  for(int ivar=0; ivar<vcModule.size();ivar++)
    vc[vcInitialPositions.at(kChannel)+ivar]=vcModule.at(ivar); 
  man.ProcessEvent();
  man.FillChannelDetectors();
  
  vc[particleTypePosition]=kRecParticle; 
  for(int ipart=0;ipart<vcRecParticle.size();ipart++)
  {
    for(int ivar=0;ivar<vcRecParticle.at(0).size();ivar++)
      vc[vcInitialPositions.at(kRecParticle)+ivar]=vcRecParticle.at(ipart).at(ivar);
    man.FillTrackingDetectors();
  }
  
  vc[particleTypePosition]=kSimParticle; 
  for(int ipart=0;ipart<vcSimParticle.size();ipart++)
  {
    for(int ivar=0;ivar<vcSimParticle.at(0).size();ivar++)
      vc[vcInitialPositions.at(kSimParticle)+ivar]=vcSimParticle.at(ipart).at(ivar);
    man.FillTrackingDetectors();
  }
  man.ProcessCorrections();
}

void init(filteredDF &d, string outFilePath, string calibFilePath)
{
  outFile = TFile::Open(outFilePath.c_str(), "recreate");
  outFile->cd();
  outTree = new TTree("tree", "tree");
  //qaFilePath=regex_replace(outFilePath, regex("qn\\.root"), "qa.root");
  qaFilePath=regex_replace(outFilePath, regex("qn\\.root"), "qa.root");
  man.SetCalibrationInputFileName(calibFilePath);
  man.SetFillOutputTree(true);
  man.SetFillCalibrationQA(true);
  man.SetFillValidationQA(true);
  man.ConnectOutputTree(outTree);
  
  InitVariables(d, varPatterns, varNames);
}

void run(filteredDF &d) {
  cout << "Number of Events: " << *(d.Count()) << endl;
  man.InitializeOnNode();
  man.SetCurrentRunName("test");
  DefineVariableFilling(d, varNames);
  d.Foreach(processEvent, {"rdfentry_", "vcEvent", "vcModule", "vcRecParticle", "vcSimParticle"});
  cout << endl;
  
  man.Finalize();
  outFile->cd();
  outTree->Write("tree");
  outFile->Close();

  TFile qaFile(qaFilePath.c_str(), "RECREATE");
  man.GetCorrectionQAList()->Write("CorrectionQAHistograms", TObject::kSingleKey);
  man.GetCorrectionList()->Write("CorrectionHistograms", TObject::kSingleKey);
  qaFile.Close();
}
