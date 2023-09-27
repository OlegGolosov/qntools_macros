const char *boldRed="\033[1;31m%s\033[0m"; // https://stackoverflow.com/questions/2616906/how-do-i-output-coloured-text-to-a-linux-terminal
                                           
map<string,string> axisTitles=
{
  {"R1","R_{1}"},
  {"v1","v_{1}"},
  {"ps","#LTQQ#GT"},
  {"pr","#LTQQ#GT"},
  {"pi","#LTQQ#GT"},
  {"trY","#it{y}"},
  {"trPt","p_{T} (GeV/#it{c})"},
  {"trY0","#it{y_{0}}"},
  {"trUt0","u_{T0}"},
  {"centrality", "Centrality (%)"},
  {"fit_fpol3", "dv_{1}/d#it{y}"},
  {"fit_fpol3star", "dv_{1}/d#it{y}"},
};

map <string, string> legendReplaceMap
{
  {"proton", "p"},
  {"pionpos", "#pi^{+}"},
  {"pionneg", "#pi^{-}"},
  {"(\\d{1,2}-\\d\\d)", "$&%"},
};

const vector <int> colors = 
{
  kBlue, kRed, kGreen+2, kMagenta+2, kBlack, kOrange+2, kPink+2,
  kTeal+2, kCyan+2, kCyan+4, kAzure, kGray+3, kOrange+7, kGreen+4
};

const vector <vector <int>> markerStyles = 
{ 
  {20, 21, 22, 23, 29, 33, 34, 39, 41, 43, 45, 47, 48, 49}, // filled
  {24, 25, 26, 32, 30, 27, 28, 37, 40, 42, 44, 46, 35, 36},  // empty
  { 1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1},  // none
  { 1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1},  // none
};

const vector <vector <int>> lineWidths = 
{
  {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
  {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
  {3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3},
  {3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3},
};

const vector <vector <int>> lineStyles = 
{
  {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
  {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
  {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
  {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
};

const vector <vector <string>> drawOptions = 
{
  {"p", "p", "p", "p", "p", "p", "p", "p", "p", "p", "p", "p", "p", "p"},
  {"p", "p", "p", "p", "p", "p", "p", "p", "p", "p", "p", "p", "p", "p"},
  {"l", "l", "l", "l", "l", "l", "l", "l", "l", "l", "l", "l", "l", "l"},
  {"l", "l", "l", "l", "l", "l", "l", "l", "l", "l", "l", "l", "l", "l"},
};

float markerSize = 1.5;
float graphShift = 0.003;
float xAxisTitleSize = 0.05;
float yAxisTitleSize = 0.05;

map <string,Qn::DataContainerStatCalculate> corrMap;
map <string,TGraphErrors*> graphMap;

template <typename T>
vector<vector<string>> findMatches(string pattern, map<string, T> tMap)
{
  if(verbose) cout << pattern << ":"  << endl;
  regex re(pattern);
  vector<vector<string>> matches;
  for(auto &p:tMap)
  {
    smatch sm;
    if(regex_match(p.first, sm, regex(pattern)))
    {
      vector<string> captures;
      for (auto c:sm) captures.push_back(c); 
      matches.push_back(captures);
      if(verbose) cout << captures.front() << endl;
    }
  }
  return matches;
}

void R1_3S(string name, string c1, string c2, string c3, Qn::AxisD rebinAxis={"",{}})
{
  vector<Qn::DataContainerStatCalculate> corrs={corrMap.at(c1), corrMap.at(c2), corrMap.at(c3)};
  for(auto &corr:corrs)
  {
    for(auto &c:corr) c.SetWeightType(Qn::Stat::WeightType::REFERENCE);
    if(!rebinAxis.Name().empty()) corr=corr.Rebin(rebinAxis);
  }
  auto result=Sqrt(2*corrs.at(0)*corrs.at(1)/corrs.at(2));
  if (!corrMap.emplace(name, result).second)
    cout << __func__ << ": ERROR! correlation already defined: " << name << endl;
}

void R1_4S(string name, string c1, string c2, float scale=1., Qn::AxisD rebinAxis={"",{}})
{
  vector<Qn::DataContainerStatCalculate> corrs={corrMap.at(c1), corrMap.at(c2)};
  for(auto &corr:corrs)
  {
    for(auto &c:corr) c.SetWeightType(Qn::Stat::WeightType::REFERENCE);
    if(!rebinAxis.Name().empty()) corr=corr.Rebin(rebinAxis);
  }
  auto result=scale*2*corrs.at(0)/corrs.at(1);
  if (!corrMap.emplace(name, result).second)
    cout << __func__ << ": ERROR! correlation already defined: " << name << endl;
}

void v1(string name, string vObs, string res)
{
  if (!corrMap.count(vObs) || !corrMap.count(res)) return;
  auto result=2.*corrMap.at(vObs)/corrMap.at(res);
  if (!corrMap.emplace(name, result).second)
    cout << __func__ << ": ERROR! correlation already defined: " << name << endl;
}

void v2(string name, string vObs, string res1, string res2)
{
  if (!corrMap.count(vObs) || !corrMap.count(res1) || !corrMap.count(res2)) return;
  auto result=4.*corrMap.at(vObs)/corrMap.at(res1)/corrMap.at(res2);
  if (!corrMap.emplace(name, result).second)
    cout << __func__ << ": ERROR! correlation already defined: " << name << endl;
}

//void addGraphsFromFile(string fileName, string postfix)
//{
//  TFile f(fileName.c_str());
//  string graphName=obj->GetName();
//  auto g=new TGraphErrors(path.c_str());
//  g->SetName(graphName.c_str());
//  g->SetTitle(graphName.c_str());
//  if (!graphMap.emplace(graphName, g).second)
//    cout << __func__ << ": ERROR! graph already defined: " << graphName << endl;
//}

void addGraphFromTxt(string pattern)
{
  auto wd=gSystem->GetWorkingDirectory();
  regex re(pattern);
  TSystemDirectory dir(wd.c_str(), wd.c_str());
  smatch sm;
  for(const auto& obj:*dir.GetListOfFiles())
  {
    string path=Form("%s%s", obj->GetTitle(), obj->GetName());
    if(regex_match(path, sm, re))
    {
      string graphName=obj->GetName();
      auto g=new TGraphErrors(path.c_str());
      g->SetName(graphName.c_str());
      g->SetTitle(graphName.c_str());
      if (!graphMap.emplace(graphName, g).second)
        cout << __func__ << ": ERROR! graph already defined: " << graphName << endl;
    }
  }
}

void scale(string pattern, double scale)
{
  if(verbose) cout << endl << __func__ << " " << scale << " ";
  auto matches=findMatches<Qn::DataContainerStatCalculate>(pattern, corrMap);
  for (auto &m:matches)
  {
    auto name=m.at(0);
    corrMap.at(name)=corrMap.at(name)*scale;
  }
}

void changeBinEdges(string pattern, string axisName, std::function<double(double)> func, string axisNameNew="")
{
  if(verbose) cout << endl << __func__ << " " << pattern << " " << axisName << "->" << axisNameNew << "";
  auto matches=findMatches<Qn::DataContainerStatCalculate>(pattern, corrMap);
  for (auto &m:matches)
  {
    auto &corr=corrMap.at(m.at(0));
    auto &axis=corr.GetAxis(axisName);
    for (auto &x:axis.GetBinEdges())
     x=func(x); 
    if(!axisNameNew.empty())
      axis.SetName(axisNameNew);
  }
}

void project(string pattern, vector<string> axes)
{
  if(verbose) cout << endl << __func__ << " ";
  auto matches=findMatches<Qn::DataContainerStatCalculate>(pattern, corrMap);
  for (auto &m:matches)
  {
    auto name=m.at(0);
    auto corr=corrMap.at(name);
    name+="_proj";
    for (auto& a:axes)
      name=name+"_"+a;
    if (!corrMap.emplace(name, corr.Projection(axes)).second)
      cout << __func__ << ": ERROR! correlation already defined: " << name << endl;
  }
}

void rebin(string pattern, Qn::AxisD axis, string comment)
{
  if(verbose) cout << endl << __func__ << " " << axis.Name() << " " << comment << " ";
  auto matches=findMatches<Qn::DataContainerStatCalculate>(pattern, corrMap);
  for (auto &m:matches)
  {
    auto name=m.at(0);
    auto corr=corrMap.at(name);
    name=name+"_rebin_"+axis.Name()+"_"+comment;
    if (!corrMap.emplace(name, corr.Rebin(axis)).second)
      cout << __func__ << ": ERROR! correlation already defined: " << name << endl;
  }
}

void unfold(string pattern, string axisName, string format="%02.0f")
{
  if(verbose) cout << endl << __func__ << " " << axisName << " ";
  auto matches=findMatches<Qn::DataContainerStatCalculate>(pattern, corrMap);
  for (auto &m:matches)
  {
    auto name=m.at(0);
    auto corr=corrMap.at(name);
    name=name+"_unfold_"+axisName+"_";
    auto axis=corr.GetAxis(axisName);
    for (int i=0;i<axis.GetNBins();i++)
    {
      auto perBinName=name;
      string range=format+"-"+format;
      perBinName.append(Form(range.c_str(), axis.GetLowerBinEdge(i), axis.GetUpperBinEdge(i)));
      Qn::AxisD perBinAxis(axis.Name(), {axis.GetLowerBinEdge(i), axis.GetUpperBinEdge(i)});
      if(verbose) cout << perBinName << "\t";
      if (!corrMap.emplace(perBinName, corr.Rebin(perBinAxis)).second)
        cout << __func__ << ": ERROR! correlation already defined: " << perBinName << endl;
    }
  }
}

void merge(string pattern) //TODO: try solving with one regexp: merge those with captures, save to different containers without captures(e.g. different centralities) 
{
  if(verbose) cout << endl << __func__ << " ";
  auto matches=findMatches<Qn::DataContainerStatCalculate>(pattern, corrMap);
  if (matches.size()==0) return;
  Qn::DataContainerStatCalculate merged;
  string mergedName=regex_replace(pattern, regex("\\(|\\)"), "");
  mergedName=regex_replace(mergedName, regex("\\|"), "_");
  for(int i=0; i<matches.size(); i++)
  {
    auto name=matches.at(i).at(0);
    if(i==0)
      merged=corrMap.at(name);
    else
      merged=merged.Apply(corrMap.at(name), [](const Qn::StatCalculate &lhs, const Qn::StatCalculate &rhs) { return Qn::Merge(lhs, rhs); });
  }
  if (!corrMap.emplace(mergedName, merged).second)
    cout << __func__ << ": ERROR! correlation already defined: " << mergedName << endl;
}

void ShiftGraphsX(TMultiGraph *mg, float dx=0.0) 
{
  if (!mg->GetListOfGraphs() || dx==0.) return;
  for (int ig=1; ig < mg->GetListOfGraphs()->GetEntries(); ++ig) 
  {
    auto g=dynamic_cast<TGraph*>(mg->GetListOfGraphs()->At(ig));
    auto *x=g->GetX();
    auto n=g->GetN();
    float shift=ig*dx*(x[n-1]-x[0]);
    for (int ix=0; ix<n; ++ix)
      x[ix]+=shift;
  }
}
     
TGraphErrors *makeGraph(string name)
{
  auto corr=corrMap.at(name);
  auto g=ToTGraph(corr);
  if(!g) return nullptr;
  auto xname=corr.GetAxes().at(0).Name();
  auto yname=name.substr(0,2);
  if(axisTitles.count(xname)) xname=axisTitles.at(xname);
  if(axisTitles.count(yname)) yname=axisTitles.at(yname);
  
  g->GetXaxis()->SetTitle(xname.c_str());
  g->GetYaxis()->SetTitle(yname.c_str());
  g->SetName(name.c_str());
  return g;
}

void fit(string pattern, const char *funcName, const char *option="qr")
{
  if(verbose) cout << endl << __func__ << " ";
  auto matches=findMatches<Qn::DataContainerStatCalculate>(pattern, corrMap);
  string postfix=Form("_fit_%s", funcName);
  for(auto &m:matches)
  {
    auto name=m.at(0);
    //auto g=ToTGraph(corrMap.at(name));
    auto g=makeGraph(name);
    g->Fit(funcName, option, "");
//    g->GetXaxis()->SetTitle(corrMap.at(name).GetAxes().at(0)->Name().c_str());
    if (!graphMap.emplace(name+postfix, g).second)
      cout << __func__ << ": ERROR! fit already defined: " << name << endl;
  }
}

void projectFitParameters(string pattern)
{
  if(verbose) cout << endl << __func__ << " ";
  auto matches=findMatches<TGraphErrors*>(pattern, graphMap);
  vector<double> xPoints;
  int nPars=0;
  string xname, yname;
  for(auto &match:matches)
  {
    auto g=graphMap.at(match.at(0));
    if(g->GetListOfFunctions()->GetSize()>0)
    {
      nPars=((TF1*)g->GetListOfFunctions()->At(0))->GetNumberFreeParameters();
      xname=match.at(1);
      yname=match.at(4);
      if(axisTitles.count(xname)) xname=axisTitles.at(xname);
      if(axisTitles.count(yname)) yname=axisTitles.at(yname);
      break;
    }
  }
  vector<vector<double>> parameters(nPars);
  vector<vector<double>> parErrors(nPars);
  for(auto &match:matches)
  {
    auto g=graphMap.at(match.at(0));
    if(g->GetListOfFunctions()->GetSize()>0)
    {
      auto f=(TF1*)g->GetListOfFunctions()->At(0);
      auto xmin=atof(match.at(2).c_str());
      auto xmax=atof(match.at(3).c_str());
      xPoints.push_back(0.5*(xmin+xmax));
      for (int i=0; i<nPars;i++)
      {
        parameters.at(i).push_back(f->GetParameters()[i]);
        parErrors.at(i).push_back(f->GetParErrors()[i]);
      }
    }
  }
  auto gname=regex_replace(pattern, regex("_\\(\\.\\*\\)-\\(\\.\\*\\)"), "");
  gname=regex_replace(gname, regex("\\(|\\||\\)"), "");
  for (int i=0; i<nPars;i++)
  {
    auto gnamei=gname+Form("_par%i", i);
    auto g=new TGraphErrors(parameters.at(i).size(), &xPoints[0], &parameters.at(i)[0], 0, &parErrors.at(i)[0]);
    g->GetXaxis()->SetTitle(xname.c_str());
    g->GetYaxis()->SetTitle(yname.c_str());
    if (!graphMap.emplace(gnamei, g).second)
      cout << __func__ << ": ERROR! graph already defined: " << gnamei << endl;
  }
}

string makeCommonName(vector<string> patterns={"R1_(psd.)_(4S)_(X)", "R1_(psd.)_(4S)_(Y)", "R1_(psd.)_(3S)_(X)", "R1_(psd.)_(3S)_(Y)"})
{
  regex r1("_?(\\([^)]*\\))");
  regex r2("(\\([^)]*\\))");
  string firstPattern;
  for (auto &p:patterns)
    if(p.size()>0)
    {
      firstPattern=p;
      break;
    }
  string common=regex_replace(firstPattern, r1, "");
  string diff="";
  for(auto &p:patterns)
  {
    diff+="_";
    smatch sm;
    while(regex_search(p, sm, r2))
    {
      string capture=sm[1];
      p.erase(p.find(capture),capture.size());
      diff+=capture;
    }
  }
  return common+diff;
}

TMultiGraph* makeMG(vector<string> patterns)
{
  auto commonName=makeCommonName(patterns);
  auto mg=new TMultiGraph(commonName.c_str(), commonName.c_str());
  for(int i=0; i<patterns.size(); i++)
  {
    string pattern=patterns.at(i);
    if(!verbose) cout << Form(boldRed, pattern.c_str()) << ":\t";
    auto matches=findMatches<Qn::DataContainerStatCalculate>(pattern, corrMap);
    bool corrFound=matches.size();
    if(!corrFound)
      matches=findMatches<TGraphErrors*>(pattern, graphMap);

    for(int j=0; j<matches.size(); j++)
    {
      TGraphErrors *g;
      auto match=matches.at(j);
      auto name=match.at(0);
      if(!verbose) cout << name << "\t";
      if(corrFound)
        g=makeGraph(name);
      else
        g=graphMap.at(name);
      if (!g) continue;
      string gtitle="";
      vector<string> captures(match.begin()+1, match.end());
      for(auto &capture:captures)
        gtitle=gtitle+capture+" ";
      for(auto &r:legendReplaceMap)
	gtitle=regex_replace(gtitle, regex(r.first), r.second);
      cout << "gtitle: " << gtitle << endl;
      g->SetTitle(gtitle.c_str());
      g->SetLineColor(colors.at(j));
      g->SetLineStyle(lineStyles.at(i).at(j));
      g->SetLineWidth(lineWidths.at(i).at(j));
      g->SetMarkerColor(colors.at(j));
      g->SetMarkerSize(markerSize);
      g->SetMarkerStyle(markerStyles.at(i).at(j));
      if(g->GetListOfFunctions()->GetSize()>0)
      {
        auto f=(TF1*)g->GetListOfFunctions()->At(0);
        f->SetLineColor(colors.at(j));
        f->SetLineStyle(lineStyles.at(i).at(j));
      }
      mg->Add(g, drawOptions.at(i).at(j).c_str());
    }
  }
  cout << endl;
  if(mg->GetListOfGraphs())
  {
    ShiftGraphsX(mg, graphShift);
    auto g=(TGraphErrors*)mg->GetListOfGraphs()->At(0);
    mg->GetXaxis()->SetTitle(g->GetXaxis()->GetTitle());
    mg->GetYaxis()->SetTitle(g->GetYaxis()->GetTitle());
  }
  return mg;
}

void drawMG(vector<string> patterns, const char* option="apl") {makeMG(patterns)->Draw(option);}
void writeMG(vector<string> patterns) {makeMG(patterns)->Write();}
void drawGraph(string name, const char* option="apl") {makeGraph(name)->Draw(option);}
void writeGraph(string name) {makeGraph(name)->Write();}
