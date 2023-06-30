TGraphErrors* flow(TProfile3D *p3orig, float min, float max, float centMin=10, float centMax=30, string profileAxis="X", int color=kBlue)
{
  TString name=Form("%s_%s_%2.1f_%2.1f_%.0f_%.0f", p3orig->GetName(), profileAxis.c_str(), min, max, centMin, centMax);
  auto p3 = (TProfile3D*)p3orig->Clone("p3_"+name);
  p3->GetZaxis()->SetRangeUser(centMin, centMax);
  TProfile2D *p2=p3->Project3DProfile("yx");
  TProfile *p;
  if(profileAxis=="X")
    p=p2->ProfileX("p_"+name, p2->GetYaxis()->FindBin(min), p2->GetYaxis()->FindBin(max));
  if(profileAxis=="Y")
    p=p2->ProfileY("p_"+name, p2->GetXaxis()->FindBin(min), p2->GetXaxis()->FindBin(max));
  vector<float> x, y, ey;
  for (int i=1;i<=p->GetNbinsX();i++)
  {
    x.push_back(p->GetXaxis()->GetBinCenter(i));
    y.push_back(p->GetBinContent(i));
    ey.push_back(p->GetBinError(i));
//    cout << i << "\t" << x.back() << "\t" << y.back() << "\t" << ey.back() << endl;
  }
  auto g=new TGraphErrors(x.size()-1, x.data(), y.data(), 0, ey.data());
  g->SetName("g_"+name);
  g->SetFillStyle(3002);
  g->SetFillColor(color);
  g->SetMarkerColor(color);
  g->SetLineColor(color);
  g->SetLineWidth(2);
  return g;  
}

void compareToModel(string dataFileName="data.root", string modelFileName="model.root", string compFileName="comp.root")
{
  smatch sm;
  regex_search(modelFileName, sm, regex("pbpb(..)_"));
  float pBeam=atof(string(sm[1]).c_str());
  map<string,vector<float>> axisMap=
  {
    {"proton_trY", {0.0, 0.6, 1.0, 1.2, 1.4, 1.6, 2.2}},
    {"proton_trPt", {0.0, 0.6, 1.0, 1.4, 1.8}},
    {"pionneg_trY", {0.0, 0.6, 1.0, 2.0}},
    {"pionneg_trPt", {0.0, 0.2, 0.4, 0.6, 1.0, 1.8}},
    {"pionpos_trY", {0.0, 0.6, 1.0, 2.0}},
    {"pionpos_trPt", {0.0, 0.2, 0.4, 0.6, 1.0, 1.8}},
  };
  TFile dataFile(dataFileName.c_str());
  TFile modelFile(modelFileName.c_str());
  TFile compFile(compFileName.c_str(), "recreate");
  vector <string> mgNames=
  {
    {"v1_pionneg_psd1_psd2_psd3_X_rebin_centrality_wide_unfold_centrality_00-10_rebin_trPt_pionneg_unfold_trPt_proj_trY_(.*)"},
    {"v1_pionneg_psd1_psd2_psd3_X_rebin_centrality_wide_unfold_centrality_00-10_rebin_trY_pionneg_unfold_trY_proj_trPt_(.*)"},
    {"v1_pionneg_psd1_psd2_psd3_X_rebin_centrality_wide_unfold_centrality_10-25_rebin_trPt_pionneg_unfold_trPt_proj_trY_(.*)"},
    {"v1_pionneg_psd1_psd2_psd3_X_rebin_centrality_wide_unfold_centrality_10-25_rebin_trY_pionneg_unfold_trY_proj_trPt_(.*)"},
    {"v1_pionneg_psd1_psd2_psd3_X_rebin_centrality_wide_unfold_centrality_25-45_rebin_trPt_pionneg_unfold_trPt_proj_trY_(.*)"},
    {"v1_pionneg_psd1_psd2_psd3_X_rebin_centrality_wide_unfold_centrality_25-45_rebin_trY_pionneg_unfold_trY_proj_trPt_(.*)"},
    {"v1_pionneg_psd1_psd2_psd3_X_rebin_centrality_wide_unfold_centrality_45-60_rebin_trPt_pionneg_unfold_trPt_proj_trY_(.*)"},
    {"v1_pionneg_psd1_psd2_psd3_X_rebin_centrality_wide_unfold_centrality_45-60_rebin_trY_pionneg_unfold_trY_proj_trPt_(.*)"},
    {"v1_pionpos_psd1_psd2_psd3_X_rebin_centrality_wide_unfold_centrality_00-10_rebin_trPt_pionpos_unfold_trPt_proj_trY_(.*)"},
    {"v1_pionpos_psd1_psd2_psd3_X_rebin_centrality_wide_unfold_centrality_00-10_rebin_trY_pionpos_unfold_trY_proj_trPt_(.*)"},
    {"v1_pionpos_psd1_psd2_psd3_X_rebin_centrality_wide_unfold_centrality_10-25_rebin_trPt_pionpos_unfold_trPt_proj_trY_(.*)"},
    {"v1_pionpos_psd1_psd2_psd3_X_rebin_centrality_wide_unfold_centrality_10-25_rebin_trY_pionpos_unfold_trY_proj_trPt_(.*)"},
    {"v1_pionpos_psd1_psd2_psd3_X_rebin_centrality_wide_unfold_centrality_25-45_rebin_trPt_pionpos_unfold_trPt_proj_trY_(.*)"},
    {"v1_pionpos_psd1_psd2_psd3_X_rebin_centrality_wide_unfold_centrality_25-45_rebin_trY_pionpos_unfold_trY_proj_trPt_(.*)"},
    {"v1_pionpos_psd1_psd2_psd3_X_rebin_centrality_wide_unfold_centrality_45-60_rebin_trPt_pionpos_unfold_trPt_proj_trY_(.*)"},
    {"v1_pionpos_psd1_psd2_psd3_X_rebin_centrality_wide_unfold_centrality_45-60_rebin_trY_pionpos_unfold_trY_proj_trPt_(.*)"},
    {"v1_proton_psd1_psd2_psd3_X_rebin_centrality_wide_unfold_centrality_00-10_rebin_trPt_proton_unfold_trPt_proj_trY_(.*)"},
    {"v1_proton_psd1_psd2_psd3_X_rebin_centrality_wide_unfold_centrality_00-10_rebin_trY_proton_unfold_trY_proj_trPt_(.*)"},
    {"v1_proton_psd1_psd2_psd3_X_rebin_centrality_wide_unfold_centrality_10-25_rebin_trPt_proton_unfold_trPt_proj_trY_(.*)"},
    {"v1_proton_psd1_psd2_psd3_X_rebin_centrality_wide_unfold_centrality_10-25_rebin_trY_proton_unfold_trY_proj_trPt_(.*)"},
    {"v1_proton_psd1_psd2_psd3_X_rebin_centrality_wide_unfold_centrality_25-45_rebin_trPt_proton_unfold_trPt_proj_trY_(.*)"},
    {"v1_proton_psd1_psd2_psd3_X_rebin_centrality_wide_unfold_centrality_25-45_rebin_trY_proton_unfold_trY_proj_trPt_(.*)"},
    {"v1_proton_psd1_psd2_psd3_X_rebin_centrality_wide_unfold_centrality_45-60_rebin_trPt_proton_unfold_trPt_proj_trY_(.*)"},
    {"v1_proton_psd1_psd2_psd3_X_rebin_centrality_wide_unfold_centrality_45-60_rebin_trY_proton_unfold_trY_proj_trPt_(.*)"},
  };

  const vector <int> colors = 
  {
    kBlue, kRed, kGreen+2, kMagenta+2, kBlack, kOrange+2, kPink+2,
    kTeal+2, kCyan+2, kCyan+4, kAzure, kGray+3, kOrange+7, kGreen+4
  };

  for (auto &name:mgNames)
  {
    auto found=regex_search(name, sm, regex("v1_(.*)_psd1_psd2_psd3_X_rebin_centrality_wide_unfold_centrality_(..)-(..)_rebin_.*_proj_(.*)_.*"));
    string part=sm[1], axis=sm[4];
    float cent1=atof(string(sm[2]).c_str()), cent2=atof(string(sm[3]).c_str());
    cout << name << endl;
//    cout << found << "\t" << part << "\t" << axis << "\t" << cent1 << "\t" << cent2 << endl;
    auto mg=(TMultiGraph*)dataFile.Get(name.c_str())->Clone();
    auto p3=(TProfile3D*)modelFile.Get(Form("v1_%s", part.c_str()));
    vector<float> bins;
    string profileAxis;
    if (axis=="trPt")
    {
      profileAxis="X";
      bins=axisMap.at(part+"_trY");
    }
    else if (axis=="trY")
    {
      profileAxis="Y";
      bins=axisMap.at(part+"_trPt");
    }
    for (int i=1;i<bins.size();i++)
    {
//      cout << bins.at(i-1) << "\t" << bins.at(i) << endl;
      mg->Add(flow(p3, bins.at(i-1), bins.at(i), cent1, cent2, profileAxis, colors.at(i-1)), "e3");
    }
    mg->Write();
  } 
  compFile.Close();
}
