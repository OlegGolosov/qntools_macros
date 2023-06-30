void shiftGraphsX(TMultiGraph *mg, float dx=0.0) 
{
  if (!mg->GetListOfGraphs() || dx==0.) return;
  for (int ig=1; ig < mg->GetListOfGraphs()->GetEntries(); ++ig) 
  {
    auto g=dynamic_cast<TGraphMultiErrors*>(mg->GetListOfGraphs()->At(ig));
    auto *x=g->GetX();
    auto n=g->GetN();
    float shift=ig*dx*(x[n-1]-x[0]);
    for (int ix=0; ix<n; ++ix)
      x[ix]+=shift;
  }
}

map<string,string> axisTitles=
{
  {"R1","R_{1}"},
  {"v1","v_{1}"},
  {"ps","#LTQQ#GT"},
  {"pr","#LTQQ#GT"},
  {"pi","#LTQQ#GT"},
  {"trPt","p_{T} (GeV/#it{c})"},
  {"trY","#it{y}"},
  {"centrality", "Centrality, %"},
  {"fit_fpol3", "dv_{1}/d#it{y}"},
  {"fit_fpol3star", "dv_{1}/d#it{y}"},
};

void plotSystematics(string fInNameList="syst.list", string fOutName="out.root")
{
  float systWidth=0.005, graphShift=0.005;
  vector<TFile*> files;
  ifstream inFileNames(fInNameList);
  string line;
  while(getline(inFileNames, line))
    if(line.at(0)!='#')
      files.push_back(new TFile(line.c_str())); 
  vector<TList*> lists;
  for (auto &f:files)
    lists.push_back(f->GetListOfKeys());
  TFile fOut(fOutName.c_str(), "recreate");
  auto nObj=lists.at(0)->GetSize();
  for(int i=0; i<nObj;i++)
  {
    auto mgIn=(TMultiGraph*)files.at(0)->Get(lists.at(0)->At(i)->GetName());
    cout << mgIn->GetName() << endl;
    auto ng=mgIn->GetListOfGraphs()->GetSize();
    auto mgOut=new TMultiGraph(mgIn->GetName(), mgIn->GetTitle());
    for (int g=0;g<ng;g++)
    {
      auto gIn=(TGraphErrors*)mgIn->GetListOfGraphs()->At(g);
      int np=gIn->GetN();
      vector<double> systErr(np, 0);
      vector<double> ex(np, systWidth*(gIn->GetPointX(np-1)-gIn->GetPointX(0)));
      for(int f=1; f<lists.size();f++)
      {
        auto mgVar=(TMultiGraph*)files.at(f)->Get(lists.at(f)->At(i)->GetName());
        auto gVar=(TGraphErrors*)mgVar->GetListOfGraphs()->At(g);
        for(int p=0;p<np;p++)
        {
          auto delta=fabs(gIn->GetPointY(p)-gVar->GetPointY(p));
          auto deltaSqErr=fabs(gIn->GetErrorY(p)*gIn->GetErrorY(p)-gVar->GetErrorY(p)*gVar->GetErrorY(p));
          auto significance=delta/sqrt(deltaSqErr);
          if(systErr.at(p)<delta)//max
            systErr.at(p)=delta;//max
//          if(significance>0.7)//rms
//            systErr.at(p)+=delta*delta;//rms
        }
      }
//      for(int p=0;p<np;p++)//rms
//        systErr.at(p)=sqrt(systErr.at(p));//rms

      auto gOut=new TGraphMultiErrors(gIn->GetN(), gIn->GetX(), gIn->GetY(), ex.data(), ex.data(), gIn->GetEY(), gIn->GetEY());
      string gTitle=gIn->GetTitle();
      gTitle=regex_replace(gTitle,regex("pionneg"),"#pi^{-}");
      gTitle=regex_replace(gTitle,regex("pionpos"),"#pi^{+}");
      gOut->SetTitle(gTitle.c_str());
      gOut->SetName(gIn->GetName());
      string titleX=gIn->GetXaxis()->GetTitle();
      string titleY=gIn->GetYaxis()->GetTitle();
      if (titleX.size()>0) gOut->GetXaxis()->SetTitle(axisTitles.at(titleX).c_str());
      if (titleY.size()>0) gOut->GetYaxis()->SetTitle(axisTitles.at(titleY).c_str());
      gOut->SetName(gIn->GetName());
      gOut->AddYError(systErr.size(), systErr.data(), systErr.data());
      gOut->SetLineColor(gIn->GetLineColor());
      gOut->SetLineStyle(gIn->GetLineStyle());
      gOut->SetLineWidth(gIn->GetLineWidth());
      gOut->SetMarkerColor(gIn->GetMarkerColor());
      gOut->SetMarkerStyle(gIn->GetMarkerStyle());
      gOut->SetMarkerSize(gIn->GetMarkerSize());
      gOut->SetFillColor(gIn->GetFillColor());
      gOut->SetFillStyle(gIn->GetFillStyle());
      gOut->GetAttLine(0)->SetLineColor(gIn->GetLineColor());
      //gOut->GetAttLine(0)->SetLineWidth(gIn->GetLineWidth());
      gOut->GetAttLine(0)->SetLineWidth(2);
      gOut->GetAttLine(1)->SetLineWidth(1);
      gOut->GetAttLine(1)->SetLineColor(gIn->GetLineColor());
      gOut->GetAttFill(1)->SetFillColor(gIn->GetLineColor());
      gOut->GetAttFill(1)->SetFillColorAlpha(gIn->GetLineColor(), 0.9);
      gOut->GetAttFill(1)->SetFillStyle(0);
      mgOut->Add(gOut, "PS ; Z; 5 s=1");
      shiftGraphsX(mgOut, graphShift);
    }
    string titleX=mgIn->GetXaxis()->GetTitle();
    string titleY=mgIn->GetYaxis()->GetTitle();
    if (titleX.size()>0) mgOut->GetXaxis()->SetTitle(axisTitles.at(titleX).c_str());
    if (titleY.size()>0) mgOut->GetYaxis()->SetTitle(axisTitles.at(titleY).c_str());
    mgOut->Write();
//    if(i==13)
//    {
//      mgOut->GetYaxis()->SetRangeUser(-0.08,0.04);
//      mgOut->Draw("a");
//      break;
//    }
  }
  fOut.Close();
}
