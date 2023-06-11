void shiftGraphsX(TMultiGraph *mg, float dx=0.0) 
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

void plotSystematics(string fInNameList="syst.list", string fOutName="out.root")
{
  float systWidth=0.005, graphShift=0.005;
  vector<TFile*> files;
  ifstream inFileNames(fInNameList);
  string line;
  while(getline(inFileNames, line))
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
          if(significance>0.7)
            systErr.at(p)+=delta*delta;
        }
      }
      for(int p=0;p<np;p++)
        systErr.at(p)=sqrt(systErr.at(p));

      auto gOut=new TGraphMultiErrors(gIn->GetN(), gIn->GetX(), gIn->GetY(), ex.data(), ex.data(), gIn->GetEY(), gIn->GetEY());
      gOut->AddYError(systErr.size(), systErr.data(), systErr.data());
      gOut->SetLineColor(gIn->GetLineColor());
      gOut->SetLineStyle(gIn->GetLineStyle());
      gOut->SetLineWidth(gIn->GetLineWidth());
      gOut->SetMarkerColor(gIn->GetMarkerColor());
      gOut->SetMarkerStyle(gIn->GetMarkerStyle());
      gOut->SetMarkerSize(gIn->GetMarkerSize());
      gOut->SetFillColor(gIn->GetFillColor());
      gOut->SetFillStyle(gIn->GetFillStyle());
      gOut->GetAttLine(0)->SetLineColor(gOut->GetLineColor());
      //gOut->GetAttLine(0)->SetLineWidth(gOut->GetLineWidth());
      gOut->GetAttLine(0)->SetLineWidth(2);
      gOut->GetAttLine(1)->SetLineWidth(0);
      gOut->GetAttFill(1)->SetFillColor(gOut->GetLineColor());
      gOut->GetAttFill(1)->SetFillColorAlpha(gOut->GetLineColor(), 0.5);
      gOut->GetAttFill(1)->SetFillStyle(3001);
      mgOut->Add(gOut, "PS ; ; 5 s=1");
      shiftGraphsX(mgOut, graphShift);
    }
    mgOut->Write();
//    if(i==10)
//    {
//      mgOut->Draw("a");
//      break;
//    }
  }
  fOut.Close();
}
