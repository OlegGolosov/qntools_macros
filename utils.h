TChain* makeChain(string& fileName, const char* treename) {
  cout << "Adding files to chain:" << endl;
  TChain *chain = new TChain(treename);
  if (fileName.rfind(".root") < fileName.size())
    chain->Add(fileName.data());
  else {
    TFileCollection fc("fc", "", fileName.c_str());
    chain->AddFileInfoList((TCollection*)fc.GetList());
  }
  chain->ls();
  return chain;
}

inline std::function<bool(double)> range(double lo, double hi) 
  {return [=](double x) { return lo <= x && x <= hi; };}

inline std::function<bool(double)> rangeStrict(double lo, double hi) 
  {return [=](double x) { return lo < x && x < hi; };}

inline std::function<bool(double)> equal(int val) 
  {return [=](double _val) { return TMath::Nint(_val) == val; };}

inline bool is (double x) 
  {return fabs(x - 1) < 1e-5;}

inline bool isNot (double x) 
  {return fabs(x) < 1e-5;}

inline bool positive (double x) 
  {return x > 0;}

inline bool negative (double x) 
  {return x < 0;}
