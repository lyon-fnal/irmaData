/*Art analyzer module to produce pileup corrected ratio histograms.
@author  Sudeshna Ganguly
@date    10/21/2019*/

#include "james_module.hh"
bool DEBUG = false;

// now the ClusterHit and PileupHit implementation
void ClusterHit::setPileupHit(PileupHit *pHit)
{
  if (pileupHit != 0)
    if (DEBUG)
      mf::LogVerbatim("RatioEast") << "hit already has a pileup hit set!";
  pileupHit = pHit;
}

void PileupHit::init()
{
  NSources = 0;
  //time = 0.;
  //energy = 0.;
  //recomputed = false;
}
void PileupHit::addSource(ClusterHit *hit)
{
  // first check if it's already in the sources
  for (ClusterHit **src = sources; src < sources + NSources; src++)
  {
    if (*src == hit)
    {
      mf::LogVerbatim("RatioEast")
          << "skipping a pileup positron source that we already added...";
      return;
    }
  }
  // if we got this far then all is good
  //if (DEBUG) mf::LogVerbatim("RatioEast") << "Adding hit to source hits...";
  sources[NSources++] = hit;
  //recomputed = false;
}

/*void PileupHit::recompute() {
  if (!recomputed) {
    time = 0.;
    energy = 0.;
    for (ClusterHit ** src=sources; src<sources+NSources; src++) {
      ClusterHit & srcref = **src;
      time += srcref.getRandomizedAnalysisTime()*srcref.getAnalysisEnergy();
      // TODO: t
      energy += srcref.getAnalysisEnergy();
    }
    time /= energy;
    recomputed = true;
  }
}
double PileupHit::getTime() {
  if (!recomputed) recompute();
  return time;
}
double PileupHit::getEnergy() {
  if (!recomputed) recompute();
  return energy;
}*/

double PileupHit::getAvgAnalysisTime()
{
  double time = 0.;
  double energy = 0.;
  for (ClusterHit **src = sources; src < sources + NSources; src++)
  {
    ClusterHit &srcref = **src;
    time += srcref.getAnalysisTime() * srcref.getAnalysisEnergy();
    energy += srcref.getAnalysisEnergy();
  }
  time /= energy;
  return time;
}
double PileupHit::getAvgRandomizedAnalysisTime()
{
  double time = 0.;
  double energy = 0.;
  for (ClusterHit **src = sources; src < sources + NSources; src++)
  {
    ClusterHit &srcref = **src;
    time += srcref.getRandomizedAnalysisTime() * srcref.getAnalysisEnergy();
    energy += srcref.getAnalysisEnergy();
  }
  time /= energy;
  return time;
}
double PileupHit::getEnergySum()
{
  double energy = 0.;
  for (ClusterHit **src = sources; src < sources + NSources; src++)
  {
    energy += (**src).getAnalysisEnergy();
  }
  return energy;
}
double PileupHit::getAvgEnergy()
{
  return getEnergySum() / ((double)NSources);
}
double PileupHit::getQuarteringConstant()
{
  // the definition of this may be subject to change!
  return sources[0]->quarteringConstant; // <== just pick the first cluster hit
}

ClusterHit allHits[2048];
unsigned int NAllHits = 0;
PileupHit allPileupHits[2048];
unsigned int NPileupHits = 0;

//class ClusterHit;
//class PileupHit;
class ratioeast;

//----------------------------

//----------------------------
///////////////////////////
class ratioeast : public art::EDAnalyzer
{

public:
  explicit ratioeast(fhicl::ParameterSet const &p);

  ratioeast(ratioeast const &) = delete;
  ratioeast(ratioeast &&) = delete;
  ratioeast &operator=(ratioeast const &) = delete;
  ratioeast &operator=(ratioeast &&) = delete;

  // Required functions.
  void analyze(art::Event const &e) override;
  //void endjob() override;

private:
  art::ServiceHandle<art::TFileService> tfs_;
  std::string readModule_, readInstance_;
  bool verbose_, debug_;
  const double Emin_;
  const double tMin_, tMax_;
  const double tBinWidth_;
  const double gapT_;
  const double DT_;
  const double DX_;
  const double DY_;
  const int randSeed1_;
  const int randSeed2_;
  //TRandom3 random1;
  //TRandom3 random2;

  const double maxTimeBin;  // = 700000.;
  const int nBins;          // = int(maxTimeBin/tBinWidth_);
  const double histMaxTime; // = nBins*tBinWidth_/1000

  double E1;
  double E2;
  double T1_up;
  double T1_um;
  double T1;
  double T2;
  double Tn, En;

  //TTree *t_;

  /* double clusterEnergy_;
  double clusterTime_;
  unsigned int caloNum_;
  unsigned int eventNum_;
  unsigned int bunchNum_;
  unsigned int midasSerialNum_;
  unsigned int subRunNum_;
  unsigned int runNum_;
  */
  std::string datasetName_reconEast_;

  //Adding histograms

  TH1F *TIME_FULL;
  TH1F *ENERGY_FULL;

  TH1F *TIME_D;
  TH1F *TIME_S1;
  TH1F *TIME_S2;

  TH1F *ENERGY_D;
  TH1F *ENERGY_S1;
  TH1F *ENERGY_S2;

  /*  TH2F *XY;
  TH1F *X;
  TH1F *Y;
  TH1F *r;
  TH1F *deltax[11];

  TH1F *deltat;
  TH1F *deltat_calo[NCalos];
  TH2F *deltaxvsdeltatallcalo;
  TH2F *deltayvsdeltatallcalo;


  TH2F *deltaxvsdeltat[NCalos];
  TH2F *deltayvsdeltat[NCalos];
  */

  TH1F *ENERGY_FULL_up;
  TH1F *ENERGY_FULL_um;
  TH1F *ENERGY_FULL_vp;
  TH1F *ENERGY_FULL_vm;

  TH1F *ENERGY_D_up;
  TH1F *ENERGY_D_um;
  TH1F *ENERGY_D_vp;
  TH1F *ENERGY_D_vm;

  TH1F *ENERGY_S1_up;
  TH1F *ENERGY_S1_um;
  TH1F *ENERGY_S1_vp;
  TH1F *ENERGY_S1_vm;

  TH1F *ENERGY_S2_up;
  TH1F *ENERGY_S2_um;
  TH1F *ENERGY_S2_vp;
  TH1F *ENERGY_S2_vm;

  TH1F *caloTimes_Up_all;
  TH1F *caloTimes_Um_all;
  TH1F *caloTimes_Vp_all;
  TH1F *caloTimes_Vm_all;

  TH1F *DTimes_Up_all;
  TH1F *DTimes_Um_all;
  TH1F *DTimes_Vp_all;
  TH1F *DTimes_Vm_all;

  TH1F *S1Times_Up_all;
  TH1F *S1Times_Um_all;
  TH1F *S1Times_Vp_all;
  TH1F *S1Times_Vm_all;

  TH1F *S2Times_Up_all;
  TH1F *S2Times_Um_all;
  TH1F *S2Times_Vp_all;
  TH1F *S2Times_Vm_all;

  TH1F *caloEnergies_All[NCalos];
  TH1F *caloTimes_All[NCalos];

  TH1F *energy_d_All[NCalos];
  TH1F *energy_s1_All[NCalos];
  TH1F *energy_s2_All[NCalos];

  TH1F *DTimes_calo_All[NCalos];
  TH1F *S1Times_calo_All[NCalos];
  TH1F *S2Times_calo_All[NCalos];

  TH1F *caloEnergies_Up[NCalos];
  TH1F *caloEnergies_Um[NCalos];
  TH1F *caloEnergies_Vp[NCalos];
  TH1F *caloEnergies_Vm[NCalos];

  TH1F *caloTimes_Up[NCalos];
  TH1F *caloTimes_Um[NCalos];
  TH1F *caloTimes_Vp[NCalos];
  TH1F *caloTimes_Vm[NCalos];

  TH1F *energy_d_Up[NCalos];
  TH1F *energy_d_Um[NCalos];
  TH1F *energy_d_Vp[NCalos];
  TH1F *energy_d_Vm[NCalos];

  TH1F *energy_s1_Up[NCalos];
  TH1F *energy_s1_Um[NCalos];
  TH1F *energy_s1_Vp[NCalos];
  TH1F *energy_s1_Vm[NCalos];

  TH1F *energy_s2_Up[NCalos];
  TH1F *energy_s2_Um[NCalos];
  TH1F *energy_s2_Vp[NCalos];
  TH1F *energy_s2_Vm[NCalos];

  TH1F *DTimes_calo_Up[NCalos];
  TH1F *DTimes_calo_Um[NCalos];
  TH1F *DTimes_calo_Vp[NCalos];
  TH1F *DTimes_calo_Vm[NCalos];

  TH1F *S1Times_calo_Up[NCalos];
  TH1F *S1Times_calo_Um[NCalos];
  TH1F *S1Times_calo_Vp[NCalos];
  TH1F *S1Times_calo_Vm[NCalos];

  TH1F *S2Times_calo_Up[NCalos];
  TH1F *S2Times_calo_Um[NCalos];
  TH1F *S2Times_calo_Vp[NCalos];
  TH1F *S2Times_calo_Vm[NCalos];

  void pileup_hist_up(double T1_up, double T2, double E1, double E2, TH1F *ENERGY_D_up, TH1F *ENERGY_S1_up, TH1F *ENERGY_S2_up, TH1F *DTimes_Up_all, TH1F *S1Times_Up_all, TH1F *S2Times_Up_all, TH1F **energy_d_Up, TH1F **energy_s1_Up, TH1F **energy_s2_Up, TH1F **DTimes_calo_Up, TH1F **S1Times_calo_Up, TH1F **S2Times_calo_Up);
  void pileup_hist_um(double T1_um, double T2, double E1, double E2, TH1F *ENERGY_D_um, TH1F *ENERGY_S1_um, TH1F *ENERGY_S2_um, TH1F *DTimes_Um_all, TH1F *S1Times_Um_all, TH1F *S2Times_Um_all, TH1F **energy_d_Um, TH1F **energy_s1_Um, TH1F **energy_s2_Um, TH1F **DTimes_calo_Um, TH1F **S1Times_calo_Um, TH1F **S2Times_calo_Um);
  void pileup_hist_vp(double T1, double T2, double E1, double E2, TH1F *ENERGY_D_vp, TH1F *ENERGY_S1_vp, TH1F *ENERGY_S2_vp, TH1F *DTimes_Vp_all, TH1F *S1Times_Vp_all, TH1F *S2Times_Vp_all, TH1F **energy_d_Vp, TH1F **energy_s1_Vp, TH1F **energy_s2_Vp, TH1F **DTimes_calo_Vp, TH1F **S1Times_calo_Vp, TH1F **S2Times_calo_Vp);
  void pileup_hist_vm(double T1, double T2, double E1, double E2, TH1F *ENERGY_D_vm, TH1F *ENERGY_S1_vm, TH1F *ENERGY_S2_vm, TH1F *DTimes_Vm_all, TH1F *S1Times_Vm_all, TH1F *S2Times_Vm_all, TH1F **energy_d_Vm, TH1F **energy_s1_Vm, TH1F **energy_s2_Vm, TH1F **DTimes_calo_Vm, TH1F **S1Times_calo_Vm, TH1F **S2Times_calo_Vm);
};

ratioeast::ratioeast(fhicl::ParameterSet const &p)

    :

      EDAnalyzer(p),
      readModule_(p.get<std::string>("readModule", "energyPartition")),
      readInstance_(p.get<std::string>("readInstance", "partition")),
      verbose_(p.get<bool>("verbose", false)),
      debug_(p.get<bool>("debug", false)),
      Emin_(p.get<double>("Emin", 1700.)),
      tMin_(p.get<double>("tmin", 18.)),
      tMax_(p.get<double>("tmax", 700.)),
      tBinWidth_(p.get<double>("tBinWidth", 0.14919)),
      gapT_(p.get<double>("gapT", 0.020)),
      DT_(p.get<double>("DeadT", 0.001)),
      DX_(p.get<double>("DeadX", 3.)),
      DY_(p.get<double>("DeadY", 3.)),
      randSeed1_(p.get<int>("randSeed1", 0)),
      randSeed2_(p.get<int>("randSeed2", 0)),
      //random1(randSeed1_),
      //random2(randSeed2_),
      maxTimeBin(p.get<double>("maxTimeBin", 700.)),
      nBins(int(maxTimeBin / tBinWidth_)),
      histMaxTime(nBins * tBinWidth_),
      E1(0),
      E2(0),
      T1_up(0),
      T1_um(0),
      T1(0),
      T2(0),
      //clusterEnergy_(0),
      //clusterTime_(0),
      //caloNum_(0),
      //eventNum_(0),
      //bunchNum_(-1),
      //midasSerialNum_(0),
      //subRunNum_(0),
      //runNum_(0),
      datasetName_reconEast_(p.get<std::string>("datasetName_reconEast", "60h"))
{
  if (debug_)
    verbose_ = true;
  //if (debug_) {
  //  LogVerbatim("RatioEast") << "[LogVerbatim]: ratioeast() ";
  //  LogInfo("RatioEast") << "[LogInfo]: ratioeast() ";
  //  if (debug_) LogVerbatim("RatioEast") << "[LogDebug]: ratioeast() ";
  //  cout << "cout...\n";
  //  cerr << "cerr...\n";
  //}
  DEBUG = debug_; // global setting for other stuff

  // global vars
  tBinWidth = tBinWidth_;
  random1.SetSeed(randSeed1_);
  random2.SetSeed(randSeed2_);

  art::ServiceHandle<art::TFileService> tfs;
  /*  t_ = tfs->make<TTree>("testTree", "testTree");
  t_->Branch("energy", &clusterEnergy_, "energy/D");
  t_->Branch("time", &clusterTime_, "time/D");
  t_->Branch("caloNum", &caloNum_, "caloNum/i");
  t_->Branch("eventNum", &eventNum_, "eventNum/i");
  t_->Branch("bunchNum", &bunchNum_, "bunchNum/i");
  t_->Branch("midasSerialNum", &midasSerialNum_, "midasSerialNum/i");
  t_->Branch("subRunNum", &subRunNum_, "subRunNum/i");
  t_->Branch("runNum", &runNum_, "runNum/i");
  */

  TIME_FULL = tfs->make<TH1F>("TIME_FULL", "TIME_FULL", nBins, 0, histMaxTime);
  ENERGY_FULL = tfs->make<TH1F>("ENERGY_FULL", "ENERGY ON ALL CALO; energy [MeV]", 500, 0, 10000);
  TIME_D = tfs->make<TH1F>("TIME_D", "TIME_D ON ALL CALO;time [#mus];counts/149ns", nBins, 0, histMaxTime);
  TIME_S1 = tfs->make<TH1F>("TIME_S1", "TIME_S1 ON ALL CALO;time [#mus];counts/149ns", nBins, 0, histMaxTime);
  TIME_S2 = tfs->make<TH1F>("TIME_S2", "TIME_S2 ON ALL CALO;time [#mus];counts/149ns", nBins, 0, histMaxTime);

  ENERGY_D = tfs->make<TH1F>("ENERGY_D", "ENERGY_D ON ALL CALO; energy [MeV]", 500, 0, 10000);
  ENERGY_S1 = tfs->make<TH1F>("ENERGY_S1", "ENERGY_S1 ON ALL CALO; energy [MeV]", 500, 0, 10000);
  ENERGY_S2 = tfs->make<TH1F>("ENERGY_S2", "ENERGY_S2 ON ALL CALO; energy [MeV]", 500, 0, 10000);

  /*  XY = tfs->make<TH2F>("XY", "XY;Horizontal position [c.w.];Vertical position [c.w.]",900,0.5,9.5,600,0.5,6.5);
  X = tfs->make<TH1F>("X", "X;Horizontal position [c.w.]", 900, 0.5, 9.5);
  Y = tfs->make<TH1F>("Y", "Y;Vertical position [c.w.]", 600, 0.5, 6.5);
  r = tfs->make<TH1F>("r", "r; sqrt(x^2+y^2)[c.w]", 1100, 0.5, 11.5);

  for(unsigned int k=0;k<11;k++){
    deltax[k] = tfs->make<TH1F>(Form("deltax_%d",k),Form("#DeltaX_%d;Horizontal position [c.w.]",k), 900, 0.5, 9.5);
  }
  deltat = tfs->make<TH1F>("DeltaT","T_{E2}-T_{E1};time [ns];",1000,0,40);
  deltaxvsdeltatallcalo = tfs->make<TH2F>("deltaxvsdeltat_allcalo", "deltax vs deltat; time separation(ns); #delta X [c.w.]",200, 0, 5, 900, 0.5, 9.5);
  deltayvsdeltatallcalo = tfs->make<TH2F>("deltayvsdeltat_allcalo", "deltay vs deltat; time separation(ns); #delta Y [c.w.]",200, 0, 5, 600, 0.5, 6.5);

  */

  ENERGY_FULL_up = tfs->make<TH1F>("ENERGY_FULL_up", "ENERGY_D ON ALL CALO; energy [MeV]", 500, 0, 10000);
  ENERGY_FULL_um = tfs->make<TH1F>("ENERGY_FULL_um", "ENERGY_D ON ALL CALO; energy [MeV]", 500, 0, 10000);
  ENERGY_FULL_vp = tfs->make<TH1F>("ENERGY_FULL_vp", "ENERGY_D ON ALL CALO; energy [MeV]", 500, 0, 10000);
  ENERGY_FULL_vm = tfs->make<TH1F>("ENERGY_FULL_vm", "ENERGY_D ON ALL CALO; energy [MeV]", 500, 0, 10000);

  ENERGY_D_up = tfs->make<TH1F>("ENERGY_D_up", "ENERGY_D_up ON ALL CALO; energy [MeV]", 500, 0, 10000);
  ENERGY_D_um = tfs->make<TH1F>("ENERGY_D_um", "ENERGY_D_um ON ALL CALO; energy [MeV]", 500, 0, 10000);
  ENERGY_D_vp = tfs->make<TH1F>("ENERGY_D_vp", "ENERGY_D_vp ON ALL CALO; energy [MeV]", 500, 0, 10000);
  ENERGY_D_vm = tfs->make<TH1F>("ENERGY_D_vm", "ENERGY_D_vm ON ALL CALO; energy [MeV]", 500, 0, 10000);

  ENERGY_S1_up = tfs->make<TH1F>("ENERGY_S1_up", "ENERGY_S1_up ON ALL CALO; energy [MeV]", 500, 0, 10000);
  ENERGY_S1_um = tfs->make<TH1F>("ENERGY_S1_um", "ENERGY_S1_um ON ALL CALO; energy [MeV]", 500, 0, 10000);
  ENERGY_S1_vp = tfs->make<TH1F>("ENERGY_S1_vp", "ENERGY_S1_vp ON ALL CALO; energy [MeV]", 500, 0, 10000);
  ENERGY_S1_vm = tfs->make<TH1F>("ENERGY_S1_vm", "ENERGY_S1_vm ON ALL CALO; energy [MeV]", 500, 0, 10000);

  ENERGY_S2_up = tfs->make<TH1F>("ENERGY_S2_up", "ENERGY_S2_up ON ALL CALO; energy [MeV]", 500, 0, 10000);
  ENERGY_S2_um = tfs->make<TH1F>("ENERGY_S2_um", "ENERGY_S2_um ON ALL CALO; energy [MeV]", 500, 0, 10000);
  ENERGY_S2_vp = tfs->make<TH1F>("ENERGY_S2_vp", "ENERGY_S2_vp ON ALL CALO; energy [MeV]", 500, 0, 10000);
  ENERGY_S2_vm = tfs->make<TH1F>("ENERGY_S2_vm", "ENERGY_S2_vm ON ALL CALO; energy [MeV]", 500, 0, 10000);

  caloTimes_Up_all = tfs->make<TH1F>("caloTimes_Up_all", "caloTimes_Up_all", nBins, 0, histMaxTime);
  caloTimes_Um_all = tfs->make<TH1F>("caloTimes_Um_all", "caloTimes_Um_all", nBins, 0, histMaxTime);
  caloTimes_Vp_all = tfs->make<TH1F>("caloTimes_Vp_all", "caloTimes_Vp_all", nBins, 0, histMaxTime);
  caloTimes_Vm_all = tfs->make<TH1F>("caloTimes_Vm_all", "caloTimes_Vm_all", nBins, 0, histMaxTime);

  DTimes_Up_all = tfs->make<TH1F>("DTimes_Up_all", "DTimes_Up_all", nBins, 0, histMaxTime);
  DTimes_Um_all = tfs->make<TH1F>("DTimes_Um_all", "DTimes_Um_all", nBins, 0, histMaxTime);
  DTimes_Vp_all = tfs->make<TH1F>("DTimes_Vp_all", "DTimes_Vp_all", nBins, 0, histMaxTime);
  DTimes_Vm_all = tfs->make<TH1F>("DTimes_Vm_all", "DTimes_Vm_all", nBins, 0, histMaxTime);

  S1Times_Up_all = tfs->make<TH1F>("S1Times_Up_all", "S1Times_Up_all", nBins, 0, histMaxTime);
  S1Times_Um_all = tfs->make<TH1F>("S1Times_Um_all", "S1Times_Um_all", nBins, 0, histMaxTime);
  S1Times_Vp_all = tfs->make<TH1F>("S1Times_Vp_all", "S1Times_Vp_all", nBins, 0, histMaxTime);
  S1Times_Vm_all = tfs->make<TH1F>("S1Times_Vm_all", "S1Times_Vm_all", nBins, 0, histMaxTime);

  S2Times_Up_all = tfs->make<TH1F>("S2Times_Up_all", "S2Times_Up_all", nBins, 0, histMaxTime);
  S2Times_Um_all = tfs->make<TH1F>("S2Times_Um_all", "S2Times_Um_all", nBins, 0, histMaxTime);
  S2Times_Vp_all = tfs->make<TH1F>("S2Times_Vp_all", "S2Times_Vp_all", nBins, 0, histMaxTime);
  S2Times_Vm_all = tfs->make<TH1F>("S2Times_Vm_all", "S2Times_Vm_all", nBins, 0, histMaxTime);

  for (unsigned int iCalo = 0; iCalo < NCalos; iCalo++)
  {
    caloEnergies_All[iCalo] = tfs->make<TH1F>(Form("EnergyOriginal_All%d", iCalo), Form("Energy_Spectrum_calo_All_%d;energy [MeV]", iCalo), 500, 0, 10000);
    caloTimes_All[iCalo] = tfs->make<TH1F>(Form("TimeOriginal_All_calo%d", iCalo), Form("Time_calo_All%d;time [#mus];counts/145.15ns", iCalo), nBins, 0, histMaxTime);

    caloEnergies_Up[iCalo] = tfs->make<TH1F>(Form("EnergyOriginal_Up%d", iCalo), Form("Energy_Spectrum_calo_Up_%d;energy [MeV]", iCalo), 500, 0, 10000);
    caloEnergies_Um[iCalo] = tfs->make<TH1F>(Form("EnergyOriginal_Um%d", iCalo), Form("Energy_Spectrum_calo_Um_%d;energy [MeV]", iCalo), 500, 0, 10000);
    caloEnergies_Vp[iCalo] = tfs->make<TH1F>(Form("EnergyOriginal_Vp%d", iCalo), Form("Energy_Spectrum_calo_Vp_%d;energy [MeV]", iCalo), 500, 0, 10000);
    caloEnergies_Vm[iCalo] = tfs->make<TH1F>(Form("EnergyOriginal_Vm%d", iCalo), Form("Energy_Spectrum_calo_Vm_%d;energy [MeV]", iCalo), 500, 0, 10000);

    caloTimes_Up[iCalo] = tfs->make<TH1F>(Form("TimeOriginal_Up_calo%d", iCalo), Form("Time_calo_%d;time [#mus];counts/147ns", iCalo), nBins, 0, histMaxTime);
    caloTimes_Um[iCalo] = tfs->make<TH1F>(Form("TimeOriginal_Um_calo%d", iCalo), Form("Time_Um_calo_%d;time [#mus];counts/147ns", iCalo), nBins, 0, histMaxTime);
    caloTimes_Vp[iCalo] = tfs->make<TH1F>(Form("TimeOriginal_Vp_calo%d", iCalo), Form("Time_Vp_calo_%d;time [#mus];counts/147ns", iCalo), nBins, 0, histMaxTime);
    caloTimes_Vm[iCalo] = tfs->make<TH1F>(Form("TimeOriginal_Vm_calo%d", iCalo), Form("Time_Vm_calo_%d;time [#mus];counts/147ns", iCalo), nBins, 0, histMaxTime);

    //deltaxvsdeltat[iCalo] = tfs->make<TH2F>(Form("deltaxvsdeltat_%d",iCalo),Form("deltax vs deltat_%d; time separation(ns); #delta X [c.w.]",iCalo),200, 0, 5, 900, 0.5, 9.5);
    //deltayvsdeltat[iCalo] = tfs->make<TH2F>(Form("deltayvsdeltat_%d",iCalo),Form("deltay vs deltat_%d; time separation(ns); #delta Y [c.w.]",iCalo),200, 0, 5, 600, 0.5, 6.5);

    //deltat_calo[iCalo] = tfs->make<TH1F>(Form("DeltaT_%d",iCalo),Form("T_{E2}-T_{E1} calo %d;time [ns];",iCalo),1000,0,40);

    energy_d_All[iCalo] = tfs->make<TH1F>(Form("energy_D_All_%d", iCalo), Form("energy_Double_calo_All_%d;energy [MeV]", iCalo), 500, 0, 10000);
    energy_s1_All[iCalo] = tfs->make<TH1F>(Form("energy_S1_All_%d", iCalo), Form("energy_Single_1_calo_All_%d;energy [MeV]", iCalo), 500, 0, 10000);
    energy_s2_All[iCalo] = tfs->make<TH1F>(Form("energy_S2_All_%d", iCalo), Form("energy_Single_2_calo_All_%d;energy [MeV]", iCalo), 500, 0, 10000);

    DTimes_calo_All[iCalo] = tfs->make<TH1F>(Form("time_D_All_calo_%d", iCalo), Form("double_pileup_time_All_calo_%d;time [#mus];counts/147ns", iCalo), nBins, 0, histMaxTime);
    S1Times_calo_All[iCalo] = tfs->make<TH1F>(Form("time_S1_All_calo_%d", iCalo), Form("single_1_pileup_time_All_calo_%d;time [#mus];counts/148ns", iCalo), nBins, 0, histMaxTime);
    S2Times_calo_All[iCalo] = tfs->make<TH1F>(Form("time_S2_All_calo_%d", iCalo), Form("single_2_pileup_time_All_calo_%d;time [#mus];counts/148ns", iCalo), nBins, 0, histMaxTime);

    energy_d_Up[iCalo] = tfs->make<TH1F>(Form("energy_D_Up_%d", iCalo), Form("energy_Double_calo_Up_%d;energy [MeV]", iCalo), 500, 0, 10000);
    energy_s1_Up[iCalo] = tfs->make<TH1F>(Form("energy_S1_Up_%d", iCalo), Form("energy_Single_1_calo_Up_%d;energy [MeV]", iCalo), 500, 0, 10000);
    energy_s2_Up[iCalo] = tfs->make<TH1F>(Form("energy_S2_Up_%d", iCalo), Form("energy_Single_2_calo_Up_%d;energy [MeV]", iCalo), 500, 0, 10000);

    energy_d_Um[iCalo] = tfs->make<TH1F>(Form("energy_D_Um_%d", iCalo), Form("energy_Double_calo_Um_%d;energy [MeV]", iCalo), 500, 0, 10000);
    energy_s1_Um[iCalo] = tfs->make<TH1F>(Form("energy_S1_Um_%d", iCalo), Form("energy_Single_1_calo_Um_%d;energy [MeV]", iCalo), 500, 0, 10000);
    energy_s2_Um[iCalo] = tfs->make<TH1F>(Form("energy_S2_Um_%d", iCalo), Form("energy_Single_2_calo_Um_%d;energy [MeV]", iCalo), 500, 0, 10000);

    energy_d_Vp[iCalo] = tfs->make<TH1F>(Form("energy_D_Vp_%d", iCalo), Form("energy_Double_calo_Vp_%d;energy [MeV]", iCalo), 500, 0, 10000);
    energy_s1_Vp[iCalo] = tfs->make<TH1F>(Form("energy_S1_Vp_%d", iCalo), Form("energy_Single_1_calo_Vp_%d;energy [MeV]", iCalo), 500, 0, 10000);
    energy_s2_Vp[iCalo] = tfs->make<TH1F>(Form("energy_S2_Vp_%d", iCalo), Form("energy_Single_2_calo_Vp_%d;energy [MeV]", iCalo), 500, 0, 10000);

    energy_d_Vm[iCalo] = tfs->make<TH1F>(Form("energy_D_Vm_%d", iCalo), Form("energy_Double_calo_%d;energy [MeV]", iCalo), 500, 0, 10000);
    energy_s1_Vm[iCalo] = tfs->make<TH1F>(Form("energy_S1_Vm_%d", iCalo), Form("energy_Single_1_calo_%d;energy [MeV]", iCalo), 500, 0, 10000);
    energy_s2_Vm[iCalo] = tfs->make<TH1F>(Form("energy_S2_Vm_%d", iCalo), Form("energy_Single_2_calo_%d;energy [MeV]", iCalo), 500, 0, 10000);

    DTimes_calo_Up[iCalo] = tfs->make<TH1F>(Form("time_D_Up_calo_%d", iCalo), Form("double_pileup_time_Up_calo_%d;time [#mus];counts/147ns", iCalo), nBins, 0, histMaxTime);
    DTimes_calo_Um[iCalo] = tfs->make<TH1F>(Form("time_D_Um_calo_%d", iCalo), Form("double pileup_time_Um_calo_%d;time [#mus];counts/147ns", iCalo), nBins, 0, histMaxTime);
    DTimes_calo_Vp[iCalo] = tfs->make<TH1F>(Form("time_D_Vp_calo_%d", iCalo), Form("double pileup time_Vp_calo_%d;time [#mus];counts/147ns", iCalo), nBins, 0, histMaxTime);
    DTimes_calo_Vm[iCalo] = tfs->make<TH1F>(Form("time_D_Vm_calo_%d", iCalo), Form("double pileup time_Vm_calo_%d;time [#mus];counts/148ns", iCalo), nBins, 0, histMaxTime);

    S1Times_calo_Up[iCalo] = tfs->make<TH1F>(Form("time_S1_Up_calo_%d", iCalo), Form("single_1_pileup_time_Up_calo_%d;time [#mus];counts/148ns", iCalo), nBins, 0, histMaxTime);
    S1Times_calo_Um[iCalo] = tfs->make<TH1F>(Form("time_S1_Um_calo_%d", iCalo), Form("single_1_pileup_time_Um_calo_%d;time [#mus];counts/148ns", iCalo), nBins, 0, histMaxTime);
    S1Times_calo_Vp[iCalo] = tfs->make<TH1F>(Form("time_S1_Vp_calo_%d", iCalo), Form("single_1_pileup_time_Vp_calo_%d;time [#mus];counts/148ns", iCalo), nBins, 0, histMaxTime);
    S1Times_calo_Vm[iCalo] = tfs->make<TH1F>(Form("time_S1_Vm_calo_%d", iCalo), Form("single_1_pileup_time_Vm_calo_%d;time [#mus];counts/148ns", iCalo), nBins, 0, histMaxTime);

    S2Times_calo_Up[iCalo] = tfs->make<TH1F>(Form("time_S2_Up_calo_%d", iCalo), Form("single_2_pileup_time_Up_calo_%d;time [#mus];counts/148ns", iCalo), nBins, 0, histMaxTime);
    S2Times_calo_Um[iCalo] = tfs->make<TH1F>(Form("time_S2_Um_calo_%d", iCalo), Form("single_2_pileup_time_Um_calo_%d;time [#mus];counts/148ns", iCalo), nBins, 0, histMaxTime);
    S2Times_calo_Vp[iCalo] = tfs->make<TH1F>(Form("time_S2_Vp_calo_%d", iCalo), Form("single_2_pileup_time_Vp_calo_%d;time [#mus];counts/148ns", iCalo), nBins, 0, histMaxTime);
    S2Times_calo_Vm[iCalo] = tfs->make<TH1F>(Form("time_S2_Vm_calo_%d", iCalo), Form("single_2_pileup_time_Vm_calo_%d;time [#mus];counts/148ns", iCalo), nBins, 0, histMaxTime);
  }
}

void ratioeast::pileup_hist_up(double T1_up, double T2, double E1, double E2, TH1F *ENERGY_D_up, TH1F *ENERGY_S1_up, TH1F *ENERGY_S2_up, TH1F *DTimes_Up_all, TH1F *S1Times_Up_all, TH1F *S2Times_Up_all, TH1F **energy_d_Up, TH1F **energy_s1_Up, TH1F **energy_s2_Up, TH1F **DTimes_calo_Up, TH1F **S1Times_calo_Up, TH1F **S2Times_calo_Up)
{

  if (E1 < Emin_ && E2 < Emin_)
  {
    if (E1 + E2 > Emin_)
    {
      ENERGY_D_up->Fill(E1 + E2);
      DTimes_Up_all->Fill(T1_up);

      for (unsigned int iCalo = 0; iCalo < NCalos; iCalo++)
      {
        energy_d_Up[iCalo]->Fill(E1 + E2);
        DTimes_calo_Up[iCalo]->Fill(T1_up);
      }
      //delete energy_d_Up;
      //delete DTimes_calo_Up;
    }
  }

  if (E1 > Emin_ && E2 < Emin_)
  {
    ENERGY_D_up->Fill(E1 + E2);
    ENERGY_S1_up->Fill(E1);
    DTimes_Up_all->Fill(T1_up);
    S1Times_Up_all->Fill(T1_up);

    for (unsigned int iCalo = 0; iCalo < NCalos; iCalo++)
    {
      energy_d_Up[iCalo]->Fill(E1 + E2);
      energy_s1_Up[iCalo]->Fill(E1);
      DTimes_calo_Up[iCalo]->Fill(T1_up);
      S1Times_calo_Up[iCalo]->Fill(T1_up);
    }
    //delete energy_d_Up;
    //delete energy_s1_Up;
    //delete DTimes_calo_Up;
    //delete S1Times_calo_Up;
  }

  if (E1 < Emin_ && E2 > Emin_)
  {
    ENERGY_D_up->Fill(E1 + E2);
    ENERGY_S2_up->Fill(E2);
    DTimes_Up_all->Fill(T1_up);
    S2Times_Up_all->Fill(T1_up);

    for (unsigned int iCalo = 0; iCalo < NCalos; iCalo++)
    {
      energy_d_Up[iCalo]->Fill(E1 + E2);
      energy_s2_Up[iCalo]->Fill(E2);
      DTimes_calo_Up[iCalo]->Fill(T1_up);
      S2Times_calo_Up[iCalo]->Fill(T1_up);
    }
    //delete energy_d_Up;
    //delete energy_s2_Up;
    //delete DTimes_calo_Up;
    //delete S2Times_calo_Up;
  }

  if (E1 > Emin_ && E2 > Emin_)
  {
    ENERGY_D_up->Fill(E1 + E2);
    ENERGY_S1_up->Fill(E1);
    ENERGY_S2_up->Fill(E2);
    DTimes_Up_all->Fill(T1_up);
    S1Times_Up_all->Fill(T1_up);
    S2Times_Up_all->Fill(T1_up);

    for (unsigned int iCalo = 0; iCalo < NCalos; iCalo++)
    {
      energy_d_Up[iCalo]->Fill(E1 + E2);
      energy_s1_Up[iCalo]->Fill(E1);
      energy_s2_Up[iCalo]->Fill(E2);
      DTimes_calo_Up[iCalo]->Fill(T1_up);
      S1Times_calo_Up[iCalo]->Fill(T1_up);
      S2Times_calo_Up[iCalo]->Fill(T1_up);
    }
    //delete energy_d_Up;
    //delete energy_s1_Up;
    //delete energy_s2_Up;
    //delete DTimes_calo_Up;
    //delete S1Times_calo_Up;
    //delete S2Times_calo_Up;
  }
}

void ratioeast::pileup_hist_um(double T1_um, double T2, double E1, double E2, TH1F *ENERGY_D_um, TH1F *ENERGY_S1_um, TH1F *ENERGY_S2_um, TH1F *DTimes_Um_all, TH1F *S1Times_Um_all, TH1F *S2Times_Um_all, TH1F **energy_d_Um, TH1F **energy_s1_Um, TH1F **energy_s2_Um, TH1F **DTimes_calo_Um, TH1F **S1Times_calo_Um, TH1F **S2Times_calo_Um)
{

  if (E1 < Emin_ && E2 < Emin_)
  {
    if (E1 + E2 > Emin_)
    {
      ENERGY_D_um->Fill(E1 + E2);
      DTimes_Um_all->Fill(T1_um);
      for (unsigned int iCalo = 0; iCalo < NCalos; iCalo++)
      {
        energy_d_Um[iCalo]->Fill(E1 + E2);
        DTimes_calo_Um[iCalo]->Fill(T1_um);
      }
      //delete energy_d_Um;
      //delete DTimes_calo_Um;
    }
  }

  if (E1 > Emin_ && E2 < Emin_)
  {
    ENERGY_D_um->Fill(E1 + E2);
    ENERGY_S1_um->Fill(E1);
    DTimes_Um_all->Fill(T1_um);
    S1Times_Um_all->Fill(T1_um);

    for (unsigned int iCalo = 0; iCalo < NCalos; iCalo++)
    {
      energy_d_Um[iCalo]->Fill(E1 + E2);
      energy_s1_Um[iCalo]->Fill(E1);
      DTimes_calo_Um[iCalo]->Fill(T1_um);
      S1Times_calo_Um[iCalo]->Fill(T1_um);
    }
    //delete energy_d_Um;
    //delete energy_s1_Um;
    //delete DTimes_calo_Um;
    //delete S1Times_calo_Um;
  }

  if (E1 < Emin_ && E2 > Emin_)
  {
    ENERGY_D_um->Fill(E1 + E2);
    ENERGY_S2_um->Fill(E2);
    DTimes_Um_all->Fill(T1_um);
    S2Times_Um_all->Fill(T1_um);

    for (unsigned int iCalo = 0; iCalo < NCalos; iCalo++)
    {
      energy_d_Um[iCalo]->Fill(E1 + E2);
      energy_s2_Um[iCalo]->Fill(E2);
      DTimes_calo_Um[iCalo]->Fill(T1_um);
      S2Times_calo_Um[iCalo]->Fill(T1_um);
    }
    //delete energy_d_Um;
    //delete energy_s1_Um;
    //delete DTimes_calo_Um;
    //delete S2Times_calo_Um;
  }

  if (E1 > Emin_ && E2 > Emin_)
  {
    ENERGY_D_um->Fill(E1 + E2);
    ENERGY_S1_um->Fill(E1);
    ENERGY_S2_um->Fill(E2);
    DTimes_Um_all->Fill(T1_um);
    S1Times_Um_all->Fill(T1_um);
    S2Times_Um_all->Fill(T1_um);

    for (unsigned int iCalo = 0; iCalo < NCalos; iCalo++)
    {
      energy_d_Um[iCalo]->Fill(E1 + E2);
      energy_s1_Um[iCalo]->Fill(E1);
      energy_s2_Um[iCalo]->Fill(E2);
      DTimes_calo_Um[iCalo]->Fill(T1_um);
      S1Times_calo_Um[iCalo]->Fill(T1_um);
      S2Times_calo_Um[iCalo]->Fill(T1_um);
    }
    //delete energy_d_Um;
    //delete energy_s1_Um;
    //delete energy_s2_Um;
    //delete DTimes_calo_Um;
    //delete S1Times_calo_Um;
    //delete S2Times_calo_Um;
  }
}

void ratioeast::pileup_hist_vp(double T1, double T2, double E1, double E2, TH1F *ENERGY_D_vp, TH1F *ENERGY_S1_vp, TH1F *ENERGY_S2_vp, TH1F *DTimes_Vp_all, TH1F *S1Times_Vp_all, TH1F *S2Times_Vp_all, TH1F **energy_d_Vp, TH1F **energy_s1_Vp, TH1F **energy_s2_Vp, TH1F **DTimes_calo_Vp, TH1F **S1Times_calo_Vp, TH1F **S2Times_calo_Vp)
{

  if (E1 < Emin_ && E2 < Emin_)
  {
    if (E1 + E2 > Emin_)
    {
      ENERGY_D_vp->Fill(E1 + E2);
      DTimes_Vp_all->Fill(T1);
      for (unsigned int iCalo = 0; iCalo < NCalos; iCalo++)
      {
        energy_d_Vp[iCalo]->Fill(E1 + E2);
        DTimes_calo_Vp[iCalo]->Fill(T1);
      }
      //delete energy_d_Vp;
      //delete DTimes_calo_Vp;
    }
  }

  if (E1 > Emin_ && E2 < Emin_)
  {
    ENERGY_D_vp->Fill(E1 + E2);
    ENERGY_S1_vp->Fill(E1);
    DTimes_Vp_all->Fill(T1);
    S1Times_Vp_all->Fill(T1);

    for (unsigned int iCalo = 0; iCalo < NCalos; iCalo++)
    {
      energy_d_Vp[iCalo]->Fill(E1 + E2);
      energy_s1_Vp[iCalo]->Fill(E1);
      DTimes_calo_Vp[iCalo]->Fill(T1);
      S1Times_calo_Vp[iCalo]->Fill(T1);
    }
    //delete energy_d_Vp;
    //delete energy_s1_Vp;
    //delete DTimes_calo_Vp;
    //delete S1Times_calo_Vp;
  }

  if (E1 < Emin_ && E2 > Emin_)
  {
    ENERGY_D_vp->Fill(E1 + E2);
    ENERGY_S2_vp->Fill(E2);
    DTimes_Vp_all->Fill(T1);
    S2Times_Vp_all->Fill(T1);

    for (unsigned int iCalo = 0; iCalo < NCalos; iCalo++)
    {
      energy_d_Vp[iCalo]->Fill(E1 + E2);
      energy_s2_Vp[iCalo]->Fill(E2);
      DTimes_calo_Vp[iCalo]->Fill(T1);
      S2Times_calo_Vp[iCalo]->Fill(T1);
    }
    //delete energy_d_Vp;
    //delete energy_s2_Vp;
    //delete DTimes_calo_Vp;
    //delete S2Times_calo_Vp;
  }

  if (E1 > Emin_ && E2 > Emin_)
  {
    ENERGY_D_vp->Fill(E1 + E2);
    ENERGY_S1_vp->Fill(E1);
    ENERGY_S2_vp->Fill(E2);
    DTimes_Vp_all->Fill(T1);
    S1Times_Vp_all->Fill(T1);
    S2Times_Vp_all->Fill(T1);

    for (unsigned int iCalo = 0; iCalo < NCalos; iCalo++)
    {
      energy_d_Vp[iCalo]->Fill(E1 + E2);
      energy_s1_Vp[iCalo]->Fill(E1);
      energy_s2_Vp[iCalo]->Fill(E2);
      DTimes_calo_Vp[iCalo]->Fill(T1);
      S1Times_calo_Vp[iCalo]->Fill(T1);
      S2Times_calo_Vp[iCalo]->Fill(T1);
    }
    //delete energy_d_Vp;
    //delete energy_s1_Vp;
    //delete energy_s2_Vp;
    //delete DTimes_calo_Vp;
    //delete S1Times_calo_Vp;
    //delete S2Times_calo_Vp;
  }
}

void ratioeast::pileup_hist_vm(double T1, double T2, double E1, double E2, TH1F *ENERGY_D_vm, TH1F *ENERGY_S1_vm, TH1F *ENERGY_S2_vm, TH1F *DTimes_Vm_all, TH1F *S1Times_Vm_all, TH1F *S2Times_Vm_all, TH1F **energy_d_Vm, TH1F **energy_s1_Vm, TH1F **energy_s2_Vm, TH1F **DTimes_calo_Vm, TH1F **S1Times_calo_Vm, TH1F **S2Times_calo_Vm)
{

  if (E1 < Emin_ && E2 < Emin_)
  {
    if (E1 + E2 > Emin_)
    {
      ENERGY_D_vm->Fill(E1 + E2);
      DTimes_Vm_all->Fill(T1);

      for (unsigned int iCalo = 0; iCalo < NCalos; iCalo++)
      {
        energy_d_Vm[iCalo]->Fill(E1 + E2);
        DTimes_calo_Vm[iCalo]->Fill(T1);
      }
      //delete energy_d_Vm;
      //delete DTimes_calo_Vm;
    }
  }

  if (E1 > Emin_ && E2 < Emin_)
  {
    ENERGY_D_vm->Fill(E1 + E2);
    ENERGY_S1_vm->Fill(E1);
    DTimes_Vm_all->Fill(T1);
    S1Times_Vm_all->Fill(T1);

    for (unsigned int iCalo = 0; iCalo < NCalos; iCalo++)
    {
      energy_d_Vm[iCalo]->Fill(E1 + E2);
      energy_s1_Vm[iCalo]->Fill(E1);
      DTimes_calo_Vm[iCalo]->Fill(T1);
      S1Times_calo_Vm[iCalo]->Fill(T1);
    }
    //delete energy_d_Vm;
    //delete energy_s1_Vm;
    //delete DTimes_calo_Vm;
    //delete S1Times_calo_Vm;
  }
  if (E1 < Emin_ && E2 > Emin_)
  {

    ENERGY_D_vm->Fill(E1 + E2);
    ENERGY_S2_vm->Fill(E2);
    DTimes_Vm_all->Fill(T1);
    S2Times_Vm_all->Fill(T1);

    for (unsigned int iCalo = 0; iCalo < NCalos; iCalo++)
    {

      energy_d_Vm[iCalo]->Fill(E1 + E2);
      energy_s2_Vm[iCalo]->Fill(E2);
      DTimes_calo_Vm[iCalo]->Fill(T1);
      S2Times_calo_Vm[iCalo]->Fill(T1);
    }
    //delete energy_d_Vm;
    //delete energy_s2_Vm;
    //delete DTimes_calo_Vm;
    //delete S2Times_calo_Vm;
  }

  if (E1 > Emin_ && E2 > Emin_)
  {
    ENERGY_D_vm->Fill(E1 + E2);
    ENERGY_S1_vm->Fill(E1);
    ENERGY_S2_vm->Fill(E2);
    DTimes_Vm_all->Fill(T1);
    S1Times_Vm_all->Fill(T1);
    S2Times_Vm_all->Fill(T1);

    for (unsigned int iCalo = 0; iCalo < NCalos; iCalo++)
    {
      energy_d_Vm[iCalo]->Fill(E1 + E2);
      energy_s1_Vm[iCalo]->Fill(E1);
      energy_s2_Vm[iCalo]->Fill(E2);
      DTimes_calo_Vm[iCalo]->Fill(T1);
      S1Times_calo_Vm[iCalo]->Fill(T1);
      S2Times_calo_Vm[iCalo]->Fill(T1);
    }
  }
  // delete [] energy_d_Vm;
  //delete energy_s1_Vm;
  //delete energy_s2_Vm;
  //delete DTimes_calo_Vm;
  //delete S1Times_calo_Vm;
  //delete S2Times_calo_Vm;
}

void ratioeast::analyze(art::Event const &e)
{

  if (debug_)
    mf::LogVerbatim("RatioEast")
        << "Enter analyze()..." << e.id();

  // get sequence and trigger indices from fc7
  //const auto &encodercol =
  //*e.getValidHandle<gm2ccc::EncoderFC7ArtRecordCollection>(
  //  {"cccUnpacker", "unpacker"});
  //  bunchNum_ = encodercol.size() ? encodercol.front().sequenceIndex : -1u;

  //runNum_ = e.run();
  //subRunNum_ = e.subRun();
  //eventNum_ = e.event();

  // get recon east data products
  /*const auto &reconEastCollection =
      *e.getValidHandle<GlobalFitArtRecordCollection>(
  {readModule_, readInstance_});

  for(const auto &cluster : reconEastCollection){
    if(cluster.time*1.25/1000 >= 29 && cluster.time*1.25/1000 <= 651 && cluster.energy > 1700 ){
      clusterTime_ = cluster.time*1.25/1000;
      clusterEnergy_ = cluster.energy;
      t_->Fill();
    }
    }*/

  // get recon east data products, using CaloGlobalFitViewArtRecords
  //  to loop over calos first
  art::Handle<CaloGlobalFitViewArtRecordCollection> h;
  if (!e.getByLabel({readModule_, readInstance_}, h))
  {
    if (debug_)
      mf::LogVerbatim("ratioeast_module")
          << "Found NO CaloGlobalFitViewArtRecords!";
    return;
  }
  CaloGlobalFitViewArtRecordCollection caloViews = *h;
  if (debug_)
    mf::LogVerbatim("RatioEast")
        << "Found " << caloViews.size() << " CaloGlobalFitViewArtRecords";

  // ### thus begins the loop over calorimeters ###
  for (auto clustersThisCalo : caloViews)
  {
    unsigned int iCalo = clustersThisCalo.calorimeterIndex - 1;
    if (debug_)
      mf::LogVerbatim("RatioEast")
          << "Loop iteration for calorimeterIndex=" << clustersThisCalo.calorimeterIndex
          << " (iCalo=" << iCalo << ")";

    // ### loop over GlobalFit clusters and put them all in ClusterHites ###
    art::PtrVector<GlobalFitArtRecord> clusters = clustersThisCalo.clusters;
    NAllHits = 0;
    NPileupHits = 0;
    for (auto thisCluster : clusters)
    {
      //if (debug_) mf::LogVerbatim("RatioEast") << "loop over clusters";
      //ClusterHit * aHitPtr = &((*thisCluster)->get());
      //const GlobalFitArtRecord * gfPtr = &(*thisCluster);
      //ClusterHit aHit(gfPtr);
      ClusterHit aHit(thisCluster.get());

      // first fill in the T method histogram (and friends) but with a
      // different set of cuts (on energy at least) than the hits destined
      // for the matrix of all cluster hits
      if (
          hitPassesTimeCuts(&aHit, tMin_, tMax_, g2HalfPeriod + 1.5 * tBinWidth_) && hitPassesEnergyCuts(&aHit, Emin_))
      {
        TIME_FULL->Fill(aHit.getRandomizedAnalysisTime());
        ENERGY_FULL->Fill(aHit.getAnalysisEnergy());
        caloEnergies_All[iCalo]->Fill(aHit.getAnalysisEnergy());
        caloTimes_All[iCalo]->Fill(aHit.getRandomizedAnalysisTime());

        // aHit.quarteringConstant should be [0.,4.)
        // which type is this?
        const int hitType = int(aHit.quarteringConstant);
        // quarter 1 (U+)
        if (hitType == 0)
        {
          ENERGY_FULL_up->Fill(aHit.getAnalysisEnergy());
          caloEnergies_Up[iCalo]->Fill(aHit.getAnalysisEnergy());
          caloTimes_Up_all->Fill(
              aHit.getRandomizedAndShiftedAnalysisTime());
          caloTimes_Up[iCalo]->Fill(
              aHit.getRandomizedAndShiftedAnalysisTime());
        }
        // quarter 2 (U-)
        else if (hitType == 1)
        {
          ENERGY_FULL_um->Fill(aHit.getAnalysisEnergy());
          caloEnergies_Um[iCalo]->Fill(aHit.getAnalysisEnergy());
          caloTimes_Um_all->Fill(
              aHit.getRandomizedAndShiftedAnalysisTime());
          caloTimes_Um[iCalo]->Fill(
              aHit.getRandomizedAndShiftedAnalysisTime());
        }
        // quarter 3 (V1)
        else if (hitType == 2)
        {
          ENERGY_FULL_vp->Fill(aHit.getAnalysisEnergy());
          caloEnergies_Vp[iCalo]->Fill(aHit.getAnalysisEnergy());
          caloTimes_Vp_all->Fill(
              aHit.getRandomizedAndShiftedAnalysisTime());
          caloTimes_Vp[iCalo]->Fill(
              aHit.getRandomizedAndShiftedAnalysisTime());
        }
        //quarter 4 (V2)
        else if (hitType == 3)
        {
          ENERGY_FULL_vm->Fill(aHit.getAnalysisEnergy());
          caloEnergies_Vm[iCalo]->Fill(aHit.getAnalysisEnergy());
          caloTimes_Vm_all->Fill(
              aHit.getRandomizedAndShiftedAnalysisTime());
          caloTimes_Vm[iCalo]->Fill(
              aHit.getRandomizedAndShiftedAnalysisTime());
        }
        // else raise some exception because quarteringConstant=4. exactly?;
      }
      // end filling T method histogram

      // put clusters with E>500MeV in allHits
      if (
          hitPassesTimeCuts(&aHit, tMin_, tMax_, g2HalfPeriod + 1.5 * tBinWidth_)
          // not here: && hitPassesEnergyCuts(&aHit,500.)
      )
        allHits[NAllHits++] = aHit;
    }
    if (debug_)
      mf::LogVerbatim("RatioEast")
          << "Finished tabulating clusters as ClusterHites (NAllHits=" << NAllHits << ")";

    // ### search through allHits and look for shadow coincidences ###
    for (ClusterHit *thisHit = allHits; thisHit < allHits + NAllHits - 1; thisHit++) //TODO NAllHits-1?
    {
      //if (debug_) mf::LogVerbatim("RatioEast") << "loop over allHits";
      ClusterHit *nextHit = thisHit + 1;

      // ### found pileup hit by deadtime coincidence ###
      if (shadowCoincidenceMatch(thisHit, nextHit, DT_))
      { // found one!
        if (debug_)
        {
          mf::LogVerbatim("RatioEast") << "Found a DTC!\n"
                                       << "  thisHit: " << thisHit->summary() << "\n"
                                       << "  nextHit: " << nextHit->summary();
        }
        PileupHit aPileupHit(thisHit, nextHit);    // create pileup hit from pair
        thisHit->setPileupHit(&aPileupHit);        // tell hits to remember their
        nextHit->setPileupHit(&aPileupHit);        //   pileup hits ('children')
        allPileupHits[NPileupHits++] = aPileupHit; // save pileup hit in global array
        thisHit++;                                 // skip nextHit (we just used it)
        // not looking for a third hit at this point
      }

      // ### *not* a pileup hit
      else
      {
        if (debug_)
          mf::LogVerbatim("RatioEast") << "no DTC: " << thisHit->summary();
        PileupHit aPileupHit(thisHit);
        allPileupHits[NPileupHits++] = aPileupHit;
        if (thisHit == allHits + NAllHits - 2)
        {
          // ### *never* a pileup hit
          if (debug_)
            mf::LogVerbatim("RatioEast") << "last hit: " << nextHit->summary();
          PileupHit aPileupHit(nextHit);
          allPileupHits[NPileupHits++] = aPileupHit;
        }
      }
    }

    // ### what did we find? ###
    if (debug_)
    {
      mf::LogVerbatim("RatioEast")
          << "Finished pileup scan!\n"
          << "  calorimeterIndex: " << clustersThisCalo.calorimeterIndex << "\n"
          << "  allHits: " << NAllHits << " entries\n"
          << "  allPileupHits: " << NPileupHits << " entries";
    }

    // ### histogram all of the things! (loop over pileup hits) ###
    for (
        PileupHit *thisPUHit = allPileupHits;
        thisPUHit < allPileupHits + NPileupHits;
        thisPUHit++)
    {
      // use these if you want names for the the next two hits
      //PileupHit * secondPUHit = thisPUHit+1;
      //PileupHit * thirdPUHit = thisPUHit+2;

      // ### code below does the gaptime+deadtime search etc ###
      // NOTE: The pileup hit doesn't have its own defined random numbers for
      // 'bin smearing' and 'quartering' because those are all associated with
      // the cluster/data product hist.  The intent here is to use the cluster
      // hit times & quartering choices for each pileup hit like this:
      //    * time cuts: use getAvgAnalysisTime(), the cluster hit analysis
      //      time (averaged over cluster hits when the pileup hit has multiple
      //      cluster hits)
      //    * search for combos of this+next+nextnext hit (as well as the gaptime
      //      and deadtime cut excluding some combos): getAvgAnalysisTime()
      //    * TIME_D/DTimes_calo_All: getAvgRandomizedAnalysisTime() plus half
      //      gap time
      //    * TIME_S1/S2 and S1/S2Times_calo_All: getAvgRandomizedAnalysisTime()
      //    *

      // skip hits that don't pass time & space (calo xtals) cuts
      if (!hitPassesTimeCuts(thisPUHit, tMin_, tMax_, g2HalfPeriod))
      {
        if (debug_)
          LogVerbatim("RatioEast") << "pileup hit failed time cuts" << endl;
        continue;
      }
      if (!hitPassesCaloXYCuts(thisPUHit))
      {
        if (debug_)
          LogVerbatim("RatioEast") << "pileup hit failed calo XY cuts" << endl;
        continue;
      }

      // skip last element because we're looking for triples
      if (thisPUHit >= allPileupHits + NPileupHits - 2)
      {
        if (debug_)
          LogVerbatim("RatioEast") << "skipping last pileup hit..." << endl;
        break;
      }

      // let's name these \Delta t_{j,i}
      double dt21 = (thisPUHit + 1)->getAvgAnalysisTime() - thisPUHit->getAvgAnalysisTime();
      double dt31 = (thisPUHit + 2)->getAvgAnalysisTime() - thisPUHit->getAvgAnalysisTime();
      double dt32 = dt31 - dt21;
      if (debug_)
        LogVerbatim("RatioEast") << "dt21=" << dt21 << " dt31=" << dt31 << " dt32=" << dt32
                                 << " gapT_=" << gapT_ << " DT_=" << DT_ << endl;
      char tmpstr[256];
      sprintf(tmpstr, "dt21=%.5f dt31=%.5f dt32=%.5f gapT_=%.5f DT_=%.5f",
              dt21, dt31, dt32, gapT_, DT_);
      if (debug_)
        LogVerbatim("RatioEast") << tmpstr;

      // skip if t2-t1 is NOT between gaptime and gaptime+deadtime
      if (dt21 >= gapT_ + DT_ || dt21 < gapT_)
      {
        if (debug_)
          LogVerbatim("RatioEast") << "pileup hit has (dt21>=gapT_+DT_ || dt21<gapT_)" << endl;
        continue;
      }

      // now let's start comparing to the following hit
      if (
          dt31 > 2. * (gapT_ + DT_)                      // category 1
          || dt32 < gapT_                                // category 2
          || (dt31 < 2. * (gapT_ + DT_) && dt32 > gapT_) // category 3
      )
      {
        double E1 = thisPUHit->getEnergySum();
        double E2 = (thisPUHit + 1)->getEnergySum();
        double T1_rand = thisPUHit->getAvgRandomizedAnalysisTime();
        double T2_rand = (thisPUHit + 1)->getAvgRandomizedAnalysisTime();
        double Td_rand = T1_rand + gapT_ / 2.; // TODO: finish post-March bugfixes & upgrades?

        //if (debug_) LogVerbatim("RatioEast")<<"T1="<<T1<<" T2="<<T2<<" E1="<<" E2="<<E2<<endl;
        sprintf(tmpstr,
                "T1=%.5f T2=%.5f Td_rand=%.5f E1=%.2f E2=%.2f",
                T1_rand, T2_rand, Td_rand, E1, E2);
        if (debug_)
          LogVerbatim("RatioEast") << tmpstr << endl;

        // full T-method pileup time histograms
        // (these get randomized times)
        if (E1 > Emin_ && E2 > Emin_)
        {
          if (debug_)
            LogVerbatim("RatioEast") << "one" << endl;
          TIME_D->Fill(Td_rand);
          TIME_S1->Fill(T1_rand);
          TIME_S2->Fill(T1_rand);
          ENERGY_D->Fill(E1 + E2);
          ENERGY_S1->Fill(E1);
          ENERGY_S2->Fill(E2);
          energy_d_All[iCalo]->Fill(E1 + E2);
          energy_s1_All[iCalo]->Fill(E1);
          energy_s2_All[iCalo]->Fill(E2);
          DTimes_calo_All[iCalo]->Fill(Td_rand);
          S1Times_calo_All[iCalo]->Fill(T1_rand);
          S2Times_calo_All[iCalo]->Fill(T1_rand);
        }
        // TODO: this 2D hist code block (and the next three)
        //if(E1 > Ethre_ && E2 > Ethre_) {
        //  ET_D->Fill(Td+randombintime_pile[i][c], E1+E2);
        //  ET_S1->Fill(T1+randombintime_pile[i][c], E1);
        //  ET_S2->Fill(T1+randombintime_pile[i][c], E2);
        //}
        else if (E2 < Emin_ && E1 > Emin_)
        {
          if (debug_)
            LogVerbatim("RatioEast") << "two" << endl;
          TIME_D->Fill(Td_rand);
          TIME_S1->Fill(T1_rand);
          ENERGY_D->Fill(E1 + E2);
          ENERGY_S1->Fill(E1);
          energy_d_All[iCalo]->Fill(E1 + E2);
          energy_s1_All[iCalo]->Fill(E1);
          DTimes_calo_All[iCalo]->Fill(Td_rand);
          S1Times_calo_All[iCalo]->Fill(T1_rand);
        }
        //if(E2 < Ethre_ && E1 > Ethre_){
        //  ET_D->Fill(Td+randombintime_pile[i][c], E1+E2);
        //  ET_S1->Fill(T1+randombintime_pile[i][c], E1);
        //}
        else if (E2 > Emin_ && E1 < Emin_)
        {
          if (debug_)
            LogVerbatim("RatioEast") << "three" << endl;
          TIME_D->Fill(Td_rand);
          TIME_S2->Fill(T1_rand);
          ENERGY_D->Fill(E1 + E2);
          ENERGY_S2->Fill(E2);
          energy_d_All[iCalo]->Fill(E1 + E2);
          energy_s2_All[iCalo]->Fill(E2);
          DTimes_calo_All[iCalo]->Fill(Td_rand);
          S2Times_calo_All[iCalo]->Fill(T1_rand);
        }
        //if(E2 > Ethre_ && E1 < Ethre_){
        //  ET_D->Fill(Td+randombintime_pile[i][c], E1+E2);
        //  ET_S2->Fill(T1+randombintime_pile[i][c], E2);
        //}
        else if (E1 < Emin_ && E2 < Emin_ && E1 + E2 > Emin_)
        {
          if (debug_)
            LogVerbatim("RatioEast") << "four" << endl;
          TIME_D->Fill(Td_rand);
          ENERGY_D->Fill(E1 + E2);
          energy_d_All[iCalo]->Fill(E1 + E2);
          DTimes_calo_All[iCalo]->Fill(Td_rand);
        }
        //if(E1 < Ethre_ && E2 < Ethre_){
        //  if((E1+E2)>Ethre_){
        //     ET_D->Fill(Td+randombintime_pile[i][c], E1+E2);
        //   }
        // }

        // ### pileup histograms for 4-spectra to create ratio ###
        //ClusterHit * firstClusterHit = thisPUHit->sources[0];
        //const int hitType = int(firstClusterHit->quarteringConstant);
        const int hitType = int(thisPUHit->getQuarteringConstant());
        // hit->quarteringConstant should be [0.,4.)
        if (hitType == 0)
        {
          pileup_hist_up(
              (thisPUHit->getAvgRandomizedAnalysisTime()) + ratioTimeShiftDefinitions[hitType],
              T2,
              E1, E2,
              ENERGY_D_up, ENERGY_S1_up, ENERGY_S2_up,
              DTimes_Up_all, S1Times_Up_all, S2Times_Up_all,
              energy_d_Up, energy_s1_Up, energy_s2_Up,
              DTimes_calo_Up, S1Times_calo_Up, S2Times_calo_Up);
        }

        else if (hitType == 1)
        {
          pileup_hist_um(
              (thisPUHit->getAvgRandomizedAnalysisTime()) + ratioTimeShiftDefinitions[hitType],
              T2,
              E1, E2,
              ENERGY_D_um, ENERGY_S1_um, ENERGY_S2_um,
              DTimes_Um_all, S1Times_Um_all, S2Times_Um_all,
              energy_d_Um, energy_s1_Um, energy_s2_Um, DTimes_calo_Um, S1Times_calo_Um, S2Times_calo_Um);
        }

        else if (hitType == 2)
        {
          pileup_hist_vp(
              (thisPUHit->getAvgRandomizedAnalysisTime()) + ratioTimeShiftDefinitions[hitType],
              T2,
              E1, E2,
              ENERGY_D_vp, ENERGY_S1_vp, ENERGY_S2_vp,
              DTimes_Vp_all, S1Times_Vp_all, S2Times_Vp_all,
              energy_d_Vp, energy_s1_Vp, energy_s2_Vp,
              DTimes_calo_Vp, S1Times_calo_Vp, S2Times_calo_Vp);
        }

        else
        {
          pileup_hist_vm(
              (thisPUHit->getAvgRandomizedAnalysisTime()) + ratioTimeShiftDefinitions[hitType],
              T2,
              E1, E2,
              ENERGY_D_vm, ENERGY_S1_vm, ENERGY_S2_vm,
              DTimes_Vm_all, S1Times_Vm_all, S2Times_Vm_all,
              energy_d_Vm, energy_s1_Vm, energy_s2_Vm,
              DTimes_calo_Vm, S1Times_calo_Vm, S2Times_calo_Vm);
        }
        /*
*/
      } // end if-statement for categories 1-3
    }   // end of histogramming loop

  } // end loop over CaloGlobalFitView records

  //////////////// REFERENCE CODE //////////////////////

  /*
    Original setup of time[]/time_pile[] and filling of histograms:
  
  // Populate randomno, randombintime, TIME_FULL, ENERGY_FULL, 
  // caloEnergies_All[iCalo], caloTimes_All[iCalo], caloEnergiesUVpm[iCalo],
  // and caloTimesUVpm[iCalo].
  // This saves every cluster (plus two random numbers each) but only puts
  // cluster time/energy in the histograms if it passes time cuts.
  // (TODO: energy cuts don't matter?)
  // Random numbers are created here & saved to randomno and randombintime,
  // and time shift is included
  vector<vector<double> > time(NCalos);
  vector<vector<double> > energy(NCalos);
  vector<vector<double> > island(NCalos);
  vector<vector<double> > xposition(NCalos);
  vector<vector<double> > yposition(NCalos);
  vector<vector<double> > randombintime(NCalos);
  vector<vector<double> > randomno(NCalos);
  time.clear();
  energy.clear();
  island.clear();
  xposition.clear();
  yposition.clear();
  randombintime.clear();
  randomno.clear();
  for(auto hit : allHits[some iCalo?]) {
    
    //for(unsigned int iCalo=0; iCalo<NCalos; iCalo++) {
        //if(iCalo==hit->iCalo) {
    // TODO: eliminate extra vars with dangerous names?
    int iCalo = hit->iCalo;
    double randomBinTime = random2.Uniform() * 0.14915 - 0.14915/2;
    double rand_time = hit->time+randomBinTime;
    double randomNum = random1.Uniform();
    double halfPeriod = g2Period/2;
    
    randomno[iCalo].push_back(randomNum);
    randombintime[iCalo].push_back(randomBinTime);
    
    // debug output
    if (debug_) {
      if(iCalo==21) {
        if (randomno.size()%20==1) {
          cout << "DEBUG: calo 21, every 20th cluster\n";
          //cout<<" random2= "<< random2.Uniform() * 0.14915 <<"random bin time test = "<< random2.Uniform() * 0.14915 - 0.14915/2 <<" randomBinTime = "<<randomBinTime<<" random1 ="<<random1.Uniform()<<" randomNum = "<<randomNum<<std::endl;
          int i = randomno.size()-1;
          cout << "  hit time: "<<hit->time<<endl;
          cout << "  xcluster: "<<hit->xcluster<<endl;
          cout << "  i: "<<i<<endl;
          cout << "  xposition[21][i]: "<<xposition[21][i]<<endl;
          cout << "  randombintime[21][i]: "<<randombintime[21][i]<<endl;
          cout << "  randomno[21][i]: "<<randomno[21][i]<<endl;
        }
      }
    }
    
    // save cluster info here (TODO: even if it doesn't pass time cuts?)
    // these times are NOT shifted for FR
    time[iCalo].push_back(hit->time);
    energy[iCalo].push_back(hit->energy);
    island[iCalo].push_back(hit->islandNum);
    xposition[iCalo].push_back(hit->xcluster);
    yposition[iCalo].push_back(hit->ycluster);
    
    if(hit->time>20 && hit->time<651 && hit->energy > Emin_) {

      TIME_FULL->Fill(rand_time);
      ENERGY_FULL->Fill(hit->energy);
      caloEnergies_All[iCalo]->Fill(hit->energy);
      caloTimes_All[iCalo]->Fill(rand_time);
      
      // quarter 1 (U+)
      if(randomNum < 0.25) {
        ENERGY_FULL_up->Fill(hit->energy);
        caloEnergies_Up[iCalo]->Fill(hit->energy);
        caloTimes_Up_all->Fill(rand_time - halfPeriod);
        caloTimes_Up[iCalo]->Fill(rand_time - halfPeriod);
      }
      
      // quarter 2 (U-)
      if(0.25 <= randomNum && randomNum < 0.5) {
        ENERGY_FULL_um->Fill(hit->energy);
        caloEnergies_Um[iCalo]->Fill(hit->energy);
        caloTimes_Um_all->Fill(rand_time + halfPeriod);
        caloTimes_Um[iCalo]->Fill(rand_time + halfPeriod);
      }
      
      // quarter 3 (V1)
      if(0.5 <= randomNum && randomNum < 0.75) {
        ENERGY_FULL_vp->Fill(hit->energy);
        caloEnergies_Vp[iCalo]->Fill(hit->energy);
        caloTimes_Vp_all->Fill(rand_time);
        caloTimes_Vp[iCalo]->Fill(rand_time);
      }
      
      // quarter 4 (V2)
      else {
        ENERGY_FULL_vm->Fill(hit->energy);
        caloEnergies_Vm[iCalo]->Fill(hit->energy);
        caloTimes_Vm_all->Fill(rand_time);
        caloTimes_Vm[iCalo]->Fill(rand_time);
      }

    }//wiggle time cut
    
      //}
    //}//loop over calonums

  }//clusterhit vector
  */

  /* 
    Original search for DT coincidences and population of time_pile[]/etc
  // create pileup histograms
  vector<vector<double> > time_pile(NCalos);
  vector<vector<double> > energy_pile(NCalos);

  for(unsigned int iCalo=0; iCalo<NCalos; iCalo++) {

    for(unsigned int jCluster=0; jCluster<time[iCalo].size()-1; jCluster++) {
      
      // skip anything outside time range
      if(!(time[iCalo][jCluster]>20 && time[iCalo][jCluster]<651)) continue;
      
      // if second is in first's dead time, and they're in the same island
      // then add a combined entry to time_pile and energy_pile and skip jCluster
      // else add to time/energy_pile individually
      // TODO: what happens if there are *two* within DT_? This loop doesn't 
      // treat the third positron as part of the pile?
      if(time[iCalo][jCluster+1]-time[iCalo][jCluster]<DT_ && island[iCalo][jCluster+1]==island[iCalo][jCluster]) {
        En=energy[iCalo][jCluster+1]+energy[iCalo][jCluster];
        energy_pile[iCalo].push_back(En);

        Tn = (time[iCalo][jCluster+1]*energy[iCalo][jCluster+1]+time[iCalo][jCluster]*energy[iCalo][jCluster])/En;
        time_pile[iCalo].push_back(Tn);

        jCluster++;

      }
      else {
        time_pile[iCalo].push_back(time[iCalo][jCluster]);
        energy_pile[iCalo].push_back(energy[iCalo][jCluster]);
        if(jCluster==time[iCalo].size()-2) {
          time_pile[iCalo].push_back(time[iCalo][jCluster+1]);
          energy_pile[iCalo].push_back(energy[iCalo][jCluster+1]);

        }
      }

      }//clusters
    }//calos
  } reference, do not uncomment
  */

  //  // loop over calos, then over entries in time/energy_pile
  //  // TODO: association between positrons in elements of randombintime,
  //  // randomno cannot be used here because time/energy_pile have fewer
  //  // entries than time/energy (because of the combined entries).
  //  for(unsigned int iCalo=0; iCalo<NCalos; iCalo++) {
  //    for(unsigned int c=0; c<time_pile[iCalo].size(); c++) {

  //      double randomBinTime = randombintime[iCalo][c];
  //      double randomNum = randomno[iCalo][c];
  //
  //      // debug output
  //      if(debug_) {
  //        if (iCalo==21) {
  //          if (c%20==1) {
  //            cout << "DEBUG: calo 21, every 20th cluster\n";
  //            //cout<<" random2= "<< random2.Uniform() * 0.14915 <<"random bin time test = "<< random2.Uniform() * 0.14915 - 0.14915/2 <<" randomBinTime = "<<randomBinTime<<" random1 ="<<random1.Uniform()<<" randomNum = "<<randomNum<<std::endl;
  //            int i = c;
  //            cout << "  i: "<<i<<endl;
  //            cout << "  xposition[21][i]: "<<xposition[21][i]<<endl;
  //            cout << "  randombintime[21][i]: "<<randombintime[21][i]<<endl;
  //            cout << "  randomno[21][i]: "<<randomno[21][i]<<endl;
  //          }
  //        }
  //      }
  //
  //      //if(energy_pile[iCalo][c] > 500 && energy_pile[iCalo][c+1] > 500 && time_pile[iCalo][c+1] > 29 && time_pile[iCalo][c+1] < 651 && time_pile[iCalo][c] > 29 && time_pile[iCalo][c] < 651){
  //      //deltat->Fill((time_pile[iCalo][c+1]-time_pile[iCalo][c])*1000);
  //      //if(energy_pile[iCalo][c] > 1700 && energy_pile[iCalo][c+1] > 1700){
  //      //deltaxvsdeltatallcalo->Fill((time_pile[iCalo][c+1]-time_pile[iCalo][c])*1000,fabs(xposition[iCalo][c+1]-xposition[iCalo][c]));
  //      //deltayvsdeltatallcalo->Fill((time_pile[iCalo][c+1]-time_pile[iCalo][c])*1000,fabs(yposition[iCalo][c+1]-yposition[iCalo][c]));
  //      //deltaxvsdeltat[iCalo]->Fill((time_pile[iCalo][c+1]-time_pile[iCalo][c])*1000,fabs(xposition[iCalo][c+1]-xposition[iCalo][c]));
  //      //deltayvsdeltat[iCalo]->Fill((time_pile[iCalo][c+1]-time_pile[iCalo][c])*1000,fabs(yposition[iCalo][c+1]-yposition[iCalo][c]));
  //      //}
  //      //}
  //
  //
  //      // if these conditions do not apply then we're done with this iteration of the loop
  //      if(
  //        !(
  //          time_pile[iCalo][c]>20
  //          && time_pile[iCalo][c]<651
  //          && fabs(xposition[iCalo][c+1]-xposition[iCalo][c]) <DX_
  //          && c<time_pile[iCalo].size()-2
  //        )
  //      ) continue;
  //
  //
  //      // do the same for all three categories
  //      if(
  //        // common condition to all three categories
  //        time_pile[iCalo][c+1]-time_pile[iCalo][c]<gapT_+DT_ && time_pile[iCalo][c+1]-time_pile[iCalo][c]>gapT_
  //        && (
  //          // category 1
  //          ( time_pile[iCalo][c+2]-time_pile[iCalo][c]>2.*(DT_+gapT_) )
  //          // category 2
  //          || ( time_pile[iCalo][c+2]-time_pile[iCalo][c+1]<gapT_ )
  //          // category 3
  //          || ( time_pile[iCalo][c+2]-time_pile[iCalo][c]<2.*(DT_+gapT_) && time_pile[iCalo][c+2]-time_pile[iCalo][c+1]>gapT_ )
  //        )
  //      ) {
  //
  //        E1 = energy_pile[iCalo][c];
  //        E2 = energy_pile[iCalo][c+1];
  //
  //        T1 = time_pile[iCalo][c]+randomBinTime;
  //        T2 = time_pile[iCalo][c+1];
  //
  //        //full T-method pileup time histograms
  //        if(E1 > Emin_ && E2 > Emin_) {
  //          TIME_D->Fill(T1);
  //          TIME_S1->Fill(T1);
  //          TIME_S2->Fill(T1);
  //          ENERGY_D->Fill(E1+E2);
  //          ENERGY_S1->Fill(E1);
  //          ENERGY_S2->Fill(E2);
  //          energy_d_All[iCalo]->Fill(E1+E2);
  //          energy_s1_All[iCalo]->Fill(E1);
  //          energy_s2_All[iCalo]->Fill(E2);
  //          DTimes_calo_All[iCalo]->Fill(T1);
  //          S1Times_calo_All[iCalo]->Fill(T1);
  //          S2Times_calo_All[iCalo]->Fill(T1);
  //        }
  //        else if(E2 < Emin_ && E1 > Emin_) {
  //          TIME_D->Fill(T1);
  //          TIME_S1->Fill(T1);
  //          ENERGY_D->Fill(E1+E2);
  //          ENERGY_S1->Fill(E1);
  //          energy_d_All[iCalo]->Fill(E1+E2);
  //          energy_s1_All[iCalo]->Fill(E1);
  //          DTimes_calo_All[iCalo]->Fill(T1);
  //          S1Times_calo_All[iCalo]->Fill(T1);
  //        }
  //        else if(E2 > Emin_ && E1 < Emin_) {
  //          TIME_D->Fill(T1);
  //          TIME_S2->Fill(T1);
  //          ENERGY_D->Fill(E1+E2);
  //          ENERGY_S2->Fill(E2);
  //          energy_d_All[iCalo]->Fill(E1+E2);
  //          energy_s2_All[iCalo]->Fill(E2);
  //          DTimes_calo_All[iCalo]->Fill(T1);
  //          S2Times_calo_All[iCalo]->Fill(T1);
  //        }
  //        else if (E1 < Emin_ && E2 < Emin_ && E1+E2>Emin_) {
  //          TIME_D->Fill(T1);
  //          ENERGY_D->Fill(E1+E2);
  //          energy_d_All[iCalo]->Fill(E1+E2);
  //          DTimes_calo_All[iCalo]->Fill(T1);
  //        }

  //        //pileup histograms for 4-spectra to create ratio
  //        if(randomNum < 0.25) {
  //          T1_up = (time_pile[iCalo][c])-g2Period/2;
  //          pileup_hist_up( T1_up, T2, E1, E2, ENERGY_D_up, ENERGY_S1_up, ENERGY_S2_up, DTimes_Up_all, S1Times_Up_all, S2Times_Up_all, energy_d_Up, energy_s1_Up, energy_s2_Up,DTimes_calo_Up,S1Times_calo_Up,S2Times_calo_Up);
  //        }

  //        else if(0.25 <= randomNum && randomNum < 0.5) {
  //          T1_um = (time_pile[iCalo][c])+g2Period/2;
  //          pileup_hist_um( T1_um, T2, E1, E2, ENERGY_D_um, ENERGY_S1_um, ENERGY_S2_um, DTimes_Um_all, S1Times_Um_all, S2Times_Um_all, energy_d_Um, energy_s1_Um, energy_s2_Um,DTimes_calo_Um,S1Times_calo_Um,S2Times_calo_Um);
  //        }

  //        else if(0.5 <= randomNum && randomNum < 0.75) {
  //          pileup_hist_vp( T1, T2, E1, E2, ENERGY_D_vp, ENERGY_S1_vp, ENERGY_S2_vp, DTimes_Vp_all, S1Times_Vp_all, S2Times_Vp_all, energy_d_Vp, energy_s1_Vp, energy_s2_Vp,DTimes_calo_Vp,S1Times_calo_Vp,S2Times_calo_Vp);
  //        }

  //        else {
  //          pileup_hist_vm( T1, T2, E1, E2, ENERGY_D_vm, ENERGY_S1_vm, ENERGY_S2_vm, DTimes_Vm_all, S1Times_Vm_all, S2Times_Vm_all, energy_d_Vm, energy_s1_Vm, energy_s2_Vm,DTimes_calo_Vm,S1Times_calo_Vm,S2Times_calo_Vm);
  //        }
  //      }
  //    }
  //  } //calo

} //analyze
DEFINE_ART_MODULE(ratioeast)
