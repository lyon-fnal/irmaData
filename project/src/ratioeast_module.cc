/*Art analyzer module to produce pileup corrected ratio histograms.
@author  Sudeshna Ganguly
@date    10/21/2019*/
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "gm2dataproducts/daq/MidasEventHeaderArtRecord.hh"
#include "gm2dataproducts/daq/CCCArtRecord.hh"
#include "TTree.h"
#include "gm2dataproducts/reconeast/GlobalFitArtRecord.hh"
#include "TH2D.h"
#include <TRandom3.h>
#include "ReconEast.hh"



class hitclass {
public:
  int caloNum;
  int islandNum;
  double time;
  double energy;
  double xcluster;
  double ycluster;
};

/*int randSeed1 = 5;
int randSeed2 = 5;

TRandom3* random1 = new TRandom3(randSeed1);
TRandom3* random2 = new TRandom3(randSeed2);
*/


///////////////////////////                                                                                                                                   
class ratioeast;
class ratioeast : public art::EDAnalyzer {

public:
  explicit ratioeast(fhicl::ParameterSet const & p);

  ratioeast(ratioeast const &) = delete;
  ratioeast(ratioeast &&) = delete;
  ratioeast & operator = (ratioeast const &) = delete;
  ratioeast & operator = (ratioeast &&) = delete;

  // Required functions.                                                                                                                                      
  void analyze(art::Event const & e) override;
  //void endJob() override;                   

private:

  art::ServiceHandle<art::TFileService> tfs_;
  std::string readModule_, readInstance_;
  const int Emin_;
  const int Ethre_;
  const double binWidth_;
  const double gapT_;
  const double DT_;
  const double DX_;
  const double DY_;
  const int c1_;
  const int c2_;
  //const int randSeed1_; const int randSeed2_;  

  //TRandom3 random1;
  //TRandom3 random2;
  
  TRandom3 random1;
  TRandom3 random2;

  const double EnergyCal_run1[24] = {
    1628.9, 1505.9, 1559.4, 1564.9, 1368.8, 1516.9, 1543.8, 1533.0,
    1518.1, 1551.7, 1582.6, 1610.8, 1604.2, 1566.5, 1528.0, 1487.0,
    1520.0, 1588.2, 1554.9, 1525.5, 1455.7, 1474.9, 1522.5, 1548.1
  };  

  const double EnergyCal_run2[24] = {
    1845.34, 1956.21, 1852.62, 1882.91, 2075.44, 1919.23, 1885.50, 1900.13,
    1893.89, 1889.55, 1880.18, 1906.99, 1920.32, 1910.87, 1913.59, 1963.29,
    1973.80, 1931.09, 1932.92, 1943.34, 1976.45, 1964.30, 1928.10, 1933.15
  };
  
  const double MaxTimeBin = 700000;
  const int nBins = int(MaxTimeBin/binWidth_);
  const double histMaxTime = nBins*binWidth_/1000;
  const double Fa = 0.2290735 * 1e6 * 1/(1e9);                                                                                                                     
  const double g2Period = (1/Fa)/1000;
  const double Lifetime = 64400; // ns 
  int Ncalo=24;
 
  /*unsigned int eventnum;
  unsigned int runnum;
  unsigned int subRunnum;
  */

 
  
  double E1;
  double E2;
  double T1_up;
  double T1_um;
  double T1_v1;
  double T1_v2;
  double T1;
  double T2;
  double Td;
  double Tn, En;

  /*TTree *t_;

  double clusterEnergy_;
  double clusterTime_;
  double clusterY_;
  double clusterX_;
  unsigned int caloNum_;*/
  //unsigned int eventNum_;
  //unsigned int bunchNum_;
  //unsigned int midasSerialNum_;
  //unsigned int subRunNum_;
  //unsigned int runNum_;
  
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

  TH2F *ETspectrum;
  TH2F *ET_D;
  TH2F *ET_S1;
  TH2F *ET_S2;

  /*  TH2F *XY;
  TH1F *X;
  TH1F *Y;
  TH1F *r;
  TH1F *deltax[11];

  TH1F *deltat;
  TH1F *deltat_calo[24];
  TH2F *deltaxvsdeltatallcalo;
  TH2F *deltayvsdeltatallcalo;
 

  TH2F *deltaxvsdeltat[24];
  TH2F *deltayvsdeltat[24];
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

  TH1F* caloEnergies_All[24];
  TH1F* caloTimes_All[24];

  TH1F *energy_d_All[24];
  TH1F *energy_s1_All[24];
  TH1F *energy_s2_All[24];

  TH1F *DTimes_calo_All[24];
  TH1F *S1Times_calo_All[24];
  TH1F *S2Times_calo_All[24];

  TH1F* caloEnergies_Up[24];
  TH1F* caloEnergies_Um[24];
  TH1F* caloEnergies_Vp[24];
  TH1F* caloEnergies_Vm[24];

  TH1F* caloTimes_Up[24];
  TH1F* caloTimes_Um[24];
  TH1F* caloTimes_Vp[24];
  TH1F* caloTimes_Vm[24];

  TH1F *energy_d_Up[24];
  TH1F *energy_d_Um[24];
  TH1F *energy_d_Vp[24];
  TH1F *energy_d_Vm[24];

  TH1F *energy_s1_Up[24];
  TH1F *energy_s1_Um[24];
  TH1F *energy_s1_Vp[24];
  TH1F *energy_s1_Vm[24];

  TH1F *energy_s2_Up[24];
  TH1F *energy_s2_Um[24];
  TH1F *energy_s2_Vp[24];
  TH1F *energy_s2_Vm[24];

  TH1F *DTimes_calo_Up[24];
  TH1F *DTimes_calo_Um[24];
  TH1F *DTimes_calo_Vp[24];
  TH1F *DTimes_calo_Vm[24];

  TH1F *S1Times_calo_Up[24];
  TH1F *S1Times_calo_Um[24];
  TH1F *S1Times_calo_Vp[24];
  TH1F *S1Times_calo_Vm[24];

  TH1F *S2Times_calo_Up[24];
  TH1F *S2Times_calo_Um[24];
  TH1F *S2Times_calo_Vp[24];
  TH1F *S2Times_calo_Vm[24];

  void  pileup_hist_up(double T1_up, double T2, double E1, double E2, TH1F *ENERGY_D_up, TH1F *ENERGY_S1_up, TH1F *ENERGY_S2_up, TH1F *DTimes_Up_all, TH1F *S1Times_Up_all, TH1F *S2Times_Up_all, TH1F **energy_d_Up, TH1F **energy_s1_Up, TH1F **energy_s2_Up,TH1F **DTimes_calo_Up,TH1F **S1Times_calo_Up,TH1F **S2Times_calo_Up); 
  void  pileup_hist_um(double T1_um, double T2, double E1, double E2, TH1F *ENERGY_D_um, TH1F *ENERGY_S1_um, TH1F *ENERGY_S2_um, TH1F *DTimes_Um_all, TH1F *S1Times_Um_all, TH1F *S2Times_Um_all, TH1F **energy_d_Um, TH1F **energy_s1_Um, TH1F **energy_s2_Um,TH1F **DTimes_calo_Um,TH1F **S1Times_calo_Um,TH1F **S2Times_calo_Um);
  void  pileup_hist_vp(double T1_v1, double T2, double E1, double E2, TH1F *ENERGY_D_vp, TH1F *ENERGY_S1_vp, TH1F *ENERGY_S2_vp, TH1F *DTimes_Vp_all, TH1F *S1Times_Vp_all, TH1F *S2Times_Vp_all, TH1F **energy_d_Vp, TH1F **energy_s1_Vp, TH1F **energy_s2_Vp,TH1F **DTimes_calo_Vp,TH1F **S1Times_calo_Vp,TH1F **S2Times_calo_Vp);
  void  pileup_hist_vm(double T1_v2, double T2, double E1, double E2, TH1F *ENERGY_D_vm, TH1F *ENERGY_S1_vm, TH1F *ENERGY_S2_vm, TH1F *DTimes_Vm_all, TH1F *S1Times_Vm_all, TH1F *S2Times_Vm_all, TH1F **energy_d_Vm, TH1F **energy_s1_Vm, TH1F **energy_s2_Vm,TH1F **DTimes_calo_Vm,TH1F **S1Times_calo_Vm,TH1F **S2Times_calo_Vm);

};

ratioeast::ratioeast(fhicl::ParameterSet const & p)

  :

  EDAnalyzer(p),
  readModule_ (p.get<std::string>("readModule", "energyPartition")),
  readInstance_ (p.get<std::string>("readInstance", "partition")),
  Emin_(p.get<int>("Emin",1700)),
  Ethre_(p.get<int>("Ethre",500)),
  binWidth_(p.get<double>("binWidth",149.15)),
  gapT_(p.get<double>("gapT",0.020)),
  DT_(p.get<double>("DeadT",0.001)),
  DX_(p.get<double>("DeadX",3.0)),
  DY_(p.get<double>("DeadY",3.0)),
  c1_(p.get<int>("c1",0)),
  c2_(p.get<int>("c2",0)),
  //randSeed1_(p.get<int>("randSeed1",c1_*fillId)),
  //randSeed2_(p.get<int>("randSeed2",c2_*fillId)),
  //randSeed1_(p.get<int>("randSeed1",c1_*eventnum)),
  //randSeed1_(p.get<int>("randSeed1",2)),
  //randSeed2_(p.get<int>("randSeed2",2)),
  random1(),
  random2(),
  E1(0),
  E2(0),
  T1_up(0),
  T1_um(0),
  T1_v1(0),
  T1_v2(0),
  T1(0),
  T2(0),
  Td(0),
  Tn(0),
  En(0),
  //clusterEnergy_(0),
  //clusterTime_(0),
  //clusterY_(0),
  //clusterX_(0),
  //caloNum_(0),
  //eventNum_(0),
  //bunchNum_(-1),
  //midasSerialNum_(0),
  //subRunNum_(0),
  //runNum_(0),
  datasetName_reconEast_(p.get<std::string>("datasetName_reconEast", "60h"))
{

  random1.SetSeed(0);
  random2.SetSeed(0);
  

  art::ServiceHandle<art::TFileService> tfs;
  /* t_ = tfs->make<TTree>("testTree", "testTree");
  t_->Branch("energy", &clusterEnergy_, "energy/D");
  t_->Branch("time", &clusterTime_, "time/D");
  t_->Branch("Y", &clusterY_, "y/D");
  t_->Branch("X", &clusterX_, "x/D");
  t_->Branch("caloNum", &caloNum_, "caloNum/i");
  */
  /*t_->Branch("eventNum", &eventNum_, "eventNum/i");
    t_->Branch("bunchNum", &bunchNum_, "bunchNum/i");
    t_->Branch("midasSerialNum", &midasSerialNum_, "midasSerialNum/i");
    t_->Branch("subRunNum", &subRunNum_, "subRunNum/i");
    t_->Branch("runNum", &runNum_, "runNum/i");*/
  
  TIME_FULL = tfs->make<TH1F>("TIME_FULL", "TIME_FULL", nBins, 0, histMaxTime);
  ENERGY_FULL =  tfs->make<TH1F>("ENERGY_FULL","ENERGY ON ALL CALO; energy [MeV]",500,0,10000);
  TIME_D =  tfs->make<TH1F>("TIME_D","TIME_D ON ALL CALO;time [#mus];counts/149ns",nBins,0,histMaxTime);
  TIME_S1 =  tfs->make<TH1F>("TIME_S1","TIME_S1 ON ALL CALO;time [#mus];counts/149ns",nBins,0,histMaxTime);
  TIME_S2 =  tfs->make<TH1F>("TIME_S2","TIME_S2 ON ALL CALO;time [#mus];counts/149ns",nBins,0,histMaxTime);

  ENERGY_D =  tfs->make<TH1F>("ENERGY_D","ENERGY_D ON ALL CALO; energy [MeV]",500,0,10000);
  ENERGY_S1 =  tfs->make<TH1F>("ENERGY_S1","ENERGY_S1 ON ALL CALO; energy [MeV]",500,0,10000);
  ENERGY_S2 =  tfs->make<TH1F>("ENERGY_S2","ENERGY_S2 ON ALL CALO; energy [MeV]",500,0,10000);

  ETspectrum = tfs->make<TH2F>("ETspectrum", "Energy and Time Spectrum;time [#mus];energy [MeV]",nBins, 0, histMaxTime,2000,0,10000);
  ET_D = tfs->make<TH2F>("ETspectrum_D", "Energy and Time Spectrum;time [#mus] for doubles; energy[MeV]",nBins, 0, histMaxTime,2000,0,10000);
  ET_S1 = tfs->make<TH2F>("ETspectrum_S1", "Energy and Time Spectrum;time [#mus] for singles 1; energy [MeV]",nBins, 0, histMaxTime,2000,0,10000);
  ET_S2 = tfs->make<TH2F>("ETspectrum_S2", "Energy and Time Spectrum;time [#mus] for singles 2; energy [MeV]",nBins, 0, histMaxTime,2000,0,10000);
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

  ENERGY_FULL_up =  tfs->make<TH1F>("ENERGY_FULL_up","ENERGY_D ON ALL CALO; energy [MeV]",500,0,10000);
  ENERGY_FULL_um =  tfs->make<TH1F>("ENERGY_FULL_um","ENERGY_D ON ALL CALO; energy [MeV]",500,0,10000);
  ENERGY_FULL_vp =  tfs->make<TH1F>("ENERGY_FULL_vp","ENERGY_D ON ALL CALO; energy [MeV]",500,0,10000);
  ENERGY_FULL_vm =  tfs->make<TH1F>("ENERGY_FULL_vm","ENERGY_D ON ALL CALO; energy [MeV]",500,0,10000);

  ENERGY_D_up =  tfs->make<TH1F>("ENERGY_D_up","ENERGY_D_up ON ALL CALO; energy [MeV]",500,0,10000);
  ENERGY_D_um =  tfs->make<TH1F>("ENERGY_D_um","ENERGY_D_um ON ALL CALO; energy [MeV]",500,0,10000);
  ENERGY_D_vp =  tfs->make<TH1F>("ENERGY_D_vp","ENERGY_D_vp ON ALL CALO; energy [MeV]",500,0,10000);
  ENERGY_D_vm =  tfs->make<TH1F>("ENERGY_D_vm","ENERGY_D_vm ON ALL CALO; energy [MeV]",500,0,10000);

  ENERGY_S1_up =  tfs->make<TH1F>("ENERGY_S1_up","ENERGY_S1_up ON ALL CALO; energy [MeV]",500,0,10000);
  ENERGY_S1_um =  tfs->make<TH1F>("ENERGY_S1_um","ENERGY_S1_um ON ALL CALO; energy [MeV]",500,0,10000);
  ENERGY_S1_vp =  tfs->make<TH1F>("ENERGY_S1_vp","ENERGY_S1_vp ON ALL CALO; energy [MeV]",500,0,10000);
  ENERGY_S1_vm =  tfs->make<TH1F>("ENERGY_S1_vm","ENERGY_S1_vm ON ALL CALO; energy [MeV]",500,0,10000);

  ENERGY_S2_up =  tfs->make<TH1F>("ENERGY_S2_up","ENERGY_S2_up ON ALL CALO; energy [MeV]",500,0,10000);
  ENERGY_S2_um =  tfs->make<TH1F>("ENERGY_S2_um","ENERGY_S2_um ON ALL CALO; energy [MeV]",500,0,10000);
  ENERGY_S2_vp =  tfs->make<TH1F>("ENERGY_S2_vp","ENERGY_S2_vp ON ALL CALO; energy [MeV]",500,0,10000);
  ENERGY_S2_vm =  tfs->make<TH1F>("ENERGY_S2_vm","ENERGY_S2_vm ON ALL CALO; energy [MeV]",500,0,10000);

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


   for(int i=0;i<24;i++){
     caloEnergies_All[i] = tfs->make<TH1F>(Form("EnergyOriginal_All%d",i),Form("Energy_Spectrum_calo_All_%d;energy [MeV]",i),500,0,10000);
     caloTimes_All[i] = tfs->make<TH1F>(Form("TimeOriginal_All_calo%d",i),Form("Time_calo_All%d;time [#mus];counts/145.15ns",i), nBins, 0, histMaxTime);
    
    caloEnergies_Up[i] = tfs->make<TH1F>(Form("EnergyOriginal_Up%d",i),Form("Energy_Spectrum_calo_Up_%d;energy [MeV]",i),500,0,10000);
    caloEnergies_Um[i] = tfs->make<TH1F>(Form("EnergyOriginal_Um%d",i),Form("Energy_Spectrum_calo_Um_%d;energy [MeV]",i),500,0,10000);
    caloEnergies_Vp[i] = tfs->make<TH1F>(Form("EnergyOriginal_Vp%d",i),Form("Energy_Spectrum_calo_Vp_%d;energy [MeV]",i),500,0,10000);
    caloEnergies_Vm[i] = tfs->make<TH1F>(Form("EnergyOriginal_Vm%d",i),Form("Energy_Spectrum_calo_Vm_%d;energy [MeV]",i),500,0,10000);

    caloTimes_Up[i] = tfs->make<TH1F>(Form("TimeOriginal_Up_calo%d",i),Form("Time_calo_%d;time [#mus];counts/147ns",i), nBins, 0, histMaxTime);
    caloTimes_Um[i] = tfs->make<TH1F>(Form("TimeOriginal_Um_calo%d",i),Form("Time_Um_calo_%d;time [#mus];counts/147ns",i), nBins, 0, histMaxTime);
    caloTimes_Vp[i] = tfs->make<TH1F>(Form("TimeOriginal_Vp_calo%d",i),Form("Time_Vp_calo_%d;time [#mus];counts/147ns",i), nBins, 0, histMaxTime);
    caloTimes_Vm[i] = tfs->make<TH1F>(Form("TimeOriginal_Vm_calo%d",i),Form("Time_Vm_calo_%d;time [#mus];counts/147ns",i), nBins, 0, histMaxTime);

    //deltaxvsdeltat[i] = tfs->make<TH2F>(Form("deltaxvsdeltat_%d",i),Form("deltax vs deltat_%d; time separation(ns); #delta X [c.w.]",i),200, 0, 5, 900, 0.5, 9.5);
    //deltayvsdeltat[i] = tfs->make<TH2F>(Form("deltayvsdeltat_%d",i),Form("deltay vs deltat_%d; time separation(ns); #delta Y [c.w.]",i),200, 0, 5, 600, 0.5, 6.5); 

    //deltat_calo[i] = tfs->make<TH1F>(Form("DeltaT_%d",i),Form("T_{E2}-T_{E1} calo %d;time [ns];",i),1000,0,40);


    energy_d_All[i] = tfs->make<TH1F>(Form("energy_D_All_%d",i),Form("energy_Double_calo_All_%d;energy [MeV]",i),500,0,10000);
    energy_s1_All[i] = tfs->make<TH1F>(Form("energy_S1_All_%d",i),Form("energy_Single_1_calo_All_%d;energy [MeV]",i),500,0,10000);
    energy_s2_All[i] = tfs->make<TH1F>(Form("energy_S2_All_%d",i),Form("energy_Single_2_calo_All_%d;energy [MeV]",i),500,0,10000);


    DTimes_calo_All[i] =  tfs->make<TH1F>(Form("time_D_All_calo_%d",i),Form("double_pileup_time_All_calo_%d;time [#mus];counts/147ns",i),nBins,0, histMaxTime);
    S1Times_calo_All[i] =  tfs->make<TH1F>(Form("time_S1_All_calo_%d",i),Form("single_1_pileup_time_All_calo_%d;time [#mus];counts/148ns",i),nBins,0, histMaxTime);
    S2Times_calo_All[i] =  tfs->make<TH1F>(Form("time_S2_All_calo_%d",i),Form("single_2_pileup_time_All_calo_%d;time [#mus];counts/148ns",i),nBins,0, histMaxTime);

    energy_d_Up[i] = tfs->make<TH1F>(Form("energy_D_Up_%d",i),Form("energy_Double_calo_Up_%d;energy [MeV]",i),500,0,10000);
    energy_s1_Up[i] = tfs->make<TH1F>(Form("energy_S1_Up_%d",i),Form("energy_Single_1_calo_Up_%d;energy [MeV]",i),500,0,10000);
    energy_s2_Up[i] = tfs->make<TH1F>(Form("energy_S2_Up_%d",i),Form("energy_Single_2_calo_Up_%d;energy [MeV]",i),500,0,10000);

    energy_d_Um[i] = tfs->make<TH1F>(Form("energy_D_Um_%d",i),Form("energy_Double_calo_Um_%d;energy [MeV]",i),500,0,10000);
    energy_s1_Um[i] = tfs->make<TH1F>(Form("energy_S1_Um_%d",i),Form("energy_Single_1_calo_Um_%d;energy [MeV]",i),500,0,10000);
    energy_s2_Um[i] = tfs->make<TH1F>(Form("energy_S2_Um_%d",i),Form("energy_Single_2_calo_Um_%d;energy [MeV]",i),500,0,10000);

    energy_d_Vp[i] = tfs->make<TH1F>(Form("energy_D_Vp_%d",i),Form("energy_Double_calo_Vp_%d;energy [MeV]",i),500,0,10000);
    energy_s1_Vp[i] = tfs->make<TH1F>(Form("energy_S1_Vp_%d",i),Form("energy_Single_1_calo_Vp_%d;energy [MeV]",i),500,0,10000);
    energy_s2_Vp[i] = tfs->make<TH1F>(Form("energy_S2_Vp_%d",i),Form("energy_Single_2_calo_Vp_%d;energy [MeV]",i),500,0,10000);

    energy_d_Vm[i] = tfs->make<TH1F>(Form("energy_D_Vm_%d",i),Form("energy_Double_calo_%d;energy [MeV]",i),500,0,10000);
    energy_s1_Vm[i] = tfs->make<TH1F>(Form("energy_S1_Vm_%d",i),Form("energy_Single_1_calo_%d;energy [MeV]",i),500,0,10000);
    energy_s2_Vm[i] = tfs->make<TH1F>(Form("energy_S2_Vm_%d",i),Form("energy_Single_2_calo_%d;energy [MeV]",i),500,0,10000);

    DTimes_calo_Up[i] =  tfs->make<TH1F>(Form("time_D_Up_calo_%d",i),Form("double_pileup_time_Up_calo_%d;time [#mus];counts/147ns",i),nBins,0, histMaxTime);
    DTimes_calo_Um[i] = tfs->make<TH1F>(Form("time_D_Um_calo_%d",i),Form("double pileup_time_Um_calo_%d;time [#mus];counts/147ns",i),nBins,0, histMaxTime);
    DTimes_calo_Vp[i] = tfs->make<TH1F>(Form("time_D_Vp_calo_%d",i),Form("double pileup time_Vp_calo_%d;time [#mus];counts/147ns",i),nBins,0, histMaxTime);
    DTimes_calo_Vm[i] = tfs->make<TH1F>(Form("time_D_Vm_calo_%d",i),Form("double pileup time_Vm_calo_%d;time [#mus];counts/148ns",i),nBins,0, histMaxTime);

    S1Times_calo_Up[i] =  tfs->make<TH1F>(Form("time_S1_Up_calo_%d",i),Form("single_1_pileup_time_Up_calo_%d;time [#mus];counts/148ns",i),nBins,0, histMaxTime);
    S1Times_calo_Um[i] = tfs->make<TH1F>(Form("time_S1_Um_calo_%d",i),Form("single_1_pileup_time_Um_calo_%d;time [#mus];counts/148ns",i),nBins,0, histMaxTime);
    S1Times_calo_Vp[i] = tfs->make<TH1F>(Form("time_S1_Vp_calo_%d",i),Form("single_1_pileup_time_Vp_calo_%d;time [#mus];counts/148ns",i),nBins,0, histMaxTime);
    S1Times_calo_Vm[i] = tfs->make<TH1F>(Form("time_S1_Vm_calo_%d",i),Form("single_1_pileup_time_Vm_calo_%d;time [#mus];counts/148ns",i),nBins,0, histMaxTime);

    S2Times_calo_Up[i] =  tfs->make<TH1F>(Form("time_S2_Up_calo_%d",i),Form("single_2_pileup_time_Up_calo_%d;time [#mus];counts/148ns",i),nBins,0, histMaxTime);
    S2Times_calo_Um[i] = tfs->make<TH1F>(Form("time_S2_Um_calo_%d",i),Form("single_2_pileup_time_Um_calo_%d;time [#mus];counts/148ns",i),nBins,0, histMaxTime);
    S2Times_calo_Vp[i] = tfs->make<TH1F>(Form("time_S2_Vp_calo_%d",i),Form("single_2_pileup_time_Vp_calo_%d;time [#mus];counts/148ns",i),nBins,0, histMaxTime);
    S2Times_calo_Vm[i] = tfs->make<TH1F>(Form("time_S2_Vm_calo_%d",i),Form("single_2_pileup_time_Vm_calo_%d;time [#mus];counts/148ns",i),nBins,0, histMaxTime);

    }
}



void  ratioeast::pileup_hist_up(double T1_up, double T2, double E1, double E2, TH1F *ENERGY_D_up, TH1F *ENERGY_S1_up, TH1F *ENERGY_S2_up, TH1F *DTimes_Up_all,TH1F *S1Times_Up_all, TH1F *S2Times_Up_all, TH1F **energy_d_Up, TH1F **energy_s1_Up,TH1F **energy_s2_Up,TH1F **DTimes_calo_Up,TH1F **S1Times_calo_Up,TH1F **S2Times_calo_Up){
 

  if(E1<Emin_ && E2<Emin_) { if(E1+E2>Emin_){
      ENERGY_D_up->Fill(E1+E2);
      DTimes_Up_all->Fill(T1_up+gapT_/2);

      for(unsigned int i=0;i<24;i++){
	energy_d_Up[i]->Fill(E1+E2);
	DTimes_calo_Up[i]->Fill(T1_up+gapT_/2);
      }
    }
  }
  
  if(E1>Emin_ && E2<Emin_) {
    ENERGY_D_up->Fill(E1+E2); 
    ENERGY_S1_up->Fill(E1);
    DTimes_Up_all->Fill(T1_up+gapT_/2);
    S1Times_Up_all->Fill(T1_up);
    
    
    for(unsigned int i=0;i<24;i++){
      energy_d_Up[i]->Fill(E1+E2);
      energy_s1_Up[i]->Fill(E1);   
      DTimes_calo_Up[i]->Fill(T1_up+gapT_/2);
      S1Times_calo_Up[i]->Fill(T1_up);
    }
  }


  if(E1<Emin_ && E2>Emin_) { 
    ENERGY_D_up->Fill(E1+E2);
    ENERGY_S2_up->Fill(E2);
    DTimes_Up_all->Fill(T1_up+gapT_/2);
    S2Times_Up_all->Fill(T1_up);

    for(unsigned int i=0;i<24;i++){
      energy_d_Up[i]->Fill(E1+E2);
      energy_s2_Up[i]->Fill(E2);
      DTimes_calo_Up[i]->Fill(T1_up+gapT_/2);
      S2Times_calo_Up[i]->Fill(T1_up);
    }
  }

  if(E1>Emin_ && E2>Emin_) {
    ENERGY_D_up->Fill(E1+E2); 
    ENERGY_S1_up->Fill(E1);
    ENERGY_S2_up->Fill(E2);
    DTimes_Up_all->Fill(T1_up+gapT_/2);
    S1Times_Up_all->Fill(T1_up);
    S2Times_Up_all->Fill(T1_up);

    for(unsigned int i=0;i<24;i++){
      energy_d_Up[i]->Fill(E1+E2);
      energy_s1_Up[i]->Fill(E1);
      energy_s2_Up[i]->Fill(E2);
      DTimes_calo_Up[i]->Fill(T1_up+gapT_/2);
      S1Times_calo_Up[i]->Fill(T1_up);
      S2Times_calo_Up[i]->Fill(T1_up);
    }
  }

}



void  ratioeast::pileup_hist_um(double T1_um, double T2, double E1, double E2, TH1F *ENERGY_D_um, TH1F *ENERGY_S1_um, TH1F *ENERGY_S2_um, TH1F *DTimes_Um_all,TH1F *S1Times_Um_all, TH1F *S2Times_Um_all, TH1F **energy_d_Um, TH1F **energy_s1_Um,TH1F **energy_s2_Um,TH1F **DTimes_calo_Um,TH1F **S1Times_calo_Um,TH1F **S2Times_calo_Um){



  if(E1<Emin_ && E2<Emin_) { if(E1+E2>Emin_){
      ENERGY_D_um->Fill(E1+E2);
      DTimes_Um_all->Fill(T1_um+gapT_/2);
      for(unsigned int i=0;i<24;i++){
        energy_d_Um[i]->Fill(E1+E2);
        DTimes_calo_Um[i]->Fill(T1_um+gapT_/2);
      }
    }
  }


  if(E1>Emin_ && E2<Emin_) {
    ENERGY_D_um->Fill(E1+E2);
    ENERGY_S1_um->Fill(E1);
    DTimes_Um_all->Fill(T1_um+gapT_/2);
    S1Times_Um_all->Fill(T1_um);


    for(unsigned int i=0;i<24;i++){
      energy_d_Um[i]->Fill(E1+E2);
      energy_s1_Um[i]->Fill(E1);
      DTimes_calo_Um[i]->Fill(T1_um+gapT_/2);
      S1Times_calo_Um[i]->Fill(T1_um);
    }
  }

  if(E1<Emin_ && E2>Emin_) {
    ENERGY_D_um->Fill(E1+E2);
    ENERGY_S2_um->Fill(E2);
    DTimes_Um_all->Fill(T1_um+gapT_/2);
    S2Times_Um_all->Fill(T1_um);

    for(unsigned int i=0;i<24;i++){
      energy_d_Um[i]->Fill(E1+E2);
      energy_s2_Um[i]->Fill(E2);
      DTimes_calo_Um[i]->Fill(T1_um+gapT_/2);
      S2Times_calo_Um[i]->Fill(T1_um);
    }
  }

  if(E1>Emin_ && E2>Emin_) {
    ENERGY_D_um->Fill(E1+E2);
    ENERGY_S1_um->Fill(E1);
    ENERGY_S2_um->Fill(E2);
    DTimes_Um_all->Fill(T1_um+gapT_/2);
    S1Times_Um_all->Fill(T1_um);
    S2Times_Um_all->Fill(T1_um);

    for(unsigned int i=0;i<24;i++){
      energy_d_Um[i]->Fill(E1+E2);
      energy_s1_Um[i]->Fill(E1);
      energy_s2_Um[i]->Fill(E2);
      DTimes_calo_Um[i]->Fill(T1_um+gapT_/2);
      S1Times_calo_Um[i]->Fill(T1_um);
      S2Times_calo_Um[i]->Fill(T1_um);
    }
  }
}


void  ratioeast::pileup_hist_vp(double T1_v1, double T2, double E1, double E2, TH1F *ENERGY_D_vp, TH1F *ENERGY_S1_vp, TH1F *ENERGY_S2_vp, TH1F *DTimes_Vp_all,TH1F *S1Times_Vp_all, TH1F *S2Times_Vp_all, TH1F **energy_d_Vp, TH1F **energy_s1_Vp,TH1F **energy_s2_Vp,TH1F **DTimes_calo_Vp,TH1F **S1Times_calo_Vp,TH1F **S2Times_calo_Vp){


  if(E1<Emin_ && E2<Emin_) { if(E1+E2>Emin_){
      ENERGY_D_vp->Fill(E1+E2);
      DTimes_Vp_all->Fill(T1_v1+gapT_/2);
      for(unsigned int i=0;i<24;i++){
        energy_d_Vp[i]->Fill(E1+E2);
        DTimes_calo_Vp[i]->Fill(T1_v1+gapT_/2);
      }
    }
  }

  if(E1>Emin_ && E2<Emin_) {
    ENERGY_D_vp->Fill(E1+E2);
    ENERGY_S1_vp->Fill(E1);
    DTimes_Vp_all->Fill(T1_v1+gapT_/2);
    S1Times_Vp_all->Fill(T1_v1);

    for(unsigned int i=0;i<24;i++){
      energy_d_Vp[i]->Fill(E1+E2);
      energy_s1_Vp[i]->Fill(E1);
      DTimes_calo_Vp[i]->Fill(T1_v1+gapT_/2);
      S1Times_calo_Vp[i]->Fill(T1_v1);
    }
  }

  if(E1<Emin_ && E2>Emin_) {
    ENERGY_D_vp->Fill(E1+E2);
    ENERGY_S2_vp->Fill(E2);
    DTimes_Vp_all->Fill(T1_v1+gapT_/2);
    S2Times_Vp_all->Fill(T1_v1);

    for(unsigned int i=0;i<24;i++){
      energy_d_Vp[i]->Fill(E1+E2);
      energy_s2_Vp[i]->Fill(E2);
      DTimes_calo_Vp[i]->Fill(T1_v1+gapT_/2);
      S2Times_calo_Vp[i]->Fill(T1_v1);
    }
  }

  if(E1>Emin_ && E2>Emin_) {
    ENERGY_D_vp->Fill(E1+E2);
    ENERGY_S1_vp->Fill(E1);
    ENERGY_S2_vp->Fill(E2);
    DTimes_Vp_all->Fill(T1_v1+gapT_/2);
    S1Times_Vp_all->Fill(T1_v1);
    S2Times_Vp_all->Fill(T1_v1);

    for(unsigned int i=0;i<24;i++){
      energy_d_Vp[i]->Fill(E1+E2);
      energy_s1_Vp[i]->Fill(E1);
      energy_s2_Vp[i]->Fill(E2);
      DTimes_calo_Vp[i]->Fill(T1_v1+gapT_/2);
      S1Times_calo_Vp[i]->Fill(T1_v1);
      S2Times_calo_Vp[i]->Fill(T1_v1);
    }
  }
}

void  ratioeast::pileup_hist_vm(double T1_v2, double T2, double E1, double E2, TH1F *ENERGY_D_vm, TH1F *ENERGY_S1_vm, TH1F *ENERGY_S2_vm, TH1F *DTimes_Vm_all,TH1F *S1Times_Vm_all, TH1F *S2Times_Vm_all, TH1F **energy_d_Vm, TH1F **energy_s1_Vm,TH1F **energy_s2_Vm,TH1F **DTimes_calo_Vm,TH1F **S1Times_calo_Vm,TH1F **S2Times_calo_Vm){

  if(E1<Emin_ && E2<Emin_) { if(E1+E2>Emin_){
      ENERGY_D_vm->Fill(E1+E2);
      DTimes_Vm_all->Fill(T1_v2+gapT_/2);

      for(unsigned int i=0;i<24;i++){
        energy_d_Vm[i]->Fill(E1+E2);
        DTimes_calo_Vm[i]->Fill(T1_v2+gapT_/2);
      }
    }
  }

  if(E1>Emin_ && E2<Emin_) {
    ENERGY_D_vm->Fill(E1+E2);
    ENERGY_S1_vm->Fill(E1);
    DTimes_Vm_all->Fill(T1_v2+gapT_/2);
    S1Times_Vm_all->Fill(T1_v2);


    for(unsigned int i=0;i<24;i++){
      energy_d_Vm[i]->Fill(E1+E2);
      energy_s1_Vm[i]->Fill(E1);
      DTimes_calo_Vm[i]->Fill(T1_v2+gapT_/2);
      S1Times_calo_Vm[i]->Fill(T1_v2);
    }
  }
  if(E1<Emin_ && E2>Emin_) {


    ENERGY_D_vm->Fill(E1+E2);
    ENERGY_S2_vm->Fill(E2);
    DTimes_Vm_all->Fill(T1_v2+gapT_/2);
    S2Times_Vm_all->Fill(T1_v2);

    for(unsigned int i=0;i<24;i++){

      energy_d_Vm[i]->Fill(E1+E2);
      energy_s2_Vm[i]->Fill(E2);
      DTimes_calo_Vm[i]->Fill(T1_v2+gapT_/2);
      S2Times_calo_Vm[i]->Fill(T1_v2);
    }
  }

  if(E1>Emin_ && E2>Emin_) {
    ENERGY_D_vm->Fill(E1+E2);
    ENERGY_S1_vm->Fill(E1);
    ENERGY_S2_vm->Fill(E2);
    DTimes_Vm_all->Fill(T1_v2+gapT_/2);
    S1Times_Vm_all->Fill(T1_v2);
    S2Times_Vm_all->Fill(T1_v2);

    for(unsigned int i=0;i<24;i++){
      energy_d_Vm[i]->Fill(E1+E2);
      energy_s1_Vm[i]->Fill(E1);
      energy_s2_Vm[i]->Fill(E2);
      DTimes_calo_Vm[i]->Fill(T1_v2+gapT_/2);
      S1Times_calo_Vm[i]->Fill(T1_v2);
      S2Times_calo_Vm[i]->Fill(T1_v2);
    }
  }
}




void ratioeast::analyze(art::Event const & e){

  mf::LogVerbatim("gm2reconeast::RatioEast") << "Enter gm2reconeast::Pileup for " << e.id();
  /*  unsigned int  eventnum = e.event(); 
  unsigned int runnum = e.run();  
  unsigned int subRunnum = e.subRun();    
  double fillId = runnum*1e6+subRunnum*1e3+eventnum;
  */
  //random1.SetSeed(c1_);                                                                       
  //random2.SetSeed(c2_);
 
  //midasSerialNum_ = gm2midastoart::getMidasEventHeader(e).eventNumber;


 

  //get sequence and trigger indices from fc7                                                          
  //const auto &encodercol =
  // *e.getValidHandle<gm2ccc::EncoderFC7ArtRecordCollection>(
  //   {"cccUnpacker", "unpacker"});
  // bunchNum_ = encodercol.size() ? encodercol.front().sequenceIndex : -1u;

  // runNum_ = e.run();
  // subRunNum_ = e.subRun();
  // eventNum_ = e.event();

  // get recon east data products                                                                                                         
  /* const auto &reconEastCollection = 
    *e.getValidHandle<gm2reconeast::GlobalFitArtRecordCollection>(
  	{readModule_, readInstance_});

   
    for(const auto &cluster : reconEastCollection){

      //if(cluster.time*1.25/1000 >= 29 && cluster.time*1.25/1000 <= 651){

      clusterTime_ = cluster.time*1.25/1000;

      //clusterEnergy_ = cluster.energy;

      //clusterEnergy_ = cluster.energy* 1700./ EnergyCal_run1[cluster.calorimeterIndex-1];

      clusterEnergy_ = (cluster.energy* 1700./ EnergyCal_run1[cluster.calorimeterIndex-1])*(1700./EnergyCal_run2[cluster.calorimeterIndex-1]);

      clusterY_ = cluster.position.second;
      clusterX_ = cluster.position.first;
      caloNum_ = cluster.calorimeterIndex;

      t_->Fill();

      //}
    }
   */
  
  art::Handle<gm2reconeast::GlobalFitArtRecordCollection> h;
  
  if(e.getByLabel({readModule_,readInstance_}, h)){
    gm2reconeast::GlobalFitArtRecordCollection clusters = *h;
    mf::LogVerbatim("gm2reconeast::RatioEast") << "Enter gm2reconeast::RatioEast for " << clusters.size();
    

    // get recon east data products                                                                                               
    
    std::vector<std::vector<double> > time(Ncalo);
    std::vector<std::vector<double> > energy(Ncalo);
    std::vector<std::vector<double> > island(Ncalo);
    std::vector<std::vector<double> > xposition(Ncalo);
    std::vector<std::vector<double> > yposition(Ncalo);
    std::vector<std::vector<double> > randombintime(Ncalo);
    std::vector<std::vector<double> > randomno(Ncalo);

    time.clear();
    energy.clear();
    island.clear();
    xposition.clear();
    yposition.clear();
    randombintime.clear();
    randomno.clear();

    std::vector<hitclass*> clusterhits;
    for(unsigned int c = 0; c < clusters.size(); ++c){
      auto aHit = new hitclass();
      aHit->time = clusters.at(c).time*1.25/1000;
      //aHit->time = ReconEast::CorrectTime( datasetName_reconEast_ , clusters.at(c).calorimeterIndex,  clusters.at(c).time, clusters.at(c).position.first,clusters.at(c).position.second); 
      //aHit->energy = clusters.at(c).energy;
      //aHit->energy = clusters.at(c).energy* 1700./ EnergyCal_run1[clusters.at(c).calorimeterIndex-1];
      aHit->energy = (clusters.at(c).energy* 1700./ EnergyCal_run1[clusters.at(c).calorimeterIndex-1])*(1700./EnergyCal_run2[clusters.at(c).calorimeterIndex-1]);   
      aHit->caloNum = clusters.at(c).calorimeterIndex;
      aHit->islandNum = clusters.at(c).islandIndex;
      aHit->xcluster = clusters.at(c).position.first;
      aHit->ycluster = clusters.at(c).position.second;
      clusterhits.push_back(aHit);
     }


    for(auto hit : clusterhits){
      for(int i=0;i<Ncalo;i++){
	if(i==hit->caloNum-1){

	  
	 double  randomBinTime = (random2.Uniform() * 0.14915) - (0.14915/2);
         double rand_time = hit->time+randomBinTime;
         double randomNum = random1.Uniform();
         double halfPeriod = g2Period/2;

	 // if(i==21){	  
	 //std::cout<<"i = "<<i<<" randomBinTime = "<<randomBinTime<<" randomNo= "<<randomNum<<std::endl;
	 // }

          randomno[i].push_back(randomNum);
          randombintime[i].push_back(randomBinTime);

	  if(hit->time>20 && hit->time<651 && hit->energy > 500){
	    ETspectrum->Fill(rand_time,hit->energy);

          }
	  if(hit->time>20 && hit->time<651 && hit->energy > Emin_){

            TIME_FULL->Fill(rand_time);
            ENERGY_FULL->Fill(hit->energy);


            caloEnergies_All[i]->Fill(hit->energy);
            caloTimes_All[i]->Fill(rand_time);


	    //std::cout<<"randomNum == "<<randomNum<<std::endl;

            if(randomNum < 0.25 && randomNum >= 0){
              ENERGY_FULL_up->Fill(hit->energy);
              caloEnergies_Up[i]->Fill(hit->energy);
              caloTimes_Up_all->Fill(rand_time - halfPeriod);
              caloTimes_Up[i]->Fill(rand_time - halfPeriod);

	      //std::cout<<" I am here 1"<<std::endl; 

            }

	    //std::cout<<"randomNum2 == "<<randomNum<<std::endl;
	    if(randomNum < 0.5 && randomNum >= 0.25){
              ENERGY_FULL_um->Fill(hit->energy);
              caloEnergies_Um[i]->Fill(hit->energy);
	      caloTimes_Um_all->Fill(rand_time + halfPeriod);
              caloTimes_Um[i]->Fill(rand_time + halfPeriod);
	      //std::cout<<" I am here 2"<<std::endl;
            }

	    //std::cout<<"randomNum3 == "<<randomNum<<std::endl;
            if(randomNum <0.75 && randomNum >= 0.50){
              ENERGY_FULL_vp->Fill(hit->energy);
              caloEnergies_Vp[i]->Fill(hit->energy);
              caloTimes_Vp_all->Fill(rand_time);
              caloTimes_Vp[i]->Fill(rand_time);
	      //std::cout<<" I am here 3"<<std::endl;
            }

	    //std::cout<<"randomNum4 == "<<randomNum<<std::endl;
            if(randomNum <1.0 && randomNum >= 0.75){
              ENERGY_FULL_vm->Fill(hit->energy);
              caloEnergies_Vm[i]->Fill(hit->energy);
              caloTimes_Vm_all->Fill(rand_time);
              caloTimes_Vm[i]->Fill(rand_time);
	      //std::cout<<" I am here 4"<<std::endl;
            }

            }//wiggle time cut     


          time[i].push_back(hit->time);
          energy[i].push_back(hit->energy);
          island[i].push_back(hit->islandNum);
          xposition[i].push_back(hit->xcluster);
	  yposition[i].push_back(hit->ycluster);
	}
      }
    }//clusterhit vector     
  
   
    
    std::vector<std::vector<double> > time_pile(Ncalo);
    std::vector<std::vector<double> > energy_pile(Ncalo);
    std::vector<std::vector<double> > randombintime_pile(Ncalo);
    std::vector<std::vector<double> > randomno_pile(Ncalo);


    time_pile.clear();
    energy_pile.clear();
    randombintime_pile.clear();
    randomno_pile.clear();

    for(int i=0;i<Ncalo;i++){
      //if(time[i].size()!=0){
      for(unsigned int j=0; j<time[i].size()-1; j++){
	if(time[i][j]>20 && time[i][j]<651){

	  if(time[i][j+1]-time[i][j]<DT_ && island[i][j+1]==island[i][j]){  
          En=energy[i][j+1]+energy[i][j];
          energy_pile[i].push_back(En);


	  Tn = (time[i][j+1]*energy[i][j+1]+time[i][j]*energy[i][j])/En;
	  time_pile[i].push_back(Tn);
	  randombintime_pile[i].push_back(randombintime[i][j]);
	  randomno_pile[i].push_back(randomno[i][j]);
       
          j++;
	 
	  }
	else{
          time_pile[i].push_back(time[i][j]);
          energy_pile[i].push_back(energy[i][j]);
	  randombintime_pile[i].push_back(randombintime[i][j]);
          randomno_pile[i].push_back(randomno[i][j]);
          if(j==time[i].size()-2){
            time_pile[i].push_back(time[i][j+1]);
            energy_pile[i].push_back(energy[i][j+1]);
	    randombintime_pile[i].push_back(randombintime[i][j+1]);
	    randomno_pile[i].push_back(randomno[i][j+1]);
	    
          }
	}
	  }//time
	
      }//clusters
      //}
    }//calos
    
        
      
    for(int i=0;i<Ncalo;i++){
      for(unsigned int c=0;c<time_pile[i].size();c++){
             
	
    /*if(energy_pile[i][c] > 500 && energy_pile[i][c+1] > 500 && time_pile[i][c+1] > 29 && time_pile[i][c+1] < 651 && time_pile[i][c] > 29 && time_pile[i][c] < 651){
	  deltat->Fill((time_pile[i][c+1]-time_pile[i][c])*1000);
	  if(energy_pile[i][c] > 1700 && energy_pile[i][c+1] > 1700){
	  deltaxvsdeltatallcalo->Fill((time_pile[i][c+1]-time_pile[i][c])*1000,fabs(xposition[i][c+1]-xposition[i][c]));
	  deltayvsdeltatallcalo->Fill((time_pile[i][c+1]-time_pile[i][c])*1000,fabs(yposition[i][c+1]-yposition[i][c])); 
	  deltaxvsdeltat[i]->Fill((time_pile[i][c+1]-time_pile[i][c])*1000,fabs(xposition[i][c+1]-xposition[i][c]));
	  deltayvsdeltat[i]->Fill((time_pile[i][c+1]-time_pile[i][c])*1000,fabs(yposition[i][c+1]-yposition[i][c]));
	}
	}*/

    
    if(time_pile[i][c]>20 && time_pile[i][c]<651 && fabs(xposition[i][c+1]-xposition[i][c]) <DX_ && fabs(yposition[i][c+1]-yposition[i][c]) <DY_ && c<time_pile[i].size()-2){

	  if(time_pile[i][c+1]-time_pile[i][c]<gapT_+DT_ && time_pile[i][c+1]-time_pile[i][c]>gapT_ && time_pile[i][c+2]-time_pile[i][c]>2*(DT_+gapT_) ) {

	    E1 = energy_pile[i][c];
            E2 = energy_pile[i][c+1];
	   
	    T1 = time_pile[i][c];	    
            T2 = time_pile[i][c+1];
	    Td = T1+gapT_/2;

	    //full T-method pileup time histograms 
	    if(E1 > Emin_ && E2 > Emin_){TIME_D->Fill(Td+randombintime_pile[i][c]);TIME_S1->Fill(T1+randombintime_pile[i][c]);TIME_S2->Fill(T1+randombintime_pile[i][c]);ENERGY_D->Fill(E1+E2); ENERGY_S1->Fill(E1); ENERGY_S2->Fill(E2); energy_d_All[i]->Fill(E1+E2); energy_s1_All[i]->Fill(E1); energy_s2_All[i]->Fill(E2); DTimes_calo_All[i]->Fill(Td+randombintime_pile[i][c]); S1Times_calo_All[i]->Fill(T1+randombintime_pile[i][c]); S2Times_calo_All[i]->Fill(T1+randombintime_pile[i][c]);} if(E1 > Ethre_ && E2 > Ethre_) {ET_D->Fill(Td+randombintime_pile[i][c], E1+E2); ET_S1->Fill(T1+randombintime_pile[i][c], E1);ET_S2->Fill(T1+randombintime_pile[i][c], E2);}

	    if(E2 < Emin_ && E1 > Emin_){TIME_D->Fill(Td+randombintime_pile[i][c]);TIME_S1->Fill(T1+randombintime_pile[i][c]);ENERGY_D->Fill(E1+E2); ENERGY_S1->Fill(E1); energy_d_All[i]->Fill(E1+E2); energy_s1_All[i]->Fill(E1); DTimes_calo_All[i]->Fill(Td+randombintime_pile[i][c]); S1Times_calo_All[i]->Fill(T1+randombintime_pile[i][c]);} if(E2 < Ethre_ && E1 > Ethre_){ ET_D->Fill(Td+randombintime_pile[i][c], E1+E2); ET_S1->Fill(T1+randombintime_pile[i][c], E1);}

	    if(E2 > Emin_ && E1 < Emin_){TIME_D->Fill(Td+randombintime_pile[i][c]);TIME_S2->Fill(T1+randombintime_pile[i][c]);ENERGY_D->Fill(E1+E2); ENERGY_S2->Fill(E2); energy_d_All[i]->Fill(E1+E2); energy_s2_All[i]->Fill(E2); DTimes_calo_All[i]->Fill(Td+randombintime_pile[i][c]); S2Times_calo_All[i]->Fill(T1+randombintime_pile[i][c]);}  if(E2 > Ethre_ && E1 < Ethre_){ ET_D->Fill(Td+randombintime_pile[i][c], E1+E2); ET_S2->Fill(T1+randombintime_pile[i][c], E2);}

	    if(E1 < Emin_ && E2 < Emin_){if((E1+E2)>Emin_){TIME_D->Fill(Td+randombintime_pile[i][c]);ENERGY_D->Fill(E1+E2); energy_d_All[i]->Fill(E1+E2); DTimes_calo_All[i]->Fill(Td+randombintime_pile[i][c]);}}if(E1 < Ethre_ && E2 < Ethre_){if((E1+E2)>Ethre_){ ET_D->Fill(Td+randombintime_pile[i][c], E1+E2);}}
	    
	    //pileup histograms for 4-spectra to create ratio	    
	    if(randomno_pile[i][c] < 0.25 && randomno_pile[i][c] > 0){
	      T1_up = (T1+randombintime_pile[i][c])-g2Period/2;
	      pileup_hist_up( T1_up, T2, E1, E2, ENERGY_D_up, ENERGY_S1_up, ENERGY_S2_up, DTimes_Up_all, S1Times_Up_all, S2Times_Up_all, energy_d_Up, energy_s1_Up, energy_s2_Up,DTimes_calo_Up,S1Times_calo_Up,S2Times_calo_Up);	   
	    }

	    if(randomno_pile[i][c] < 0.50 && randomno_pile[i][c] > 0.25){
	      T1_um = (T1+randombintime_pile[i][c])+g2Period/2;
              pileup_hist_um( T1_um, T2, E1, E2, ENERGY_D_um, ENERGY_S1_um, ENERGY_S2_um, DTimes_Um_all, S1Times_Um_all, S2Times_Um_all, energy_d_Um, energy_s1_Um, energy_s2_Um,DTimes_calo_Um,S1Times_calo_Um,S2Times_calo_Um);

            }


            if(randomno_pile[i][c] < 0.75 && randomno_pile[i][c] > 0.50){
	      T1_v1 = T1+randombintime_pile[i][c];
              pileup_hist_vp( T1_v1, T2, E1, E2, ENERGY_D_vp, ENERGY_S1_vp, ENERGY_S2_vp, DTimes_Vp_all, S1Times_Vp_all, S2Times_Vp_all, energy_d_Vp, energy_s1_Vp, energy_s2_Vp,DTimes_calo_Vp,S1Times_calo_Vp,S2Times_calo_Vp);
	    }


	    if(randomno_pile[i][c] < 1.00 && randomno_pile[i][c] > 0.75){
	      T1_v2 = T1+randombintime_pile[i][c];
              pileup_hist_vm( T1_v2, T2, E1, E2, ENERGY_D_vm, ENERGY_S1_vm, ENERGY_S2_vm, DTimes_Vm_all, S1Times_Vm_all, S2Times_Vm_all, energy_d_Vm, energy_s1_Vm, energy_s2_Vm,DTimes_calo_Vm,S1Times_calo_Vm,S2Times_calo_Vm);
	    }
	  }//category 1

	  
	  if(time_pile[i][c+1]-time_pile[i][c]<gapT_+DT_ && time_pile[i][c+1]-time_pile[i][c]>gapT_ && time_pile[i][c+2]-time_pile[i][c+1]<gapT_) {
	    E1 = energy_pile[i][c];
            E2 = energy_pile[i][c+1];
            T1 = time_pile[i][c];
            T2 = time_pile[i][c+1];
	    Td = T1+gapT_/2;

	    //full T-method pileup time histograms                                                      
	    if(E1 > Emin_ && E2 > Emin_){TIME_D->Fill(Td+randombintime_pile[i][c]);TIME_S1->Fill(T1+randombintime_pile[i][c]);TIME_S2->Fill(T1+randombintime_pile[i][c]);ENERGY_D->Fill(E1+E2); ENERGY_S1->Fill(E1); ENERGY_S2->Fill(E2);  energy_d_All[i]->Fill(E1+E2); energy_s1_All[i]->Fill(E1); energy_s2_All[i]->Fill(E2); DTimes_calo_All[i]->Fill(Td+randombintime_pile[i][c]); S1Times_calo_All[i]->Fill(T1+randombintime_pile[i][c]); S2Times_calo_All[i]->Fill(T1+randombintime_pile[i][c]);} if(E1 > Ethre_ && E2 > Ethre_){ET_D->Fill(Td+randombintime_pile[i][c], E1+E2); ET_S1->Fill(T1+randombintime_pile[i][c], E1); ET_S2->Fill(T1+randombintime_pile[i][c], E2);}

            if(E2 < Emin_ && E1 > Emin_){TIME_D->Fill(Td+randombintime_pile[i][c]);TIME_S1->Fill(T1+randombintime_pile[i][c]);ENERGY_D->Fill(E1+E2); ENERGY_S1->Fill(E1); energy_d_All[i]->Fill(E1+E2); energy_s1_All[i]->Fill(E1); DTimes_calo_All[i]->Fill(Td+randombintime_pile[i][c]); S1Times_calo_All[i]->Fill(T1+randombintime_pile[i][c]);} if(E2 < Ethre_ && E1 > Ethre_){ET_D->Fill(Td+randombintime_pile[i][c], E1+E2); ET_S1->Fill(T1+randombintime_pile[i][c], E1);}

            if(E2 > Emin_ && E1 < Emin_){TIME_D->Fill(Td+randombintime_pile[i][c]);TIME_S2->Fill(T1+randombintime_pile[i][c]);ENERGY_D->Fill(E1+E2); ENERGY_S2->Fill(E2); energy_d_All[i]->Fill(E1+E2); energy_s2_All[i]->Fill(E2); DTimes_calo_All[i]->Fill(Td+randombintime_pile[i][c]); S2Times_calo_All[i]->Fill(T1+randombintime_pile[i][c]);} if(E2 > Ethre_ && E1 < Ethre_){ET_D->Fill(Td+randombintime_pile[i][c], E1+E2); ET_S2->Fill(T1+randombintime_pile[i][c], E2);}

            if(E1 < Emin_ && E2 < Emin_){if((E1+E2)>Emin_){TIME_D->Fill(Td+randombintime_pile[i][c]);ENERGY_D->Fill(E1+E2);  energy_d_All[i]->Fill(E1+E2); DTimes_calo_All[i]->Fill(Td+randombintime_pile[i][c]);}}if(E1 < Ethre_ && E2 < Ethre_){if((E1+E2)>Ethre_){ET_D->Fill(Td+randombintime_pile[i][c], E1+E2);}}
	    
	    //pileup histograms for 4-spectra to create ratio

	    if(randomno_pile[i][c] < 0.25 && randomno_pile[i][c] > 0){
	      T1_up = (T1+randombintime_pile[i][c])-g2Period/2;
              pileup_hist_up( T1_up, T2, E1, E2, ENERGY_D_up, ENERGY_S1_up, ENERGY_S2_up, DTimes_Up_all, S1Times_Up_all, S2Times_Up_all, energy_d_Up, energy_s1_Up, energy_s2_Up,DTimes_calo_Up,S1Times_calo_Up,S2Times_calo_Up);

            }


            if(randomno_pile[i][c] < 0.50 && randomno_pile[i][c] > 0.25){
	      T1_um = (T1+randombintime_pile[i][c])+g2Period/2;
              pileup_hist_um( T1_um, T2, E1, E2, ENERGY_D_um, ENERGY_S1_um, ENERGY_S2_um, DTimes_Um_all, S1Times_Um_all, S2Times_Um_all, energy_d_Um, energy_s1_Um, energy_s2_Um,DTimes_calo_Um,S1Times_calo_Um,S2Times_calo_Um);

            }


            if(randomno_pile[i][c] < 0.75 && randomno_pile[i][c] > 0.50){
	      T1_v1 = T1+randombintime_pile[i][c];
              pileup_hist_vp( T1_v1, T2, E1, E2, ENERGY_D_vp, ENERGY_S1_vp, ENERGY_S2_vp, DTimes_Vp_all, S1Times_Vp_all, S2Times_Vp_all, energy_d_Vp, energy_s1_Vp, energy_s2_Vp,DTimes_calo_Vp,S1Times_calo_Vp,S2Times_calo_Vp);
            }


            if(randomno_pile[i][c] < 1.00 && randomno_pile[i][c] > 0.75){
	      T1_v2 = T1+randombintime_pile[i][c];
              pileup_hist_vm( T1_v2, T2, E1, E2, ENERGY_D_vm, ENERGY_S1_vm, ENERGY_S2_vm, DTimes_Vm_all, S1Times_Vm_all, S2Times_Vm_all, energy_d_Vm, energy_s1_Vm, energy_s2_Vm,DTimes_calo_Vm,S1Times_calo_Vm,S2Times_calo_Vm);
            }
	  }//category2	  


	  if(time_pile[i][c+1]-time_pile[i][c]<gapT_+DT_ && time_pile[i][c+1]-time_pile[i][c]>gapT_ && time_pile[i][c+2]-time_pile[i][c]<2*(DT_+gapT_) && time_pile[i][c+2]-time_pile[i][c+1]>gapT_) {

	    E1 = energy_pile[i][c];
            E2 = energy_pile[i][c+1];
            T1 = time_pile[i][c];
	    T2 = time_pile[i][c+1];
	    Td = T1+gapT_/2;

	    //full T-method pileup time histograms                                                      
	    if(E1 > Emin_ && E2 > Emin_){TIME_D->Fill(Td+randombintime_pile[i][c]);TIME_S1->Fill(T1+randombintime_pile[i][c]);TIME_S2->Fill(T1+randombintime_pile[i][c]); ENERGY_D->Fill(E1+E2); ENERGY_S1->Fill(E1); ENERGY_S2->Fill(E2); energy_d_All[i]->Fill(E1+E2); energy_s1_All[i]->Fill(E1); energy_s2_All[i]->Fill(E2); DTimes_calo_All[i]->Fill(Td+randombintime_pile[i][c]); S1Times_calo_All[i]->Fill(T1+randombintime_pile[i][c]); S2Times_calo_All[i]->Fill(T1+randombintime_pile[i][c]);}  if(E1 > Ethre_ && E2 > Ethre_){ET_D->Fill(Td+randombintime_pile[i][c], E1+E2); ET_S1->Fill(T1+randombintime_pile[i][c], E1); ET_S2->Fill(T1+randombintime_pile[i][c], E2);}

            if(E2 < Emin_ && E1 > Emin_){TIME_D->Fill(Td+randombintime_pile[i][c]);TIME_S1->Fill(T1+randombintime_pile[i][c]);ENERGY_D->Fill(E1+E2); ENERGY_S1->Fill(E1);  energy_d_All[i]->Fill(E1+E2); energy_s1_All[i]->Fill(E1); DTimes_calo_All[i]->Fill(Td+randombintime_pile[i][c]); S1Times_calo_All[i]->Fill(T1+randombintime_pile[i][c]);} if(E2 < Ethre_ && E1 > Ethre_){ET_D->Fill(Td+randombintime_pile[i][c], E1+E2); ET_S1->Fill(T1+randombintime_pile[i][c], E1);}

            if(E2 > Emin_ && E1 < Emin_){TIME_D->Fill(Td+randombintime_pile[i][c]);TIME_S2->Fill(T1+randombintime_pile[i][c]);ENERGY_D->Fill(E1+E2); ENERGY_S2->Fill(E2); energy_d_All[i]->Fill(E1+E2); energy_s2_All[i]->Fill(E2); DTimes_calo_All[i]->Fill(Td+randombintime_pile[i][c]); S2Times_calo_All[i]->Fill(T1+randombintime_pile[i][c]);} if(E2 > Ethre_ && E1 < Ethre_){ET_D->Fill(Td+randombintime_pile[i][c], E1+E2); ET_S2->Fill(T1+randombintime_pile[i][c], E2);}

            if(E1 < Emin_ && E2 < Emin_){if((E1+E2)>Emin_){TIME_D->Fill(Td+randombintime_pile[i][c]);ENERGY_D->Fill(E1+E2);  energy_d_All[i]->Fill(E1+E2); DTimes_calo_All[i]->Fill(Td+randombintime_pile[i][c]);}}if(E1 < Ethre_ && E2 < Ethre_){if((E1+E2)>Ethre_){ET_D->Fill(Td+randombintime_pile[i][c], E1+E2);}}
	    

	    //pileup histograms for 4-spectra to create ratio
	    if(randomno_pile[i][c] < 0.25 && randomno_pile[i][c] > 0){
	      T1_up = (T1+randombintime_pile[i][c])-g2Period/2;  
              pileup_hist_up( T1_up, T2, E1, E2, ENERGY_D_up, ENERGY_S1_up, ENERGY_S2_up, DTimes_Up_all, S1Times_Up_all, S2Times_Up_all, energy_d_Up, energy_s1_Up, energy_s2_Up,DTimes_calo_Up,S1Times_calo_Up,S2Times_calo_Up);

            }

            if(randomno_pile[i][c] < 0.50 && randomno_pile[i][c] > 0.25){
	      T1_um = (T1+randombintime_pile[i][c])+g2Period/2; 
              pileup_hist_um( T1_um, T2, E1, E2, ENERGY_D_um, ENERGY_S1_um, ENERGY_S2_um, DTimes_Um_all, S1Times_Um_all, S2Times_Um_all, energy_d_Um, energy_s1_Um, energy_s2_Um,DTimes_calo_Um,S1Times_calo_Um,S2Times_calo_Um);

            }

            if(randomno_pile[i][c] < 0.75 && randomno_pile[i][c] > 0.50){
	      T1_v1= T1+randombintime_pile[i][c];
              pileup_hist_vp( T1_v1, T2, E1, E2, ENERGY_D_vp, ENERGY_S1_vp, ENERGY_S2_vp, DTimes_Vp_all, S1Times_Vp_all, S2Times_Vp_all, energy_d_Vp, energy_s1_Vp, energy_s2_Vp,DTimes_calo_Vp,S1Times_calo_Vp,S2Times_calo_Vp);
            }

            if(randomno_pile[i][c] < 1.00 && randomno_pile[i][c] > 0.75){
	      T1_v2= T1+randombintime_pile[i][c];
              pileup_hist_vm( T1_v2, T2, E1, E2, ENERGY_D_vm, ENERGY_S1_vm, ENERGY_S2_vm, DTimes_Vm_all, S1Times_Vm_all, S2Times_Vm_all, energy_d_Vm, energy_s1_Vm, energy_s2_Vm,DTimes_calo_Vm,S1Times_calo_Vm,S2Times_calo_Vm);
	    }
	      
	  }//category3
	  
	  
	}
	
    }//pileup vector
    }//time vector 
} //calo 
            
  
}//analyze
DEFINE_ART_MODULE(ratioeast)
