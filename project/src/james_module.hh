
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

using std::cerr;
using std::cout;
using std::endl;
using std::string;
//using std::unique_ptr;
using gm2reconeast::CaloGlobalFitViewArtRecord;
using gm2reconeast::CaloGlobalFitViewArtRecordCollection;
using gm2reconeast::GlobalFitArtRecord;
using gm2reconeast::GlobalFitArtRecordCollection;
using mf::LogDebug;
using mf::LogInfo;
using mf::LogVerbatim;
using std::vector;

class ClusterHit;
class PileupHit;

//int randSeed1 = 0;
//int randSeed2 = 0;
TRandom3 random1; // seeds set in ratioeast constructor
TRandom3 random2;
double tBinWidth; // will be set to ratioeast.tBinWidth_

const double Fa = 0.2290735 * 1e6 * 1 / (1e9); // TODO: extract from Blinders library
const double g2Period = (1. / Fa) / 1000.;
const double g2HalfPeriod = g2Period / 2.;
const double ratioTimeShiftDefinitions[] = {
    // easy convention indexed by 0,1,2,3
    -g2HalfPeriod, // U+
    g2HalfPeriod,  // U-
    0.,            // V1
    0.             // V2
};
const double muonLifetime = 64400.;    // ns
static const unsigned int NCalos = 24; // NOTE: zero-indexed (0 through 23)

const double EnergyCal_run1[24] = {
    1628.9, 1505.9, 1559.4, 1564.9, 1368.8, 1516.9, 1543.8, 1533.0,
    1518.1, 1551.7, 1582.6, 1610.8, 1604.2, 1566.5, 1528.0, 1487.0,
    1520.0, 1588.2, 1554.9, 1525.5, 1455.7, 1474.9, 1522.5, 1548.1};

const double EnergyCal_run2[24] = {
    1845.34, 1956.21, 1852.62, 1882.91, 2075.44, 1919.23, 1885.50, 1900.13,
    1893.89, 1889.55, 1880.18, 1906.99, 1920.32, 1910.87, 1913.59, 1963.29,
    1973.80, 1931.09, 1932.92, 1943.34, 1976.45, 1964.30, 1928.10, 1933.15};

/*
  simple implementation: do conversions & calibrations immediately in the
  constructor, save every piece of the sum(s), and define 'getter' 
  functions which explicitly name the adjustments that have been done
*/
class PileupHit
{
public:
  ClusterHit *sources[5];
  unsigned int NSources;
  // There are multiple types of 'average time' (and for all we know there
  // will be multiple types of 'average energy' at some point).  And the
  // data-product/cluster hits already store all of the time/energy
  // information.  So instead of storing information here and having a
  // 'state watcher' (bool recomputed) for newly-added cluster hits (and
  // having to store multiple time values for the different types of times
  // in the cluster hit) let's just recompute it every time we need it
  // (which seems like will always be less than 10 recomputations or so).
  //double time;
  //double energy;
  //bool recomputed;

  void init();
  PileupHit() { init(); }
  // TODO: variadic template functions (someday)
  PileupHit(ClusterHit *hit)
  {
    init();
    addSource(hit);
  }
  PileupHit(ClusterHit *hit1, ClusterHit *hit2)
  {
    init();
    addSource(hit1);
    addSource(hit2);
  }
  PileupHit(ClusterHit *hit1, ClusterHit *hit2, ClusterHit *hit3)
  {
    init();
    addSource(hit1);
    addSource(hit2);
    addSource(hit3);
  }
  void addSource(ClusterHit *hit);
  void recompute();
  //double getTime();
  double getAvgAnalysisTime();
  double getAvgRandomizedAnalysisTime();
  double getEnergySum();
  double getAvgEnergy();
  double getQuarteringConstant(); // definition may change
};
class ClusterHit
{
public:
  const GlobalFitArtRecord *gfRecord;
  unsigned int NShadows;
  PileupHit *pileupHit;

  // I would just point to (or reference) the value in the art record to
  // save space, but that gets complicted fast.  Reference members make a
  // classes default constructor 'ineligible', so we can't allocate arrays
  // of them.  Pointers would work but xCluster and yCluster are in a STL
  // pair<double,double>, and we shouldn't raw-point to something in an STL
  // container.
  double rawTime;
  double rawEnergy;
  unsigned int iCalo;
  unsigned int islandNum;
  double xCluster;
  double yCluster;

  // values after unit conversion & calibration (no randomization or ratio shift)
  double analysisTime;
  double analysisEnergy;

  // values used to quarter dataset and shift times for FR
  double quarteringConstant; // [0,4) so just do int(quarteringConstant) for an index [0,3]
  double ratioTimeShift;     // from ratioTimeShiftDefinition[]
  double randomTimeShift;    // uniform between +-0.5*tBinWidth

  ClusterHit(){}; // default constructible (but this gets deleted if the class
  // owns any references)
  ClusterHit(const GlobalFitArtRecord *globalFitArtRecord) : gfRecord(globalFitArtRecord),
                                                             //NShadows(0),
                                                             pileupHit(0),
                                                             rawTime(globalFitArtRecord->time),
                                                             rawEnergy(globalFitArtRecord->energy),
                                                             iCalo(globalFitArtRecord->calorimeterIndex - 1), // art record has 1-24, we're using 0-23
                                                             islandNum(globalFitArtRecord->islandIndex),
                                                             xCluster(globalFitArtRecord->position.first),
                                                             yCluster(globalFitArtRecord->position.second),
                                                             analysisTime(0.),
                                                             analysisEnergy(0.),
                                                             quarteringConstant(4. * random1.Uniform()),
                                                             ratioTimeShift(ratioTimeShiftDefinitions[int(quarteringConstant)]),
                                                             randomTimeShift(tBinWidth * (random2.Uniform() - 0.5))
  {

    // convert [1.25ns] to [us] and apply time calib (if there is one)
    analysisTime = rawTime * 1.25 / 1000. + 0.;

    // no unit conversion needed, but there *is* an energy calibration
    //analysisEnergy = rawEnergy * 1700. / energyCalibRun2[iCalo];
    analysisEnergy = rawEnergy * 1700. / EnergyCal_run1[iCalo] * 1700. / EnergyCal_run2[iCalo];
  };
  //void addShadows(ClusterHit * shadowHit) {}

  // values after unit conversion & calibration:
  void setPileupHit(PileupHit *pHit);

  // time values after unit conversion & calibration
  double getAnalysisTime() const { return analysisTime; }

  // values after unit conversion, calibration, & randomization:
  double getRandomizedAnalysisTime() const { return getAnalysisTime() + randomTimeShift; }

  // values after unit conversion, calibration, randomization, & ratio time shift:
  double getRandomizedAndShiftedAnalysisTime() const
  {
    return getRandomizedAnalysisTime() + ratioTimeShift;
  }

  // energy in MeV after per-calo calibration
  double getAnalysisEnergy() const { return analysisEnergy; }

  // quick description mostly for debugging
  string summary() const
  {
    char tmpstr[128];
    sprintf(tmpstr, "ta=%.2f Ea=%.2f c=%d i=%d r1=%.6f/4 r2=%.6f/binwidth[us]",
            analysisTime, analysisEnergy,
            iCalo, islandNum,
            quarteringConstant, randomTimeShift);
    return string(tmpstr);
  }
};

// TODO: don't have two separate functions for hits & pileup hits, it's just goofy
// TODO: well actually they might need separate criteria...
bool hitPassesCaloXYCuts(ClusterHit *Hit)
{
  return true; // TODO: implement (flexibly? fcl param?)
}
bool hitPassesCaloXYCuts(PileupHit *puHit)
{
  return true; // TODO: implement (flexibly? fcl param?)
}

bool hitPassesTimeCuts(ClusterHit *hit, double tmin, double tmax, double margin = 0.)
{
  /*
    e.g. margin=g2HalfPeriod+1.5*binWidth_
  */
  return tmin - margin < hit->getAnalysisTime() && hit->getAnalysisTime() < tmax + margin;
}
bool hitPassesTimeCuts(PileupHit *puHit, double tmin, double tmax, double margin = 0.)
{
  return tmin - margin < puHit->getAvgAnalysisTime() && puHit->getAvgAnalysisTime() < tmax + margin;
}

bool shadowCoincidenceMatch(ClusterHit *h1, ClusterHit *h2, double deadTime)
{
  //cerr<<"t1:"<<h1->getAnalysisTime()<<"  t2:"<<h2->getAnalysisTime()<<"\n";
  return h2->getAnalysisTime() - h1->getAnalysisTime() < deadTime && h1->islandNum == h2->islandNum;
}

bool hitPassesEnergyCuts(ClusterHit *hit, double EMin)
{
  return hit->getAnalysisEnergy() >= EMin;
}
