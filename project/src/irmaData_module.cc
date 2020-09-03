////////////////////////////////////////////////////////////////////////
// Class:       irmaData
// Plugin Type: analyzer (art v2_10_03)
// File:        irmaData_module.cc
//
// Generated at Wed Sep  2 16:48:27 2020 by root using cetskelgen
// from cetlib version v3_02_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

class irmaData;


class irmaData : public art::EDAnalyzer {
public:
  explicit irmaData(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  irmaData(irmaData const &) = delete;
  irmaData(irmaData &&) = delete;
  irmaData & operator = (irmaData const &) = delete;
  irmaData & operator = (irmaData &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

private:

  // Declare member data here.

};


irmaData::irmaData(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{}

void irmaData::analyze(art::Event const & e)
{
  LOG_INFO("HI") << "Bonjour!";
}

DEFINE_ART_MODULE(irmaData)
