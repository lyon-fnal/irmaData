/**
 * @file IrmaData_module.cc
 * @author Adam Lyon (lyon@fnal.gov)
 * @brief 
 * 
 * See https://redmine.fnal.gov/redmine/projects/dunetpc/repository/revisions/13e0e8d239b9ea7701f4416594a96e21b28d6eb7/entry/dune/CVN/art/GCNH5_module.cc for an example of HDF5 in art using hep_hpc
 * 
 * @date 2020-09-09
 */

// Standard art includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Specifc art includes
#include "gm2dataproducts/reconeast/GlobalFitArtRecord.hh"

// HDF5 includes
#include "hep_hpc/hdf5/make_ntuple.hpp"

// Boost includes (for the file name)
#include <boost/uuid/uuid.hpp>            // uuid class
#include <boost/uuid/uuid_generators.hpp> // generators
#include <boost/uuid/uuid_io.hpp>         // streaming

#include <sstream>
#include <iomanip>

using hep_hpc::hdf5::Column;
using hep_hpc::hdf5::make_scalar_column;

class IrmaData;

/**
 * @brief Analyzer to convert IRMA data from art to HDF5
 * 
 */
class IrmaData : public art::EDAnalyzer
{
public:
  explicit IrmaData(fhicl::ParameterSet const &p);
  ~IrmaData() noexcept {};

  // Plugins should not be copied or assigned.
  IrmaData(IrmaData const &) = delete;
  IrmaData(IrmaData &&) = delete;
  IrmaData &operator=(IrmaData const &) = delete;
  IrmaData &operator=(IrmaData &&) = delete;

  // Required functions.
  void analyze(art::Event const &e) override;

private:
  // The output file name
  std::string fileName_;

  // The read module
  std::string readModule_;

  // The read instance
  std::string readInstance_;

  // The output HDF5 file
  hep_hpc::hdf5::File h5File_;

  // Cluster ntuple
  hep_hpc::hdf5::Ntuple<Column<int, 1>,   // Run #
                        Column<int, 1>,   // SubRun #
                        Column<int, 1>,   // Event #
                        Column<int, 1>,   // Calorimeter index
                        Column<int, 1>,   // Island index
                        Column<float, 1>, // time
                        Column<float, 1>, // energy
                        Column<float, 1>, // x
                        Column<float, 1>> // y
      clusterTuple_;                      
};

IrmaData::IrmaData(fhicl::ParameterSet const &p)
    : EDAnalyzer(p),
      fileName_(p.get<std::string>("fileName")),
      readModule_(p.get<std::string>("readModule", "energyPartition")),
      readInstance_(p.get<std::string>("readInstance", "partition")),
      h5File_(fileName_, H5F_ACC_TRUNC),
      clusterTuple_(hep_hpc::hdf5::make_ntuple({h5File_, "ReconEastClusters", 1000},
                                               make_scalar_column<int>("run"),
                                               make_scalar_column<int>("subrun"),
                                               make_scalar_column<int>("event"),
                                               make_scalar_column<int>("caloIndex"),
                                               make_scalar_column<int>("islandIndex"),
                                               make_scalar_column<float>("time"),
                                               make_scalar_column<float>("energy"),
                                               make_scalar_column<float>("x"),
                                               make_scalar_column<float>("y")))
{
}

/**
 * @brief Fill the HDF5 Ntuple for this event
 * 
 * @param e The event from art
 */
void IrmaData::analyze(art::Event const &e)
{
  int run = e.id().run();
  int subrun = e.id().subRun();
  int event = e.id().event();

  // Get the calorimter information
  art::Handle<gm2reconeast::CaloGlobalFitViewArtRecordCollection> caloViews;
  e.getByLabel({readModule_, readInstance_}, caloViews);

  // Loop over the calorimeters with clusters
  for (auto &clustersThisCaloPtrVtr : *caloViews)
  {

    // Loop over clusters
    for (auto &aClusterPtr : clustersThisCaloPtrVtr.clusters)
    {

      const auto &aCluster = aClusterPtr.get();

      // Fill the ntuple
      clusterTuple_.insert(
          run, subrun, event,
          aCluster->calorimeterIndex,
          aCluster->islandIndex,
          aCluster->time,
          aCluster->energy,
          aCluster->position.first,
          aCluster->position.second);
    }
  }
}

DEFINE_ART_MODULE(IrmaData)
