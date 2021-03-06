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
#include "gm2dataproducts/daq/CCCArtRecord.hh"

// HDF5 includes
#include "hep_hpc/hdf5/make_ntuple.hpp"

using hep_hpc::hdf5::Column;
using hep_hpc::hdf5::make_scalar_column;
using hep_hpc::hdf5::PropertyList;

const int chunkSize = 1024 * 1024 / sizeof(int);

// A convenience function to make columns in the Ntuple with the correct chunksize and compression
template <typename T>
inline
hep_hpc::hdf5::Column<T, 1ull>
make_scalar_column_with_my_params(std::string name) {
  return make_scalar_column<T>(name, chunkSize, {PropertyList{H5P_DATASET_CREATE}(&H5Pset_shuffle)(&H5Pset_deflate,6u)} );
}

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

  // module and instance for fc7
  std::string fc7Module_;
  std::string fc7Instance_;

  // The output HDF5 file
  hep_hpc::hdf5::File h5File_;

  // Cluster ntuple
  hep_hpc::hdf5::Ntuple<Column<int, 1>,   // Run #
                        Column<int, 1>,   // SubRun #
                        Column<int, 1>,   // Event #
                        Column<int, 1>,   // Calorimeter index
                        Column<int, 1>,   // Island index
                        Column<int, 1>,   // Bunch number
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
      fc7Module_(p.get<std::string>("fc7Module", "cccUnpacker")),
      fc7Instance_(p.get<std::string>("fc7Instance", "unpacker")),
      h5File_(fileName_, H5F_ACC_TRUNC),
      clusterTuple_(hep_hpc::hdf5::make_ntuple({h5File_, "ReconEastClusters", 1000},
                                               make_scalar_column_with_my_params<int>("run"),
                                               make_scalar_column_with_my_params<int>("subrun"),
                                               make_scalar_column_with_my_params<int>("event"),
                                               make_scalar_column_with_my_params<int>("caloIndex"),
                                               make_scalar_column_with_my_params<int>("islandIndex"),
                                               make_scalar_column_with_my_params<int>("bunchNum"),
                                               make_scalar_column_with_my_params<float>("time"),
                                               make_scalar_column_with_my_params<float>("energy"),
                                               make_scalar_column_with_my_params<float>("x"),
                                               make_scalar_column_with_my_params<float>("y")))
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

  // Get the bunch number
  art::Handle<gm2ccc::EncoderFC7ArtRecordCollection> fc7H;
  e.getByLabel({fc7Module_, fc7Instance_}, fc7H);
  const auto &fc7 = *fc7H;
  int bunchNum = fc7.size() ? fc7.front().sequenceIndex : -1;

  // Get the calorimter information
  art::Handle<gm2reconeast::CaloGlobalFitViewArtRecordCollection> caloViews;
  e.getByLabel({readModule_, readInstance_}, caloViews);

  // Loop over the calorimeters with clusters
  for (auto &clustersThisCaloPtrVtr : *caloViews)
  {

    // Loop over clusters
    for (auto &aClusterPtr : clustersThisCaloPtrVtr.clusters)
    {
      const auto &aCluster = *(aClusterPtr.get());

      // Fill the ntuple
      clusterTuple_.insert(
          run, subrun, event,
          aCluster.calorimeterIndex,
          aCluster.islandIndex,
          bunchNum,
          aCluster.time,
          aCluster.energy,
          aCluster.position.first,
          aCluster.position.second);
    }
  }
}

DEFINE_ART_MODULE(IrmaData)
