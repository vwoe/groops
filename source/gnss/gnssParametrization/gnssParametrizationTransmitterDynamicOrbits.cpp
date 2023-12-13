/***********************************************/
/**
* @file gnssParametrizationTransmitterDynamicOrbits.cpp
*
* @brief Orbits by variational equations.
* @see GnssParametrization
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2021-01-23
*
*/
/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "files/fileInstrument.h"
#include "files/fileMatrix.h"
#include "classes/ephemerides/ephemerides.h"
#include "classes/parametrizationAcceleration/parametrizationAcceleration.h"
#include "classes/timeSeries/timeSeries.h"
#include "classes/timeSeries/timeSeries.h"
#include "classes/platformSelector/platformSelector.h"
#include "misc/observation/variationalEquationFromFile.h"
#include "gnss/gnss.h"
#include "gnss/gnssParametrization/gnssParametrizationTransmitterDynamicOrbits.h"

/***********************************************/

GnssParametrizationTransmitterDynamicOrbits::GnssParametrizationTransmitterDynamicOrbits(Config &config)
{
  try
  {
    TimeSeriesPtr stochasticPulse;

    readConfig(config, "name",                        name,                        Config::OPTIONAL, "parameter.transmitterDynamicOrbits", "used for parameter selection");
    readConfig(config, "selectTransmitters",          selectTransmitters,          Config::MUSTSET,  "",     "");
    readConfig(config, "outputfileOrbit",             fileNameOrbit,               Config::OPTIONAL, "",     "variable {prn} available");
    readConfig(config, "outputfileParameters",        fileNameParameter,           Config::OPTIONAL, "",     "variable {prn} available");
    readConfig(config, "inputfileVariational",        fileNameVariational,         Config::MUSTSET,  "variational_{loopTime:%D}.{prn}.dat", "variable {prn} available");
    readConfig(config, "stochasticPulse",             stochasticPulse,             Config::DEFAULT,  "",     "[mu/s] parametrization of stochastic pulses");
    readConfig(config, "parametrizationAcceleration", parametrizationAcceleration, Config::DEFAULT,  "",     "orbit force parameters");
    readConfig(config, "ephemerides",                 ephemerides,                 Config::MUSTSET,  "",     "");
    readConfig(config, "minEstimableEpochsRatio",     minEstimableEpochsRatio,     Config::DEFAULT,  "0.75", "drop satellites with lower ratio of estimable epochs to total epochs");
    readConfig(config, "integrationDegree",           integrationDegree,           Config::DEFAULT,  "7",    "integration of forces by polynomial approximation of degree n");
    readConfig(config, "interpolationDegree",         interpolationDegree,         Config::DEFAULT,  "7",    "for orbit interpolation and velocity calculation");
    if(isCreateSchema(config)) return;

    pulses = stochasticPulse->times();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

GnssParametrizationTransmitterDynamicOrbits::~GnssParametrizationTransmitterDynamicOrbits()
{
  for(Parameter *para : parameters)
    delete para;
}

/***********************************************/

void GnssParametrizationTransmitterDynamicOrbits::init(Gnss *gnss, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    this->gnss = gnss;
    auto selectedTransmitters = gnss->selectTransmitters(selectTransmitters);

    VariableList fileNameVariableList;
    fileNameVariableList.setVariable("prn", "***");
    parameters.resize(gnss->transmitters.size(), nullptr);
    for(UInt idTrans=0; idTrans<gnss->transmitters.size(); idTrans++)
      if(selectedTransmitters.at(idTrans) && gnss->transmitters.at(idTrans)->useable())
      {
        auto para = new Parameter();
        parameters.at(idTrans) = para;
        para->trans = gnss->transmitters.at(idTrans);

        // find first and last valid epoch
        const Time timeStart = para->trans->timesPosVel.front();
        const Time timeEnd   = para->trans->timesPosVel.back();

        // stochastic pulses in interval
        std::vector<Time> pulsesInterval;
        std::copy_if(pulses.begin(), pulses.end(), std::back_inserter(pulsesInterval), [&](auto &p) {return (timeStart < p) && (p < timeEnd);});

        fileNameVariableList.setVariable("prn", para->trans->name());
        VariationalEquationFromFile file;
        file.open(fileNameVariational(fileNameVariableList), nullptr/*parametrizationGravity*/, parametrizationAcceleration, pulsesInterval, ephemerides, integrationDegree);

        auto variationalEquation = file.integrateArc(timeStart, timeEnd, TRUE/*computePosition*/, TRUE/*computeVelocity*/);
        para->times     = variationalEquation.times;
        para->PosDesign = variationalEquation.PosDesign;
        para->VelDesign = variationalEquation.VelDesign;
        para->x         = Vector(para->PosDesign.columns());
        para->polynomial.init(para->times, interpolationDegree);

        // parameter names
        file.parameterNameSatellite(para->parameterNames);
        file.parameterNameSatelliteArc(para->parameterNames);
        for(auto &name : para->parameterNames)
          name.object = para->trans->name();
      }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationTransmitterDynamicOrbits::requirements(GnssNormalEquationInfo &normalEquationInfo,
                                                               std::vector<UInt> &transCount, std::vector<UInt> &/*transCountEpoch*/,
                                                               std::vector<UInt> &/*recvCount*/,  std::vector<UInt> &/*recvCountEpoch*/)
{
  try
  {
    if(isEnabled(normalEquationInfo, name) && !normalEquationInfo.isEachReceiverSeparately)
      for(auto para : parameters)
        if(para && para->trans->useable())
         transCount.at(para->trans->idTrans()) += static_cast<UInt>(minEstimableEpochsRatio * normalEquationInfo.idEpochs.size());
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationTransmitterDynamicOrbits::initParameter(GnssNormalEquationInfo &normalEquationInfo)
{
  try
  {
   for(auto para : parameters)
      if(para)
        para->index = GnssParameterIndex();
    if(!isEnabled(normalEquationInfo, name) || normalEquationInfo.isEachReceiverSeparately)
      return;

    UInt countPara = 0;
    for(auto para : parameters)
      if(para && para->trans->useable())
      {
        para->index = normalEquationInfo.parameterNamesTransmitter(para->trans->idTrans(), para->parameterNames);
        countPara += para->parameterNames.size();
      }
    if(countPara)
      logInfo<<countPara%"%9i transmitter dynamic orbit parameters"s<<Log::endl;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationTransmitterDynamicOrbits::aprioriParameter(const GnssNormalEquationInfo &normalEquationInfo, MatrixSliceRef x0) const
{
  try
  {
    if(Parallel::isMaster(normalEquationInfo.comm))
      for(auto para : parameters)
        if(para && para->index)
          copy(para->x, x0.row(normalEquationInfo.index(para->index), para->x.rows()));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationTransmitterDynamicOrbits::designMatrix(const GnssNormalEquationInfo &/*normalEquationInfo*/, const GnssObservationEquation &eqn, GnssDesignMatrix &A) const
{
  try
  {
    auto para = parameters.at(eqn.transmitter->idTrans());
    if(para && para->index)
      matMult(1., eqn.A.column(GnssObservationEquation::idxPosTrans, 3),
              para->polynomial.interpolate({eqn.timeTrans}, para->PosDesign, 3),
              A.column(para->index));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double GnssParametrizationTransmitterDynamicOrbits::updateParameter(const GnssNormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef /*Wz*/)
{
  try
  {
    Double maxChange = 0;
    Gnss::InfoParameterChange info("mm");
    for(auto para : parameters)
      if(para && para->index)
      {
        const Vector dx = x.row(normalEquationInfo.index(para->index), para->x.rows());
        para->x += dx;
        const Vector dpos = para->polynomial.interpolate(para->trans->timesPosVel, para->PosDesign * dx, 3);
        const Vector dvel = para->polynomial.interpolate(para->trans->timesPosVel, para->VelDesign * dx, 3);
        para->trans->pos += reshape(dpos, 3, para->trans->pos.rows()).trans();
        para->trans->vel += reshape(dvel, 3, para->trans->vel.rows()).trans();

        for(UInt i=0; i<para->trans->timesPosVel.size(); i++)
          if(info.update(1e3*norm(dpos.row(3*i, 3))))
            info.info = "position transmitter ("+para->trans->name()+", "+para->trans->timesPosVel.at(i).dateTimeStr()+")";
      }
    info.synchronizeAndPrint(normalEquationInfo.comm, 1e-3, maxChange);
    return maxChange;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationTransmitterDynamicOrbits::writeResults(const GnssNormalEquationInfo &normalEquationInfo, const std::string &suffix) const
{
  try
  {
    if(!isEnabled(normalEquationInfo, name) || !Parallel::isMaster(normalEquationInfo.comm))
      return;

    if(!fileNameOrbit.empty())
    {
      VariableList fileNameVariableList;
      fileNameVariableList.setVariable("prn", "***");
      logStatus<<"write transmitter orbits to files <"<<fileNameOrbit(fileNameVariableList).appendBaseName(suffix)<<">"<<Log::endl;
      for(auto para : parameters)
        if(para && para->index)
        {
          OrbitArc arc;
          for(UInt idEpoch=0; idEpoch<para->times.size(); idEpoch++)
          {
            OrbitEpoch epoch;
            epoch.time     = para->times.at(idEpoch);
            epoch.position = para->trans->positionCoM(para->times.at(idEpoch));
            epoch.velocity = para->trans->velocity(para->times.at(idEpoch));
            arc.push_back(epoch);
          }
          fileNameVariableList.setVariable("prn", para->trans->name());
          InstrumentFile::write(fileNameOrbit(fileNameVariableList).appendBaseName(suffix), arc);
        }
    }

    if(!fileNameParameter.empty())
    {
      VariableList fileNameVariableList;
      fileNameVariableList.setVariable("prn", "***");
      logStatus<<"write estimated transmitter parameters to files <"<<fileNameParameter(fileNameVariableList).appendBaseName(suffix)<<">"<<Log::endl;
      for(auto para : parameters)
        if(para && para->index)
        {
          fileNameVariableList.setVariable("prn", para->trans->name());
          writeFileMatrix(fileNameParameter(fileNameVariableList).appendBaseName(suffix), para->x);
        }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
