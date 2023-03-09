/***********************************************/
/**
* @file gnssParametrizationLeoDynamicOrbits.cpp
*
* @brief Orbits by variational equations.
* @see GnssParametrization
*
* @author Torsten Mayer-Guerr
* @date 2021-01-23
*
*/
/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "files/fileInstrument.h"
#include "files/fileMatrix.h"
#include "files/fileParameterName.h"
#include "gnss/gnss.h"
#include "classes/ephemerides/ephemerides.h"
#include "classes/parametrizationAcceleration/parametrizationAcceleration.h"
#include "classes/timeSeries/timeSeries.h"
#include "classes/timeSeries/timeSeries.h"
#include "classes/platformSelector/platformSelector.h"
#include "misc/observation/variationalEquationFromFile.h"
#include "gnss/gnssParametrization/gnssParametrizationLeoDynamicOrbits.h"

/***********************************************/

GnssParametrizationLeoDynamicOrbits::GnssParametrizationLeoDynamicOrbits(Config &config)
{
  try
  {
    TimeSeriesPtr stochasticPulse;

    readConfig(config, "name",                        name,                        Config::OPTIONAL, "parameter.dynamicOrbits", "used for parameter selection");
    readConfig(config, "selectReceivers",             selectReceivers,             Config::MUSTSET,  "",     "");
    readConfig(config, "outputfileOrbit",             fileNameOrbit,               Config::OPTIONAL, "",     "variable {station} available");
    readConfig(config, "outputfileParameters",        fileNameParameter,           Config::OPTIONAL, "",     "variable {station} available");
    readConfig(config, "inputfileVariational",        fileNameVariational,         Config::MUSTSET,  "variational_{loopTime:%D}.{station}.dat", "variable {station} available");
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

GnssParametrizationLeoDynamicOrbits::~GnssParametrizationLeoDynamicOrbits()
{
  for(Parameter *para : parameters)
    delete para;
}

/***********************************************/

void GnssParametrizationLeoDynamicOrbits::init(Gnss *gnss, Parallel::CommunicatorPtr comm)
{
  try
  {
    this->gnss = gnss;
    auto selectedReceivers = gnss->selectReceivers(selectReceivers);

    VariableList fileNameVariableList;
    addVariable("station", "****", fileNameVariableList);
    parameters.resize(gnss->receivers.size(), nullptr);
    Vector recvProcess(gnss->receivers.size());
    for(UInt idRecv=0; idRecv<gnss->receivers.size(); idRecv++)
      if(selectedReceivers.at(idRecv) && gnss->receivers.at(idRecv)->useable())
      {
        recvProcess(idRecv) = Parallel::myRank(comm)+1;
        auto para = new Parameter();
        parameters.at(idRecv) = para;
        para->recv = gnss->receivers.at(idRecv);

        if(!para->recv->isMyRank())
          continue;

        // find first and last valid epoch
        Time timeStart, timeEnd;
        for(UInt idEpoch=0; idEpoch<gnss->times.size(); idEpoch++)
          if(para->recv->useable(idEpoch))
          {
            timeStart = gnss->times.at(idEpoch);
            break;
          }
        for(UInt idEpoch=gnss->times.size(); idEpoch-->0;)
          if(para->recv->useable(idEpoch))
          {
            timeEnd = gnss->times.at(idEpoch);
            break;
          }

        // stochastic pulses in interval
        std::vector<Time> pulsesInterval;
        std::copy_if(pulses.begin(), pulses.end(), std::back_inserter(pulsesInterval), [&](auto &p) {return (timeStart < p) && (p < timeEnd);});

        fileNameVariableList["station"]->setValue(para->recv->name());
        VariationalEquationFromFile file;
        file.open(fileNameVariational(fileNameVariableList), nullptr/*parametrizationGravity*/, parametrizationAcceleration, pulsesInterval, ephemerides, integrationDegree);

        auto variationalEquation = file.integrateArc(timeStart, timeEnd, TRUE/*computePosition*/, TRUE/*computeVelocity*/);
        para->times     = variationalEquation.times;
        para->PosDesign = variationalEquation.PosDesign;
        para->VelDesign = variationalEquation.VelDesign;
        para->pos       = variationalEquation.pos0;
        para->vel       = variationalEquation.vel0;
        para->x         = Vector(para->PosDesign.columns());
        para->polynomial.init(para->times, interpolationDegree);

        std::vector<Time> times;
        for(UInt idEpoch=0; idEpoch<gnss->times.size(); idEpoch++)
          if(para->recv->useable(idEpoch))
            times.push_back(para->recv->timeCorrected(idEpoch));

        const Vector pos = para->polynomial.interpolate(times, para->pos, 3);
        const Vector vel = para->polynomial.interpolate(times, para->vel, 3);
        UInt i=0;
        for(UInt idEpoch=0; idEpoch<gnss->times.size(); idEpoch++)
          if(para->recv->useable(idEpoch))
          {
            para->recv->pos.at(idEpoch) = Vector3d(pos.row(3*i, 3));
            para->recv->vel.at(idEpoch) = Vector3d(vel.row(3*i, 3));
            i++;
          }

        // parameter names
        file.parameterNameSatellite(para->parameterNames);
        file.parameterNameSatelliteArc(para->parameterNames);
        for(auto &name : para->parameterNames)
          name.object = para->recv->name();
      }

    // distribute process id of receivers
    Parallel::reduceSum(recvProcess, 0, comm);
    Parallel::broadCast(recvProcess, 0, comm);

    // synchronize parameter names
    for(UInt idRecv=0; idRecv<gnss->receivers.size(); idRecv++)
      if(recvProcess(idRecv))
        Parallel::broadCast(parameters.at(idRecv)->parameterNames, static_cast<UInt>(recvProcess(idRecv)-1), comm);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationLeoDynamicOrbits::requirements(GnssNormalEquationInfo &normalEquationInfo,
                                                       std::vector<UInt> &/*transCount*/, std::vector<UInt> &/*transCountEpoch*/,
                                                       std::vector<UInt> &recvCount,  std::vector<UInt> &/*recvCountEpoch*/)
{
  try
  {
    if(isEnabled(normalEquationInfo, name))
      for(auto para : parameters)
        if(para && para->recv->isMyRank())
         recvCount.at(para->recv->idRecv()) += static_cast<UInt>(minEstimableEpochsRatio * normalEquationInfo.idEpochs.size());
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationLeoDynamicOrbits::initParameter(GnssNormalEquationInfo &normalEquationInfo)
{
  try
  {
   for(auto para : parameters)
      if(para)
        para->index = GnssParameterIndex();
    if(!isEnabled(normalEquationInfo, name))
      return;

    UInt countPara = 0;
    for(auto para : parameters)
      if(para && para->recv->useable() && normalEquationInfo.estimateReceiver.at(para->recv->idRecv()))
      {
        para->index = normalEquationInfo.parameterNamesReceiver(para->recv->idRecv(), para->parameterNames);
        countPara += para->parameterNames.size();
      }
    if(countPara)
      logInfo<<countPara%"%9i dynamic orbit parameters"s<<Log::endl;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationLeoDynamicOrbits::aprioriParameter(const GnssNormalEquationInfo &normalEquationInfo, MatrixSliceRef x0) const
{
  try
  {
    for(auto para : parameters)
      if(para && para->index && para->recv->isMyRank())
        copy(para->x, x0.row(normalEquationInfo.index(para->index), para->x.rows()));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationLeoDynamicOrbits::designMatrix(const GnssNormalEquationInfo &/*normalEquationInfo*/, const GnssObservationEquation &eqn, GnssDesignMatrix &A) const
{
  try
  {
    auto para = parameters.at(eqn.receiver->idRecv());
    if(para && para->index)
      matMult(1., eqn.A.column(GnssObservationEquation::idxPosRecv, 3),
              para->polynomial.interpolate({eqn.timeRecv}, para->PosDesign, 3),
              A.column(para->index));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double GnssParametrizationLeoDynamicOrbits::updateParameter(const GnssNormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef /*Wz*/)
{
  try
  {
    Double maxChange = 0;
    Gnss::InfoParameterChange info("mm");
    for(auto para : parameters)
      if(para && para->index && para->recv->isMyRank())
      {
        const Vector dx = x.row(normalEquationInfo.index(para->index), para->x.rows());
        para->x   += dx;
        para->pos += para->PosDesign * dx;
        para->vel += para->VelDesign * dx;

        std::vector<Time> times;
        for(UInt idEpoch=0; idEpoch<gnss->times.size(); idEpoch++)
          if(para->recv->useable(idEpoch))
            times.push_back(para->recv->timeCorrected(idEpoch));

        const Vector pos = para->polynomial.interpolate(times, para->pos, 3);
        const Vector vel = para->polynomial.interpolate(times, para->vel, 3);
        UInt i=0;
        for(UInt idEpoch=0; idEpoch<gnss->times.size(); idEpoch++)
          if(para->recv->useable(idEpoch))
          {
            const Vector3d dpos = Vector3d(pos.row(3*i, 3)) - para->recv->pos.at(idEpoch);
            para->recv->pos.at(idEpoch) = Vector3d(pos.row(3*i, 3));
            para->recv->vel.at(idEpoch) = Vector3d(vel.row(3*i, 3));
            if(info.update(1e3*dpos.r()))
              info.info = "position receiver ("+para->recv->name()+", "+times.at(i).dateTimeStr()+")";
            i++;
          }
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

void GnssParametrizationLeoDynamicOrbits::writeResults(const GnssNormalEquationInfo &normalEquationInfo, const std::string &suffix) const
{
  try
  {
    if(!isEnabled(normalEquationInfo, name))
      return;

    if(!fileNameOrbit.empty())
    {
      VariableList fileNameVariableList;
      addVariable("station", "***", fileNameVariableList);
      logStatus<<"write receiver orbits to files <"<<fileNameOrbit(fileNameVariableList).appendBaseName(suffix)<<">"<<Log::endl;
      for(auto para : parameters)
        if(para && para->index && para->recv->isMyRank())
        {

          OrbitArc arc;
          for(UInt idEpoch=0; idEpoch<para->times.size(); idEpoch++)
          {
            OrbitEpoch epoch;
            epoch.time     = para->times.at(idEpoch);
            epoch.position = Vector3d(para->pos.row(3*idEpoch, 3));
            epoch.velocity = Vector3d(para->vel.row(3*idEpoch, 3));
            arc.push_back(epoch);
          }
          fileNameVariableList["station"]->setValue(para->recv->name());
          InstrumentFile::write(fileNameOrbit(fileNameVariableList).appendBaseName(suffix), arc);
        }
    }

    if(!fileNameParameter.empty())
    {
      VariableList fileNameVariableList;
      addVariable("station", "***", fileNameVariableList);
      logStatus<<"write estimated receiver parameters to files <"<<fileNameParameter(fileNameVariableList).appendBaseName(suffix)<<">"<<Log::endl;
      for(auto para : parameters)
        if(para && para->index && para->recv->isMyRank())
        {
          fileNameVariableList["station"]->setValue(para->recv->name());
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
