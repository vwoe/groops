/***********************************************/
/**
* @file orbitPropagatorEuler.h
*
* @brief Propagate a dynamic orbit using Euler's method.
* @see orbitPropagator
*
* @author Matthias Ellmer
* @date 2017-01-19
*
*/
/***********************************************/

#ifndef __GROOPS_ORBITPROPAGATOREULER__
#define __GROOPS_ORBITPROPAGATOREULER__

// Latex documentation
#ifdef DOCSTRING_OrbitPropagator
static const char *docstringOrbitPropagatorEuler = R"(
\subsection{Euler}
This class implements Euler's method to propagate a satellite orbit under the influence of \configClass{Forces}{forcesType}.
Satellite is assumed to be oriented along-track.
)";
#endif

/***********************************************/

#include "classes/orbitPropagator/orbitPropagator.h"

/***** CLASS ***********************************/

/** @brief Propagate orbit using Euler's method.
* @ingroup orbitPropagatorGroup
* @see orbitPropagator */
class OrbitPropagatorEuler : public OrbitPropagator
{
public:
  OrbitPropagatorEuler(Config &/*config*/) {}

  OrbitArc integrateArc(const OrbitEpoch &startEpoch, const Time &sampling, UInt posCount, ForcesPtr forces, SatelliteModelPtr satellite, EarthRotationPtr earthRotation,
                        EphemeridesPtr ephemerides, Bool timing) const override;
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline OrbitArc OrbitPropagatorEuler::integrateArc(const OrbitEpoch &startEpoch, const Time &sampling, UInt posCount, ForcesPtr forces,
                                                   SatelliteModelPtr satellite, EarthRotationPtr earthRotation, EphemeridesPtr ephemerides, Bool timing) const
{
  try
  {
    OrbitArc orbit;
    orbit.push_back(startEpoch);
    orbit.back().acceleration = acceleration(startEpoch, forces, satellite, earthRotation, ephemerides);
    const Double dt = sampling.seconds();

    Single::forEach(posCount-1, [&](UInt k)
    {
      OrbitEpoch epoch;
      epoch.time         = startEpoch.time + (k+1) * sampling;
      epoch.position     = orbit.at(k).position + dt * orbit.at(k).velocity;
      epoch.velocity     = orbit.at(k).velocity + dt * orbit.at(k).acceleration;
      epoch.acceleration = acceleration(epoch, forces, satellite, earthRotation, ephemerides);
      orbit.push_back(epoch);
    }, timing);

    return orbit;
  }
  catch (std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif /* __GROOPS_ORBITPROPAGATOREULER__ */
