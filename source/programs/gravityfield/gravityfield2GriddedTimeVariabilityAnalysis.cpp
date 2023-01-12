/***********************************************/
/**
* @file gravityfield2GriddedTimeVariabilityAnalysis.cpp
*
* @brief Analyse time variability of gravity field with respect to reference signal.
*
* @author Torsten Mayer-Guerr
* @date 2009-05-19
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This programm computes a time series of functionals of a \class{gravityfield}{gravityfieldType}
on each point at a given \class{grid}{gridType}. The type of functional (e.g gravity anomalies or geoid heights)
can be choosen with \class{kernel}{kernelType}. At each grid point the error root mean square (RMS),
the correlation coefficient or the Nash-Sutcliffe coefficient (NSC) will be computed frome the time series
with respect to the \config{referencefield}. The functional will be saved together with points expressed
as ellipsoidal coordinates (longitude, latitude, height) based on a reference ellipsoid with parameters \config{R}
and \config{inverseFlattening}. To speed up the computation the gravity field will be converted to
spherical harmonics before the computation.

To visualize the results use \program{PlotMap}.

This programm is parallelized.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileGriddedData.h"
#include "classes/grid/grid.h"
#include "classes/timeSeries/timeSeries.h"
#include "classes/kernel/kernel.h"
#include "classes/gravityfield/gravityfield.h"
#include "misc/miscGriddedData.h"

/***** CLASS ***********************************/

/** @brief Analyse time variability of gravity field with respect to reference signal.
* @ingroup programGroup */
class Gravityfield2GriddedTimeVariabilityAnalysis
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(Gravityfield2GriddedTimeVariabilityAnalysis, PARALLEL, "Analyse time variability at grid points", Gravityfield, Grid)

/***********************************************/

void Gravityfield2GriddedTimeVariabilityAnalysis::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    enum functionalType {CORRELATION, ERRORRMS, NSC, PERCENTAGE};
    std::string     choice;
    FileName        outName;
    GravityfieldPtr gravityfield,referencefield;
    KernelPtr       kernel;
    GridPtr         grid;
    TimeSeriesPtr   timeSeries;
    Double          a, f;
    functionalType functional;

    readConfig(config, "outputfileGriddedData",  outName,      Config::MUSTSET,  "", "");
    readConfig(config, "grid",               grid,         Config::MUSTSET, "", "");
    readConfig(config, "kernel",             kernel,       Config::MUSTSET, "", "");
    readConfig(config, "gravityfield",       gravityfield, Config::MUSTSET, "", "");
    readConfig(config, "referencefield",     referencefield, Config::MUSTSET, "", "gravity field to compare to");

    readConfigChoice(config, "functional", choice, Config::MUSTSET, "", "correlation coefficient with respect to reference field");
    if(readConfigChoiceElement(config, "correlationCoefficient", choice, ""))                                                           functional = CORRELATION;
    if(readConfigChoiceElement(config, "errorRMS",               choice, "error rms with respect to reference field"))                  functional = ERRORRMS;
    if(readConfigChoiceElement(config, "NSC",                    choice, "Nash-Sutcliffe-Coefficient with respect to reference field")) functional = NSC;
    if(readConfigChoiceElement(config, "percentage",             choice, "Percentage of reference field explained by signal"))          functional = PERCENTAGE;
    endChoice(config);

    readConfig(config, "timeSeries",         timeSeries,   Config::MUSTSET, "", "");
    readConfig(config, "R",                  a,            Config::OPTIONAL, STRING_DEFAULT_GRS80_a, "reference radius for ellipsoidal coordinates on output");
    readConfig(config, "inverseFlattening",  f,            Config::OPTIONAL, STRING_DEFAULT_GRS80_f, "reference flattening for ellipsoidal coordinates on output, 0: spherical coordinates");
    if(isCreateSchema(config)) return;

    std::vector<Time>     times  = timeSeries->times();
    std::vector<Vector3d> points = grid->points();
    std::vector<Double>   areas  = grid->areas();

    // Create gridded data
    // -------------------
    logStatus<<"create gravity functionals on grid"<<Log::endl;
    UInt count = times.size();
    std::vector< std::vector<Double> > field(count);
    std::vector< std::vector<Double> > reffield(count);
    std::vector<Double> mean(points.size(), 0);
    std::vector<Double> refmean(points.size(), 0);

    //logTimerStart;
    for(UInt i=0; i<count; i++)
    {
      //logTimerLoop(i,count);
      SphericalHarmonics harm = gravityfield->sphericalHarmonics(times.at(i));
      field.at(i) = MiscGriddedData::synthesisSphericalHarmonics(harm, points, kernel, comm);

      harm = referencefield->sphericalHarmonics(times.at(i));
      reffield.at(i) = MiscGriddedData::synthesisSphericalHarmonics(harm, points, kernel, comm);
    }
    //logTimerLoopEnd(count);

    if(Parallel::isMaster(comm))
    {
      // Compute mean
      // ------------
      logStatus<<"compute mean"<<Log::endl;
      std::vector<Double> mean(points.size(), 0);
      std::vector<Double> refmean(points.size(), 0);
      //logTimerStart;
      for(UInt i=0; i<count; i++)
      {
        //logTimerLoop(i,count);
        for(UInt k=0; k<points.size(); k++)
          mean.at(k) += field.at(i).at(k)/count;
        for(UInt k=0; k<points.size(); k++)
          refmean.at(k) += reffield.at(i).at(k)/count;
        }
      }
      //logTimerLoopEnd(count);

      std::vector<Double> values(points.size(), 0);
      switch(functional)
      {
        case CORRELATION:
        {
          logStatus<<"compute correlation coefficients"<<Log::endl;
          std::vector<Double> rms(points.size(), 0);
          std::vector<Double> refrms(points.size(), 0);
          std::vector<Double> crossrms(points.size(), 0);

          //logTimerStart;
          for(UInt i=0; i<count; i++) // Alle Zeitpunkt)
          {
            //logTimerLoop(i,count);
            for(UInt k=0; k<points.size(); k++)
            {
              rms.at(k) += pow(field.at(i).at(k)-mean.at(k), 2)/count;
              refrms.at(k) += pow(reffield.at(i).at(k)-refmean.at(k), 2)/count;
              crossrms.at(k) += (field.at(i).at(k)-mean.at(k))*(reffield.at(i).at(k)-refmean.at(k))/count;
            }
          }
          //logTimerLoopEnd(count);

          for(UInt k=0; k<points.size(); k++)
            values.at(k) = crossrms.at(k)/sqrt(rms.at(k)*refrms.at(k));
          break;
        }
        case ERRORRMS:
        {
          logStatus<<"compute error rms"<<Log::endl;
          std::vector<Double> rms(points.size(), 0);

          //logTimerStart;
          for(UInt i=0; i<count; i++) // Alle Zeitpunkt)
          {
            //logTimerLoop(i,count);
            for(UInt k=0; k<points.size(); k++)
              rms.at(k) += pow((field.at(i).at(k)-mean.at(k))-(reffield.at(i).at(k)-refmean.at(k)), 2)/count;
          }
          //logTimerLoopEnd(count);

          for(UInt k=0; k<points.size(); k++)
            values.at(k) = sqrt(rms.at(k));
          break;
        }
        case NSC :
        {
          logStatus<<"compute Nash Sutcliffe coefficients"<<Log::endl;
          std::vector<Double> sum1(points.size(), 0);
          std::vector<Double> sum2(points.size(), 0);

          //logTimerStart;
          for(UInt i=0; i<count; i++) // Alle Zeitpunkte
          {
            //logTimerLoop(i,count);
            for(UInt k=0; k<points.size(); k++)
            {
              sum1.at(k) += pow(field.at(i).at(k)-reffield.at(i).at(k),2);
              sum2.at(k) +=  pow(field.at(i).at(k)-mean.at(k), 2);
            }
          }
          //logTimerLoopEnd(count);

          for(UInt k=0; k<points.size(); k++)
            values.at(k) = 1 - sum1.at(k)/sum2.at(k);
          break;
        }
        case PERCENTAGE:
        {
          logStatus<<"compute percentage of reference signal"<<Log::endl;
          std::vector<Double> diffrms(points.size(), 0);
          std::vector<Double> refrms(points.size(), 0);

          //logTimerStart;
          for(UInt i=0; i<count; i++)
          {
            //logTimerLoop(i,count);
            for(UInt k=0; k<points.size(); k++)
            {
              diffrms.at(k) += pow((field.at(i).at(k)-mean.at(k)) - (reffield.at(i).at(k)-refmean.at(k)), 2)/count;
              refrms.at(k) += pow(reffield.at(i).at(k)-refmean.at(k), 2)/count;
            }
          }
          //logTimerLoopEnd(count);

          for(UInt k=0; k<points.size(); k++)
            values.at(k) = 1 - (diffrms.at(k)/refrms.at(k));
          break;
        }
      }

      logStatus<<"save values to file <"<<outName<<">"<<Log::endl;
	  GriddedData griddedData(Ellipsoid(a,f), points, areas, {values});
      writeFileGriddedData(outName, griddedData);
	  
      //writeFileValueGrid(outName, ValueGrid(Ellipsoid(a,f), points, areas, values));
    }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
