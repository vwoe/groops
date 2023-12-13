/***********************************************/
/**
* @file griddedData2Matrix.cpp
*
* @brief Write grid to matrix file.
*
* @author Torsten Mayer-Guerr
* @date 2010-06-22
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program converts \configFile{inputfileGriddedData}{griddedData}
to \configFile{outputfileMatrix}{matrix} with data columns.
The grid is expressed as ellipsoidal coordinates (longitude, latitude, height)
based on a reference ellipsoid with parameters \config{R} and \config{inverseFlattening}.
The content of the output matrix can be controlled by \config{outColumn} expressions
applied to every grid point. The common data variables for grids are available,
see \reference{dataVariables}{general.parser:dataVariables}.
)";

/***********************************************/

#include "programs/program.h"
#include "parser/dataVariables.h"
#include "files/fileGriddedData.h"
#include "files/fileMatrix.h"
#include "misc/miscGriddedData.h"

/***** CLASS ***********************************/

/** @brief Write grid to matrix file.
* @ingroup programsGroup */
class GriddedData2Matrix
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GriddedData2Matrix, SINGLEPROCESS, "write grid to matrix file", Grid, Matrix)
GROOPS_RENAMED_PROGRAM(Grid2Matrix, GriddedData2Matrix, date2time(2020, 02, 02))

/***********************************************/

void GriddedData2Matrix::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName fileNameOut, fileNameGrid;
    Double   a, f;
    std::vector<ExpressionVariablePtr> expression;

    readConfig(config, "outputfileMatrix",     fileNameOut,  Config::MUSTSET,  "", "point list as matrix with longitude and latitude values in columns and possible additional columns");
    readConfig(config, "inputfileGriddedData", fileNameGrid, Config::MUSTSET,  "",      "");
    readConfig(config, "R",                    a,            Config::DEFAULT,  STRING_DEFAULT_GRS80_a, "reference radius for ellipsoidal coordinates on output");
    readConfig(config, "inverseFlattening",    f,            Config::DEFAULT,  STRING_DEFAULT_GRS80_f, "reference flattening for ellipsoidal coordinates on output, 0: spherical coordinates");
    readConfig(config, "outColumn",            expression,   Config::OPTIONAL, R"(["longitude", "latitude", "height", "data0"])", "expression (variables: longitude, latitude, height, area, data0, data1, ...)");
    if(isCreateSchema(config)) return;

    // create grid
    // -----------
    logStatus<<"read grid file <"<<fileNameGrid<<">"<<Log::endl;
    GriddedData grid;
    readFileGriddedData(fileNameGrid, grid);
    grid.ellipsoid = Ellipsoid(a,f);
    MiscGriddedData::printStatistics(grid);

    // create data variables
    // ---------------------
    VariableList varList;
    addDataVariables(grid, varList);
    std::for_each(expression.begin(), expression.end(), [&](auto expr) {expr->simplify(varList);});

    // calculate output grid
    // ---------------------
    logStatus<<"calculate output matrix"<<Log::endl;
    Matrix A(grid.points.size(), expression.size());
    Single::forEach(A.rows(), [&](UInt i)
    {
      evaluateDataVariables(grid, i, varList);
      for(UInt k=0; k<A.columns(); k++)
        A(i,k) = expression.at(k)->evaluate(varList);
    });

    // write
    // -----
    logStatus<<"write matrix to file <"<<fileNameOut<<">"<<Log::endl;
    writeFileMatrix(fileNameOut, A);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
