/***********************************************/
/**
* @file loopMatrix.h
*
* @brief Loop over rows of a matrix.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2017-01-27
*
*/
/***********************************************/

#ifndef __GROOPS_LOOPMATRIX__
#define __GROOPS_LOOPMATRIX__

// Latex documentation
#ifdef DOCSTRING_Loop
static const char *docstringLoopMatrix = R"(
\subsection{Matrix}
Loop over rows of a matrix. To define the loop variables the standard
data variables of the matrix are available, see~\reference{dataVariables}{general.parser:dataVariables}.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "parser/dataVariables.h"
#include "files/fileMatrix.h"
#include "classes/loop/loop.h"

/***** CLASS ***********************************/

/** @brief Loop over rows of a matrix.
* @ingroup loopGroup
* @see Loop */
class LoopMatrix : public Loop
{
  Matrix          A;
  std::vector<ExpressionVariablePtr> variables;
  std::string     nameIndex, nameCount;
  VariableList    varListMatrix;

public:
  LoopMatrix(Config &config);

  UInt count() const override {return A.rows();}
  Bool iteration(VariableList &varList) override;
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline LoopMatrix::LoopMatrix(Config &config)
{
  try
  {
    FileName              fileName;
    ExpressionVariablePtr startRow, countRows;
    Bool                  transpose;

    readConfig(config, "inputfile",          fileName,  Config::MUSTSET,  "",                 "");
    readConfig(config, "transpose",          transpose, Config::DEFAULT,  "0",                "effectively loop over columns");
    readConfig(config, "startRow",           startRow,  Config::DEFAULT,  "0",                "start at this row (variable: rows)");
    readConfig(config, "countRows",          countRows, Config::DEFAULT,  "rows",             "use this many rows (variable: rows)");
    readConfig(config, "variableLoop",       variables, Config::MUSTSET,  "loopNumber=data0", "define a variable by name = expression (input columns are named data0, data1, ...)");
    readConfig(config, "variableLoopIndex",  nameIndex, Config::OPTIONAL, "",                 "variable with index of current iteration (starts with zero)");
    readConfig(config, "variableLoopCount",  nameCount, Config::OPTIONAL, "",                 "variable with total number of iterations");
    readConfigCondition(config);
    if(isCreateSchema(config)) return;

    readFileMatrix(fileName, A);
    if(transpose)
      A = A.trans();

    varListMatrix.setVariable("rows", A.rows());
    A = A.row(static_cast<UInt>(startRow->evaluate(varListMatrix)), static_cast<UInt>(countRows->evaluate(varListMatrix)));

    for(auto &variable : variables)
      variable->parseVariableName(); // get real variable names, otherwise all named 'variableLoop'
    addDataVariables(A, varListMatrix);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Bool LoopMatrix::iteration(VariableList &varList)
{
  if(index() >= count())
    return FALSE;

  evaluateDataVariables(A, index(), varListMatrix);
  for(auto &variable : variables)
    varList.setVariable(variable->name(), variable->evaluate(varListMatrix));
  if(!nameIndex.empty()) varList.setVariable(nameIndex, index());
  if(!nameCount.empty()) varList.setVariable(nameCount, count());

  return checkCondition(varList);
}

/***********************************************/

#endif
