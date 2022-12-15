/***********************************************/
/**
* @file matrix2EmpiricalCovariance.cpp
*
* @brief Estimate the empirical covariance matrix of a matrix file.
**
* @author Viviana Woehnke
* @date 2022-12-11
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program estimates a spatial and temporal covariance matrix from a time series of \configFile{matrix}{matrix} files according to
\begin{equation}
\M\Sigma(\Delta i) = \frac{1}{N}\sum_{i=1}^N \M x_i \M x_{i+\Delta i}^T
\end{equation}
where $\Delta i$ is given by \config{differenceStep}.

)";

/***********************************************/

#include "programs/program.h"
#include "files/fileMatrix.h"
#include "classes/timeSeries/timeSeries.h"

/***** CLASS ***********************************/

/** @brief Estimate the empirical covariance matrix of an instrument file.
* @ingroup programsGroup */
class Matrix2EmpiricalCovariance
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(Matrix2EmpiricalCovariance, SINGLEPROCESS, "Estimate the empirical covariance matrix of a matrix file", Matrix, Covariance)

/***********************************************/

void Matrix2EmpiricalCovariance::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    FileName 			  	fileNameOut;
    std::vector<FileName> 	fileNameIn;
	Bool            		removeMean, regularize, factorOfTrace;
	UInt            		delay;
	Double					regularizationFactor;
	
	regularize = FALSE;

    readConfig(config, "outputfileCovarianceMatrix", fileNameOut, Config::MUSTSET,  "",  "");
    readConfig(config, "inputfileMatrix",        	 fileNameIn,  Config::MUSTSET,  "",  "");
    readConfig(config, "removeMean",                 removeMean,  Config::DEFAULT,  "1",   "");
    readConfig(config, "differenceStep",             delay,       Config::DEFAULT,  "0",   "choose dt for: x,i(t) - x,j(t+dt)");
	if(readConfigSequence(config, "regularize",      			  Config::OPTIONAL,  "",   "add value to main diagonal"))
	{
	regularize = TRUE;
    readConfig(config, "factorOfTrace",              factorOfTrace,  Config::DEFAULT,  "0",   "multiply factor with trace of matrix");		
	readConfig(config, "regularizationFactor",       regularizationFactor,  Config::DEFAULT,  "0.0",   "add to main diagonal");
	endSequence(config);
	}
	
    if(isCreateSchema(config)) return;

	// ============================

	Matrix data_tmp;
	readFileMatrix(fileNameIn.at(0), data_tmp);
	UInt dim = data_tmp.rows();
	
	UInt count = fileNameIn.size();
	
	Matrix covarianceMatrix(dim, dim);
	std::vector<Matrix> data(count);
	

	// read matrix data: loop over number of files
	for(UInt t=0; t<count; t++)
	{
		Matrix data_t;
		readFileMatrix(fileNameIn.at(t), data_t);
		data.at(t) = data_t;
	};	
		
	// remove temporal mean
	// --------------------
	if (removeMean)
	{
		Vector dataMean(dim);
		for (UInt i=0; i<count; i++)
			dataMean += Vector(data.at(i));
		dataMean *= 1./count;
		for (UInt i=0; i<count; i++)
			data.at(i) -= dataMean;
	}

	UInt count_delay = data.size()-delay;
	// compute cov. mat: loop over number of files - delay
	// --------------------

	Single::forEach(count_delay, [&](UInt t)
	{
		Vector x1(dim), x2(dim);
		x1 = data.at(t);
		x2 = data.at(t+delay);
		
		matMult(1., x1, x2.trans(), covarianceMatrix);
	});
		
	if(delay==0)
      covarianceMatrix.setType(Matrix::SYMMETRIC);
	covarianceMatrix *= 1./count_delay;
	
	// regularize matrix
	if(regularize)
	{	
		// multiply with trace of matrix
		if(factorOfTrace)
		{
			Double r = trace(covarianceMatrix);
			regularizationFactor *= r;
		}
		
		for(UInt i=0; i<dim; i++)
			covarianceMatrix(i,i) += regularizationFactor;
	}
		
	logStatus << "save covariance matrix to <"<<fileNameOut<<">"<< Log::endl;
	writeFileMatrix(fileNameOut, covarianceMatrix);

  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/