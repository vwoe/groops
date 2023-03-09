/***********************************************/
/**
* @file gnssGlonassFrequencyNumberUpdate.cpp
*
* @brief Update/set GLONASS frequency number in transmitter info files.
*
* @author Sebastian Strasser
* @date 2019-08-29
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Update/set GLONASS frequency number in \configFile{inputfileTransmitterInfo}{platform} files.

PRN/SVN to frequency number source: \url{http://semisys.gfz-potsdam.de/semisys/api/?symname=2002&format=json&satellite=GLO}.

See also \program{GnssAntex2AntennaDefinition}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileGnssReceiverDefinition.h"
#include "files/fileGnssStationInfo.h"
#include "files/fileMatrix.h"

/***** CLASS ***********************************/

/** @brief Update/set GLONASS frequency number in transmitter info files.
* @ingroup programsGroup */
class GnssGlonassFrequencyNumberUpdate
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GnssGlonassFrequencyNumberUpdate, SINGLEPROCESS, "Update/set GLONASS frequency number in transmitter info files.", Gnss)

/***********************************************/

void GnssGlonassFrequencyNumberUpdate::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName fileNameOutTransmitterInfo, fileNameInTransmitterInfo, fileNameInPrnSvn2FrequencyNumber;
    std::vector<std::string> prnList;
    std::string variableNamePrn;

    readConfig(config, "outputfileTransmitterInfo",    fileNameOutTransmitterInfo,       Config::OPTIONAL, "",    "templated for PRN list (variableNamePrn)");
    readConfig(config, "inputfileTransmitterInfo",     fileNameInTransmitterInfo,        Config::MUSTSET,  "",    "templated for PRN list (variableNamePrn)");
    readConfig(config, "inputfilePrn2FrequencyNumber", fileNameInPrnSvn2FrequencyNumber, Config::MUSTSET,  "",    "GROOPS matrix with columns: GLONASS PRN, SVN, mjdStart, mjdEnd, frequencyNumber");
    readConfig(config, "prn",                          prnList,                          Config::OPTIONAL, "",    "PRN (e.g. G01) for transmitter info files");
    readConfig(config, "variableNamePrn",              variableNamePrn,                  Config::OPTIONAL, "prn", "variable name for PRN in transmitter info files");
    if(isCreateSchema(config)) return;

    VariableList varList;
    std::vector<GnssStationInfo> transmitterInfos(prnList.size());
    if(!fileNameInTransmitterInfo.empty())
    {
      logStatus<<"read transmitter infos from <"<<fileNameInTransmitterInfo<<">"<< Log::endl;
      for(UInt idPrn = 0; idPrn < prnList.size(); idPrn++)
      {
        varList[variableNamePrn]->setValue(prnList.at(idPrn));
        readFileGnssStationInfo(fileNameInTransmitterInfo(varList), transmitterInfos.at(idPrn));
      }
    }

    std::vector<GnssReceiverDefinitionPtr> receiverDefList;

    logStatus<<"read GLONASS PRN/SVN to frequency number matrix from <"<<fileNameInPrnSvn2FrequencyNumber<<">"<< Log::endl;
    Matrix prnSvn2FreqNo;
    readFileMatrix(fileNameInPrnSvn2FrequencyNumber, prnSvn2FreqNo);

    //---------------------------------------------------------------

    std::vector<std::vector<GnssReceiverInfo>> receiverInfos(prnList.size());
    std::vector<GnssReceiverDefinitionPtr> receiverDefListNew;
    for(UInt i = 0; i < prnSvn2FreqNo.rows(); i++)
    {
      const std::string prn    = prnSvn2FreqNo(i,0)%"R%02i"s;
      const std::string svn    = prnSvn2FreqNo(i,1)%"R%03i"s;
      const Time timeStart     = mjd2time(prnSvn2FreqNo(i,2));
      const Time timeEnd       = mjd2time(prnSvn2FreqNo(i,3));
      const std::string freqNo = prnSvn2FreqNo(i,4)%"%i"s;

      Bool found = FALSE;
      for(UInt idTrans = 0; idTrans < transmitterInfos.size(); idTrans++)
        if(transmitterInfos.at(idTrans).markerNumber == prn)
        {
          // new receiver list entry
          UInt idRecv = transmitterInfos.at(idTrans).findReceiver(0.5*(timeStart+timeEnd));
          if(idRecv == NULLINDEX)
            continue;
          found = TRUE;

          GnssReceiverInfo info(transmitterInfos.at(idTrans).receiver.at(idRecv));
          info.timeStart = timeStart;
          info.timeEnd   = timeEnd;
          info.version   = freqNo;
          receiverInfos.at(idTrans).push_back(info);
        }

      if(!found)
        logWarning<<prn<<": no transmitter info entry found for time period "<<timeStart.dateTimeStr()<<" to "<<timeEnd.dateTimeStr()<<Log::endl;
    }

    //---------------------------------------------------------------

    if(!fileNameOutTransmitterInfo.empty())
    {
      logStatus<<"write transmitter infos to <"<<fileNameOutTransmitterInfo<<">"<< Log::endl;
      for(UInt idPrn = 0; idPrn < prnList.size(); idPrn++)
      {
        varList[variableNamePrn]->setValue(prnList.at(idPrn));
        std::sort(receiverInfos.at(idPrn).begin(), receiverInfos.at(idPrn).end(), [](GnssReceiverInfo &info1, GnssReceiverInfo &info2){ return info1.timeStart < info2.timeStart; });
        transmitterInfos.at(idPrn).receiver = receiverInfos.at(idPrn);
        writeFileGnssStationInfo(fileNameOutTransmitterInfo(varList), transmitterInfos.at(idPrn));
      }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
