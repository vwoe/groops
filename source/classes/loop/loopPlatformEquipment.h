/***********************************************/
/**
* @file loopPlatformEquipment.h
*
* @brief Loop over equipment of a platform file.
*
* @author Torsten Mayer-Guerr
* @date 2022-11-11
*
*/
/***********************************************/

#ifndef __GROOPS_LOOPPLATFORMEQUIPMENT__
#define __GROOPS_LOOPPLATFORMEQUIPMENT__

// Latex documentation
#ifdef DOCSTRING_Loop
static const char *docstringLoopPlatformEquipment = R"(
\subsection{PlatformEquipment}
Loop over specific equipment of a \file{platform file}{platform}.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "classes/loop/loop.h"
#include "files/filePlatform.h"

/***** CLASS ***********************************/

/** @brief Loop over equipment of a platform file.
* @ingroup LoopGroup
* @see Loop */
class LoopPlatformEquipment : public Loop
{
  Platform                platform;
  PlatformEquipment::Type type;
  std::string             nameName, nameSerial, nameInfo, nameTimeStart, nameTimeEnd;
  std::string             nameIndex, nameCount;
  std::string             namePositionX, namePositionY, namePositionZ;
  UInt                    index;

public:
  LoopPlatformEquipment(Config &config);

  UInt count() const override {return std::count_if(platform.equipments.begin(), platform.equipments.end(), [&](const auto &x)
                                                   {return (type == PlatformEquipment::UNDEFINED) || (x->getType() == type);});}
  Bool iteration(VariableList &varList) override;
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline LoopPlatformEquipment::LoopPlatformEquipment(Config &config)
{
  try
  {
    FileName    fileName;
    std::string choice;

    readConfig(config, "inputfilePlatform", fileName, Config::MUSTSET, "", "platform info file");
    if(readConfigChoice(config, "equipmentType", choice, Config::MUSTSET, "", "equipment type to loop over"))
    {
      if(readConfigChoiceElement(config, "all",          choice, "loop over all types"))    type = PlatformEquipment::UNDEFINED;
      if(readConfigChoiceElement(config, "gnssAntenna",  choice, "loop over antennas"))     type = PlatformEquipment::GNSSANTENNA;
      if(readConfigChoiceElement(config, "gnssReceiver", choice, "loop over receivers"))    type = PlatformEquipment::GNSSRECEIVER;
      if(readConfigChoiceElement(config, "other",        choice, "loop over other types"))  type = PlatformEquipment::OTHER;
      endChoice(config);
    }
    readConfig(config, "variableLoopName",      nameName,            Config::OPTIONAL,  "loopName",      "variable with name");
    readConfig(config, "variableLoopSerial",    nameSerial,          Config::OPTIONAL,  "loopSerial",    "variable with serial");
    readConfig(config, "variableLoopInfo",      nameInfo,            Config::OPTIONAL,  "loopInfo",      "variable with radome (antenna) or version (receiver)");
    readConfig(config, "variableLoopTimeStart", nameTimeStart,       Config::OPTIONAL,  "loopTimeStart", "variable with start time");
    readConfig(config, "variableLoopTimeEnd",   nameTimeEnd,         Config::OPTIONAL,  "loopTimeEnd",   "variable with end time");
    readConfig(config, "variableLoopPositionX", namePositionX,       Config::OPTIONAL,  "loopPositionX", "variable with position x");
    readConfig(config, "variableLoopPositionY", namePositionY,       Config::OPTIONAL,  "loopPositionY", "variable with position y");
    readConfig(config, "variableLoopPositionY", namePositionZ,       Config::OPTIONAL,  "loopPositionZ", "variable with position z");
    readConfig(config, "variableLoopIndex",     nameIndex,           Config::OPTIONAL,  "",              "variable with index of current iteration (starts with zero)");
    readConfig(config, "variableLoopCount",     nameCount,           Config::OPTIONAL,  "",              "variable with total number of iterations");
    if(isCreateSchema(config)) return;

    readFilePlatform(fileName, platform);
    index = 0;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Bool LoopPlatformEquipment::iteration(VariableList &varList)
{
  if(index >= count())
    return FALSE;

  UInt idx = 0;
  for(const auto &eq : platform.equipments)
    if(((type == PlatformEquipment::UNDEFINED) || (eq->getType() == type)) && (idx++ == index))
    {
      if(!nameIndex.empty())     addVariable(nameIndex,     index,               varList);
      if(!nameCount.empty())     addVariable(nameCount,     count(),             varList);
      if(!nameName.empty())      addVariable(nameName,      eq->name,            varList);
      if(!nameSerial.empty())    addVariable(nameSerial,    eq->serial,          varList);
      if(!nameTimeStart.empty()) addVariable(nameTimeStart, eq->timeStart.mjd(), varList);
      if(!nameTimeEnd.empty())   addVariable(nameTimeEnd,   eq->timeEnd.mjd(),   varList);
      if(!namePositionX.empty()) addVariable(namePositionX, eq->position.x(),    varList);
      if(!namePositionY.empty()) addVariable(namePositionY, eq->position.y(),    varList);
      if(!namePositionZ.empty()) addVariable(namePositionZ, eq->position.z(),    varList);
      if(!nameInfo.empty())
      {
        switch(eq->getType())
        {
          case PlatformEquipment::GNSSANTENNA:  addVariable(nameInfo, std::dynamic_pointer_cast<PlatformGnssAntenna>(eq)->radome,   varList); break;
          case PlatformEquipment::GNSSRECEIVER: addVariable(nameInfo, std::dynamic_pointer_cast<PlatformGnssReceiver>(eq)->version, varList); break;
          case PlatformEquipment::OTHER:        break;
          case PlatformEquipment::UNDEFINED:    break;
          default:                              break;
        }
      }
      break;
    }

  index++;
  return TRUE;
}

/***********************************************/

#endif
