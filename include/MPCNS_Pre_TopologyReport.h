#pragma once

#include "MPCNS_Pre_Data.h"
#include "MPCNS_Pre_Group.h"
#include "MPCNS_Pre_Parameter.h"

class TopologyReport
{
public:
    static void Generate(Preprocess_Data_Structure *ptr, Preprocess_Group *group, Param *param);
};
