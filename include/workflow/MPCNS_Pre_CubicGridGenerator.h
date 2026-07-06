#pragma once

#include "core/MPCNS_Pre_Parameter.h"

class CubicPeriodicGridGenerator
{
public:
    static bool IsEnabled(Param &param);
    static void GenerateIfEnabled(Param &param);
};
