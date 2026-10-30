#pragma once
#include "core/MPCNS_Pre_Data.h"
#include <cmath>
#include <cstdint>
#include <string>

// Relative computing cost of one cell of a block, used by both the automatic
// split and the block-to-process assignment.
//   solid blocks (name == solid_block_name, default "Solid") cost 1.0 per cell;
//   every other block (fluid) costs fluid_cost_weight per cell (default 1.0,
//   which reproduces the plain cell-count balance).
inline double BlockCostFactor(const block_info &blk, Param *par)
{
    if (par == nullptr || !par->HasDou("fluid_cost_weight"))
        return 1.0;
    const double fluid_weight = par->GetDou("fluid_cost_weight");
    const std::string solid_name =
        par->HasStr("solid_block_name") ? par->GetStr("solid_block_name") : "Solid";
    if (!blk.hiera.empty() && blk.hiera[0].blk_name == solid_name)
        return 1.0;
    return fluid_weight;
}

inline int64_t BlockCells(const block_info &blk)
{
    return static_cast<int64_t>(blk.ijkmax[0]) * blk.ijkmax[1] * blk.ijkmax[2];
}

inline int64_t BlockWeightedCells(const block_info &blk, Param *par)
{
    return static_cast<int64_t>(std::llround(BlockCells(blk) * BlockCostFactor(blk, par)));
}

// Ghost-inclusive cell count of a block: prod_d (cells_d + 2*g) over the active
// dimensions, where cells_d = ijkmax[d]-1 and g = cost_ghost_layers.  The solver
// trial run shows per-rank compute follows this allocated volume, which also
// explains the apparent fixed per-block cost.  With g = 0 (default) the legacy
// plain ijkmax product is used.
inline double BlockAllocCells(const block_info &blk, Param *par)
{
    const int g = (par != nullptr && par->HasInt("cost_ghost_layers")) ? par->GetInt("cost_ghost_layers") : 0;
    if (g <= 0)
        return static_cast<double>(BlockCells(blk));
    const int dim = (par != nullptr && par->HasInt("dimension")) ? par->GetInt("dimension") : 3;
    double cells = 1.0;
    for (int d = 0; d < 3; ++d)
        cells *= (d < dim) ? static_cast<double>(blk.ijkmax[d] - 1 + 2 * g) : 1.0;
    return cells;
}

// Weight used when assigning blocks to processes: cost factor times the
// (ghost-inclusive) cell count plus a fixed per-block overhead, in solid-cell
// equivalents.  block_cost_overhead defaults to 0.
inline int64_t BlockAssignWeight(const block_info &blk, Param *par)
{
    const double overhead =
        (par != nullptr && par->HasDou("block_cost_overhead")) ? par->GetDou("block_cost_overhead") : 0.0;
    return static_cast<int64_t>(std::llround(BlockAllocCells(blk, par) * BlockCostFactor(blk, par) + overhead));
}
