#include "partition/MPCNS_Pre_Group.h"
#include "core/0_CONST_DEFINE.h"
#include "partition/MPCNS_Pre_BlockWeight.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <iostream>
#include <iomanip>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

#if MPCNS_IF_METIS == 1
#include "metis.h"
#endif

namespace
{
// Set by Metis_allocate_group; BlockWeight() is the weighted cost used for balance.
Param *g_weight_par = nullptr;

struct BlockPartitionInfo
{
    int32_t id;
    int64_t weight;
    double center[3];
};

int64_t BlockWeight(const block_info &blk_info)
{
    return BlockAssignWeight(blk_info, g_weight_par);
}

BlockPartitionInfo MakeBlockInfo(Preprocess_Data_Structure *ptr, int32_t id)
{
    block_info &blk_info = ptr->blk(id);
    BlockPartitionInfo info;
    info.id = id;
    info.weight = BlockWeight(blk_info);
    info.center[0] = 0.5 * (blk_info.ijkmax[0] + 1);
    info.center[1] = 0.5 * (blk_info.ijkmax[1] + 1);
    info.center[2] = 0.5 * (blk_info.ijkmax[2] + 1);

    double total_volume = 0.0;
    double weighted_center[3] = {0.0, 0.0, 0.0};
    for (int32_t ihier = 0; ihier < blk_info.hiera.size(); ihier++)
    {
        coordinate_blk &hier = blk_info.hiera[ihier];
        double volume = 1.0;
        double center[3];
        for (int32_t d = 0; d < 3; d++)
        {
            volume *= std::abs(hier.oring_sup[d] - hier.oring_sub[d]) + 1;
            center[d] = 0.5 * (hier.oring_sub[d] + hier.oring_sup[d]);
        }
        total_volume += volume;
        for (int32_t d = 0; d < 3; d++)
            weighted_center[d] += center[d] * volume;
    }

    if (total_volume > 0.0)
    {
        for (int32_t d = 0; d < 3; d++)
            info.center[d] = weighted_center[d] / total_volume;
    }
    return info;
}

std::vector<BlockPartitionInfo> MakeBlockInfos(Preprocess_Data_Structure *ptr)
{
    std::vector<BlockPartitionInfo> blocks;
    blocks.reserve(ptr->blk.Getsize1());
    for (int32_t i = 0; i < ptr->blk.Getsize1(); i++)
        blocks.push_back(MakeBlockInfo(ptr, i));
    return blocks;
}

void PrintAndValidatePartition(Preprocess_Group *group, Preprocess_Data_Structure *ptr, int32_t processnum)
{
    for (int32_t i = 0; i < ptr->blk.Getsize1(); i++)
    {
        std::cout << "\tThe Block No. " << i << "\t assigned into Proc_\t" << group->block_proc_number(i) << std::endl;
    }
    std::cout << "\t------------------------------------------------------------------------------------------\n";
    for (int32_t i = 0; i < processnum; i++)
        std::cout << "\tThe Proc No. " << i << "\t has \t" << group->proc_num(i) << "\t Blocks" << std::endl;
    std::cout << "\t------------------------------------------------------------------------------------------\n";
    {
        std::vector<int64_t> weighted(processnum, 0);
        int64_t total = 0;
        for (int32_t i = 0; i < ptr->blk.Getsize1(); i++)
        {
            const int64_t w = BlockWeight(ptr->blk(i));
            weighted[group->block_proc_number(i)] += w;
            total += w;
        }
        const int64_t max_w = *std::max_element(weighted.begin(), weighted.end());
        const int64_t min_w = *std::min_element(weighted.begin(), weighted.end());
        const double avg_w = static_cast<double>(total) / processnum;
        std::cout << "\tWeighted load (fluid_cost_weight applied): avg = " << avg_w << ", max = " << max_w
                  << ", min = " << min_w << ", max/avg = " << (avg_w > 0.0 ? max_w / avg_w : 0.0) << std::endl;
        std::cout << "\t------------------------------------------------------------------------------------------\n";
    }
    for (int32_t i = 0; i < processnum; i++)
    {
        if (group->proc_num(i) == 0)
        {
            std::cout << "\tWarning!\tThe Proc No. " << i << "\t has no Block ! ! ! Refine the splitting scale!\n";
            exit(-1);
        }
    }
}

void ApplyPartition(Preprocess_Group *group, Preprocess_Data_Structure *ptr, const std::vector<int32_t> &part, int32_t processnum)
{
    for (int32_t i = 0; i < processnum; i++)
    {
        group->pro_grid_number(i) = 0;
        group->proc_num(i) = 0;
        group->pro_block_index(i).clear();
    }

    for (int32_t i = 0; i < ptr->blk.Getsize1(); i++)
    {
        const int32_t proc = part[i];
        group->block_proc_number(i) = proc;
        group->pro_block_index(proc).push_back(i);
        group->block_proc_index(i) = group->pro_block_index(proc).size() - 1;
        group->proc_num(proc) = group->pro_block_index(proc).size();
        group->pro_grid_number(proc) += static_cast<int32_t>(BlockCells(ptr->blk(i))); // raw cells, as stored in the group file
    }

    PrintAndValidatePartition(group, ptr, processnum);
}

std::vector<int32_t> GreedyPartition(Preprocess_Data_Structure *ptr, int32_t processnum)
{
    std::vector<int32_t> part(ptr->blk.Getsize1(), 0);
    std::vector<int64_t> load(processnum, 0);
    for (int32_t blk_index = 0; blk_index < ptr->blk.Getsize1(); blk_index++)
    {
        int32_t proc = 0;
        for (int32_t proc_index = 1; proc_index < processnum; proc_index++)
        {
            if (load[proc] > load[proc_index])
                proc = proc_index;
        }
        part[blk_index] = proc;
        load[proc] += BlockWeight(ptr->blk(blk_index));
    }
    return part;
}

// Longest-processing-time-first: blocks in descending weight, each to the least
// loaded process.  Best balance for a given block set; ignores block locality.
std::vector<int32_t> LptPartition(Preprocess_Data_Structure *ptr, int32_t processnum)
{
    std::vector<int32_t> order(ptr->blk.Getsize1());
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(), [&](int32_t a, int32_t b) {
        return BlockWeight(ptr->blk(a)) > BlockWeight(ptr->blk(b));
    });
    std::vector<int32_t> part(ptr->blk.Getsize1(), 0);
    std::vector<int64_t> load(processnum, 0);
    std::vector<int32_t> count(processnum, 0);
    for (int32_t blk_index : order)
    {
        int32_t proc = 0;
        for (int32_t p = 1; p < processnum; p++)
            if (load[p] < load[proc])
                proc = p;
        part[blk_index] = proc;
        load[proc] += BlockWeight(ptr->blk(blk_index));
        count[proc]++;
    }
    return part;
}

uint64_t ExpandMortonBits(uint32_t value)
{
    uint64_t x = value & 0x1fffff;
    x = (x | x << 32) & 0x1f00000000ffff;
    x = (x | x << 16) & 0x1f0000ff0000ff;
    x = (x | x << 8) & 0x100f00f00f00f00f;
    x = (x | x << 4) & 0x10c30c30c30c30c3;
    x = (x | x << 2) & 0x1249249249249249;
    return x;
}

std::vector<int32_t> MortonPartition(Preprocess_Data_Structure *ptr, int32_t processnum)
{
    std::vector<BlockPartitionInfo> blocks = MakeBlockInfos(ptr);
    double min_center[3] = {blocks[0].center[0], blocks[0].center[1], blocks[0].center[2]};
    double max_center[3] = {blocks[0].center[0], blocks[0].center[1], blocks[0].center[2]};
    for (int32_t i = 0; i < blocks.size(); i++)
    {
        for (int32_t d = 0; d < 3; d++)
        {
            min_center[d] = std::min(min_center[d], blocks[i].center[d]);
            max_center[d] = std::max(max_center[d], blocks[i].center[d]);
        }
    }

    std::sort(blocks.begin(), blocks.end(), [&](const BlockPartitionInfo &a, const BlockPartitionInfo &b) {
        uint64_t code_a = 0;
        uint64_t code_b = 0;
        for (int32_t d = 0; d < 3; d++)
        {
            const double range = max_center[d] - min_center[d];
            const double normalized_a = (range > 0.0) ? (a.center[d] - min_center[d]) / range : 0.0;
            const double normalized_b = (range > 0.0) ? (b.center[d] - min_center[d]) / range : 0.0;
            const uint32_t ia = static_cast<uint32_t>(std::max(0.0, std::min(2097151.0, normalized_a * 2097151.0)));
            const uint32_t ib = static_cast<uint32_t>(std::max(0.0, std::min(2097151.0, normalized_b * 2097151.0)));
            code_a |= ExpandMortonBits(ia) << d;
            code_b |= ExpandMortonBits(ib) << d;
        }
        if (code_a == code_b)
            return a.id < b.id;
        return code_a < code_b;
    });

    int64_t total_weight = 0;
    for (int32_t i = 0; i < blocks.size(); i++)
        total_weight += blocks[i].weight;

    std::vector<int32_t> part(ptr->blk.Getsize1(), 0);
    int64_t accumulated = 0;
    int32_t proc = 0;
    for (int32_t i = 0; i < blocks.size(); i++)
    {
        part[blocks[i].id] = proc;
        accumulated += blocks[i].weight;
        const int32_t remaining_blocks = blocks.size() - i - 1;
        while (proc + 1 < processnum)
        {
            const int32_t remaining_procs_after_advance = processnum - (proc + 1);
            const int64_t target = (total_weight * (proc + 1)) / processnum;
            if ((accumulated >= target || remaining_blocks == remaining_procs_after_advance) &&
                remaining_blocks >= remaining_procs_after_advance)
            {
                proc++;
            }
            else
            {
                break;
            }
        }
    }
    return part;
}

void RcbPartitionRange(const std::vector<BlockPartitionInfo> &input_blocks, int32_t proc_begin, int32_t proc_end, std::vector<int32_t> &part)
{
    if (input_blocks.empty())
        return;
    if (proc_end - proc_begin <= 1 || input_blocks.size() == 1)
    {
        for (int32_t i = 0; i < input_blocks.size(); i++)
            part[input_blocks[i].id] = proc_begin;
        return;
    }

    double min_center[3] = {input_blocks[0].center[0], input_blocks[0].center[1], input_blocks[0].center[2]};
    double max_center[3] = {input_blocks[0].center[0], input_blocks[0].center[1], input_blocks[0].center[2]};
    for (int32_t i = 0; i < input_blocks.size(); i++)
    {
        for (int32_t d = 0; d < 3; d++)
        {
            min_center[d] = std::min(min_center[d], input_blocks[i].center[d]);
            max_center[d] = std::max(max_center[d], input_blocks[i].center[d]);
        }
    }

    int32_t axis = 0;
    for (int32_t d = 1; d < 3; d++)
    {
        if (max_center[d] - min_center[d] > max_center[axis] - min_center[axis])
            axis = d;
    }

    std::vector<BlockPartitionInfo> blocks(input_blocks);
    std::sort(blocks.begin(), blocks.end(), [axis](const BlockPartitionInfo &a, const BlockPartitionInfo &b) {
        if (a.center[axis] == b.center[axis])
            return a.id < b.id;
        return a.center[axis] < b.center[axis];
    });

    const int32_t proc_count = proc_end - proc_begin;
    const int32_t left_proc_count = proc_count / 2;
    int64_t total_weight = 0;
    for (int32_t i = 0; i < blocks.size(); i++)
        total_weight += blocks[i].weight;
    const int64_t target_left_weight = total_weight * left_proc_count / proc_count;

    int32_t split = 1;
    int64_t accumulated = 0;
    int64_t best_diff = std::numeric_limits<int64_t>::max();
    for (int32_t i = 0; i < blocks.size() - 1; i++)
    {
        accumulated += blocks[i].weight;
        const int64_t diff = std::llabs(accumulated - target_left_weight);
        if (diff < best_diff)
        {
            best_diff = diff;
            split = i + 1;
        }
    }

    std::vector<BlockPartitionInfo> left(blocks.begin(), blocks.begin() + split);
    std::vector<BlockPartitionInfo> right(blocks.begin() + split, blocks.end());
    RcbPartitionRange(left, proc_begin, proc_begin + left_proc_count, part);
    RcbPartitionRange(right, proc_begin + left_proc_count, proc_end, part);
}

std::vector<int32_t> RcbPartition(Preprocess_Data_Structure *ptr, int32_t processnum)
{
    std::vector<BlockPartitionInfo> blocks = MakeBlockInfos(ptr);
    std::vector<int32_t> part(ptr->blk.Getsize1(), 0);
    RcbPartitionRange(blocks, 0, processnum, part);
    return part;
}

std::string GetPartitionMethod(Param *par)
{
    if (par->HasStr("partition_method"))
        return par->GetStr("partition_method");
#if MPCNS_IF_METIS == 1
    return "metis";
#else
    return "morton";
#endif
}
} // namespace

void Preprocess_Group::Metis_allocate_group(Preprocess_Data_Structure *ptr, Param *par)
{
    int32_t processnum = par->GetInt("proc_num");
    g_weight_par = par;
    if (processnum == 1)
    {
        ApplyPartition(this, ptr, std::vector<int32_t>(ptr->blk.Getsize1(), 0), processnum);
        return;
    }

    std::string partition_method = GetPartitionMethod(par);
    if (partition_method == "greedy")
    {
        ApplyPartition(this, ptr, GreedyPartition(ptr, processnum), processnum);
        return;
    }
    if (partition_method == "lpt")
    {
        ApplyPartition(this, ptr, LptPartition(ptr, processnum), processnum);
        return;
    }
    if (partition_method == "morton")
    {
        ApplyPartition(this, ptr, MortonPartition(ptr, processnum), processnum);
        return;
    }
    if (partition_method == "rcb")
    {
        ApplyPartition(this, ptr, RcbPartition(ptr, processnum), processnum);
        return;
    }
    if (partition_method != "metis")
    {
        std::cout << "\tWarning!\tUnknown partition_method = " << partition_method << ", use morton instead.\n";
        ApplyPartition(this, ptr, MortonPartition(ptr, processnum), processnum);
        return;
    }

#if MPCNS_IF_METIS == 0
    std::cout << "\tWarning!\tpartition_method = metis, but METIS is disabled. Use morton instead.\n";
    ApplyPartition(this, ptr, MortonPartition(ptr, processnum), processnum);
    return;
#endif

#if MPCNS_IF_METIS == 1
    std::vector<int32_t> metis_part(ptr->blk.Getsize1(), 0);
    {
        int32_t nvtxs, ncon, nparts;
        int32_t *xadj, *adjncy, *vwgt, *vsize, *adjwgt;
        double *tpwgts, *ubvec;
        int32_t *objval, *part;

        objval = new int32_t;
        part = new int32_t[ptr->blk.Getsize1()];

        // The number of vertices in the graph.
        nvtxs = ptr->blk.Getsize1();

        // The number of balancing constraints. It should be at least 1.
        ncon = 1;

        // The number of parts to partition the graph.
        nparts = processnum;

        // This is an array of size nparts×ncon that specifies the desired weight for each partition and constraint.
        tpwgts = NULL;

        // This is an array of size ncon that specifies the allowed load imbalance tolerance for each constraint.
        ubvec = NULL;

        // The size of the vertices for computing the total communication volume
        vsize = NULL;

        // !xadj adjncy 用于表示一个无向图 xadj代表该顶点相连的顶点数 adjncy代表相连的顶点ID
        // !part返回顶点的分区向量——数组大小为顶点数，part(i)的值为其所在的分区
        xadj = new int32_t[ptr->blk.Getsize1() + 1];
        vwgt = new int32_t[ptr->blk.Getsize1()];
        int32_t count = 0;
        for (int32_t i = 0; i < ptr->blk.Getsize1(); i++)
        {
            for (int32_t j = 0; j < ptr->blk.Getsize1(); j++)
            {
                if (j == i)
                    continue;
                count += ptr->inp(i, j).inner.size();
            }
        }
        adjncy = new int32_t[count];
        adjwgt = new int32_t[count];
        // 设置xadj adjncy
        count = 0;
        for (int32_t i = 0; i < ptr->blk.Getsize1(); i++)
        {
            part[i] = 0;
            xadj[i] = count;
            vwgt[i] = static_cast<int32_t>(BlockWeight(ptr->blk(i)));
            for (int32_t j = 0; j < ptr->blk.Getsize1(); j++)
            {
                if (j == i)
                    continue;
                for (int32_t k = 0; k < ptr->inp(i, j).inner.size(); k++)
                {
                    adjncy[count] = j;
                    adjwgt[count] = ptr->inp(i, j).inner[k].volume;
                    count++;
                }
            }
        }
        xadj[ptr->blk.Getsize1()] = count;
        //-----------------------------------------------
        // 设置options
        int32_t options[METIS_NOPTIONS];
        for (int i = 0; i < 40; i++)
            options[i] = 0;
        METIS_SetDefaultOptions(options);
        //! 0 C-style 数组 1 fortran-style 数组
        options[METIS_OPTION_NUMBERING] = 0;
        // 0 通信条数最小化 1 通信量最小化
        options[METIS_OPTION_OBJTYPE] = 1;

        int32_t return_M = METIS_PartGraphKway(&nvtxs, &ncon, xadj, adjncy, vwgt, NULL, adjwgt, &nparts, NULL, NULL, options, objval, part);
        // int32_t return_M = METIS_PartGraphKway(&nvtxs, &ncon, xadj, adjncy, vwgt, vsize, adjwgt, &nparts, tpwgts, NULL, options, objval, part);

        for (int32_t i = 0; i < ptr->blk.Getsize1(); i++)
            metis_part[i] = part[i];

        delete[] xadj;
        delete[] adjncy;
        delete[] vwgt;
        delete[] adjwgt;
        delete objval;
        delete[] part;
        (void)return_M;
    }
    ApplyPartition(this, ptr, metis_part, processnum);
    return;
#endif
}

void Preprocess_Group::Load_Balance_Output()
{
    std::ofstream file;
    file.open("./OUTPUT/load_balance.dat", std::ios_base::out);
    if (!file.is_open())
    {
        std::cout << "Error, cannot write the file:\t" << "./OUTPUT/load_balance.dat, Please double check! ! !\n";
        std::cout << "Press any key to exit...\n";
        std::string command;
        std::cin >> command;
        exit(-1);
    }
    file << "variables=Number_of_Process, Grid_Load\n";
    for (int32_t i = 0; i < pro_grid_number.Getsize1(); i++)
    {
        file << i << "\t" << pro_grid_number(i) << "\n";
    }
    file.close();
}
