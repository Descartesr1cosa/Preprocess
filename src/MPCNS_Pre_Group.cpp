#include "MPCNS_Pre_Group.h"
#include "0_CONST_DEFINE.h"
#include <iostream>
#include <iomanip>

#if MPCNS_IF_METIS == 1
#include "metis.h"
#endif

void Preprocess_Group::Metis_allocate_group(Preprocess_Data_Structure *ptr, Param *par)
{
    int32_t processnum = par->GetInt("proc_num");
    bool ifsplit = par->GetBoo("ifsplit");
    if (processnum == 1)
    {
        int32_t *part = new int32_t[ptr->blk.Getsize1()];
        for (int32_t i = 0; i < ptr->blk.Getsize1(); i++)
            part[i] = 0;
        for (int32_t i = 0; i < ptr->blk.Getsize1(); i++)
        {
            block_proc_number(i) = part[i];
            std::cout << "\tThe Block No. " << i << "\t assigned into Proc_\t" << part[i] << std::endl;
        }

        for (int32_t i = 0; i < ptr->blk.Getsize1(); i++)
        {
            int32_t blk_num = i;
            pro_block_index(block_proc_number(i)).push_back(blk_num);
            block_proc_index(i) = pro_block_index(block_proc_number(i)).size() - 1;
        }

        for (int32_t i = 0; i < processnum; i++)
        {
            proc_num(i) = pro_block_index(i).size();
            int32_t temp = 0;
            for (int32_t j = 0; j < pro_block_index(i).size(); j++)
            {
                block_info &blk_info = ptr->blk(pro_block_index(i)[j]);
                temp += blk_info.ijkmax[0] * blk_info.ijkmax[1] * blk_info.ijkmax[2];
            }
            pro_grid_number(i) = temp;
        }
        return;
    }
    else if (!ifsplit || MPCNS_IF_METIS == 0)
    {
        // Greedy Algorithms
        int32_t grid_number = ptr->blk.Getsize1();
        int32_t proc_number = processnum;
        int temp_grid_point = 0;

        // 初始化，进程中分到的块数为0，分到的网格点数为0
        for (int32_t proc_index = 0; proc_index < proc_number; proc_index++)
        {
            pro_grid_number(proc_index) = 0;
            proc_num(proc_index) = 0;
        }

        for (int32_t blk_index = 0; blk_index < grid_number; blk_index++)
        {
            int temp_grid_point = pro_grid_number(0); // 暂存0进程网格点数
            int set_proc_num = 0;                     // 首先默认分到0进程
            for (int32_t proc_index = 1; proc_index < proc_number; proc_index++)
            {
                // Find minimal proc
                if (temp_grid_point > pro_grid_number(proc_index))
                {
                    set_proc_num = proc_index;
                    temp_grid_point = pro_grid_number(set_proc_num);
                }
            }
            // 将网格块blk_index放入进程号set_proc_num中
            block_proc_number(blk_index) = set_proc_num;
            pro_block_index(set_proc_num).push_back(blk_index);
            block_proc_index(blk_index) = pro_block_index(set_proc_num).size() - 1;
            proc_num(set_proc_num) = pro_block_index(set_proc_num).size();

            block_info &blk_info = ptr->blk(blk_index);
            pro_grid_number(set_proc_num) += blk_info.ijkmax[0] * blk_info.ijkmax[1] * blk_info.ijkmax[2];
        }
        for (int32_t i = 0; i < ptr->blk.Getsize1(); i++)
        {
            std::cout << "\tThe Block No. " << i << "\t assigned into Proc_\t" << block_proc_number(i) << std::endl;
        }
        std::cout << "\t------------------------------------------------------------------------------------------\n";
        for (int32_t i = 0; i < processnum; i++)
            std::cout << "\tThe Proc No. " << i << "\t has \t" << proc_num(i) << "\t Blocks" << std::endl;
        std::cout << "\t------------------------------------------------------------------------------------------\n";
        for (int32_t i = 0; i < processnum; i++)
        {
            if (proc_num(i) == 0)
            {
                std::cout << "\tWarning!\tThe Proc No. " << i << "\t has no Block ! ! ! Refine the splitting scale!\n";
                exit(-1);
            }
        }
        return;
    }

#if MPCNS_IF_METIS == 1

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
        vwgt[i] = ptr->blk(i).ijkmax[0] * ptr->blk(i).ijkmax[1] * ptr->blk(i).ijkmax[2];
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
    {
        block_proc_number(i) = part[i];
        std::cout << "\tThe Block No. " << i << "\t assigned into Proc_\t" << part[i] << std::endl;
    }

    for (int32_t i = 0; i < ptr->blk.Getsize1(); i++)
    {
        int32_t blk_num = i;
        pro_block_index(block_proc_number(i)).push_back(blk_num);
        block_proc_index(i) = pro_block_index(block_proc_number(i)).size() - 1;
    }

    for (int32_t i = 0; i < processnum; i++)
    {
        proc_num(i) = pro_block_index(i).size();
        int32_t temp = 0;
        for (int32_t j = 0; j < pro_block_index(i).size(); j++)
        {
            block_info &blk_info = ptr->blk(pro_block_index(i)[j]);
            temp += blk_info.ijkmax[0] * blk_info.ijkmax[1] * blk_info.ijkmax[2];
        }
        pro_grid_number(i) = temp;
    }

    std::cout << "\t------------------------------------------------------------------------------------------\n";
    for (int32_t i = 0; i < processnum; i++)
        std::cout << "\tThe Proc No. " << i << "\t has \t" << proc_num(i) << "\t Blocks" << std::endl;
    std::cout << "\t------------------------------------------------------------------------------------------\n";
    for (int32_t i = 0; i < processnum; i++)
    {
        if (proc_num(i) == 0)
        {
            std::cout << "\tWarning!\tThe Proc No. " << i << "\t has no Block ! ! ! Refine the splitting scale!\n";
            exit(-1);
        }
    }

    delete[] xadj;
    delete[] adjncy;
    delete[] vwgt;
    delete[] adjwgt;
    delete objval;
    delete[] part;
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
