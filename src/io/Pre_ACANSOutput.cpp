#include "io/Pre_ACANSOutput.h"
#include <iostream>
#include <iomanip>

void OUTPUT_ACANS::output_inner_file()
{
    std::ofstream file;
    int32_t procnum = pr->GetInt("proc_num");
    file.open("./OUTPUT/geometry/boundary_condition/inner   0.txt", std::ios_base::out);
    if (!file)
    {
        std::string command;
        command = "mkdir OUTPUT/geometry/boundary_condition";
        std::system(command.c_str());
        command = "mkdir OUTPUT\\geometry\\boundary_condition";
        std::system(command.c_str());
    }
    file.close();

    //=============================================================================================
    for (int32_t my_id = 0; my_id < procnum; my_id++)
    {
        std::string _my_id_s;
        if (my_id < 10)
        {
            _my_id_s = "   " + std::to_string(my_id);
        }
        else if (my_id < 100)
        {
            _my_id_s = "  " + std::to_string(my_id);
        }
        else if (my_id < 1000)
        {
            _my_id_s = " " + std::to_string(my_id);
        }
        else // 这说明并行进程数不得超过9999
        {
            _my_id_s = std::to_string(my_id);
        }
        file.open("./OUTPUT/geometry/boundary_condition/inner" + _my_id_s + ".txt", std::ios_base::out);
        // 首先统计一共有多少内部关系
        int32_t count = 0;
        int32_t *count_blk = new int32_t[gp->pro_block_index(my_id).size()];
        for (int32_t iiblock = 0; iiblock < gp->pro_block_index(my_id).size(); iiblock++)
        {
            count_blk[iiblock] = 0;
            for (int32_t jjblock = 0; jjblock < gp->pro_block_index(my_id).size(); jjblock++)
            {
                Inner &inner_info = pt->inp(gp->pro_block_index(my_id)[iiblock], gp->pro_block_index(my_id)[jjblock]);
                count += inner_info.inner.size();
                count_blk[iiblock] += inner_info.inner.size();
            }
        }
        file << "\t" << count << std::endl;
        for (int32_t iiblock = 0; iiblock < gp->pro_block_index(my_id).size(); iiblock++)
            file << "\t" << count_blk[iiblock];
        file << std::endl;
        file << "imin\timax\tjmin\tjmax\tkmin\tkmax\tid\n";

        for (int32_t iiblock = 0; iiblock < gp->pro_block_index(my_id).size(); iiblock++)
        {
            for (int32_t jjblock = 0; jjblock < gp->pro_block_index(my_id).size(); jjblock++)
            {
                Inner &inner_info = pt->inp(gp->pro_block_index(my_id)[iiblock], gp->pro_block_index(my_id)[jjblock]);
                for (int32_t iface = 0; iface < inner_info.inner.size(); iface++)
                {
                    file << inner_info.inner[iface].sub[0] << "\t" << inner_info.inner[iface].sup[0] << "\t"
                         << inner_info.inner[iface].sub[1] << "\t" << inner_info.inner[iface].sup[1] << "\t"
                         << inner_info.inner[iface].sub[2] << "\t" << inner_info.inner[iface].sup[2] << "\t"
                         << iiblock + 1 << std::endl;
                    file << inner_info.inner[iface].tar_sub[0] << "\t" << inner_info.inner[iface].tar_sup[0] << "\t"
                         << inner_info.inner[iface].tar_sub[1] << "\t" << inner_info.inner[iface].tar_sup[1] << "\t"
                         << inner_info.inner[iface].tar_sub[2] << "\t" << inner_info.inner[iface].tar_sup[2] << "\t"
                         << jjblock + 1 << std::endl;
                }
            }
        }

        file.close();
        delete[] count_blk;
        std::cout << "\t\tFinish output the Inner file\t"
                  << "inner" + _my_id_s + ".txt\n";
    }
}

void OUTPUT_ACANS::output_boundary_file()
{
    std::ofstream file;
    int32_t procnum = pr->GetInt("proc_num");
    file.open("./OUTPUT/geometry/boundary_condition/boundary   0.txt", std::ios_base::out);
    if (!file)
    {
        std::string command;
        command = "mkdir OUTPUT/geometry/boundary_condition";
        std::system(command.c_str());
        command = "mkdir OUTPUT\\geometry\\boundary_condition";
        std::system(command.c_str());
    }
    file.close();

    //=============================================================================================
    for (int32_t my_id = 0; my_id < procnum; my_id++)
    {
        std::string _my_id_s;
        if (my_id < 10)
        {
            _my_id_s = "   " + std::to_string(my_id);
        }
        else if (my_id < 100)
        {
            _my_id_s = "  " + std::to_string(my_id);
        }
        else if (my_id < 1000)
        {
            _my_id_s = " " + std::to_string(my_id);
        }
        else // 这说明并行进程数不得超过9999
        {
            _my_id_s = std::to_string(my_id);
        }
        file.open("./OUTPUT/geometry/boundary_condition/boundary" + _my_id_s + ".txt", std::ios_base::out);

        for (int32_t iiblock = 0; iiblock < gp->pro_block_index(my_id).size(); iiblock++)
        {
            block_info &blk_info = pt->blk(gp->pro_block_index(my_id)[iiblock]);
            file << "\t" << blk_info.bound.size() << std::endl;
            for (int32_t iface = 0; iface < blk_info.bound.size(); iface++)
            {
                if (blk_info.bound[iface].phy_kind == 2 || (blk_info.bound[iface].phy_kind >= 200) && (blk_info.bound[iface].phy_kind <= 299))
                {
                    file << std::setw(6) << abs(blk_info.bound[iface].sub[0]) - 1 << std::setw(6) << abs(blk_info.bound[iface].sup[0]) - 1 << std::setw(6)
                         << abs(blk_info.bound[iface].sub[1]) - 1 << std::setw(6) << abs(blk_info.bound[iface].sup[1]) - 1 << std::setw(6)
                         << abs(blk_info.bound[iface].sub[2]) - 1 << std::setw(6) << abs(blk_info.bound[iface].sup[2]) - 1 << std::setw(6)
                         << blk_info.bound[iface].phy_kind << "  " << pt->phy_name[blk_info.bound[iface].phy_kind] << std::endl;
                }
                else
                {
                    file << std::setw(6) << abs(blk_info.bound[iface].sub[0]) - 1 << std::setw(6) << abs(blk_info.bound[iface].sup[0]) - 1 << std::setw(6)
                         << abs(blk_info.bound[iface].sub[1]) - 1 << std::setw(6) << abs(blk_info.bound[iface].sup[1]) - 1 << std::setw(6)
                         << abs(blk_info.bound[iface].sub[2]) - 1 << std::setw(6) << abs(blk_info.bound[iface].sup[2]) - 1 << std::setw(6)
                         << blk_info.bound[iface].phy_kind << std::endl;
                    // file << abs(blk_info.bound[iface].sub[0]) - 1 << "\t" << abs(blk_info.bound[iface].sup[0]) - 1 << "\t"
                    //  << abs(blk_info.bound[iface].sub[1]) - 1 << "\t" << abs(blk_info.bound[iface].sup[1]) - 1 << "\t"
                    //  << abs(blk_info.bound[iface].sub[2]) - 1 << "\t" << abs(blk_info.bound[iface].sup[2]) - 1 << "\t"
                    //  << blk_info.bound[iface].phy_kind << std::endl;
                }
            }
        }

        file.close();
        std::cout << "\t\tFinish output the Boundary file\t"
                  << "boundary" + _my_id_s + ".txt\n";
    }
}

void OUTPUT_ACANS::process_para_id()
{
    int32_t procnum = pr->GetInt("proc_num");
    parallel_id.SetSize(pt->blk.Getsize1(), pt->blk.Getsize1());
    for (int32_t iblock = 0; iblock < pt->blk.Getsize1(); iblock++)
    {
        for (int32_t jblock = 0; jblock < pt->blk.Getsize1(); jblock++)
        {
            parallel_id(iblock, jblock) = pt->inp(iblock, jblock);
            for (int32_t kface = 0; kface < parallel_id(iblock, jblock).inner.size(); kface++)
            {
                parallel_id(iblock, jblock).inner[kface].volume = -1;
            }
        }
    }

    // for (int32_t myid = 0; myid < pr->proc_num; myid++)
    // {
    //     for (int32_t tar_myid = 0; tar_myid < pr->proc_num; tar_myid++)
    //     {
    //         if (tar_myid == myid)
    //             continue;
    //         int32_t count = 1; //计算send_flag数量的时候选择从1开始

    //         for (int32_t iiblock = 0; iiblock < gp->proc_num(myid); iiblock++)
    //         {
    //             int32_t iblock = gp->pro_block_index(myid)[iiblock];

    //             //对于目标进程
    //             for (int32_t jjblock = 0; jjblock < gp->proc_num(tar_myid); jjblock++)
    //             {
    //                 int32_t jblock = gp->pro_block_index(tar_myid)[jjblock];
    //                 for (int32_t kface = 0; kface < parallel_id(iblock, jblock).inner.size(); kface++)
    //                 {
    //                     //计算volume作为send_flag的id
    //                     parallel_id(iblock, jblock).inner[kface].volume = myid * 1000000 + tar_myid * 1000 + count;
    //                     // //如果是周期边界条件产生的内部边界条件，flag标记为负数
    //                     // if (parallel_id(iblock, jblock).inner[kface].is_parallel == 1)
    //                     //     parallel_id(iblock, jblock).inner[kface].volume = -parallel_id(iblock, jblock).inner[kface].volume;
    //                     count++;
    //                 }
    //             }
    //         }
    //     }
    // }

    int32_t count = 1; // 计算send_flag数量的时候选择从1开始
    for (int32_t myid = 0; myid < procnum; myid++)
    {
        for (int32_t tar_myid = 0; tar_myid < procnum; tar_myid++)
        {
            if (tar_myid == myid)
                continue;

            for (int32_t iiblock = 0; iiblock < gp->proc_num(myid); iiblock++)
            {
                int32_t iblock = gp->pro_block_index(myid)[iiblock];

                // 对于目标进程
                for (int32_t jjblock = 0; jjblock < gp->proc_num(tar_myid); jjblock++)
                {
                    int32_t jblock = gp->pro_block_index(tar_myid)[jjblock];
                    for (int32_t kface = 0; kface < parallel_id(iblock, jblock).inner.size(); kface++)
                    {
                        // 计算volume作为send_flag的id
                        parallel_id(iblock, jblock).inner[kface].volume = count;
                        // //如果是周期边界条件产生的内部边界条件，flag标记为负数
                        // if (parallel_id(iblock, jblock).inner[kface].is_parallel == 1)
                        //     parallel_id(iblock, jblock).inner[kface].volume = -parallel_id(iblock, jblock).inner[kface].volume;
                        count++;
                    }
                }
            }
        }
    }
}

void OUTPUT_ACANS::output_parallel_file()
{
    int32_t procnum = pr->GetInt("proc_num");
    std::ofstream file;
    file.open("./OUTPUT/geometry/boundary_condition/parallel   0.txt", std::ios_base::out);
    if (!file)
    {
        std::string command;
        command = "mkdir OUTPUT/geometry/boundary_condition";
        std::system(command.c_str());
        command = "mkdir OUTPUT\\geometry\\boundary_condition";
        std::system(command.c_str());
    }
    file.close();
    //=============================================================================================
    process_para_id();
    //=============================================================================================
    for (int32_t my_id = 0; my_id < procnum; my_id++)
    {
        std::string _my_id_s;
        if (my_id < 10)
        {
            _my_id_s = "   " + std::to_string(my_id);
        }
        else if (my_id < 100)
        {
            _my_id_s = "  " + std::to_string(my_id);
        }
        else if (my_id < 1000)
        {
            _my_id_s = " " + std::to_string(my_id);
        }
        else // 这说明并行进程数不得超过9999
        {
            _my_id_s = std::to_string(my_id);
        }
        file.open("./OUTPUT/geometry/boundary_condition/parallel" + _my_id_s + ".txt", std::ios_base::out);
        // 首先统计一共有多少外部通信关系
        int32_t count = 0;
        int32_t *count_blk = new int32_t[gp->pro_block_index(my_id).size()];
        for (int32_t iiblock = 0; iiblock < gp->pro_block_index(my_id).size(); iiblock++)
        {
            count_blk[iiblock] = 0;
            // 针对进程中的每一块iiblock(块号为gp->pro_block_index(my_id)[iiblock])
            // 考察其他块（块号jjblock）与该块的连接关系
            for (int32_t jjblock = 0; jjblock < pt->blk.Getsize1(); jjblock++)
            {
                // 检测本块、jjblock是否被分配到同一个进程中
                if (gp->block_proc_number(jjblock) == my_id)
                {
                    // 内部通信，忽略
                    continue;
                }
                else
                {
                    Inner &inner_info = pt->inp(gp->pro_block_index(my_id)[iiblock], jjblock);
                    count += inner_info.inner.size();
                    count_blk[iiblock] += inner_info.inner.size();
                }
            }
        }
        file << "\t" << gp->proc_num(my_id) << "\t" << count << std::endl;
        for (int32_t iiblock = 0; iiblock < gp->pro_block_index(my_id).size(); iiblock++)
            file << "\t" << count_blk[iiblock];
        file << std::endl;
        file << "imin\timax\tjmin\tjmax\tkmin\tkmax\ts_r_id\ts_flag\t\tr_flag\n";

        for (int32_t iiblock = 0; iiblock < gp->pro_block_index(my_id).size(); iiblock++)
        {
            // 针对进程中的每一块iiblock(块号为gp->pro_block_index(my_id)[iiblock])
            // 考察其他块（块号jjblock）与该块的连接关系
            for (int32_t jjblock = 0; jjblock < pt->blk.Getsize1(); jjblock++)
            {
                // 检测本块、jjblock是否被分配到同一个进程中
                if (gp->block_proc_number(jjblock) == my_id)
                {
                    // 内部通信，忽略
                    continue;
                }
                else
                {
                    // 外部通信，需要输出
                    Inner &inner_info = parallel_id(gp->pro_block_index(my_id)[iiblock], jjblock);

                    for (int32_t iface = 0; iface < inner_info.inner.size(); iface++)
                    {
                        int32_t tar_face = find_index_inner_face(&parallel_id(jjblock, gp->pro_block_index(my_id)[iiblock]), inner_info.inner[iface].tar_sub, inner_info.inner[iface].tar_sup);
                        file << inner_info.inner[iface].sub[0] << "\t" << inner_info.inner[iface].sup[0] << "\t"
                             << inner_info.inner[iface].sub[1] << "\t" << inner_info.inner[iface].sup[1] << "\t"
                             << inner_info.inner[iface].sub[2] << "\t" << inner_info.inner[iface].sup[2] << "\t"
                             << gp->block_proc_number(jjblock) << "\t" << inner_info.inner[iface].volume << "\t\t"
                             << parallel_id(jjblock, gp->pro_block_index(my_id)[iiblock]).inner[tar_face].volume
                             << std::endl;
                    }
                }
            }
        }

        file.close();
        delete[] count_blk;
        std::cout << "\t\tFinish output the Parallel file\t"
                  << "parallel" + _my_id_s + ".txt\n";
    }
}

void OUTPUT_ACANS::output_gridfile()
{
    std::ofstream file;
    std::ifstream file_temp;
    int32_t procnum = pr->GetInt("proc_num");
    int32_t grd_readtype = pr->GetInt("grd_readtype");
    int32_t grd_writetype = pr->GetInt("grd_writetype");
    file.open("./OUTPUT/geometry/grid_file/grid   0.grd", std::ios_base::out);
    if (!file)
    {
        std::string command;
        command = "mkdir OUTPUT/geometry/grid_file";
        std::system(command.c_str());
        command = "mkdir OUTPUT\\geometry\\grid_file";
        std::system(command.c_str());
    }
    file.close();
    //=============================================================================================

    for (int32_t my_id = 0; my_id < procnum; my_id++)
    {
        std::string _my_id_s;
        if (my_id < 10)
        {
            _my_id_s = "   " + std::to_string(my_id);
        }
        else if (my_id < 100)
        {
            _my_id_s = "  " + std::to_string(my_id);
        }
        else if (my_id < 1000)
        {
            _my_id_s = " " + std::to_string(my_id);
        }
        else // 这说明并行进程数不得超过9999
        {
            _my_id_s = std::to_string(my_id);
        }

        if (grd_writetype == 0)
        {
            file.open("./OUTPUT/geometry/grid_file/grid" + _my_id_s + ".grd", std::ios_base::out);
            /******输出网格块总数******/
            file << gp->proc_num(my_id) << std::endl;
            // /******输出网格块名称******/
            // // Remark: 原始ACANS以及plot3D格式的网格文件中没有这一部分
            // for (int32_t i = 0; i < gp->proc_num(my_id); i++)
            //     file << pt->blk(gp->pro_block_index(my_id)[i]).hiera[0].blk_name << std::endl;
            /******输出网格块点数******/
            for (int32_t i = 0; i < gp->proc_num(my_id); i++)
                file << pt->blk(gp->pro_block_index(my_id)[i]).ijkmax[0] << "\t"
                     << pt->blk(gp->pro_block_index(my_id)[i]).ijkmax[1] << "\t"
                     << pt->blk(gp->pro_block_index(my_id)[i]).ijkmax[2] << std::endl;
            /******输出网格坐标******/
        }
        else
        {
            file.open("./OUTPUT/geometry/grid_file/grid" + _my_id_s + ".grd", std::ios::out | std::ios::binary);
            int32_t temp;
            temp = 4;

            if (grd_writetype == 1)
            {
                file.write((char *)&temp, sizeof(temp));
            }

            file.write((char *)&gp->proc_num(my_id), sizeof(gp->proc_num(my_id)));

            if (grd_writetype == 1)
            {
                file.write((char *)&temp, sizeof(temp));
                temp = gp->proc_num(my_id) * 12;
                file.write((char *)&temp, sizeof(temp));
            }

            for (int i = 0; i < gp->proc_num(my_id); i++)
            {

                int idim = pt->blk(gp->pro_block_index(my_id)[i]).ijkmax[0];
                int jdim = pt->blk(gp->pro_block_index(my_id)[i]).ijkmax[1];
                int kdim = pt->blk(gp->pro_block_index(my_id)[i]).ijkmax[2];

                file.write((char *)&idim, sizeof(idim));
                file.write((char *)&jdim, sizeof(jdim));
                file.write((char *)&kdim, sizeof(kdim));
            }

            if (grd_writetype == 1)
            {
                file.write((char *)&temp, sizeof(temp));
            }
        }
        for (int32_t iblock = 0; iblock < gp->proc_num(my_id); iblock++)
        {
            block_info &iblk = pt->blk(gp->pro_block_index(my_id)[iblock]);
            double3D x(iblk.ijkmax[0], iblk.ijkmax[1], iblk.ijkmax[2], -1);
            double3D y(iblk.ijkmax[0], iblk.ijkmax[1], iblk.ijkmax[2], -1);
            double3D z(iblk.ijkmax[0], iblk.ijkmax[1], iblk.ijkmax[2], -1);
            for (int32_t sub_blk = 0; sub_blk < iblk.hiera.size(); sub_blk++)
            {
                int32_t origin_blk_num = iblk.hiera[sub_blk].oring_num;
                int32_t numofblock, imax, jmax, kmax;
                if (grd_readtype == 1 || grd_readtype == 2)
                {

                    //==================================================
                    // 读入原始网格文件,尝试读入gridgen的unformatted格式
                    file_temp.open(pr->GetStr("gfilename"), std::ios::in | std::ios::binary);
                    int32_t temp;
                    if (grd_readtype == 1)
                    {
                        file_temp.read((char *)&temp, sizeof(temp));
                    }
                    file_temp.read((char *)&temp, sizeof(temp));
                    numofblock = temp;
                    if (grd_readtype == 1)
                    {
                        file_temp.read((char *)&temp, sizeof(temp));
                        file_temp.read((char *)&temp, sizeof(temp));
                    }
                    // file_temp.close();
                    // std::cout << numofblock<<std::endl;

                    int32_t *mx = new int32_t[numofblock];
                    int32_t *my = new int32_t[numofblock];
                    int32_t *mz = new int32_t[numofblock];
                    for (int iiblock = 0; iiblock < numofblock; iiblock++)
                    {
                        file_temp.read((char *)&temp, sizeof(temp));
                        mx[iiblock] = temp;
                        file_temp.read((char *)&temp, sizeof(temp));
                        my[iiblock] = temp;
                        file_temp.read((char *)&temp, sizeof(temp));
                        mz[iiblock] = temp;
                        // std::cout << " "<<mx[iiblock] << " " << my[iiblock] << " " << mz[iiblock]<<std::endl;
                    }
                    imax = mx[origin_blk_num];
                    jmax = my[origin_blk_num];
                    kmax = mz[origin_blk_num];

                    if (grd_readtype == 1)
                    {
                        file_temp.read((char *)&temp, sizeof(temp));
                    }
                    // std::cout << std::endl;
                    double3D origin_x(imax, jmax, kmax, -1);
                    double3D origin_y(imax, jmax, kmax, -1);
                    double3D origin_z(imax, jmax, kmax, -1);
                    for (int iiblock = 0; iiblock < numofblock; iiblock++)
                    {
                        double temp_double;
                        double temp_binary;
                        if (grd_readtype == 1)
                        {
                            file_temp.read((char *)&temp, sizeof(temp));
                        }
                        if (iiblock == origin_blk_num)
                        {
                            for (int k = 0; k < mz[iiblock]; ++k)
                                for (int j = 0; j < my[iiblock]; ++j)
                                    for (int i = 0; i < mx[iiblock]; ++i)
                                    {
                                        file_temp.read((char *)&temp_binary, sizeof(temp_binary));
                                        origin_x(i + 1, j + 1, k + 1) = temp_binary;
                                        // std::cout << " " << temp_binary << std::endl;
                                    }
                            for (int k = 0; k < mz[iiblock]; ++k)
                                for (int j = 0; j < my[iiblock]; ++j)
                                    for (int i = 0; i < mx[iiblock]; ++i)
                                    {
                                        file_temp.read((char *)&temp_binary, sizeof(temp_binary));
                                        origin_y(i + 1, j + 1, k + 1) = temp_binary;
                                        // std::cout << " " << temp_binary << std::endl;
                                    }
                            for (int k = 0; k < mz[iiblock]; ++k)
                                for (int j = 0; j < my[iiblock]; ++j)
                                    for (int i = 0; i < mx[iiblock]; ++i)
                                    {
                                        file_temp.read((char *)&temp_binary, sizeof(temp_binary));
                                        origin_z(i + 1, j + 1, k + 1) = temp_binary;
                                        // std::cout << " " << temp_binary << std::endl;
                                    }
                            break;
                        }
                        else
                        {
                            for (int k = 0; k < mz[iiblock]; ++k)
                                for (int j = 0; j < my[iiblock]; ++j)
                                    for (int i = 0; i < mx[iiblock]; ++i)
                                    {
                                        file_temp.read((char *)&temp_binary, sizeof(temp_binary));
                                        temp_double = temp_binary;
                                        // std::cout << " " << temp_binary << std::endl;
                                    }
                            for (int k = 0; k < mz[iiblock]; ++k)
                                for (int j = 0; j < my[iiblock]; ++j)
                                    for (int i = 0; i < mx[iiblock]; ++i)
                                    {
                                        file_temp.read((char *)&temp_binary, sizeof(temp_binary));
                                        temp_double = temp_binary;
                                        // std::cout << " " << temp_binary << std::endl;
                                    }
                            for (int k = 0; k < mz[iiblock]; ++k)
                                for (int j = 0; j < my[iiblock]; ++j)
                                    for (int i = 0; i < mx[iiblock]; ++i)
                                    {
                                        file_temp.read((char *)&temp_binary, sizeof(temp_binary));
                                        temp_double = temp_binary;
                                        // std::cout << " " << temp_binary << std::endl;
                                    }
                        }
                        if (grd_readtype == 1)
                        {
                            file_temp.read((char *)&temp, sizeof(temp));
                        }
                    }
                    file_temp.close();
                    delete[] mx;
                    delete[] my;
                    delete[] mz;
                    // 读入原始网格文件结束
                    //==================================================
                    int32_t ijk[3], sub[3], sup[3], tar_sub[3], tar_sup[3], Transform[3];
                    for (int i = 0; i < 3; ++i)
                    {
                        sub[i] = iblk.hiera[sub_blk].sub[i];
                        sup[i] = iblk.hiera[sub_blk].sup[i];
                        tar_sub[i] = iblk.hiera[sub_blk].oring_sub[i];
                        tar_sup[i] = iblk.hiera[sub_blk].oring_sup[i];
                        Transform[i] = iblk.hiera[sub_blk].Transform[i];
                    }
                    for (int k = sub[2]; k <= sup[2]; ++k)
                    {
                        for (int j = sub[1]; j <= sup[1]; ++j)
                        {
                            for (int i = sub[0]; i <= sup[0]; ++i)
                            {
                                ijk[abs(Transform[0]) - 1] = (Transform[0] > 0) ? (i - sub[0]) + tar_sub[abs(Transform[0]) - 1] : tar_sup[abs(Transform[0]) - 1] - (i - sub[0]);
                                ijk[abs(Transform[1]) - 1] = (Transform[1] > 0) ? (j - sub[1]) + tar_sub[abs(Transform[1]) - 1] : tar_sup[abs(Transform[1]) - 1] - (j - sub[1]);
                                ijk[abs(Transform[2]) - 1] = (Transform[2] > 0) ? (k - sub[2]) + tar_sub[abs(Transform[2]) - 1] : tar_sup[abs(Transform[2]) - 1] - (k - sub[2]);
                                x(i, j, k) = origin_x(ijk[0], ijk[1], ijk[2]);
                                y(i, j, k) = origin_y(ijk[0], ijk[1], ijk[2]);
                                z(i, j, k) = origin_z(ijk[0], ijk[1], ijk[2]);
                            }
                        }
                    }
                }
                else if (grd_readtype == 0)
                {
                    //==================================================
                    //    读入原始网格文件
                    file_temp.open(pr->GetStr("gfilename"), std::ios_base::in);
                    file_temp >> numofblock;
                    int32_t *mx = new int32_t[numofblock];
                    int32_t *my = new int32_t[numofblock];
                    int32_t *mz = new int32_t[numofblock];
                    for (int iiblock = 0; iiblock < numofblock; iiblock++)
                        file_temp >> mx[iiblock] >> my[iiblock] >> mz[iiblock];
                    imax = mx[origin_blk_num];
                    jmax = my[origin_blk_num];
                    kmax = mz[origin_blk_num];

                    double3D origin_x(imax, jmax, kmax, -1);
                    double3D origin_y(imax, jmax, kmax, -1);
                    double3D origin_z(imax, jmax, kmax, -1);
                    for (int iiblock = 0; iiblock < numofblock; iiblock++)
                    {
                        double temp_double;
                        if (iiblock == origin_blk_num)
                        {
                            for (int k = 0; k < mz[iiblock]; ++k)
                                for (int j = 0; j < my[iiblock]; ++j)
                                    for (int i = 0; i < mx[iiblock]; ++i)
                                        file_temp >> origin_x(i + 1, j + 1, k + 1);
                            for (int k = 0; k < mz[iiblock]; ++k)
                                for (int j = 0; j < my[iiblock]; ++j)
                                    for (int i = 0; i < mx[iiblock]; ++i)
                                        file_temp >> origin_y(i + 1, j + 1, k + 1);
                            for (int k = 0; k < mz[iiblock]; ++k)
                                for (int j = 0; j < my[iiblock]; ++j)
                                    for (int i = 0; i < mx[iiblock]; ++i)
                                        file_temp >> origin_z(i + 1, j + 1, k + 1);
                            break;
                        }
                        else
                        {
                            for (int k = 0; k < mz[iiblock]; ++k)
                                for (int j = 0; j < my[iiblock]; ++j)
                                    for (int i = 0; i < mx[iiblock]; ++i)
                                        file_temp >> temp_double;
                            for (int k = 0; k < mz[iiblock]; ++k)
                                for (int j = 0; j < my[iiblock]; ++j)
                                    for (int i = 0; i < mx[iiblock]; ++i)
                                        file_temp >> temp_double;
                            for (int k = 0; k < mz[iiblock]; ++k)
                                for (int j = 0; j < my[iiblock]; ++j)
                                    for (int i = 0; i < mx[iiblock]; ++i)
                                        file_temp >> temp_double;
                        }
                    }
                    file_temp.close();
                    delete[] mx;
                    delete[] my;
                    delete[] mz;
                    // 读入原始网格文件结束
                    //==================================================
                    int32_t ijk[3], sub[3], sup[3], tar_sub[3], tar_sup[3], Transform[3];
                    for (int i = 0; i < 3; ++i)
                    {
                        sub[i] = iblk.hiera[sub_blk].sub[i];
                        sup[i] = iblk.hiera[sub_blk].sup[i];
                        tar_sub[i] = iblk.hiera[sub_blk].oring_sub[i];
                        tar_sup[i] = iblk.hiera[sub_blk].oring_sup[i];
                        Transform[i] = iblk.hiera[sub_blk].Transform[i];
                    }

                    for (int k = sub[2]; k <= sup[2]; ++k)
                    {
                        for (int j = sub[1]; j <= sup[1]; ++j)
                        {
                            for (int i = sub[0]; i <= sup[0]; ++i)
                            {
                                ijk[abs(Transform[0]) - 1] = (Transform[0] > 0) ? (i - sub[0]) + tar_sub[abs(Transform[0]) - 1] : tar_sup[abs(Transform[0]) - 1] - (i - sub[0]);
                                ijk[abs(Transform[1]) - 1] = (Transform[1] > 0) ? (j - sub[1]) + tar_sub[abs(Transform[1]) - 1] : tar_sup[abs(Transform[1]) - 1] - (j - sub[1]);
                                ijk[abs(Transform[2]) - 1] = (Transform[2] > 0) ? (k - sub[2]) + tar_sub[abs(Transform[2]) - 1] : tar_sup[abs(Transform[2]) - 1] - (k - sub[2]);
                                x(i, j, k) = origin_x(ijk[0], ijk[1], ijk[2]);
                                y(i, j, k) = origin_y(ijk[0], ijk[1], ijk[2]);
                                z(i, j, k) = origin_z(ijk[0], ijk[1], ijk[2]);
                            }
                        }
                    }
                }
                //
                //
                ////将数组origin_xyz中的数据存入xyz中
                // int32_t ijk[3], sub[3], sup[3], tar_sub[3], tar_sup[3], Transform[3];
                // for (int i = 0; i < 3; ++i)
                //{
                //     sub[i] = iblk.hiera[sub_blk].sub[i];
                //     sup[i] = iblk.hiera[sub_blk].sup[i];
                //     tar_sub[i] = iblk.hiera[sub_blk].oring_sub[i];
                //     tar_sup[i] = iblk.hiera[sub_blk].oring_sup[i];
                //     Transform[i] = iblk.hiera[sub_blk].Transform[i];
                // }
                // for (int k = sub[2]; k <= sup[2]; ++k)
                //{
                //     for (int j = sub[1]; j <= sup[1]; ++j)
                //     {
                //         for (int i = sub[0]; i <= sup[0]; ++i)
                //         {
                //             ijk[abs(Transform[0]) - 1] = (Transform[0] > 0) ? (i - sub[0]) + tar_sub[abs(Transform[0]) - 1] : tar_sup[abs(Transform[0]) - 1] - (i - sub[0]);
                //             ijk[abs(Transform[1]) - 1] = (Transform[1] > 0) ? (j - sub[1]) + tar_sub[abs(Transform[1]) - 1] : tar_sup[abs(Transform[1]) - 1] - (j - sub[1]);
                //             ijk[abs(Transform[2]) - 1] = (Transform[2] > 0) ? (k - sub[2]) + tar_sub[abs(Transform[2]) - 1] : tar_sup[abs(Transform[2]) - 1] - (k - sub[2]);
                //             x(i, j, k) = origin_x(ijk[0], ijk[1], ijk[2]);
                //             y(i, j, k) = origin_y(ijk[0], ijk[1], ijk[2]);
                //             z(i, j, k) = origin_z(ijk[0], ijk[1], ijk[2]);
                //         }
                //     }
                // }
                ////将本子网格块的origin_xyz中的数据存入了xyz中
                //==================================================
            }
            //=======================================================
            // 输出
            if (grd_writetype == 0)
            {
                int32_t count = 0;
                file << std::setiosflags(std::ios::scientific) << std::setiosflags(std::ios::right) << std::setprecision(17) << std::setw(30);
                for (int k = 1; k <= iblk.ijkmax[2]; ++k)
                {
                    for (int j = 1; j <= iblk.ijkmax[1]; ++j)
                    {
                        for (int i = 1; i <= iblk.ijkmax[0]; ++i)
                        {

                            if (count == 3)
                            {
                                file << std::setiosflags(std::ios::scientific) << std::setw(28) << std::setprecision(18) << x(i, j, k) << std::endl;
                                count = 0;
                                continue;
                            }
                            file << std::setiosflags(std::ios::scientific) << std::setw(28) << std::setprecision(18) << x(i, j, k);
                            count++;
                        }
                    }
                }

                for (int k = 1; k <= iblk.ijkmax[2]; ++k)
                {
                    for (int j = 1; j <= iblk.ijkmax[1]; ++j)
                    {
                        for (int i = 1; i <= iblk.ijkmax[0]; ++i)
                        {

                            if (count == 3)
                            {
                                file << std::setiosflags(std::ios::scientific) << std::setw(28) << std::setprecision(18) << y(i, j, k) << std::endl;
                                count = 0;
                                continue;
                            }
                            file << std::setiosflags(std::ios::scientific) << std::setw(28) << std::setprecision(18) << y(i, j, k);
                            count++;
                        }
                    }
                }

                for (int k = 1; k <= iblk.ijkmax[2]; ++k)
                {
                    for (int j = 1; j <= iblk.ijkmax[1]; ++j)
                    {
                        for (int i = 1; i <= iblk.ijkmax[0]; ++i)
                        {

                            if (count == 3)
                            {
                                file << std::setiosflags(std::ios::scientific) << std::setw(28) << std::setprecision(18) << z(i, j, k) << std::endl;
                                count = 0;
                                continue;
                            }
                            file << std::setiosflags(std::ios::scientific) << std::setw(28) << std::setprecision(18) << z(i, j, k);
                            count++;
                        }
                    }
                }
                file << std::endl;
            }
            else
            {

                int32_t dtemp;
                dtemp = iblk.ijkmax[2] * iblk.ijkmax[1] * iblk.ijkmax[0] * 24;
                if (grd_writetype == 1)
                {
                    file.write((char *)&dtemp, sizeof(dtemp));
                }

                for (int k = 1; k <= iblk.ijkmax[2]; ++k)
                {
                    for (int j = 1; j <= iblk.ijkmax[1]; ++j)
                    {
                        for (int i = 1; i <= iblk.ijkmax[0]; ++i)
                        {
                            file.write((char *)&x(i, j, k), sizeof(x(i, j, k)));
                        }
                    }
                }
                for (int k = 1; k <= iblk.ijkmax[2]; ++k)
                {
                    for (int j = 1; j <= iblk.ijkmax[1]; ++j)
                    {
                        for (int i = 1; i <= iblk.ijkmax[0]; ++i)
                        {
                            file.write((char *)&y(i, j, k), sizeof(y(i, j, k)));
                        }
                    }
                }

                for (int k = 1; k <= iblk.ijkmax[2]; ++k)
                {
                    for (int j = 1; j <= iblk.ijkmax[1]; ++j)
                    {
                        for (int i = 1; i <= iblk.ijkmax[0]; ++i)
                        {
                            file.write((char *)&z(i, j, k), sizeof(z(i, j, k)));
                        }
                    }
                }

                if (grd_writetype == 1)
                {
                    file.write((char *)&dtemp, sizeof(dtemp));
                }
            }
        }
        file.close();
        std::cout << "\t\tFinish output the Grid file\t"
                  << "grid" + _my_id_s + ".grd\n";
    }
}
