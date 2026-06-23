#include "MPCNS_Pre_Split.h"
#include <math.h>
// #include "preprocess_output.h"
void Split(Preprocess_Data_Structure *ptr, Param *par)
{
    if (!par->GetBoo("ifsplit"))
        return;
    //=============================================================================================
    // 开始自动剖分
    //---------------------------------------------------------------------------------------------
    bool stop_the_loop = true;
    int32_t proc_num = par->GetInt("proc_num");
    int32_t dimension = par->GetInt("dimension");
    double split_scale = par->GetDou("split_scale");
    int32_t averge_grid = 0, max_grid = 0, this_block_number = 0;
    do
    {
        stop_the_loop = true;
        // 统计此时的平均网格数量
        averge_grid = 0;
        for (int32_t iblock = 0; iblock < ptr->blk.Getsize1(); iblock++)
            averge_grid += ptr->blk(iblock).ijkmax[0] * ptr->blk(iblock).ijkmax[1] * ptr->blk(iblock).ijkmax[2];
        averge_grid = (int32_t)(averge_grid / proc_num);
        // 允许的最大网格数量，凡是大于该网格量的都剖分
        max_grid = averge_grid * (split_scale + 1.0);
        //---------------------------------------------------------------------------
        int32_t split_index;
        for (int32_t iblock = 0; iblock < ptr->blk.Getsize1(); iblock++)
        {
            if (ptr->blk(iblock).ijkmax[0] * ptr->blk(iblock).ijkmax[1] * ptr->blk(iblock).ijkmax[2] > max_grid)
            {
                this_block_number = ptr->blk(iblock).ijkmax[0] * ptr->blk(iblock).ijkmax[1] * ptr->blk(iblock).ijkmax[2];
                split_index = iblock;
                stop_the_loop = false;
                break;
            }
        }
        if (stop_the_loop)
            break;
        //---------------------------------------------------------------------------
        if (this_block_number < 2 * max_grid)
        {
            // 1~2倍平均网格点数则《尝试》采用最大维中点二分
            split_mid_location(ptr, par, split_index);
        }
        else if (this_block_number < 4 * max_grid)
        {
            // 2~4倍平均网格点数则《尝试》采用一维剖分，返回剖分方向和位置
            split_1D_location(ptr, par, split_index, max_grid);
        }
        else if (this_block_number < 8 * max_grid || dimension == 2)
        {
            // 4~8倍平均网格点数则《尝试》采用二维剖分，子进程中直接修改剖分方向与位置
            split_2D_location(ptr, par, split_index, max_grid);
        }
        else if (dimension == 3)
        {
            // 8倍以上的平均网格点数切维数为3则《尝试》采用三维剖分，子进程中直接修改剖分方向与位置
            split_3D_location(ptr, par, split_index, max_grid);
        }
        std::cout << "\t\t->We split the grid " << split_index << "\tNow the number is " << ptr->blk.Getsize1() << std::endl;
        ptr->check_read_info(par);
    } while (!stop_the_loop);

    std::cout << "\t\tThe average grid number is: " << averge_grid << std::endl;
    ptr->check_read_info(par);
    output_split_info_inp(ptr, par);
}

void split_mid_location(Preprocess_Data_Structure *ptr, Param *par, int32_t split_index)
{
    int32_t ijkmax = 0, direction = 0, location = 0;
    for (int32_t i = 0; i < 3; i++)
    {
        if (ptr->blk(split_index).ijkmax[i] > ijkmax)
        {
            ijkmax = ptr->blk(split_index).ijkmax[i];
            direction = i;
        }
    }

    location = (int)(ptr->blk(split_index).ijkmax[direction] / 2);
    direction++;
    ptr->split_block(par, split_index, direction, location);
}

void split_1D_location(Preprocess_Data_Structure *ptr, Param *par, int32_t split_index, int32_t average_grid)
{
    int32_t ijkmax = 0, direction = 0, location = 0;
    for (int32_t i = 0; i < 3; i++)
    {
        if (ptr->blk(split_index).ijkmax[i] > ijkmax)
        {
            ijkmax = ptr->blk(split_index).ijkmax[i];
            direction = i;
        }
    }

    location = (int)(ptr->blk(split_index).ijkmax[direction] * average_grid / ptr->blk(split_index).ijkmax[0] / ptr->blk(split_index).ijkmax[1] / ptr->blk(split_index).ijkmax[2]);
    direction++;
    ptr->split_block(par, split_index, direction, location);
}

void split_2D_location(Preprocess_Data_Structure *ptr, Param *par, int32_t split_index, int32_t average_grid)
{
    int32_t ijkmax = 1e7, direction = 0, location = 0;
    int32_t split_limit = par->GetInt("split_limit");
    for (int32_t i = 0; i < 3; i++)
    {
        if (ptr->blk(split_index).ijkmax[i] < ijkmax)
        {
            ijkmax = ptr->blk(split_index).ijkmax[i];
            direction = i;
        }
    }
    location = (int)sqrt(average_grid / ptr->blk(split_index).ijkmax[direction]);

    if (location < split_limit)
    {
        direction++;
        split_1D_location(ptr, par, split_index, average_grid);
        return;
    }

    for (int32_t i = 0; i < 3; i++)
    {
        if (i == direction)
            continue;
        if (abs(ptr->blk(split_index).ijkmax[i] - location) < split_limit || location > ptr->blk(split_index).ijkmax[i])
        {
            direction++;
            split_1D_location(ptr, par, split_index, average_grid);
            return;
        }
    }

    // 一剖为4
    int32_t old_number_max = ptr->blk.Getsize1();
    direction++;
    if (direction == 1)
    {
        ptr->split_block(par, split_index, 2, location);
        ptr->split_block(par, old_number_max, 3, location);
        ptr->split_block(par, split_index, 3, location);
    }
    else if (direction == 2)
    {
        ptr->split_block(par, split_index, 1, location);
        ptr->split_block(par, old_number_max, 3, location);
        ptr->split_block(par, split_index, 3, location);
    }
    else if (direction == 3)
    {
        ptr->split_block(par, split_index, 1, location);
        ptr->split_block(par, old_number_max, 2, location);
        ptr->split_block(par, split_index, 2, location);
    }
}

void split_3D_location(Preprocess_Data_Structure *ptr, Param *par, int32_t split_index, int32_t average_grid)
{
    int32_t direction = 0, location = 0;
    int32_t split_limit = par->GetInt("split_limit");

    location = (int)pow(average_grid, 1.0 / 3.0);

    if (location < split_limit)
    {
        split_2D_location(ptr, par, split_index, average_grid);
        return;
    }

    for (int32_t i = 0; i < 3; i++)
    {
        if (abs(ptr->blk(split_index).ijkmax[i] - location) < split_limit || location > ptr->blk(split_index).ijkmax[i])
        {
            split_2D_location(ptr, par, split_index, average_grid);
            return;
        }
    }

    // 一剖为8
    int32_t old_number_max = ptr->blk.Getsize1();
    ptr->split_block(par, split_index, 1, location);
    ptr->split_block(par, split_index, 2, location);
    ptr->split_block(par, split_index, 3, location);

    ptr->split_block(par, old_number_max, 2, location);
    ptr->split_block(par, old_number_max, 3, location);

    ptr->split_block(par, old_number_max + 1, 3, location);
    ptr->split_block(par, old_number_max + 3, 3, location);
}

void output_split_info_inp(Preprocess_Data_Structure *ptr, Param *par)
{
    std::ofstream file;
    file.open("./OUTPUT/split_info.inp", std::ios_base::out);
    if (!file)
    {
        std::string command;
        command = "mkdir OUTPUT&";
        std::system(command.c_str());
        file.open("./OUTPUT/split_info.inp", std::ios_base::out);
    }
    file << "\t1\n";
    file << "\t" << ptr->blk.Getsize1() << "\n";

    for (int32_t i = 0; i < ptr->blk.Getsize1(); i++)
    {
        file << "\t" << ptr->blk(i).ijkmax[0] << "\t" << ptr->blk(i).ijkmax[1] << "\t" << ptr->blk(i).ijkmax[2] << std::endl;
        file << ptr->blk(i).hiera[0].blk_name << "\n";
        int32_t count = 0;
        count = ptr->blk(i).bound.size();
        for (int32_t j = 0; j < ptr->blk.Getsize1(); j++)
            count += ptr->inp(i, j).inner.size();

        file << "\t" << count << "\n";
        for (int32_t j = 0; j < ptr->blk(i).bound.size(); j++)
        {
            file << "\t" << ptr->blk(i).bound[j].sub[0] << "\t" << ptr->blk(i).bound[j].sup[0]
                 << "\t" << ptr->blk(i).bound[j].sub[1] << "\t" << ptr->blk(i).bound[j].sup[1]
                 << "\t" << ptr->blk(i).bound[j].sub[2] << "\t" << ptr->blk(i).bound[j].sup[2]
                 << "\t" << ptr->blk(i).bound[j].phy_kind << "\n";
        }

        for (int32_t j = 0; j < ptr->blk.Getsize1(); j++)
        {
            for (int32_t k = 0; k < ptr->inp(i, j).inner.size(); k++)
            {
                coordinate_inp &temp_inp = ptr->inp(i, j).inner[k];
                file << "\t" << temp_inp.sub[0] << "\t" << temp_inp.sup[0]
                     << "\t" << temp_inp.sub[1] << "\t" << temp_inp.sup[1]
                     << "\t" << temp_inp.sub[2] << "\t" << temp_inp.sup[2] << "\t-1\n";
                if (temp_inp.is_parallel == 0) // 一般的内部边界条件
                {
                    file << "\t" << temp_inp.tar_sub[0] << "\t" << temp_inp.tar_sup[0]
                         << "\t" << temp_inp.tar_sub[1] << "\t" << temp_inp.tar_sup[1]
                         << "\t" << temp_inp.tar_sub[2] << "\t" << temp_inp.tar_sup[2]
                         << "\t" << j + 1 << std::endl;
                }
                else // 周期边界条件
                {
                    file << "\t" << temp_inp.tar_sub[0] << "\t" << temp_inp.tar_sup[0]
                         << "\t" << temp_inp.tar_sub[1] << "\t" << temp_inp.tar_sup[1]
                         << "\t" << temp_inp.tar_sub[2] << "\t" << temp_inp.tar_sup[2]
                         << "\t" << -j - 1 << std::endl;
                }
            }
        }
    }
    file.close();
}