#include "MPCNS_Pre_Split.h"
#include <math.h>
// #include "preprocess_output.h"

namespace
{
int32_t GetPoleAzimuthDirection(Param *par)
{
    if (!par->HasInt("pole_azimuth_direction"))
        return 3;
    return par->GetInt("pole_azimuth_direction");
}

std::string GetPoleBoundaryName(Param *par)
{
    if (!par->HasStr("pole_boundary_name"))
        return "Pole";
    return par->GetStr("pole_boundary_name");
}

bool ForbidPoleAzimuthSplit(Param *par)
{
    if (!par->HasBoo("forbid_pole_azimuth_split"))
        return true;
    return par->GetBoo("forbid_pole_azimuth_split");
}

int32_t FaceNormalDirection(const coordinate_phy &face, int32_t dimension)
{
    int32_t normal = 0;
    for (int32_t i = 0; i < dimension; i++)
    {
        if (face.sub[i] == face.sup[i])
        {
            if (normal != 0)
                return 0;
            normal = i + 1;
        }
    }
    return normal;
}

bool HasPoleFaceCutByAzimuth(Preprocess_Data_Structure *ptr, Param *par, int32_t split_index)
{
    if (!ForbidPoleAzimuthSplit(par))
        return false;

    const std::string pole_boundary_name = GetPoleBoundaryName(par);
    const std::map<std::string, int32_t>::iterator it = ptr->phy_index.find(pole_boundary_name);
    if (it == ptr->phy_index.end())
        return false;

    const int32_t pole_kind = it->second;
    const int32_t azimuth_direction = GetPoleAzimuthDirection(par);
    const int32_t dimension = par->GetInt("dimension");

    for (int32_t iface = 0; iface < ptr->blk(split_index).bound.size(); iface++)
    {
        const coordinate_phy &face = ptr->blk(split_index).bound[iface];
        if (face.phy_kind != pole_kind)
            continue;

        const int32_t normal_direction = FaceNormalDirection(face, dimension);
        if (normal_direction == azimuth_direction)
        {
            std::cout << "\tWarning!\tPole boundary '" << pole_boundary_name
                      << "' has normal direction equal to pole_azimuth_direction = "
                      << azimuth_direction << " on block " << split_index << ".\n";
            continue;
        }
        return true;
    }
    return false;
}

bool IsSplitDirectionForbidden(Preprocess_Data_Structure *ptr, Param *par, int32_t split_index, int32_t direction)
{
    return direction == GetPoleAzimuthDirection(par) && HasPoleFaceCutByAzimuth(ptr, par, split_index);
}

void ReportPoleBoundaryMapping(Preprocess_Data_Structure *ptr, Param *par)
{
    if (!ForbidPoleAzimuthSplit(par))
        return;

    const std::string pole_boundary_name = GetPoleBoundaryName(par);
    const std::map<std::string, int32_t>::iterator it = ptr->phy_index.find(pole_boundary_name);
    if (it == ptr->phy_index.end())
    {
        std::cout << "\tWarning!\tCannot find pole_boundary_name = " << pole_boundary_name
                  << " in fvbnd boundary names. Pole azimuth split constraint is inactive.\n";
        return;
    }

    std::cout << "\t\tPole boundary '" << pole_boundary_name << "' corresponds to physical kind "
              << it->second << ".\n";
}

bool SplitBlockIfAllowed(Preprocess_Data_Structure *ptr, Param *par, int32_t split_index, int32_t direction, int32_t location)
{
    if (IsSplitDirectionForbidden(ptr, par, split_index, direction))
        return false;
    ptr->split_block(par, split_index, direction, location);
    return true;
}

int32_t LargestAllowedDirection(Preprocess_Data_Structure *ptr, Param *par, int32_t split_index)
{
    int32_t direction = 0;
    int32_t ijkmax = 0;
    for (int32_t i = 0; i < 3; i++)
    {
        if (IsSplitDirectionForbidden(ptr, par, split_index, i + 1))
            continue;
        if (ptr->blk(split_index).ijkmax[i] > ijkmax)
        {
            ijkmax = ptr->blk(split_index).ijkmax[i];
            direction = i + 1;
        }
    }
    return direction;
}
}

void Split(Preprocess_Data_Structure *ptr, Param *par)
{
    if (!par->GetBoo("ifsplit"))
        return;
    ReportPoleBoundaryMapping(ptr, par);
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
        bool split_success = false;
        if (this_block_number < 2 * max_grid)
        {
            // 1~2倍平均网格点数则《尝试》采用最大维中点二分
            split_success = split_mid_location(ptr, par, split_index);
        }
        else if (this_block_number < 4 * max_grid)
        {
            // 2~4倍平均网格点数则《尝试》采用一维剖分，返回剖分方向和位置
            split_success = split_1D_location(ptr, par, split_index, max_grid);
        }
        else if (this_block_number < 8 * max_grid || dimension == 2)
        {
            // 4~8倍平均网格点数则《尝试》采用二维剖分，子进程中直接修改剖分方向与位置
            split_success = split_2D_location(ptr, par, split_index, max_grid);
        }
        else if (dimension == 3)
        {
            // 8倍以上的平均网格点数切维数为3则《尝试》采用三维剖分，子进程中直接修改剖分方向与位置
            split_success = split_3D_location(ptr, par, split_index, max_grid);
        }
        if (!split_success)
        {
            std::cout << "\tWarning!\tNo legal split direction for block " << split_index
                      << " under Pole azimuth split constraint. Stop splitting.\n";
            break;
        }
        std::cout << "\t\t->We split the grid " << split_index << "\tNow the number is " << ptr->blk.Getsize1() << std::endl;
        ptr->check_read_info(par);
    } while (!stop_the_loop);

    std::cout << "\t\tThe average grid number is: " << averge_grid << std::endl;
    ptr->check_read_info(par);
    output_split_info_inp(ptr, par);
}

bool split_mid_location(Preprocess_Data_Structure *ptr, Param *par, int32_t split_index)
{
    int32_t direction = LargestAllowedDirection(ptr, par, split_index);
    if (direction == 0)
        return false;
    int32_t location = 0;
    direction--;
    location = (int)(ptr->blk(split_index).ijkmax[direction] / 2);
    direction++;
    return SplitBlockIfAllowed(ptr, par, split_index, direction, location);
}

bool split_1D_location(Preprocess_Data_Structure *ptr, Param *par, int32_t split_index, int32_t average_grid)
{
    int32_t direction = LargestAllowedDirection(ptr, par, split_index);
    if (direction == 0)
        return false;
    int32_t location = 0;
    direction--;
    location = (int)(ptr->blk(split_index).ijkmax[direction] * average_grid / ptr->blk(split_index).ijkmax[0] / ptr->blk(split_index).ijkmax[1] / ptr->blk(split_index).ijkmax[2]);
    direction++;
    return SplitBlockIfAllowed(ptr, par, split_index, direction, location);
}

bool split_2D_location(Preprocess_Data_Structure *ptr, Param *par, int32_t split_index, int32_t average_grid)
{
    int32_t ijkmax = 1e7, direction = 0, location = 0;
    int32_t split_limit = par->GetInt("split_limit");
    if (HasPoleFaceCutByAzimuth(ptr, par, split_index))
    {
        direction = GetPoleAzimuthDirection(par) - 1;
        if (direction < 0 || direction >= 3)
            return split_1D_location(ptr, par, split_index, average_grid);
    }
    else
    {
        for (int32_t i = 0; i < 3; i++)
        {
            if (ptr->blk(split_index).ijkmax[i] < ijkmax)
            {
                ijkmax = ptr->blk(split_index).ijkmax[i];
                direction = i;
            }
        }
    }
    location = (int)sqrt(average_grid / ptr->blk(split_index).ijkmax[direction]);

    if (location < split_limit)
    {
        direction++;
        return split_1D_location(ptr, par, split_index, average_grid);
    }

    for (int32_t i = 0; i < 3; i++)
    {
        if (i == direction)
            continue;
        if (abs(ptr->blk(split_index).ijkmax[i] - location) < split_limit || location > ptr->blk(split_index).ijkmax[i])
        {
            direction++;
            return split_1D_location(ptr, par, split_index, average_grid);
        }
    }

    // 一剖为4
    int32_t old_number_max = ptr->blk.Getsize1();
    direction++;
    if (direction == 1)
    {
        if (!SplitBlockIfAllowed(ptr, par, split_index, 2, location))
            return false;
        if (!SplitBlockIfAllowed(ptr, par, old_number_max, 3, location))
            return false;
        if (!SplitBlockIfAllowed(ptr, par, split_index, 3, location))
            return false;
    }
    else if (direction == 2)
    {
        if (!SplitBlockIfAllowed(ptr, par, split_index, 1, location))
            return false;
        if (!SplitBlockIfAllowed(ptr, par, old_number_max, 3, location))
            return false;
        if (!SplitBlockIfAllowed(ptr, par, split_index, 3, location))
            return false;
    }
    else if (direction == 3)
    {
        if (!SplitBlockIfAllowed(ptr, par, split_index, 1, location))
            return false;
        if (!SplitBlockIfAllowed(ptr, par, old_number_max, 2, location))
            return false;
        if (!SplitBlockIfAllowed(ptr, par, split_index, 2, location))
            return false;
    }
    return true;
}

bool split_3D_location(Preprocess_Data_Structure *ptr, Param *par, int32_t split_index, int32_t average_grid)
{
    int32_t location = 0;
    int32_t split_limit = par->GetInt("split_limit");

    if (HasPoleFaceCutByAzimuth(ptr, par, split_index))
        return split_2D_location(ptr, par, split_index, average_grid);

    location = (int)pow(average_grid, 1.0 / 3.0);

    if (location < split_limit)
    {
        return split_2D_location(ptr, par, split_index, average_grid);
    }

    for (int32_t i = 0; i < 3; i++)
    {
        if (abs(ptr->blk(split_index).ijkmax[i] - location) < split_limit || location > ptr->blk(split_index).ijkmax[i])
        {
            return split_2D_location(ptr, par, split_index, average_grid);
        }
    }

    // 一剖为8
    int32_t old_number_max = ptr->blk.Getsize1();
    if (!SplitBlockIfAllowed(ptr, par, split_index, 1, location))
        return false;
    if (!SplitBlockIfAllowed(ptr, par, split_index, 2, location))
        return false;
    if (!SplitBlockIfAllowed(ptr, par, split_index, 3, location))
        return false;

    if (!SplitBlockIfAllowed(ptr, par, old_number_max, 2, location))
        return false;
    if (!SplitBlockIfAllowed(ptr, par, old_number_max, 3, location))
        return false;

    if (!SplitBlockIfAllowed(ptr, par, old_number_max + 1, 3, location))
        return false;
    if (!SplitBlockIfAllowed(ptr, par, old_number_max + 3, 3, location))
        return false;
    return true;
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
