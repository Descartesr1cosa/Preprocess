#include "MPCNS_Pre_InputRefinement.h"

#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <vector>

namespace
{
struct BlockInfo
{
    int old_n[3];
    int new_n[3];
    std::string name;
};

struct GridBlock
{
    int n[3];
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;
};

int RefinedIndex(int q, int ratio)
{
    if (q == 0)
        return 0;

    const int abs_q = std::abs(q);
    const int refined = (abs_q - 1) * ratio + 1;
    return q < 0 ? -refined : refined;
}

int RefinedSize(int old_size, int ratio)
{
    return ratio * (old_size - 1) + 1;
}

std::string BaseName(const std::string &path)
{
    const std::string::size_type pos = path.find_last_of("/\\");
    return pos == std::string::npos ? path : path.substr(pos + 1);
}

std::string ToString(int value)
{
    std::ostringstream stream;
    stream << value;
    return stream.str();
}

void EnsureDirectory(const std::string &path)
{
    if (mkdir(path.c_str(), 0755) != 0 && errno != EEXIST)
    {
        std::cout << "#Fatal Error: Cannot create directory:\t" << path << std::endl;
        exit(-1);
    }
}

void CopyFile(const std::string &source, const std::string &target)
{
    std::ifstream in(source.c_str(), std::ios::binary);
    if (!in)
    {
        std::cout << "#Fatal Error: Cannot open source file:\t" << source << std::endl;
        exit(-1);
    }

    std::ofstream out(target.c_str(), std::ios::binary);
    if (!out)
    {
        std::cout << "#Fatal Error: Cannot open target file:\t" << target << std::endl;
        exit(-1);
    }

    out << in.rdbuf();
}

int GridIndex(int i, int j, int k, const int n[3])
{
    return (k * n[1] + j) * n[0] + i;
}

double InterpolateValue(const std::vector<double> &v, const int old_n[3], int i, int j, int k, int ratio)
{
    const double pos[3] = {
        old_n[0] == 1 ? 0.0 : static_cast<double>(i) / static_cast<double>(ratio),
        old_n[1] == 1 ? 0.0 : static_cast<double>(j) / static_cast<double>(ratio),
        old_n[2] == 1 ? 0.0 : static_cast<double>(k) / static_cast<double>(ratio)};

    int lo[3];
    int hi[3];
    double t[3];
    for (int d = 0; d < 3; ++d)
    {
        lo[d] = static_cast<int>(std::floor(pos[d]));
        if (lo[d] >= old_n[d] - 1)
        {
            lo[d] = old_n[d] - 1;
            hi[d] = lo[d];
            t[d] = 0.0;
        }
        else
        {
            hi[d] = lo[d] + 1;
            t[d] = pos[d] - static_cast<double>(lo[d]);
        }
    }

    double value = 0.0;
    for (int dk = 0; dk <= 1; ++dk)
        for (int dj = 0; dj <= 1; ++dj)
            for (int di = 0; di <= 1; ++di)
            {
                const double wi = di == 0 ? 1.0 - t[0] : t[0];
                const double wj = dj == 0 ? 1.0 - t[1] : t[1];
                const double wk = dk == 0 ? 1.0 - t[2] : t[2];
                const int oi = di == 0 ? lo[0] : hi[0];
                const int oj = dj == 0 ? lo[1] : hi[1];
                const int ok = dk == 0 ? lo[2] : hi[2];
                value += wi * wj * wk * v[GridIndex(oi, oj, ok, old_n)];
            }

    return value;
}

std::vector<GridBlock> ReadGridFile(const std::string &filename, int read_type)
{
    std::vector<GridBlock> blocks;
    if (read_type == 1 || read_type == 2)
    {
        std::ifstream file(filename.c_str(), std::ios::binary);
        if (!file)
        {
            std::cout << "#Fatal Error: Cannot open grid file:\t" << filename << std::endl;
            exit(-1);
        }

        int temp = 0;
        if (read_type == 1)
            file.read(reinterpret_cast<char *>(&temp), sizeof(temp));
        file.read(reinterpret_cast<char *>(&temp), sizeof(temp));
        const int block_count = temp;
        if (read_type == 1)
        {
            file.read(reinterpret_cast<char *>(&temp), sizeof(temp));
            file.read(reinterpret_cast<char *>(&temp), sizeof(temp));
        }

        blocks.resize(block_count);
        for (int iblock = 0; iblock < block_count; ++iblock)
            for (int d = 0; d < 3; ++d)
                file.read(reinterpret_cast<char *>(&blocks[iblock].n[d]), sizeof(int));

        if (read_type == 1)
            file.read(reinterpret_cast<char *>(&temp), sizeof(temp));

        for (int iblock = 0; iblock < block_count; ++iblock)
        {
            GridBlock &block = blocks[iblock];
            const int count = block.n[0] * block.n[1] * block.n[2];
            block.x.resize(count);
            block.y.resize(count);
            block.z.resize(count);

            if (read_type == 1)
                file.read(reinterpret_cast<char *>(&temp), sizeof(temp));
            file.read(reinterpret_cast<char *>(&block.x[0]), sizeof(double) * count);
            file.read(reinterpret_cast<char *>(&block.y[0]), sizeof(double) * count);
            file.read(reinterpret_cast<char *>(&block.z[0]), sizeof(double) * count);
            if (read_type == 1)
                file.read(reinterpret_cast<char *>(&temp), sizeof(temp));
        }
    }
    else
    {
        std::ifstream file(filename.c_str());
        if (!file)
        {
            std::cout << "#Fatal Error: Cannot open grid file:\t" << filename << std::endl;
            exit(-1);
        }

        int block_count = 0;
        file >> block_count;
        blocks.resize(block_count);
        for (int iblock = 0; iblock < block_count; ++iblock)
            file >> blocks[iblock].n[0] >> blocks[iblock].n[1] >> blocks[iblock].n[2];

        for (int iblock = 0; iblock < block_count; ++iblock)
        {
            GridBlock &block = blocks[iblock];
            const int count = block.n[0] * block.n[1] * block.n[2];
            block.x.resize(count);
            block.y.resize(count);
            block.z.resize(count);
            for (int i = 0; i < count; ++i)
                file >> block.x[i];
            for (int i = 0; i < count; ++i)
                file >> block.y[i];
            for (int i = 0; i < count; ++i)
                file >> block.z[i];
        }
    }

    return blocks;
}

void WriteBinaryGridFile(const std::string &filename, const std::vector<GridBlock> &blocks)
{
    std::ofstream file(filename.c_str(), std::ios::binary);
    if (!file)
    {
        std::cout << "#Fatal Error: Cannot write grid file:\t" << filename << std::endl;
        exit(-1);
    }

    const int block_count = static_cast<int>(blocks.size());
    file.write(reinterpret_cast<const char *>(&block_count), sizeof(block_count));
    for (int iblock = 0; iblock < block_count; ++iblock)
        for (int d = 0; d < 3; ++d)
            file.write(reinterpret_cast<const char *>(&blocks[iblock].n[d]), sizeof(int));

    for (int iblock = 0; iblock < block_count; ++iblock)
    {
        const GridBlock &block = blocks[iblock];
        const int count = block.n[0] * block.n[1] * block.n[2];
        file.write(reinterpret_cast<const char *>(&block.x[0]), sizeof(double) * count);
        file.write(reinterpret_cast<const char *>(&block.y[0]), sizeof(double) * count);
        file.write(reinterpret_cast<const char *>(&block.z[0]), sizeof(double) * count);
    }
}

std::vector<GridBlock> RefineGridBlocks(const std::vector<GridBlock> &old_blocks,
                                        const std::vector<BlockInfo> &block_info,
                                        int ratio)
{
    if (old_blocks.size() != block_info.size())
    {
        std::cout << "#Fatal Error: inp block count does not match grid block count." << std::endl;
        exit(-1);
    }

    std::vector<GridBlock> refined(old_blocks.size());
    for (std::size_t iblock = 0; iblock < old_blocks.size(); ++iblock)
    {
        const GridBlock &old_block = old_blocks[iblock];
        GridBlock &new_block = refined[iblock];
        for (int d = 0; d < 3; ++d)
        {
            if (old_block.n[d] != block_info[iblock].old_n[d])
            {
                std::cout << "#Fatal Error: inp/grid size mismatch at block " << iblock + 1 << std::endl;
                exit(-1);
            }
            new_block.n[d] = block_info[iblock].new_n[d];
        }

        const int count = new_block.n[0] * new_block.n[1] * new_block.n[2];
        new_block.x.resize(count);
        new_block.y.resize(count);
        new_block.z.resize(count);

        for (int k = 0; k < new_block.n[2]; ++k)
            for (int j = 0; j < new_block.n[1]; ++j)
                for (int i = 0; i < new_block.n[0]; ++i)
                {
                    const int index = GridIndex(i, j, k, new_block.n);
                    new_block.x[index] = InterpolateValue(old_block.x, old_block.n, i, j, k, ratio);
                    new_block.y[index] = InterpolateValue(old_block.y, old_block.n, i, j, k, ratio);
                    new_block.z[index] = InterpolateValue(old_block.z, old_block.n, i, j, k, ratio);
                }
    }

    return refined;
}

std::vector<BlockInfo> RewriteInpFile(const std::string &source, const std::string &target, int dimension, int ratio)
{
    std::ifstream in(source.c_str());
    if (!in)
    {
        std::cout << "#Fatal Error: Cannot open inp file:\t" << source << std::endl;
        exit(-1);
    }

    std::ofstream out(target.c_str());
    if (!out)
    {
        std::cout << "#Fatal Error: Cannot write inp file:\t" << target << std::endl;
        exit(-1);
    }

    std::string header;
    std::getline(in, header);
    out << "# refined from " << BaseName(source) << " by MPCNS preprocess\n";

    int block_count = 0;
    in >> block_count;
    out << block_count << "\n";

    std::vector<BlockInfo> blocks(block_count);
    for (int iblock = 0; iblock < block_count; ++iblock)
    {
        BlockInfo &block = blocks[iblock];
        if (dimension == 2)
        {
            in >> block.old_n[0] >> block.old_n[1];
            block.old_n[2] = 1;
        }
        else
        {
            in >> block.old_n[0] >> block.old_n[1] >> block.old_n[2];
        }
        in >> block.name;

        for (int d = 0; d < 3; ++d)
            block.new_n[d] = RefinedSize(block.old_n[d], ratio);

        if (dimension == 2)
            out << block.new_n[0] << "\t" << block.new_n[1] << "\t" << block.name << "\n";
        else
            out << block.new_n[0] << "\t" << block.new_n[1] << "\t" << block.new_n[2] << "\t" << block.name << "\n";

        int face_count = 0;
        in >> face_count;
        out << face_count << "\n";

        for (int iface = 0; iface < face_count; ++iface)
        {
            int sub[3] = {1, 1, 1};
            int sup[3] = {1, 1, 1};
            int first_index = 0;
            for (int d = 0; d < dimension; ++d)
                in >> sub[d] >> sup[d];
            in >> first_index;

            for (int d = 0; d < dimension; ++d)
                out << RefinedIndex(sub[d], ratio) << "\t" << RefinedIndex(sup[d], ratio) << "\t";
            out << first_index;

            if (first_index <= 0)
            {
                int tar_sub[3] = {1, 1, 1};
                int tar_sup[3] = {1, 1, 1};
                int target_block = 0;
                for (int d = 0; d < dimension; ++d)
                    in >> tar_sub[d] >> tar_sup[d];
                in >> target_block;

                for (int d = 0; d < dimension; ++d)
                    out << "\t" << RefinedIndex(tar_sub[d], ratio) << "\t" << RefinedIndex(tar_sup[d], ratio);
                out << "\t" << target_block;
            }
            out << "\n";
        }
    }

    return blocks;
}
} // namespace

void InputRefinement::RefineIfEnabled(Param &param)
{
    const int refine_level = param.HasInt("refine_level") ? param.GetInt("refine_level") : 0;
    if (refine_level <= 0)
        return;

    if (refine_level >= static_cast<int>(sizeof(int) * 8 - 2))
    {
        std::cout << "#Fatal Error: refine_level is too large:\t" << refine_level << std::endl;
        exit(-1);
    }

    const int ratio = 1 << refine_level;
    const int dimension = param.GetInt("dimension");
    const std::string out_dir = "./INPUT/refine_level_" + ToString(refine_level);
    const std::string old_inp = param.GetStr("bfilename");
    const std::string old_grid = param.GetStr("gfilename");
    const std::string old_fvbnd = param.GetStr("ffilename");
    const std::string new_inp = out_dir + "/" + BaseName(old_inp);
    const std::string new_grid = out_dir + "/" + BaseName(old_grid);
    const std::string new_fvbnd = out_dir + "/" + BaseName(old_fvbnd);

    if (param.GetInt("myid") != 0)
    {
        param.AddParam("bfilename", new_inp);
        param.AddParam("gfilename", new_grid);
        param.AddParam("ffilename", new_fvbnd);
        param.AddParam("grd_readtype", 2);
        param.AddParam("if_split_group_info", false);
        return;
    }

    EnsureDirectory("./INPUT");
    EnsureDirectory(out_dir);

    std::cout << "---->(1.5) Refining input grid by refine_level=" << refine_level
              << " (R=" << ratio << ")...\n";

    const std::vector<BlockInfo> blocks = RewriteInpFile(old_inp, new_inp, dimension, ratio);
    CopyFile(old_fvbnd, new_fvbnd);

    const std::vector<GridBlock> old_grid_blocks = ReadGridFile(old_grid, param.GetInt("grd_readtype"));
    const std::vector<GridBlock> new_grid_blocks = RefineGridBlocks(old_grid_blocks, blocks, ratio);
    WriteBinaryGridFile(new_grid, new_grid_blocks);

    param.AddParam("bfilename", new_inp);
    param.AddParam("gfilename", new_grid);
    param.AddParam("ffilename", new_fvbnd);
    param.AddParam("grd_readtype", 2);
    param.AddParam("if_split_group_info", false);

    std::cout << "\t-->Refined INPUT files have been generated in " << out_dir << "\n";
}
