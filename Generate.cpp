#include "Generate.h"

void Make_grid(std::string path)
{
    std::ofstream out;
    out.precision(15);
    std::vector<double> all_X, all_Y, all_Z;
    std::ifstream in(path + "grid.txt");
    double X, Y, Z, kx, ky, kz;
    int Nx, Ny, Nz;
    int count_x, count_y, count_z;
    in >> count_x >> count_y >> count_z;
    all_X.resize(count_x);
    all_Y.resize(count_y);
    all_Z.resize(count_z);

    in >> all_X[0] >> all_Y[0] >> all_Z[0];
    for (int curr_count_x = 0; curr_count_x < count_x - 1; )
    {
        in >> X >> Nx >> kx;
        double hx;
        if (kx == 1)
        {
            hx = (X - all_X[curr_count_x]) / Nx;
            for (int p = 1; p < Nx; p++)
            {
                all_X[curr_count_x + p] = all_X[curr_count_x] + hx * p;
            }
            curr_count_x += Nx;
        }
        else
        {
            hx = (X - all_X[curr_count_x]) * (kx - 1) / (pow(kx, Nx) - 1);
            for (int p = 0; p < Nx - 1; curr_count_x++, p++)
            {
                all_X[curr_count_x + 1] = all_X[curr_count_x] + hx * pow(kx, p);
            }
            curr_count_x++;
        }
        all_X[curr_count_x] = X;
    }
    for (int curr_count_y = 0; curr_count_y < count_y - 1; )
    {
        in >> Y >> Ny >> ky;
        double hy;
        if (ky == 1)
        {
            hy = (Y - all_Y[curr_count_y]) / Ny;
            for (int p = 1; p < Ny; p++)
            {
                all_Y[curr_count_y + p] = all_Y[curr_count_y] + hy * p;
            }
            curr_count_y += Ny;
        }
        else
        {
            hy = (Y - all_Y[curr_count_y]) * (ky - 1) / (pow(ky, Ny) - 1);
            for (int p = 0; p < Ny - 1; curr_count_y++, p++)
            {
                all_Y[curr_count_y + 1] = all_Y[curr_count_y] + hy * pow(ky, p);
            }
            curr_count_y++;
        }
        all_Y[curr_count_y] = Y;
    }
    for (int curr_count_z = 0; curr_count_z < count_z - 1; )
    {
        in >> Z >> Nz >> kz;
        double hy;
        if (kz == 1)
        {
            hy = (Z - all_Z[curr_count_z]) / Nz;
            for (int p = 1; p < Nz; p++)
            {
                all_Z[curr_count_z + p] = all_Z[curr_count_z] + hy * p;
            }
            curr_count_z += Nz;
        }
        else
        {
            hy = (Z - all_Z[curr_count_z]) * (kz - 1) / (pow(kz, Nz) - 1);
            for (int p = 0; p < Nz - 1; curr_count_z++, p++)
            {
                all_Z[curr_count_z + 1] = all_Z[curr_count_z] + hy * pow(kz, p);
            }
            curr_count_z++;
        }
        all_Z[curr_count_z] = Z;
    }
    in.close();

    out.open("xyz.txt");
    // ребра параллельные х
    for (int i = 0; i < count_z; i++)
    {
        for (int j = 0; j < count_x - 1; j++)
        {
            for (int k = 0; k < count_y; k++)
            {
                out << all_X[j] << "\t" << all_Y[k] << "\t" << all_Z[i] << "\t"
                    << all_X[j + 1] << "\t" << all_Y[k] << "\t" << all_Z[i] << "\n";
            }
        }
    }
    // ребра параллельные y
    for (int i = 0; i < count_z; i++)
    {
        for (int j = 0; j < count_y - 1; j++)
        {
            for (int k = 0; k < count_x; k++)
            {
                out << all_X[k] << "\t" << all_Y[j] << "\t" << all_Z[i] << "\t"
                    << all_X[k] << "\t" << all_Y[j + 1] << "\t" << all_Z[i] << "\n";
            }
        }
    }
    // ребра параллельные z
    for (int i = 0; i < count_z - 1; i++)
    {
        for (int k = 0; k < count_y; k++)
        {
            for (int j = 0; j < count_x; j++)
            {
                out << all_X[j] << "\t" << all_Y[k] << "\t" << all_Z[i] << "\t"
                    << all_X[j] << "\t" << all_Y[k] << "\t" << all_Z[i + 1] << "\n";
            }
        }
    }
    out.close();

    // input area
    // mu, sigma same in all area

    int shift_y = count_z * count_y * (count_x - 1),
        shift_z = shift_y + count_z * count_x * (count_y - 1);
    out.open("elem.txt");
    for (int i = 0; i < count_z - 1; i++)
    {
        for (int j = 0; j < count_x - 1; j++)
        {
            for (int k = 0; k < count_y - 1; k++)
            {
                out << i * count_y * (count_x - 1) + j * count_y + k << " "
                    << i * count_y * (count_x - 1) + j * count_y + k + 1 << " "
                    << (i + 1) * count_y * (count_x - 1) + j * count_y + k << " "
                    << (i + 1) * count_y * (count_x - 1) + j * count_y + k + 1 << " "

                    << shift_y + i * count_x * (count_y - 1) + j + k * count_x << " "
                    << shift_y + i * count_x * (count_y - 1) + j + 1 + k * count_x << " "
                    << shift_y + (i + 1) * count_x * (count_y - 1) + j + k * count_x << " "
                    << shift_y + (i + 1) * count_x * (count_y - 1) + j + 1 + k * count_x << " "

                    << shift_z + i * count_y * count_x + k * count_x + j << " "
                    << shift_z + i * count_y * count_x + k * count_x + 1 + j << " "
                    << shift_z + i * count_y * count_x + (k + 1) * count_x + j << " "
                    << shift_z + i * count_y * count_x + (k + 1) * count_x + 1 + j

                    << " 0 0 \n";
            }
        }
    }
    out.close();

    // bounder
    out.open("S1.txt");
    // параллельно х
    out << count_y * (count_x - 1) * 2 + (count_z - 2) * (count_x - 1) * 2 << "\n";
    // верх - низ
    for (int k = 0; k < count_y; k++)
    {
        for (int j = 0; j < count_x - 1; j++)
        {
            out << (count_z - 1) * count_y * (count_x - 1) + j * count_y + k << " " << j * count_y + k << " ";
        }
    }
    // лево - право
    for (int i = 1; i < count_z - 1; i++)
    {
        for (int j = 0; j < count_x - 1; j++)
        {
            out << i * count_y * (count_x - 1) + j * count_y << " " << i * count_y * (count_x - 1) + (j + 1) * count_y - 1 << " ";
        }
    }

    // параллельно у
    out << "\n" << 2 * count_x * (count_y - 1) + 2 * (count_z - 2) * (count_y - 1) << "\n";
    // верх - низ
    for (int k = 0; k < count_y - 1; k++)
    {
        for (int j = 0; j < count_x; j++)
        {
            out << shift_y + (count_z - 1) * count_x * (count_y - 1) + k * count_x + j << " "
                << shift_y + k * count_x + j << " ";
        }
    }
    // перед - зад
    for (int i = 1; i < count_z - 1; i++)
    {
        for (int k = 0; k < count_y - 1; k++)
        {
            out << shift_y + i * count_x * (count_y - 1) + k * count_x << " "
                << shift_y + i * count_x * (count_y - 1) + (k + 1) * count_x - 1 << " ";
        }
    }

    // параллельно z
    out << "\n" << 2 * (count_z - 1) * count_x + 2 * (count_z - 1) * (count_y - 2) << "\n";
    // лево - право
    for (int i = 0; i < count_z - 1; i++)
    {
        for (int j = 0; j < count_x; j++)
        {
            out << shift_z + i * count_y * count_x + j << " "
                << shift_z + i * count_y * count_x + count_x * (count_y - 1) + j << " ";
        }
    }
    // перед - зад
    for (int i = 0; i < count_z - 1; i++)
    {
        for (int k = 1; k < count_y - 1; k++)
        {
            out << shift_z + i * count_y * count_x + k * count_x << " "
                << shift_z + i * count_y * count_x + (k + 1) * count_x - 1 << " ";
        }
    }
    out.close();

    out.open("info.txt");
    out << (count_x - 1) * count_y * count_z + count_x * (count_y - 1) * count_z + count_x * count_y * (count_z - 1) << " 1 "
        << (count_x - 1) * (count_y - 1) * (count_z - 1) << " 3";
    out.close();

    out.open("material.txt");
    out << "1 1";
    out.close();
}

void Create_time_grid()
{
    std::ofstream out;
    out.precision(15);
    std::ifstream in;
    // time grid
    in.open("time_grid.txt");
    std::vector<double> time_grid;
    double T, kt;
    int Nt;
    int count_t;
    in >> count_t;
    time_grid.resize(count_t);

    in >> time_grid[0];
    for (int curr_count_t = 0; curr_count_t < count_t - 1; )
    {
        in >> T >> Nt >> kt;
        double ht;
        if (kt == 1)
        {
            ht = (T - time_grid[curr_count_t]) / Nt;
            for (int p = 1; p < Nt; p++)
            {
                time_grid[curr_count_t + p] = time_grid[curr_count_t] + ht * p;
            }
            curr_count_t += Nt;
        }
        else
        {
            ht = (T - time_grid[curr_count_t]) * (kt - 1) / (pow(kt, Nt) - 1);
            double pow_kt = 1;
            for (int p = 0; p < Nt - 1; curr_count_t++, p++)
            {
                time_grid[curr_count_t + 1] = time_grid[curr_count_t] + ht * pow_kt;
                pow_kt *= kt;
            }
            curr_count_t++;
        }
        time_grid[curr_count_t] = T;
    }
    in.close();
    out.open("time.txt");
    out << time_grid.size() << "\n";
    for (int i = 0; i < count_t; i++)
    {
        out << time_grid[i] << " ";
    }
    out.close();
}
