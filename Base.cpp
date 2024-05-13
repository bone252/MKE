#include "Base.h"

double FEM::func_f(double x, double y, double z, int f_id, int num) // значение f по индексу f_id 
{
    double t = time_grid[i_t];

    switch (f_id)
    {
    case 0:
        switch (num)
        {
        case 0: // F[0]
            return   3 * t * y * t
                ;
            break;
        case 1: // F[1]
            return  3 * t * z * t
                ;
            break;
        case 2: // F[2]
            return 3 * t * x * t
                ;
            break;
        default:
            break;
        }
    default:
        std::cout << "can't find f № " << f_id << "\n";
        break;
    }
    return 0;
}

double FEM::func_S1(double x, double y, double z, int s1_id) // значение краевого S1 по индексу f_id
{
    double t = time_grid[i_t];

    switch (s1_id)
    {
    case 0:
        return  y * t * t * t
            ;
    case 1:
        return   z * t * t * t
            ;
    case 2:
        return  x * t * t * t
            ;

    default:
        std::cout << "can't find s1 № " << s1_id << "\n";
        break;
    }
}

double FEM::true_func(double x, double y, double z, int num)
{
    double t = time_grid[i_t];
    switch (num)
    {
    case 0:
        return y * t * t * t
            ;
    case 1:
        return  z * t * t * t
            ;
    case 2:
        return x * t * t * t
            ;
    default:
        std::cout << "can't find true_func № " << num << "\n";
        break;
    }
}

int FEM::Input() // чтение данных
{
    int //N, // всего ребер
        //Nx, Ny, Nz,
        Nmat, // всего материалов
        //Kel, // всего элементов
        NS1; // всего 1ых краевых
    std::ifstream in;

    in.open("info.txt");
    in >> N >> Nmat >> Kel >> NS1;
    in.close();

    in.open("xyz.txt");

    all_edges.resize(N);
    for (int i = 0; i < N; i++)
    {
        in >> all_edges[i]._1.x >>
            all_edges[i]._1.y >>
            all_edges[i]._1.z >>
            all_edges[i]._2.x >>
            all_edges[i]._2.y >>
            all_edges[i]._2.z;
    }
    in.close();

    in.open("S1.txt");
    S1.resize(NS1);
    for (int i = 0; i < NS1; i++)
    {
        int size;
        in >> size;
        S1[i].resize(size);
        for (int j = 0; j < size; j++)
        {
            in >> S1[i][j];
        }
    }
    in.close();

    in.open("material.txt");
    all_materials.resize(Nmat);
    for (int i = 0; i < Nmat; i++)
    {
        in >> all_materials[i].mu >> all_materials[i].sigma;
    }
    in.close();

    in.open("elem.txt");
    all_elems.resize(Kel);
    for (int i = 0; i < Kel; i++)
    {
        all_elems[i].edge_loc.resize(12);
        in >> all_elems[i].edge_loc[0] >> all_elems[i].edge_loc[1]
            >> all_elems[i].edge_loc[2] >> all_elems[i].edge_loc[3]
            >> all_elems[i].edge_loc[4] >> all_elems[i].edge_loc[5]
            >> all_elems[i].edge_loc[6] >> all_elems[i].edge_loc[7]
            >> all_elems[i].edge_loc[8] >> all_elems[i].edge_loc[9]
            >> all_elems[i].edge_loc[10] >> all_elems[i].edge_loc[11]
            >> all_elems[i].mater >> all_elems[i].f_id;
    }
    in.close();

    q.resize(3);
    q[0].resize(N);
    q[1].resize(N);
    q[2].resize(N);
    in.open("time.txt");
    int Nt;
    in >> Nt;
    time_grid.resize(Nt);
    for (int i = 0; i < Nt; i++)
    {
        in >> time_grid[i];
    }
    // заполним 0 и 1 слои
    for (int k = 0; k < 2; k++)
    {
        i_t = k;
        for (int i = 0; i < N; i++)
        {
            double x_mid = (all_edges[i]._2.x + all_edges[i]._1.x) / 2,
                y_mid = (all_edges[i]._2.y + all_edges[i]._1.y) / 2,
                z_mid = (all_edges[i]._2.z + all_edges[i]._1.z) / 2;
            int num;
            if (all_edges[i]._2.x != all_edges[i]._1.x)
                num = 0;
            else if (all_edges[i]._2.y != all_edges[i]._1.y)
                num = 1;
            else
                num = 2;
            q[k][i] = true_func(x_mid, y_mid, z_mid, num);
        }
    }
    in.close();


    return 0;
}

void FEM::Init_const_matrix()
{
    D.matrix[0][0] = 4;
    D.matrix[0][1] = 2;
    D.matrix[0][2] = 2;
    D.matrix[0][3] = 1;
    D.matrix[1][0] = 2;
    D.matrix[1][1] = 4;
    D.matrix[1][2] = 1;
    D.matrix[1][3] = 2;
    D.matrix[2][0] = 2;
    D.matrix[2][1] = 1;
    D.matrix[2][2] = 4;
    D.matrix[2][3] = 2;
    D.matrix[3][0] = 1;
    D.matrix[3][1] = 2;
    D.matrix[3][2] = 2;
    D.matrix[3][3] = 4;

    G1.matrix[0][0] = 2;
    G1.matrix[0][1] = 1;
    G1.matrix[0][2] = -2;
    G1.matrix[0][3] = -1;
    G1.matrix[1][0] = 1;
    G1.matrix[1][1] = 2;
    G1.matrix[1][2] = -1;
    G1.matrix[1][3] = -2;
    G1.matrix[2][0] = -2;
    G1.matrix[2][1] = -1;
    G1.matrix[2][2] = 2;
    G1.matrix[2][3] = 1;
    G1.matrix[3][0] = -1;
    G1.matrix[3][1] = -2;
    G1.matrix[3][2] = 1;
    G1.matrix[3][3] = 2;

    G2.matrix[0][0] = 2;
    G2.matrix[0][1] = -2;
    G2.matrix[0][2] = 1;
    G2.matrix[0][3] = -1;
    G2.matrix[1][0] = -2;
    G2.matrix[1][1] = 2;
    G2.matrix[1][2] = -1;
    G2.matrix[1][3] = 1;
    G2.matrix[2][0] = 1;
    G2.matrix[2][1] = -1;
    G2.matrix[2][2] = 2;
    G2.matrix[2][3] = -2;
    G2.matrix[3][0] = -1;
    G2.matrix[3][1] = 1;
    G2.matrix[3][2] = -2;
    G2.matrix[3][3] = 2;

    G3.matrix[0][0] = -2;
    G3.matrix[0][1] = 2;
    G3.matrix[0][2] = -1;
    G3.matrix[0][3] = 1;
    G3.matrix[1][0] = -1;
    G3.matrix[1][1] = 1;
    G3.matrix[1][2] = -2;
    G3.matrix[1][3] = 2;
    G3.matrix[2][0] = 2;
    G3.matrix[2][1] = -2;
    G3.matrix[2][2] = 1;
    G3.matrix[2][3] = -1;
    G3.matrix[3][0] = 1;
    G3.matrix[3][1] = -1;
    G3.matrix[3][2] = 2;
    G3.matrix[3][3] = -2;

    GT3.matrix[0][0] = -2;
    GT3.matrix[1][0] = 2;
    GT3.matrix[2][0] = -1;
    GT3.matrix[3][0] = 1;
    GT3.matrix[0][1] = -1;
    GT3.matrix[1][1] = 1;
    GT3.matrix[2][1] = -2;
    GT3.matrix[3][1] = 2;
    GT3.matrix[0][2] = 2;
    GT3.matrix[1][2] = -2;
    GT3.matrix[2][2] = 1;
    GT3.matrix[3][2] = -1;
    GT3.matrix[0][3] = 1;
    GT3.matrix[1][3] = -1;
    GT3.matrix[2][3] = 2;
    GT3.matrix[3][3] = -2;
}

int FEM::Loc_matrix(int el_id, std::vector<Matrix>& A_loc, std::vector<double>& b_loc)
{
    double x1 = all_edges[all_elems[el_id].edge_loc[0]]._1.x,
        x2 = all_edges[all_elems[el_id].edge_loc[0]]._2.x,
        y1 = all_edges[all_elems[el_id].edge_loc[4]]._1.y,
        y2 = all_edges[all_elems[el_id].edge_loc[4]]._2.y,
        z1 = all_edges[all_elems[el_id].edge_loc[8]]._1.z,
        z2 = all_edges[all_elems[el_id].edge_loc[8]]._2.z,
        hx = x2 - x1,
        hy = y2 - y1,
        hz = z2 - z1,
        mu = all_materials[all_elems[el_id].mater].mu,
        gam = all_materials[all_elems[el_id].mater].sigma,
        koef_M = hx * hy * hz / 36;

    //локальная матрица масс для gam = 1
    A_loc[0] = D * koef_M;
    A_loc[4] = D * koef_M;
    A_loc[8] = D * koef_M;

    //b = C*F
    std::vector<double> F(12);
    double x_mid = (x2 + x1) / 2,
        y_mid = (y2 + y1) / 2,
        z_mid = (z2 + z1) / 2;
    F[0] = func_f(x_mid, y1, z1, all_elems[el_id].f_id, 0);
    F[1] = func_f(x_mid, y2, z1, all_elems[el_id].f_id, 0);
    F[2] = func_f(x_mid, y1, z2, all_elems[el_id].f_id, 0);
    F[3] = func_f(x_mid, y2, z2, all_elems[el_id].f_id, 0);
    F[4] = func_f(x1, y_mid, z1, all_elems[el_id].f_id, 1);
    F[5] = func_f(x2, y_mid, z1, all_elems[el_id].f_id, 1);
    F[6] = func_f(x1, y_mid, z2, all_elems[el_id].f_id, 1);
    F[7] = func_f(x2, y_mid, z2, all_elems[el_id].f_id, 1);
    F[8] = func_f(x1, y1, z_mid, all_elems[el_id].f_id, 2);
    F[9] = func_f(x2, y1, z_mid, all_elems[el_id].f_id, 2);
    F[10] = func_f(x1, y2, z_mid, all_elems[el_id].f_id, 2);
    F[11] = func_f(x2, y2, z_mid, all_elems[el_id].f_id, 2);

    for (int k = 0; k < 12; k += 4)
    {
        for (int i = k; i < k + 4; i++)
        {
            b_loc[i] = 0;
            for (int j = 0; j < 4; j++)
            {
                b_loc[i] += A_loc[k].matrix[i % 4][j] * F[k + j];
            }
        }
    }

    // собираем локальную матрицу с учетом гаммы и матрицы жесткости
    A_loc[0] = A_loc[0] * gam + (G1 * (hx * hy / (6 * hz)) + G2 * (hx * hz / (6 * hy))) * (1 / mu);
    A_loc[1] = G2 * (-hz / 6 / mu);
    A_loc[2] = G3 * (hy / 6 / mu);
    A_loc[3] = G2 * (-hz / 6 / mu);
    A_loc[4] = A_loc[4] * gam + (G1 * (hx * hy / (6 * hz)) + G2 * (hy * hz / (6 * hx))) * (1 / mu);
    A_loc[5] = G1 * (-hx / 6 / mu);
    A_loc[6] = GT3 * (hy / 6 / mu);
    A_loc[7] = G1 * (-hx / 6 / mu);
    A_loc[8] = A_loc[8] * gam + (G1 * (hx * hz / (6 * hy)) + G2 * (hy * hz / (6 * hx))) * (1 / mu);

    return 0;
}

int FEM::Loc_matrix_time(int el_id, std::vector<Matrix>& M_loc, std::vector<Matrix>& G_loc/*, std::vector<double>& b_loc*/)
{
    double x1 = all_edges[all_elems[el_id].edge_loc[0]]._1.x,
        x2 = all_edges[all_elems[el_id].edge_loc[0]]._2.x,
        y1 = all_edges[all_elems[el_id].edge_loc[4]]._1.y,
        y2 = all_edges[all_elems[el_id].edge_loc[4]]._2.y,
        z1 = all_edges[all_elems[el_id].edge_loc[8]]._1.z,
        z2 = all_edges[all_elems[el_id].edge_loc[8]]._2.z,
        hx = x2 - x1,
        hy = y2 - y1,
        hz = z2 - z1,
        mu = all_materials[all_elems[el_id].mater].mu,
        gam = all_materials[all_elems[el_id].mater].sigma, // это сигма
        koef_M = hx * hy * hz / 36;

    // локальная матрица масс
    M_loc[0] = D * koef_M * gam;
    M_loc[4] = D * koef_M * gam;
    M_loc[8] = D * koef_M * gam;

    // собираем локальную матрицу жесткости
    G_loc[0] = (G1 * (hx * hy / (6 * hz)) + G2 * (hx * hz / (6 * hy))) * (1 / mu);
    G_loc[1] = G2 * (-hz / 6 / mu);
    G_loc[2] = G3 * (hy / 6 / mu);
    G_loc[3] = G2 * (-hz / 6 / mu);
    G_loc[4] = (G1 * (hx * hy / (6 * hz)) + G2 * (hy * hz / (6 * hx))) * (1 / mu);
    G_loc[5] = G1 * (-hx / 6 / mu);
    G_loc[6] = GT3 * (hy / 6 / mu);
    G_loc[7] = G1 * (-hx / 6 / mu);
    G_loc[8] = (G1 * (hx * hz / (6 * hy)) + G2 * (hy * hz / (6 * hx))) * (1 / mu);

    return 0;
}

int FEM::Loc_b_time(int el_id, std::vector<double>& b_loc)
{
    std::vector<Matrix> M_loc(9);
    double x1 = all_edges[all_elems[el_id].edge_loc[0]]._1.x,
        x2 = all_edges[all_elems[el_id].edge_loc[0]]._2.x,
        y1 = all_edges[all_elems[el_id].edge_loc[4]]._1.y,
        y2 = all_edges[all_elems[el_id].edge_loc[4]]._2.y,
        z1 = all_edges[all_elems[el_id].edge_loc[8]]._1.z,
        z2 = all_edges[all_elems[el_id].edge_loc[8]]._2.z,
        hx = x2 - x1,
        hy = y2 - y1,
        hz = z2 - z1,
        koef_M = hx * hy * hz / 36;

    //локальная матрица масс для гамма = 1
    M_loc[0] = D * koef_M;
    M_loc[4] = D * koef_M;
    M_loc[8] = D * koef_M;

    //b = C*F
    std::vector<double> F(12);
    double x_mid = (x2 + x1) / 2,
        y_mid = (y2 + y1) / 2,
        z_mid = (z2 + z1) / 2;
    F[0] = func_f(x_mid, y1, z1, all_elems[el_id].f_id, 0);
    F[1] = func_f(x_mid, y2, z1, all_elems[el_id].f_id, 0);
    F[2] = func_f(x_mid, y1, z2, all_elems[el_id].f_id, 0);
    F[3] = func_f(x_mid, y2, z2, all_elems[el_id].f_id, 0);
    F[4] = func_f(x1, y_mid, z1, all_elems[el_id].f_id, 1);
    F[5] = func_f(x2, y_mid, z1, all_elems[el_id].f_id, 1);
    F[6] = func_f(x1, y_mid, z2, all_elems[el_id].f_id, 1);
    F[7] = func_f(x2, y_mid, z2, all_elems[el_id].f_id, 1);
    F[8] = func_f(x1, y1, z_mid, all_elems[el_id].f_id, 2);
    F[9] = func_f(x2, y1, z_mid, all_elems[el_id].f_id, 2);
    F[10] = func_f(x1, y2, z_mid, all_elems[el_id].f_id, 2);
    F[11] = func_f(x2, y2, z_mid, all_elems[el_id].f_id, 2);

    for (int k = 0; k < 12; k += 4)
    {
        for (int i = k; i < k + 4; i++)
        {
            b_loc[i] = 0;
            for (int j = 0; j < 4; j++)
            {
                b_loc[i] += M_loc[k].matrix[i % 4][j] * F[k + j];
            }
        }
    }

    return 0;
}

int FEM::GeneratePortrait() // генерация портрета
{
    di.resize(N);
    diM.resize(N);
    diG.resize(N);
    b.resize(N);
    ia.resize(N + 1);
    ja.resize(144 * Kel);
    std::vector<int> temp_list1(144 * Kel),
        temp_list2(144 * Kel);
    std::vector<int> listbeg(N);
    int listsize = 0;
    for (int i = 0; i < N; i++)
    {
        listbeg[i] = 0;
    }
    for (int ielem = 0; ielem < Kel; ielem++)
    {
        for (int i = 0; i < 12; i++)
        {
            int k = all_elems[ielem].edge_loc[i];
            for (int j = i + 1; j < 12; j++)
            {
                int ind1 = k;
                int ind2 = all_elems[ielem].edge_loc[j];
                if (ind2 < ind1)
                {
                    ind1 = ind2;
                    ind2 = k;
                }
                int iaddr = listbeg[ind2];
                if (iaddr == 0)
                {
                    listsize++;
                    listbeg[ind2] = listsize;
                    temp_list1[listsize] = ind1;
                    temp_list2[listsize] = 0;
                }
                else
                {
                    while (temp_list1[iaddr] < ind1 && temp_list2[iaddr] > 0)
                    {
                        iaddr = temp_list2[iaddr];
                    }
                    if (temp_list1[iaddr] > ind1)
                    {
                        listsize++;
                        temp_list1[listsize] = temp_list1[iaddr];
                        temp_list2[listsize] = temp_list2[iaddr];
                        temp_list1[iaddr] = ind1;
                        temp_list2[iaddr] = listsize;
                    }
                    else if (temp_list1[iaddr] < ind1)
                    {
                        listsize++;
                        temp_list2[iaddr] = listsize;
                        temp_list1[listsize] = ind1;
                        temp_list2[listsize] = 0;
                    }
                }
            }
        }
    }
    ia[0] = 0;
    for (int i = 0; i < N; i++)
    {
        ia[i + 1] = ia[i];
        int iaddr = listbeg[i];
        while (iaddr != 0)
        {
            ja[ia[i + 1]] = temp_list1[iaddr];
            ia[i + 1]++;
            iaddr = temp_list2[iaddr];
        }
    }

    ja.resize(ia[N]);
    al.resize(ia[N]);
    alM.resize(ia[N]);
    alG.resize(ia[N]);

    return 0;
}

int FEM::AddLocal(int el_id, std::vector<Matrix>& A_loc)
// внесение локальных A в глобальную СЛАУ
{
    std::vector<int> L = all_elems[el_id].edge_loc;
    int k = all_elems[el_id].edge_loc.size(); // размерность локальной матрицы
    for (int i = 0; i < k / 3; i++)
    {
        di[L[i]] += A_loc[0].matrix[i][i];
        di[L[i + 4]] += A_loc[4].matrix[i][i];
        di[L[i + 8]] += A_loc[8].matrix[i][i];
    }

    for (int i = 0; i < 12; i++)
    {
        int temp = ia[L[i]];
        for (int j = 0; j < i; j++)
        {
            int loc_pos = (i / 4) * 3 + j / 4,
                loc_i = i % 4,
                loc_j = j % 4;
            for (int k = temp; k < ia[L[i] + 1]; k++)
            {
                if (ja[k] == L[j])
                {
                    al[k] += A_loc[loc_pos].matrix[loc_i][loc_j];
                    //std::cout << loc_pos << " " << loc_i << " " << loc_j << "\n";
                    k++;
                    break;
                }
            }
        }
    }

    return 0;
}

int FEM::AddLocalM(int el_id, std::vector<Matrix>& M_loc)
// внесение локальных A в глобальную СЛАУ
{
    std::vector<int> L = all_elems[el_id].edge_loc;
    int k = all_elems[el_id].edge_loc.size(); // размерность локальной матрицы
    for (int i = 0; i < k / 3; i++)
    {
        diM[L[i]] += M_loc[0].matrix[i][i];
        diM[L[i + 4]] += M_loc[4].matrix[i][i];
        diM[L[i + 8]] += M_loc[8].matrix[i][i];
    }

    for (int i = 0; i < 12; i++)
    {
        int temp = ia[L[i]];
        for (int j = 0; j < i; j++)
        {
            int loc_pos = (i / 4) * 3 + j / 4,
                loc_i = i % 4,
                loc_j = j % 4;
            for (int k = temp; k < ia[L[i] + 1]; k++)
            {
                if (ja[k] == L[j])
                {
                    alM[k] += M_loc[loc_pos].matrix[loc_i][loc_j];
                    k++;
                    break;
                }
            }
        }
    }

    return 0;
}
int FEM::AddLocalG(int el_id, std::vector<Matrix>& G_loc)
// внесение локальных A в глобальную СЛАУ
{
    std::vector<int> L = all_elems[el_id].edge_loc;
    int k = all_elems[el_id].edge_loc.size(); // размерность локальной матрицы
    for (int i = 0; i < k / 3; i++)
    {
        diG[L[i]] += G_loc[0].matrix[i][i];
        diG[L[i + 4]] += G_loc[4].matrix[i][i];
        diG[L[i + 8]] += G_loc[8].matrix[i][i];
    }
    for (int i = 0; i < 12; i++)
    {
        int temp = ia[L[i]];
        for (int j = 0; j < i; j++)
        {
            int loc_pos = (i / 4) * 3 + j / 4,
                loc_i = i % 4,
                loc_j = j % 4;
            for (int k = temp; k < ia[L[i] + 1]; k++)
            {
                if (ja[k] == L[j])
                {
                    alG[k] += G_loc[loc_pos].matrix[loc_i][loc_j];
                    k++;
                    break;
                }
            }
        }
    }

    return 0;
}

int FEM::AddLocal_b(int el_id, std::vector<double>& b_loc)
// внесение локальных b  в глобальную СЛАУ
{
    std::vector<int> L = all_elems[el_id].edge_loc;
    int k = all_elems[el_id].edge_loc.size(); // размерность локальной матрицы
    for (int i = 0; i < k; i++)
    {
        b[L[i]] += b_loc[i];
    }
    return 0;
}

int FEM::SetS1() // учет первых краевых
{
    int NS1 = S1.size();
    for (int i = 0; i < NS1; i++)
    {
        int s1_id = i;
        for (int j = 0; j < S1[i].size(); j++)
        {
            int node_id = S1[i][j];
            double
                x_mid = (all_edges[node_id]._2.x + all_edges[node_id]._1.x) / 2,
                y_mid = (all_edges[node_id]._2.y + all_edges[node_id]._1.y) / 2,
                z_mid = (all_edges[node_id]._2.z + all_edges[node_id]._1.z) / 2;
            di[node_id] = M;
            b[node_id] = M * func_S1(x_mid, y_mid, z_mid, s1_id);
        }
    }
    return 0;
}

int FEM::SolveTask(std::vector<double>& res)
{
    Init_const_matrix();
    Input();
    GeneratePortrait();
    std::vector<Matrix> A_loc(9);
    std::vector<double> b_loc(12);

    for (int i = 0; i < all_elems.size(); i++)
    {
        Loc_matrix(i, A_loc, b_loc);
        AddLocal(i, A_loc);
        AddLocal_b(i, b_loc);
    }
    SetS1();

    cgm_solver t(ia, ja, di, al, b);
    t.llt_preconditioning(q[2]);
    res = q[2];

    std::vector<double> U_in_point(3);
    ValueInPointXYZ(U_in_point, 1.2, 1.3, 1.4, 2);
    std::cout.precision(15);
    std::cout.imbue(std::locale("Russian"));
    std::cout << "u*\n" << true_func(1.2, 1.3, 1.4, 0) << "\n"
        << true_func(1.2, 1.3, 1.4, 1) << "\n"
        << true_func(1.2, 1.3, 1.4, 2) << "\n";
    std::cout << "u\n" << U_in_point[0] << "\n"
        << U_in_point[1] << "\n"
        << U_in_point[2] << "\n";
    return 0;
}

int FEM::SolveTask_time()
{
    std::ofstream out("result.txt");
    out.precision(15);
    out.imbue(std::locale("Russian"));

    Init_const_matrix();
    Input();
    GeneratePortrait();
    std::vector<Matrix> M_loc(9), G_loc(9);
    std::vector<double> b_loc(12);
    //----------------------------------------
    std::vector<node> test_points;
    std::ifstream input("test points.txt");
    int n_test;
    input >> n_test;
    test_points.resize(n_test);
    for (int i = 0; i < n_test; i++)
    {
        input >> test_points[i].x >> test_points[i].y >> test_points[i].z;
    }
    //-----------------------------------------
    std::cout << "start make M, G\n";
    for (int i = 0; i < all_elems.size(); i++)
    {
        std::cout << "on elem " << i << "\n";
        Loc_matrix_time(i, M_loc, G_loc);
        AddLocalM(i, M_loc);
        AddLocalG(i, G_loc);
    }

    double koef1, koef2, koef3;
    double dt01 = 0, dt02 = 0, dt12 = 0;


    for (i_t = 2; i_t < time_grid.size(); i_t++)
    {
        if (dt01 != time_grid[i_t] - time_grid[i_t - 1] ||
            dt02 != time_grid[i_t] - time_grid[i_t - 2] ||
            dt12 != time_grid[i_t - 1] - time_grid[i_t - 2])
        {
            dt01 = time_grid[i_t] - time_grid[i_t - 1];
            dt02 = time_grid[i_t] - time_grid[i_t - 2];
            dt12 = time_grid[i_t - 1] - time_grid[i_t - 2];
            koef1 = (dt01 + dt02) / dt01 / dt02;
            koef2 = dt02 / dt01 / dt12;
            koef3 = -dt01 / dt02 / dt12;

            plus_matrix_matrix(diM, alM, diG, alG, di, al, koef1); // A = koef1 * M + G
        }
        // (koef1 * M + G) * q_(2) = b_(j) + koef2 * M * q_(1) + koef3 * M * q_(0) 

        // пересобрать вектор правой части
        std::fill(b.begin(), b.end(), 0);
        for (int i = 0; i < all_elems.size(); i++)
        {
            Loc_b_time(i, b_loc);
            AddLocal_b(i, b_loc);
        }
        std::vector<double> temp(N);
        // koef2 * M * q_(1)
        mult_matrix_vector(ia, ja, diM, alM, q[1], temp);
        mult_vector_double(temp, koef2);
        plus_vector_vector(b, temp);

        std::fill(temp.begin(), temp.end(), 0);
        // koef3 * M * q_(0)
        mult_matrix_vector(ia, ja, diM, alM, q[0], temp);
        mult_vector_double(temp, koef3);
        plus_vector_vector(b, temp);

        SetS1();

        cgm_solver t(ia, ja, di, al, b);
        t.llt_preconditioning(q[2]);

        if ((int)time_grid[i_t] > 1 && (int)time_grid[i_t] != (int)time_grid[i_t - 1])
        {
            out << "time = " << "\t" << time_grid[i_t] << "\n"
                << "x" << "\t" << "y" << "\t" << "z" << "\t" << "u*" << "\t" << "u" << "\t" << "|u* - u|\n";
            for (int i = 0; i < test_points.size(); i++)
            {
                std::vector<double> U_in_point(3);
                double tr_u0 = true_func(test_points[i].x, test_points[i].y, test_points[i].z, 0),
                    tr_u1 = true_func(test_points[i].x, test_points[i].y, test_points[i].z, 1),
                    tr_u2 = true_func(test_points[i].x, test_points[i].y, test_points[i].z, 2);
                ValueInPointXYZ(U_in_point, test_points[i].x, test_points[i].y, test_points[i].z, 2);

                out
                    << test_points[i].x << "\t" << test_points[i].y << "\t" << test_points[i].z << "\t"
                    << tr_u0 << "\t" << U_in_point[0] << "\t" << abs(U_in_point[0] - tr_u0) << "\n\t\t\t"
                    << tr_u1 << "\t" << U_in_point[1] << "\t" << abs(U_in_point[1] - tr_u1) << "\n\t\t\t"
                    << tr_u2 << "\t" << U_in_point[2] << "\t" << abs(U_in_point[2] - tr_u2) << "\n";
            }
        }
        q[0].swap(q[1]);
        q[1].swap(q[2]);
    }
    return 0;
}

void FEM::ValueInPointXYZ(std::vector<double>& U_in_point, double x, double y, double z, int time_id)
{
    // find elem
    double x1, x2, y1, y2, z1, z2;
    int cur_el = -1;
    for (int i = 0; i < all_elems.size() && cur_el == -1; i++)
    {

        x1 = all_edges[all_elems[i].edge_loc[0]]._1.x;
        x2 = all_edges[all_elems[i].edge_loc[0]]._2.x;
        y1 = all_edges[all_elems[i].edge_loc[0]]._1.y;
        y2 = all_edges[all_elems[i].edge_loc[1]]._1.y;
        z1 = all_edges[all_elems[i].edge_loc[0]]._1.z;
        z2 = all_edges[all_elems[i].edge_loc[3]]._1.z;
        if (x1 < x && x < x2 &&
            y1 < y && y < y2 &&
            z1 < z && z < z2)
            cur_el = i;
    }
    if (cur_el == -1)
    {
        std::cout << "can't find (" << x << " " << y << " " << z << ")\n";
        return;
    }
    double hx = x2 - x1, hy = y2 - y1, hz = z2 - z1;

    U_in_point[0] =
        q[time_id][all_elems[cur_el].edge_loc[0]] * (y2 - y) * (z2 - z) / hy / hz +
        q[time_id][all_elems[cur_el].edge_loc[1]] * (y - y1) * (z2 - z) / hy / hz +
        q[time_id][all_elems[cur_el].edge_loc[2]] * (y2 - y) * (z - z1) / hy / hz +
        q[time_id][all_elems[cur_el].edge_loc[3]] * (y - y1) * (z - z1) / hy / hz;
    U_in_point[1] =
        q[time_id][all_elems[cur_el].edge_loc[4]] * (x2 - x) * (z2 - z) / hx / hz +
        q[time_id][all_elems[cur_el].edge_loc[5]] * (x - x1) * (z2 - z) / hx / hz +
        q[time_id][all_elems[cur_el].edge_loc[6]] * (x2 - x) * (z - z1) / hx / hz +
        q[time_id][all_elems[cur_el].edge_loc[7]] * (x - x1) * (z - z1) / hx / hz;
    U_in_point[2] =
        q[time_id][all_elems[cur_el].edge_loc[8]] * (x2 - x) * (y2 - y) / hx / hy +
        q[time_id][all_elems[cur_el].edge_loc[9]] * (x - x1) * (y2 - y) / hx / hy +
        q[time_id][all_elems[cur_el].edge_loc[10]] * (x2 - x) * (y - y1) / hx / hy +
        q[time_id][all_elems[cur_el].edge_loc[11]] * (x - x1) * (y - y1) / hx / hy;

}

void mult_matrix_vector(std::vector<int>& ia, std::vector<int>& ja, std::vector<double>& di, std::vector<double>& al, std::vector<double>& x, std::vector<double>& y)
{
    for (int i = 0; i < di.size(); i++)
    {
        y[i] = di[i] * x[i];
        for (int j = ia[i]; j < ia[i + 1]; j++)
        {
            int k = ja[j];
            y[i] += al[j] * x[k];
            y[k] += al[j] * x[i];
        }
    }
}
void mult_matrix_double(std::vector<double>& di, std::vector<double>& al, double x)
{
    for (int i = 0; i < di.size(); i++)
    {
        di[i] *= x;
    }
    for (int i = 0; i < al.size(); i++)
    {
        al[i] *= x;
    }
}

void mult_vector_double(std::vector<double>& a, double x)
{
    for (int i = 0; i < a.size(); i++)
    {
        a[i] *= x;
    }
}
void plus_vector_vector(std::vector<double>& a, std::vector<double>& b)
{
    for (int i = 0; i < a.size(); i++)
    {
        a[i] += b[i];
    }
}
void plus_matrix_matrix(
    std::vector<double>& diM, std::vector<double>& alM, std::vector<double>& diG, std::vector<double>& alG,
    std::vector<double>& di, std::vector<double>& al, double koef1)
{
    for (int i = 0; i < di.size(); i++)
    {
        di[i] = diM[i];
    }
    for (int i = 0; i < al.size(); i++)
    {
        al[i] = alM[i];
    }
    mult_matrix_double(di, al, koef1);
    for (int i = 0; i < di.size(); i++)
    {
        di[i] += diG[i];
    }
    for (int i = 0; i < al.size(); i++)
    {
        al[i] += alG[i];
    }
}