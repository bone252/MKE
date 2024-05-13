#pragma once
#include <fstream>
#include <iostream>
#include "Matrix.h"
#include "cgm.h"

struct node // структура узла
{
    double x;
    double y;
    double z;
};

struct edge // структура ребра
{
    node _1, _2;
};

struct material // структура материала
{
    double mu;
    double sigma;
};

struct element // структура элемента
{
    std::vector<int> edge_loc;
    int mater;
    int f_id;
};

void mult_matrix_vector(std::vector<int>& ia, std::vector<int>& ja, std::vector<double>& di, std::vector<double>& al, std::vector<double>& x, std::vector<double>& y); // умножение матрицы на вектор х, результат в у
void mult_matrix_double(std::vector<double>& di, std::vector<double>& al, double x); // умножение на месте матрицы на число
void mult_vector_double(std::vector<double>& a, double x); // умножение на месте вектора на число
void plus_vector_vector(std::vector<double>& a, std::vector<double>& b); // сложение векторов, в вектор а
void plus_matrix_matrix(std::vector<double>& diM, std::vector<double>& alM, std::vector<double>& diG, std::vector<double>& alG, std::vector<double>& di, std::vector<double>& al, double koef1); // A = koef1 * M + G


class FEM
{
private:
    std::vector<edge> all_edges;         // все ребра в порядке глобальной нумерации
    std::vector<element> all_elems;      // все элементы в прядке глобальной нумерации
    std::vector<material> all_materials; // все материалы по индексам
    std::vector<std::vector<int>> S1;    // S1[i][j] на j-ом узле заданы i-ые краевые 1 рода

    std::vector<std::vector<double>> q;  // значения q на трех слоях: q0, q1 известны, q2 ищем
    int i_t = 0;                         // текущий временной слой
    std::vector<double> time_grid;       // сетка по времени

    // матрица A
    std::vector<int> ia, ja;       // портрет
    std::vector<double> di, al, b; // симметричная

    // М
    std::vector<double> diM, alM;
    // G
    std::vector<double> diG, alG;

    int N, Kel; // кол-во ребер и элементов

    const double M = 1E+15; // большое число, используется для учета краевых первого рода

    Matrix D, G1, G2, G3, GT3; // константные матрицы, используемые для сборки локальных
    void Init_const_matrix();  // инициализация константных матриц

    double func_f(double x, double y, double z, int f_id, int num); // значение f по индексу f_id 
    double func_S1(double x, double y, double z, int s1_id);        // значение краевого S1 по индексу s1_id
    double true_func(double x, double y, double z, int num);        // истинная функция (используется для определения погрешности решения)
    int Input(); // чтение данных
    int Loc_matrix(int el_id, std::vector<Matrix>& A_loc, std::vector<double>& b_loc); // сборка локальной матрицы A = M + G
    int Loc_matrix_time(int el_id, std::vector<Matrix>& M_loc, std::vector<Matrix>& G_loc); // сборка локальных матриц M и G
    int Loc_b_time(int el_id, std::vector<double>& b_loc); // сборка локального вектора b
    int GeneratePortrait(); // генерация портрета матрицы А
    int AddLocal(int el_id, std::vector<Matrix>& A_loc);  // внесение локальной матрицы в глобальную А
    int AddLocalM(int el_id, std::vector<Matrix>& M_loc); // внесение локальной матрицы в глобальную M
    int AddLocalG(int el_id, std::vector<Matrix>& G_loc); // внесение локальной матрицы в глобальную G
    int AddLocal_b(int el_id, std::vector<double>& b_loc); // внесение локального вектора правой части в глобальный
    int SetS1(); // учет первых краевых условий

public:
    int SolveTask(std::vector<double>& res); // решение стационарной задачи
    int SolveTask_time(); // решение нестационарной задачи

    void ValueInPointXYZ(std::vector<double>& A, double x, double y, double z, int time_id); // значение численно найденной функции в произвольной точке пространства в одном из трех последних временных узлов

};
