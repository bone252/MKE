#pragma once
#include <fstream>
#include <iostream>
#include "Matrix.h"
#include "cgm.h"

struct node // ��������� ����
{
    double x;
    double y;
    double z;
};

struct edge // ��������� �����
{
    node _1, _2;
};

struct material // ��������� ���������
{
    double mu;
    double sigma;
};

struct element // ��������� ��������
{
    std::vector<int> edge_loc;
    int mater;
    int f_id;
};

void mult_matrix_vector(std::vector<int>& ia, std::vector<int>& ja, std::vector<double>& di, std::vector<double>& al, std::vector<double>& x, std::vector<double>& y); // ��������� ������� �� ������ �, ��������� � �
void mult_matrix_double(std::vector<double>& di, std::vector<double>& al, double x); // ��������� �� ����� ������� �� �����
void mult_vector_double(std::vector<double>& a, double x); // ��������� �� ����� ������� �� �����
void plus_vector_vector(std::vector<double>& a, std::vector<double>& b); // �������� ��������, � ������ �
void plus_matrix_matrix(std::vector<double>& diM, std::vector<double>& alM, std::vector<double>& diG, std::vector<double>& alG, std::vector<double>& di, std::vector<double>& al, double koef1); // A = koef1 * M + G


class FEM
{
private:
    std::vector<edge> all_edges;         // ��� ����� � ������� ���������� ���������
    std::vector<element> all_elems;      // ��� �������� � ������ ���������� ���������
    std::vector<material> all_materials; // ��� ��������� �� ��������
    std::vector<std::vector<int>> S1;    // S1[i][j] �� j-�� ���� ������ i-�� ������� 1 ����

    std::vector<std::vector<double>> q;  // �������� q �� ���� �����: q0, q1 ��������, q2 ����
    int i_t = 0;                         // ������� ��������� ����
    std::vector<double> time_grid;       // ����� �� �������

    // ������� A
    std::vector<int> ia, ja;       // �������
    std::vector<double> di, al, b; // ������������

    // �
    std::vector<double> diM, alM;
    // G
    std::vector<double> diG, alG;

    int N, Kel; // ���-�� ����� � ���������

    const double M = 1E+15; // ������� �����, ������������ ��� ����� ������� ������� ����

    Matrix D, G1, G2, G3, GT3; // ����������� �������, ������������ ��� ������ ���������
    void Init_const_matrix();  // ������������� ����������� ������

    double func_f(double x, double y, double z, int f_id, int num); // �������� f �� ������� f_id 
    double func_S1(double x, double y, double z, int s1_id);        // �������� �������� S1 �� ������� s1_id
    double true_func(double x, double y, double z, int num);        // �������� ������� (������������ ��� ����������� ����������� �������)
    int Input(); // ������ ������
    int Loc_matrix(int el_id, std::vector<Matrix>& A_loc, std::vector<double>& b_loc); // ������ ��������� ������� A = M + G
    int Loc_matrix_time(int el_id, std::vector<Matrix>& M_loc, std::vector<Matrix>& G_loc); // ������ ��������� ������ M � G
    int Loc_b_time(int el_id, std::vector<double>& b_loc); // ������ ���������� ������� b
    int GeneratePortrait(); // ��������� �������� ������� �
    int AddLocal(int el_id, std::vector<Matrix>& A_loc);  // �������� ��������� ������� � ���������� �
    int AddLocalM(int el_id, std::vector<Matrix>& M_loc); // �������� ��������� ������� � ���������� M
    int AddLocalG(int el_id, std::vector<Matrix>& G_loc); // �������� ��������� ������� � ���������� G
    int AddLocal_b(int el_id, std::vector<double>& b_loc); // �������� ���������� ������� ������ ����� � ����������
    int SetS1(); // ���� ������ ������� �������

public:
    int SolveTask(std::vector<double>& res); // ������� ������������ ������
    int SolveTask_time(); // ������� �������������� ������

    void ValueInPointXYZ(std::vector<double>& A, double x, double y, double z, int time_id); // �������� �������� ��������� ������� � ������������ ����� ������������ � ����� �� ���� ��������� ��������� �����

};
