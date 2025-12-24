#pragma once
#define _CRT_SECURE_NO_WARNINGS
/*
 *
 *  ���ڲ�ʹ�õ���ƽ���Ĳ��֣����Է������˴������һ�����Ľ��ƫ��
 *
 */

#include "Checked.h"
#include "roots.h"
#include "utility.h"
#include "GMRES.h"

#include <iostream>

void _get_plane_n(const enum Axis axi, const double *p1, const double *p2, const double *p3, double *n)
{
    double p12[3], p32[3];

    _get_normalized(p1, p2, p12);
    _get_normalized(p3, p2, p32);
    _multcross(p12, p32, n);

    switch (axi)
    {
    case Axis::XPOSITIVE:
        if (n[0] < 0)
            _get_reverse(n);
        break;
    case Axis::XNEGATIVE:
        if (n[0] > 0)
            _get_reverse(n);
        break;
    case Axis::YPOSITIVE:
        if (n[1] < 0)
            _get_reverse(n);
        break;
    case Axis::YNEGATIVE:
        if (n[1] > 0)
            _get_reverse(n);
        break;
    case Axis::ZPOSITIVE:
        if (n[2] < 0)
            _get_reverse(n);
        break;
    case Axis::ZNEGATIVE:
        if (n[2] > 0)
            _get_reverse(n);
        break;
    }
    _get_normalized(n);
}

void _get_plane_n(const double *p1, const double *p2, const double *p3, double *n)
{
    double p12[3], p32[3];

    _get_normalized(p1, p2, p12);
    _get_normalized(p3, p2, p32);
    _multcross(p12, p32, n);
    _get_normalized(n);
}

void _get_matK(Ktensor K, Eigen::Matrix<double, 3, 3> &mat)
{
    mat(0, 0) = K.x;
    mat(1, 1) = K.y;
    mat(2, 2) = K.z;
    mat(0, 1) = K.xy;
    mat(1, 0) = K.xy;
    mat(0, 2) = K.xz;
    mat(2, 0) = K.xz;
    mat(1, 2) = K.yz;
    mat(2, 1) = K.yz;
}

void _get_w(double *w, const Ktensor *K, double pc[][3], double p[][3],
            int *idx, int *idx_8, int ii, double *n, double S)
{
    Eigen::RowVector3d v1, v2, v3;
    Eigen::RowVector3d n_vec;
    Eigen::Matrix<double, 3, 3> matX;
    Eigen::Matrix<double, 3, 3> matK;
    Eigen::Matrix<double, 3, 3> X_ivs;
    Eigen::Vector3d x0, x1, x2;

    n_vec << n[0], n[1], n[2];
    for (int l = 0; l < 3; l++)
    {
        v1[l] = p[idx[0]][l] - pc[ii][l];
        v2[l] = p[idx[1]][l] - pc[ii][l];
        v3[l] = p[idx[2]][l] - pc[ii][l];
    }
    //_get_matK(K[idx_8[ii]], matK);
    matK(0, 0) = K[idx_8[ii]].x;
    matK(1, 1) = K[idx_8[ii]].y;
    matK(2, 2) = K[idx_8[ii]].z;
    matK(0, 1) = K[idx_8[ii]].xy;
    matK(1, 0) = K[idx_8[ii]].xy;
    matK(0, 2) = K[idx_8[ii]].xz;
    matK(2, 0) = K[idx_8[ii]].xz;
    matK(1, 2) = K[idx_8[ii]].yz;
    matK(2, 1) = K[idx_8[ii]].yz;

    matX.row(0) = v1;
    matX.row(1) = v2;
    matX.row(2) = v3;
    X_ivs = matX.inverse();
    x0 = X_ivs.col(0);
    x1 = X_ivs.col(1);
    x2 = X_ivs.col(2);
    w[0] = -S * n_vec * matK * x0; // 0'-6
    w[1] = -S * n_vec * matK * x1; // 4'-6
    w[2] = -S * n_vec * matK * x2; // 8'-6
}

void _get_vec(int i, int j, int k, const Ktensor *perm,
              const double *verts, int idx_row, double para_cell[8], int idx_cell[8], int divn)
{
    int nx = nxx * divn;
    int ny = nyy * divn;
    int nz = nzz * divn;
    std::mdspan<const double, std::extents<size_t, -1, -1, -1, 8, 3>> _verts(verts, nz, ny, nx);

    double c[8][3];      // 8���������ĵ�����
    double p[12][3];     // 12������������
    double n[12][3];     // 12���������Ӧ���񽻽���ķ�����
    double pp[6][3];     // 12�����������������6������
    double S[12];        // 12�����������
    double w[12][3];     // 12�������湹�ɷ��̵���߸���ϵ��
    double w_rvs[12][3]; // 12�������湹�ɷ��̵��ұ߸���ϵ��
    int idx_c[8][3];     // 8�������Ӧ������������������p[i][3]��������i
    int idx_w[12];       // 12�������湹�ɷ��̵���ߵ�c[i][3]������i
    int idx_w_rvs[12];   // 12�������湹�ɷ��̵��ұߵ�c[i][3]������i

    int i0 = k * nx * ny + j * nx + i;
    int i1 = i0 + nx;
    _get_pcood(nx, ny, nz, i, j, k, verts, c[0]);
    _get_pcood(nx, ny, nz, i + 1, j, k, verts, c[1]);
    _get_pcood(nx, ny, nz, i + 1, j + 1, k, verts, c[2]);
    _get_pcood(nx, ny, nz, i, j + 1, k, verts, c[3]);
    _get_pcood(nx, ny, nz, i, j, k + 1, verts, c[4]);
    _get_pcood(nx, ny, nz, i + 1, j, k + 1, verts, c[5]);
    _get_pcood(nx, ny, nz, i + 1, j + 1, k + 1, verts, c[6]);
    _get_pcood(nx, ny, nz, i, j + 1, k + 1, verts, c[7]);

    for (int l = 0; l < 3; l++)
    {
        p[0][l] = (_verts(k + 1, j + 1, i + 1, 0, l) + _verts(k + 1, j + 1, i + 1, 2, l) + _verts(k + 1, j + 1, i + 1, 4, l) + _verts(k + 1, j + 1, i + 1, 6, l)) / 4;
        p[1][l] = (_verts(k + 1, j, i + 1, 0, l) + _verts(k + 1, j, i + 1, 2, l) + _verts(k + 1, j, i + 1, 4, l) + _verts(k + 1, j, i + 1, 6, l)) / 4;
        p[2][l] = (_verts(k, j, i + 1, 0, l) + _verts(k, j, i + 1, 2, l) + _verts(k, j, i + 1, 4, l) + _verts(k, j, i + 1, 6, l)) / 4;
        p[3][l] = (_verts(k, j + 1, i + 1, 0, l) + _verts(k, j + 1, i + 1, 2, l) + _verts(k, j + 1, i + 1, 4, l) + _verts(k, j + 1, i + 1, 6, l)) / 4;

        p[4][l] = (_verts(k + 1, j, i + 1, 2, l) + _verts(k + 1, j, i + 1, 3, l) + _verts(k + 1, j, i + 1, 6, l) + _verts(k + 1, j, i + 1, 7, l)) / 4;
        p[5][l] = (_verts(k, j, i + 1, 2, l) + _verts(k, j, i + 1, 3, l) + _verts(k, j, i + 1, 6, l) + _verts(k, j, i + 1, 7, l)) / 4;
        p[6][l] = (_verts(k, j, i, 2, l) + _verts(k, j, i, 3, l) + _verts(k, j, i, 6, l) + _verts(k, j, i, 7, l)) / 4;
        p[7][l] = (_verts(k + 1, j, i, 2, l) + _verts(k + 1, j, i, 3, l) + _verts(k + 1, j, i, 6, l) + _verts(k + 1, j, i, 7, l)) / 4;

        p[8][l] = (_verts(k, j + 1, i + 1, 4, l) + _verts(k, j + 1, i + 1, 5, l) + _verts(k, j + 1, i + 1, 6, l) + _verts(k, j + 1, i + 1, 7, l)) / 4;
        p[9][l] = (_verts(k, j + 1, i, 4, l) + _verts(k, j + 1, i, 5, l) + _verts(k, j + 1, i, 6, l) + _verts(k, j + 1, i, 7, l)) / 4;
        p[10][l] = (_verts(k, j, i, 4, l) + _verts(k, j, i, 5, l) + _verts(k, j, i, 6, l) + _verts(k, j, i, 7, l)) / 4;
        p[11][l] = (_verts(k, j, i + 1, 4, l) + _verts(k, j, i + 1, 5, l) + _verts(k, j, i + 1, 6, l) + _verts(k, j, i + 1, 7, l)) / 4;

        pp[0][l] = (_verts(k, j, i, 5, l) + _verts(k, j, i, 7, l)) / 2;
        pp[1][l] = (_verts(k, j, i + 1, 6, l) + _verts(k, j, i + 1, 7, l)) / 2;
        pp[2][l] = (_verts(k, j + 1, i, 5, l) + _verts(k, j + 1, i, 7, l)) / 2;
        pp[3][l] = (_verts(k, j, i, 6, l) + _verts(k, j, i, 7, l)) / 2;
        pp[4][l] = (_verts(k, j, i, 3, l) + _verts(k, j, i, 7, l)) / 2;
        pp[5][l] = (_verts(k + 1, j, i, 3, l) + _verts(k + 1, j, i, 7, l)) / 2;
    }
    S[0] = 2 * _get_area(pp[2], pp[5], p[0]);
    S[1] = 2 * _get_area(pp[0], pp[5], p[1]);
    S[2] = 2 * _get_area(pp[0], pp[4], p[2]);
    S[3] = 2 * _get_area(pp[2], pp[4], p[3]);
    S[4] = 2 * _get_area(pp[1], pp[5], p[4]);
    S[5] = 2 * _get_area(pp[1], pp[4], p[5]);
    S[6] = 2 * _get_area(pp[3], pp[4], p[6]);
    S[7] = 2 * _get_area(pp[3], pp[5], p[7]);
    S[8] = 2 * _get_area(pp[1], pp[2], p[8]);
    S[9] = 2 * _get_area(pp[2], pp[3], p[9]);
    S[10] = 2 * _get_area(pp[0], pp[3], p[10]);
    S[11] = 2 * _get_area(pp[0], pp[1], p[11]);

    // ���������Ϊ��
    _get_plane_n(Axis::XPOSITIVE, pp[2], pp[5], p[0], n[0]);
    _get_plane_n(Axis::XPOSITIVE, pp[0], pp[5], p[1], n[1]);
    _get_plane_n(Axis::XPOSITIVE, pp[0], pp[4], p[2], n[2]);
    _get_plane_n(Axis::XPOSITIVE, pp[2], pp[4], p[3], n[3]);
    _get_plane_n(Axis::YPOSITIVE, pp[1], pp[5], p[4], n[4]);
    _get_plane_n(Axis::YPOSITIVE, pp[1], pp[4], p[5], n[5]);
    _get_plane_n(Axis::YPOSITIVE, pp[3], pp[4], p[6], n[6]);
    _get_plane_n(Axis::YPOSITIVE, pp[3], pp[5], p[7], n[7]);
    _get_plane_n(Axis::ZPOSITIVE, pp[1], pp[2], p[8], n[8]);
    _get_plane_n(Axis::ZPOSITIVE, pp[2], pp[3], p[9], n[9]);
    _get_plane_n(Axis::ZPOSITIVE, pp[0], pp[3], p[10], n[10]);
    _get_plane_n(Axis::ZPOSITIVE, pp[0], pp[1], p[11], n[11]);

    idx_cell[0] = k * nx * ny + j * nx + i;
    idx_c[0][0] = 2;
    idx_c[0][1] = 6;
    idx_c[0][2] = 10;
    idx_cell[1] = k * nx * ny + j * nx + i + 1;
    idx_c[1][0] = 2;
    idx_c[1][1] = 5;
    idx_c[1][2] = 11;
    idx_cell[2] = k * nx * ny + (j + 1) * nx + i + 1;
    idx_c[2][0] = 3;
    idx_c[2][1] = 5;
    idx_c[2][2] = 8;
    idx_cell[3] = k * nx * ny + (j + 1) * nx + i;
    idx_c[3][0] = 3;
    idx_c[3][1] = 6;
    idx_c[3][2] = 9;

    idx_cell[4] = (k + 1) * nx * ny + j * nx + i;
    idx_c[4][0] = 1;
    idx_c[4][1] = 7;
    idx_c[4][2] = 10;
    idx_cell[5] = (k + 1) * nx * ny + j * nx + i + 1;
    idx_c[5][0] = 1;
    idx_c[5][1] = 4;
    idx_c[5][2] = 11;
    idx_cell[6] = (k + 1) * nx * ny + (j + 1) * nx + i + 1;
    idx_c[6][0] = 0;
    idx_c[6][1] = 4;
    idx_c[6][2] = 8;
    idx_cell[7] = (k + 1) * nx * ny + (j + 1) * nx + i;
    idx_c[7][0] = 0;
    idx_c[7][1] = 7;
    idx_c[7][2] = 9; // 8������ļ���

    idx_w[0] = 6;
    idx_w_rvs[0] = 7;
    idx_w[1] = 5;
    idx_w_rvs[1] = 4;
    idx_w[2] = 1;
    idx_w_rvs[2] = 0;
    idx_w[3] = 2;
    idx_w_rvs[3] = 3;

    idx_w[4] = 6;
    idx_w_rvs[4] = 5;
    idx_w[5] = 2;
    idx_w_rvs[5] = 1;
    idx_w[6] = 3;
    idx_w_rvs[6] = 0;
    idx_w[7] = 7;
    idx_w_rvs[7] = 4;

    idx_w[8] = 6;
    idx_w_rvs[8] = 2;
    idx_w[9] = 7;
    idx_w_rvs[9] = 3;
    idx_w[10] = 4;
    idx_w_rvs[10] = 0;
    idx_w[11] = 5;
    idx_w_rvs[11] = 1; // 12����������������

    for (int l = 0; l < 12; l++)
    {
        _get_w(w[l], perm, c, p, idx_c[idx_w[l]], idx_cell, idx_w[l], n[l], S[l]);
        _get_w(w_rvs[l], perm, c, p, idx_c[idx_w_rvs[l]], idx_cell, idx_w_rvs[l], n[l], S[l]);
    }

    Eigen::Matrix<double, 12, 12, Eigen::RowMajor> matA;
    Eigen::Matrix<double, 12, 12, Eigen::RowMajor> matC;
    Eigen::Matrix<double, 12, 8, Eigen::RowMajor> matB;
    ;
    Eigen::Matrix<double, 12, 8, Eigen::RowMajor> matD;
    Eigen::Matrix<double, 12, 8, Eigen::RowMajor> matT;

    for (int l = 0; l < 12; l++)
    { // 12���������
        for (int m = 0; m < 12; m++)
        { // 12����
            matA(l, m) = 0.0;
            for (int n = 0; n < 3; n++)
            {
                if (m == idx_c[idx_w[l]][n])
                {
                    matA(l, m) += w[l][n];
                }
                if (m == idx_c[idx_w_rvs[l]][n])
                {
                    matA(l, m) -= w_rvs[l][n];
                }
            }
        }
    }

    for (int l = 0; l < 12; l++)
    {
        for (int m = 0; m < 8; m++)
        {
            matB(l, m) = 0.0;
            if (m == idx_w[l])
            {
                matB(l, m) += w[l][0] + w[l][1] + w[l][2];
            }
            if (m == idx_w_rvs[l])
            {
                matB(l, m) -= w_rvs[l][0] + w_rvs[l][1] + w_rvs[l][2];
            }
        }
    }

    for (int l = 0; l < 12; l++)
    {
        for (int m = 0; m < 12; m++)
        {
            matC(l, m) = 0.0;
            for (int n = 0; n < 3; n++)
            {
                if (m == idx_c[idx_w[l]][n])
                {
                    matC(l, m) += w[l][n];
                }
            }
        }
    }

    for (int l = 0; l < 12; l++)
    {
        for (int m = 0; m < 8; m++)
        {
            matD(l, m) = 0.0;
            if (m == idx_w[l])
            {
                matD(l, m) += w[l][0] + w[l][1] + w[l][2];
            }
        }
    }

    matT = matC * matA.inverse() * matB - matD;
    for (int l = 0; l < 8; ++l)
    {
        para_cell[l] = matT(idx_row, l);
    }
}

void _get_tolQx(const Ktensor *perm, const double *p, int i, int j, int k,
                const double *verts, double &Qx, int divn)
{
    int nx = nxx * divn;
    int ny = nyy * divn;
    int nz = nzz * divn;
    std::mdspan<const double, std::extents<size_t, -1, -1, -1, 8, 3>> _verts(verts, nz, ny, nx);

    const double *pp[4];
    double idxbeta, idx, idxempty;
    double pc[3], pc1[3], fc[3];
    double S = 0.0, D = 0.0;
    double T[8];
    int idx_8[8];

    int i0 = k * nx * ny + j * nx + i;
    int i1 = i0 + 1;
    _get_pcood(nx, ny, nz, i, j, k, verts, pc);
    _get_pcood(nx, ny, nz, i + 1, j, k, verts, pc1);
    pp[0] = &_verts(k, j, i, 1, 0);
    pp[1] = &_verts(k, j, i, 3, 0);
    pp[2] = &_verts(k, j, i, 5, 0);
    pp[3] = &_verts(k, j, i, 7, 0);
    for (int i = 0; i < 3; i++)
    {
        fc[i] = (pp[0][i] + pp[1][i] + pp[2][i] + pp[3][i]) / 4;
    }
    D = _get_distance(pc, pp[0], pp[1], pp[2]) + _get_distance(pc1, pp[0], pp[1], pp[2]);
    S = _get_area(pp[0], pp[1], pp[2], pp[3]) / 4;

    Qx = 0.0;
    // 1
    if (k > 0 && j > 0 && i < nx - 1)
    {
        _get_vec(i, j - 1, k - 1, perm, verts, 0, T, idx_8, divn);
        for (int l = 0; l < 8; l++)
        {
            Qx -= T[l] * p[idx_8[l]];
        }
    }
    else
    {
        Qx += (p[i1] - p[i0]) * 2.0 * perm[i0].x * perm[i1].x / (perm[i0].x + perm[i1].x) * S / D;
    }
    // 3
    if (k > 0 && j < ny - 1 && i < nx - 1)
    {
        _get_vec(i, j, k - 1, perm, verts, 1, T, idx_8, divn);
        for (int l = 0; l < 8; l++)
        {
            Qx -= T[l] * p[idx_8[l]];
        }
    }
    else
    {
        Qx += (p[i1] - p[i0]) * 2.0 * perm[i0].x * perm[i1].x / (perm[i0].x + perm[i1].x) * S / D;
    }
    // 5
    if (k < nz - 1 && j > 0 && i < nx - 1)
    {
        _get_vec(i, j - 1, k, perm, verts, 3, T, idx_8, divn);
        for (int l = 0; l < 8; l++)
        {
            Qx -= T[l] * p[idx_8[l]];
        }
    }
    else
    {
        Qx += (p[i1] - p[i0]) * 2.0 * perm[i0].x * perm[i1].x / (perm[i0].x + perm[i1].x) * S / D;
    }
    // 7
    if (k < nz - 1 && j < ny - 1 && i < nx - 1)
    {
        _get_vec(i, j, k, perm, verts, 2, T, idx_8, divn);
        for (int l = 0; l < 8; l++)
        {
            Qx -= T[l] * p[idx_8[l]];
        }
    }
    else
    {
        Qx += (p[i1] - p[i0]) * 2.0 * perm[i0].x * perm[i1].x / (perm[i0].x + perm[i1].x) * S / D;
    }
}

void _get_tolQy(const Ktensor *perm, const double *p, int i, int j, int k,
                const double *verts, double &Qy, int divn)
{
    int nx = nxx * divn;
    int ny = nyy * divn;
    int nz = nzz * divn;
    std::mdspan<const double, std::extents<size_t, -1, -1, -1, 8, 3>> _verts(verts, nz, ny, nx);

    const double *pp[4];
    double idxbeta, idx, idxempty;
    double pc[3], pc1[3], fc[3];
    double S = 0.0, D = 0.0;
    double T[8];
    int idx_8[8];

    int i0 = k * nx * ny + j * nx + i;
    int i1 = i0 + nx;
    _get_pcood(nx, ny, nz, i, j, k, verts, pc);
    _get_pcood(nx, ny, nz, i, j + 1, k, verts, pc1);
    pp[0] = &_verts(k, j, i, 2, 0);
    pp[1] = &_verts(k, j, i, 3, 0);
    pp[2] = &_verts(k, j, i, 6, 0);
    pp[3] = &_verts(k, j, i, 7, 0);
    for (int i = 0; i < 3; i++)
    {
        fc[i] = (pp[0][i] + pp[1][i] + pp[2][i] + pp[3][i]) / 4;
    }
    D = _get_distance(pc, pp[0], pp[1], pp[2]) + _get_distance(pc1, pp[0], pp[1], pp[2]);
    S = _get_area(pp[0], pp[1], pp[2], pp[3]) / 4;

    Qy = 0.0;
    // 2
    if (k > 0 && j < ny - 1 && i > 0)
    {
        _get_vec(i - 1, j, k - 1, perm, verts, 4, T, idx_8, divn);
        for (int l = 0; l < 8; l++)
        {
            Qy -= T[l] * p[idx_8[l]];
        }
    }
    else
    {
        Qy += (p[i1] - p[i0]) * 2.0 * perm[i0].y * perm[i1].y / (perm[i0].y + perm[i1].y) * S / D;
    }
    // 3
    if (k > 0 && j < ny - 1 && i < nx - 1)
    {
        _get_vec(i, j, k - 1, perm, verts, 7, T, idx_8, divn);
        for (int l = 0; l < 8; l++)
        {
            Qy -= T[l] * p[idx_8[l]];
        }
    }
    else
    {
        Qy += (p[i1] - p[i0]) * 2.0 * perm[i0].y * perm[i1].y / (perm[i0].y + perm[i1].y) * S / D;
    }
    // 6
    if (k < nz - 1 && j < ny - 1 && i > 0)
    {
        _get_vec(i - 1, j, k, perm, verts, 5, T, idx_8, divn);
        for (int l = 0; l < 8; l++)
        {
            Qy -= T[l] * p[idx_8[l]];
        }
    }
    else
    {
        Qy += (p[i1] - p[i0]) * 2.0 * perm[i0].y * perm[i1].y / (perm[i0].y + perm[i1].y) * S / D;
    }
    // 7
    if (k < nz - 1 && j < ny - 1 && i < nx - 1)
    {
        _get_vec(i, j, k, perm, verts, 6, T, idx_8, divn);
        for (int l = 0; l < 8; l++)
        {
            Qy -= T[l] * p[idx_8[l]];
        }
    }
    else
    {
        Qy += (p[i1] - p[i0]) * 2.0 * perm[i0].y * perm[i1].y / (perm[i0].y + perm[i1].y) * S / D;
    }
}

void _get_tolQz(const Ktensor *perm, const double *p, int i, int j, int k,
                const double *verts, double &Qz, int divn)
{
    int nx = nxx * divn;
    int ny = nyy * divn;
    int nz = nzz * divn;
    std::mdspan<const double, std::extents<size_t, -1, -1, -1, 8, 3>> _verts(verts, nz, ny, nx);

    const double *pp[4];
    double idxbeta, idx, idxempty;
    double pc[3], pc1[3], fc[3];
    double S = 0.0, D = 0.0;
    double T[8];
    int idx_8[8];

    int i0 = k * nx * ny + j * nx + i;
    int i1 = i0 + nx * ny;
    _get_pcood(nx, ny, nz, i, j, k, verts, pc);
    _get_pcood(nx, ny, nz, i, j, k + 1, verts, pc1);
    pp[0] = &_verts(k, j, i, 4, 0);
    pp[1] = &_verts(k, j, i, 5, 0);
    pp[2] = &_verts(k, j, i, 6, 0);
    pp[3] = &_verts(k, j, i, 7, 0);
    for (int i = 0; i < 3; i++)
    {
        fc[i] = (pp[0][i] + pp[1][i] + pp[2][i] + pp[3][i]) / 4;
    }
    D = _get_distance(pc, pp[0], pp[1], pp[2]) + _get_distance(pc1, pp[0], pp[1], pp[2]);
    S = _get_area(pp[0], pp[1], pp[2], pp[3]) / 4;

    Qz = 0.0;
    // 4
    if (k < nz - 1 && j > 0 && i > 0)
    {
        _get_vec(i - 1, j - 1, k, perm, verts, 8, T, idx_8, divn);
        for (int l = 0; l < 8; l++)
        {
            Qz -= T[l] * p[idx_8[l]];
        }
    }
    else
    {
        Qz += (p[i1] - p[i0]) * 2.0 * perm[i0].z * perm[i1].z / (perm[i0].z + perm[i1].z) * S / D;
    }
    // 5
    if (k < nz - 1 && j > 0 && i < nx - 1)
    {
        _get_vec(i, j - 1, k, perm, verts, 9, T, idx_8, divn);
        for (int l = 0; l < 8; l++)
        {
            Qz -= T[l] * p[idx_8[l]];
        }
    }
    else
    {
        Qz += (p[i1] - p[i0]) * 2.0 * perm[i0].z * perm[i1].z / (perm[i0].z + perm[i1].z) * S / D;
    }
    // 6
    if (k < nz - 1 && j < ny - 1 && i > 0)
    {
        _get_vec(i - 1, j, k, perm, verts, 11, T, idx_8, divn);
        for (int l = 0; l < 8; l++)
        {
            Qz -= T[l] * p[idx_8[l]];
        }
    }
    else
    {
        Qz += (p[i1] - p[i0]) * 2.0 * perm[i0].z * perm[i1].z / (perm[i0].z + perm[i1].z) * S / D;
    }
    // 7
    if (k < nz - 1 && j < ny - 1 && i < nx - 1)
    {
        _get_vec(i, j, k, perm, verts, 10, T, idx_8, divn);
        for (int l = 0; l < 8; l++)
        {
            Qz -= T[l] * p[idx_8[l]];
        }
    }
    else
    {
        Qz += (p[i1] - p[i0]) * 2.0 * perm[i0].z * perm[i1].z / (perm[i0].z + perm[i1].z) * S / D;
    }
}

void _get_F(const Ktensor *perm, const double *p, int i, int j, int k,
            const double *verts, double &F, int divn)
{
    int nx = nxx * divn;
    int ny = nyy * divn;
    int nz = nzz * divn;
    std::mdspan<const double, std::extents<size_t, -1, -1, -1, 8, 3>> _verts(verts, nz, ny, nx);

    const double *pp[4];
    double idxbeta, idx;
    double idxempty;
    double pc[3], pc1[3];
    double fc[3];
    double S = 0.0, D = 0.0;
    double Q;
    F = 0.0;
    int i0 = k * nx * ny + j * nx + i;
    _get_pcood(nx, ny, nz, i, j, k, verts, pc);
    // Qx
    if (i > 0)
    {
        int i1 = i0 - 1;
        _get_tolQx(perm, p, i - 1, j, k, verts, Q, divn);
        F -= Q;
    }
    else
    {
        pp[0] = &_verts(k, j, i, 0, 0);
        pp[1] = &_verts(k, j, i, 2, 0);
        pp[2] = &_verts(k, j, i, 4, 0);
        pp[3] = &_verts(k, j, i, 6, 0);
        S = _get_area(pp[0], pp[1], pp[2]) + _get_area(pp[3], pp[1], pp[2]);
        D = _get_distance(pc, pp[0], pp[1], pp[2]);
        F += S / D * perm[i0].x * (Plow - p[i0]);
    }

    if (i < nx - 1)
    {
        _get_tolQx(perm, p, i, j, k, verts, Q, divn);
        F += Q;
    }
    else
    {
        pp[0] = &_verts(k, j, i, 1, 0);
        pp[1] = &_verts(k, j, i, 3, 0);
        pp[2] = &_verts(k, j, i, 5, 0);
        pp[3] = &_verts(k, j, i, 7, 0);
        S = _get_area(pp[0], pp[1], pp[2]) + _get_area(pp[3], pp[1], pp[2]);
        D = _get_distance(pc, pp[0], pp[1], pp[2]);
        F += S / D * perm[i0].x * (Phigh - p[i0]);
    }

    // Qy
    if (j > 0)
    {
        int i1 = i0 - nx;
        _get_tolQy(perm, p, i, j - 1, k, verts, Q, divn);
        F -= Q;
    }
    else
    {
        F += 0.0;
    }

    if (j < ny - 1)
    {
        _get_tolQy(perm, p, i, j, k, verts, Q, divn);
        F += Q;
    }
    else
    {
        F += 0.0;
    }

    // Qz
    if (k > 0)
    {
        int i1 = i0 - nx * ny;
        _get_tolQz(perm, p, i, j, k - 1, verts, Q, divn);
        F -= Q;
    }
    else
    {
        F += 0.0;
    }

    if (k < nx - 1)
    {
        _get_tolQz(perm, p, i, j, k, verts, Q, divn);
        F += Q;
    }
    else
    {
        F += 0.0;
    }
    assert(!std::isnan(F));
}

void _get_csr(const Ktensor *perm, double *p, int i0, int j0, int k0,
              const double *verts, double F, int i1, int j1, int k1,
              int *Idx, double *Val, int &nnz, int &count, int divn)
{
    int nx = nxx * divn;
    int ny = nyy * divn;
    int nz = nzz * divn;
    double F1 = 0.0;
    double val = 0.0;
    int ii = k1 * nx * ny + j1 * nx + i1;
    p[ii] += epsilonP;
    _get_F(perm, p, i0, j0, k0, verts, F1, divn);
    // std::cout << F1 << "\n";
    val = (F1 - F) / epsilonP;
    p[ii] -= epsilonP;

    if (val != 0.0)
    {
        Val[nnz] = val;
        Idx[nnz] = ii;
        nnz = nnz + 1;
        count = count + 1;
    }
}
double &C(int n, int m)
{
    for (int idx = Ptr[n]; idx < Ptr[n + 1]; ++idx)
    {
        if (Idx[idx] == m)
        {
            return Val[idx];
        }
    }
    assert(false);
    static double zero = 0.0;
    return zero;
};
void MPFA(int i, int j, int k, const enum Axis axi)
{
    double t[4][3]{};
    double ratio[4]{};
    double pb = 0.0;
    int cur = k * nx * ny + j * nx + i;
    // 边界网格的边界压力
    if (i == 0)
    {
        pb = Plow;
    }
    else if (i == nx)
    {
        pb = Phigh;
    }
    switch (axi)
    {
    case Axis::XPOSITIVE:
        assert(i < nx);

        // 3个面的点，2个未知数,8个点

        // 12条边上其余点的处理
        if (i == 0 && k == 0 && j != 0 && j != ny)
        {
            // 两个网格中心下标
            int i6 = cur, i5 = cur - nx;
            // 0,1,4,8,11交接面中点下标
            double p0[3];
            double p1[3];
            double p4[3];
            double p8[3];
            double p11[3];
            _get_centroid(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 2, 0), &verts1(k, j, i, 4, 0), &verts1(k, j, i, 6, 0), p0);
            _get_centroid(&verts1(k, j - 1, i, 0, 0), &verts1(k, j - 1, i, 2, 0), &verts1(k, j - 1, i, 4, 0), &verts1(k, j - 1, i, 6, 0), p1);
            _get_centroid(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 1, 0), &verts1(k, j, i, 4, 0), &verts1(k, j, i, 5, 0), p4);
            _get_centroid(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 1, 0), &verts1(k, j, i, 2, 0), &verts1(k, j, i, 3, 0), p8);
            _get_centroid(&verts1(k, j - 1, i, 0, 0), &verts1(k, j - 1, i, 1, 0), &verts1(k, j - 1, i, 2, 0), &verts1(k, j - 1, i, 3, 0), p11);
            // 0,1,2,5交接边点下标
            double pp0[3];
            double pp1[3];
            double pp2[3];
            double pp5[3];
            _get_midpoint(&verts1(k, j - 1, i, 0, 0), &verts1(k, j - 1, i, 2, 0), pp0);
            _get_midpoint(&verts1(k, j - 1, i, 2, 0), &verts1(k, j - 1, i, 3, 0), pp1);
            _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 2, 0), pp2);
            _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 4, 0), pp5);
            double n0[3];
            double n1[3];
            double n4[3];
            double n8[3];
            double n11[3];
            _get_surface_normal(&verts1(k, j, i, 0, 0), pp2, p0, pp5, n0, Axis::XPOSITIVE);
            _get_surface_normal(&verts1(k, j - 1, i, 2, 0), pp0, p1, pp5, n1, Axis::XPOSITIVE);
            _get_surface_normal(&verts1(k, j, i, 0, 0), pp5, p4, pp1, n4, Axis::YPOSITIVE);
            _get_surface_normal(&verts1(k, j, i, 0, 0), pp1, p8, pp2, n8, Axis::ZPOSITIVE);
            _get_surface_normal(&verts1(k, j - 1, i, 2, 0), pp0, p11, pp1, n11, Axis::ZPOSITIVE);
            _get_triangle_normal(&verts1(k, j, i, 0, 0), p4, pp1, t[2], Axis::YPOSITIVE);
            ratio[2] = _get_ratio(n4, t[2]);
            Eigen::Matrix<double, 6, 6> A;
            A.setZero();
            A.row(0) << 0, 0, 0, p0[0] - _cx(k, j, i), p0[1] - _cy(k, j, i), p0[2] - _cz(k, j, i);
            A.row(1) << p1[0] - _cx(k, j - 1, i), p1[1] - _cy(k, j - 1, i), p1[2] - _cz(k, j - 1, i), 0, 0, 0;
            A.row(2) << _cx(k, j - 1, i) - p4[0], _cy(k, j - 1, i) - p4[1], _cz(k, j - 1, i) - p4[2], p4[0] - _cx(k, j, i), p4[1] - _cy(k, j, i), p4[2] - _cz(k, j, i);
            Eigen::Matrix3d matK6, matK5;
            _get_matK(pem1[i6], matK6);
            _get_matK(pem1[i5], matK5);
            Eigen::RowVector3d v0, v1, v4, v8, v11;
            v0 << n0[0], n0[1], n0[2];
            v1 << n1[0], n1[1], n1[2];
            v4 << n4[0], n4[1], n4[2];
            v8 << n8[0], n8[1], n8[2];
            v11 << n11[0], n11[1], n11[2];
            A.block(3, 0, 1, 3) = v4 * matK5;
            A.block(3, 3, 1, 3) = -v4 * matK6;
            A.block(4, 3, 1, 3) = v8 * matK6;
            A.block(5, 0, 1, 3) = v11 * matK5;
            A = A.inverse();
            Eigen::RowVectorXd r = ratio[2] * v4 * matK5 * A.topRows(3);
            B[i5] -= pb * (r[0] + r[1]);
            C(i5, i5) += -r[1] + r[2];
            C(i5, i6) += -r[0] - r[2];
            r = -ratio[2] * v4 * matK6 * A.bottomRows(3);
            B[i6] -= pb * (r[0] + r[1]);
            C(i6, i5) += -r[1] + r[2];
            C(i6, i6) += -r[0] - r[2];
        }
        else if (i == 0 && k == nz && j != 0 && j != ny)
        {
            // 两个网格中心下标
            int i2 = cur - nx * ny, i1 = i2 - nx;
            // 2,3,5,8,11交接面中点下标
            double p2[3];
            double p3[3];
            double p5[3];
            double p8[3];
            double p11[3];
            _get_centroid(&verts1(k - 1, j - 1, i, 0, 0), &verts1(k - 1, j - 1, i, 2, 0), &verts1(k - 1, j - 1, i, 4, 0), &verts1(k - 1, j - 1, i, 6, 0), p2);
            _get_centroid(&verts1(k - 1, j, i, 0, 0), &verts1(k - 1, j, i, 2, 0), &verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 6, 0), p3);
            _get_centroid(&verts1(k - 1, j, i, 0, 0), &verts1(k - 1, j, i, 1, 0), &verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 5, 0), p5);
            _get_centroid(&verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 5, 0), &verts1(k - 1, j, i, 6, 0), &verts1(k - 1, j, i, 7, 0), p8);
            _get_centroid(&verts1(k - 1, j - 1, i, 4, 0), &verts1(k - 1, j - 1, i, 5, 0), &verts1(k - 1, j - 1, i, 6, 0), &verts1(k - 1, j - 1, i, 7, 0), p11);
            // 0,1,2,4交接边点下标
            double pp0[3];
            double pp1[3];
            double pp2[3];
            double pp4[3];
            _get_midpoint(&verts1(k - 1, j - 1, i, 4, 0), &verts1(k - 1, j - 1, i, 6, 0), pp0);
            _get_midpoint(&verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 6, 0), pp2);
            _get_midpoint(&verts1(k - 1, j - 1, i, 6, 0), &verts1(k - 1, j - 1, i, 7, 0), pp1);
            _get_midpoint(&verts1(k - 1, j, i, 0, 0), &verts1(k - 1, j, i, 4, 0), pp4);
            double n2[3];
            double n3[3];
            double n5[3];
            double n8[3];
            double n11[3];
            _get_surface_normal(&verts1(k - 1, j - 1, i, 6, 0), pp0, p2, pp4, n2, Axis::XPOSITIVE);
            _get_surface_normal(&verts1(k - 1, j, i, 4, 0), pp2, p3, pp4, n3, Axis::XPOSITIVE);
            _get_surface_normal(&verts1(k - 1, j, i, 4, 0), pp4, p5, pp1, n5, Axis::YPOSITIVE);
            _get_surface_normal(&verts1(k - 1, j, i, 4, 0), pp1, p8, pp2, n8, Axis::ZPOSITIVE);
            _get_surface_normal(&verts1(k - 1, j - 1, i, 6, 0), pp0, p11, pp1, n11, Axis::ZPOSITIVE);
            _get_triangle_normal(&verts1(k - 1, j, i, 4, 0), p5, pp1, t[0], Axis::YPOSITIVE);
            ratio[0] = _get_ratio(n5, t[0]);
            Eigen::Matrix<double, 6, 6> A;
            A.setZero();
            A.row(0) << 0, 0, 0, p3[0] - _cx(k - 1, j, i), p3[1] - _cy(k - 1, j, i), p3[2] - _cz(k - 1, j, i);
            A.row(1) << p2[0] - _cx(k - 1, j - 1, i), p2[1] - _cy(k - 1, j - 1, i), p2[2] - _cz(k - 1, j - 1, i), 0, 0, 0;
            A.row(2) << _cx(k - 1, j - 1, i) - p5[0], _cy(k - 1, j - 1, i) - p5[1], _cz(k - 1, j - 1, i) - p5[2], p5[0] - _cx(k - 1, j, i), p5[1] - _cy(k - 1, j, i), p5[2] - _cz(k - 1, j, i);
            Eigen::Matrix3d matK2, matK1;
            _get_matK(pem1[i2], matK2);
            _get_matK(pem1[i1], matK1);
            Eigen::RowVector3d v2, v3, v5, v8, v11;
            v2 << n2[0], n2[1], n2[2];
            v3 << n3[0], n3[1], n3[2];
            v5 << n5[0], n5[1], n5[2];
            v8 << n8[0], n8[1], n8[2];
            v11 << n11[0], n11[1], n11[2];
            A.block(3, 0, 1, 3) = v5 * matK1;
            A.block(3, 3, 1, 3) = -v5 * matK2;
            A.block(4, 3, 1, 3) = v8 * matK2;
            A.block(5, 0, 1, 3) = v11 * matK1;
            A = A.inverse();
            Eigen::RowVectorXd r = ratio[0] * v5 * matK1 * A.topRows(3);
            B[i1] -= pb * (r[0] + r[1]);
            C(i1, i1) += -r[1] + r[2];
            C(i1, i2) += -r[0] - r[2];
            r = -ratio[0] * v5 * matK2 * A.bottomRows(3);
            B[i2] -= pb * (r[0] + r[1]);
            C(i2, i1) += -r[1] + r[2];
            C(i2, i2) += -r[0] - r[2];
        }
        else if (i == 0 && j == 0 && k != 0 && k != nz)
        {
            int i6 = cur, i2 = cur - nx * ny;
            // 0,3,4,5,8交接面中点下标
            double p0[3];
            double p3[3];
            double p4[3];
            double p5[3];
            double p8[3];
            _get_centroid(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 2, 0), &verts1(k, j, i, 4, 0), &verts1(k, j, i, 6, 0), p0);
            _get_centroid(&verts1(k - 1, j, i, 0, 0), &verts1(k - 1, j, i, 2, 0), &verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 6, 0), p3);
            _get_centroid(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 1, 0), &verts1(k, j, i, 4, 0), &verts1(k, j, i, 5, 0), p4);
            _get_centroid(&verts1(k - 1, j, i, 0, 0), &verts1(k - 1, j, i, 1, 0), &verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 5, 0), p5);
            _get_centroid(&verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 5, 0), &verts1(k - 1, j, i, 6, 0), &verts1(k - 1, j, i, 7, 0), p8);
            // 1,2,4,5交接边点下标
            double pp1[3];
            double pp2[3];
            double pp4[3];
            double pp5[3];
            _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 1, 0), pp1);
            _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 2, 0), pp2);
            _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 4, 0), pp5);
            _get_midpoint(&verts1(k - 1, j, i, 0, 0), &verts1(k - 1, j, i, 4, 0), pp4);
            double n0[3];
            double n3[3];
            double n4[3];
            double n5[3];
            double n8[3];
            _get_surface_normal(&verts1(k, j, i, 0, 0), pp2, p0, pp5, n0, Axis::XPOSITIVE);
            _get_surface_normal(&verts1(k - 1, j, i, 4, 0), pp2, p3, pp4, n3, Axis::XPOSITIVE);
            _get_surface_normal(&verts1(k, j, i, 0, 0), pp5, p4, pp1, n4, Axis::YPOSITIVE);
            _get_surface_normal(&verts1(k - 1, j, i, 4, 0), pp4, p5, pp1, n5, Axis::YPOSITIVE);
            _get_surface_normal(&verts1(k - 1, j, i, 4, 0), pp1, p8, pp2, n8, Axis::ZPOSITIVE);
            _get_triangle_normal(&verts1(k - 1, j, i, 4, 0), pp1, p8, t[1], Axis::ZPOSITIVE);
            ratio[1] = _get_ratio(n8, t[1]);
            Eigen::Matrix<double, 6, 6> A;
            A.setZero();
            A.row(0) << 0, 0, 0, p0[0] - _cx(k, j, i), p0[1] - _cy(k, j, i), p0[2] - _cz(k, j, i);
            A.row(1) << p3[0] - _cx(k - 1, j, i), p3[1] - _cy(k - 1, j, i), p3[2] - _cz(k - 1, j, i), 0, 0, 0;
            A.row(2) << _cx(k - 1, j, i) - p8[0], _cy(k - 1, j, i) - p8[1], _cz(k - 1, j, i) - p8[2], p8[0] - _cx(k, j, i), p8[1] - _cy(k, j, i), p8[2] - _cz(k, j, i);
            Eigen::Matrix3d matK6, matK2;
            _get_matK(pem1[i6], matK6);
            _get_matK(pem1[i2], matK2);
            Eigen::RowVector3d v0, v3, v4, v5, v8;
            v0 << n0[0], n0[1], n0[2];
            v3 << n3[0], n3[1], n3[2];
            v4 << n4[0], n4[1], n4[2];
            v5 << n5[0], n5[1], n5[2];
            v8 << n8[0], n8[1], n8[2];
            A.block(3, 0, 1, 3) = v8 * matK2;
            A.block(3, 3, 1, 3) = -v8 * matK6;
            A.block(4, 3, 1, 3) = v4 * matK6;
            A.block(5, 0, 1, 3) = v5 * matK2;
            A = A.inverse();
            Eigen::RowVectorXd r = ratio[1] * v8 * matK2 * A.topRows(3);
            B[i2] -= pb * (r[0] + r[1]);
            C(i2, i2) += -r[1] + r[2];
            C(i2, i6) += -r[0] - r[2];
            r = -ratio[1] * v8 * matK6 * A.bottomRows(3);
            B[i6] -= pb * (r[0] + r[1]);
            C(i6, i2) += -r[1] + r[2];
            C(i6, i6) += -r[0] - r[2];
            // To be implemented
        }
        else if (i == 0 && j == ny && k != 0 && k != nz)
        {
            int i5 = cur - nx, i1 = i5 - nx * ny;
            // 1,2,4,5,11交接面中点下标
            double p1[3];
            double p2[3];
            double p4[3];
            double p5[3];
            double p11[3];
            _get_centroid(&verts1(k, j - 1, i, 0, 0), &verts1(k, j - 1, i, 2, 0), &verts1(k, j - 1, i, 4, 0), &verts1(k, j - 1, i, 6, 0), p1);
            _get_centroid(&verts1(k - 1, j - 1, i, 0, 0), &verts1(k - 1, j - 1, i, 2, 0), &verts1(k - 1, j - 1, i, 4, 0), &verts1(k - 1, j - 1, i, 6, 0), p2);
            _get_centroid(&verts1(k, j - 1, i, 2, 0), &verts1(k, j - 1, i, 3, 0), &verts1(k, j - 1, i, 6, 0), &verts1(k, j - 1, i, 7, 0), p4);
            _get_centroid(&verts1(k - 1, j - 1, i, 2, 0), &verts1(k - 1, j - 1, i, 3, 0), &verts1(k - 1, j - 1, i, 6, 0), &verts1(k - 1, j - 1, i, 7, 0), p5);
            _get_centroid(&verts1(k - 1, j - 1, i, 4, 0), &verts1(k - 1, j - 1, i, 5, 0), &verts1(k - 1, j - 1, i, 6, 0), &verts1(k - 1, j - 1, i, 7, 0), p11);
            // 0,1,4,5交接边点下标
            double pp0[3];
            double pp1[3];
            double pp4[3];
            double pp5[3];
            _get_midpoint(&verts1(k, j - 1, i, 0, 0), &verts1(k, j - 1, i, 2, 0), pp0);
            _get_midpoint(&verts1(k, j - 1, i, 3, 0), &verts1(k, j - 1, i, 2, 0), pp1);
            _get_midpoint(&verts1(k, j - 1, i, 6, 0), &verts1(k, j - 1, i, 2, 0), pp5);
            _get_midpoint(&verts1(k - 1, j - 1, i, 2, 0), &verts1(k - 1, j - 1, i, 6, 0), pp4);
            double n1[3];
            double n2[3];
            double n4[3];
            double n5[3];
            double n11[3];
            _get_surface_normal(&verts1(k, j - 1, i, 2, 0), pp0, p1, pp5, n1, Axis::XPOSITIVE);
            _get_surface_normal(&verts1(k - 1, j - 1, i, 6, 0), pp4, p2, pp0, n2, Axis::XPOSITIVE);
            _get_surface_normal(&verts1(k, j - 1, i, 2, 0), pp1, p4, pp5, n4, Axis::YPOSITIVE);
            _get_surface_normal(&verts1(k - 1, j - 1, i, 6, 0), pp4, p5, pp1, n5, Axis::YPOSITIVE);
            _get_surface_normal(&verts1(k - 1, j - 1, i, 6, 0), pp0, p11, pp1, n11, Axis::ZPOSITIVE);
            _get_triangle_normal(&verts1(k - 1, j - 1, i, 6, 0), p11, pp1, t[3], Axis::ZPOSITIVE);
            ratio[3] = _get_ratio(n11, t[3]);
            Eigen::Matrix<double, 6, 6> A;
            A.setZero();
            A.row(0) << 0, 0, 0, p1[0] - _cx(k, j - 1, i), p1[1] - _cy(k, j - 1, i), p1[2] - _cz(k, j - 1, i);
            A.row(1) << p2[0] - _cx(k - 1, j - 1, i), p2[1] - _cy(k - 1, j - 1, i), p2[2] - _cz(k - 1, j - 1, i), 0, 0, 0;
            A.row(2) << _cx(k - 1, j - 1, i) - p11[0], _cy(k - 1, j - 1, i) - p11[1], _cz(k - 1, j - 1, i) - p11[2], p11[0] - _cx(k, j - 1, i), p11[1] - _cy(k, j - 1, i), p11[2] - _cz(k, j - 1, i);
            Eigen::Matrix3d matK5, matK1;
            _get_matK(pem1[i5], matK5);
            _get_matK(pem1[i1], matK1);
            Eigen::RowVector3d v1, v2, v4, v5, v11;
            v1 << n1[0], n1[1], n1[2];
            v2 << n2[0], n2[1], n2[2];
            v4 << n4[0], n4[1], n4[2];
            v5 << n5[0], n5[1], n5[2];
            v11 << n11[0], n11[1], n11[2];
            A.block(3, 0, 1, 3) = v11 * matK1;
            A.block(3, 3, 1, 3) = -v11 * matK5;
            A.block(4, 3, 1, 3) = v4 * matK5;
            A.block(5, 0, 1, 3) = v5 * matK1;
            A = A.inverse();
            Eigen::RowVectorXd r = ratio[3] * v11 * matK1 * A.topRows(3);
            B[i1] -= pb * (r[0] + r[1]);
            C(i1, i1) += -r[1] + r[2];
            C(i1, i5) += -r[0] - r[2];
            r = -ratio[3] * v11 * matK5 * A.bottomRows(3);
            B[i5] -= pb * (r[0] + r[1]);
            C(i5, i1) += -r[1] + r[2];
            C(i5, i5) += -r[0] - r[2];
            // To be implemented
        }
        // 六个表面其余点
        else if (i == 0 && j != 0 && j != ny && k != 0 && k != nz)
        {
            int i6 = cur, i5 = i6 - nx, i2 = i6 - nx * ny, i1 = i5 - nx * ny;
            // 0,1,2,3,4,5,8,11交接面中点下标
            double p0[3];
            double p1[3];
            double p2[3];
            double p3[3];
            double p4[3];
            double p5[3];
            double p8[3];
            double p11[3];
            _get_centroid(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 1, 0), &verts1(k, j, i, 2, 0), &verts1(k, j, i, 3, 0), p8);
            _get_centroid(&verts1(k, j - 1, i, 0, 0), &verts1(k, j - 1, i, 1, 0), &verts1(k, j - 1, i, 2, 0), &verts1(k, j - 1, i, 3, 0), p11);
            _get_centroid(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 1, 0), &verts1(k, j, i, 4, 0), &verts1(k, j, i, 5, 0), p4);
            _get_centroid(&verts1(k - 1, j, i, 0, 0), &verts1(k - 1, j, i, 1, 0), &verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 5, 0), p5);
            _get_centroid(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 2, 0), &verts1(k, j, i, 4, 0), &verts1(k, j, i, 6, 0), p0);
            _get_centroid(&verts1(k - 1, j, i, 0, 0), &verts1(k - 1, j, i, 2, 0), &verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 6, 0), p3);
            _get_centroid(&verts1(k, j - 1, i, 0, 0), &verts1(k, j - 1, i, 2, 0), &verts1(k, j - 1, i, 4, 0), &verts1(k, j - 1, i, 6, 0), p1);
            _get_centroid(&verts1(k - 1, j - 1, i, 0, 0), &verts1(k - 1, j - 1, i, 2, 0), &verts1(k - 1, j - 1, i, 4, 0), &verts1(k - 1, j - 1, i, 6, 0), p2);
            // 交接边点下标0,1,2,4,5
            double pp0[3];
            double pp1[3];
            double pp2[3];
            double pp4[3];
            double pp5[3];
            _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 1, 0), pp1);
            _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 2, 0), pp2);
            _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 4, 0), pp5);
            _get_midpoint(&verts1(k - 1, j, i, 0, 0), &verts1(k - 1, j, i, 4, 0), pp4);
            _get_midpoint(&verts1(k, j - 1, i, 0, 0), &verts1(k, j - 1, i, 2, 0), pp0);
            double n0[3];
            double n1[3];
            double n2[3];
            double n3[3];
            double n4[3];
            double n5[3];
            double n8[3];
            double n11[3];
            _get_surface_normal(&verts1(k, j, i, 0, 0), pp2, p0, pp5, n0, Axis::XPOSITIVE);
            _get_surface_normal(&verts1(k, j, i, 0, 0), pp0, p1, pp5, n1, Axis::XPOSITIVE);
            _get_surface_normal(&verts1(k, j, i, 0, 0), pp0, p2, pp4, n2, Axis::XPOSITIVE);
            _get_surface_normal(&verts1(k, j, i, 0, 0), pp4, p3, pp2, n3, Axis::XPOSITIVE);
            _get_surface_normal(&verts1(k, j, i, 0, 0), pp1, p4, pp5, n4, Axis::YPOSITIVE);
            _get_surface_normal(&verts1(k, j, i, 0, 0), pp4, p5, pp1, n5, Axis::YPOSITIVE);
            _get_surface_normal(&verts1(k, j, i, 0, 0), pp2, p8, pp1, n8, Axis::ZPOSITIVE);
            _get_surface_normal(&verts1(k, j, i, 0, 0), pp0, p11, pp1, n11, Axis::ZPOSITIVE);
            _get_triangle_normal(&verts1(k, j, i, 0, 0), p5, pp1, t[0], Axis::YPOSITIVE);
            _get_triangle_normal(&verts1(k, j, i, 0, 0), p8, pp1, t[1], Axis::ZPOSITIVE);
            _get_triangle_normal(&verts1(k, j, i, 0, 0), pp1, p4, t[2], Axis::YPOSITIVE);
            _get_triangle_normal(&verts1(k, j, i, 0, 0), p11, pp1, t[3], Axis::ZPOSITIVE);
            ratio[0] = _get_ratio(n5, t[0]);
            ratio[1] = _get_ratio(n8, t[1]);
            ratio[2] = _get_ratio(n4, t[2]);
            ratio[3] = _get_ratio(n11, t[3]);
            Eigen::Matrix<double, 12, 12> A;
            A.setZero();
            A.block(0, 9, 1, 3) << p0[0] - _cx(k, j, i), p0[1] - _cy(k, j, i), p0[2] - _cz(k, j, i);
            A.block(1, 6, 1, 3) << p1[0] - _cx(k, j - 1, i), p1[1] - _cy(k, j - 1, i), p1[2] - _cz(k, j - 1, i);
            A.block(2, 0, 1, 3) << p2[0] - _cx(k - 1, j - 1, i), p2[1] - _cy(k - 1, j - 1, i), p2[2] - _cz(k - 1, j - 1, i);
            A.block(3, 3, 1, 3) << p3[0] - _cx(k - 1, j, i), p3[1] - _cy(k - 1, j, i), p3[2] - _cz(k - 1, j, i);
            Eigen::Matrix3d matK6, matK5, matK2, matK1;
            _get_matK(pem1[i6], matK6);
            _get_matK(pem1[i5], matK5);
            _get_matK(pem1[i2], matK2);
            _get_matK(pem1[i1], matK1);
            Eigen::RowVector3d v0, v1, v2, v3, v4, v5, v8, v11;
            v0 << n0[0], n0[1], n0[2];
            v1 << n1[0], n1[1], n1[2];
            v2 << n2[0], n2[1], n2[2];
            v3 << n3[0], n3[1], n3[2];
            v4 << n4[0], n4[1], n4[2];
            v5 << n5[0], n5[1], n5[2];
            v8 << n8[0], n8[1], n8[2];
            v11 << n11[0], n11[1], n11[2];
            A.block(4, 6, 1, 6) << _cx(k, j - 1, i) - p4[0], _cy(k, j - 1, i) - p4[1], _cz(k, j - 1, i) - p4[2], p4[0] - _cx(k, j, i), p4[1] - _cy(k, j, i), p4[2] - _cz(k, j, i);
            A.block(5, 6, 1, 3) = v4 * matK5;
            A.block(5, 9, 1, 3) = -v4 * matK6;
            A.block(6, 0, 1, 6) << _cx(k - 1, j - 1, i) - p5[0], _cy(k - 1, j - 1, i) - p5[1], _cz(k - 1, j - 1, i) - p5[2], p5[0] - _cx(k - 1, j, i), p5[1] - _cy(k - 1, j, i), p5[2] - _cz(k - 1, j, i);
            A.block(7, 0, 1, 3) = v5 * matK1;
            A.block(7, 3, 1, 3) = -v5 * matK2;
            A.block(8, 3, 1, 3) << _cx(k - 1, j, i) - p8[0], _cy(k - 1, j, i) - p8[1], _cz(k - 1, j, i) - p8[2];
            A.block(8, 9, 1, 3) << p8[0] - _cx(k, j, i), p8[1] - _cy(k, j, i), p8[2] - _cz(k, j, i);
            A.block(9, 3, 1, 3) = v8 * matK2;
            A.block(9, 9, 1, 3) = -v8 * matK6;
            A.block(10, 0, 1, 3) << _cx(k - 1, j - 1, i) - p11[0], _cy(k - 1, j - 1, i) - p11[1], _cz(k - 1, j - 1, i) - p11[2];
            A.block(10, 6, 1, 3) << p11[0] - _cx(k, j - 1, i), p11[1] - _cy(k, j - 1, i), p11[2] - _cz(k, j - 1, i);
            A.block(11, 0, 1, 3) = v11 * matK1;
            A.block(11, 6, 1, 3) = -v11 * matK5;
            A = A.inverse();
            Eigen::RowVectorXd r = (ratio[0] * v5 + ratio[3] * v11) * matK1 * A.topRows(3);
            B[i1] -= pb * (r[0] + r[1] + r[2] + r[3]);
            C(i1, i1) += -r[2] + r[6] + r[10];
            C(i1, i2) += -r[3] - r[6] + r[8];
            C(i1, i5) += -r[1] + r[4] - r[10];
            C(i1, i6) += -r[0] - r[4] - r[8];
            r = (-ratio[0] * v5 + ratio[1] * v8) * matK2 * A.block(3, 0, 3, 12);
            B[i2] -= pb * (r[0] + r[1] + r[2] + r[3]);
            C(i2, i1) += -r[2] + r[6] + r[10];
            C(i2, i2) += -r[3] - r[6] + r[8];
            C(i2, i5) += -r[1] + r[4] - r[10];
            C(i2, i6) += -r[0] - r[4] - r[8];
            r = (ratio[2] * v4 - ratio[3] * v11) * matK5 * A.block(6, 0, 3, 12);
            B[i5] -= pb * (r[0] + r[1] + r[2] + r[3]);
            C(i5, i1) += -r[2] + r[6] + r[10];
            C(i5, i2) += -r[3] - r[6] + r[8];
            C(i5, i5) += -r[1] + r[4] - r[10];
            C(i5, i6) += -r[0] - r[4] - r[8];
            r = -(ratio[2] * v4 + ratio[1] * v8) * matK6 * A.bottomRows(3);
            B[i6] -= pb * (r[0] + r[1] + r[2] + r[3]);
            C(i6, i1) += -r[2] + r[6] + r[10];
            C(i6, i2) += -r[3] - r[6] + r[8];
            C(i6, i5) += -r[1] + r[4] - r[10];
            C(i6, i6) += -r[0] - r[4] - r[8];
            //  To be implemented
        }
        else if (j == 0 && i != 0 && k != 0 && k != nz)
        {
            int i6 = cur, i7 = i6 - 1, i2 = i6 - nx * ny, i3 = i7 - nx * ny;
            // 0,3,4,5,6,7,8,9交接面中点下标
            double p0[3];
            double p3[3];
            double p4[3];
            double p5[3];
            double p6[3];
            double p7[3];
            double p8[3];
            double p9[3];
            _get_centroid(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 2, 0), &verts1(k, j, i, 4, 0), &verts1(k, j, i, 6, 0), p0);
            _get_centroid(&verts1(k - 1, j, i, 0, 0), &verts1(k - 1, j, i, 2, 0), &verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 6, 0), p3);
            _get_centroid(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 1, 0), &verts1(k, j, i, 2, 0), &verts1(k, j, i, 3, 0), p8);
            _get_centroid(&verts1(k, j, i - 1, 0, 0), &verts1(k, j, i - 1, 1, 0), &verts1(k, j, i - 1, 2, 0), &verts1(k, j, i - 1, 3, 0), p9);
            _get_centroid(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 1, 0), &verts1(k, j, i, 4, 0), &verts1(k, j, i, 5, 0), p4);
            _get_centroid(&verts1(k - 1, j, i, 0, 0), &verts1(k - 1, j, i, 1, 0), &verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 5, 0), p5);
            _get_centroid(&verts1(k - 1, j, i - 1, 0, 0), &verts1(k - 1, j, i - 1, 1, 0), &verts1(k - 1, j, i - 1, 4, 0), &verts1(k - 1, j, i - 1, 5, 0), p6);
            _get_centroid(&verts1(k, j, i - 1, 0, 0), &verts1(k, j, i - 1, 1, 0), &verts1(k, j, i - 1, 4, 0), &verts1(k, j, i - 1, 5, 0), p7);
            // 1,2,3,4,5交接边点下
            double pp1[3];
            double pp2[3];
            double pp3[3];
            double pp4[3];
            double pp5[3];
            _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 2, 0), pp2);
            _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 1, 0), pp1);
            _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 4, 0), pp5);
            _get_midpoint(&verts1(k - 1, j, i, 0, 0), &verts1(k - 1, j, i, 4, 0), pp4);
            _get_midpoint(&verts1(k, j, i - 1, 0, 0), &verts1(k, j, i - 1, 1, 0), pp3);
            double n0[3];
            double n3[3];
            double n4[3];
            double n5[3];
            double n6[3];
            double n7[3];
            double n8[3];
            double n9[3];
            _get_surface_normal(&verts1(k, j, i, 0, 0), pp2, p0, pp5, n0, Axis::XPOSITIVE);
            _get_surface_normal(&verts1(k, j, i, 0, 0), pp2, p3, pp4, n3, Axis::XPOSITIVE);
            _get_surface_normal(&verts1(k, j, i, 0, 0), pp1, p4, pp5, n4, Axis::YPOSITIVE);
            _get_surface_normal(&verts1(k, j, i, 0, 0), pp4, p5, pp1, n5, Axis::YPOSITIVE);
            _get_surface_normal(&verts1(k, j, i, 0, 0), pp4, p6, pp3, n6, Axis::YPOSITIVE);
            _get_surface_normal(&verts1(k, j, i, 0, 0), pp5, p7, pp3, n7, Axis::YPOSITIVE);
            _get_surface_normal(&verts1(k, j, i, 0, 0), pp2, p8, pp1, n8, Axis::ZPOSITIVE);
            _get_surface_normal(&verts1(k, j, i, 0, 0), pp2, p9, pp3, n9, Axis::ZPOSITIVE);
            _get_triangle_normal(&verts1(k, j, i, 0, 0), p8, pp1, t[1], Axis::ZPOSITIVE);
            ratio[1] = _get_ratio(n8, t[1]);
            Eigen::Matrix<double, 12, 12> A;
            A.setZero();
            Eigen::Matrix3d matK6, matK7, matK2, matK3;
            _get_matK(pem1[i6], matK6);
            _get_matK(pem1[i7], matK7);
            _get_matK(pem1[i2], matK2);
            _get_matK(pem1[i3], matK3);
            Eigen::RowVector3d v0, v3, v4, v5, v6, v7, v8, v9;
            v0 << n0[0], n0[1], n0[2];
            v3 << n3[0], n3[1], n3[2];
            v4 << n4[0], n4[1], n4[2];
            v5 << n5[0], n5[1], n5[2];
            v6 << n6[0], n6[1], n6[2];
            v7 << n7[0], n7[1], n7[2];
            v8 << n8[0], n8[1], n8[2];
            v9 << n9[0], n9[1], n9[2];
            A.block(0, 6, 1, 6) << _cx(k, j, i) - p0[0], _cy(k, j, i) - p0[1], _cz(k, j, i) - p0[2], p0[0] - _cx(k, j, i - 1), p0[1] - _cy(k, j, i - 1), p0[2] - _cz(k, j, i - 1);
            A.block(1, 6, 1, 3) = v0 * matK6;
            A.block(1, 9, 1, 3) = -v0 * matK7;
            A.block(2, 0, 1, 6) << _cx(k - 1, j, i) - p3[0], _cy(k - 1, j, i) - p3[1], _cz(k - 1, j, i) - p3[2], p3[0] - _cx(k - 1, j, i - 1), p3[1] - _cy(k - 1, j, i - 1), p3[2] - _cz(k - 1, j, i - 1);
            A.block(3, 0, 1, 3) = v3 * matK2;
            A.block(3, 3, 1, 3) = -v3 * matK3;
            A.block(4, 6, 1, 3) = v4 * matK6;
            A.block(5, 0, 1, 3) = v5 * matK2;
            A.block(6, 3, 1, 3) = v6 * matK3;
            A.block(7, 9, 1, 3) = v7 * matK7;
            A.block(8, 0, 1, 3) << _cx(k - 1, j, i) - p8[0], _cy(k - 1, j, i) - p8[1], _cz(k - 1, j, i) - p8[2];
            A.block(8, 6, 1, 3) << p8[0] - _cx(k, j, i), p8[1] - _cy(k, j, i), p8[2] - _cz(k, j, i);
            A.block(9, 0, 1, 3) = v8 * matK2;
            A.block(9, 6, 1, 3) = -v8 * matK6;
            A.block(10, 3, 1, 3) << _cx(k - 1, j, i - 1) - p9[0], _cy(k - 1, j, i - 1) - p9[1], _cz(k - 1, j, i - 1) - p9[2];
            A.block(10, 9, 1, 3) << p9[0] - _cx(k, j, i - 1), p9[1] - _cy(k, j, i - 1), p9[2] - _cz(k, j, i - 1);
            A.block(11, 3, 1, 3) = v9 * matK3;
            A.block(11, 9, 1, 3) = -v9 * matK7;
            A = A.inverse();
            Eigen::RowVectorXd r = (ratio[1] * v8) * matK2 * A.topRows(3);
            C(i2, i2) += r[2] + r[8];
            C(i2, i3) += -r[2] + r[10];
            C(i2, i6) += r[0] - r[8];
            C(i2, i7) += -r[0] - r[10];
            r = (-ratio[1] * v8) * matK6 * A.block(6, 0, 3, 12);
            C(i6, i2) += r[2] + r[8];
            C(i6, i3) += -r[2] + r[10];
            C(i6, i6) += r[0] - r[8];
            C(i6, i7) += -r[0] - r[10];
            // To be implemented
        }
        else if (j == ny && i != 0 && k != 0 && k != nz)
        {
            int i5 = cur - nx, i4 = i5 - 1, i1 = i5 - nx * ny, i0 = i4 - nx * ny;
            // 1,2,4,5,6,7,10,11交接面中点下标
            double p1[3];
            double p2[3];
            double p4[3];
            double p5[3];
            double p6[3];
            double p7[3];
            double p10[3];
            double p11[3];
            _get_centroid(&verts1(k, j - 1, i, 2, 0), &verts1(k, j - 1, i, 0, 0), &verts1(k, j - 1, i, 4, 0), &verts1(k, j - 1, i, 6, 0), p1);
            _get_centroid(&verts1(k - 1, j - 1, i, 2, 0), &verts1(k - 1, j - 1, i, 0, 0), &verts1(k - 1, j - 1, i, 4, 0), &verts1(k - 1, j - 1, i, 6, 0), p2);
            _get_centroid(&verts1(k, j - 1, i, 2, 0), &verts1(k, j - 1, i, 3, 0), &verts1(k, j - 1, i, 6, 0), &verts1(k, j - 1, i, 7, 0), p4);
            _get_centroid(&verts1(k, j - 1, i - 1, 2, 0), &verts1(k, j - 1, i - 1, 3, 0), &verts1(k, j - 1, i - 1, 6, 0), &verts1(k, j - 1, i - 1, 7, 0), p7);
            _get_centroid(&verts1(k - 1, j - 1, i, 2, 0), &verts1(k - 1, j - 1, i, 3, 0), &verts1(k - 1, j - 1, i, 6, 0), &verts1(k - 1, j - 1, i, 7, 0), p5);
            _get_centroid(&verts1(k - 1, j - 1, i - 1, 2, 0), &verts1(k - 1, j - 1, i - 1, 3, 0), &verts1(k - 1, j - 1, i - 1, 6, 0), &verts1(k - 1, j - 1, i - 1, 7, 0), p6);
            _get_centroid(&verts1(k, j - 1, i, 2, 0), &verts1(k, j - 1, i, 0, 0), &verts1(k, j - 1, i, 1, 0), &verts1(k, j - 1, i, 3, 0), p11);
            _get_centroid(&verts1(k, j - 1, i - 1, 2, 0), &verts1(k, j - 1, i - 1, 0, 0), &verts1(k, j - 1, i - 1, 1, 0), &verts1(k, j - 1, i - 1, 3, 0), p10);
            // 0,1,3,4,5交接边点下标
            double pp0[3];
            double pp1[3];
            double pp3[3];
            double pp4[3];
            double pp5[3];
            _get_midpoint(&verts1(k, j - 1, i, 2, 0), &verts1(k, j - 1, i, 0, 0), pp0);
            _get_midpoint(&verts1(k, j - 1, i, 2, 0), &verts1(k, j - 1, i, 3, 0), pp1);
            _get_midpoint(&verts1(k, j - 1, i, 2, 0), &verts1(k, j - 1, i, 6, 0), pp5);
            _get_midpoint(&verts1(k, j - 1, i - 1, 2, 0), &verts1(k, j - 1, i - 1, 3, 0), pp3);
            _get_midpoint(&verts1(k - 1, j - 1, i, 2, 0), &verts1(k - 1, j - 1, i, 6, 0), pp4);
            double n1[3];
            double n2[3];
            double n4[3];
            double n5[3];
            double n6[3];
            double n7[3];
            double n10[3];
            double n11[3];
            _get_surface_normal(&verts1(k, j - 1, i, 2, 0), pp0, p1, pp5, n1, Axis::XPOSITIVE);
            _get_surface_normal(&verts1(k, j - 1, i, 2, 0), pp0, p2, pp4, n2, Axis::XPOSITIVE);
            _get_surface_normal(&verts1(k, j - 1, i, 2, 0), pp1, p4, pp5, n4, Axis::YPOSITIVE);
            _get_surface_normal(&verts1(k, j - 1, i, 2, 0), pp4, p5, pp1, n5, Axis::YPOSITIVE);
            _get_surface_normal(&verts1(k, j - 1, i, 2, 0), pp4, p6, pp3, n6, Axis::YPOSITIVE);
            _get_surface_normal(&verts1(k, j - 1, i, 2, 0), pp5, p7, pp3, n7, Axis::YPOSITIVE);
            _get_surface_normal(&verts1(k, j - 1, i, 2, 0), pp0, p10, pp3, n10, Axis::ZPOSITIVE);
            _get_surface_normal(&verts1(k, j - 1, i, 2, 0), pp0, p11, pp1, n11, Axis::ZPOSITIVE);
            _get_triangle_normal(&verts1(k, j - 1, i, 2, 0), p11, pp1, t[3], Axis::ZPOSITIVE);
            ratio[3] = _get_ratio(n11, t[3]);
            Eigen::Matrix<double, 12, 12> A;
            A.setZero();
            Eigen::Matrix3d matK5, matK4, matK1, matK0;
            _get_matK(pem1[i5], matK5);
            _get_matK(pem1[i4], matK4);
            _get_matK(pem1[i1], matK1);
            _get_matK(pem1[i0], matK0);
            Eigen::RowVector3d v1, v2, v4, v5, v6, v7, v10, v11;
            v1 << n1[0], n1[1], n1[2];
            v2 << n2[0], n2[1], n2[2];
            v4 << n4[0], n4[1], n4[2];
            v5 << n5[0], n5[1], n5[2];
            v6 << n6[0], n6[1], n6[2];
            v7 << n7[0], n7[1], n7[2];
            v10 << n10[0], n10[1], n10[2];
            v11 << n11[0], n11[1], n11[2];
            A.block(0, 6, 1, 6) << _cx(k, j - 1, i - 1) - p1[0], _cy(k, j - 1, i - 1) - p1[1], _cz(k, j - 1, i - 1) - p1[2], p1[0] - _cx(k, j - 1, i), p1[1] - _cy(k, j - 1, i), p1[2] - _cz(k, j - 1, i);
            A.block(1, 6, 1, 3) = v1 * matK4;
            A.block(1, 9, 1, 3) = -v1 * matK5;
            A.block(2, 0, 1, 6) << _cx(k - 1, j - 1, i - 1) - p2[0], _cy(k - 1, j - 1, i - 1) - p2[1], _cz(k - 1, j - 1, i - 1) - p2[2], p2[0] - _cx(k - 1, j - 1, i), p2[1] - _cy(k - 1, j - 1, i), p2[2] - _cz(k - 1, j - 1, i);
            A.block(3, 0, 1, 3) = v2 * matK0;
            A.block(3, 3, 1, 3) = -v2 * matK1;
            A.block(4, 9, 1, 3) = v4 * matK5;
            A.block(5, 3, 1, 3) = v5 * matK1;
            A.block(6, 0, 1, 3) = v6 * matK0;
            A.block(7, 6, 1, 3) = v7 * matK4;
            A.block(8, 0, 1, 3) << _cx(k - 1, j - 1, i - 1) - p10[0], _cy(k - 1, j - 1, i - 1) - p10[1], _cz(k - 1, j - 1, i - 1) - p10[2];
            A.block(8, 6, 1, 3) << p10[0] - _cx(k, j - 1, i - 1), p10[1] - _cy(k, j - 1, i - 1), p10[2] - _cz(k, j - 1, i - 1);
            A.block(9, 0, 1, 3) = v10 * matK0;
            A.block(9, 6, 1, 3) = -v10 * matK4;
            A.block(10, 3, 1, 3) << _cx(k - 1, j - 1, i) - p11[0], _cy(k - 1, j - 1, i) - p11[1], _cz(k - 1, j - 1, i) - p11[2];
            A.block(10, 9, 1, 3) << p11[0] - _cx(k, j - 1, i), p11[1] - _cy(k, j - 1, i), p11[2] - _cz(k, j - 1, i);
            A.block(11, 3, 1, 3) = v11 * matK1;
            A.block(11, 9, 1, 3) = -v11 * matK5;
            A = A.inverse();
            Eigen::RowVectorXd r = (ratio[3] * v11) * matK1 * A.block(3, 0, 3, 12);
            C(i1, i0) += r[2] + r[8];
            C(i1, i1) += -r[2] + r[10];
            C(i1, i4) += r[0] - r[8];
            C(i1, i5) += -r[0] - r[10];
            r = (-ratio[3] * v11) * matK5 * A.bottomRows(3);
            C(i5, i0) += r[2] + r[8];
            C(i5, i1) += -r[2] + r[10];
            C(i5, i4) += r[0] - r[8];
            C(i5, i5) += -r[0] - r[10];
            // To be implemented
        }
        else if (k == 0 && i != 0 && j != 0 && j != ny)
        {
            int i6 = cur, i7 = cur - 1, i5 = cur - nx, i4 = i5 - 1;
            // 0,1,4,7,8,9,10,11交接面中点下标
            double p0[3];
            double p1[3];
            double p4[3];
            double p7[3];
            double p8[3];
            double p9[3];
            double p10[3];
            double p11[3];
            _get_centroid(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 2, 0), &verts1(k, j, i, 4, 0), &verts1(k, j, i, 6, 0), p0);
            _get_centroid(&verts1(k, j - 1, i, 0, 0), &verts1(k, j - 1, i, 2, 0), &verts1(k, j - 1, i, 4, 0), &verts1(k, j - 1, i, 6, 0), p1);
            _get_centroid(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 1, 0), &verts1(k, j, i, 4, 0), &verts1(k, j, i, 5, 0), p4);
            _get_centroid(&verts1(k, j, i - 1, 0, 0), &verts1(k, j, i - 1, 1, 0), &verts1(k, j, i - 1, 4, 0), &verts1(k, j, i - 1, 5, 0), p7);
            _get_centroid(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 1, 0), &verts1(k, j, i, 2, 0), &verts1(k, j, i, 3, 0), p8);
            _get_centroid(&verts1(k, j - 1, i, 0, 0), &verts1(k, j - 1, i, 1, 0), &verts1(k, j - 1, i, 2, 0), &verts1(k, j - 1, i, 3, 0), p11);
            _get_centroid(&verts1(k, j, i - 1, 0, 0), &verts1(k, j, i - 1, 1, 0), &verts1(k, j, i - 1, 2, 0), &verts1(k, j, i - 1, 3, 0), p9);
            _get_centroid(&verts1(k, j - 1, i - 1, 0, 0), &verts1(k, j - 1, i - 1, 1, 0), &verts1(k, j - 1, i - 1, 2, 0), &verts1(k, j - 1, i - 1, 3, 0), p10);
            // 0,1,2,3,5 交接边点下标
            double pp0[3];
            double pp1[3];
            double pp2[3];
            double pp3[3];
            double pp5[3];
            _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 2, 0), pp2);
            _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 1, 0), pp1);
            _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 4, 0), pp5);
            _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j - 1, i, 0, 0), pp0);
            _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j, i - 1, 0, 0), pp3);
            double n0[3];
            double n1[3];
            double n4[3];
            double n7[3];
            double n8[3];
            double n9[3];
            double n10[3];
            double n11[3];
            _get_surface_normal(&verts1(k, j, i, 0, 0), pp2, p0, pp5, n0, Axis::XPOSITIVE);
            _get_surface_normal(&verts1(k, j, i, 0, 0), pp5, p1, pp0, n1, Axis::XPOSITIVE);
            _get_surface_normal(&verts1(k, j, i, 0, 0), pp1, p4, pp5, n4, Axis::YPOSITIVE);
            _get_surface_normal(&verts1(k, j, i, 0, 0), pp5, p7, pp3, n7, Axis::YPOSITIVE);
            _get_surface_normal(&verts1(k, j, i, 0, 0), pp2, p8, pp1, n8, Axis::ZPOSITIVE);
            _get_surface_normal(&verts1(k, j, i, 0, 0), pp0, p11, pp1, n11, Axis::ZPOSITIVE);
            _get_surface_normal(&verts1(k, j, i, 0, 0), pp0, p10, pp3, n10, Axis::ZPOSITIVE);
            _get_surface_normal(&verts1(k, j, i, 0, 0), pp3, p9, pp2, n9, Axis::ZPOSITIVE);
            _get_triangle_normal(&verts1(k, j, i, 0, 0), pp1, p4, t[2], Axis::YPOSITIVE);
            ratio[2] = _get_ratio(n4, t[2]);
            Eigen::Matrix<double, 12, 12> A;
            A.setZero();
            Eigen::Matrix3d matK6, matK7, matK5, matK4;
            _get_matK(pem1[i6], matK6);
            _get_matK(pem1[i7], matK7);
            _get_matK(pem1[i5], matK5);
            _get_matK(pem1[i4], matK4);
            Eigen::RowVector3d v0, v1, v4, v7, v8, v9, v10, v11;
            v0 << n0[0], n0[1], n0[2];
            v1 << n1[0], n1[1], n1[2];
            v4 << n4[0], n4[1], n4[2];
            v7 << n7[0], n7[1], n7[2];
            v8 << n8[0], n8[1], n8[2];
            v9 << n9[0], n9[1], n9[2];
            v10 << n10[0], n10[1], n10[2];
            v11 << n11[0], n11[1], n11[2];
            A.block(0, 6, 1, 6) << _cx(k, j, i) - p0[0], _cy(k, j, i) - p0[1], _cz(k, j, i) - p0[2], p0[0] - _cx(k, j, i - 1), p0[1] - _cy(k, j, i - 1), p0[2] - _cz(k, j, i - 1);
            A.block(1, 6, 1, 3) = v0 * matK6;
            A.block(1, 9, 1, 3) = -v0 * matK7;
            A.block(2, 0, 1, 6) << p1[0] - _cx(k, j - 1, i - 1), p1[1] - _cy(k, j - 1, i - 1), p1[2] - _cz(k, j - 1, i - 1), _cx(k, j - 1, i) - p1[0], _cy(k, j - 1, i) - p1[1], _cz(k, j - 1, i) - p1[2];
            A.block(3, 0, 1, 3) = v1 * matK4;
            A.block(3, 3, 1, 3) = -v1 * matK5;
            A.block(4, 3, 1, 6) << _cx(k, j - 1, i) - p4[0], _cy(k, j - 1, i) - p4[1], _cz(k, j - 1, i) - p4[2], p4[0] - _cx(k, j, i), p4[1] - _cy(k, j, i), p4[2] - _cz(k, j, i);
            A.block(5, 3, 1, 3) = v4 * matK5;
            A.block(5, 6, 1, 3) = -v4 * matK6;
            A.block(6, 0, 1, 3) << _cx(k, j - 1, i - 1) - p7[0], _cy(k, j - 1, i - 1) - p7[1], _cz(k, j - 1, i - 1) - p7[2];
            A.block(6, 9, 1, 3) << p7[0] - _cx(k, j, i - 1), p7[1] - _cy(k, j, i - 1), p7[2] - _cz(k, j, i - 1);
            A.block(7, 0, 1, 3) = v7 * matK4;
            A.block(7, 9, 1, 3) = -v7 * matK7;
            A.block(8, 6, 1, 3) = v8 * matK6;
            A.block(9, 9, 1, 3) = v9 * matK7;
            A.block(10, 0, 1, 3) = v10 * matK4;
            A.block(11, 3, 1, 3) = v11 * matK5;
            A = A.inverse();
            Eigen::RowVectorXd r = (ratio[2] * v4) * matK5 * A.block(3, 0, 3, 12);
            C(i5, i4) += -r[2] + r[6];
            C(i5, i5) += r[2] + r[4];
            C(i5, i6) += r[0] - r[4];
            C(i5, i7) += -r[0] - r[6];
            r = (-ratio[2] * v4) * matK6 * A.block(6, 0, 3, 12);
            C(i6, i4) += -r[2] + r[6];
            C(i6, i5) += r[2] + r[4];
            C(i6, i6) += r[0] - r[4];
            C(i6, i7) += -r[0] - r[6];
            // To be implemented
        }
        else if (k == nz && i != 0 && j != 0 && j != ny)
        {
            int i2 = cur - nx * ny, i3 = i2 - 1, i1 = i2 - nx, i0 = i1 - 1;
            // 2,3,5,6,8,9,10,11交接面中点下标
            double p2[3];
            double p3[3];
            double p5[3];
            double p6[3];
            double p8[3];
            double p9[3];
            double p10[3];
            double p11[3];
            _get_centroid(&verts1(k - 1, j, i, 2, 0), &verts1(k - 1, j, i, 0, 0), &verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 6, 0), p3);
            _get_centroid(&verts1(k - 1, j - 1, i, 2, 0), &verts1(k - 1, j - 1, i, 0, 0), &verts1(k - 1, j - 1, i, 4, 0), &verts1(k - 1, j - 1, i, 6, 0), p2);
            _get_centroid(&verts1(k - 1, j, i, 0, 0), &verts1(k - 1, j, i, 1, 0), &verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 5, 0), p5);
            _get_centroid(&verts1(k - 1, j, i - 1, 0, 0), &verts1(k - 1, j, i - 1, 1, 0), &verts1(k - 1, j, i - 1, 4, 0), &verts1(k - 1, j, i - 1, 5, 0), p6);
            _get_centroid(&verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 5, 0), &verts1(k - 1, j, i, 6, 0), &verts1(k - 1, j, i, 7, 0), p8);
            _get_centroid(&verts1(k - 1, j, i - 1, 4, 0), &verts1(k - 1, j, i - 1, 5, 0), &verts1(k - 1, j, i - 1, 6, 0), &verts1(k - 1, j, i - 1, 7, 0), p9);
            _get_centroid(&verts1(k - 1, j - 1, i, 4, 0), &verts1(k - 1, j - 1, i, 5, 0), &verts1(k - 1, j - 1, i, 6, 0), &verts1(k - 1, j - 1, i, 7, 0), p11);
            _get_centroid(&verts1(k - 1, j - 1, i - 1, 4, 0), &verts1(k - 1, j - 1, i - 1, 5, 0), &verts1(k - 1, j - 1, i - 1, 6, 0), &verts1(k - 1, j - 1, i - 1, 7, 0), p10);
            // 0,1,2,3,4 交接边点下标
            double pp0[3];
            double pp1[3];
            double pp2[3];
            double pp3[3];
            double pp4[3];
            _get_midpoint(&verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 6, 0), pp2);
            _get_midpoint(&verts1(k - 1, j - 1, i, 4, 0), &verts1(k - 1, j - 1, i, 6, 0), pp0);
            _get_midpoint(&verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 5, 0), pp1);
            _get_midpoint(&verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i - 1, 4, 0), pp3);
            _get_midpoint(&verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 0, 0), pp4);
            double n2[3];
            double n3[3];
            double n5[3];
            double n6[3];
            double n8[3];
            double n9[3];
            double n10[3];
            double n11[3];
            _get_surface_normal(&verts1(k - 1, j, i, 4, 0), pp0, p2, pp4, n2, Axis::XPOSITIVE);
            _get_surface_normal(&verts1(k - 1, j, i, 4, 0), pp4, p3, pp2, n3, Axis::XPOSITIVE);
            _get_surface_normal(&verts1(k - 1, j, i, 4, 0), pp1, p5, pp4, n5, Axis::YPOSITIVE);
            _get_surface_normal(&verts1(k - 1, j, i, 4, 0), pp4, p6, pp3, n6, Axis::YPOSITIVE);
            _get_surface_normal(&verts1(k - 1, j, i, 4, 0), pp2, p8, pp1, n8, Axis::ZPOSITIVE);
            _get_surface_normal(&verts1(k - 1, j, i, 4, 0), pp3, p9, pp2, n9, Axis::ZPOSITIVE);
            _get_surface_normal(&verts1(k - 1, j, i, 4, 0), pp0, p11, pp1, n11, Axis::ZPOSITIVE);
            _get_surface_normal(&verts1(k - 1, j, i, 4, 0), pp0, p10, pp3, n10, Axis::ZPOSITIVE);
            _get_triangle_normal(&verts1(k - 1, j, i, 4, 0), p5, pp1, t[0], Axis::YPOSITIVE);
            ratio[0] = _get_ratio(n5, t[0]);
            Eigen::Matrix<double, 12, 12> A;
            A.setZero();
            Eigen::Matrix3d matK2, matK3, matK1, matK0;
            _get_matK(pem1[i2], matK2);
            _get_matK(pem1[i3], matK3);
            _get_matK(pem1[i1], matK1);
            _get_matK(pem1[i0], matK0);
            Eigen::RowVector3d v2, v3, v5, v6, v8, v9, v10, v11;
            v2 << n2[0], n2[1], n2[2];
            v3 << n3[0], n3[1], n3[2];
            v5 << n5[0], n5[1], n5[2];
            v6 << n6[0], n6[1], n6[2];
            v8 << n8[0], n8[1], n8[2];
            v9 << n9[0], n9[1], n9[2];
            v10 << n10[0], n10[1], n10[2];
            v11 << n11[0], n11[1], n11[2];
            A.block(0, 0, 1, 6) << _cx(k - 1, j - 1, i - 1) - p2[0], _cy(k - 1, j - 1, i - 1) - p2[1], _cz(k - 1, j - 1, i - 1) - p2[2], p2[0] - _cx(k - 1, j - 1, i), p2[1] - _cy(k - 1, j - 1, i), p2[2] - _cz(k - 1, j - 1, i);
            A.block(1, 0, 1, 3) = v2 * matK0;
            A.block(1, 3, 1, 3) = -v2 * matK1;
            A.block(2, 6, 1, 6) << p3[0] - _cx(k - 1, j, i), p3[1] - _cy(k - 1, j, i), p3[2] - _cz(k - 1, j, i), _cx(k - 1, j, i - 1) - p3[0], _cy(k - 1, j, i - 1) - p3[1], _cz(k - 1, j, i - 1) - p3[2];
            A.block(3, 6, 1, 3) = v3 * matK2;
            A.block(3, 9, 1, 3) = -v3 * matK3;
            A.block(4, 3, 1, 6) << _cx(k - 1, j - 1, i) - p5[0], _cy(k - 1, j - 1, i) - p5[1], _cz(k - 1, j - 1, i) - p5[2], p5[0] - _cx(k - 1, j, i), p5[1] - _cy(k - 1, j, i), p5[2] - _cz(k - 1, j, i);
            A.block(5, 3, 1, 3) = v5 * matK1;
            A.block(5, 6, 1, 3) = -v5 * matK2;
            A.block(6, 0, 1, 3) << _cx(k - 1, j - 1, i - 1) - p6[0], _cy(k - 1, j - 1, i - 1) - p6[1], _cz(k - 1, j - 1, i - 1) - p6[2];
            A.block(6, 9, 1, 3) << p6[0] - _cx(k - 1, j, i - 1), p6[1] - _cy(k - 1, j, i - 1), p6[2] - _cz(k - 1, j, i - 1);
            A.block(7, 0, 1, 3) = v6 * matK0;
            A.block(7, 9, 1, 3) = -v6 * matK3;
            A.block(8, 6, 1, 3) = v8 * matK2;
            A.block(9, 9, 1, 3) = v9 * matK3;
            A.block(10, 0, 1, 3) = v10 * matK0;
            A.block(11, 3, 1, 3) = v11 * matK1;
            A = A.inverse();
            Eigen::RowVectorXd r = (ratio[0] * v5) * matK1 * A.block(3, 0, 3, 12);
            C(i1, i0) += r[0] + r[6];
            C(i1, i1) += -r[0] + r[4];
            C(i1, i2) += -r[2] - r[4];
            C(i1, i3) += r[2] - r[6];
            r = (-ratio[0] * v5) * matK2 * A.block(6, 0, 3, 12);
            C(i2, i0) += r[0] + r[6];
            C(i2, i1) += -r[0] + r[4];
            C(i2, i2) += -r[2] - r[4];
            C(i2, i3) += r[2] - r[6];
            //  To be implemented
        }
        else if (i != 0 && j != 0 && j != ny && k != 0 && k != nz)
        {
            int idx[8];
            idx[6] = cur;         // i6
            idx[7] = idx[6] - 1;  // i7
            idx[4] = idx[7] - nx; // i4
            idx[5] = idx[6] - nx; // i5

            int layerOffset = nx * ny;
            idx[0] = idx[4] - layerOffset; // i0
            idx[1] = idx[5] - layerOffset; // i1
            idx[2] = idx[6] - layerOffset; // i2
            idx[3] = idx[7] - layerOffset; // i3
            double p[12][3];
            _get_centroid(&verts1(k, j, i, 2, 0), &verts1(k, j, i, 0, 0), &verts1(k, j, i, 4, 0), &verts1(k, j, i, 6, 0), p[0]);
            _get_centroid(&verts1(k, j - 1, i, 2, 0), &verts1(k, j - 1, i, 0, 0), &verts1(k, j - 1, i, 4, 0), &verts1(k, j - 1, i, 6, 0), p[1]);
            _get_centroid(&verts1(k - 1, j - 1, i, 2, 0), &verts1(k - 1, j - 1, i, 0, 0), &verts1(k - 1, j - 1, i, 4, 0), &verts1(k - 1, j - 1, i, 6, 0), p[2]);
            _get_centroid(&verts1(k - 1, j, i, 2, 0), &verts1(k - 1, j, i, 0, 0), &verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 6, 0), p[3]);
            _get_centroid(&verts1(k, j, i, 1, 0), &verts1(k, j, i, 0, 0), &verts1(k, j, i, 4, 0), &verts1(k, j, i, 5, 0), p[4]);
            _get_centroid(&verts1(k - 1, j, i, 1, 0), &verts1(k - 1, j, i, 0, 0), &verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 5, 0), p[5]);
            _get_centroid(&verts1(k - 1, j, i - 1, 1, 0), &verts1(k - 1, j, i - 1, 0, 0), &verts1(k - 1, j, i - 1, 4, 0), &verts1(k - 1, j, i - 1, 5, 0), p[6]);
            _get_centroid(&verts1(k, j, i - 1, 1, 0), &verts1(k, j, i - 1, 0, 0), &verts1(k, j, i - 1, 4, 0), &verts1(k, j, i - 1, 5, 0), p[7]);
            _get_centroid(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 1, 0), &verts1(k, j, i, 2, 0), &verts1(k, j, i, 3, 0), p[8]);
            _get_centroid(&verts1(k, j, i - 1, 0, 0), &verts1(k, j, i - 1, 1, 0), &verts1(k, j, i - 1, 2, 0), &verts1(k, j, i - 1, 3, 0), p[9]);
            _get_centroid(&verts1(k, j - 1, i - 1, 0, 0), &verts1(k, j - 1, i - 1, 1, 0), &verts1(k, j - 1, i - 1, 2, 0), &verts1(k, j - 1, i - 1, 3, 0), p[10]);
            _get_centroid(&verts1(k, j - 1, i, 0, 0), &verts1(k, j - 1, i, 1, 0), &verts1(k, j - 1, i, 2, 0), &verts1(k, j - 1, i, 3, 0), p[11]);
            double pp[6][3];
            _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 1, 0), pp[1]);
            _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j, i - 1, 0, 0), pp[3]);
            _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 2, 0), pp[2]);
            _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j - 1, i, 0, 0), pp[0]);
            _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 4, 0), pp[5]);
            _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k - 1, j, i, 0, 0), pp[4]);
            double n[12][3];
            _get_surface_normal(&verts1(k, j, i, 0, 0), pp[2], p[0], pp[5], n[0], Axis::XPOSITIVE);
            _get_surface_normal(&verts1(k, j, i, 0, 0), pp[0], p[1], pp[5], n[1], Axis::XPOSITIVE);
            _get_surface_normal(&verts1(k, j, i, 0, 0), pp[0], p[2], pp[4], n[2], Axis::XPOSITIVE);
            _get_surface_normal(&verts1(k, j, i, 0, 0), pp[4], p[3], pp[2], n[3], Axis::XPOSITIVE);
            _get_surface_normal(&verts1(k, j, i, 0, 0), pp[1], p[4], pp[5], n[4], Axis::YPOSITIVE);
            _get_surface_normal(&verts1(k, j, i, 0, 0), pp[1], p[5], pp[4], n[5], Axis::YPOSITIVE);
            _get_surface_normal(&verts1(k, j, i, 0, 0), pp[3], p[6], pp[4], n[6], Axis::YPOSITIVE);
            _get_surface_normal(&verts1(k, j, i, 0, 0), pp[3], p[7], pp[5], n[7], Axis::YPOSITIVE);
            _get_surface_normal(&verts1(k, j, i, 0, 0), pp[1], p[8], pp[2], n[8], Axis::ZPOSITIVE);
            _get_surface_normal(&verts1(k, j, i, 0, 0), pp[3], p[9], pp[2], n[9], Axis::ZPOSITIVE);
            _get_surface_normal(&verts1(k, j, i, 0, 0), pp[3], p[10], pp[0], n[10], Axis::ZPOSITIVE);
            _get_surface_normal(&verts1(k, j, i, 0, 0), pp[0], p[11], pp[1], n[11], Axis::ZPOSITIVE);
            _get_triangle_normal(&verts1(k, j, i, 0, 0), p[5], pp[1], t[0], Axis::YPOSITIVE);
            _get_triangle_normal(&verts1(k, j, i, 0, 0), p[8], pp[1], t[1], Axis::ZPOSITIVE);
            _get_triangle_normal(&verts1(k, j, i, 0, 0), pp[1], p[4], t[2], Axis::YPOSITIVE);
            _get_triangle_normal(&verts1(k, j, i, 0, 0), p[11], pp[1], t[3], Axis::ZPOSITIVE);
            ratio[0] = _get_ratio(n[5], t[0]);
            ratio[1] = _get_ratio(n[8], t[1]);
            ratio[2] = _get_ratio(n[4], t[2]);
            ratio[3] = _get_ratio(n[11], t[3]);
            Eigen::Matrix<double, 24, 24> A;
            A.setZero();
            Eigen::Matrix3d matK[8];
            for (int l = 0; l < 8; ++l)
                _get_matK(pem1[idx[l]], matK[l]);
            Eigen::RowVector3d v[12];
            for (int l = 0; l < 12; ++l)
                v[l] << n[l][0], n[l][1], n[l][2];
            A.block(0, 18, 1, 6) << _cx(k, j, i) - p[0][0], _cy(k, j, i) - p[0][1], _cz(k, j, i) - p[0][2], p[0][0] - _cx(k, j, i - 1), p[0][1] - _cy(k, j, i - 1), p[0][2] - _cz(k, j, i - 1);
            A.block(1, 18, 1, 3) = v[0] * matK[6];
            A.block(1, 21, 1, 3) = -v[0] * matK[7];
            A.block(2, 12, 1, 6) << _cx(k, j - 1, i - 1) - p[1][0], _cy(k, j - 1, i - 1) - p[1][1], _cz(k, j - 1, i - 1) - p[1][2], p[1][0] - _cx(k, j - 1, i), p[1][1] - _cy(k, j - 1, i), p[1][2] - _cz(k, j - 1, i);
            A.block(3, 12, 1, 3) = v[1] * matK[4];
            A.block(3, 15, 1, 3) = -v[1] * matK[5];
            A.block(4, 0, 1, 6) << _cx(k - 1, j - 1, i - 1) - p[2][0], _cy(k - 1, j - 1, i - 1) - p[2][1], _cz(k - 1, j - 1, i - 1) - p[2][2], p[2][0] - _cx(k - 1, j - 1, i), p[2][1] - _cy(k - 1, j - 1, i), p[2][2] - _cz(k - 1, j - 1, i);
            A.block(5, 0, 1, 3) = v[2] * matK[0];
            A.block(5, 3, 1, 3) = -v[2] * matK[1];
            A.block(6, 6, 1, 6) << _cx(k - 1, j, i) - p[3][0], _cy(k - 1, j, i) - p[3][1], _cz(k - 1, j, i) - p[3][2], p[3][0] - _cx(k - 1, j, i - 1), p[3][1] - _cy(k - 1, j, i - 1), p[3][2] - _cz(k - 1, j, i - 1);
            A.block(7, 6, 1, 3) = v[3] * matK[2];
            A.block(7, 9, 1, 3) = -v[3] * matK[3];
            A.block(8, 15, 1, 6) << _cx(k, j - 1, i) - p[4][0], _cy(k, j - 1, i) - p[4][1], _cz(k, j - 1, i) - p[4][2], p[4][0] - _cx(k, j, i), p[4][1] - _cy(k, j, i), p[4][2] - _cz(k, j, i);
            A.block(9, 15, 1, 3) = v[4] * matK[5];
            A.block(9, 18, 1, 3) = -v[4] * matK[6];
            A.block(10, 3, 1, 6) << _cx(k - 1, j - 1, i) - p[5][0], _cy(k - 1, j - 1, i) - p[5][1], _cz(k - 1, j - 1, i) - p[5][2], p[5][0] - _cx(k - 1, j, i), p[5][1] - _cy(k - 1, j, i), p[5][2] - _cz(k - 1, j, i);
            A.block(11, 3, 1, 3) = v[5] * matK[1];
            A.block(11, 6, 1, 3) = -v[5] * matK[2];
            A.block(12, 0, 1, 3) << _cx(k - 1, j - 1, i - 1) - p[6][0], _cy(k - 1, j - 1, i - 1) - p[6][1], _cz(k - 1, j - 1, i - 1) - p[6][2];
            A.block(12, 9, 1, 3) << p[6][0] - _cx(k - 1, j, i - 1), p[6][1] - _cy(k - 1, j, i - 1), p[6][2] - _cz(k - 1, j, i - 1);
            A.block(13, 0, 1, 3) = v[6] * matK[0];
            A.block(13, 9, 1, 3) = -v[6] * matK[3];
            A.block(14, 12, 1, 3) << _cx(k, j - 1, i - 1) - p[7][0], _cy(k, j - 1, i - 1) - p[7][1], _cz(k, j - 1, i - 1) - p[7][2];
            A.block(14, 21, 1, 3) << p[7][0] - _cx(k, j, i - 1), p[7][1] - _cy(k, j, i - 1), p[7][2] - _cz(k, j, i - 1);
            A.block(15, 12, 1, 3) = v[7] * matK[4];
            A.block(15, 21, 1, 3) = -v[7] * matK[7];
            A.block(16, 6, 1, 3) << _cx(k - 1, j, i) - p[8][0], _cy(k - 1, j, i) - p[8][1], _cz(k - 1, j, i) - p[8][2];
            A.block(16, 18, 1, 3) << p[8][0] - _cx(k, j, i), p[8][1] - _cy(k, j, i), p[8][2] - _cz(k, j, i);
            A.block(17, 6, 1, 3) = v[8] * matK[2];
            A.block(17, 18, 1, 3) = -v[8] * matK[6];
            A.block(18, 9, 1, 3) << _cx(k - 1, j, i - 1) - p[9][0], _cy(k - 1, j, i - 1) - p[9][1], _cz(k - 1, j, i - 1) - p[9][2];
            A.block(18, 21, 1, 3) << p[9][0] - _cx(k, j, i - 1), p[9][1] - _cy(k, j, i - 1), p[9][2] - _cz(k, j, i - 1);
            A.block(19, 9, 1, 3) = v[9] * matK[3];
            A.block(19, 21, 1, 3) = -v[9] * matK[7];
            A.block(20, 0, 1, 3) << _cx(k - 1, j - 1, i - 1) - p[10][0], _cy(k - 1, j - 1, i - 1) - p[10][1], _cz(k - 1, j - 1, i - 1) - p[10][2];
            A.block(20, 12, 1, 3) << p[10][0] - _cx(k, j - 1, i - 1), p[10][1] - _cy(k, j - 1, i - 1), p[10][2] - _cz(k, j - 1, i - 1);
            A.block(21, 0, 1, 3) = v[10] * matK[0];
            A.block(21, 12, 1, 3) = -v[10] * matK[4];
            A.block(22, 3, 1, 3) << _cx(k - 1, j - 1, i) - p[11][0], _cy(k - 1, j - 1, i) - p[11][1], _cz(k - 1, j - 1, i) - p[11][2];
            A.block(22, 15, 1, 3) << p[11][0] - _cx(k, j - 1, i), p[11][1] - _cy(k, j - 1, i), p[11][2] - _cz(k, j - 1, i);
            A.block(23, 3, 1, 3) = v[11] * matK[1];
            A.block(23, 15, 1, 3) = -v[11] * matK[5];
            A = A.inverse();
            Eigen::RowVectorXd r = (ratio[0] * v[5] + ratio[3] * v[11]) * matK[1] * A.block(3, 0, 3, 24);
            C(idx[1], idx[0]) += r[4] + r[12] + r[20];
            C(idx[1], idx[1]) += -r[4] + r[10] + r[22];
            C(idx[1], idx[2]) += r[6] - r[10] + r[16];
            C(idx[1], idx[3]) += -r[6] - r[12] + r[18];
            C(idx[1], idx[4]) += r[2] + r[14] - r[20];
            C(idx[1], idx[5]) += -r[2] + r[8] - r[22];
            C(idx[1], idx[6]) += r[0] - r[8] - r[16];
            C(idx[1], idx[7]) += -r[0] - r[14] - r[18];
            r = (-ratio[0] * v[5] + ratio[1] * v[8]) * matK[2] * A.block(6, 0, 3, 24);
            C(idx[2], idx[0]) += r[4] + r[12] + r[20];
            C(idx[2], idx[1]) += -r[4] + r[10] + r[22];
            C(idx[2], idx[2]) += r[6] - r[10] + r[16];
            C(idx[2], idx[3]) += -r[6] - r[12] + r[18];
            C(idx[2], idx[4]) += r[2] + r[14] - r[20];
            C(idx[2], idx[5]) += -r[2] + r[8] - r[22];
            C(idx[2], idx[6]) += r[0] - r[8] - r[16];
            C(idx[2], idx[7]) += -r[0] - r[14] - r[18];
            r = (ratio[2] * v[4] - ratio[3] * v[11]) * matK[5] * A.block(15, 0, 3, 24);
            C(idx[5], idx[0]) += r[4] + r[12] + r[20];
            C(idx[5], idx[1]) += -r[4] + r[10] + r[22];
            C(idx[5], idx[2]) += r[6] - r[10] + r[16];
            C(idx[5], idx[3]) += -r[6] - r[12] + r[18];
            C(idx[5], idx[4]) += r[2] + r[14] - r[20];
            C(idx[5], idx[5]) += -r[2] + r[8] - r[22];
            C(idx[5], idx[6]) += r[0] - r[8] - r[16];
            C(idx[5], idx[7]) += -r[0] - r[14] - r[18];
            r = (-ratio[2] * v[4] - ratio[1] * v[8]) * matK[6] * A.block(18, 0, 3, 24);
            C(idx[6], idx[0]) += r[4] + r[12] + r[20];
            C(idx[6], idx[1]) += -r[4] + r[10] + r[22];
            C(idx[6], idx[2]) += r[6] - r[10] + r[16];
            C(idx[6], idx[3]) += -r[6] - r[12] + r[18];
            C(idx[6], idx[4]) += r[2] + r[14] - r[20];
            C(idx[6], idx[5]) += -r[2] + r[8] - r[22];
            C(idx[6], idx[6]) += r[0] - r[8] - r[16];
            C(idx[6], idx[7]) += -r[0] - r[14] - r[18];
            // 内部角点
        }
        break;
    case Axis::XNEGATIVE:
        assert(i > 0);
        if (i == nx && k == 0 && j != 0 && j != ny)
        {
            // 两个网格中心下标
            int i7 = cur - 1, i4 = cur - nx - 1;
            // 0,1,7,9,10交接面中点下标
            double p0[3];
            double p1[3];
            double p7[3];
            double p9[3];
            double p10[3];
            _get_centroid(&verts1(k, j, i - 1, 1, 0), &verts1(k, j, i - 1, 3, 0), &verts1(k, j, i - 1, 5, 0), &verts1(k, j, i - 1, 7, 0), p0);
            _get_centroid(&verts1(k, j - 1, i - 1, 1, 0), &verts1(k, j - 1, i - 1, 3, 0), &verts1(k, j - 1, i - 1, 5, 0), &verts1(k, j - 1, i - 1, 7, 0), p1);
            _get_centroid(&verts1(k, j, i - 1, 1, 0), &verts1(k, j, i - 1, 0, 0), &verts1(k, j, i - 1, 4, 0), &verts1(k, j, i - 1, 5, 0), p7);
            _get_centroid(&verts1(k, j, i - 1, 0, 0), &verts1(k, j, i - 1, 1, 0), &verts1(k, j, i - 1, 2, 0), &verts1(k, j, i - 1, 3, 0), p9);
            _get_centroid(&verts1(k, j - 1, i - 1, 0, 0), &verts1(k, j - 1, i - 1, 1, 0), &verts1(k, j - 1, i - 1, 2, 0), &verts1(k, j - 1, i - 1, 3, 0), p10);
            // 2,3,0,5交接边点下标
            double pp2[3];
            double pp3[3];
            double pp0[3];
            double pp5[3];
            _get_midpoint(&verts1(k, j - 1, i - 1, 3, 0), &verts1(k, j - 1, i - 1, 1, 0), pp0);
            _get_midpoint(&verts1(k, j - 1, i - 1, 2, 0), &verts1(k, j - 1, i - 1, 3, 0), pp3);
            _get_midpoint(&verts1(k, j, i - 1, 1, 0), &verts1(k, j, i - 1, 3, 0), pp2);
            _get_midpoint(&verts1(k, j, i - 1, 1, 0), &verts1(k, j, i - 1, 5, 0), pp5);
            double n0[3];
            double n1[3];
            double n7[3];
            double n9[3];
            double n10[3];
            _get_surface_normal(&verts1(k, j, i - 1, 1, 0), pp2, p0, pp5, n0, Axis::XPOSITIVE);
            _get_surface_normal(&verts1(k, j - 1, i - 1, 3, 0), pp0, p1, pp5, n1, Axis::XPOSITIVE);
            _get_surface_normal(&verts1(k, j, i - 1, 1, 0), pp5, p7, pp3, n7, Axis::YPOSITIVE);
            _get_surface_normal(&verts1(k, j, i - 1, 1, 0), pp3, p9, pp2, n9, Axis::ZPOSITIVE);
            _get_surface_normal(&verts1(k, j - 1, i - 1, 3, 0), pp0, p10, pp3, n10, Axis::ZPOSITIVE);
            _get_triangle_normal(&verts1(k, j, i - 1, 1, 0), p7, pp3, t[2], Axis::YPOSITIVE);
            ratio[2] = _get_ratio(n7, t[2]);
            Eigen::Matrix<double, 6, 6> A;
            A.setZero();
            A.row(0) << 0, 0, 0, p0[0] - _cx(k, j, i - 1), p0[1] - _cy(k, j, i - 1), p0[2] - _cz(k, j, i - 1);
            A.row(1) << p1[0] - _cx(k, j - 1, i - 1), p1[1] - _cy(k, j - 1, i - 1), p1[2] - _cz(k, j - 1, i - 1), 0, 0, 0;
            A.row(2) << _cx(k, j - 1, i - 1) - p7[0], _cy(k, j - 1, i - 1) - p7[1], _cz(k, j - 1, i - 1) - p7[2], p7[0] - _cx(k, j, i - 1), p7[1] - _cy(k, j, i - 1), p7[2] - _cz(k, j, i - 1);
            Eigen::Matrix3d matK7, matK4;
            _get_matK(pem1[i7], matK7);
            _get_matK(pem1[i4], matK4);
            Eigen::RowVector3d v0, v1, v7, v9, v10;
            v0 << n0[0], n0[1], n0[2];
            v1 << n1[0], n1[1], n1[2];
            v7 << n7[0], n7[1], n7[2];
            v9 << n9[0], n9[1], n9[2];
            v10 << n10[0], n10[1], n10[2];
            A.block(3, 0, 1, 3) = v7 * matK4;
            A.block(3, 3, 1, 3) = -v7 * matK7;
            A.block(4, 3, 1, 3) = v9 * matK7;
            A.block(5, 0, 1, 3) = v10 * matK4;
            A = A.inverse();
            Eigen::RowVectorXd r = (ratio[2] * v7) * matK4 * A.topRows(3);
            B[i4] -= pb * (r[0] + r[1]);
            C(i4, i4) += -r[1] + r[2];
            C(i4, i7) += -r[0] - r[2];
            r = (-ratio[2] * v7) * matK7 * A.bottomRows(3);
            B[i7] -= pb * (r[0] + r[1]);
            C(i7, i4) += -r[1] + r[2];
            C(i7, i7) += -r[0] - r[2];
        }
        else if (i == nx && k == nz && j != 0 && j != ny)
        {
            int i3 = cur - nx * ny - 1, i0 = i3 - nx;
            // 2,3,6,9,10交接面中点下标
            double p2[3];
            double p3[3];
            double p6[3];
            double p9[3];
            double p10[3];
            _get_centroid(&verts1(k - 1, j - 1, i - 1, 1, 0), &verts1(k - 1, j - 1, i - 1, 3, 0), &verts1(k - 1, j - 1, i - 1, 5, 0), &verts1(k - 1, j - 1, i - 1, 7, 0), p2);
            _get_centroid(&verts1(k - 1, j, i - 1, 1, 0), &verts1(k - 1, j, i - 1, 3, 0), &verts1(k - 1, j, i - 1, 5, 0), &verts1(k - 1, j, i - 1, 7, 0), p3);
            _get_centroid(&verts1(k - 1, j, i - 1, 1, 0), &verts1(k - 1, j, i - 1, 0, 0), &verts1(k - 1, j, i - 1, 4, 0), &verts1(k - 1, j, i - 1, 5, 0), p6);
            _get_centroid(&verts1(k - 1, j, i - 1, 4, 0), &verts1(k - 1, j, i - 1, 5, 0), &verts1(k - 1, j, i - 1, 6, 0), &verts1(k - 1, j, i - 1, 7, 0), p9);
            _get_centroid(&verts1(k - 1, j - 1, i - 1, 4, 0), &verts1(k - 1, j - 1, i - 1, 5, 0), &verts1(k - 1, j - 1, i - 1, 6, 0), &verts1(k - 1, j - 1, i - 1, 7, 0), p10);
            // 0,2,3,4交接边点下标
            double pp0[3];
            double pp2[3];
            double pp3[3];
            double pp4[3];
            _get_midpoint(&verts1(k - 1, j - 1, i - 1, 5, 0), &verts1(k - 1, j - 1, i - 1, 7, 0), pp0);
            _get_midpoint(&verts1(k - 1, j, i - 1, 5, 0), &verts1(k - 1, j, i - 1, 7, 0), pp2);
            _get_midpoint(&verts1(k - 1, j, i - 1, 4, 0), &verts1(k - 1, j, i - 1, 5, 0), pp3);
            _get_midpoint(&verts1(k - 1, j - 1, i - 1, 7, 0), &verts1(k - 1, j - 1, i - 1, 3, 0), pp4);
            double n2[3];
            double n3[3];
            double n6[3];
            double n9[3];
            double n10[3];
            _get_surface_normal(&verts1(k - 1, j - 1, i - 1, 7, 0), pp0, p2, pp4, n2, Axis::XPOSITIVE);
            _get_surface_normal(&verts1(k - 1, j, i - 1, 5, 0), pp2, p3, pp4, n3, Axis::XPOSITIVE);
            _get_surface_normal(&verts1(k - 1, j, i - 1, 5, 0), pp4, p6, pp3, n6, Axis::YPOSITIVE);
            _get_surface_normal(&verts1(k - 1, j, i - 1, 5, 0), pp3, p9, pp2, n9, Axis::ZPOSITIVE);
            _get_surface_normal(&verts1(k - 1, j - 1, i - 1, 7, 0), pp0, p10, pp3, n10, Axis::ZPOSITIVE);
            _get_triangle_normal(&verts1(k - 1, j, i - 1, 5, 0), p6, pp3, t[0], Axis::YPOSITIVE);
            ratio[0] = _get_ratio(n6, t[0]);
            Eigen::Matrix<double, 6, 6> A;
            A.setZero();
            A.row(0) << 0, 0, 0, p3[0] - _cx(k - 1, j, i - 1), p3[1] - _cy(k - 1, j, i - 1), p3[2] - _cz(k - 1, j, i - 1);
            A.row(1) << p2[0] - _cx(k - 1, j - 1, i - 1), p2[1] - _cy(k - 1, j - 1, i - 1), p2[2] - _cz(k - 1, j - 1, i - 1), 0, 0, 0;
            A.row(2) << _cx(k - 1, j - 1, i - 1) - p6[0], _cy(k - 1, j - 1, i - 1) - p6[1], _cz(k - 1, j - 1, i - 1) - p6[2], p6[0] - _cx(k - 1, j, i - 1), p6[1] - _cy(k - 1, j, i - 1), p6[2] - _cz(k - 1, j, i - 1);
            Eigen::Matrix3d matK3, matK0;
            _get_matK(pem1[i3], matK3);
            _get_matK(pem1[i0], matK0);
            Eigen::RowVector3d v2, v3, v6, v9, v10;
            v2 << n2[0], n2[1], n2[2];
            v3 << n3[0], n3[1], n3[2];
            v6 << n6[0], n6[1], n6[2];
            v9 << n9[0], n9[1], n9[2];
            v10 << n10[0], n10[1], n10[2];
            A.block(3, 0, 1, 3) = v6 * matK0;
            A.block(3, 3, 1, 3) = -v6 * matK3;
            A.block(4, 3, 1, 3) = v9 * matK3;
            A.block(5, 0, 1, 3) = v10 * matK0;
            A = A.inverse();
            Eigen::RowVectorXd r = (ratio[0] * v6) * matK0 * A.topRows(3);
            B[i0] -= pb * (r[0] + r[1]);
            C(i0, i0) += -r[1] + r[2];
            C(i0, i3) += -r[0] - r[2];
            r = (-ratio[0] * v6) * matK3 * A.bottomRows(3);
            B[i3] -= pb * (r[0] + r[1]);
            C(i3, i0) += -r[1] + r[2];
            C(i3, i3) += -r[0] - r[2];
        }
        else if (i == nx && j == 0 && k != 0 && k != nz)
        {
            int i7 = cur - 1, i3 = i7 - nx * ny;
            // 0,3,6,7,9交接面中点下标
            double p0[3];
            double p3[3];
            double p6[3];
            double p7[3];
            double p9[3];
            _get_centroid(&verts1(k, j, i - 1, 1, 0), &verts1(k, j, i - 1, 3, 0), &verts1(k, j, i - 1, 5, 0), &verts1(k, j, i - 1, 7, 0), p0);
            _get_centroid(&verts1(k - 1, j, i - 1, 1, 0), &verts1(k - 1, j, i - 1, 3, 0), &verts1(k - 1, j, i - 1, 5, 0), &verts1(k - 1, j, i - 1, 7, 0), p3);
            _get_centroid(&verts1(k, j, i - 1, 1, 0), &verts1(k, j, i - 1, 0, 0), &verts1(k, j, i - 1, 4, 0), &verts1(k, j, i - 1, 5, 0), p7);
            _get_centroid(&verts1(k - 1, j, i - 1, 1, 0), &verts1(k - 1, j, i - 1, 0, 0), &verts1(k - 1, j, i - 1, 4, 0), &verts1(k - 1, j, i - 1, 5, 0), p6);
            _get_centroid(&verts1(k - 1, j, i - 1, 4, 0), &verts1(k - 1, j, i - 1, 5, 0), &verts1(k - 1, j, i - 1, 6, 0), &verts1(k - 1, j, i - 1, 7, 0), p9);
            // 2,3,4,5交接边点下标
            double pp2[3];
            double pp3[3];
            double pp4[3];
            double pp5[3];
            _get_midpoint(&verts1(k, j, i - 1, 1, 0), &verts1(k, j, i - 1, 3, 0), pp2);
            _get_midpoint(&verts1(k, j, i - 1, 1, 0), &verts1(k, j, i - 1, 0, 0), pp3);
            _get_midpoint(&verts1(k, j, i - 1, 1, 0), &verts1(k, j, i - 1, 5, 0), pp5);
            _get_midpoint(&verts1(k - 1, j, i - 1, 1, 0), &verts1(k - 1, j, i - 1, 5, 0), pp4);
            double n0[3];
            double n3[3];
            double n6[3];
            double n7[3];
            double n9[3];
            _get_surface_normal(&verts1(k, j, i - 1, 1, 0), pp2, p0, pp5, n0, Axis::XPOSITIVE);
            _get_surface_normal(&verts1(k - 1, j, i - 1, 5, 0), pp4, p3, pp2, n3, Axis::XPOSITIVE);
            _get_surface_normal(&verts1(k, j, i - 1, 1, 0), pp3, p7, pp5, n7, Axis::YPOSITIVE);
            _get_surface_normal(&verts1(k - 1, j, i - 1, 5, 0), pp4, p6, pp3, n6, Axis::YPOSITIVE);
            _get_surface_normal(&verts1(k - 1, j, i - 1, 5, 0), pp2, p9, pp3, n9, Axis::ZPOSITIVE);
            _get_triangle_normal(&verts1(k - 1, j, i - 1, 5, 0), p9, pp3, t[1], Axis::ZPOSITIVE);
            ratio[1] = _get_ratio(n9, t[1]);
            Eigen::Matrix<double, 6, 6> A;
            A.setZero();
            A.row(0) << 0, 0, 0, p0[0] - _cx(k, j, i - 1), p0[1] - _cy(k, j, i - 1), p0[2] - _cz(k, j, i - 1);
            A.row(1) << p3[0] - _cx(k - 1, j, i - 1), p3[1] - _cy(k - 1, j, i - 1), p3[2] - _cz(k - 1, j, i - 1), 0, 0, 0;
            A.row(2) << _cx(k - 1, j, i - 1) - p9[0], _cy(k - 1, j, i - 1) - p9[1], _cz(k - 1, j, i - 1) - p9[2], p9[0] - _cx(k, j, i - 1), p9[1] - _cy(k, j, i - 1), p9[2] - _cz(k, j, i - 1);
            Eigen::Matrix3d matK7, matK3;
            _get_matK(pem1[i7], matK7);
            _get_matK(pem1[i3], matK3);
            Eigen::RowVector3d v0, v3, v6, v7, v9;
            v0 << n0[0], n0[1], n0[2];
            v3 << n3[0], n3[1], n3[2];
            v6 << n6[0], n6[1], n6[2];
            v7 << n7[0], n7[1], n7[2];
            v9 << n9[0], n9[1], n9[2];
            A.block(3, 0, 1, 3) = v9 * matK3;
            A.block(3, 3, 1, 3) = -v9 * matK7;
            A.block(4, 3, 1, 3) = v7 * matK7;
            A.block(5, 0, 1, 3) = v6 * matK3;
            A = A.inverse();
            Eigen::RowVectorXd r = (ratio[1] * v9) * matK3 * A.topRows(3);
            B[i3] -= pb * (r[0] + r[1]);
            C(i3, i3) += -r[1] + r[2];
            C(i3, i7) += -r[0] - r[2];
            r = (-ratio[1] * v9) * matK7 * A.bottomRows(3);
            B[i7] -= pb * (r[0] + r[1]);
            C(i7, i3) += -r[1] + r[2];
            C(i7, i7) += -r[0] - r[2];
            //  To be implemented
        }
        else if (i == nx && j == ny && k != 0 && k != nz)
        {
            int i4 = cur - 1 - nx, i0 = i4 - nx * ny;
            // 1,2,6,7,10交接面中点下标
            double p1[3];
            double p2[3];
            double p6[3];
            double p7[3];
            double p10[3];
            _get_centroid(&verts1(k, j - 1, i - 1, 1, 0), &verts1(k, j - 1, i - 1, 3, 0), &verts1(k, j - 1, i - 1, 5, 0), &verts1(k, j - 1, i - 1, 7, 0), p1);
            _get_centroid(&verts1(k - 1, j - 1, i - 1, 1, 0), &verts1(k - 1, j - 1, i - 1, 3, 0), &verts1(k - 1, j - 1, i - 1, 5, 0), &verts1(k - 1, j - 1, i - 1, 7, 0), p2);
            _get_centroid(&verts1(k, j - 1, i - 1, 3, 0), &verts1(k, j - 1, i - 1, 2, 0), &verts1(k, j - 1, i - 1, 6, 0), &verts1(k, j - 1, i - 1, 7, 0), p7);
            _get_centroid(&verts1(k - 1, j - 1, i - 1, 3, 0), &verts1(k - 1, j - 1, i - 1, 2, 0), &verts1(k - 1, j - 1, i - 1, 6, 0), &verts1(k - 1, j - 1, i - 1, 7, 0), p6);
            _get_centroid(&verts1(k - 1, j - 1, i - 1, 5, 0), &verts1(k - 1, j - 1, i - 1, 4, 0), &verts1(k - 1, j - 1, i - 1, 6, 0), &verts1(k - 1, j - 1, i - 1, 7, 0), p10);
            // 0,3,4,5交接边点下标
            double pp0[3];
            double pp3[3];
            double pp4[3];
            double pp5[3];
            _get_midpoint(&verts1(k, j - 1, i - 1, 1, 0), &verts1(k, j - 1, i - 1, 3, 0), pp0);
            _get_midpoint(&verts1(k, j - 1, i - 1, 2, 0), &verts1(k, j - 1, i - 1, 3, 0), pp3);
            _get_midpoint(&verts1(k, j - 1, i - 1, 7, 0), &verts1(k, j - 1, i - 1, 3, 0), pp5);
            _get_midpoint(&verts1(k - 1, j - 1, i - 1, 3, 0), &verts1(k - 1, j - 1, i - 1, 7, 0), pp4);
            double n1[3];
            double n2[3];
            double n6[3];
            double n7[3];
            double n10[3];
            _get_surface_normal(&verts1(k, j - 1, i - 1, 3, 0), pp0, p1, pp5, n1, Axis::XPOSITIVE);
            _get_surface_normal(&verts1(k - 1, j - 1, i - 1, 7, 0), pp4, p2, pp0, n2, Axis::XPOSITIVE);
            _get_surface_normal(&verts1(k, j - 1, i - 1, 3, 0), pp3, p7, pp5, n7, Axis::YPOSITIVE);
            _get_surface_normal(&verts1(k - 1, j - 1, i - 1, 7, 0), pp4, p6, pp3, n6, Axis::YPOSITIVE);
            _get_surface_normal(&verts1(k - 1, j - 1, i - 1, 7, 0), pp0, p10, pp3, n10, Axis::ZPOSITIVE);
            _get_triangle_normal(&verts1(k - 1, j - 1, i - 1, 7, 0), p10, pp3, t[3], Axis::ZPOSITIVE);
            ratio[3] = _get_ratio(n10, t[3]);
            Eigen::Matrix<double, 6, 6> A;
            A.setZero();
            A.row(0) << 0, 0, 0, p1[0] - _cx(k, j - 1, i - 1), p1[1] - _cy(k, j - 1, i - 1), p1[2] - _cz(k, j - 1, i - 1);
            A.row(1) << p2[0] - _cx(k - 1, j - 1, i - 1), p2[1] - _cy(k - 1, j - 1, i - 1), p2[2] - _cz(k - 1, j - 1, i - 1), 0, 0, 0;
            A.row(2) << _cx(k - 1, j - 1, i - 1) - p10[0], _cy(k - 1, j - 1, i - 1) - p10[1], _cz(k - 1, j - 1, i - 1) - p10[2], p10[0] - _cx(k, j - 1, i - 1), p10[1] - _cy(k, j - 1, i - 1), p10[2] - _cz(k, j - 1, i - 1);
            Eigen::Matrix3d matK4, matK0;
            _get_matK(pem1[i4], matK4);
            _get_matK(pem1[i0], matK0);
            Eigen::RowVector3d v1, v2, v6, v7, v10;
            v1 << n1[0], n1[1], n1[2];
            v2 << n2[0], n2[1], n2[2];
            v6 << n6[0], n6[1], n6[2];
            v7 << n7[0], n7[1], n7[2];
            v10 << n10[0], n10[1], n10[2];
            A.block(3, 0, 1, 3) = v10 * matK0;
            A.block(3, 3, 1, 3) = -v10 * matK4;
            A.block(4, 3, 1, 3) = v7 * matK4;
            A.block(5, 0, 1, 3) = v6 * matK0;
            A = A.inverse();
            Eigen::RowVectorXd r = (ratio[3] * v10) * matK0 * A.topRows(3);
            B[i0] -= pb * (r[0] + r[1]);
            C(i0, i0) += -r[1] + r[2];
            C(i0, i4) += -r[0] - r[2];
            r = (-ratio[3] * v10) * matK4 * A.bottomRows(3);
            B[i4] -= pb * (r[0] + r[1]);
            C(i4, i0) += -r[1] + r[2];
            C(i4, i4) += -r[0] - r[2];
            // To be implemented
        }
        // 六个表面其余点
        else if (i == nx && j != 0 && j != ny && k != 0 && k != nz)
        {
            int i7 = cur - 1, i4 = i7 - nx, i3 = i7 - nx * ny, i0 = i4 - nx * ny;
            // 0,1,2,3,6,7,9,10交接面中点下标
            double p0[3];
            double p1[3];
            double p2[3];
            double p3[3];
            double p6[3];
            double p7[3];
            double p9[3];
            double p10[3];
            _get_centroid(&verts1(k, j, i - 1, 1, 0), &verts1(k, j, i - 1, 3, 0), &verts1(k, j, i - 1, 5, 0), &verts1(k, j, i - 1, 7, 0), p0);
            _get_centroid(&verts1(k, j - 1, i - 1, 1, 0), &verts1(k, j - 1, i - 1, 3, 0), &verts1(k, j - 1, i - 1, 5, 0), &verts1(k, j - 1, i - 1, 7, 0), p1);
            _get_centroid(&verts1(k - 1, j - 1, i - 1, 1, 0), &verts1(k - 1, j - 1, i - 1, 3, 0), &verts1(k - 1, j - 1, i - 1, 5, 0), &verts1(k - 1, j - 1, i - 1, 7, 0), p2);
            _get_centroid(&verts1(k - 1, j, i - 1, 1, 0), &verts1(k - 1, j, i - 1, 3, 0), &verts1(k - 1, j, i - 1, 5, 0), &verts1(k - 1, j, i - 1, 7, 0), p3);
            _get_centroid(&verts1(k - 1, j, i - 1, 1, 0), &verts1(k - 1, j, i - 1, 0, 0), &verts1(k - 1, j, i - 1, 5, 0), &verts1(k - 1, j, i - 1, 4, 0), p6);
            _get_centroid(&verts1(k, j, i - 1, 1, 0), &verts1(k, j, i - 1, 0, 0), &verts1(k, j, i - 1, 5, 0), &verts1(k, j, i - 1, 4, 0), p7);
            _get_centroid(&verts1(k, j, i - 1, 1, 0), &verts1(k, j, i - 1, 3, 0), &verts1(k, j, i - 1, 2, 0), &verts1(k, j, i - 1, 0, 0), p9);
            _get_centroid(&verts1(k, j - 1, i - 1, 1, 0), &verts1(k, j - 1, i - 1, 3, 0), &verts1(k, j - 1, i - 1, 2, 0), &verts1(k, j - 1, i - 1, 0, 0), p10);
            // 0,2,3,4,5交接边点下标
            double pp0[3];
            double pp2[3];
            double pp3[3];
            double pp4[3];
            double pp5[3];
            _get_midpoint(&verts1(k, j, i - 1, 1, 0), &verts1(k, j, i - 1, 3, 0), pp2);
            _get_midpoint(&verts1(k, j, i - 1, 1, 0), &verts1(k, j, i - 1, 0, 0), pp3);
            _get_midpoint(&verts1(k, j, i - 1, 1, 0), &verts1(k, j, i - 1, 5, 0), pp5);
            _get_midpoint(&verts1(k, j - 1, i - 1, 1, 0), &verts1(k, j - 1, i - 1, 3, 0), pp0);
            _get_midpoint(&verts1(k - 1, j, i - 1, 1, 0), &verts1(k - 1, j, i - 1, 5, 0), pp4);
            double n0[3];
            double n1[3];
            double n2[3];
            double n3[3];
            double n6[3];
            double n7[3];
            double n9[3];
            double n10[3];
            _get_surface_normal(&verts1(k, j, i - 1, 1, 0), pp2, p0, pp5, n0, Axis::XPOSITIVE);
            _get_surface_normal(&verts1(k, j, i - 1, 1, 0), pp0, p1, pp5, n1, Axis::XPOSITIVE);
            _get_surface_normal(&verts1(k, j, i - 1, 1, 0), pp0, p2, pp4, n2, Axis::XPOSITIVE);
            _get_surface_normal(&verts1(k, j, i - 1, 1, 0), pp4, p3, pp2, n3, Axis::XPOSITIVE);
            _get_surface_normal(&verts1(k, j, i - 1, 1, 0), pp3, p6, pp4, n6, Axis::YPOSITIVE);
            _get_surface_normal(&verts1(k, j, i - 1, 1, 0), pp5, p7, pp3, n7, Axis::YPOSITIVE);
            _get_surface_normal(&verts1(k, j, i - 1, 1, 0), pp2, p9, pp3, n9, Axis::ZPOSITIVE);
            _get_surface_normal(&verts1(k, j, i - 1, 1, 0), pp0, p10, pp3, n10, Axis::ZPOSITIVE);
            _get_triangle_normal(&verts1(k, j, i - 1, 1, 0), pp3, p6, t[0], Axis::YPOSITIVE);
            _get_triangle_normal(&verts1(k, j, i - 1, 1, 0), p7, pp3, t[2], Axis::YPOSITIVE);
            _get_triangle_normal(&verts1(k, j, i - 1, 1, 0), p9, pp3, t[1], Axis::ZPOSITIVE);
            _get_triangle_normal(&verts1(k, j, i - 1, 1, 0), p10, pp3, t[3], Axis::ZPOSITIVE);
            ratio[0] = _get_ratio(n6, t[0]);
            ratio[1] = _get_ratio(n9, t[1]);
            ratio[2] = _get_ratio(n7, t[2]);
            ratio[3] = _get_ratio(n10, t[3]);
            Eigen::Matrix<double, 12, 12> A;
            A.setZero();
            A.block(0, 9, 1, 3) << p0[0] - _cx(k, j, i - 1), p0[1] - _cy(k, j, i - 1), p0[2] - _cz(k, j, i - 1);
            A.block(1, 6, 1, 3) << p1[0] - _cx(k, j - 1, i - 1), p1[1] - _cy(k, j - 1, i - 1), p1[2] - _cz(k, j - 1, i - 1);
            A.block(2, 0, 1, 3) << p2[0] - _cx(k - 1, j - 1, i - 1), p2[1] - _cy(k - 1, j - 1, i - 1), p2[2] - _cz(k - 1, j - 1, i - 1);
            A.block(3, 3, 1, 3) << p3[0] - _cx(k - 1, j, i - 1), p3[1] - _cy(k - 1, j, i - 1), p3[2] - _cz(k - 1, j, i - 1);
            Eigen::Matrix3d matK7, matK4, matK3, matK0;
            _get_matK(pem1[i7], matK7);
            _get_matK(pem1[i4], matK4);
            _get_matK(pem1[i3], matK3);
            _get_matK(pem1[i0], matK0);
            Eigen::RowVector3d v0, v1, v2, v3, v6, v7, v9, v10;
            v0 << n0[0], n0[1], n0[2];
            v1 << n1[0], n1[1], n1[2];
            v2 << n2[0], n2[1], n2[2];
            v3 << n3[0], n3[1], n3[2];
            v6 << n6[0], n6[1], n6[2];
            v7 << n7[0], n7[1], n7[2];
            v9 << n9[0], n9[1], n9[2];
            v10 << n10[0], n10[1], n10[2];
            A.block(4, 0, 1, 6) << _cx(k - 1, j - 1, i - 1) - p6[0], _cy(k - 1, j - 1, i - 1) - p6[1], _cz(k - 1, j - 1, i - 1) - p6[2], p6[0] - _cx(k - 1, j, i - 1), p6[1] - _cy(k - 1, j, i - 1), p6[2] - _cz(k - 1, j, i - 1);
            A.block(5, 0, 1, 3) = v6 * matK0;
            A.block(5, 3, 1, 3) = -v6 * matK3;
            A.block(6, 6, 1, 6) << _cx(k, j - 1, i - 1) - p7[0], _cy(k, j - 1, i - 1) - p7[1], _cz(k, j - 1, i - 1) - p7[2], p7[0] - _cx(k, j, i - 1), p7[1] - _cy(k, j, i - 1), p7[2] - _cz(k, j, i - 1);
            A.block(7, 6, 1, 3) = v7 * matK4;
            A.block(7, 9, 1, 3) = -v7 * matK7;
            A.block(8, 3, 1, 3) << _cx(k - 1, j, i - 1) - p9[0], _cy(k - 1, j, i - 1) - p9[1], _cz(k - 1, j, i - 1) - p9[2];
            A.block(8, 9, 1, 3) << p9[0] - _cx(k, j, i - 1), p9[1] - _cy(k, j, i - 1), p9[2] - _cz(k, j, i - 1);
            A.block(9, 3, 1, 3) = v9 * matK3;
            A.block(9, 9, 1, 3) = -v9 * matK7;
            A.block(10, 0, 1, 3) << _cx(k - 1, j - 1, i - 1) - p10[0], _cy(k - 1, j - 1, i - 1) - p10[1], _cz(k - 1, j - 1, i - 1) - p10[2];
            A.block(10, 6, 1, 3) << p10[0] - _cx(k, j - 1, i - 1), p10[1] - _cy(k, j - 1, i - 1), p10[2] - _cz(k, j - 1, i - 1);
            A.block(11, 0, 1, 3) = v10 * matK0;
            A.block(11, 6, 1, 3) = -v10 * matK4;
            A = A.inverse();
            Eigen::RowVectorXd r = (ratio[0] * v6 + ratio[3] * v10) * matK0 * A.topRows(3);
            B[i0] -= pb * (r[0] + r[1] + r[2] + r[3]);
            C(i0, i0) += -r[2] + r[4] + r[10];
            C(i0, i3) += -r[3] - r[4] + r[8];
            C(i0, i4) += -r[1] + r[6] - r[10];
            C(i0, i7) += -r[0] - r[6] - r[8];
            r = (-ratio[0] * v6 + ratio[1] * v9) * matK3 * A.block(3, 0, 3, 12);
            B[i3] -= pb * (r[0] + r[1] + r[2] + r[3]);
            C(i3, i0) += -r[2] + r[4] + r[10];
            C(i3, i3) += -r[3] - r[4] + r[8];
            C(i3, i4) += -r[1] + r[6] - r[10];
            C(i3, i7) += -r[0] - r[6] - r[8];
            r = (-ratio[3] * v10 + ratio[2] * v7) * matK4 * A.block(6, 0, 3, 12);
            B[i4] -= pb * (r[0] + r[1] + r[2] + r[3]);
            C(i4, i0) += -r[2] + r[4] + r[10];
            C(i4, i3) += -r[3] - r[4] + r[8];
            C(i4, i4) += -r[1] + r[6] - r[10];
            C(i4, i7) += -r[0] - r[6] - r[8];
            r = (-ratio[2] * v7 - ratio[1] * v9) * matK7 * A.bottomRows(3);
            B[i7] -= pb * (r[0] + r[1] + r[2] + r[3]);
            C(i7, i0) += -r[2] + r[4] + r[10];
            C(i7, i3) += -r[3] - r[4] + r[8];
            C(i7, i4) += -r[1] + r[6] - r[10];
            C(i7, i7) += -r[0] - r[6] - r[8];

            // To be implemented
        }
        else if (j == 0 && i != nx && k != 0 && k != nz)
        {
            int i6 = cur, i7 = i6 - 1, i2 = i6 - nx * ny, i3 = i7 - nx * ny;
            // 0,3,4,5,6,7,8,9交接面中点下标
            double p0[3];
            double p3[3];
            double p4[3];
            double p5[3];
            double p6[3];
            double p7[3];
            double p8[3];
            double p9[3];
            _get_centroid(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 2, 0), &verts1(k, j, i, 4, 0), &verts1(k, j, i, 6, 0), p0);
            _get_centroid(&verts1(k - 1, j, i, 0, 0), &verts1(k - 1, j, i, 2, 0), &verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 6, 0), p3);
            _get_centroid(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 1, 0), &verts1(k, j, i, 2, 0), &verts1(k, j, i, 3, 0), p8);
            _get_centroid(&verts1(k, j, i - 1, 0, 0), &verts1(k, j, i - 1, 1, 0), &verts1(k, j, i - 1, 2, 0), &verts1(k, j, i - 1, 3, 0), p9);
            _get_centroid(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 1, 0), &verts1(k, j, i, 4, 0), &verts1(k, j, i, 5, 0), p4);
            _get_centroid(&verts1(k - 1, j, i, 0, 0), &verts1(k - 1, j, i, 1, 0), &verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 5, 0), p5);
            _get_centroid(&verts1(k - 1, j, i - 1, 0, 0), &verts1(k - 1, j, i - 1, 1, 0), &verts1(k - 1, j, i - 1, 4, 0), &verts1(k - 1, j, i - 1, 5, 0), p6);
            _get_centroid(&verts1(k, j, i - 1, 0, 0), &verts1(k, j, i - 1, 1, 0), &verts1(k, j, i - 1, 4, 0), &verts1(k, j, i - 1, 5, 0), p7);
            // 1,2,3,4,5交接边点下
            double pp1[3];
            double pp2[3];
            double pp3[3];
            double pp4[3];
            double pp5[3];
            _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 2, 0), pp2);
            _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 1, 0), pp1);
            _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 4, 0), pp5);
            _get_midpoint(&verts1(k - 1, j, i, 0, 0), &verts1(k - 1, j, i, 4, 0), pp4);
            _get_midpoint(&verts1(k, j, i - 1, 0, 0), &verts1(k, j, i - 1, 1, 0), pp3);
            double n0[3];
            double n3[3];
            double n4[3];
            double n5[3];
            double n6[3];
            double n7[3];
            double n8[3];
            double n9[3];
            _get_surface_normal(&verts1(k, j, i, 0, 0), pp2, p0, pp5, n0, Axis::XPOSITIVE);
            _get_surface_normal(&verts1(k, j, i, 0, 0), pp2, p3, pp4, n3, Axis::XPOSITIVE);
            _get_surface_normal(&verts1(k, j, i, 0, 0), pp1, p4, pp5, n4, Axis::YPOSITIVE);
            _get_surface_normal(&verts1(k, j, i, 0, 0), pp4, p5, pp1, n5, Axis::YPOSITIVE);
            _get_surface_normal(&verts1(k, j, i, 0, 0), pp4, p6, pp3, n6, Axis::YPOSITIVE);
            _get_surface_normal(&verts1(k, j, i, 0, 0), pp5, p7, pp3, n7, Axis::YPOSITIVE);
            _get_surface_normal(&verts1(k, j, i, 0, 0), pp2, p8, pp1, n8, Axis::ZPOSITIVE);
            _get_surface_normal(&verts1(k, j, i, 0, 0), pp2, p9, pp3, n9, Axis::ZPOSITIVE);
            _get_triangle_normal(&verts1(k, j, i, 0, 0), p9, pp3, t[1], Axis::ZPOSITIVE);
            ratio[1] = _get_ratio(n9, t[1]);
            Eigen::Matrix<double, 12, 12> A;
            A.setZero();
            Eigen::Matrix3d matK6, matK7, matK2, matK3;
            _get_matK(pem1[i6], matK6);
            _get_matK(pem1[i7], matK7);
            _get_matK(pem1[i2], matK2);
            _get_matK(pem1[i3], matK3);
            Eigen::RowVector3d v0, v3, v4, v5, v6, v7, v8, v9;
            v0 << n0[0], n0[1], n0[2];
            v3 << n3[0], n3[1], n3[2];
            v4 << n4[0], n4[1], n4[2];
            v5 << n5[0], n5[1], n5[2];
            v6 << n6[0], n6[1], n6[2];
            v7 << n7[0], n7[1], n7[2];
            v8 << n8[0], n8[1], n8[2];
            v9 << n9[0], n9[1], n9[2];
            A.block(0, 6, 1, 6) << _cx(k, j, i) - p0[0], _cy(k, j, i) - p0[1], _cz(k, j, i) - p0[2], p0[0] - _cx(k, j, i - 1), p0[1] - _cy(k, j, i - 1), p0[2] - _cz(k, j, i - 1);
            A.block(1, 6, 1, 3) = v0 * matK6;
            A.block(1, 9, 1, 3) = -v0 * matK7;
            A.block(2, 0, 1, 6) << _cx(k - 1, j, i) - p3[0], _cy(k - 1, j, i) - p3[1], _cz(k - 1, j, i) - p3[2], p3[0] - _cx(k - 1, j, i - 1), p3[1] - _cy(k - 1, j, i - 1), p3[2] - _cz(k - 1, j, i - 1);
            A.block(3, 0, 1, 3) = v3 * matK2;
            A.block(3, 3, 1, 3) = -v3 * matK3;
            A.block(4, 6, 1, 3) = v4 * matK6;
            A.block(5, 0, 1, 3) = v5 * matK2;
            A.block(6, 3, 1, 3) = v6 * matK3;
            A.block(7, 9, 1, 3) = v7 * matK7;
            A.block(8, 0, 1, 3) << _cx(k - 1, j, i) - p8[0], _cy(k - 1, j, i) - p8[1], _cz(k - 1, j, i) - p8[2];
            A.block(8, 6, 1, 3) << p8[0] - _cx(k, j, i), p8[1] - _cy(k, j, i), p8[2] - _cz(k, j, i);
            A.block(9, 0, 1, 3) = v8 * matK2;
            A.block(9, 6, 1, 3) = -v8 * matK6;
            A.block(10, 3, 1, 3) << _cx(k - 1, j, i - 1) - p9[0], _cy(k - 1, j, i - 1) - p9[1], _cz(k - 1, j, i - 1) - p9[2];
            A.block(10, 9, 1, 3) << p9[0] - _cx(k, j, i - 1), p9[1] - _cy(k, j, i - 1), p9[2] - _cz(k, j, i - 1);
            A.block(11, 3, 1, 3) = v9 * matK3;
            A.block(11, 9, 1, 3) = -v9 * matK7;
            A = A.inverse();
            Eigen::RowVectorXd r = (ratio[1] * v9) * matK3 * A.block(3, 0, 3, 12);
            C(i3, i2) += r[2] + r[8];
            C(i3, i3) += -r[2] + r[10];
            C(i3, i6) += r[0] - r[8];
            C(i3, i7) += -r[0] - r[10];
            r = (-ratio[1] * v9) * matK7 * A.bottomRows(3);
            C(i7, i2) += r[2] + r[8];
            C(i7, i3) += -r[2] + r[10];
            C(i7, i6) += r[0] - r[8];
            C(i7, i7) += -r[0] - r[10];
            // To be implemented
        }
        else if (j == ny && i != nx && k != 0 && k != nz)
        {
            int i5 = cur - nx, i4 = i5 - 1, i1 = i5 - nx * ny, i0 = i4 - nx * ny;
            // 1,2,4,5,6,7,10,11交接面中点下标
            double p1[3];
            double p2[3];
            double p4[3];
            double p5[3];
            double p6[3];
            double p7[3];
            double p10[3];
            double p11[3];
            _get_centroid(&verts1(k, j - 1, i, 2, 0), &verts1(k, j - 1, i, 0, 0), &verts1(k, j - 1, i, 4, 0), &verts1(k, j - 1, i, 6, 0), p1);
            _get_centroid(&verts1(k - 1, j - 1, i, 2, 0), &verts1(k - 1, j - 1, i, 0, 0), &verts1(k - 1, j - 1, i, 4, 0), &verts1(k - 1, j - 1, i, 6, 0), p2);
            _get_centroid(&verts1(k, j - 1, i, 2, 0), &verts1(k, j - 1, i, 3, 0), &verts1(k, j - 1, i, 6, 0), &verts1(k, j - 1, i, 7, 0), p4);
            _get_centroid(&verts1(k, j - 1, i - 1, 2, 0), &verts1(k, j - 1, i - 1, 3, 0), &verts1(k, j - 1, i - 1, 6, 0), &verts1(k, j - 1, i - 1, 7, 0), p7);
            _get_centroid(&verts1(k - 1, j - 1, i, 2, 0), &verts1(k - 1, j - 1, i, 3, 0), &verts1(k - 1, j - 1, i, 6, 0), &verts1(k - 1, j - 1, i, 7, 0), p5);
            _get_centroid(&verts1(k - 1, j - 1, i - 1, 2, 0), &verts1(k - 1, j - 1, i - 1, 3, 0), &verts1(k - 1, j - 1, i - 1, 6, 0), &verts1(k - 1, j - 1, i - 1, 7, 0), p6);
            _get_centroid(&verts1(k, j - 1, i, 2, 0), &verts1(k, j - 1, i, 0, 0), &verts1(k, j - 1, i, 1, 0), &verts1(k, j - 1, i, 3, 0), p11);
            _get_centroid(&verts1(k, j - 1, i - 1, 2, 0), &verts1(k, j - 1, i - 1, 0, 0), &verts1(k, j - 1, i - 1, 1, 0), &verts1(k, j - 1, i - 1, 3, 0), p10);
            // 0,1,3,4,5交接边点下标
            double pp0[3];
            double pp1[3];
            double pp3[3];
            double pp4[3];
            double pp5[3];
            _get_midpoint(&verts1(k, j - 1, i, 2, 0), &verts1(k, j - 1, i, 0, 0), pp0);
            _get_midpoint(&verts1(k, j - 1, i, 2, 0), &verts1(k, j - 1, i, 3, 0), pp1);
            _get_midpoint(&verts1(k, j - 1, i, 2, 0), &verts1(k, j - 1, i, 6, 0), pp5);
            _get_midpoint(&verts1(k, j - 1, i - 1, 2, 0), &verts1(k, j - 1, i - 1, 3, 0), pp3);
            _get_midpoint(&verts1(k - 1, j - 1, i, 2, 0), &verts1(k - 1, j - 1, i, 6, 0), pp4);
            double n1[3];
            double n2[3];
            double n4[3];
            double n5[3];
            double n6[3];
            double n7[3];
            double n10[3];
            double n11[3];
            _get_surface_normal(&verts1(k, j - 1, i, 2, 0), pp0, p1, pp5, n1, Axis::XPOSITIVE);
            _get_surface_normal(&verts1(k, j - 1, i, 2, 0), pp0, p2, pp4, n2, Axis::XPOSITIVE);
            _get_surface_normal(&verts1(k, j - 1, i, 2, 0), pp1, p4, pp5, n4, Axis::YPOSITIVE);
            _get_surface_normal(&verts1(k, j - 1, i, 2, 0), pp4, p5, pp1, n5, Axis::YPOSITIVE);
            _get_surface_normal(&verts1(k, j - 1, i, 2, 0), pp4, p6, pp3, n6, Axis::YPOSITIVE);
            _get_surface_normal(&verts1(k, j - 1, i, 2, 0), pp5, p7, pp3, n7, Axis::YPOSITIVE);
            _get_surface_normal(&verts1(k, j - 1, i, 2, 0), pp0, p10, pp3, n10, Axis::ZPOSITIVE);
            _get_surface_normal(&verts1(k, j - 1, i, 2, 0), pp0, p11, pp1, n11, Axis::ZPOSITIVE);
            _get_triangle_normal(&verts1(k, j - 1, i, 2, 0), p10, pp3, t[3], Axis::ZPOSITIVE);
            ratio[3] = _get_ratio(n10, t[3]);
            Eigen::Matrix<double, 12, 12>
                A;
            A.setZero();
            Eigen::Matrix3d matK5, matK4, matK1, matK0;
            _get_matK(pem1[i5], matK5);
            _get_matK(pem1[i4], matK4);
            _get_matK(pem1[i1], matK1);
            _get_matK(pem1[i0], matK0);
            Eigen::RowVector3d v1, v2, v4, v5, v6, v7, v10, v11;
            v1 << n1[0], n1[1], n1[2];
            v2 << n2[0], n2[1], n2[2];
            v4 << n4[0], n4[1], n4[2];
            v5 << n5[0], n5[1], n5[2];
            v6 << n6[0], n6[1], n6[2];
            v7 << n7[0], n7[1], n7[2];
            v10 << n10[0], n10[1], n10[2];
            v11 << n11[0], n11[1], n11[2];
            A.block(0, 6, 1, 6) << _cx(k, j - 1, i - 1) - p1[0], _cy(k, j - 1, i - 1) - p1[1], _cz(k, j - 1, i - 1) - p1[2], p1[0] - _cx(k, j - 1, i), p1[1] - _cy(k, j - 1, i), p1[2] - _cz(k, j - 1, i);
            A.block(1, 6, 1, 3) = v1 * matK4;
            A.block(1, 9, 1, 3) = -v1 * matK5;
            A.block(2, 0, 1, 6) << _cx(k - 1, j - 1, i - 1) - p2[0], _cy(k - 1, j - 1, i - 1) - p2[1], _cz(k - 1, j - 1, i - 1) - p2[2], p2[0] - _cx(k - 1, j - 1, i), p2[1] - _cy(k - 1, j - 1, i), p2[2] - _cz(k - 1, j - 1, i);
            A.block(3, 0, 1, 3) = v2 * matK0;
            A.block(3, 3, 1, 3) = -v2 * matK1;
            A.block(4, 9, 1, 3) = v4 * matK5;
            A.block(5, 3, 1, 3) = v5 * matK1;
            A.block(6, 0, 1, 3) = v6 * matK0;
            A.block(7, 6, 1, 3) = v7 * matK4;
            A.block(8, 0, 1, 3) << _cx(k - 1, j - 1, i - 1) - p10[0], _cy(k - 1, j - 1, i - 1) - p10[1], _cz(k - 1, j - 1, i - 1) - p10[2];
            A.block(8, 6, 1, 3) << p10[0] - _cx(k, j - 1, i - 1), p10[1] - _cy(k, j - 1, i - 1), p10[2] - _cz(k, j - 1, i - 1);
            A.block(9, 0, 1, 3) = v10 * matK0;
            A.block(9, 6, 1, 3) = -v10 * matK4;
            A.block(10, 3, 1, 3) << _cx(k - 1, j - 1, i) - p11[0], _cy(k - 1, j - 1, i) - p11[1], _cz(k - 1, j - 1, i) - p11[2];
            A.block(10, 9, 1, 3) << p11[0] - _cx(k, j - 1, i), p11[1] - _cy(k, j - 1, i), p11[2] - _cz(k, j - 1, i);
            A.block(11, 3, 1, 3) = v11 * matK1;
            A.block(11, 9, 1, 3) = -v11 * matK5;
            A = A.inverse();
            Eigen::RowVectorXd r = (ratio[3] * v10) * matK0 * A.topRows(3);
            C(i0, i0) += r[2] + r[8];
            C(i0, i1) += -r[2] + r[10];
            C(i0, i4) += r[0] - r[8];
            C(i0, i5) += -r[0] - r[10];
            r = (-ratio[3] * v10) * matK4 * A.block(6, 0, 3, 12);
            C(i4, i0) += r[2] + r[8];
            C(i4, i1) += -r[2] + r[10];
            C(i4, i4) += r[0] - r[8];
            C(i4, i5) += -r[0] - r[10];
            // To be implemented
        }
        else if (k == 0 && i != nx && j != 0 && j != ny)
        {
            int i6 = cur, i7 = cur - 1, i5 = cur - nx, i4 = i5 - 1;
            // 0,1,4,7,8,9,10,11交接面中点下标
            double p0[3];
            double p1[3];
            double p4[3];
            double p7[3];
            double p8[3];
            double p9[3];
            double p10[3];
            double p11[3];
            _get_centroid(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 2, 0), &verts1(k, j, i, 4, 0), &verts1(k, j, i, 6, 0), p0);
            _get_centroid(&verts1(k, j - 1, i, 0, 0), &verts1(k, j - 1, i, 2, 0), &verts1(k, j - 1, i, 4, 0), &verts1(k, j - 1, i, 6, 0), p1);
            _get_centroid(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 1, 0), &verts1(k, j, i, 4, 0), &verts1(k, j, i, 5, 0), p4);
            _get_centroid(&verts1(k, j, i - 1, 0, 0), &verts1(k, j, i - 1, 1, 0), &verts1(k, j, i - 1, 4, 0), &verts1(k, j, i - 1, 5, 0), p7);
            _get_centroid(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 1, 0), &verts1(k, j, i, 2, 0), &verts1(k, j, i, 3, 0), p8);
            _get_centroid(&verts1(k, j - 1, i, 0, 0), &verts1(k, j - 1, i, 1, 0), &verts1(k, j - 1, i, 2, 0), &verts1(k, j - 1, i, 3, 0), p11);
            _get_centroid(&verts1(k, j, i - 1, 0, 0), &verts1(k, j, i - 1, 1, 0), &verts1(k, j, i - 1, 2, 0), &verts1(k, j, i - 1, 3, 0), p9);
            _get_centroid(&verts1(k, j - 1, i - 1, 0, 0), &verts1(k, j - 1, i - 1, 1, 0), &verts1(k, j - 1, i - 1, 2, 0), &verts1(k, j - 1, i - 1, 3, 0), p10);
            // 0,1,2,3,5 交接边点下标
            double pp0[3];
            double pp1[3];
            double pp2[3];
            double pp3[3];
            double pp5[3];
            _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 2, 0), pp2);
            _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 1, 0), pp1);
            _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 4, 0), pp5);
            _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j - 1, i, 0, 0), pp0);
            _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j, i - 1, 0, 0), pp3);
            double n0[3];
            double n1[3];
            double n4[3];
            double n7[3];
            double n8[3];
            double n9[3];
            double n10[3];
            double n11[3];
            _get_surface_normal(&verts1(k, j, i, 0, 0), pp2, p0, pp5, n0, Axis::XPOSITIVE);
            _get_surface_normal(&verts1(k, j, i, 0, 0), pp5, p1, pp0, n1, Axis::XPOSITIVE);
            _get_surface_normal(&verts1(k, j, i, 0, 0), pp1, p4, pp5, n4, Axis::YPOSITIVE);
            _get_surface_normal(&verts1(k, j, i, 0, 0), pp5, p7, pp3, n7, Axis::YPOSITIVE);
            _get_surface_normal(&verts1(k, j, i, 0, 0), pp2, p8, pp1, n8, Axis::ZPOSITIVE);
            _get_surface_normal(&verts1(k, j, i, 0, 0), pp0, p11, pp1, n11, Axis::ZPOSITIVE);
            _get_surface_normal(&verts1(k, j, i, 0, 0), pp0, p10, pp3, n10, Axis::ZPOSITIVE);
            _get_surface_normal(&verts1(k, j, i, 0, 0), pp3, p9, pp2, n9, Axis::ZPOSITIVE);
            _get_triangle_normal(&verts1(k, j, i, 0, 0), p7, pp3, t[2], Axis::YPOSITIVE);
            ratio[2] = _get_ratio(n7, t[2]);
            Eigen::Matrix<double, 12, 12> A;
            A.setZero();
            Eigen::Matrix3d matK6, matK7, matK5, matK4;
            _get_matK(pem1[i6], matK6);
            _get_matK(pem1[i7], matK7);
            _get_matK(pem1[i5], matK5);
            _get_matK(pem1[i4], matK4);
            Eigen::RowVector3d v0, v1, v4, v7, v8, v9, v10, v11;
            v0 << n0[0], n0[1], n0[2];
            v1 << n1[0], n1[1], n1[2];
            v4 << n4[0], n4[1], n4[2];
            v7 << n7[0], n7[1], n7[2];
            v8 << n8[0], n8[1], n8[2];
            v9 << n9[0], n9[1], n9[2];
            v10 << n10[0], n10[1], n10[2];
            v11 << n11[0], n11[1], n11[2];
            A.block(0, 6, 1, 6) << _cx(k, j, i) - p0[0], _cy(k, j, i) - p0[1], _cz(k, j, i) - p0[2], p0[0] - _cx(k, j, i - 1), p0[1] - _cy(k, j, i - 1), p0[2] - _cz(k, j, i - 1);
            A.block(1, 6, 1, 3) = v0 * matK6;
            A.block(1, 9, 1, 3) = -v0 * matK7;
            A.block(2, 0, 1, 6) << p1[0] - _cx(k, j - 1, i - 1), p1[1] - _cy(k, j - 1, i - 1), p1[2] - _cz(k, j - 1, i - 1), _cx(k, j - 1, i) - p1[0], _cy(k, j - 1, i) - p1[1], _cz(k, j - 1, i) - p1[2];
            A.block(3, 0, 1, 3) = v1 * matK4;
            A.block(3, 3, 1, 3) = -v1 * matK5;
            A.block(4, 3, 1, 6) << _cx(k, j - 1, i) - p4[0], _cy(k, j - 1, i) - p4[1], _cz(k, j - 1, i) - p4[2], p4[0] - _cx(k, j, i), p4[1] - _cy(k, j, i), p4[2] - _cz(k, j, i);
            A.block(5, 3, 1, 3) = v4 * matK5;
            A.block(5, 6, 1, 3) = -v4 * matK6;
            A.block(6, 0, 1, 3) << _cx(k, j - 1, i - 1) - p7[0], _cy(k, j - 1, i - 1) - p7[1], _cz(k, j - 1, i - 1) - p7[2];
            A.block(6, 9, 1, 3) << p7[0] - _cx(k, j, i - 1), p7[1] - _cy(k, j, i - 1), p7[2] - _cz(k, j, i - 1);
            A.block(7, 0, 1, 3) = v7 * matK4;
            A.block(7, 9, 1, 3) = -v7 * matK7;
            A.block(8, 6, 1, 3) = v8 * matK6;
            A.block(9, 9, 1, 3) = v9 * matK7;
            A.block(10, 0, 1, 3) = v10 * matK4;
            A.block(11, 3, 1, 3) = v11 * matK5;
            A = A.inverse();
            Eigen::RowVectorXd r = (ratio[2] * v7) * matK4 * A.topRows(3);
            C(i4, i4) += -r[2] + r[6];
            C(i4, i5) += r[2] + r[4];
            C(i4, i6) += r[0] - r[4];
            C(i4, i7) += -r[0] - r[6];
            r = (-ratio[2] * v7) * matK7 * A.bottomRows(3);
            C(i7, i4) += -r[2] + r[6];
            C(i7, i5) += r[2] + r[4];
            C(i7, i6) += r[0] - r[4];
            C(i7, i7) += -r[0] - r[6];
            // To be implemented
        }
        else if (k == nz && i != nx && j != 0 && j != ny)
        {
            int i2 = cur - nx * ny, i3 = i2 - 1, i1 = i2 - nx, i0 = i1 - 1;
            // 2,3,5,6,8,9,10,11交接面中点下标
            double p2[3];
            double p3[3];
            double p5[3];
            double p6[3];
            double p8[3];
            double p9[3];
            double p10[3];
            double p11[3];
            _get_centroid(&verts1(k - 1, j, i, 2, 0), &verts1(k - 1, j, i, 0, 0), &verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 6, 0), p3);
            _get_centroid(&verts1(k - 1, j - 1, i, 2, 0), &verts1(k - 1, j - 1, i, 0, 0), &verts1(k - 1, j - 1, i, 4, 0), &verts1(k - 1, j - 1, i, 6, 0), p2);
            _get_centroid(&verts1(k - 1, j, i, 0, 0), &verts1(k - 1, j, i, 1, 0), &verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 5, 0), p5);
            _get_centroid(&verts1(k - 1, j, i - 1, 0, 0), &verts1(k - 1, j, i - 1, 1, 0), &verts1(k - 1, j, i - 1, 4, 0), &verts1(k - 1, j, i - 1, 5, 0), p6);
            _get_centroid(&verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 5, 0), &verts1(k - 1, j, i, 6, 0), &verts1(k - 1, j, i, 7, 0), p8);
            _get_centroid(&verts1(k - 1, j, i - 1, 4, 0), &verts1(k - 1, j, i - 1, 5, 0), &verts1(k - 1, j, i - 1, 6, 0), &verts1(k - 1, j, i - 1, 7, 0), p9);
            _get_centroid(&verts1(k - 1, j - 1, i, 4, 0), &verts1(k - 1, j - 1, i, 5, 0), &verts1(k - 1, j - 1, i, 6, 0), &verts1(k - 1, j - 1, i, 7, 0), p11);
            _get_centroid(&verts1(k - 1, j - 1, i - 1, 4, 0), &verts1(k - 1, j - 1, i - 1, 5, 0), &verts1(k - 1, j - 1, i - 1, 6, 0), &verts1(k - 1, j - 1, i - 1, 7, 0), p10);
            // 0,1,2,3,4 交接边点下标
            double pp0[3];
            double pp1[3];
            double pp2[3];
            double pp3[3];
            double pp4[3];
            _get_midpoint(&verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 6, 0), pp2);
            _get_midpoint(&verts1(k - 1, j - 1, i, 4, 0), &verts1(k - 1, j - 1, i, 6, 0), pp0);
            _get_midpoint(&verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 5, 0), pp1);
            _get_midpoint(&verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i - 1, 4, 0), pp3);
            _get_midpoint(&verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 0, 0), pp4);
            double n2[3];
            double n3[3];
            double n5[3];
            double n6[3];
            double n8[3];
            double n9[3];
            double n10[3];
            double n11[3];
            _get_surface_normal(&verts1(k - 1, j, i, 4, 0), pp0, p2, pp4, n2, Axis::XPOSITIVE);
            _get_surface_normal(&verts1(k - 1, j, i, 4, 0), pp4, p3, pp2, n3, Axis::XPOSITIVE);
            _get_surface_normal(&verts1(k - 1, j, i, 4, 0), pp1, p5, pp4, n5, Axis::YPOSITIVE);
            _get_surface_normal(&verts1(k - 1, j, i, 4, 0), pp4, p6, pp3, n6, Axis::YPOSITIVE);
            _get_surface_normal(&verts1(k - 1, j, i, 4, 0), pp2, p8, pp1, n8, Axis::ZPOSITIVE);
            _get_surface_normal(&verts1(k - 1, j, i, 4, 0), pp3, p9, pp2, n9, Axis::ZPOSITIVE);
            _get_surface_normal(&verts1(k - 1, j, i, 4, 0), pp0, p11, pp1, n11, Axis::ZPOSITIVE);
            _get_surface_normal(&verts1(k - 1, j, i, 4, 0), pp0, p10, pp3, n10, Axis::ZPOSITIVE);
            _get_triangle_normal(&verts1(k - 1, j, i, 4, 0), p6, pp3, t[0], Axis::YPOSITIVE);
            ratio[0] = _get_ratio(n6, t[0]);
            Eigen::Matrix<double, 12, 12> A;
            A.setZero();
            Eigen::Matrix3d matK2, matK3, matK1, matK0;
            _get_matK(pem1[i2], matK2);
            _get_matK(pem1[i3], matK3);
            _get_matK(pem1[i1], matK1);
            _get_matK(pem1[i0], matK0);
            Eigen::RowVector3d v2, v3, v5, v6, v8, v9, v10, v11;
            v2 << n2[0], n2[1], n2[2];
            v3 << n3[0], n3[1], n3[2];
            v5 << n5[0], n5[1], n5[2];
            v6 << n6[0], n6[1], n6[2];
            v8 << n8[0], n8[1], n8[2];
            v9 << n9[0], n9[1], n9[2];
            v10 << n10[0], n10[1], n10[2];
            v11 << n11[0], n11[1], n11[2];
            A.block(0, 0, 1, 6) << _cx(k - 1, j - 1, i - 1) - p2[0], _cy(k - 1, j - 1, i - 1) - p2[1], _cz(k - 1, j - 1, i - 1) - p2[2], p2[0] - _cx(k - 1, j - 1, i), p2[1] - _cy(k - 1, j - 1, i), p2[2] - _cz(k - 1, j - 1, i);
            A.block(1, 0, 1, 3) = v2 * matK0;
            A.block(1, 3, 1, 3) = -v2 * matK1;
            A.block(2, 6, 1, 6) << p3[0] - _cx(k - 1, j, i), p3[1] - _cy(k - 1, j, i), p3[2] - _cz(k - 1, j, i), _cx(k - 1, j, i - 1) - p3[0], _cy(k - 1, j, i - 1) - p3[1], _cz(k - 1, j, i - 1) - p3[2];
            A.block(3, 6, 1, 3) = v3 * matK2;
            A.block(3, 9, 1, 3) = -v3 * matK3;
            A.block(4, 3, 1, 6) << _cx(k - 1, j - 1, i) - p5[0], _cy(k - 1, j - 1, i) - p5[1], _cz(k - 1, j - 1, i) - p5[2], p5[0] - _cx(k - 1, j, i), p5[1] - _cy(k - 1, j, i), p5[2] - _cz(k - 1, j, i);
            A.block(5, 3, 1, 3) = v5 * matK1;
            A.block(5, 6, 1, 3) = -v5 * matK2;
            A.block(6, 0, 1, 3) << _cx(k - 1, j - 1, i - 1) - p6[0], _cy(k - 1, j - 1, i - 1) - p6[1], _cz(k - 1, j - 1, i - 1) - p6[2];
            A.block(6, 9, 1, 3) << p6[0] - _cx(k - 1, j, i - 1), p6[1] - _cy(k - 1, j, i - 1), p6[2] - _cz(k - 1, j, i - 1);
            A.block(7, 0, 1, 3) = v6 * matK0;
            A.block(7, 9, 1, 3) = -v6 * matK3;
            A.block(8, 6, 1, 3) = v8 * matK2;
            A.block(9, 9, 1, 3) = v9 * matK3;
            A.block(10, 0, 1, 3) = v10 * matK0;
            A.block(11, 3, 1, 3) = v11 * matK1;
            A = A.inverse();
            Eigen::RowVectorXd r = (ratio[0] * v6) * matK0 * A.topRows(3);
            C(i0, i0) += r[0] + r[6];
            C(i0, i1) += -r[0] + r[4];
            C(i0, i2) += -r[2] - r[4];
            C(i0, i3) += r[2] - r[6];
            r = (-ratio[0] * v6) * matK3 * A.bottomRows(3);
            C(i3, i0) += r[0] + r[6];
            C(i3, i1) += -r[0] + r[4];
            C(i3, i2) += -r[2] - r[4];
            C(i3, i3) += r[2] - r[6];
            //  To be implemented
        }
        else if (i != nx && j != 0 && j != ny && k != 0 && k != nz)
        {
            int idx[8];
            idx[6] = cur;         // i6
            idx[7] = idx[6] - 1;  // i7
            idx[4] = idx[7] - nx; // i4
            idx[5] = idx[6] - nx; // i5

            int layerOffset = nx * ny;
            idx[0] = idx[4] - layerOffset; // i0
            idx[1] = idx[5] - layerOffset; // i1
            idx[2] = idx[6] - layerOffset; // i2
            idx[3] = idx[7] - layerOffset; // i3
            double p[12][3];
            _get_centroid(&verts1(k, j, i, 2, 0), &verts1(k, j, i, 0, 0), &verts1(k, j, i, 4, 0), &verts1(k, j, i, 6, 0), p[0]);
            _get_centroid(&verts1(k, j - 1, i, 2, 0), &verts1(k, j - 1, i, 0, 0), &verts1(k, j - 1, i, 4, 0), &verts1(k, j - 1, i, 6, 0), p[1]);
            _get_centroid(&verts1(k - 1, j - 1, i, 2, 0), &verts1(k - 1, j - 1, i, 0, 0), &verts1(k - 1, j - 1, i, 4, 0), &verts1(k - 1, j - 1, i, 6, 0), p[2]);
            _get_centroid(&verts1(k - 1, j, i, 2, 0), &verts1(k - 1, j, i, 0, 0), &verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 6, 0), p[3]);
            _get_centroid(&verts1(k, j, i, 1, 0), &verts1(k, j, i, 0, 0), &verts1(k, j, i, 4, 0), &verts1(k, j, i, 5, 0), p[4]);
            _get_centroid(&verts1(k - 1, j, i, 1, 0), &verts1(k - 1, j, i, 0, 0), &verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 5, 0), p[5]);
            _get_centroid(&verts1(k - 1, j, i - 1, 1, 0), &verts1(k - 1, j, i - 1, 0, 0), &verts1(k - 1, j, i - 1, 4, 0), &verts1(k - 1, j, i - 1, 5, 0), p[6]);
            _get_centroid(&verts1(k, j, i - 1, 1, 0), &verts1(k, j, i - 1, 0, 0), &verts1(k, j, i - 1, 4, 0), &verts1(k, j, i - 1, 5, 0), p[7]);
            _get_centroid(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 1, 0), &verts1(k, j, i, 2, 0), &verts1(k, j, i, 3, 0), p[8]);
            _get_centroid(&verts1(k, j, i - 1, 0, 0), &verts1(k, j, i - 1, 1, 0), &verts1(k, j, i - 1, 2, 0), &verts1(k, j, i - 1, 3, 0), p[9]);
            _get_centroid(&verts1(k, j - 1, i - 1, 0, 0), &verts1(k, j - 1, i - 1, 1, 0), &verts1(k, j - 1, i - 1, 2, 0), &verts1(k, j - 1, i - 1, 3, 0), p[10]);
            _get_centroid(&verts1(k, j - 1, i, 0, 0), &verts1(k, j - 1, i, 1, 0), &verts1(k, j - 1, i, 2, 0), &verts1(k, j - 1, i, 3, 0), p[11]);
            double pp[6][3];
            _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 1, 0), pp[1]);
            _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j, i - 1, 0, 0), pp[3]);
            _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 2, 0), pp[2]);
            _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j - 1, i, 0, 0), pp[0]);
            _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 4, 0), pp[5]);
            _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k - 1, j, i, 0, 0), pp[4]);
            double n[12][3];
            _get_surface_normal(&verts1(k, j, i, 0, 0), pp[2], p[0], pp[5], n[0], Axis::XPOSITIVE);
            _get_surface_normal(&verts1(k, j, i, 0, 0), pp[0], p[1], pp[5], n[1], Axis::XPOSITIVE);
            _get_surface_normal(&verts1(k, j, i, 0, 0), pp[0], p[2], pp[4], n[2], Axis::XPOSITIVE);
            _get_surface_normal(&verts1(k, j, i, 0, 0), pp[4], p[3], pp[2], n[3], Axis::XPOSITIVE);
            _get_surface_normal(&verts1(k, j, i, 0, 0), pp[1], p[4], pp[5], n[4], Axis::YPOSITIVE);
            _get_surface_normal(&verts1(k, j, i, 0, 0), pp[1], p[5], pp[4], n[5], Axis::YPOSITIVE);
            _get_surface_normal(&verts1(k, j, i, 0, 0), pp[3], p[6], pp[4], n[6], Axis::YPOSITIVE);
            _get_surface_normal(&verts1(k, j, i, 0, 0), pp[3], p[7], pp[5], n[7], Axis::YPOSITIVE);
            _get_surface_normal(&verts1(k, j, i, 0, 0), pp[1], p[8], pp[2], n[8], Axis::ZPOSITIVE);
            _get_surface_normal(&verts1(k, j, i, 0, 0), pp[3], p[9], pp[2], n[9], Axis::ZPOSITIVE);
            _get_surface_normal(&verts1(k, j, i, 0, 0), pp[3], p[10], pp[0], n[10], Axis::ZPOSITIVE);
            _get_surface_normal(&verts1(k, j, i, 0, 0), pp[0], p[11], pp[1], n[11], Axis::ZPOSITIVE);
            _get_triangle_normal(&verts1(k, j, i, 0, 0), pp[3], p[6], t[0], Axis::YPOSITIVE);
            _get_triangle_normal(&verts1(k, j, i, 0, 0), pp[3], p[7], t[2], Axis::YPOSITIVE);
            _get_triangle_normal(&verts1(k, j, i, 0, 0), pp[3], p[9], t[1], Axis::ZPOSITIVE);
            _get_triangle_normal(&verts1(k, j, i, 0, 0), pp[3], p[10], t[3], Axis::ZPOSITIVE);
            ratio[0] = _get_ratio(p[6], t[0]);
            ratio[1] = _get_ratio(p[9], t[1]);
            ratio[2] = _get_ratio(p[7], t[2]);
            ratio[3] = _get_ratio(p[10], t[3]);
            Eigen::Matrix<double, 24, 24> A;
            A.setZero();
            Eigen::Matrix3d matK[8];
            for (int l = 0; l < 8; ++l)
                _get_matK(pem1[idx[l]], matK[l]);
            Eigen::RowVector3d v[12];
            for (int l = 0; l < 12; ++l)
                v[l] << n[l][0], n[l][1], n[l][2];
            A.block(0, 18, 1, 6) << _cx(k, j, i) - p[0][0], _cy(k, j, i) - p[0][1], _cz(k, j, i) - p[0][2], p[0][0] - _cx(k, j, i - 1), p[0][1] - _cy(k, j, i - 1), p[0][2] - _cz(k, j, i - 1);
            A.block(1, 18, 1, 3) = v[0] * matK[6];
            A.block(1, 21, 1, 3) = -v[0] * matK[7];
            A.block(2, 12, 1, 6) << _cx(k, j - 1, i - 1) - p[1][0], _cy(k, j - 1, i - 1) - p[1][1], _cz(k, j - 1, i - 1) - p[1][2], p[1][0] - _cx(k, j - 1, i), p[1][1] - _cy(k, j - 1, i), p[1][2] - _cz(k, j - 1, i);
            A.block(3, 12, 1, 3) = v[1] * matK[4];
            A.block(3, 15, 1, 3) = -v[1] * matK[5];
            A.block(4, 0, 1, 6) << _cx(k - 1, j - 1, i - 1) - p[2][0], _cy(k - 1, j - 1, i - 1) - p[2][1], _cz(k - 1, j - 1, i - 1) - p[2][2], p[2][0] - _cx(k - 1, j - 1, i), p[2][1] - _cy(k - 1, j - 1, i), p[2][2] - _cz(k - 1, j - 1, i);
            A.block(5, 0, 1, 3) = v[2] * matK[0];
            A.block(5, 3, 1, 3) = -v[2] * matK[1];
            A.block(6, 6, 1, 6) << _cx(k - 1, j, i) - p[3][0], _cy(k - 1, j, i) - p[3][1], _cz(k - 1, j, i) - p[3][2], p[3][0] - _cx(k - 1, j, i - 1), p[3][1] - _cy(k - 1, j, i - 1), p[3][2] - _cz(k - 1, j, i - 1);
            A.block(7, 6, 1, 3) = v[3] * matK[2];
            A.block(7, 9, 1, 3) = -v[3] * matK[3];
            A.block(8, 15, 1, 6) << _cx(k, j - 1, i) - p[4][0], _cy(k, j - 1, i) - p[4][1], _cz(k, j - 1, i) - p[4][2], p[4][0] - _cx(k, j, i), p[4][1] - _cy(k, j, i), p[4][2] - _cz(k, j, i);
            A.block(9, 15, 1, 3) = v[4] * matK[5];
            A.block(9, 18, 1, 3) = -v[4] * matK[6];
            A.block(10, 3, 1, 6) << _cx(k - 1, j - 1, i) - p[5][0], _cy(k - 1, j - 1, i) - p[5][1], _cz(k - 1, j - 1, i) - p[5][2], p[5][0] - _cx(k - 1, j, i), p[5][1] - _cy(k - 1, j, i), p[5][2] - _cz(k - 1, j, i);
            A.block(11, 3, 1, 3) = v[5] * matK[1];
            A.block(11, 6, 1, 3) = -v[5] * matK[2];
            A.block(12, 0, 1, 3) << _cx(k - 1, j - 1, i - 1) - p[6][0], _cy(k - 1, j - 1, i - 1) - p[6][1], _cz(k - 1, j - 1, i - 1) - p[6][2];
            A.block(12, 9, 1, 3) << p[6][0] - _cx(k - 1, j, i - 1), p[6][1] - _cy(k - 1, j, i - 1), p[6][2] - _cz(k - 1, j, i - 1);
            A.block(13, 0, 1, 3) = v[6] * matK[0];
            A.block(13, 9, 1, 3) = -v[6] * matK[3];
            A.block(14, 12, 1, 3) << _cx(k, j - 1, i - 1) - p[7][0], _cy(k, j - 1, i - 1) - p[7][1], _cz(k, j - 1, i - 1) - p[7][2];
            A.block(14, 21, 1, 3) << p[7][0] - _cx(k, j, i - 1), p[7][1] - _cy(k, j, i - 1), p[7][2] - _cz(k, j, i - 1);
            A.block(15, 12, 1, 3) = v[7] * matK[4];
            A.block(15, 21, 1, 3) = -v[7] * matK[7];
            A.block(16, 6, 1, 3) << _cx(k - 1, j, i) - p[8][0], _cy(k - 1, j, i) - p[8][1], _cz(k - 1, j, i) - p[8][2];
            A.block(16, 18, 1, 3) << p[8][0] - _cx(k, j, i), p[8][1] - _cy(k, j, i), p[8][2] - _cz(k, j, i);
            A.block(17, 6, 1, 3) = v[8] * matK[2];
            A.block(17, 18, 1, 3) = -v[8] * matK[6];
            A.block(18, 9, 1, 3) << _cx(k - 1, j, i - 1) - p[9][0], _cy(k - 1, j, i - 1) - p[9][1], _cz(k - 1, j, i - 1) - p[9][2];
            A.block(18, 21, 1, 3) << p[9][0] - _cx(k, j, i - 1), p[9][1] - _cy(k, j, i - 1), p[9][2] - _cz(k, j, i - 1);
            A.block(19, 9, 1, 3) = v[9] * matK[3];
            A.block(19, 21, 1, 3) = -v[9] * matK[7];
            A.block(20, 0, 1, 3) << _cx(k - 1, j - 1, i - 1) - p[10][0], _cy(k - 1, j - 1, i - 1) - p[10][1], _cz(k - 1, j - 1, i - 1) - p[10][2];
            A.block(20, 12, 1, 3) << p[10][0] - _cx(k, j - 1, i - 1), p[10][1] - _cy(k, j - 1, i - 1), p[10][2] - _cz(k, j - 1, i - 1);
            A.block(21, 0, 1, 3) = v[10] * matK[0];
            A.block(21, 12, 1, 3) = -v[10] * matK[4];
            A.block(22, 3, 1, 3) << _cx(k - 1, j - 1, i) - p[11][0], _cy(k - 1, j - 1, i) - p[11][1], _cz(k - 1, j - 1, i) - p[11][2];
            A.block(22, 15, 1, 3) << p[11][0] - _cx(k, j - 1, i), p[11][1] - _cy(k, j - 1, i), p[11][2] - _cz(k, j - 1, i);
            A.block(23, 3, 1, 3) = v[11] * matK[1];
            A.block(23, 15, 1, 3) = -v[11] * matK[5];
            A = A.inverse();
            Eigen::RowVectorXd r = (ratio[0] * v[6] + ratio[3] * v[10]) * matK[0] * A.topRows(3);
            C(idx[0], idx[0]) += r[4] + r[12] + r[20];
            C(idx[0], idx[1]) += -r[4] + r[10] + r[22];
            C(idx[0], idx[2]) += r[6] - r[10] + r[16];
            C(idx[0], idx[3]) += -r[6] - r[12] + r[18];
            C(idx[0], idx[4]) += r[2] + r[14] - r[20];
            C(idx[0], idx[5]) += -r[2] + r[8] - r[22];
            C(idx[0], idx[6]) += r[0] - r[8] - r[16];
            C(idx[0], idx[7]) += -r[0] - r[14] - r[18];
            r = (-ratio[0] * v[6] + ratio[1] * v[9]) * matK[3] * A.block(9, 0, 3, 24);
            C(idx[3], idx[0]) += r[4] + r[12] + r[20];
            C(idx[3], idx[1]) += -r[4] + r[10] + r[22];
            C(idx[3], idx[2]) += r[6] - r[10] + r[16];
            C(idx[3], idx[3]) += -r[6] - r[12] + r[18];
            C(idx[3], idx[4]) += r[2] + r[14] - r[20];
            C(idx[3], idx[5]) += -r[2] + r[8] - r[22];
            C(idx[3], idx[6]) += r[0] - r[8] - r[16];
            C(idx[3], idx[7]) += -r[0] - r[14] - r[18];
            r = (ratio[2] * v[7] - ratio[3] * v[10]) * matK[4] * A.block(12, 0, 3, 24);
            C(idx[4], idx[0]) += r[4] + r[12] + r[20];
            C(idx[4], idx[1]) += -r[4] + r[10] + r[22];
            C(idx[4], idx[2]) += r[6] - r[10] + r[16];
            C(idx[4], idx[3]) += -r[6] - r[12] + r[18];
            C(idx[4], idx[4]) += r[2] + r[14] - r[20];
            C(idx[4], idx[5]) += -r[2] + r[8] - r[22];
            C(idx[4], idx[6]) += r[0] - r[8] - r[16];
            C(idx[4], idx[7]) += -r[0] - r[14] - r[18];
            r = (-ratio[2] * v[7] - ratio[1] * v[9]) * matK[7] * A.block(21, 0, 3, 24);
            C(idx[7], idx[0]) += r[4] + r[12] + r[20];
            C(idx[7], idx[1]) += -r[4] + r[10] + r[22];
            C(idx[7], idx[2]) += r[6] - r[10] + r[16];
            C(idx[7], idx[3]) += -r[6] - r[12] + r[18];
            C(idx[7], idx[4]) += r[2] + r[14] - r[20];
            C(idx[7], idx[5]) += -r[2] + r[8] - r[22];
            C(idx[7], idx[6]) += r[0] - r[8] - r[16];
            C(idx[7], idx[7]) += -r[0] - r[14] - r[18];
            // 内部角点
        }
        break;
    case Axis::YPOSITIVE:
        assert(j < ny);
        break;
    case Axis::YNEGATIVE:
        assert(j > 0);
        break;
    case Axis::ZPOSITIVE:
        assert(k < nz);
        break;
    case Axis::ZNEGATIVE:
        assert(k > 0);
        break;
    }
}

void _get_Qx(const Ktensor *perm, const double *verts, int divn, int &Iter, double &Q)
{
    std::cout << "divn = " << divn << std::endl;
    nx = nxx * divn;
    ny = nyy * divn;
    nz = nzz * divn;
    int n, nnz, count;

    _verts1 = new double[nx * ny * nz * 8 * 3];
    verts1 = std::mdspan<double, std::extents<size_t, -1, -1, -1, 8, 3>>(_verts1, nz, ny, nx);
    pem1 = new Ktensor[nx * ny * nz]();
    double *p = new double[nx * ny * nz]();
    Ptr = new int[nx * ny * nz + 1]();
    Idx = new int[27 * nx * ny * nz]();
    Val = new double[27 * nx * ny * nz]();
    B = new double[nx * ny * nz]();
    _get_divide(verts, _verts1, perm, pem1, divn);
    cx = new double[nx * ny * nz];
    cy = new double[nx * ny * nz];
    cz = new double[nx * ny * nz];

    _cx = std::mdspan(cx, nz, ny, nx);
    _cy = std::mdspan(cy, nz, ny, nx);
    _cz = std::mdspan(cz, nz, ny, nx);

    /*double para_8[8];
    int idx_8[8];
    _get_vec(0, 0, 0, pem1, _verts1, 8, para_8, idx_8, 2);*/

    /**
     * Calculate centroids.
     */
    for (int k = 0; k < nz; ++k)
    {
        for (int j = 0; j < ny; ++j)
        {
            for (int i = 0; i < nx; ++i)
            {
                _cx(k, j, i) = 0.0;
                _cy(k, j, i) = 0.0;
                _cz(k, j, i) = 0.0;

                for (int p = 0; p < 8; ++p)
                {
                    _cx(k, j, i) += verts1(k, j, i, p, 0);
                    _cy(k, j, i) += verts1(k, j, i, p, 1);
                    _cz(k, j, i) += verts1(k, j, i, p, 2);
                }

                _cx(k, j, i) /= 8;
                _cy(k, j, i) /= 8;
                _cz(k, j, i) /= 8;
            }
        }
    }

    // int *actnum = new int[nx * ny * nz];
    // std::fill(actnum, actnum + nx * ny * nz, 1);

    //_get_tran(nx, ny, nz, actnum, pem1, _verts1, cx, cy, cz, tranx, trany, tranz);

    double epsilonP = 1.0e-6, residual = 0.0001;
    n = 0;
    nnz = 0;
    Ptr[0] = 0;
    for (int k = 0; k < nz; k++)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int i = 0; i < nx; i++)
            {
                count = 0;
                for (int z = -1; z <= 1; ++z)
                    for (int y = -1; y <= 1; ++y)
                        for (int x = -1; x <= 1; ++x)
                        {
                            int ii = k + z;
                            int jj = j + y;
                            int ii_ = i + x;
                            if (ii >= 0 && ii < nz && jj >= 0 && jj < ny && ii_ >= 0 && ii_ < nx)
                            {
                                Idx[nnz] = ii * nx * ny + jj * nx + ii_;
                                nnz++;
                                count++;
                            }
                        }
                Ptr[n + 1] = Ptr[n] + count;
                ++n;
            }
        }
    }
    double pb = 0.0;
    for (int k = 0; k <= nz; k++)
        for (int j = 0; j <= ny; j++)
            for (int i = 0; i <= nx; i++)
            {
                int cur = k * nx * ny + j * nx + i;
                // 边界网格的边界压力

                if (i == 0)
                {
                    pb = Plow;
                }
                else if (i == nx)
                {
                    pb = Phigh;
                }
                // 3个面的点，2个未知数,8个点
                if (i == 0 && j == 0 && k == 0)
                {
                    // 0,4,8,交接面中点下标
                    double p0[3];
                    double p4[3];
                    double p8[3];
                    _get_centroid(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 2, 0), &verts1(k, j, i, 4, 0), &verts1(k, j, i, 6, 0), p0);
                    _get_centroid(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 1, 0), &verts1(k, j, i, 4, 0), &verts1(k, j, i, 5, 0), p4);
                    _get_centroid(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 1, 0), &verts1(k, j, i, 2, 0), &verts1(k, j, i, 3, 0), p8);
                    // 1,2,5,交接边点下标
                    double pp1[3];
                    double pp2[3];
                    double pp5[3];
                    _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 1, 0), pp1);
                    _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 2, 0), pp2);
                    _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 4, 0), pp5);
                    double n0[3];
                    double n4[3];
                    double n8[3];
                    double n48[3];
                    // 通过两个绝热面的流量为0,流量的方向可根据其法向量确定
                    _get_surface_normal(&verts1(k, j, i, 0, 0), pp2, p0, pp5, n0, Axis::XPOSITIVE);
                    _get_surface_normal(&verts1(k, j, i, 0, 0), pp1, p4, pp5, n4, Axis::YPOSITIVE);
                    _get_surface_normal(&verts1(k, j, i, 0, 0), pp1, p8, pp2, n8, Axis::ZPOSITIVE);
                    _cross_product(n4, n8, n48);
                    Eigen::Matrix3d matK;
                    _get_matK(pem1[cur], matK);
                    Eigen::RowVector3d Dx;
                    Eigen::Vector3d Dn;
                    Dx << p0[0] - _cx(k, j, i), p0[1] - _cy(k, j, i), p0[2] - _cz(k, j, i);
                    Dn << n48[0], n48[1], n48[2];
                    double cA = Dx * matK.inverse() * Dn;
                    cA = _dot_product(n0, n48) / cA;
                    B[cur] += pb * cA;
                    C(cur, cur) += cA;
                }
                else if (i == 0 && j == ny && k == 0)
                {
                    cur -= nx;
                    // 1,4,11交接面中点下标
                    double p1[3];
                    double p4[3];
                    double p11[3];
                    _get_centroid(&verts1(k, j - 1, i, 2, 0), &verts1(k, j - 1, i, 3, 0), &verts1(k, j - 1, i, 6, 0), &verts1(k, j - 1, i, 7, 0), p4);
                    _get_centroid(&verts1(k, j - 1, i, 0, 0), &verts1(k, j - 1, i, 2, 0), &verts1(k, j - 1, i, 4, 0), &verts1(k, j - 1, i, 6, 0), p1);
                    _get_centroid(&verts1(k, j - 1, i, 0, 0), &verts1(k, j - 1, i, 1, 0), &verts1(k, j - 1, i, 2, 0), &verts1(k, j - 1, i, 3, 0), p11);
                    // 0,1,5,交接边点下标
                    double pp0[3];
                    double pp1[3];
                    double pp5[3];
                    _get_midpoint(&verts1(k, j - 1, i, 0, 0), &verts1(k, j - 1, i, 2, 0), pp0);
                    _get_midpoint(&verts1(k, j - 1, i, 2, 0), &verts1(k, j - 1, i, 3, 0), pp1);
                    _get_midpoint(&verts1(k, j - 1, i, 2, 0), &verts1(k, j - 1, i, 6, 0), pp5);
                    double n1[3];
                    double n4[3];
                    double n11[3];
                    double n411[3];
                    // 通过两个绝热面的流量为0,流量的方向可根据其法向量确定
                    _get_surface_normal(&verts1(k, j - 1, i, 2, 0), pp0, p1, pp5, n1, Axis::XPOSITIVE);
                    _get_surface_normal(&verts1(k, j - 1, i, 2, 0), pp1, p4, pp5, n4, Axis::YPOSITIVE);
                    _get_surface_normal(&verts1(k, j - 1, i, 2, 0), pp1, p11, pp0, n11, Axis::ZPOSITIVE);
                    _cross_product(n4, n11, n411);
                    Eigen::Matrix3d matK;
                    _get_matK(pem1[cur], matK);
                    Eigen::RowVector3d Dx;
                    Eigen::Vector3d Dn;
                    Dx << p1[0] - _cx(k, j - 1, i), p1[1] - _cy(k, j - 1, i), p1[2] - _cz(k, j - 1, i);
                    Dn << n411[0], n411[1], n411[2];
                    double cA = Dx * matK.inverse() * Dn;
                    cA = _dot_product(n1, n411) / cA;
                    B[cur] += pb * cA;
                    C(cur, cur) += cA;
                }
                else if (i == 0 && j == 0 && k == nz)
                {
                    cur -= nx * ny;
                    // 3,5,8交接面中点下标
                    double p3[3];
                    double p5[3];
                    double p8[3];
                    _get_centroid(&verts1(k - 1, j, i, 0, 0), &verts1(k - 1, j, i, 2, 0), &verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 6, 0), p3);
                    _get_centroid(&verts1(k - 1, j, i, 0, 0), &verts1(k - 1, j, i, 1, 0), &verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 5, 0), p5);
                    _get_centroid(&verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 5, 0), &verts1(k - 1, j, i, 6, 0), &verts1(k - 1, j, i, 7, 0), p8);
                    // 1,2,4交接边点下标
                    double pp1[3];
                    double pp2[3];
                    double pp4[3];
                    _get_midpoint(&verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 5, 0), pp1);
                    _get_midpoint(&verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 6, 0), pp2);
                    _get_midpoint(&verts1(k - 1, j, i, 0, 0), &verts1(k - 1, j, i, 4, 0), pp4);
                    double n3[3];
                    double n5[3];
                    double n8[3];
                    double n58[3];
                    // 通过两个绝热面的流量为0,流量的方向可根据其法向量确定
                    _get_surface_normal(&verts1(k - 1, j, i, 4, 0), pp4, p3, pp2, n3, Axis::XPOSITIVE);
                    _get_surface_normal(&verts1(k - 1, j, i, 4, 0), pp1, p5, pp4, n5, Axis::YPOSITIVE);
                    _get_surface_normal(&verts1(k - 1, j, i, 4, 0), pp1, p8, pp2, n8, Axis::ZPOSITIVE);
                    _cross_product(n5, n8, n58);
                    Eigen::Matrix3d matK;
                    _get_matK(pem1[cur], matK);
                    Eigen::RowVector3d Dx;
                    Eigen::Vector3d Dn;
                    Dx << p3[0] - _cx(k - 1, j, i), p3[1] - _cy(k - 1, j, i), p3[2] - _cz(k - 1, j, i);
                    Dn << n58[0], n58[1], n58[2];
                    double cA = Dx * matK.inverse() * Dn;
                    cA = _dot_product(n3, n58) / cA;
                    B[cur] += pb * cA;
                    C(cur, cur) += cA;
                }
                else if (i == 0 && j == ny && k == nz)
                {
                    cur -= nx * ny + nx;
                    // 2,5,11交接面中点下标
                    double p2[3];
                    double p5[3];
                    double p11[3];
                    _get_centroid(&verts1(k - 1, j - 1, i, 0, 0), &verts1(k - 1, j - 1, i, 2, 0), &verts1(k - 1, j - 1, i, 4, 0), &verts1(k - 1, j - 1, i, 6, 0), p2);
                    _get_centroid(&verts1(k - 1, j - 1, i, 2, 0), &verts1(k - 1, j - 1, i, 3, 0), &verts1(k - 1, j - 1, i, 6, 0), &verts1(k - 1, j - 1, i, 7, 0), p5);
                    _get_centroid(&verts1(k - 1, j - 1, i, 4, 0), &verts1(k - 1, j - 1, i, 5, 0), &verts1(k - 1, j - 1, i, 6, 0), &verts1(k - 1, j - 1, i, 7, 0), p11);
                    // 0,1,4交接边点下标
                    double pp0[3];
                    double pp1[3];
                    double pp4[3];
                    _get_midpoint(&verts1(k - 1, j - 1, i, 6, 0), &verts1(k - 1, j - 1, i, 4, 0), pp0);
                    _get_midpoint(&verts1(k - 1, j - 1, i, 6, 0), &verts1(k - 1, j - 1, i, 7, 0), pp1);
                    _get_midpoint(&verts1(k - 1, j - 1, i, 2, 0), &verts1(k - 1, j - 1, i, 6, 0), pp4);
                    double n2[3];
                    double n5[3];
                    double n11[3];
                    double n511[3];
                    // 通过两个绝热面的流量为0,流量的方向可根据其法向量确定
                    _get_surface_normal(&verts1(k - 1, j - 1, i, 6, 0), pp0, p2, pp4, n2, Axis::XPOSITIVE);
                    _get_surface_normal(&verts1(k - 1, j - 1, i, 6, 0), pp1, p5, pp4, n5, Axis::YPOSITIVE);
                    _get_surface_normal(&verts1(k - 1, j - 1, i, 6, 0), pp1, p11, pp0, n11, Axis::ZPOSITIVE);
                    _cross_product(n5, n11, n511);
                    Eigen::Matrix3d matK;
                    _get_matK(pem1[cur], matK);
                    Eigen::RowVector3d Dx;
                    Eigen::Vector3d Dn;
                    Dx << p2[0] - _cx(k - 1, j - 1, i), p2[1] - _cy(k - 1, j - 1, i), p2[2] - _cz(k - 1, j - 1, i);
                    Dn << n511[0], n511[1], n511[2];
                    double cA = Dx * matK.inverse() * Dn;
                    cA = _dot_product(n2, n511) / cA;
                    B[cur] += pb * cA;
                    C(cur, cur) += cA;
                }
                else if (i == nx && j == 0 && k == 0)
                {
                    cur -= 1;
                    // 0,7,9交接面中点下标
                    double p0[3];
                    double p7[3];
                    double p9[3];
                    _get_centroid(&verts1(k, j, i - 1, 1, 0), &verts1(k, j, i - 1, 3, 0), &verts1(k, j, i - 1, 5, 0), &verts1(k, j, i - 1, 7, 0), p0);
                    _get_centroid(&verts1(k, j, i - 1, 0, 0), &verts1(k, j, i - 1, 1, 0), &verts1(k, j, i - 1, 4, 0), &verts1(k, j, i - 1, 5, 0), p7);
                    _get_centroid(&verts1(k, j, i - 1, 0, 0), &verts1(k, j, i - 1, 1, 0), &verts1(k, j, i - 1, 2, 0), &verts1(k, j, i - 1, 3, 0), p9);
                    // 2,3,5交接边点下标
                    double pp2[3];
                    double pp3[3];
                    double pp5[3];
                    _get_midpoint(&verts1(k, j, i - 1, 1, 0), &verts1(k, j, i - 1, 3, 0), pp2);
                    _get_midpoint(&verts1(k, j, i - 1, 0, 0), &verts1(k, j, i - 1, 1, 0), pp3);
                    _get_midpoint(&verts1(k, j, i - 1, 1, 0), &verts1(k, j, i - 1, 5, 0), pp5);
                    double n0[3];
                    double n7[3];
                    double n9[3];
                    double n79[3];
                    // 通过两个绝热面的流量为0,流量的方向可根据其法向量确定
                    _get_surface_normal(&verts1(k, j, i - 1, 1, 0), pp2, p0, pp5, n0, Axis::XPOSITIVE);
                    _get_surface_normal(&verts1(k, j, i - 1, 1, 0), pp3, p7, pp5, n7, Axis::YPOSITIVE);
                    _get_surface_normal(&verts1(k, j, i - 1, 1, 0), pp3, p9, pp2, n9, Axis::ZPOSITIVE);
                    _cross_product(n7, n9, n79);
                    Eigen::Matrix3d matK;
                    _get_matK(pem1[cur], matK);
                    Eigen::RowVector3d Dx;
                    Eigen::Vector3d Dn;
                    Dx << p0[0] - _cx(k, j, i - 1), p0[1] - _cy(k, j, i - 1), p0[2] - _cz(k, j, i - 1);
                    Dn << n79[0], n79[1], n79[2];
                    double cA = Dx * matK.inverse() * Dn;
                    cA = _dot_product(n0, n79) / cA;
                    B[cur] -= pb * cA;
                    C(cur, cur) -= cA;
                }
                else if (i == nx && j == ny && k == 0)
                {
                    cur -= nx + 1;
                    // 1,7,10交接面中点下标
                    double p1[3];
                    double p7[3];
                    double p10[3];
                    _get_centroid(&verts1(k, j - 1, i - 1, 2, 0), &verts1(k, j - 1, i - 1, 3, 0), &verts1(k, j - 1, i - 1, 6, 0), &verts1(k, j - 1, i - 1, 7, 0), p7);
                    _get_centroid(&verts1(k, j - 1, i - 1, 1, 0), &verts1(k, j - 1, i - 1, 3, 0), &verts1(k, j - 1, i - 1, 5, 0), &verts1(k, j - 1, i - 1, 7, 0), p1);
                    _get_centroid(&verts1(k, j - 1, i - 1, 1, 0), &verts1(k, j - 1, i - 1, 0, 0), &verts1(k, j - 1, i - 1, 2, 0), &verts1(k, j - 1, i - 1, 3, 0), p10);
                    // 0,3,5交接边点下标
                    double pp0[3];
                    double pp3[3];
                    double pp5[3];
                    _get_midpoint(&verts1(k, j - 1, i - 1, 3, 0), &verts1(k, j - 1, i - 1, 1, 0), pp0);
                    _get_midpoint(&verts1(k, j - 1, i - 1, 2, 0), &verts1(k, j - 1, i - 1, 3, 0), pp3);
                    _get_midpoint(&verts1(k, j - 1, i - 1, 3, 0), &verts1(k, j - 1, i - 1, 7, 0), pp5);
                    double n1[3];
                    double n7[3];
                    double n10[3];
                    double n710[3];
                    // 通过两个绝热面的流量为0,流量的方向可根据其法向量确定
                    _get_surface_normal(&verts1(k, j - 1, i - 1, 3, 0), pp0, p1, pp5, n1, Axis::XPOSITIVE);
                    _get_surface_normal(&verts1(k, j - 1, i - 1, 3, 0), pp3, p7, pp5, n7, Axis::YPOSITIVE);
                    _get_surface_normal(&verts1(k, j - 1, i - 1, 3, 0), pp3, p10, pp0, n10, Axis::ZPOSITIVE);
                    _cross_product(n7, n10, n710);
                    Eigen::Matrix3d matK;
                    _get_matK(pem1[cur], matK);
                    Eigen::RowVector3d Dx;
                    Eigen::Vector3d Dn;
                    Dx << p1[0] - _cx(k, j - 1, i - 1), p1[1] - _cy(k, j - 1, i - 1), p1[2] - _cz(k, j - 1, i - 1);
                    Dn << n710[0], n710[1], n710[2];
                    double cA = Dx * matK.inverse() * Dn;
                    cA = _dot_product(n1, n710) / cA;
                    B[cur] -= pb * cA;
                    C(cur, cur) -= cA;
                }
                else if (i == nx && j == 0 && k == nz)
                {
                    cur -= nx * ny + 1;
                    // 3,6,9交接面中点下标
                    double p3[3];
                    double p6[3];
                    double p9[3];
                    _get_centroid(&verts1(k - 1, j, i - 1, 1, 0), &verts1(k - 1, j, i - 1, 3, 0), &verts1(k - 1, j, i - 1, 5, 0), &verts1(k - 1, j, i - 1, 7, 0), p3);
                    _get_centroid(&verts1(k - 1, j, i - 1, 0, 0), &verts1(k - 1, j, i - 1, 1, 0), &verts1(k - 1, j, i - 1, 4, 0), &verts1(k - 1, j, i - 1, 5, 0), p6);
                    _get_centroid(&verts1(k - 1, j, i - 1, 4, 0), &verts1(k - 1, j, i - 1, 5, 0), &verts1(k - 1, j, i - 1, 6, 0), &verts1(k - 1, j, i - 1, 7, 0), p9);
                    // 2,3,4交接边点下标
                    double pp2[3];
                    double pp3[3];
                    double pp4[3];
                    _get_midpoint(&verts1(k - 1, j, i - 1, 5, 0), &verts1(k - 1, j, i - 1, 7, 0), pp2);
                    _get_midpoint(&verts1(k - 1, j, i - 1, 5, 0), &verts1(k - 1, j, i - 1, 4, 0), pp3);
                    _get_midpoint(&verts1(k - 1, j, i - 1, 5, 0), &verts1(k - 1, j, i - 1, 1, 0), pp4);
                    double n3[3];
                    double n6[3];
                    double n9[3];
                    double n69[3];
                    // 通过两个绝热面的流量为0,流量的方向可根据其法向量确定
                    _get_surface_normal(&verts1(k - 1, j, i - 1, 5, 0), pp4, p3, pp2, n3, Axis::XPOSITIVE);
                    _get_surface_normal(&verts1(k - 1, j, i - 1, 5, 0), pp3, p6, pp4, n6, Axis::YPOSITIVE);
                    _get_surface_normal(&verts1(k - 1, j, i - 1, 5, 0), pp3, p9, pp2, n9, Axis::ZPOSITIVE);
                    _cross_product(n6, n9, n69);
                    Eigen::Matrix3d matK;
                    _get_matK(pem1[cur], matK);
                    Eigen::RowVector3d Dx;
                    Eigen::Vector3d Dn;
                    Dx << p3[0] - _cx(k - 1, j, i - 1), p3[1] - _cy(k - 1, j, i - 1), p3[2] - _cz(k - 1, j, i - 1);
                    Dn << n69[0], n69[1], n69[2];
                    double cA = Dx * matK.inverse() * Dn;
                    cA = _dot_product(n3, n69) / cA;
                    B[cur] -= pb * cA;
                    C(cur, cur) -= cA;
                }
                else if (i == nx && j == ny && k == nz)
                {
                    cur -= nx * ny + nx + 1;
                    // 2,6,10交接面中点下标
                    double p2[3];
                    double p6[3];
                    double p10[3];
                    _get_centroid(&verts1(k - 1, j - 1, i - 1, 1, 0), &verts1(k - 1, j - 1, i - 1, 3, 0), &verts1(k - 1, j - 1, i - 1, 5, 0), &verts1(k - 1, j - 1, i - 1, 7, 0), p2);
                    _get_centroid(&verts1(k - 1, j - 1, i - 1, 2, 0), &verts1(k - 1, j - 1, i - 1, 3, 0), &verts1(k - 1, j - 1, i - 1, 6, 0), &verts1(k - 1, j - 1, i - 1, 7, 0), p6);
                    _get_centroid(&verts1(k - 1, j - 1, i - 1, 4, 0), &verts1(k - 1, j - 1, i - 1, 5, 0), &verts1(k - 1, j - 1, i - 1, 6, 0), &verts1(k - 1, j - 1, i - 1, 7, 0), p10);
                    // 0,3,4交接边点下标
                    double pp0[3];
                    double pp3[3];
                    double pp4[3];
                    _get_midpoint(&verts1(k - 1, j - 1, i - 1, 5, 0), &verts1(k - 1, j - 1, i - 1, 7, 0), pp0);
                    _get_midpoint(&verts1(k - 1, j - 1, i - 1, 6, 0), &verts1(k - 1, j - 1, i - 1, 7, 0), pp3);
                    _get_midpoint(&verts1(k - 1, j - 1, i - 1, 3, 0), &verts1(k - 1, j - 1, i - 1, 7, 0), pp4);
                    double n2[3];
                    double n6[3];
                    double n10[3];
                    double n610[3];
                    // 通过两个绝热面的流量为0,流量的方向可根据其法向量确定
                    _get_surface_normal(&verts1(k - 1, j - 1, i - 1, 7, 0), pp0, p2, pp4, n2, Axis::XPOSITIVE);
                    _get_surface_normal(&verts1(k - 1, j - 1, i - 1, 7, 0), pp3, p6, pp4, n6, Axis::YPOSITIVE);
                    _get_surface_normal(&verts1(k - 1, j - 1, i - 1, 7, 0), pp3, p10, pp0, n10, Axis::ZPOSITIVE);
                    _cross_product(n6, n10, n610);
                    Eigen::Matrix3d matK;
                    _get_matK(pem1[cur], matK);
                    Eigen::RowVector3d Dx;
                    Eigen::Vector3d Dn;
                    Dx << p2[0] - _cx(k - 1, j - 1, i - 1), p2[1] - _cy(k - 1, j - 1, i - 1), p2[2] - _cz(k - 1, j - 1, i - 1);
                    Dn << n610[0], n610[1], n610[2];
                    double cA = Dx * matK.inverse() * Dn;
                    cA = _dot_product(n2, n610) / cA;
                    B[cur] -= pb * cA;
                    C(cur, cur) -= cA;
                }
                // 12条边上其余点的处理
                else if (j == 0 && k == 0)
                {
                    // 两个网格中心下标
                    int i6 = cur, i7 = cur - 1;
                    // 0,4,7,8,9交接面中点下标
                    double p0[3];
                    double p4[3];
                    double p7[3];
                    double p8[3];
                    double p9[3];
                    _get_centroid(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 2, 0), &verts1(k, j, i, 4, 0), &verts1(k, j, i, 6, 0), p0);
                    _get_centroid(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 1, 0), &verts1(k, j, i, 4, 0), &verts1(k, j, i, 5, 0), p4);
                    _get_centroid(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 1, 0), &verts1(k, j, i, 2, 0), &verts1(k, j, i, 3, 0), p8);
                    _get_centroid(&verts1(k, j, i - 1, 0, 0), &verts1(k, j, i - 1, 1, 0), &verts1(k, j, i - 1, 4, 0), &verts1(k, j, i - 1, 5, 0), p7);
                    _get_centroid(&verts1(k, j, i - 1, 0, 0), &verts1(k, j, i - 1, 1, 0), &verts1(k, j, i - 1, 2, 0), &verts1(k, j, i - 1, 3, 0), p9);
                    // 1,2,3,5,交接边点下标
                    double pp1[3];
                    double pp2[3];
                    double pp3[3];
                    double pp5[3];
                    _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 1, 0), pp1);
                    _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 2, 0), pp2);
                    _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 4, 0), pp5);
                    _get_midpoint(&verts1(k, j, i - 1, 0, 0), &verts1(k, j, i - 1, 1, 0), pp3);
                    double n0[3];
                    double n4[3];
                    double n8[3];
                    double n48[3];
                    double n7[3];
                    double n9[3];
                    double n79[3];
                    // 通过两个绝热面的流量为0,流量的方向可根据其法向量确定
                    _get_surface_normal(&verts1(k, j, i, 0, 0), pp2, p0, pp5, n0, Axis::XPOSITIVE);
                    _get_surface_normal(&verts1(k, j, i, 0, 0), pp1, p4, pp5, n4, Axis::YPOSITIVE);
                    _get_surface_normal(&verts1(k, j, i, 0, 0), pp1, p8, pp2, n8, Axis::ZPOSITIVE);
                    _get_surface_normal(&verts1(k, j, i - 1, 1, 0), pp3, p7, pp5, n7, Axis::YPOSITIVE);
                    _get_surface_normal(&verts1(k, j, i - 1, 1, 0), pp3, p9, pp2, n9, Axis::ZPOSITIVE);
                    _cross_product(n4, n8, n48);
                    _cross_product(n7, n9, n79);
                    Eigen::Matrix3d matK6, matK7;
                    _get_matK(pem1[i6], matK6);
                    _get_matK(pem1[i7], matK7);
                    Eigen::RowVector3d Dx6, Dx7;
                    Eigen::Vector3d Dn48, Dn79;
                    Dx6 << p0[0] - _cx(k, j, i), p0[1] - _cy(k, j, i), p0[2] - _cz(k, j, i);
                    Dn48 << n48[0], n48[1], n48[2];
                    double cA6 = Dx6 * matK6.inverse() * Dn48;
                    cA6 = _dot_product(n0, n48) / cA6;
                    Dx7 << p0[0] - _cx(k, j, i - 1), p0[1] - _cy(k, j, i - 1), p0[2] - _cz(k, j, i - 1);
                    Dn79 << n79[0], n79[1], n79[2];
                    double cA7 = Dx7 * matK7.inverse() * Dn79;
                    cA7 = _dot_product(n0, n79) / cA7;
                    double cA = 0.5 * _harmonic(cA6, -cA7);
                    C(i6, i6) += cA;
                    C(i7, i7) += cA;
                    C(i6, i7) -= cA;
                    C(i7, i6) -= cA;
                }
                else if (j == ny && k == 0)
                {
                    // 两个网格中心下标
                    int i5 = cur - nx, i4 = cur - nx - 1;
                    // 1,4,7,10,11交接面中点下标
                    double p1[3];
                    double p4[3];
                    double p7[3];
                    double p10[3];
                    double p11[3];
                    _get_centroid(&verts1(k, j - 1, i, 2, 0), &verts1(k, j - 1, i, 3, 0), &verts1(k, j - 1, i, 6, 0), &verts1(k, j - 1, i, 7, 0), p4);
                    _get_centroid(&verts1(k, j - 1, i, 0, 0), &verts1(k, j - 1, i, 2, 0), &verts1(k, j - 1, i, 4, 0), &verts1(k, j - 1, i, 6, 0), p1);
                    _get_centroid(&verts1(k, j - 1, i, 0, 0), &verts1(k, j - 1, i, 1, 0), &verts1(k, j - 1, i, 2, 0), &verts1(k, j - 1, i, 3, 0), p11);
                    _get_centroid(&verts1(k, j - 1, i - 1, 2, 0), &verts1(k, j - 1, i - 1, 3, 0), &verts1(k, j - 1, i - 1, 6, 0), &verts1(k, j - 1, i - 1, 7, 0), p7);
                    _get_centroid(&verts1(k, j - 1, i - 1, 1, 0), &verts1(k, j - 1, i - 1, 0, 0), &verts1(k, j - 1, i - 1, 2, 0), &verts1(k, j - 1, i - 1, 3, 0), p10);
                    // 0,1,3,5交接边点下标
                    double pp0[3];
                    double pp1[3];
                    double pp3[3];
                    double pp5[3];
                    _get_midpoint(&verts1(k, j - 1, i - 1, 3, 0), &verts1(k, j - 1, i - 1, 1, 0), pp0);
                    _get_midpoint(&verts1(k, j - 1, i - 1, 2, 0), &verts1(k, j - 1, i - 1, 3, 0), pp3);
                    _get_midpoint(&verts1(k, j - 1, i - 1, 3, 0), &verts1(k, j - 1, i - 1, 7, 0), pp5);
                    _get_midpoint(&verts1(k, j - 1, i, 2, 0), &verts1(k, j - 1, i, 3, 0), pp1);
                    double n1[3];
                    double n4[3];
                    double n11[3];
                    double n411[3];
                    double n7[3];
                    double n10[3];
                    double n710[3];
                    // 通过两个绝热面的流量为0,流量的方向可根据其法向量确定
                    _get_surface_normal(&verts1(k, j - 1, i - 1, 3, 0), pp0, p1, pp5, n1, Axis::XPOSITIVE);
                    _get_surface_normal(&verts1(k, j - 1, i - 1, 3, 0), pp3, p7, pp5, n7, Axis::YPOSITIVE);
                    _get_surface_normal(&verts1(k, j - 1, i - 1, 3, 0), pp3, p10, pp0, n10, Axis::ZPOSITIVE);
                    _get_surface_normal(&verts1(k, j - 1, i, 2, 0), pp1, p4, pp5, n4, Axis::YPOSITIVE);
                    _get_surface_normal(&verts1(k, j - 1, i, 2, 0), pp1, p11, pp0, n11, Axis::ZPOSITIVE);
                    _cross_product(n4, n11, n411);
                    _cross_product(n7, n10, n710);
                    Eigen::Matrix3d matK4, matK5;
                    _get_matK(pem1[i4], matK4);
                    _get_matK(pem1[i5], matK5);
                    Eigen::RowVector3d Dx4, Dx5;
                    Eigen::Vector3d Dn411, Dn710;
                    Dx5 << p1[0] - _cx(k, j - 1, i), p1[1] - _cy(k, j - 1, i), p1[2] - _cz(k, j - 1, i);
                    Dn411 << n411[0], n411[1], n411[2];
                    double cA5 = Dx5 * matK5.inverse() * Dn411;
                    cA5 = _dot_product(n1, n411) / cA5;
                    Dx4 << p1[0] - _cx(k, j - 1, i - 1), p1[1] - _cy(k, j - 1, i - 1), p1[2] - _cz(k, j - 1, i - 1);
                    Dn710 << n710[0], n710[1], n710[2];
                    double cA4 = Dx4 * matK4.inverse() * Dn710;
                    cA4 = _dot_product(n1, n710) / cA4;
                    double cA = 0.5 * _harmonic(cA5, -cA4);
                    C(i5, i5) += cA;
                    C(i4, i4) += cA;
                    C(i5, i4) -= cA;
                    C(i4, i5) -= cA;
                }
                else if (j == 0 && k == nz)
                {
                    //  两个网格中心下标
                    int i2 = cur - nx * ny, i3 = cur - nx * ny - 1;
                    // 3,5,6,8,9,交接面中点下标
                    double p3[3];
                    double p5[3];
                    double p6[3];
                    double p8[3];
                    double p9[3];
                    _get_centroid(&verts1(k - 1, j, i, 0, 0), &verts1(k - 1, j, i, 2, 0), &verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 6, 0), p3);
                    _get_centroid(&verts1(k - 1, j, i, 0, 0), &verts1(k - 1, j, i, 1, 0), &verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 5, 0), p5);
                    _get_centroid(&verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 5, 0), &verts1(k - 1, j, i, 6, 0), &verts1(k - 1, j, i, 7, 0), p8);
                    _get_centroid(&verts1(k - 1, j, i - 1, 0, 0), &verts1(k - 1, j, i - 1, 1, 0), &verts1(k - 1, j, i - 1, 4, 0), &verts1(k - 1, j, i - 1, 5, 0), p6);
                    _get_centroid(&verts1(k - 1, j, i - 1, 4, 0), &verts1(k - 1, j, i - 1, 5, 0), &verts1(k - 1, j, i - 1, 6, 0), &verts1(k - 1, j, i - 1, 7, 0), p9);
                    // 1,2,3,4,交接边点下标
                    double pp1[3];
                    double pp2[3];
                    double pp3[3];
                    double pp4[3];
                    _get_midpoint(&verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 5, 0), pp1);
                    _get_midpoint(&verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 6, 0), pp2);
                    _get_midpoint(&verts1(k - 1, j, i - 1, 5, 0), &verts1(k - 1, j, i - 1, 4, 0), pp3);
                    _get_midpoint(&verts1(k - 1, j, i, 0, 0), &verts1(k - 1, j, i, 4, 0), pp4);
                    double n3[3];
                    double n5[3];
                    double n8[3];
                    double n58[3];
                    double n6[3];
                    double n9[3];
                    double n69[3];
                    // 通过两个绝热面的流量为0,流量的方向可根据其法向量确定
                    _get_surface_normal(&verts1(k - 1, j, i, 4, 0), pp4, p3, pp2, n3, Axis::XPOSITIVE);
                    _get_surface_normal(&verts1(k - 1, j, i, 4, 0), pp1, p5, pp4, n5, Axis::YPOSITIVE);
                    _get_surface_normal(&verts1(k - 1, j, i, 4, 0), pp1, p8, pp2, n8, Axis::ZPOSITIVE);
                    _get_surface_normal(&verts1(k - 1, j, i - 1, 5, 0), pp3, p6, pp4, n6, Axis::YPOSITIVE);
                    _get_surface_normal(&verts1(k - 1, j, i - 1, 5, 0), pp3, p9, pp2, n9, Axis::ZPOSITIVE);
                    _cross_product(n5, n8, n58);
                    _cross_product(n6, n9, n69);
                    Eigen::Matrix3d matK2, matK3;
                    _get_matK(pem1[i2], matK2);
                    _get_matK(pem1[i3], matK3);
                    Eigen::RowVector3d Dx2, Dx3;
                    Eigen::Vector3d Dn58, Dn69;
                    Dx2 << p3[0] - _cx(k - 1, j, i), p3[1] - _cy(k - 1, j, i), p3[2] - _cz(k - 1, j, i);
                    Dn58 << n58[0], n58[1], n58[2];
                    double cA2 = Dx2 * matK2.inverse() * Dn58;
                    cA2 = _dot_product(n3, n58) / cA2;
                    Dx3 << p3[0] - _cx(k - 1, j, i - 1), p3[1] - _cy(k - 1, j, i - 1), p3[2] - _cz(k - 1, j, i - 1);
                    Dn69 << n69[0], n69[1], n69[2];
                    double cA3 = Dx3 * matK3.inverse() * Dn69;
                    cA3 = _dot_product(n3, n69) / cA3;
                    double cA = 0.5 * _harmonic(cA2, -cA3);
                    C(i2, i2) += cA;
                    C(i3, i3) += cA;
                    C(i2, i3) -= cA;
                    C(i3, i2) -= cA;
                }
                else if (j == ny && k == nz)
                {
                    // 两个网格中心下标
                    int i1 = cur - nx * ny - nx, i0 = cur - nx * ny - nx - 1;
                    // 2,5,6,10,11交接面中点下标
                    double p2[3];
                    double p5[3];
                    double p6[3];
                    double p10[3];
                    double p11[3];
                    _get_centroid(&verts1(k - 1, j - 1, i, 2, 0), &verts1(k - 1, j - 1, i, 3, 0), &verts1(k - 1, j - 1, i, 6, 0), &verts1(k - 1, j - 1, i, 7, 0), p5);
                    _get_centroid(&verts1(k - 1, j - 1, i, 0, 0), &verts1(k - 1, j - 1, i, 2, 0), &verts1(k - 1, j - 1, i, 4, 0), &verts1(k - 1, j - 1, i, 6, 0), p2);
                    _get_centroid(&verts1(k - 1, j - 1, i, 4, 0), &verts1(k - 1, j - 1, i, 5, 0), &verts1(k - 1, j - 1, i, 6, 0), &verts1(k - 1, j - 1, i, 7, 0), p11);
                    _get_centroid(&verts1(k - 1, j - 1, i - 1, 2, 0), &verts1(k - 1, j - 1, i - 1, 3, 0), &verts1(k - 1, j - 1, i - 1, 6, 0), &verts1(k - 1, j - 1, i - 1, 7, 0), p6);
                    _get_centroid(&verts1(k - 1, j - 1, i - 1, 4, 0), &verts1(k - 1, j - 1, i - 1, 5, 0), &verts1(k - 1, j - 1, i - 1, 6, 0), &verts1(k - 1, j - 1, i - 1, 7, 0), p10);
                    // 0,1,3,4交接边点下标
                    double pp0[3];
                    double pp1[3];
                    double pp3[3];
                    double pp4[3];
                    _get_midpoint(&verts1(k - 1, j - 1, i, 6, 0), &verts1(k - 1, j - 1, i, 4, 0), pp0);
                    _get_midpoint(&verts1(k - 1, j - 1, i, 6, 0), &verts1(k - 1, j - 1, i, 7, 0), pp1);
                    _get_midpoint(&verts1(k - 1, j - 1, i - 1, 6, 0), &verts1(k - 1, j - 1, i - 1, 7, 0), pp3);
                    _get_midpoint(&verts1(k - 1, j - 1, i, 2, 0), &verts1(k - 1, j - 1, i, 6, 0), pp4);
                    double n2[3];
                    double n5[3];
                    double n11[3];
                    double n511[3];
                    double n6[3];
                    double n10[3];
                    double n610[3];
                    // 通过两个绝热面的流量为0,流量的方向可根据其法向量确定
                    _get_surface_normal(&verts1(k - 1, j - 1, i, 6, 0), pp0, p2, pp4, n2, Axis::XPOSITIVE);
                    _get_surface_normal(&verts1(k - 1, j - 1, i, 6, 0), pp1, p5, pp4, n5, Axis::YPOSITIVE);
                    _get_surface_normal(&verts1(k - 1, j - 1, i, 6, 0), pp1, p11, pp0, n11, Axis::ZPOSITIVE);
                    _get_surface_normal(&verts1(k - 1, j - 1, i - 1, 7, 0), pp3, p6, pp4, n6, Axis::YPOSITIVE);
                    _get_surface_normal(&verts1(k - 1, j - 1, i - 1, 7, 0), pp3, p10, pp0, n10, Axis::ZPOSITIVE);
                    _cross_product(n5, n11, n511);
                    _cross_product(n6, n10, n610);
                    Eigen::Matrix3d matK0, matK1;
                    _get_matK(pem1[i0], matK0);
                    _get_matK(pem1[i1], matK1);
                    Eigen::RowVector3d Dx0, Dx1;
                    Eigen::Vector3d Dn511, Dn610;
                    Dx1 << p2[0] - _cx(k - 1, j - 1, i), p2[1] - _cy(k - 1, j - 1, i), p2[2] - _cz(k - 1, j - 1, i);
                    Dn511 << n511[0], n511[1], n511[2];
                    double cA1 = Dx1 * matK1.inverse() * Dn511;
                    cA1 = _dot_product(n2, n511) / cA1;
                    Dx0 << p2[0] - _cx(k - 1, j - 1, i - 1), p2[1] - _cy(k - 1, j - 1, i - 1), p2[2] - _cz(k - 1, j - 1, i - 1);
                    Dn610 << n610[0], n610[1], n610[2];
                    double cA0 = Dx0 * matK0.inverse() * Dn610;
                    cA0 = _dot_product(n2, n610) / cA0;
                    double cA = 0.5 * _harmonic(cA1, -cA0);
                    C(i1, i1) += cA;
                    C(i0, i0) += cA;
                    C(i1, i0) -= cA;
                    C(i0, i1) -= cA;
                }
                else if (i == 0 && k == 0)
                {
                    // 两个网格中心下标
                    int i6 = cur, i5 = cur - nx;
                    // 0,1,4,8,11交接面中点下标
                    double p0[3];
                    double p1[3];
                    double p4[3];
                    double p8[3];
                    double p11[3];
                    _get_centroid(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 2, 0), &verts1(k, j, i, 4, 0), &verts1(k, j, i, 6, 0), p0);
                    _get_centroid(&verts1(k, j - 1, i, 0, 0), &verts1(k, j - 1, i, 2, 0), &verts1(k, j - 1, i, 4, 0), &verts1(k, j - 1, i, 6, 0), p1);
                    _get_centroid(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 1, 0), &verts1(k, j, i, 4, 0), &verts1(k, j, i, 5, 0), p4);
                    _get_centroid(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 1, 0), &verts1(k, j, i, 2, 0), &verts1(k, j, i, 3, 0), p8);
                    _get_centroid(&verts1(k, j - 1, i, 0, 0), &verts1(k, j - 1, i, 1, 0), &verts1(k, j - 1, i, 2, 0), &verts1(k, j - 1, i, 3, 0), p11);
                    // 0,1,2,5交接边点下标
                    double pp0[3];
                    double pp1[3];
                    double pp2[3];
                    double pp5[3];
                    _get_midpoint(&verts1(k, j - 1, i, 0, 0), &verts1(k, j - 1, i, 2, 0), pp0);
                    _get_midpoint(&verts1(k, j - 1, i, 2, 0), &verts1(k, j - 1, i, 3, 0), pp1);
                    _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 2, 0), pp2);
                    _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 4, 0), pp5);
                    double n0[3];
                    double n1[3];
                    double n4[3];
                    double n8[3];
                    double n11[3];
                    _get_surface_normal(&verts1(k, j, i, 0, 0), pp2, p0, pp5, n0, Axis::XPOSITIVE);
                    _get_surface_normal(&verts1(k, j - 1, i, 2, 0), pp0, p1, pp5, n1, Axis::XPOSITIVE);
                    _get_surface_normal(&verts1(k, j, i, 0, 0), pp5, p4, pp1, n4, Axis::YPOSITIVE);
                    _get_surface_normal(&verts1(k, j, i, 0, 0), pp1, p8, pp2, n8, Axis::ZPOSITIVE);
                    _get_surface_normal(&verts1(k, j - 1, i, 2, 0), pp0, p11, pp1, n11, Axis::ZPOSITIVE);
                    Eigen::Matrix<double, 6, 6> A;
                    A.setZero();
                    A.row(0) << 0, 0, 0, p0[0] - _cx(k, j, i), p0[1] - _cy(k, j, i), p0[2] - _cz(k, j, i);
                    A.row(1) << p1[0] - _cx(k, j - 1, i), p1[1] - _cy(k, j - 1, i), p1[2] - _cz(k, j - 1, i), 0, 0, 0;
                    A.row(2) << _cx(k, j - 1, i) - p4[0], _cy(k, j - 1, i) - p4[1], _cz(k, j - 1, i) - p4[2], p4[0] - _cx(k, j, i), p4[1] - _cy(k, j, i), p4[2] - _cz(k, j, i);
                    Eigen::Matrix3d matK6, matK5;
                    _get_matK(pem1[i6], matK6);
                    _get_matK(pem1[i5], matK5);
                    Eigen::RowVector3d v0, v1, v4, v8, v11;
                    v0 << n0[0], n0[1], n0[2];
                    v1 << n1[0], n1[1], n1[2];
                    v4 << n4[0], n4[1], n4[2];
                    v8 << n8[0], n8[1], n8[2];
                    v11 << n11[0], n11[1], n11[2];
                    A.block(3, 0, 1, 3) = v4 * matK5;
                    A.block(3, 3, 1, 3) = -v4 * matK6;
                    A.block(4, 3, 1, 3) = v8 * matK6;
                    A.block(5, 0, 1, 3) = v11 * matK5;
                    A = A.inverse();
                    Eigen::RowVectorXd r = (v4 - v1) * matK5 * A.topRows(3);
                    B[i5] -= pb * (r[0] + r[1]);
                    C(i5, i5) += -r[1] + r[2];
                    C(i5, i6) += -r[0] - r[2];
                    r = -(v0 + v4) * matK6 * A.bottomRows(3);
                    B[i6] -= pb * (r[0] + r[1]);
                    C(i6, i5) += -r[1] + r[2];
                    C(i6, i6) += -r[0] - r[2];
                }
                else if (i == nx && k == 0)
                {
                    // 两个网格中心下标
                    int i7 = cur - 1, i4 = cur - nx - 1;
                    // 0,1,7,9,10交接面中点下标
                    double p0[3];
                    double p1[3];
                    double p7[3];
                    double p9[3];
                    double p10[3];
                    _get_centroid(&verts1(k, j, i - 1, 1, 0), &verts1(k, j, i - 1, 3, 0), &verts1(k, j, i - 1, 5, 0), &verts1(k, j, i - 1, 7, 0), p0);
                    _get_centroid(&verts1(k, j - 1, i - 1, 1, 0), &verts1(k, j - 1, i - 1, 3, 0), &verts1(k, j - 1, i - 1, 5, 0), &verts1(k, j - 1, i - 1, 7, 0), p1);
                    _get_centroid(&verts1(k, j, i - 1, 1, 0), &verts1(k, j, i - 1, 0, 0), &verts1(k, j, i - 1, 4, 0), &verts1(k, j, i - 1, 5, 0), p7);
                    _get_centroid(&verts1(k, j, i - 1, 0, 0), &verts1(k, j, i - 1, 1, 0), &verts1(k, j, i - 1, 2, 0), &verts1(k, j, i - 1, 3, 0), p9);
                    _get_centroid(&verts1(k, j - 1, i - 1, 0, 0), &verts1(k, j - 1, i - 1, 1, 0), &verts1(k, j - 1, i - 1, 2, 0), &verts1(k, j - 1, i - 1, 3, 0), p10);
                    // 2,3,0,5交接边点下标
                    double pp2[3];
                    double pp3[3];
                    double pp0[3];
                    double pp5[3];
                    _get_midpoint(&verts1(k, j - 1, i - 1, 3, 0), &verts1(k, j - 1, i - 1, 1, 0), pp0);
                    _get_midpoint(&verts1(k, j - 1, i - 1, 2, 0), &verts1(k, j - 1, i - 1, 3, 0), pp3);
                    _get_midpoint(&verts1(k, j, i - 1, 1, 0), &verts1(k, j, i - 1, 3, 0), pp2);
                    _get_midpoint(&verts1(k, j, i - 1, 1, 0), &verts1(k, j, i - 1, 5, 0), pp5);
                    double n0[3];
                    double n1[3];
                    double n7[3];
                    double n9[3];
                    double n10[3];
                    _get_surface_normal(&verts1(k, j, i - 1, 1, 0), pp2, p0, pp5, n0, Axis::XPOSITIVE);
                    _get_surface_normal(&verts1(k, j - 1, i - 1, 3, 0), pp0, p1, pp5, n1, Axis::XPOSITIVE);
                    _get_surface_normal(&verts1(k, j, i - 1, 1, 0), pp5, p7, pp3, n7, Axis::YPOSITIVE);
                    _get_surface_normal(&verts1(k, j, i - 1, 1, 0), pp3, p9, pp2, n9, Axis::ZPOSITIVE);
                    _get_surface_normal(&verts1(k, j - 1, i - 1, 3, 0), pp0, p10, pp3, n10, Axis::ZPOSITIVE);
                    Eigen::Matrix<double, 6, 6> A;
                    A.setZero();
                    A.row(0) << 0, 0, 0, p0[0] - _cx(k, j, i - 1), p0[1] - _cy(k, j, i - 1), p0[2] - _cz(k, j, i - 1);
                    A.row(1) << p1[0] - _cx(k, j - 1, i - 1), p1[1] - _cy(k, j - 1, i - 1), p1[2] - _cz(k, j - 1, i - 1), 0, 0, 0;
                    A.row(2) << _cx(k, j - 1, i - 1) - p7[0], _cy(k, j - 1, i - 1) - p7[1], _cz(k, j - 1, i - 1) - p7[2], p7[0] - _cx(k, j, i - 1), p7[1] - _cy(k, j, i - 1), p7[2] - _cz(k, j, i - 1);
                    Eigen::Matrix3d matK7, matK4;
                    _get_matK(pem1[i7], matK7);
                    _get_matK(pem1[i4], matK4);
                    Eigen::RowVector3d v0, v1, v7, v9, v10;
                    v0 << n0[0], n0[1], n0[2];
                    v1 << n1[0], n1[1], n1[2];
                    v7 << n7[0], n7[1], n7[2];
                    v9 << n9[0], n9[1], n9[2];
                    v10 << n10[0], n10[1], n10[2];
                    A.block(3, 0, 1, 3) = v7 * matK4;
                    A.block(3, 3, 1, 3) = -v7 * matK7;
                    A.block(4, 3, 1, 3) = v9 * matK7;
                    A.block(5, 0, 1, 3) = v10 * matK4;
                    A = A.inverse();
                    Eigen::RowVectorXd r = (v7 + v1 - v10) * matK4 * A.topRows(3);
                    B[i4] -= pb * (r[0] + r[1]);
                    C(i4, i4) += -r[1] + r[2];
                    C(i4, i7) += -r[0] - r[2];
                    r = (v0 - v7 - v9) * matK7 * A.bottomRows(3);
                    B[i7] -= pb * (r[0] + r[1]);
                    C(i7, i4) += -r[1] + r[2];
                    C(i7, i7) += -r[0] - r[2];
                }
                else if (i == 0 && k == nz)
                {
                    // 两个网格中心下标
                    int i2 = cur - nx * ny, i1 = i2 - nx;
                    // 2,3,5,8,11交接面中点下标
                    double p2[3];
                    double p3[3];
                    double p5[3];
                    double p8[3];
                    double p11[3];
                    _get_centroid(&verts1(k - 1, j - 1, i, 0, 0), &verts1(k - 1, j - 1, i, 2, 0), &verts1(k - 1, j - 1, i, 4, 0), &verts1(k - 1, j - 1, i, 6, 0), p2);
                    _get_centroid(&verts1(k - 1, j, i, 0, 0), &verts1(k - 1, j, i, 2, 0), &verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 6, 0), p3);
                    _get_centroid(&verts1(k - 1, j, i, 0, 0), &verts1(k - 1, j, i, 1, 0), &verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 5, 0), p5);
                    _get_centroid(&verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 5, 0), &verts1(k - 1, j, i, 6, 0), &verts1(k - 1, j, i, 7, 0), p8);
                    _get_centroid(&verts1(k - 1, j - 1, i, 4, 0), &verts1(k - 1, j - 1, i, 5, 0), &verts1(k - 1, j - 1, i, 6, 0), &verts1(k - 1, j - 1, i, 7, 0), p11);
                    // 0,1,2,4交接边点下标
                    double pp0[3];
                    double pp1[3];
                    double pp2[3];
                    double pp4[3];
                    _get_midpoint(&verts1(k - 1, j - 1, i, 4, 0), &verts1(k - 1, j - 1, i, 6, 0), pp0);
                    _get_midpoint(&verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 6, 0), pp2);
                    _get_midpoint(&verts1(k - 1, j - 1, i, 6, 0), &verts1(k - 1, j - 1, i, 7, 0), pp1);
                    _get_midpoint(&verts1(k - 1, j, i, 0, 0), &verts1(k - 1, j, i, 4, 0), pp4);
                    double n2[3];
                    double n3[3];
                    double n5[3];
                    double n8[3];
                    double n11[3];
                    _get_surface_normal(&verts1(k - 1, j - 1, i, 6, 0), pp0, p2, pp4, n2, Axis::XPOSITIVE);
                    _get_surface_normal(&verts1(k - 1, j, i, 4, 0), pp2, p3, pp4, n3, Axis::XPOSITIVE);
                    _get_surface_normal(&verts1(k - 1, j, i, 4, 0), pp4, p5, pp1, n5, Axis::YPOSITIVE);
                    _get_surface_normal(&verts1(k - 1, j, i, 4, 0), pp1, p8, pp2, n8, Axis::ZPOSITIVE);
                    _get_surface_normal(&verts1(k - 1, j - 1, i, 6, 0), pp0, p11, pp1, n11, Axis::ZPOSITIVE);
                    Eigen::Matrix<double, 6, 6> A;
                    A.setZero();
                    A.row(0) << 0, 0, 0, p3[0] - _cx(k - 1, j, i), p3[1] - _cy(k - 1, j, i), p3[2] - _cz(k - 1, j, i);
                    A.row(1) << p2[0] - _cx(k - 1, j - 1, i), p2[1] - _cy(k - 1, j - 1, i), p2[2] - _cz(k - 1, j - 1, i), 0, 0, 0;
                    A.row(2) << _cx(k - 1, j - 1, i) - p5[0], _cy(k - 1, j - 1, i) - p5[1], _cz(k - 1, j - 1, i) - p5[2], p5[0] - _cx(k - 1, j, i), p5[1] - _cy(k - 1, j, i), p5[2] - _cz(k - 1, j, i);
                    Eigen::Matrix3d matK2, matK1;
                    _get_matK(pem1[i2], matK2);
                    _get_matK(pem1[i1], matK1);
                    Eigen::RowVector3d v2, v3, v5, v8, v11;
                    v2 << n2[0], n2[1], n2[2];
                    v3 << n3[0], n3[1], n3[2];
                    v5 << n5[0], n5[1], n5[2];
                    v8 << n8[0], n8[1], n8[2];
                    v11 << n11[0], n11[1], n11[2];
                    A.block(3, 0, 1, 3) = v5 * matK1;
                    A.block(3, 3, 1, 3) = -v5 * matK2;
                    A.block(4, 3, 1, 3) = v8 * matK2;
                    A.block(5, 0, 1, 3) = v11 * matK1;
                    A = A.inverse();
                    Eigen::RowVectorXd r = (v5 - v2) * matK1 * A.topRows(3);
                    B[i1] -= pb * (r[0] + r[1]);
                    C(i1, i1) += -r[1] + r[2];
                    C(i1, i2) += -r[0] - r[2];
                    r = -(v3 + v5) * matK2 * A.bottomRows(3);
                    B[i2] -= pb * (r[0] + r[1]);
                    C(i2, i1) += -r[1] + r[2];
                    C(i2, i2) += -r[0] - r[2];
                }
                else if (i == nx && k == nz)
                {
                    int i3 = cur - nx * ny - 1, i0 = i3 - nx;
                    // 2,3,6,9,10交接面中点下标
                    double p2[3];
                    double p3[3];
                    double p6[3];
                    double p9[3];
                    double p10[3];
                    _get_centroid(&verts1(k - 1, j - 1, i - 1, 1, 0), &verts1(k - 1, j - 1, i - 1, 3, 0), &verts1(k - 1, j - 1, i - 1, 5, 0), &verts1(k - 1, j - 1, i - 1, 7, 0), p2);
                    _get_centroid(&verts1(k - 1, j, i - 1, 1, 0), &verts1(k - 1, j, i - 1, 3, 0), &verts1(k - 1, j, i - 1, 5, 0), &verts1(k - 1, j, i - 1, 7, 0), p3);
                    _get_centroid(&verts1(k - 1, j, i - 1, 1, 0), &verts1(k - 1, j, i - 1, 0, 0), &verts1(k - 1, j, i - 1, 4, 0), &verts1(k - 1, j, i - 1, 5, 0), p6);
                    _get_centroid(&verts1(k - 1, j, i - 1, 4, 0), &verts1(k - 1, j, i - 1, 5, 0), &verts1(k - 1, j, i - 1, 6, 0), &verts1(k - 1, j, i - 1, 7, 0), p9);
                    _get_centroid(&verts1(k - 1, j - 1, i - 1, 4, 0), &verts1(k - 1, j - 1, i - 1, 5, 0), &verts1(k - 1, j - 1, i - 1, 6, 0), &verts1(k - 1, j - 1, i - 1, 7, 0), p10);
                    // 0,2,3,4交接边点下标
                    double pp0[3];
                    double pp2[3];
                    double pp3[3];
                    double pp4[3];
                    _get_midpoint(&verts1(k - 1, j - 1, i - 1, 5, 0), &verts1(k - 1, j - 1, i - 1, 7, 0), pp0);
                    _get_midpoint(&verts1(k - 1, j, i - 1, 5, 0), &verts1(k - 1, j, i - 1, 7, 0), pp2);
                    _get_midpoint(&verts1(k - 1, j, i - 1, 4, 0), &verts1(k - 1, j, i - 1, 5, 0), pp3);
                    _get_midpoint(&verts1(k - 1, j - 1, i - 1, 7, 0), &verts1(k - 1, j - 1, i - 1, 3, 0), pp4);
                    double n2[3];
                    double n3[3];
                    double n6[3];
                    double n9[3];
                    double n10[3];
                    _get_surface_normal(&verts1(k - 1, j - 1, i - 1, 7, 0), pp0, p2, pp4, n2, Axis::XPOSITIVE);
                    _get_surface_normal(&verts1(k - 1, j, i - 1, 5, 0), pp2, p3, pp4, n3, Axis::XPOSITIVE);
                    _get_surface_normal(&verts1(k - 1, j, i - 1, 5, 0), pp4, p6, pp3, n6, Axis::YPOSITIVE);
                    _get_surface_normal(&verts1(k - 1, j, i - 1, 5, 0), pp3, p9, pp2, n9, Axis::ZPOSITIVE);
                    _get_surface_normal(&verts1(k - 1, j - 1, i - 1, 7, 0), pp0, p10, pp3, n10, Axis::ZPOSITIVE);
                    Eigen::Matrix<double, 6, 6> A;
                    A.setZero();
                    A.row(0) << 0, 0, 0, p3[0] - _cx(k - 1, j, i - 1), p3[1] - _cy(k - 1, j, i - 1), p3[2] - _cz(k - 1, j, i - 1);
                    A.row(1) << p2[0] - _cx(k - 1, j - 1, i - 1), p2[1] - _cy(k - 1, j - 1, i - 1), p2[2] - _cz(k - 1, j - 1, i - 1), 0, 0, 0;
                    A.row(2) << _cx(k - 1, j - 1, i - 1) - p6[0], _cy(k - 1, j - 1, i - 1) - p6[1], _cz(k - 1, j - 1, i - 1) - p6[2], p6[0] - _cx(k - 1, j, i - 1), p6[1] - _cy(k - 1, j, i - 1), p6[2] - _cz(k - 1, j, i - 1);
                    Eigen::Matrix3d matK3, matK0;
                    _get_matK(pem1[i3], matK3);
                    _get_matK(pem1[i0], matK0);
                    Eigen::RowVector3d v2, v3, v6, v9, v10;
                    v2 << n2[0], n2[1], n2[2];
                    v3 << n3[0], n3[1], n3[2];
                    v6 << n6[0], n6[1], n6[2];
                    v9 << n9[0], n9[1], n9[2];
                    v10 << n10[0], n10[1], n10[2];
                    A.block(3, 0, 1, 3) = v6 * matK0;
                    A.block(3, 3, 1, 3) = -v6 * matK3;
                    A.block(4, 3, 1, 3) = v9 * matK3;
                    A.block(5, 0, 1, 3) = v10 * matK0;
                    A = A.inverse();
                    Eigen::RowVectorXd r = (v6 + v2) * matK0 * A.topRows(3);
                    B[i0] -= pb * (r[0] + r[1]);
                    C(i0, i0) += -r[1] + r[2];
                    C(i0, i3) += -r[0] - r[2];
                    r = (v3 - v6) * matK3 * A.bottomRows(3);
                    B[i3] -= pb * (r[0] + r[1]);
                    C(i3, i0) += -r[1] + r[2];
                    C(i3, i3) += -r[0] - r[2];
                }
                else if (i == 0 && j == 0)
                {
                    int i6 = cur, i2 = cur - nx * ny;
                    // 0,3,4,5,8交接面中点下标
                    double p0[3];
                    double p3[3];
                    double p4[3];
                    double p5[3];
                    double p8[3];
                    _get_centroid(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 2, 0), &verts1(k, j, i, 4, 0), &verts1(k, j, i, 6, 0), p0);
                    _get_centroid(&verts1(k - 1, j, i, 0, 0), &verts1(k - 1, j, i, 2, 0), &verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 6, 0), p3);
                    _get_centroid(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 1, 0), &verts1(k, j, i, 4, 0), &verts1(k, j, i, 5, 0), p4);
                    _get_centroid(&verts1(k - 1, j, i, 0, 0), &verts1(k - 1, j, i, 1, 0), &verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 5, 0), p5);
                    _get_centroid(&verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 5, 0), &verts1(k - 1, j, i, 6, 0), &verts1(k - 1, j, i, 7, 0), p8);
                    // 1,2,4,5交接边点下标
                    double pp1[3];
                    double pp2[3];
                    double pp4[3];
                    double pp5[3];
                    _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 1, 0), pp1);
                    _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 2, 0), pp2);
                    _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 4, 0), pp5);
                    _get_midpoint(&verts1(k - 1, j, i, 0, 0), &verts1(k - 1, j, i, 4, 0), pp4);
                    double n0[3];
                    double n3[3];
                    double n4[3];
                    double n5[3];
                    double n8[3];
                    _get_surface_normal(&verts1(k, j, i, 0, 0), pp2, p0, pp5, n0, Axis::XPOSITIVE);
                    _get_surface_normal(&verts1(k - 1, j, i, 4, 0), pp2, p3, pp4, n3, Axis::XPOSITIVE);
                    _get_surface_normal(&verts1(k, j, i, 0, 0), pp5, p4, pp1, n4, Axis::YPOSITIVE);
                    _get_surface_normal(&verts1(k - 1, j, i, 4, 0), pp4, p5, pp1, n5, Axis::YPOSITIVE);
                    _get_surface_normal(&verts1(k - 1, j, i, 4, 0), pp1, p8, pp2, n8, Axis::ZPOSITIVE);
                    Eigen::Matrix<double, 6, 6> A;
                    A.setZero();
                    A.row(0) << 0, 0, 0, p0[0] - _cx(k, j, i), p0[1] - _cy(k, j, i), p0[2] - _cz(k, j, i);
                    A.row(1) << p3[0] - _cx(k - 1, j, i), p3[1] - _cy(k - 1, j, i), p3[2] - _cz(k - 1, j, i), 0, 0, 0;
                    A.row(2) << _cx(k - 1, j, i) - p8[0], _cy(k - 1, j, i) - p8[1], _cz(k - 1, j, i) - p8[2], p8[0] - _cx(k, j, i), p8[1] - _cy(k, j, i), p8[2] - _cz(k, j, i);
                    Eigen::Matrix3d matK6, matK2;
                    _get_matK(pem1[i6], matK6);
                    _get_matK(pem1[i2], matK2);
                    Eigen::RowVector3d v0, v3, v4, v5, v8;
                    v0 << n0[0], n0[1], n0[2];
                    v3 << n3[0], n3[1], n3[2];
                    v4 << n4[0], n4[1], n4[2];
                    v5 << n5[0], n5[1], n5[2];
                    v8 << n8[0], n8[1], n8[2];
                    A.block(3, 0, 1, 3) = v8 * matK2;
                    A.block(3, 3, 1, 3) = -v8 * matK6;
                    A.block(4, 3, 1, 3) = v4 * matK6;
                    A.block(5, 0, 1, 3) = v5 * matK2;
                    A = A.inverse();
                    Eigen::RowVectorXd r = (v8 - v3) * matK2 * A.topRows(3);
                    B[i2] -= pb * (r[0] + r[1]);
                    C(i2, i2) += -r[1] + r[2];
                    C(i2, i6) += -r[0] - r[2];
                    r = -(v0 + v8) * matK6 * A.bottomRows(3);
                    B[i6] -= pb * (r[0] + r[1]);
                    C(i6, i2) += -r[1] + r[2];
                    C(i6, i6) += -r[0] - r[2];
                    // To be implemented
                }
                else if (i == nx && j == 0)
                {
                    int i7 = cur - 1, i3 = i7 - nx * ny;
                    // 0,3,6,7,9交接面中点下标
                    double p0[3];
                    double p3[3];
                    double p6[3];
                    double p7[3];
                    double p9[3];
                    _get_centroid(&verts1(k, j, i - 1, 1, 0), &verts1(k, j, i - 1, 3, 0), &verts1(k, j, i - 1, 5, 0), &verts1(k, j, i - 1, 7, 0), p0);
                    _get_centroid(&verts1(k - 1, j, i - 1, 1, 0), &verts1(k - 1, j, i - 1, 3, 0), &verts1(k - 1, j, i - 1, 5, 0), &verts1(k - 1, j, i - 1, 7, 0), p3);
                    _get_centroid(&verts1(k, j, i - 1, 1, 0), &verts1(k, j, i - 1, 0, 0), &verts1(k, j, i - 1, 4, 0), &verts1(k, j, i - 1, 5, 0), p7);
                    _get_centroid(&verts1(k - 1, j, i - 1, 1, 0), &verts1(k - 1, j, i - 1, 0, 0), &verts1(k - 1, j, i - 1, 4, 0), &verts1(k - 1, j, i - 1, 5, 0), p6);
                    _get_centroid(&verts1(k - 1, j, i - 1, 4, 0), &verts1(k - 1, j, i - 1, 5, 0), &verts1(k - 1, j, i - 1, 6, 0), &verts1(k - 1, j, i - 1, 7, 0), p9);
                    // 2,3,4,5交接边点下标
                    double pp2[3];
                    double pp3[3];
                    double pp4[3];
                    double pp5[3];
                    _get_midpoint(&verts1(k, j, i - 1, 1, 0), &verts1(k, j, i - 1, 3, 0), pp2);
                    _get_midpoint(&verts1(k, j, i - 1, 1, 0), &verts1(k, j, i - 1, 0, 0), pp3);
                    _get_midpoint(&verts1(k, j, i - 1, 1, 0), &verts1(k, j, i - 1, 5, 0), pp5);
                    _get_midpoint(&verts1(k - 1, j, i - 1, 1, 0), &verts1(k - 1, j, i - 1, 5, 0), pp4);
                    double n0[3];
                    double n3[3];
                    double n6[3];
                    double n7[3];
                    double n9[3];
                    _get_surface_normal(&verts1(k, j, i - 1, 1, 0), pp2, p0, pp5, n0, Axis::XPOSITIVE);
                    _get_surface_normal(&verts1(k - 1, j, i - 1, 5, 0), pp4, p3, pp2, n3, Axis::XPOSITIVE);
                    _get_surface_normal(&verts1(k, j, i - 1, 1, 0), pp3, p7, pp5, n7, Axis::YPOSITIVE);
                    _get_surface_normal(&verts1(k - 1, j, i - 1, 5, 0), pp4, p6, pp3, n6, Axis::YPOSITIVE);
                    _get_surface_normal(&verts1(k - 1, j, i - 1, 5, 0), pp2, p9, pp3, n9, Axis::ZPOSITIVE);
                    Eigen::Matrix<double, 6, 6> A;
                    A.setZero();
                    A.row(0) << 0, 0, 0, p0[0] - _cx(k, j, i - 1), p0[1] - _cy(k, j, i - 1), p0[2] - _cz(k, j, i - 1);
                    A.row(1) << p3[0] - _cx(k - 1, j, i - 1), p3[1] - _cy(k - 1, j, i - 1), p3[2] - _cz(k - 1, j, i - 1), 0, 0, 0;
                    A.row(2) << _cx(k - 1, j, i - 1) - p9[0], _cy(k - 1, j, i - 1) - p9[1], _cz(k - 1, j, i - 1) - p9[2], p9[0] - _cx(k, j, i - 1), p9[1] - _cy(k, j, i - 1), p9[2] - _cz(k, j, i - 1);
                    Eigen::Matrix3d matK7, matK3;
                    _get_matK(pem1[i7], matK7);
                    _get_matK(pem1[i3], matK3);
                    Eigen::RowVector3d v0, v3, v6, v7, v9;
                    v0 << n0[0], n0[1], n0[2];
                    v3 << n3[0], n3[1], n3[2];
                    v6 << n6[0], n6[1], n6[2];
                    v7 << n7[0], n7[1], n7[2];
                    v9 << n9[0], n9[1], n9[2];
                    A.block(3, 0, 1, 3) = v9 * matK3;
                    A.block(3, 3, 1, 3) = -v9 * matK7;
                    A.block(4, 3, 1, 3) = v7 * matK7;
                    A.block(5, 0, 1, 3) = v6 * matK3;
                    A = A.inverse();
                    Eigen::RowVectorXd r = (v9 + v3) * matK3 * A.topRows(3);
                    B[i3] -= pb * (r[0] + r[1]);
                    C(i3, i3) += -r[1] + r[2];
                    C(i3, i7) += -r[0] - r[2];
                    r = (v0 - v9) * matK7 * A.bottomRows(3);
                    B[i7] -= pb * (r[0] + r[1]);
                    C(i7, i3) += -r[1] + r[2];
                    C(i7, i7) += -r[0] - r[2];
                    //  To be implemented
                }
                else if (i == 0 && j == ny)
                {
                    int i5 = cur - nx, i1 = i5 - nx * ny;
                    // 1,2,4,5,11交接面中点下标
                    double p1[3];
                    double p2[3];
                    double p4[3];
                    double p5[3];
                    double p11[3];
                    _get_centroid(&verts1(k, j - 1, i, 0, 0), &verts1(k, j - 1, i, 2, 0), &verts1(k, j - 1, i, 4, 0), &verts1(k, j - 1, i, 6, 0), p1);
                    _get_centroid(&verts1(k - 1, j - 1, i, 0, 0), &verts1(k - 1, j - 1, i, 2, 0), &verts1(k - 1, j - 1, i, 4, 0), &verts1(k - 1, j - 1, i, 6, 0), p2);
                    _get_centroid(&verts1(k, j - 1, i, 2, 0), &verts1(k, j - 1, i, 3, 0), &verts1(k, j - 1, i, 6, 0), &verts1(k, j - 1, i, 7, 0), p4);
                    _get_centroid(&verts1(k - 1, j - 1, i, 2, 0), &verts1(k - 1, j - 1, i, 3, 0), &verts1(k - 1, j - 1, i, 6, 0), &verts1(k - 1, j - 1, i, 7, 0), p5);
                    _get_centroid(&verts1(k - 1, j - 1, i, 4, 0), &verts1(k - 1, j - 1, i, 5, 0), &verts1(k - 1, j - 1, i, 6, 0), &verts1(k - 1, j - 1, i, 7, 0), p11);
                    // 0,1,4,5交接边点下标
                    double pp0[3];
                    double pp1[3];
                    double pp4[3];
                    double pp5[3];
                    _get_midpoint(&verts1(k, j - 1, i, 0, 0), &verts1(k, j - 1, i, 2, 0), pp0);
                    _get_midpoint(&verts1(k, j - 1, i, 3, 0), &verts1(k, j - 1, i, 2, 0), pp1);
                    _get_midpoint(&verts1(k, j - 1, i, 6, 0), &verts1(k, j - 1, i, 2, 0), pp5);
                    _get_midpoint(&verts1(k - 1, j - 1, i, 2, 0), &verts1(k - 1, j - 1, i, 6, 0), pp4);
                    double n1[3];
                    double n2[3];
                    double n4[3];
                    double n5[3];
                    double n11[3];
                    _get_surface_normal(&verts1(k, j - 1, i, 2, 0), pp0, p1, pp5, n1, Axis::XPOSITIVE);
                    _get_surface_normal(&verts1(k - 1, j - 1, i, 6, 0), pp4, p2, pp0, n2, Axis::XPOSITIVE);
                    _get_surface_normal(&verts1(k, j - 1, i, 2, 0), pp1, p4, pp5, n4, Axis::YPOSITIVE);
                    _get_surface_normal(&verts1(k - 1, j - 1, i, 6, 0), pp4, p5, pp1, n5, Axis::YPOSITIVE);
                    _get_surface_normal(&verts1(k - 1, j - 1, i, 6, 0), pp0, p11, pp1, n11, Axis::ZPOSITIVE);
                    Eigen::Matrix<double, 6, 6> A;
                    A.setZero();
                    A.row(0) << 0, 0, 0, p1[0] - _cx(k, j - 1, i), p1[1] - _cy(k, j - 1, i), p1[2] - _cz(k, j - 1, i);
                    A.row(1) << p2[0] - _cx(k - 1, j - 1, i), p2[1] - _cy(k - 1, j - 1, i), p2[2] - _cz(k - 1, j - 1, i), 0, 0, 0;
                    A.row(2) << _cx(k - 1, j - 1, i) - p11[0], _cy(k - 1, j - 1, i) - p11[1], _cz(k - 1, j - 1, i) - p11[2], p11[0] - _cx(k, j - 1, i), p11[1] - _cy(k, j - 1, i), p11[2] - _cz(k, j - 1, i);
                    Eigen::Matrix3d matK5, matK1;
                    _get_matK(pem1[i5], matK5);
                    _get_matK(pem1[i1], matK1);
                    Eigen::RowVector3d v1, v2, v4, v5, v11;
                    v1 << n1[0], n1[1], n1[2];
                    v2 << n2[0], n2[1], n2[2];
                    v4 << n4[0], n4[1], n4[2];
                    v5 << n5[0], n5[1], n5[2];
                    v11 << n11[0], n11[1], n11[2];
                    A.block(3, 0, 1, 3) = v11 * matK1;
                    A.block(3, 3, 1, 3) = -v11 * matK5;
                    A.block(4, 3, 1, 3) = v4 * matK5;
                    A.block(5, 0, 1, 3) = v5 * matK1;
                    A = A.inverse();
                    Eigen::RowVectorXd r = (v11 - v2) * matK1 * A.topRows(3);
                    B[i1] -= pb * (r[0] + r[1]);
                    C(i1, i1) += -r[1] + r[2];
                    C(i1, i5) += -r[0] - r[2];
                    r = -(v1 + v11) * matK5 * A.bottomRows(3);
                    B[i5] -= pb * (r[0] + r[1]);
                    C(i5, i1) += -r[1] + r[2];
                    C(i5, i5) += -r[0] - r[2];
                    // To be implemented
                }
                else if (i == nx && j == ny)
                {
                    int i4 = cur - 1 - nx, i0 = i4 - nx * ny;
                    // 1,2,6,7,10交接面中点下标
                    double p1[3];
                    double p2[3];
                    double p6[3];
                    double p7[3];
                    double p10[3];
                    _get_centroid(&verts1(k, j - 1, i - 1, 1, 0), &verts1(k, j - 1, i - 1, 3, 0), &verts1(k, j - 1, i - 1, 5, 0), &verts1(k, j - 1, i - 1, 7, 0), p1);
                    _get_centroid(&verts1(k - 1, j - 1, i - 1, 1, 0), &verts1(k - 1, j - 1, i - 1, 3, 0), &verts1(k - 1, j - 1, i - 1, 5, 0), &verts1(k - 1, j - 1, i - 1, 7, 0), p2);
                    _get_centroid(&verts1(k, j - 1, i - 1, 3, 0), &verts1(k, j - 1, i - 1, 2, 0), &verts1(k, j - 1, i - 1, 6, 0), &verts1(k, j - 1, i - 1, 7, 0), p7);
                    _get_centroid(&verts1(k - 1, j - 1, i - 1, 3, 0), &verts1(k - 1, j - 1, i - 1, 2, 0), &verts1(k - 1, j - 1, i - 1, 6, 0), &verts1(k - 1, j - 1, i - 1, 7, 0), p6);
                    _get_centroid(&verts1(k - 1, j - 1, i - 1, 5, 0), &verts1(k - 1, j - 1, i - 1, 4, 0), &verts1(k - 1, j - 1, i - 1, 6, 0), &verts1(k - 1, j - 1, i - 1, 7, 0), p10);
                    // 0,3,4,5交接边点下标
                    double pp0[3];
                    double pp3[3];
                    double pp4[3];
                    double pp5[3];
                    _get_midpoint(&verts1(k, j - 1, i - 1, 1, 0), &verts1(k, j - 1, i - 1, 3, 0), pp0);
                    _get_midpoint(&verts1(k, j - 1, i - 1, 2, 0), &verts1(k, j - 1, i - 1, 3, 0), pp3);
                    _get_midpoint(&verts1(k, j - 1, i - 1, 7, 0), &verts1(k, j - 1, i - 1, 3, 0), pp5);
                    _get_midpoint(&verts1(k - 1, j - 1, i - 1, 3, 0), &verts1(k - 1, j - 1, i - 1, 7, 0), pp4);
                    double n1[3];
                    double n2[3];
                    double n6[3];
                    double n7[3];
                    double n10[3];
                    _get_surface_normal(&verts1(k, j - 1, i - 1, 3, 0), pp0, p1, pp5, n1, Axis::XPOSITIVE);
                    _get_surface_normal(&verts1(k - 1, j - 1, i - 1, 7, 0), pp4, p2, pp0, n2, Axis::XPOSITIVE);
                    _get_surface_normal(&verts1(k, j - 1, i - 1, 3, 0), pp3, p7, pp5, n7, Axis::YPOSITIVE);
                    _get_surface_normal(&verts1(k - 1, j - 1, i - 1, 7, 0), pp4, p6, pp3, n6, Axis::YPOSITIVE);
                    _get_surface_normal(&verts1(k - 1, j - 1, i - 1, 7, 0), pp0, p10, pp3, n10, Axis::ZPOSITIVE);
                    Eigen::Matrix<double, 6, 6> A;
                    A.setZero();
                    A.row(0) << 0, 0, 0, p1[0] - _cx(k, j - 1, i - 1), p1[1] - _cy(k, j - 1, i - 1), p1[2] - _cz(k, j - 1, i - 1);
                    A.row(1) << p2[0] - _cx(k - 1, j - 1, i - 1), p2[1] - _cy(k - 1, j - 1, i - 1), p2[2] - _cz(k - 1, j - 1, i - 1), 0, 0, 0;
                    A.row(2) << _cx(k - 1, j - 1, i - 1) - p10[0], _cy(k - 1, j - 1, i - 1) - p10[1], _cz(k - 1, j - 1, i - 1) - p10[2], p10[0] - _cx(k, j - 1, i - 1), p10[1] - _cy(k, j - 1, i - 1), p10[2] - _cz(k, j - 1, i - 1);
                    Eigen::Matrix3d matK4, matK0;
                    _get_matK(pem1[i4], matK4);
                    _get_matK(pem1[i0], matK0);
                    Eigen::RowVector3d v1, v2, v6, v7, v10;
                    v1 << n1[0], n1[1], n1[2];
                    v2 << n2[0], n2[1], n2[2];
                    v6 << n6[0], n6[1], n6[2];
                    v7 << n7[0], n7[1], n7[2];
                    v10 << n10[0], n10[1], n10[2];
                    A.block(3, 0, 1, 3) = v10 * matK0;
                    A.block(3, 3, 1, 3) = -v10 * matK4;
                    A.block(4, 3, 1, 3) = v7 * matK4;
                    A.block(5, 0, 1, 3) = v6 * matK0;
                    A = A.inverse();
                    Eigen::RowVectorXd r = (v10 + v2) * matK0 * A.topRows(3);
                    B[i0] -= pb * (r[0] + r[1]);
                    C(i0, i0) += -r[1] + r[2];
                    C(i0, i4) += -r[0] - r[2];
                    r = (v1 - v10) * matK4 * A.bottomRows(3);
                    B[i4] -= pb * (r[0] + r[1]);
                    C(i4, i0) += -r[1] + r[2];
                    C(i4, i4) += -r[0] - r[2];
                    // To be implemented
                }
                // 六个表面其余点
                else if (i == 0)
                {
                    int i6 = cur, i5 = i6 - nx, i2 = i6 - nx * ny, i1 = i5 - nx * ny;
                    // 0,1,2,3,4,5,8,11交接面中点下标
                    double p0[3];
                    double p1[3];
                    double p2[3];
                    double p3[3];
                    double p4[3];
                    double p5[3];
                    double p8[3];
                    double p11[3];
                    _get_centroid(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 1, 0), &verts1(k, j, i, 2, 0), &verts1(k, j, i, 3, 0), p8);
                    _get_centroid(&verts1(k, j - 1, i, 0, 0), &verts1(k, j - 1, i, 1, 0), &verts1(k, j - 1, i, 2, 0), &verts1(k, j - 1, i, 3, 0), p11);
                    _get_centroid(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 1, 0), &verts1(k, j, i, 4, 0), &verts1(k, j, i, 5, 0), p4);
                    _get_centroid(&verts1(k - 1, j, i, 0, 0), &verts1(k - 1, j, i, 1, 0), &verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 5, 0), p5);
                    _get_centroid(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 2, 0), &verts1(k, j, i, 4, 0), &verts1(k, j, i, 6, 0), p0);
                    _get_centroid(&verts1(k - 1, j, i, 0, 0), &verts1(k - 1, j, i, 2, 0), &verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 6, 0), p3);
                    _get_centroid(&verts1(k, j - 1, i, 0, 0), &verts1(k, j - 1, i, 2, 0), &verts1(k, j - 1, i, 4, 0), &verts1(k, j - 1, i, 6, 0), p1);
                    _get_centroid(&verts1(k - 1, j - 1, i, 0, 0), &verts1(k - 1, j - 1, i, 2, 0), &verts1(k - 1, j - 1, i, 4, 0), &verts1(k - 1, j - 1, i, 6, 0), p2);
                    // 交接边点下标0,1,2,4,5
                    double pp0[3];
                    double pp1[3];
                    double pp2[3];
                    double pp4[3];
                    double pp5[3];
                    _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 1, 0), pp1);
                    _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 2, 0), pp2);
                    _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 4, 0), pp5);
                    _get_midpoint(&verts1(k - 1, j, i, 0, 0), &verts1(k - 1, j, i, 4, 0), pp4);
                    _get_midpoint(&verts1(k, j - 1, i, 0, 0), &verts1(k, j - 1, i, 2, 0), pp0);
                    double n0[3];
                    double n1[3];
                    double n2[3];
                    double n3[3];
                    double n4[3];
                    double n5[3];
                    double n8[3];
                    double n11[3];
                    _get_surface_normal(&verts1(k, j, i, 0, 0), pp2, p0, pp5, n0, Axis::XPOSITIVE);
                    _get_surface_normal(&verts1(k, j, i, 0, 0), pp0, p1, pp5, n1, Axis::XPOSITIVE);
                    _get_surface_normal(&verts1(k, j, i, 0, 0), pp0, p2, pp4, n2, Axis::XPOSITIVE);
                    _get_surface_normal(&verts1(k, j, i, 0, 0), pp4, p3, pp2, n3, Axis::XPOSITIVE);
                    _get_surface_normal(&verts1(k, j, i, 0, 0), pp1, p4, pp5, n4, Axis::YPOSITIVE);
                    _get_surface_normal(&verts1(k, j, i, 0, 0), pp4, p5, pp1, n5, Axis::YPOSITIVE);
                    _get_surface_normal(&verts1(k, j, i, 0, 0), pp2, p8, pp1, n8, Axis::ZPOSITIVE);
                    _get_surface_normal(&verts1(k, j, i, 0, 0), pp0, p11, pp1, n11, Axis::ZPOSITIVE);
                    Eigen::Matrix<double, 12, 12> A;
                    A.setZero();
                    A.block(0, 9, 1, 3) << p0[0] - _cx(k, j, i), p0[1] - _cy(k, j, i), p0[2] - _cz(k, j, i);
                    A.block(1, 6, 1, 3) << p1[0] - _cx(k, j - 1, i), p1[1] - _cy(k, j - 1, i), p1[2] - _cz(k, j - 1, i);
                    A.block(2, 0, 1, 3) << p2[0] - _cx(k - 1, j - 1, i), p2[1] - _cy(k - 1, j - 1, i), p2[2] - _cz(k - 1, j - 1, i);
                    A.block(3, 3, 1, 3) << p3[0] - _cx(k - 1, j, i), p3[1] - _cy(k - 1, j, i), p3[2] - _cz(k - 1, j, i);
                    Eigen::Matrix3d matK6, matK5, matK2, matK1;
                    _get_matK(pem1[i6], matK6);
                    _get_matK(pem1[i5], matK5);
                    _get_matK(pem1[i2], matK2);
                    _get_matK(pem1[i1], matK1);
                    Eigen::RowVector3d v0, v1, v2, v3, v4, v5, v8, v11;
                    v0 << n0[0], n0[1], n0[2];
                    v1 << n1[0], n1[1], n1[2];
                    v2 << n2[0], n2[1], n2[2];
                    v3 << n3[0], n3[1], n3[2];
                    v4 << n4[0], n4[1], n4[2];
                    v5 << n5[0], n5[1], n5[2];
                    v8 << n8[0], n8[1], n8[2];
                    v11 << n11[0], n11[1], n11[2];
                    A.block(4, 6, 1, 6) << _cx(k, j - 1, i) - p4[0], _cy(k, j - 1, i) - p4[1], _cz(k, j - 1, i) - p4[2], p4[0] - _cx(k, j, i), p4[1] - _cy(k, j, i), p4[2] - _cz(k, j, i);
                    A.block(5, 6, 1, 3) = v4 * matK5;
                    A.block(5, 9, 1, 3) = -v4 * matK6;
                    A.block(6, 0, 1, 6) << _cx(k - 1, j - 1, i) - p5[0], _cy(k - 1, j - 1, i) - p5[1], _cz(k - 1, j - 1, i) - p5[2], p5[0] - _cx(k - 1, j, i), p5[1] - _cy(k - 1, j, i), p5[2] - _cz(k - 1, j, i);
                    A.block(7, 0, 1, 3) = v5 * matK1;
                    A.block(7, 3, 1, 3) = -v5 * matK2;
                    A.block(8, 3, 1, 3) << _cx(k - 1, j, i) - p8[0], _cy(k - 1, j, i) - p8[1], _cz(k - 1, j, i) - p8[2];
                    A.block(8, 9, 1, 3) << p8[0] - _cx(k, j, i), p8[1] - _cy(k, j, i), p8[2] - _cz(k, j, i);
                    A.block(9, 3, 1, 3) = v8 * matK2;
                    A.block(9, 9, 1, 3) = -v8 * matK6;
                    A.block(10, 0, 1, 3) << _cx(k - 1, j - 1, i) - p11[0], _cy(k - 1, j - 1, i) - p11[1], _cz(k - 1, j - 1, i) - p11[2];
                    A.block(10, 6, 1, 3) << p11[0] - _cx(k, j - 1, i), p11[1] - _cy(k, j - 1, i), p11[2] - _cz(k, j - 1, i);
                    A.block(11, 0, 1, 3) = v11 * matK1;
                    A.block(11, 6, 1, 3) = -v11 * matK5;
                    A = A.inverse();
                    Eigen::RowVectorXd r = (-v2 + v5 + v11) * matK1 * A.topRows(3);
                    B[i1] -= pb * (r[0] + r[1] + r[2] + r[3]);
                    C(i1, i1) += -r[2] + r[6] + r[10];
                    C(i1, i2) += -r[3] - r[6] + r[8];
                    C(i1, i5) += -r[1] + r[4] - r[10];
                    C(i1, i6) += -r[0] - r[4] - r[8];
                    r = (-v3 - v5 + v8) * matK2 * A.block(3, 0, 3, 12);
                    B[i2] -= pb * (r[0] + r[1] + r[2] + r[3]);
                    C(i2, i1) += -r[2] + r[6] + r[10];
                    C(i2, i2) += -r[3] - r[6] + r[8];
                    C(i2, i5) += -r[1] + r[4] - r[10];
                    C(i2, i6) += -r[0] - r[4] - r[8];
                    r = (-v1 + v4 - v11) * matK5 * A.block(6, 0, 3, 12);
                    B[i5] -= pb * (r[0] + r[1] + r[2] + r[3]);
                    C(i5, i1) += -r[2] + r[6] + r[10];
                    C(i5, i2) += -r[3] - r[6] + r[8];
                    C(i5, i5) += -r[1] + r[4] - r[10];
                    C(i5, i6) += -r[0] - r[4] - r[8];
                    r = -(v0 + v4 + v8) * matK6 * A.bottomRows(3);
                    B[i6] -= pb * (r[0] + r[1] + r[2] + r[3]);
                    C(i6, i1) += -r[2] + r[6] + r[10];
                    C(i6, i2) += -r[3] - r[6] + r[8];
                    C(i6, i5) += -r[1] + r[4] - r[10];
                    C(i6, i6) += -r[0] - r[4] - r[8];
                    //  To be implemented
                }
                else if (i == nx)
                {
                    int i7 = cur - 1, i4 = i7 - nx, i3 = i7 - nx * ny, i0 = i4 - nx * ny;
                    // 0,1,2,3,6,7,9,10交接面中点下标
                    double p0[3];
                    double p1[3];
                    double p2[3];
                    double p3[3];
                    double p6[3];
                    double p7[3];
                    double p9[3];
                    double p10[3];
                    _get_centroid(&verts1(k, j, i - 1, 1, 0), &verts1(k, j, i - 1, 3, 0), &verts1(k, j, i - 1, 5, 0), &verts1(k, j, i - 1, 7, 0), p0);
                    _get_centroid(&verts1(k, j - 1, i - 1, 1, 0), &verts1(k, j - 1, i - 1, 3, 0), &verts1(k, j - 1, i - 1, 5, 0), &verts1(k, j - 1, i - 1, 7, 0), p1);
                    _get_centroid(&verts1(k - 1, j - 1, i - 1, 1, 0), &verts1(k - 1, j - 1, i - 1, 3, 0), &verts1(k - 1, j - 1, i - 1, 5, 0), &verts1(k - 1, j - 1, i - 1, 7, 0), p2);
                    _get_centroid(&verts1(k - 1, j, i - 1, 1, 0), &verts1(k - 1, j, i - 1, 3, 0), &verts1(k - 1, j, i - 1, 5, 0), &verts1(k - 1, j, i - 1, 7, 0), p3);
                    _get_centroid(&verts1(k - 1, j, i - 1, 1, 0), &verts1(k - 1, j, i - 1, 0, 0), &verts1(k - 1, j, i - 1, 5, 0), &verts1(k - 1, j, i - 1, 4, 0), p6);
                    _get_centroid(&verts1(k, j, i - 1, 1, 0), &verts1(k, j, i - 1, 0, 0), &verts1(k, j, i - 1, 5, 0), &verts1(k, j, i - 1, 4, 0), p7);
                    _get_centroid(&verts1(k, j, i - 1, 1, 0), &verts1(k, j, i - 1, 3, 0), &verts1(k, j, i - 1, 2, 0), &verts1(k, j, i - 1, 0, 0), p9);
                    _get_centroid(&verts1(k, j - 1, i - 1, 1, 0), &verts1(k, j - 1, i - 1, 3, 0), &verts1(k, j - 1, i - 1, 2, 0), &verts1(k, j - 1, i - 1, 0, 0), p10);
                    // 0,2,3,4,5交接边点下标
                    double pp0[3];
                    double pp2[3];
                    double pp3[3];
                    double pp4[3];
                    double pp5[3];
                    _get_midpoint(&verts1(k, j, i - 1, 1, 0), &verts1(k, j, i - 1, 3, 0), pp2);
                    _get_midpoint(&verts1(k, j, i - 1, 1, 0), &verts1(k, j, i - 1, 0, 0), pp3);
                    _get_midpoint(&verts1(k, j, i - 1, 1, 0), &verts1(k, j, i - 1, 5, 0), pp5);
                    _get_midpoint(&verts1(k, j - 1, i - 1, 1, 0), &verts1(k, j - 1, i - 1, 3, 0), pp0);
                    _get_midpoint(&verts1(k - 1, j, i - 1, 1, 0), &verts1(k - 1, j, i - 1, 5, 0), pp4);
                    double n0[3];
                    double n1[3];
                    double n2[3];
                    double n3[3];
                    double n6[3];
                    double n7[3];
                    double n9[3];
                    double n10[3];
                    _get_surface_normal(&verts1(k, j, i - 1, 1, 0), pp2, p0, pp5, n0, Axis::XPOSITIVE);
                    _get_surface_normal(&verts1(k, j, i - 1, 1, 0), pp0, p1, pp5, n1, Axis::XPOSITIVE);
                    _get_surface_normal(&verts1(k, j, i - 1, 1, 0), pp0, p2, pp4, n2, Axis::XPOSITIVE);
                    _get_surface_normal(&verts1(k, j, i - 1, 1, 0), pp4, p3, pp2, n3, Axis::XPOSITIVE);
                    _get_surface_normal(&verts1(k, j, i - 1, 1, 0), pp3, p6, pp4, n6, Axis::YPOSITIVE);
                    _get_surface_normal(&verts1(k, j, i - 1, 1, 0), pp5, p7, pp3, n7, Axis::YPOSITIVE);
                    _get_surface_normal(&verts1(k, j, i - 1, 1, 0), pp2, p9, pp3, n9, Axis::ZPOSITIVE);
                    _get_surface_normal(&verts1(k, j, i - 1, 1, 0), pp0, p10, pp3, n10, Axis::ZPOSITIVE);
                    Eigen::Matrix<double, 12, 12> A;
                    A.setZero();
                    A.block(0, 9, 1, 3) << p0[0] - _cx(k, j, i - 1), p0[1] - _cy(k, j, i - 1), p0[2] - _cz(k, j, i - 1);
                    A.block(1, 6, 1, 3) << p1[0] - _cx(k, j - 1, i - 1), p1[1] - _cy(k, j - 1, i - 1), p1[2] - _cz(k, j - 1, i - 1);
                    A.block(2, 0, 1, 3) << p2[0] - _cx(k - 1, j - 1, i - 1), p2[1] - _cy(k - 1, j - 1, i - 1), p2[2] - _cz(k - 1, j - 1, i - 1);
                    A.block(3, 3, 1, 3) << p3[0] - _cx(k - 1, j, i - 1), p3[1] - _cy(k - 1, j, i - 1), p3[2] - _cz(k - 1, j, i - 1);
                    Eigen::Matrix3d matK7, matK4, matK3, matK0;
                    _get_matK(pem1[i7], matK7);
                    _get_matK(pem1[i4], matK4);
                    _get_matK(pem1[i3], matK3);
                    _get_matK(pem1[i0], matK0);
                    Eigen::RowVector3d v0, v1, v2, v3, v6, v7, v9, v10;
                    v0 << n0[0], n0[1], n0[2];
                    v1 << n1[0], n1[1], n1[2];
                    v2 << n2[0], n2[1], n2[2];
                    v3 << n3[0], n3[1], n3[2];
                    v6 << n6[0], n6[1], n6[2];
                    v7 << n7[0], n7[1], n7[2];
                    v9 << n9[0], n9[1], n9[2];
                    v10 << n10[0], n10[1], n10[2];
                    A.block(4, 0, 1, 6) << _cx(k - 1, j - 1, i - 1) - p6[0], _cy(k - 1, j - 1, i - 1) - p6[1], _cz(k - 1, j - 1, i - 1) - p6[2], p6[0] - _cx(k - 1, j, i - 1), p6[1] - _cy(k - 1, j, i - 1), p6[2] - _cz(k - 1, j, i - 1);
                    A.block(5, 0, 1, 3) = v6 * matK0;
                    A.block(5, 3, 1, 3) = -v6 * matK3;
                    A.block(6, 6, 1, 6) << _cx(k, j - 1, i - 1) - p7[0], _cy(k, j - 1, i - 1) - p7[1], _cz(k, j - 1, i - 1) - p7[2], p7[0] - _cx(k, j, i - 1), p7[1] - _cy(k, j, i - 1), p7[2] - _cz(k, j, i - 1);
                    A.block(7, 6, 1, 3) = v7 * matK4;
                    A.block(7, 9, 1, 3) = -v7 * matK7;
                    A.block(8, 3, 1, 3) << _cx(k - 1, j, i - 1) - p9[0], _cy(k - 1, j, i - 1) - p9[1], _cz(k - 1, j, i - 1) - p9[2];
                    A.block(8, 9, 1, 3) << p9[0] - _cx(k, j, i - 1), p9[1] - _cy(k, j, i - 1), p9[2] - _cz(k, j, i - 1);
                    A.block(9, 3, 1, 3) = v9 * matK3;
                    A.block(9, 9, 1, 3) = -v9 * matK7;
                    A.block(10, 0, 1, 3) << _cx(k - 1, j - 1, i - 1) - p10[0], _cy(k - 1, j - 1, i - 1) - p10[1], _cz(k - 1, j - 1, i - 1) - p10[2];
                    A.block(10, 6, 1, 3) << p10[0] - _cx(k, j - 1, i - 1), p10[1] - _cy(k, j - 1, i - 1), p10[2] - _cz(k, j - 1, i - 1);
                    A.block(11, 0, 1, 3) = v10 * matK0;
                    A.block(11, 6, 1, 3) = -v10 * matK4;
                    A = A.inverse();
                    Eigen::RowVectorXd r = (v2 + v6 + v10) * matK0 * A.topRows(3);
                    B[i0] -= pb * (r[0] + r[1] + r[2] + r[3]);
                    C(i0, i0) += -r[2] + r[4] + r[10];
                    C(i0, i3) += -r[3] - r[4] + r[8];
                    C(i0, i4) += -r[1] + r[6] - r[10];
                    C(i0, i7) += -r[0] - r[6] - r[8];
                    r = (v3 - v6 + v9) * matK3 * A.block(3, 0, 3, 12);
                    B[i3] -= pb * (r[0] + r[1] + r[2] + r[3]);
                    C(i3, i0) += -r[2] + r[4] + r[10];
                    C(i3, i3) += -r[3] - r[4] + r[8];
                    C(i3, i4) += -r[1] + r[6] - r[10];
                    C(i3, i7) += -r[0] - r[6] - r[8];
                    r = (v1 - v10 + v7) * matK4 * A.block(6, 0, 3, 12);
                    B[i4] -= pb * (r[0] + r[1] + r[2] + r[3]);
                    C(i4, i0) += -r[2] + r[4] + r[10];
                    C(i4, i3) += -r[3] - r[4] + r[8];
                    C(i4, i4) += -r[1] + r[6] - r[10];
                    C(i4, i7) += -r[0] - r[6] - r[8];
                    r = (v0 - v7 - v9) * matK7 * A.bottomRows(3);
                    B[i7] -= pb * (r[0] + r[1] + r[2] + r[3]);
                    C(i7, i0) += -r[2] + r[4] + r[10];
                    C(i7, i3) += -r[3] - r[4] + r[8];
                    C(i7, i4) += -r[1] + r[6] - r[10];
                    C(i7, i7) += -r[0] - r[6] - r[8];

                    // To be implemented
                }
                else if (j == 0)
                {
                    int i6 = cur, i7 = i6 - 1, i2 = i6 - nx * ny, i3 = i7 - nx * ny;
                    // 0,3,4,5,6,7,8,9交接面中点下标
                    double p0[3];
                    double p3[3];
                    double p4[3];
                    double p5[3];
                    double p6[3];
                    double p7[3];
                    double p8[3];
                    double p9[3];
                    _get_centroid(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 2, 0), &verts1(k, j, i, 4, 0), &verts1(k, j, i, 6, 0), p0);
                    _get_centroid(&verts1(k - 1, j, i, 0, 0), &verts1(k - 1, j, i, 2, 0), &verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 6, 0), p3);
                    _get_centroid(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 1, 0), &verts1(k, j, i, 2, 0), &verts1(k, j, i, 3, 0), p8);
                    _get_centroid(&verts1(k, j, i - 1, 0, 0), &verts1(k, j, i - 1, 1, 0), &verts1(k, j, i - 1, 2, 0), &verts1(k, j, i - 1, 3, 0), p9);
                    _get_centroid(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 1, 0), &verts1(k, j, i, 4, 0), &verts1(k, j, i, 5, 0), p4);
                    _get_centroid(&verts1(k - 1, j, i, 0, 0), &verts1(k - 1, j, i, 1, 0), &verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 5, 0), p5);
                    _get_centroid(&verts1(k - 1, j, i - 1, 0, 0), &verts1(k - 1, j, i - 1, 1, 0), &verts1(k - 1, j, i - 1, 4, 0), &verts1(k - 1, j, i - 1, 5, 0), p6);
                    _get_centroid(&verts1(k, j, i - 1, 0, 0), &verts1(k, j, i - 1, 1, 0), &verts1(k, j, i - 1, 4, 0), &verts1(k, j, i - 1, 5, 0), p7);
                    // 1,2,3,4,5交接边点下
                    double pp1[3];
                    double pp2[3];
                    double pp3[3];
                    double pp4[3];
                    double pp5[3];
                    _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 2, 0), pp2);
                    _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 1, 0), pp1);
                    _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 4, 0), pp5);
                    _get_midpoint(&verts1(k - 1, j, i, 0, 0), &verts1(k - 1, j, i, 4, 0), pp4);
                    _get_midpoint(&verts1(k, j, i - 1, 0, 0), &verts1(k, j, i - 1, 1, 0), pp3);
                    double n0[3];
                    double n3[3];
                    double n4[3];
                    double n5[3];
                    double n6[3];
                    double n7[3];
                    double n8[3];
                    double n9[3];
                    _get_surface_normal(&verts1(k, j, i, 0, 0), pp2, p0, pp5, n0, Axis::XPOSITIVE);
                    _get_surface_normal(&verts1(k, j, i, 0, 0), pp2, p3, pp4, n3, Axis::XPOSITIVE);
                    _get_surface_normal(&verts1(k, j, i, 0, 0), pp1, p4, pp5, n4, Axis::YPOSITIVE);
                    _get_surface_normal(&verts1(k, j, i, 0, 0), pp4, p5, pp1, n5, Axis::YPOSITIVE);
                    _get_surface_normal(&verts1(k, j, i, 0, 0), pp4, p6, pp3, n6, Axis::YPOSITIVE);
                    _get_surface_normal(&verts1(k, j, i, 0, 0), pp5, p7, pp3, n7, Axis::YPOSITIVE);
                    _get_surface_normal(&verts1(k, j, i, 0, 0), pp2, p8, pp1, n8, Axis::ZPOSITIVE);
                    _get_surface_normal(&verts1(k, j, i, 0, 0), pp2, p9, pp3, n9, Axis::ZPOSITIVE);
                    Eigen::Matrix<double, 12, 12> A;
                    A.setZero();
                    Eigen::Matrix3d matK6, matK7, matK2, matK3;
                    _get_matK(pem1[i6], matK6);
                    _get_matK(pem1[i7], matK7);
                    _get_matK(pem1[i2], matK2);
                    _get_matK(pem1[i3], matK3);
                    Eigen::RowVector3d v0, v3, v4, v5, v6, v7, v8, v9;
                    v0 << n0[0], n0[1], n0[2];
                    v3 << n3[0], n3[1], n3[2];
                    v4 << n4[0], n4[1], n4[2];
                    v5 << n5[0], n5[1], n5[2];
                    v6 << n6[0], n6[1], n6[2];
                    v7 << n7[0], n7[1], n7[2];
                    v8 << n8[0], n8[1], n8[2];
                    v9 << n9[0], n9[1], n9[2];
                    A.block(0, 6, 1, 6) << _cx(k, j, i) - p0[0], _cy(k, j, i) - p0[1], _cz(k, j, i) - p0[2], p0[0] - _cx(k, j, i - 1), p0[1] - _cy(k, j, i - 1), p0[2] - _cz(k, j, i - 1);
                    A.block(1, 6, 1, 3) = v0 * matK6;
                    A.block(1, 9, 1, 3) = -v0 * matK7;
                    A.block(2, 0, 1, 6) << _cx(k - 1, j, i) - p3[0], _cy(k - 1, j, i) - p3[1], _cz(k - 1, j, i) - p3[2], p3[0] - _cx(k - 1, j, i - 1), p3[1] - _cy(k - 1, j, i - 1), p3[2] - _cz(k - 1, j, i - 1);
                    A.block(3, 0, 1, 3) = v3 * matK2;
                    A.block(3, 3, 1, 3) = -v3 * matK3;
                    A.block(4, 6, 1, 3) = v4 * matK6;
                    A.block(5, 0, 1, 3) = v5 * matK2;
                    A.block(6, 3, 1, 3) = v6 * matK3;
                    A.block(7, 9, 1, 3) = v7 * matK7;
                    A.block(8, 0, 1, 3) << _cx(k - 1, j, i) - p8[0], _cy(k - 1, j, i) - p8[1], _cz(k - 1, j, i) - p8[2];
                    A.block(8, 6, 1, 3) << p8[0] - _cx(k, j, i), p8[1] - _cy(k, j, i), p8[2] - _cz(k, j, i);
                    A.block(9, 0, 1, 3) = v8 * matK2;
                    A.block(9, 6, 1, 3) = -v8 * matK6;
                    A.block(10, 3, 1, 3) << _cx(k - 1, j, i - 1) - p9[0], _cy(k - 1, j, i - 1) - p9[1], _cz(k - 1, j, i - 1) - p9[2];
                    A.block(10, 9, 1, 3) << p9[0] - _cx(k, j, i - 1), p9[1] - _cy(k, j, i - 1), p9[2] - _cz(k, j, i - 1);
                    A.block(11, 3, 1, 3) = v9 * matK3;
                    A.block(11, 9, 1, 3) = -v9 * matK7;
                    A = A.inverse();
                    Eigen::RowVectorXd r = (-v3 - v5 + v8) * matK2 * A.topRows(3);
                    C(i2, i2) += r[2] + r[8];
                    C(i2, i3) += -r[2] + r[10];
                    C(i2, i6) += r[0] - r[8];
                    C(i2, i7) += -r[0] - r[10];
                    r = (v3 + v9 - v6) * matK3 * A.block(3, 0, 3, 12);
                    C(i3, i2) += r[2] + r[8];
                    C(i3, i3) += -r[2] + r[10];
                    C(i3, i6) += r[0] - r[8];
                    C(i3, i7) += -r[0] - r[10];
                    r = (-v0 - v4 - v8) * matK6 * A.block(6, 0, 3, 12);
                    C(i6, i2) += r[2] + r[8];
                    C(i6, i3) += -r[2] + r[10];
                    C(i6, i6) += r[0] - r[8];
                    C(i6, i7) += -r[0] - r[10];
                    r = (v0 - v7 - v9) * matK7 * A.bottomRows(3);
                    C(i7, i2) += r[2] + r[8];
                    C(i7, i3) += -r[2] + r[10];
                    C(i7, i6) += r[0] - r[8];
                    C(i7, i7) += -r[0] - r[10];
                    // To be implemented
                }
                else if (j == ny)
                {
                    int i5 = cur - nx, i4 = i5 - 1, i1 = i5 - nx * ny, i0 = i4 - nx * ny;
                    // 1,2,4,5,6,7,10,11交接面中点下标
                    double p1[3];
                    double p2[3];
                    double p4[3];
                    double p5[3];
                    double p6[3];
                    double p7[3];
                    double p10[3];
                    double p11[3];
                    _get_centroid(&verts1(k, j - 1, i, 2, 0), &verts1(k, j - 1, i, 0, 0), &verts1(k, j - 1, i, 4, 0), &verts1(k, j - 1, i, 6, 0), p1);
                    _get_centroid(&verts1(k - 1, j - 1, i, 2, 0), &verts1(k - 1, j - 1, i, 0, 0), &verts1(k - 1, j - 1, i, 4, 0), &verts1(k - 1, j - 1, i, 6, 0), p2);
                    _get_centroid(&verts1(k, j - 1, i, 2, 0), &verts1(k, j - 1, i, 3, 0), &verts1(k, j - 1, i, 6, 0), &verts1(k, j - 1, i, 7, 0), p4);
                    _get_centroid(&verts1(k, j - 1, i - 1, 2, 0), &verts1(k, j - 1, i - 1, 3, 0), &verts1(k, j - 1, i - 1, 6, 0), &verts1(k, j - 1, i - 1, 7, 0), p7);
                    _get_centroid(&verts1(k - 1, j - 1, i, 2, 0), &verts1(k - 1, j - 1, i, 3, 0), &verts1(k - 1, j - 1, i, 6, 0), &verts1(k - 1, j - 1, i, 7, 0), p5);
                    _get_centroid(&verts1(k - 1, j - 1, i - 1, 2, 0), &verts1(k - 1, j - 1, i - 1, 3, 0), &verts1(k - 1, j - 1, i - 1, 6, 0), &verts1(k - 1, j - 1, i - 1, 7, 0), p6);
                    _get_centroid(&verts1(k, j - 1, i, 2, 0), &verts1(k, j - 1, i, 0, 0), &verts1(k, j - 1, i, 1, 0), &verts1(k, j - 1, i, 3, 0), p11);
                    _get_centroid(&verts1(k, j - 1, i - 1, 2, 0), &verts1(k, j - 1, i - 1, 0, 0), &verts1(k, j - 1, i - 1, 1, 0), &verts1(k, j - 1, i - 1, 3, 0), p10);
                    // 0,1,3,4,5交接边点下标
                    double pp0[3];
                    double pp1[3];
                    double pp3[3];
                    double pp4[3];
                    double pp5[3];
                    _get_midpoint(&verts1(k, j - 1, i, 2, 0), &verts1(k, j - 1, i, 0, 0), pp0);
                    _get_midpoint(&verts1(k, j - 1, i, 2, 0), &verts1(k, j - 1, i, 3, 0), pp1);
                    _get_midpoint(&verts1(k, j - 1, i, 2, 0), &verts1(k, j - 1, i, 6, 0), pp5);
                    _get_midpoint(&verts1(k, j - 1, i - 1, 2, 0), &verts1(k, j - 1, i - 1, 3, 0), pp3);
                    _get_midpoint(&verts1(k - 1, j - 1, i, 2, 0), &verts1(k - 1, j - 1, i, 6, 0), pp4);
                    double n1[3];
                    double n2[3];
                    double n4[3];
                    double n5[3];
                    double n6[3];
                    double n7[3];
                    double n10[3];
                    double n11[3];
                    _get_surface_normal(&verts1(k, j - 1, i, 2, 0), pp0, p1, pp5, n1, Axis::XPOSITIVE);
                    _get_surface_normal(&verts1(k, j - 1, i, 2, 0), pp0, p2, pp4, n2, Axis::XPOSITIVE);
                    _get_surface_normal(&verts1(k, j - 1, i, 2, 0), pp1, p4, pp5, n4, Axis::YPOSITIVE);
                    _get_surface_normal(&verts1(k, j - 1, i, 2, 0), pp4, p5, pp1, n5, Axis::YPOSITIVE);
                    _get_surface_normal(&verts1(k, j - 1, i, 2, 0), pp4, p6, pp3, n6, Axis::YPOSITIVE);
                    _get_surface_normal(&verts1(k, j - 1, i, 2, 0), pp5, p7, pp3, n7, Axis::YPOSITIVE);
                    _get_surface_normal(&verts1(k, j - 1, i, 2, 0), pp0, p10, pp3, n10, Axis::ZPOSITIVE);
                    _get_surface_normal(&verts1(k, j - 1, i, 2, 0), pp0, p11, pp1, n11, Axis::ZPOSITIVE);
                    Eigen::Matrix<double, 12, 12> A;
                    A.setZero();
                    Eigen::Matrix3d matK5, matK4, matK1, matK0;
                    _get_matK(pem1[i5], matK5);
                    _get_matK(pem1[i4], matK4);
                    _get_matK(pem1[i1], matK1);
                    _get_matK(pem1[i0], matK0);
                    Eigen::RowVector3d v1, v2, v4, v5, v6, v7, v10, v11;
                    v1 << n1[0], n1[1], n1[2];
                    v2 << n2[0], n2[1], n2[2];
                    v4 << n4[0], n4[1], n4[2];
                    v5 << n5[0], n5[1], n5[2];
                    v6 << n6[0], n6[1], n6[2];
                    v7 << n7[0], n7[1], n7[2];
                    v10 << n10[0], n10[1], n10[2];
                    v11 << n11[0], n11[1], n11[2];
                    A.block(0, 6, 1, 6) << _cx(k, j - 1, i - 1) - p1[0], _cy(k, j - 1, i - 1) - p1[1], _cz(k, j - 1, i - 1) - p1[2], p1[0] - _cx(k, j - 1, i), p1[1] - _cy(k, j - 1, i), p1[2] - _cz(k, j - 1, i);
                    A.block(1, 6, 1, 3) = v1 * matK4;
                    A.block(1, 9, 1, 3) = -v1 * matK5;
                    A.block(2, 0, 1, 6) << _cx(k - 1, j - 1, i - 1) - p2[0], _cy(k - 1, j - 1, i - 1) - p2[1], _cz(k - 1, j - 1, i - 1) - p2[2], p2[0] - _cx(k - 1, j - 1, i), p2[1] - _cy(k - 1, j - 1, i), p2[2] - _cz(k - 1, j - 1, i);
                    A.block(3, 0, 1, 3) = v2 * matK0;
                    A.block(3, 3, 1, 3) = -v2 * matK1;
                    A.block(4, 9, 1, 3) = v4 * matK5;
                    A.block(5, 3, 1, 3) = v5 * matK1;
                    A.block(6, 0, 1, 3) = v6 * matK0;
                    A.block(7, 6, 1, 3) = v7 * matK4;
                    A.block(8, 0, 1, 3) << _cx(k - 1, j - 1, i - 1) - p10[0], _cy(k - 1, j - 1, i - 1) - p10[1], _cz(k - 1, j - 1, i - 1) - p10[2];
                    A.block(8, 6, 1, 3) << p10[0] - _cx(k, j - 1, i - 1), p10[1] - _cy(k, j - 1, i - 1), p10[2] - _cz(k, j - 1, i - 1);
                    A.block(9, 0, 1, 3) = v10 * matK0;
                    A.block(9, 6, 1, 3) = -v10 * matK4;
                    A.block(10, 3, 1, 3) << _cx(k - 1, j - 1, i) - p11[0], _cy(k - 1, j - 1, i) - p11[1], _cz(k - 1, j - 1, i) - p11[2];
                    A.block(10, 9, 1, 3) << p11[0] - _cx(k, j - 1, i), p11[1] - _cy(k, j - 1, i), p11[2] - _cz(k, j - 1, i);
                    A.block(11, 3, 1, 3) = v11 * matK1;
                    A.block(11, 9, 1, 3) = -v11 * matK5;
                    A = A.inverse();
                    Eigen::RowVectorXd r = (v2 + v6 + v10) * matK0 * A.topRows(3);
                    C(i0, i0) += r[2] + r[8];
                    C(i0, i1) += -r[2] + r[10];
                    C(i0, i4) += r[0] - r[8];
                    C(i0, i5) += -r[0] - r[10];
                    r = (-v2 + v5 + v11) * matK1 * A.block(3, 0, 3, 12);
                    C(i1, i0) += r[2] + r[8];
                    C(i1, i1) += -r[2] + r[10];
                    C(i1, i4) += r[0] - r[8];
                    C(i1, i5) += -r[0] - r[10];
                    r = (v1 + v7 - v10) * matK4 * A.block(6, 0, 3, 12);
                    C(i4, i0) += r[2] + r[8];
                    C(i4, i1) += -r[2] + r[10];
                    C(i4, i4) += r[0] - r[8];
                    C(i4, i5) += -r[0] - r[10];
                    r = (-v1 + v4 - v11) * matK5 * A.bottomRows(3);
                    C(i5, i0) += r[2] + r[8];
                    C(i5, i1) += -r[2] + r[10];
                    C(i5, i4) += r[0] - r[8];
                    C(i5, i5) += -r[0] - r[10];
                    // To be implemented
                }
                else if (k == 0)
                {
                    int i6 = cur, i7 = cur - 1, i5 = cur - nx, i4 = i5 - 1;
                    // 0,1,4,7,8,9,10,11交接面中点下标
                    double p0[3];
                    double p1[3];
                    double p4[3];
                    double p7[3];
                    double p8[3];
                    double p9[3];
                    double p10[3];
                    double p11[3];
                    _get_centroid(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 2, 0), &verts1(k, j, i, 4, 0), &verts1(k, j, i, 6, 0), p0);
                    _get_centroid(&verts1(k, j - 1, i, 0, 0), &verts1(k, j - 1, i, 2, 0), &verts1(k, j - 1, i, 4, 0), &verts1(k, j - 1, i, 6, 0), p1);
                    _get_centroid(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 1, 0), &verts1(k, j, i, 4, 0), &verts1(k, j, i, 5, 0), p4);
                    _get_centroid(&verts1(k, j, i - 1, 0, 0), &verts1(k, j, i - 1, 1, 0), &verts1(k, j, i - 1, 4, 0), &verts1(k, j, i - 1, 5, 0), p7);
                    _get_centroid(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 1, 0), &verts1(k, j, i, 2, 0), &verts1(k, j, i, 3, 0), p8);
                    _get_centroid(&verts1(k, j - 1, i, 0, 0), &verts1(k, j - 1, i, 1, 0), &verts1(k, j - 1, i, 2, 0), &verts1(k, j - 1, i, 3, 0), p11);
                    _get_centroid(&verts1(k, j, i - 1, 0, 0), &verts1(k, j, i - 1, 1, 0), &verts1(k, j, i - 1, 2, 0), &verts1(k, j, i - 1, 3, 0), p9);
                    _get_centroid(&verts1(k, j - 1, i - 1, 0, 0), &verts1(k, j - 1, i - 1, 1, 0), &verts1(k, j - 1, i - 1, 2, 0), &verts1(k, j - 1, i - 1, 3, 0), p10);
                    // 0,1,2,3,5 交接边点下标
                    double pp0[3];
                    double pp1[3];
                    double pp2[3];
                    double pp3[3];
                    double pp5[3];
                    _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 2, 0), pp2);
                    _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 1, 0), pp1);
                    _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 4, 0), pp5);
                    _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j - 1, i, 0, 0), pp0);
                    _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j, i - 1, 0, 0), pp3);
                    double n0[3];
                    double n1[3];
                    double n4[3];
                    double n7[3];
                    double n8[3];
                    double n9[3];
                    double n10[3];
                    double n11[3];
                    _get_surface_normal(&verts1(k, j, i, 0, 0), pp2, p0, pp5, n0, Axis::XPOSITIVE);
                    _get_surface_normal(&verts1(k, j, i, 0, 0), pp5, p1, pp0, n1, Axis::XPOSITIVE);
                    _get_surface_normal(&verts1(k, j, i, 0, 0), pp1, p4, pp5, n4, Axis::YPOSITIVE);
                    _get_surface_normal(&verts1(k, j, i, 0, 0), pp5, p7, pp3, n7, Axis::YPOSITIVE);
                    _get_surface_normal(&verts1(k, j, i, 0, 0), pp2, p8, pp1, n8, Axis::ZPOSITIVE);
                    _get_surface_normal(&verts1(k, j, i, 0, 0), pp0, p11, pp1, n11, Axis::ZPOSITIVE);
                    _get_surface_normal(&verts1(k, j, i, 0, 0), pp0, p10, pp3, n10, Axis::ZPOSITIVE);
                    _get_surface_normal(&verts1(k, j, i, 0, 0), pp3, p9, pp2, n9, Axis::ZPOSITIVE);
                    Eigen::Matrix<double, 12, 12> A;
                    A.setZero();
                    Eigen::Matrix3d matK6, matK7, matK5, matK4;
                    _get_matK(pem1[i6], matK6);
                    _get_matK(pem1[i7], matK7);
                    _get_matK(pem1[i5], matK5);
                    _get_matK(pem1[i4], matK4);
                    Eigen::RowVector3d v0, v1, v4, v7, v8, v9, v10, v11;
                    v0 << n0[0], n0[1], n0[2];
                    v1 << n1[0], n1[1], n1[2];
                    v4 << n4[0], n4[1], n4[2];
                    v7 << n7[0], n7[1], n7[2];
                    v8 << n8[0], n8[1], n8[2];
                    v9 << n9[0], n9[1], n9[2];
                    v10 << n10[0], n10[1], n10[2];
                    v11 << n11[0], n11[1], n11[2];
                    A.block(0, 6, 1, 6) << _cx(k, j, i) - p0[0], _cy(k, j, i) - p0[1], _cz(k, j, i) - p0[2], p0[0] - _cx(k, j, i - 1), p0[1] - _cy(k, j, i - 1), p0[2] - _cz(k, j, i - 1);
                    A.block(1, 6, 1, 3) = v0 * matK6;
                    A.block(1, 9, 1, 3) = -v0 * matK7;
                    A.block(2, 0, 1, 6) << p1[0] - _cx(k, j - 1, i - 1), p1[1] - _cy(k, j - 1, i - 1), p1[2] - _cz(k, j - 1, i - 1), _cx(k, j - 1, i) - p1[0], _cy(k, j - 1, i) - p1[1], _cz(k, j - 1, i) - p1[2];
                    A.block(3, 0, 1, 3) = v1 * matK4;
                    A.block(3, 3, 1, 3) = -v1 * matK5;
                    A.block(4, 3, 1, 6) << _cx(k, j - 1, i) - p4[0], _cy(k, j - 1, i) - p4[1], _cz(k, j - 1, i) - p4[2], p4[0] - _cx(k, j, i), p4[1] - _cy(k, j, i), p4[2] - _cz(k, j, i);
                    A.block(5, 3, 1, 3) = v4 * matK5;
                    A.block(5, 6, 1, 3) = -v4 * matK6;
                    A.block(6, 0, 1, 3) << _cx(k, j - 1, i - 1) - p7[0], _cy(k, j - 1, i - 1) - p7[1], _cz(k, j - 1, i - 1) - p7[2];
                    A.block(6, 9, 1, 3) << p7[0] - _cx(k, j, i - 1), p7[1] - _cy(k, j, i - 1), p7[2] - _cz(k, j, i - 1);
                    A.block(7, 0, 1, 3) = v7 * matK4;
                    A.block(7, 9, 1, 3) = -v7 * matK7;
                    A.block(8, 6, 1, 3) = v8 * matK6;
                    A.block(9, 9, 1, 3) = v9 * matK7;
                    A.block(10, 0, 1, 3) = v10 * matK4;
                    A.block(11, 3, 1, 3) = v11 * matK5;
                    A = A.inverse();
                    Eigen::RowVectorXd r = (v1 + v7 - v10) * matK4 * A.topRows(3);
                    C(i4, i4) += -r[2] + r[6];
                    C(i4, i5) += r[2] + r[4];
                    C(i4, i6) += r[0] - r[4];
                    C(i4, i7) += -r[0] - r[6];
                    r = (-v1 + v4 - v11) * matK5 * A.block(3, 0, 3, 12);
                    C(i5, i4) += -r[2] + r[6];
                    C(i5, i5) += r[2] + r[4];
                    C(i5, i6) += r[0] - r[4];
                    C(i5, i7) += -r[0] - r[6];
                    r = (-v0 - v8 - v4) * matK6 * A.block(6, 0, 3, 12);
                    C(i6, i4) += -r[2] + r[6];
                    C(i6, i5) += r[2] + r[4];
                    C(i6, i6) += r[0] - r[4];
                    C(i6, i7) += -r[0] - r[6];
                    r = (v0 - v9 - v7) * matK7 * A.bottomRows(3);
                    C(i7, i4) += -r[2] + r[6];
                    C(i7, i5) += r[2] + r[4];
                    C(i7, i6) += r[0] - r[4];
                    C(i7, i7) += -r[0] - r[6];
                    // To be implemented
                }
                else if (k == nz)
                {
                    int i2 = cur - nx * ny, i3 = i2 - 1, i1 = i2 - nx, i0 = i1 - 1;
                    // 2,3,5,6,8,9,10,11交接面中点下标
                    double p2[3];
                    double p3[3];
                    double p5[3];
                    double p6[3];
                    double p8[3];
                    double p9[3];
                    double p10[3];
                    double p11[3];
                    _get_centroid(&verts1(k - 1, j, i, 2, 0), &verts1(k - 1, j, i, 0, 0), &verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 6, 0), p3);
                    _get_centroid(&verts1(k - 1, j - 1, i, 2, 0), &verts1(k - 1, j - 1, i, 0, 0), &verts1(k - 1, j - 1, i, 4, 0), &verts1(k - 1, j - 1, i, 6, 0), p2);
                    _get_centroid(&verts1(k - 1, j, i, 0, 0), &verts1(k - 1, j, i, 1, 0), &verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 5, 0), p5);
                    _get_centroid(&verts1(k - 1, j, i - 1, 0, 0), &verts1(k - 1, j, i - 1, 1, 0), &verts1(k - 1, j, i - 1, 4, 0), &verts1(k - 1, j, i - 1, 5, 0), p6);
                    _get_centroid(&verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 5, 0), &verts1(k - 1, j, i, 6, 0), &verts1(k - 1, j, i, 7, 0), p8);
                    _get_centroid(&verts1(k - 1, j, i - 1, 4, 0), &verts1(k - 1, j, i - 1, 5, 0), &verts1(k - 1, j, i - 1, 6, 0), &verts1(k - 1, j, i - 1, 7, 0), p9);
                    _get_centroid(&verts1(k - 1, j - 1, i, 4, 0), &verts1(k - 1, j - 1, i, 5, 0), &verts1(k - 1, j - 1, i, 6, 0), &verts1(k - 1, j - 1, i, 7, 0), p11);
                    _get_centroid(&verts1(k - 1, j - 1, i - 1, 4, 0), &verts1(k - 1, j - 1, i - 1, 5, 0), &verts1(k - 1, j - 1, i - 1, 6, 0), &verts1(k - 1, j - 1, i - 1, 7, 0), p10);
                    // 0,1,2,3,4 交接边点下标
                    double pp0[3];
                    double pp1[3];
                    double pp2[3];
                    double pp3[3];
                    double pp4[3];
                    _get_midpoint(&verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 6, 0), pp2);
                    _get_midpoint(&verts1(k - 1, j - 1, i, 4, 0), &verts1(k - 1, j - 1, i, 6, 0), pp0);
                    _get_midpoint(&verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 5, 0), pp1);
                    _get_midpoint(&verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i - 1, 4, 0), pp3);
                    _get_midpoint(&verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 0, 0), pp4);
                    double n2[3];
                    double n3[3];
                    double n5[3];
                    double n6[3];
                    double n8[3];
                    double n9[3];
                    double n10[3];
                    double n11[3];
                    _get_surface_normal(&verts1(k - 1, j, i, 4, 0), pp0, p2, pp4, n2, Axis::XPOSITIVE);
                    _get_surface_normal(&verts1(k - 1, j, i, 4, 0), pp4, p3, pp2, n3, Axis::XPOSITIVE);
                    _get_surface_normal(&verts1(k - 1, j, i, 4, 0), pp1, p5, pp4, n5, Axis::YPOSITIVE);
                    _get_surface_normal(&verts1(k - 1, j, i, 4, 0), pp4, p6, pp3, n6, Axis::YPOSITIVE);
                    _get_surface_normal(&verts1(k - 1, j, i, 4, 0), pp2, p8, pp1, n8, Axis::ZPOSITIVE);
                    _get_surface_normal(&verts1(k - 1, j, i, 4, 0), pp3, p9, pp2, n9, Axis::ZPOSITIVE);
                    _get_surface_normal(&verts1(k - 1, j, i, 4, 0), pp0, p11, pp1, n11, Axis::ZPOSITIVE);
                    _get_surface_normal(&verts1(k - 1, j, i, 4, 0), pp0, p10, pp3, n10, Axis::ZPOSITIVE);
                    Eigen::Matrix<double, 12, 12> A;
                    A.setZero();
                    Eigen::Matrix3d matK2, matK3, matK1, matK0;
                    _get_matK(pem1[i2], matK2);
                    _get_matK(pem1[i3], matK3);
                    _get_matK(pem1[i1], matK1);
                    _get_matK(pem1[i0], matK0);
                    Eigen::RowVector3d v2, v3, v5, v6, v8, v9, v10, v11;
                    v2 << n2[0], n2[1], n2[2];
                    v3 << n3[0], n3[1], n3[2];
                    v5 << n5[0], n5[1], n5[2];
                    v6 << n6[0], n6[1], n6[2];
                    v8 << n8[0], n8[1], n8[2];
                    v9 << n9[0], n9[1], n9[2];
                    v10 << n10[0], n10[1], n10[2];
                    v11 << n11[0], n11[1], n11[2];
                    A.block(0, 0, 1, 6) << _cx(k - 1, j - 1, i - 1) - p2[0], _cy(k - 1, j - 1, i - 1) - p2[1], _cz(k - 1, j - 1, i - 1) - p2[2], p2[0] - _cx(k - 1, j - 1, i), p2[1] - _cy(k - 1, j - 1, i), p2[2] - _cz(k - 1, j - 1, i);
                    A.block(1, 0, 1, 3) = v2 * matK0;
                    A.block(1, 3, 1, 3) = -v2 * matK1;
                    A.block(2, 6, 1, 6) << p3[0] - _cx(k - 1, j, i), p3[1] - _cy(k - 1, j, i), p3[2] - _cz(k - 1, j, i), _cx(k - 1, j, i - 1) - p3[0], _cy(k - 1, j, i - 1) - p3[1], _cz(k - 1, j, i - 1) - p3[2];
                    A.block(3, 6, 1, 3) = v3 * matK2;
                    A.block(3, 9, 1, 3) = -v3 * matK3;
                    A.block(4, 3, 1, 6) << _cx(k - 1, j - 1, i) - p5[0], _cy(k - 1, j - 1, i) - p5[1], _cz(k - 1, j - 1, i) - p5[2], p5[0] - _cx(k - 1, j, i), p5[1] - _cy(k - 1, j, i), p5[2] - _cz(k - 1, j, i);
                    A.block(5, 3, 1, 3) = v5 * matK1;
                    A.block(5, 6, 1, 3) = -v5 * matK2;
                    A.block(6, 0, 1, 3) << _cx(k - 1, j - 1, i - 1) - p6[0], _cy(k - 1, j - 1, i - 1) - p6[1], _cz(k - 1, j - 1, i - 1) - p6[2];
                    A.block(6, 9, 1, 3) << p6[0] - _cx(k - 1, j, i - 1), p6[1] - _cy(k - 1, j, i - 1), p6[2] - _cz(k - 1, j, i - 1);
                    A.block(7, 0, 1, 3) = v6 * matK0;
                    A.block(7, 9, 1, 3) = -v6 * matK3;
                    A.block(8, 6, 1, 3) = v8 * matK2;
                    A.block(9, 9, 1, 3) = v9 * matK3;
                    A.block(10, 0, 1, 3) = v10 * matK0;
                    A.block(11, 3, 1, 3) = v11 * matK1;
                    A = A.inverse();
                    Eigen::RowVectorXd r = (v2 + v6 + v10) * matK0 * A.topRows(3);
                    C(i0, i0) += r[0] + r[6];
                    C(i0, i1) += -r[0] + r[4];
                    C(i0, i2) += -r[2] - r[4];
                    C(i0, i3) += r[2] - r[6];
                    r = (-v2 + v5 + v11) * matK1 * A.block(3, 0, 3, 12);
                    C(i1, i0) += r[0] + r[6];
                    C(i1, i1) += -r[0] + r[4];
                    C(i1, i2) += -r[2] - r[4];
                    C(i1, i3) += r[2] - r[6];
                    r = (-v3 - v5 + v8) * matK2 * A.block(6, 0, 3, 12);
                    C(i2, i0) += r[0] + r[6];
                    C(i2, i1) += -r[0] + r[4];
                    C(i2, i2) += -r[2] - r[4];
                    C(i2, i3) += r[2] - r[6];
                    r = (v3 - v6 + v9) * matK3 * A.bottomRows(3);
                    C(i3, i0) += r[0] + r[6];
                    C(i3, i1) += -r[0] + r[4];
                    C(i3, i2) += -r[2] - r[4];
                    C(i3, i3) += r[2] - r[6];
                    //  To be implemented
                }
                else
                {
                    int idx[8];
                    idx[6] = cur;         // i6
                    idx[7] = idx[6] - 1;  // i7
                    idx[4] = idx[7] - nx; // i4
                    idx[5] = idx[6] - nx; // i5

                    int layerOffset = nx * ny;
                    idx[0] = idx[4] - layerOffset; // i0
                    idx[1] = idx[5] - layerOffset; // i1
                    idx[2] = idx[6] - layerOffset; // i2
                    idx[3] = idx[7] - layerOffset; // i3
                    double p[12][3];
                    _get_centroid(&verts1(k, j, i, 2, 0), &verts1(k, j, i, 0, 0), &verts1(k, j, i, 4, 0), &verts1(k, j, i, 6, 0), p[0]);
                    _get_centroid(&verts1(k, j - 1, i, 2, 0), &verts1(k, j - 1, i, 0, 0), &verts1(k, j - 1, i, 4, 0), &verts1(k, j - 1, i, 6, 0), p[1]);
                    _get_centroid(&verts1(k - 1, j - 1, i, 2, 0), &verts1(k - 1, j - 1, i, 0, 0), &verts1(k - 1, j - 1, i, 4, 0), &verts1(k - 1, j - 1, i, 6, 0), p[2]);
                    _get_centroid(&verts1(k - 1, j, i, 2, 0), &verts1(k - 1, j, i, 0, 0), &verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 6, 0), p[3]);
                    _get_centroid(&verts1(k, j, i, 1, 0), &verts1(k, j, i, 0, 0), &verts1(k, j, i, 4, 0), &verts1(k, j, i, 5, 0), p[4]);
                    _get_centroid(&verts1(k - 1, j, i, 1, 0), &verts1(k - 1, j, i, 0, 0), &verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 5, 0), p[5]);
                    _get_centroid(&verts1(k - 1, j, i - 1, 1, 0), &verts1(k - 1, j, i - 1, 0, 0), &verts1(k - 1, j, i - 1, 4, 0), &verts1(k - 1, j, i - 1, 5, 0), p[6]);
                    _get_centroid(&verts1(k, j, i - 1, 1, 0), &verts1(k, j, i - 1, 0, 0), &verts1(k, j, i - 1, 4, 0), &verts1(k, j, i - 1, 5, 0), p[7]);
                    _get_centroid(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 1, 0), &verts1(k, j, i, 2, 0), &verts1(k, j, i, 3, 0), p[8]);
                    _get_centroid(&verts1(k, j, i - 1, 0, 0), &verts1(k, j, i - 1, 1, 0), &verts1(k, j, i - 1, 2, 0), &verts1(k, j, i - 1, 3, 0), p[9]);
                    _get_centroid(&verts1(k, j - 1, i - 1, 0, 0), &verts1(k, j - 1, i - 1, 1, 0), &verts1(k, j - 1, i - 1, 2, 0), &verts1(k, j - 1, i - 1, 3, 0), p[10]);
                    _get_centroid(&verts1(k, j - 1, i, 0, 0), &verts1(k, j - 1, i, 1, 0), &verts1(k, j - 1, i, 2, 0), &verts1(k, j - 1, i, 3, 0), p[11]);
                    double pp[6][3];
                    _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 1, 0), pp[1]);
                    _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j, i - 1, 0, 0), pp[3]);
                    _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 2, 0), pp[2]);
                    _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j - 1, i, 0, 0), pp[0]);
                    _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 4, 0), pp[5]);
                    _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k - 1, j, i, 0, 0), pp[4]);
                    double n[12][3];
                    _get_surface_normal(&verts1(k, j, i, 0, 0), pp[2], p[0], pp[5], n[0], Axis::XPOSITIVE);
                    _get_surface_normal(&verts1(k, j, i, 0, 0), pp[0], p[1], pp[5], n[1], Axis::XPOSITIVE);
                    _get_surface_normal(&verts1(k, j, i, 0, 0), pp[0], p[2], pp[4], n[2], Axis::XPOSITIVE);
                    _get_surface_normal(&verts1(k, j, i, 0, 0), pp[4], p[3], pp[2], n[3], Axis::XPOSITIVE);
                    _get_surface_normal(&verts1(k, j, i, 0, 0), pp[1], p[4], pp[5], n[4], Axis::YPOSITIVE);
                    _get_surface_normal(&verts1(k, j, i, 0, 0), pp[1], p[5], pp[4], n[5], Axis::YPOSITIVE);
                    _get_surface_normal(&verts1(k, j, i, 0, 0), pp[3], p[6], pp[4], n[6], Axis::YPOSITIVE);
                    _get_surface_normal(&verts1(k, j, i, 0, 0), pp[3], p[7], pp[5], n[7], Axis::YPOSITIVE);
                    _get_surface_normal(&verts1(k, j, i, 0, 0), pp[1], p[8], pp[2], n[8], Axis::ZPOSITIVE);
                    _get_surface_normal(&verts1(k, j, i, 0, 0), pp[3], p[9], pp[2], n[9], Axis::ZPOSITIVE);
                    _get_surface_normal(&verts1(k, j, i, 0, 0), pp[3], p[10], pp[0], n[10], Axis::ZPOSITIVE);
                    _get_surface_normal(&verts1(k, j, i, 0, 0), pp[0], p[11], pp[1], n[11], Axis::ZPOSITIVE);
                    Eigen::Matrix<double, 24, 24> A;
                    A.setZero();
                    Eigen::Matrix3d matK[8];
                    for (int l = 0; l < 8; ++l)
                        _get_matK(pem1[idx[l]], matK[l]);
                    Eigen::RowVector3d v[12];
                    for (int l = 0; l < 12; ++l)
                        v[l] << n[l][0], n[l][1], n[l][2];
                    A.block(0, 18, 1, 6) << _cx(k, j, i) - p[0][0], _cy(k, j, i) - p[0][1], _cz(k, j, i) - p[0][2], p[0][0] - _cx(k, j, i - 1), p[0][1] - _cy(k, j, i - 1), p[0][2] - _cz(k, j, i - 1);
                    A.block(1, 18, 1, 3) = v[0] * matK[6];
                    A.block(1, 21, 1, 3) = -v[0] * matK[7];
                    A.block(2, 12, 1, 6) << _cx(k, j - 1, i - 1) - p[1][0], _cy(k, j - 1, i - 1) - p[1][1], _cz(k, j - 1, i - 1) - p[1][2], p[1][0] - _cx(k, j - 1, i), p[1][1] - _cy(k, j - 1, i), p[1][2] - _cz(k, j - 1, i);
                    A.block(3, 12, 1, 3) = v[1] * matK[4];
                    A.block(3, 15, 1, 3) = -v[1] * matK[5];
                    A.block(4, 0, 1, 6) << _cx(k - 1, j - 1, i - 1) - p[2][0], _cy(k - 1, j - 1, i - 1) - p[2][1], _cz(k - 1, j - 1, i - 1) - p[2][2], p[2][0] - _cx(k - 1, j - 1, i), p[2][1] - _cy(k - 1, j - 1, i), p[2][2] - _cz(k - 1, j - 1, i);
                    A.block(5, 0, 1, 3) = v[2] * matK[0];
                    A.block(5, 3, 1, 3) = -v[2] * matK[1];
                    A.block(6, 6, 1, 6) << _cx(k - 1, j, i) - p[3][0], _cy(k - 1, j, i) - p[3][1], _cz(k - 1, j, i) - p[3][2], p[3][0] - _cx(k - 1, j, i - 1), p[3][1] - _cy(k - 1, j, i - 1), p[3][2] - _cz(k - 1, j, i - 1);
                    A.block(7, 6, 1, 3) = v[3] * matK[2];
                    A.block(7, 9, 1, 3) = -v[3] * matK[3];
                    A.block(8, 15, 1, 6) << _cx(k, j - 1, i) - p[4][0], _cy(k, j - 1, i) - p[4][1], _cz(k, j - 1, i) - p[4][2], p[4][0] - _cx(k, j, i), p[4][1] - _cy(k, j, i), p[4][2] - _cz(k, j, i);
                    A.block(9, 15, 1, 3) = v[4] * matK[5];
                    A.block(9, 18, 1, 3) = -v[4] * matK[6];
                    A.block(10, 3, 1, 6) << _cx(k - 1, j - 1, i) - p[5][0], _cy(k - 1, j - 1, i) - p[5][1], _cz(k - 1, j - 1, i) - p[5][2], p[5][0] - _cx(k - 1, j, i), p[5][1] - _cy(k - 1, j, i), p[5][2] - _cz(k - 1, j, i);
                    A.block(11, 3, 1, 3) = v[5] * matK[1];
                    A.block(11, 6, 1, 3) = -v[5] * matK[2];
                    A.block(12, 0, 1, 3) << _cx(k - 1, j - 1, i - 1) - p[6][0], _cy(k - 1, j - 1, i - 1) - p[6][1], _cz(k - 1, j - 1, i - 1) - p[6][2];
                    A.block(12, 9, 1, 3) << p[6][0] - _cx(k - 1, j, i - 1), p[6][1] - _cy(k - 1, j, i - 1), p[6][2] - _cz(k - 1, j, i - 1);
                    A.block(13, 0, 1, 3) = v[6] * matK[0];
                    A.block(13, 9, 1, 3) = -v[6] * matK[3];
                    A.block(14, 12, 1, 3) << _cx(k, j - 1, i - 1) - p[7][0], _cy(k, j - 1, i - 1) - p[7][1], _cz(k, j - 1, i - 1) - p[7][2];
                    A.block(14, 21, 1, 3) << p[7][0] - _cx(k, j, i - 1), p[7][1] - _cy(k, j, i - 1), p[7][2] - _cz(k, j, i - 1);
                    A.block(15, 12, 1, 3) = v[7] * matK[4];
                    A.block(15, 21, 1, 3) = -v[7] * matK[7];
                    A.block(16, 6, 1, 3) << _cx(k - 1, j, i) - p[8][0], _cy(k - 1, j, i) - p[8][1], _cz(k - 1, j, i) - p[8][2];
                    A.block(16, 18, 1, 3) << p[8][0] - _cx(k, j, i), p[8][1] - _cy(k, j, i), p[8][2] - _cz(k, j, i);
                    A.block(17, 6, 1, 3) = v[8] * matK[2];
                    A.block(17, 18, 1, 3) = -v[8] * matK[6];
                    A.block(18, 9, 1, 3) << _cx(k - 1, j, i - 1) - p[9][0], _cy(k - 1, j, i - 1) - p[9][1], _cz(k - 1, j, i - 1) - p[9][2];
                    A.block(18, 21, 1, 3) << p[9][0] - _cx(k, j, i - 1), p[9][1] - _cy(k, j, i - 1), p[9][2] - _cz(k, j, i - 1);
                    A.block(19, 9, 1, 3) = v[9] * matK[3];
                    A.block(19, 21, 1, 3) = -v[9] * matK[7];
                    A.block(20, 0, 1, 3) << _cx(k - 1, j - 1, i - 1) - p[10][0], _cy(k - 1, j - 1, i - 1) - p[10][1], _cz(k - 1, j - 1, i - 1) - p[10][2];
                    A.block(20, 12, 1, 3) << p[10][0] - _cx(k, j - 1, i - 1), p[10][1] - _cy(k, j - 1, i - 1), p[10][2] - _cz(k, j - 1, i - 1);
                    A.block(21, 0, 1, 3) = v[10] * matK[0];
                    A.block(21, 12, 1, 3) = -v[10] * matK[4];
                    A.block(22, 3, 1, 3) << _cx(k - 1, j - 1, i) - p[11][0], _cy(k - 1, j - 1, i) - p[11][1], _cz(k - 1, j - 1, i) - p[11][2];
                    A.block(22, 15, 1, 3) << p[11][0] - _cx(k, j - 1, i), p[11][1] - _cy(k, j - 1, i), p[11][2] - _cz(k, j - 1, i);
                    A.block(23, 3, 1, 3) = v[11] * matK[1];
                    A.block(23, 15, 1, 3) = -v[11] * matK[5];
                    A = A.inverse();
                    Eigen::RowVectorXd r = (v[2] + v[6] + v[10]) * matK[0] * A.topRows(3);
                    C(idx[0], idx[0]) += r[4] + r[12] + r[20];
                    C(idx[0], idx[1]) += -r[4] + r[10] + r[22];
                    C(idx[0], idx[2]) += r[6] - r[10] + r[16];
                    C(idx[0], idx[3]) += -r[6] - r[12] + r[18];
                    C(idx[0], idx[4]) += r[2] + r[14] - r[20];
                    C(idx[0], idx[5]) += -r[2] + r[8] - r[22];
                    C(idx[0], idx[6]) += r[0] - r[8] - r[16];
                    C(idx[0], idx[7]) += -r[0] - r[14] - r[18];
                    r = (-v[2] + v[5] + v[11]) * matK[1] * A.block(3, 0, 3, 24);
                    C(idx[1], idx[0]) += r[4] + r[12] + r[20];
                    C(idx[1], idx[1]) += -r[4] + r[10] + r[22];
                    C(idx[1], idx[2]) += r[6] - r[10] + r[16];
                    C(idx[1], idx[3]) += -r[6] - r[12] + r[18];
                    C(idx[1], idx[4]) += r[2] + r[14] - r[20];
                    C(idx[1], idx[5]) += -r[2] + r[8] - r[22];
                    C(idx[1], idx[6]) += r[0] - r[8] - r[16];
                    C(idx[1], idx[7]) += -r[0] - r[14] - r[18];
                    r = (-v[3] - v[5] + v[8]) * matK[2] * A.block(6, 0, 3, 24);
                    C(idx[2], idx[0]) += r[4] + r[12] + r[20];
                    C(idx[2], idx[1]) += -r[4] + r[10] + r[22];
                    C(idx[2], idx[2]) += r[6] - r[10] + r[16];
                    C(idx[2], idx[3]) += -r[6] - r[12] + r[18];
                    C(idx[2], idx[4]) += r[2] + r[14] - r[20];
                    C(idx[2], idx[5]) += -r[2] + r[8] - r[22];
                    C(idx[2], idx[6]) += r[0] - r[8] - r[16];
                    C(idx[2], idx[7]) += -r[0] - r[14] - r[18];
                    r = (v[3] - v[6] + v[9]) * matK[3] * A.block(9, 0, 3, 24);
                    C(idx[3], idx[0]) += r[4] + r[12] + r[20];
                    C(idx[3], idx[1]) += -r[4] + r[10] + r[22];
                    C(idx[3], idx[2]) += r[6] - r[10] + r[16];
                    C(idx[3], idx[3]) += -r[6] - r[12] + r[18];
                    C(idx[3], idx[4]) += r[2] + r[14] - r[20];
                    C(idx[3], idx[5]) += -r[2] + r[8] - r[22];
                    C(idx[3], idx[6]) += r[0] - r[8] - r[16];
                    C(idx[3], idx[7]) += -r[0] - r[14] - r[18];
                    r = (v[1] + v[7] - v[10]) * matK[4] * A.block(12, 0, 3, 24);
                    C(idx[4], idx[0]) += r[4] + r[12] + r[20];
                    C(idx[4], idx[1]) += -r[4] + r[10] + r[22];
                    C(idx[4], idx[2]) += r[6] - r[10] + r[16];
                    C(idx[4], idx[3]) += -r[6] - r[12] + r[18];
                    C(idx[4], idx[4]) += r[2] + r[14] - r[20];
                    C(idx[4], idx[5]) += -r[2] + r[8] - r[22];
                    C(idx[4], idx[6]) += r[0] - r[8] - r[16];
                    C(idx[4], idx[7]) += -r[0] - r[14] - r[18];
                    r = (-v[1] + v[4] - v[11]) * matK[5] * A.block(15, 0, 3, 24);
                    C(idx[5], idx[0]) += r[4] + r[12] + r[20];
                    C(idx[5], idx[1]) += -r[4] + r[10] + r[22];
                    C(idx[5], idx[2]) += r[6] - r[10] + r[16];
                    C(idx[5], idx[3]) += -r[6] - r[12] + r[18];
                    C(idx[5], idx[4]) += r[2] + r[14] - r[20];
                    C(idx[5], idx[5]) += -r[2] + r[8] - r[22];
                    C(idx[5], idx[6]) += r[0] - r[8] - r[16];
                    C(idx[5], idx[7]) += -r[0] - r[14] - r[18];
                    r = (-v[0] - v[4] - v[8]) * matK[6] * A.block(18, 0, 3, 24);
                    C(idx[6], idx[0]) += r[4] + r[12] + r[20];
                    C(idx[6], idx[1]) += -r[4] + r[10] + r[22];
                    C(idx[6], idx[2]) += r[6] - r[10] + r[16];
                    C(idx[6], idx[3]) += -r[6] - r[12] + r[18];
                    C(idx[6], idx[4]) += r[2] + r[14] - r[20];
                    C(idx[6], idx[5]) += -r[2] + r[8] - r[22];
                    C(idx[6], idx[6]) += r[0] - r[8] - r[16];
                    C(idx[6], idx[7]) += -r[0] - r[14] - r[18];
                    r = (v[0] - v[7] - v[9]) * matK[7] * A.block(21, 0, 3, 24);
                    C(idx[7], idx[0]) += r[4] + r[12] + r[20];
                    C(idx[7], idx[1]) += -r[4] + r[10] + r[22];
                    C(idx[7], idx[2]) += r[6] - r[10] + r[16];
                    C(idx[7], idx[3]) += -r[6] - r[12] + r[18];
                    C(idx[7], idx[4]) += r[2] + r[14] - r[20];
                    C(idx[7], idx[5]) += -r[2] + r[8] - r[22];
                    C(idx[7], idx[6]) += r[0] - r[8] - r[16];
                    C(idx[7], idx[7]) += -r[0] - r[14] - r[18];
                    // 内部角点
                }
            }
    pmgmres_ilu_cr(nx * ny * nz, nnz, Ptr, Idx, Val, p, B, 1000, 1000, 1e-5, 1e-5);
    Q = 0;
    for (int i = 0; i <= 0; ++i)
        for (int k = 0; k <= nz; ++k)
            for (int j = 0; j <= ny; ++j)
            {
                int cur = k * nx * ny + j * nx;
                if (j == 0 && k == 0)
                {
                    // 0,4,8,交接面中点下标
                    double p0[3];
                    double p4[3];
                    double p8[3];
                    _get_centroid(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 2, 0), &verts1(k, j, i, 4, 0), &verts1(k, j, i, 6, 0), p0);
                    _get_centroid(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 1, 0), &verts1(k, j, i, 4, 0), &verts1(k, j, i, 5, 0), p4);
                    _get_centroid(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 1, 0), &verts1(k, j, i, 2, 0), &verts1(k, j, i, 3, 0), p8);
                    // 1,2,5,交接边点下标
                    double pp1[3];
                    double pp2[3];
                    double pp5[3];
                    _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 1, 0), pp1);
                    _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 2, 0), pp2);
                    _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 4, 0), pp5);
                    double n0[3];
                    double n4[3];
                    double n8[3];
                    double n48[3];
                    // 通过两个绝热面的流量为0,流量的方向可根据其法向量确定
                    _get_surface_normal(&verts1(k, j, i, 0, 0), pp2, p0, pp5, n0, Axis::XPOSITIVE);
                    _get_surface_normal(&verts1(k, j, i, 0, 0), pp1, p4, pp5, n4, Axis::YPOSITIVE);
                    _get_surface_normal(&verts1(k, j, i, 0, 0), pp1, p8, pp2, n8, Axis::ZPOSITIVE);
                    _cross_product(n4, n8, n48);
                    Eigen::Matrix3d matK;
                    _get_matK(pem1[cur], matK);
                    Eigen::RowVector3d Dx;
                    Eigen::Vector3d Dn;
                    Dx << p0[0] - _cx(k, j, i), p0[1] - _cy(k, j, i), p0[2] - _cz(k, j, i);
                    Dn << n48[0], n48[1], n48[2];
                    double cA = Dx * matK.inverse() * Dn;
                    cA = _dot_product(n0, n48) / cA;
                    Q -= cA * (p[cur] - Plow);
                }
                else if (j == ny && k == 0)
                {
                    cur -= nx;
                    // 1,4,11交接面中点下标
                    double p1[3];
                    double p4[3];
                    double p11[3];
                    _get_centroid(&verts1(k, j - 1, i, 2, 0), &verts1(k, j - 1, i, 3, 0), &verts1(k, j - 1, i, 6, 0), &verts1(k, j - 1, i, 7, 0), p4);
                    _get_centroid(&verts1(k, j - 1, i, 0, 0), &verts1(k, j - 1, i, 2, 0), &verts1(k, j - 1, i, 4, 0), &verts1(k, j - 1, i, 6, 0), p1);
                    _get_centroid(&verts1(k, j - 1, i, 0, 0), &verts1(k, j - 1, i, 1, 0), &verts1(k, j - 1, i, 2, 0), &verts1(k, j - 1, i, 3, 0), p11);
                    // 0,1,5,交接边点下标
                    double pp0[3];
                    double pp1[3];
                    double pp5[3];
                    _get_midpoint(&verts1(k, j - 1, i, 0, 0), &verts1(k, j - 1, i, 2, 0), pp0);
                    _get_midpoint(&verts1(k, j - 1, i, 2, 0), &verts1(k, j - 1, i, 3, 0), pp1);
                    _get_midpoint(&verts1(k, j - 1, i, 2, 0), &verts1(k, j - 1, i, 6, 0), pp5);
                    double n1[3];
                    double n4[3];
                    double n11[3];
                    double n411[3];
                    // 通过两个绝热面的流量为0,流量的方向可根据其法向量确定
                    _get_surface_normal(&verts1(k, j - 1, i, 2, 0), pp0, p1, pp5, n1, Axis::XPOSITIVE);
                    _get_surface_normal(&verts1(k, j - 1, i, 2, 0), pp1, p4, pp5, n4, Axis::YPOSITIVE);
                    _get_surface_normal(&verts1(k, j - 1, i, 2, 0), pp1, p11, pp0, n11, Axis::ZPOSITIVE);
                    _cross_product(n4, n11, n411);
                    Eigen::Matrix3d matK;
                    _get_matK(pem1[cur], matK);
                    Eigen::RowVector3d Dx;
                    Eigen::Vector3d Dn;
                    Dx << p1[0] - _cx(k, j - 1, i), p1[1] - _cy(k, j - 1, i), p1[2] - _cz(k, j - 1, i);
                    Dn << n411[0], n411[1], n411[2];
                    double cA = Dx * matK.inverse() * Dn;
                    cA = _dot_product(n1, n411) / cA;
                    Q -= cA * (p[cur] - Plow);
                }
                else if (j == 0 && k == nz)
                {
                    cur -= nx * ny;
                    // 3,5,8交接面中点下标
                    double p3[3];
                    double p5[3];
                    double p8[3];
                    _get_centroid(&verts1(k - 1, j, i, 0, 0), &verts1(k - 1, j, i, 2, 0), &verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 6, 0), p3);
                    _get_centroid(&verts1(k - 1, j, i, 0, 0), &verts1(k - 1, j, i, 1, 0), &verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 5, 0), p5);
                    _get_centroid(&verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 5, 0), &verts1(k - 1, j, i, 6, 0), &verts1(k - 1, j, i, 7, 0), p8);
                    // 1,2,4交接边点下标
                    double pp1[3];
                    double pp2[3];
                    double pp4[3];
                    _get_midpoint(&verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 5, 0), pp1);
                    _get_midpoint(&verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 6, 0), pp2);
                    _get_midpoint(&verts1(k - 1, j, i, 0, 0), &verts1(k - 1, j, i, 4, 0), pp4);
                    double n3[3];
                    double n5[3];
                    double n8[3];
                    double n58[3];
                    // 通过两个绝热面的流量为0,流量的方向可根据其法向量确定
                    _get_surface_normal(&verts1(k - 1, j, i, 4, 0), pp4, p3, pp2, n3, Axis::XPOSITIVE);
                    _get_surface_normal(&verts1(k - 1, j, i, 4, 0), pp1, p5, pp4, n5, Axis::YPOSITIVE);
                    _get_surface_normal(&verts1(k - 1, j, i, 4, 0), pp1, p8, pp2, n8, Axis::ZPOSITIVE);
                    _cross_product(n5, n8, n58);
                    Eigen::Matrix3d matK;
                    _get_matK(pem1[cur], matK);
                    Eigen::RowVector3d Dx;
                    Eigen::Vector3d Dn;
                    Dx << p3[0] - _cx(k - 1, j, i), p3[1] - _cy(k - 1, j, i), p3[2] - _cz(k - 1, j, i);
                    Dn << n58[0], n58[1], n58[2];
                    double cA = Dx * matK.inverse() * Dn;
                    cA = _dot_product(n3, n58) / cA;
                    Q -= cA * (p[cur] - Plow);
                }
                else if (j == ny && k == nz)
                {
                    cur -= nx * ny + nx;
                    // 2,5,11交接面中点下标
                    double p2[3];
                    double p5[3];
                    double p11[3];
                    _get_centroid(&verts1(k - 1, j - 1, i, 0, 0), &verts1(k - 1, j - 1, i, 2, 0), &verts1(k - 1, j - 1, i, 4, 0), &verts1(k - 1, j - 1, i, 6, 0), p2);
                    _get_centroid(&verts1(k - 1, j - 1, i, 2, 0), &verts1(k - 1, j - 1, i, 3, 0), &verts1(k - 1, j - 1, i, 6, 0), &verts1(k - 1, j - 1, i, 7, 0), p5);
                    _get_centroid(&verts1(k - 1, j - 1, i, 4, 0), &verts1(k - 1, j - 1, i, 5, 0), &verts1(k - 1, j - 1, i, 6, 0), &verts1(k - 1, j - 1, i, 7, 0), p11);
                    // 0,1,4交接边点下标
                    double pp0[3];
                    double pp1[3];
                    double pp4[3];
                    _get_midpoint(&verts1(k - 1, j - 1, i, 6, 0), &verts1(k - 1, j - 1, i, 4, 0), pp0);
                    _get_midpoint(&verts1(k - 1, j - 1, i, 6, 0), &verts1(k - 1, j - 1, i, 7, 0), pp1);
                    _get_midpoint(&verts1(k - 1, j - 1, i, 2, 0), &verts1(k - 1, j - 1, i, 6, 0), pp4);
                    double n2[3];
                    double n5[3];
                    double n11[3];
                    double n511[3];
                    // 通过两个绝热面的流量为0,流量的方向可根据其法向量确定
                    _get_surface_normal(&verts1(k - 1, j - 1, i, 6, 0), pp0, p2, pp4, n2, Axis::XPOSITIVE);
                    _get_surface_normal(&verts1(k - 1, j - 1, i, 6, 0), pp1, p5, pp4, n5, Axis::YPOSITIVE);
                    _get_surface_normal(&verts1(k - 1, j - 1, i, 6, 0), pp1, p11, pp0, n11, Axis::ZPOSITIVE);
                    _cross_product(n5, n11, n511);
                    Eigen::Matrix3d matK;
                    _get_matK(pem1[cur], matK);
                    Eigen::RowVector3d Dx;
                    Eigen::Vector3d Dn;
                    Dx << p2[0] - _cx(k - 1, j - 1, i), p2[1] - _cy(k - 1, j - 1, i), p2[2] - _cz(k - 1, j - 1, i);
                    Dn << n511[0], n511[1], n511[2];
                    double cA = Dx * matK.inverse() * Dn;
                    cA = _dot_product(n2, n511) / cA;
                    Q -= cA * (p[cur] - Plow);
                }
                else if (j == 0)
                {
                    int i6 = cur, i2 = cur - nx * ny;
                    // 0,3,4,5,8交接面中点下标
                    double p0[3];
                    double p3[3];
                    double p4[3];
                    double p5[3];
                    double p8[3];
                    _get_centroid(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 2, 0), &verts1(k, j, i, 4, 0), &verts1(k, j, i, 6, 0), p0);
                    _get_centroid(&verts1(k - 1, j, i, 0, 0), &verts1(k - 1, j, i, 2, 0), &verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 6, 0), p3);
                    _get_centroid(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 1, 0), &verts1(k, j, i, 4, 0), &verts1(k, j, i, 5, 0), p4);
                    _get_centroid(&verts1(k - 1, j, i, 0, 0), &verts1(k - 1, j, i, 1, 0), &verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 5, 0), p5);
                    _get_centroid(&verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 5, 0), &verts1(k - 1, j, i, 6, 0), &verts1(k - 1, j, i, 7, 0), p8);
                    // 1,2,4,5交接边点下标
                    double pp1[3];
                    double pp2[3];
                    double pp4[3];
                    double pp5[3];
                    _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 1, 0), pp1);
                    _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 2, 0), pp2);
                    _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 4, 0), pp5);
                    _get_midpoint(&verts1(k - 1, j, i, 0, 0), &verts1(k - 1, j, i, 4, 0), pp4);
                    double n0[3];
                    double n3[3];
                    double n4[3];
                    double n5[3];
                    double n8[3];
                    _get_surface_normal(&verts1(k, j, i, 0, 0), pp2, p0, pp5, n0, Axis::XPOSITIVE);
                    _get_surface_normal(&verts1(k - 1, j, i, 4, 0), pp2, p3, pp4, n3, Axis::XPOSITIVE);
                    _get_surface_normal(&verts1(k, j, i, 0, 0), pp5, p4, pp1, n4, Axis::YPOSITIVE);
                    _get_surface_normal(&verts1(k - 1, j, i, 4, 0), pp4, p5, pp1, n5, Axis::YPOSITIVE);
                    _get_surface_normal(&verts1(k - 1, j, i, 4, 0), pp1, p8, pp2, n8, Axis::ZPOSITIVE);
                    Eigen::Matrix<double, 6, 6> A;
                    A.setZero();
                    A.row(0) << 0, 0, 0, p0[0] - _cx(k, j, i), p0[1] - _cy(k, j, i), p0[2] - _cz(k, j, i);
                    A.row(1) << p3[0] - _cx(k - 1, j, i), p3[1] - _cy(k - 1, j, i), p3[2] - _cz(k - 1, j, i), 0, 0, 0;
                    A.row(2) << _cx(k - 1, j, i) - p8[0], _cy(k - 1, j, i) - p8[1], _cz(k - 1, j, i) - p8[2], p8[0] - _cx(k, j, i), p8[1] - _cy(k, j, i), p8[2] - _cz(k, j, i);
                    Eigen::Matrix3d matK6, matK2;
                    _get_matK(pem1[i6], matK6);
                    _get_matK(pem1[i2], matK2);
                    Eigen::RowVector3d v0, v3, v4, v5, v8;
                    v0 << n0[0], n0[1], n0[2];
                    v3 << n3[0], n3[1], n3[2];
                    v4 << n4[0], n4[1], n4[2];
                    v5 << n5[0], n5[1], n5[2];
                    v8 << n8[0], n8[1], n8[2];
                    A.block(3, 0, 1, 3) = v8 * matK2;
                    A.block(3, 3, 1, 3) = -v8 * matK6;
                    A.block(4, 3, 1, 3) = v4 * matK6;
                    A.block(5, 0, 1, 3) = v5 * matK2;
                    A = A.inverse();
                    Eigen::RowVectorXd r = (v8 - v3) * matK2 * A.topRows(3);
                    Q -= Plow * (r[0] + r[1]);
                    Q -= (-r[1] + r[2]) * p[i2];
                    Q -= (-r[0] - r[2]) * p[i6];
                    r = -(v0 + v8) * matK6 * A.bottomRows(3);
                    Q -= Plow * (r[0] + r[1]);
                    Q -= (-r[1] + r[2]) * p[i2];
                    Q -= (-r[0] - r[2]) * p[i6];
                }
                else if (j == ny)
                {
                    int i5 = cur - nx, i1 = i5 - nx * ny;
                    // 1,2,4,5,11交接面中点下标
                    double p1[3];
                    double p2[3];
                    double p4[3];
                    double p5[3];
                    double p11[3];
                    _get_centroid(&verts1(k, j - 1, i, 0, 0), &verts1(k, j - 1, i, 2, 0), &verts1(k, j - 1, i, 4, 0), &verts1(k, j - 1, i, 6, 0), p1);
                    _get_centroid(&verts1(k - 1, j - 1, i, 0, 0), &verts1(k - 1, j - 1, i, 2, 0), &verts1(k - 1, j - 1, i, 4, 0), &verts1(k - 1, j - 1, i, 6, 0), p2);
                    _get_centroid(&verts1(k, j - 1, i, 2, 0), &verts1(k, j - 1, i, 3, 0), &verts1(k, j - 1, i, 6, 0), &verts1(k, j - 1, i, 7, 0), p4);
                    _get_centroid(&verts1(k - 1, j - 1, i, 2, 0), &verts1(k - 1, j - 1, i, 3, 0), &verts1(k - 1, j - 1, i, 6, 0), &verts1(k - 1, j - 1, i, 7, 0), p5);
                    _get_centroid(&verts1(k - 1, j - 1, i, 4, 0), &verts1(k - 1, j - 1, i, 5, 0), &verts1(k - 1, j - 1, i, 6, 0), &verts1(k - 1, j - 1, i, 7, 0), p11);
                    // 0,1,4,5交接边点下标
                    double pp0[3];
                    double pp1[3];
                    double pp4[3];
                    double pp5[3];
                    _get_midpoint(&verts1(k, j - 1, i, 0, 0), &verts1(k, j - 1, i, 2, 0), pp0);
                    _get_midpoint(&verts1(k, j - 1, i, 3, 0), &verts1(k, j - 1, i, 2, 0), pp1);
                    _get_midpoint(&verts1(k, j - 1, i, 6, 0), &verts1(k, j - 1, i, 2, 0), pp5);
                    _get_midpoint(&verts1(k - 1, j - 1, i, 2, 0), &verts1(k - 1, j - 1, i, 6, 0), pp4);
                    double n1[3];
                    double n2[3];
                    double n4[3];
                    double n5[3];
                    double n11[3];
                    _get_surface_normal(&verts1(k, j - 1, i, 2, 0), pp0, p1, pp5, n1, Axis::XPOSITIVE);
                    _get_surface_normal(&verts1(k - 1, j - 1, i, 6, 0), pp4, p2, pp0, n2, Axis::XPOSITIVE);
                    _get_surface_normal(&verts1(k, j - 1, i, 2, 0), pp1, p4, pp5, n4, Axis::YPOSITIVE);
                    _get_surface_normal(&verts1(k - 1, j - 1, i, 6, 0), pp4, p5, pp1, n5, Axis::YPOSITIVE);
                    _get_surface_normal(&verts1(k - 1, j - 1, i, 6, 0), pp0, p11, pp1, n11, Axis::ZPOSITIVE);
                    Eigen::Matrix<double, 6, 6> A;
                    A.setZero();
                    A.row(0) << 0, 0, 0, p1[0] - _cx(k, j - 1, i), p1[1] - _cy(k, j - 1, i), p1[2] - _cz(k, j - 1, i);
                    A.row(1) << p2[0] - _cx(k - 1, j - 1, i), p2[1] - _cy(k - 1, j - 1, i), p2[2] - _cz(k - 1, j - 1, i), 0, 0, 0;
                    A.row(2) << _cx(k - 1, j - 1, i) - p11[0], _cy(k - 1, j - 1, i) - p11[1], _cz(k - 1, j - 1, i) - p11[2], p11[0] - _cx(k, j - 1, i), p11[1] - _cy(k, j - 1, i), p11[2] - _cz(k, j - 1, i);
                    Eigen::Matrix3d matK5, matK1;
                    _get_matK(pem1[i5], matK5);
                    _get_matK(pem1[i1], matK1);
                    Eigen::RowVector3d v1, v2, v4, v5, v11;
                    v1 << n1[0], n1[1], n1[2];
                    v2 << n2[0], n2[1], n2[2];
                    v4 << n4[0], n4[1], n4[2];
                    v5 << n5[0], n5[1], n5[2];
                    v11 << n11[0], n11[1], n11[2];
                    A.block(3, 0, 1, 3) = v11 * matK1;
                    A.block(3, 3, 1, 3) = -v11 * matK5;
                    A.block(4, 3, 1, 3) = v4 * matK5;
                    A.block(5, 0, 1, 3) = v5 * matK1;
                    A = A.inverse();
                    Eigen::RowVectorXd r = (v11 - v2) * matK1 * A.topRows(3);
                    Q -= Plow * (r[0] + r[1]);
                    Q -= (-r[1] + r[2]) * p[i1];
                    Q -= (-r[0] - r[2]) * p[i5];
                    r = -(v1 + v11) * matK5 * A.bottomRows(3);
                    Q -= Plow * (r[0] + r[1]);
                    Q -= (-r[1] + r[2]) * p[i1];
                    Q -= (-r[0] - r[2]) * p[i5];
                }
                else if (k == 0)
                {
                    // 两个网格中心下标
                    int i6 = cur, i5 = cur - nx;
                    // 0,1,4,8,11交接面中点下标
                    double p0[3];
                    double p1[3];
                    double p4[3];
                    double p8[3];
                    double p11[3];
                    _get_centroid(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 2, 0), &verts1(k, j, i, 4, 0), &verts1(k, j, i, 6, 0), p0);
                    _get_centroid(&verts1(k, j - 1, i, 0, 0), &verts1(k, j - 1, i, 2, 0), &verts1(k, j - 1, i, 4, 0), &verts1(k, j - 1, i, 6, 0), p1);
                    _get_centroid(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 1, 0), &verts1(k, j, i, 4, 0), &verts1(k, j, i, 5, 0), p4);
                    _get_centroid(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 1, 0), &verts1(k, j, i, 2, 0), &verts1(k, j, i, 3, 0), p8);
                    _get_centroid(&verts1(k, j - 1, i, 0, 0), &verts1(k, j - 1, i, 1, 0), &verts1(k, j - 1, i, 2, 0), &verts1(k, j - 1, i, 3, 0), p11);
                    // 0,1,2,5交接边点下标
                    double pp0[3];
                    double pp1[3];
                    double pp2[3];
                    double pp5[3];
                    _get_midpoint(&verts1(k, j - 1, i, 0, 0), &verts1(k, j - 1, i, 2, 0), pp0);
                    _get_midpoint(&verts1(k, j - 1, i, 2, 0), &verts1(k, j - 1, i, 3, 0), pp1);
                    _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 2, 0), pp2);
                    _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 4, 0), pp5);
                    double n0[3];
                    double n1[3];
                    double n4[3];
                    double n8[3];
                    double n11[3];
                    _get_surface_normal(&verts1(k, j, i, 0, 0), pp2, p0, pp5, n0, Axis::XPOSITIVE);
                    _get_surface_normal(&verts1(k, j - 1, i, 2, 0), pp0, p1, pp5, n1, Axis::XPOSITIVE);
                    _get_surface_normal(&verts1(k, j, i, 0, 0), pp5, p4, pp1, n4, Axis::YPOSITIVE);
                    _get_surface_normal(&verts1(k, j, i, 0, 0), pp1, p8, pp2, n8, Axis::ZPOSITIVE);
                    _get_surface_normal(&verts1(k, j - 1, i, 2, 0), pp0, p11, pp1, n11, Axis::ZPOSITIVE);
                    Eigen::Matrix<double, 6, 6> A;
                    A.setZero();
                    A.row(0) << 0, 0, 0, p0[0] - _cx(k, j, i), p0[1] - _cy(k, j, i), p0[2] - _cz(k, j, i);
                    A.row(1) << p1[0] - _cx(k, j - 1, i), p1[1] - _cy(k, j - 1, i), p1[2] - _cz(k, j - 1, i), 0, 0, 0;
                    A.row(2) << _cx(k, j - 1, i) - p4[0], _cy(k, j - 1, i) - p4[1], _cz(k, j - 1, i) - p4[2], p4[0] - _cx(k, j, i), p4[1] - _cy(k, j, i), p4[2] - _cz(k, j, i);
                    Eigen::Matrix3d matK6, matK5;
                    _get_matK(pem1[i6], matK6);
                    _get_matK(pem1[i5], matK5);
                    Eigen::RowVector3d v0, v1, v4, v8, v11;
                    v0 << n0[0], n0[1], n0[2];
                    v1 << n1[0], n1[1], n1[2];
                    v4 << n4[0], n4[1], n4[2];
                    v8 << n8[0], n8[1], n8[2];
                    v11 << n11[0], n11[1], n11[2];
                    A.block(3, 0, 1, 3) = v4 * matK5;
                    A.block(3, 3, 1, 3) = -v4 * matK6;
                    A.block(4, 3, 1, 3) = v8 * matK6;
                    A.block(5, 0, 1, 3) = v11 * matK5;
                    A = A.inverse();
                    Eigen::RowVectorXd r = (v4 - v1) * matK5 * A.topRows(3);
                    Q -= Plow * (r[0] + r[1]);
                    Q -= (-r[1] + r[2]) * p[i5];
                    Q -= (-r[0] - r[2]) * p[i6];
                    r = -(v0 + v4) * matK6 * A.bottomRows(3);
                    Q -= Plow * (r[0] + r[1]);
                    Q -= (-r[1] + r[2]) * p[i5];
                    Q -= (-r[0] - r[2]) * p[i6];
                }
                else if (k == nz)
                {
                    // 两个网格中心下标
                    int i2 = cur - nx * ny, i1 = i2 - nx;
                    // 2,3,5,8,11交接面中点下标
                    double p2[3];
                    double p3[3];
                    double p5[3];
                    double p8[3];
                    double p11[3];
                    _get_centroid(&verts1(k - 1, j - 1, i, 0, 0), &verts1(k - 1, j - 1, i, 2, 0), &verts1(k - 1, j - 1, i, 4, 0), &verts1(k - 1, j - 1, i, 6, 0), p2);
                    _get_centroid(&verts1(k - 1, j, i, 0, 0), &verts1(k - 1, j, i, 2, 0), &verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 6, 0), p3);
                    _get_centroid(&verts1(k - 1, j, i, 0, 0), &verts1(k - 1, j, i, 1, 0), &verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 5, 0), p5);
                    _get_centroid(&verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 5, 0), &verts1(k - 1, j, i, 6, 0), &verts1(k - 1, j, i, 7, 0), p8);
                    _get_centroid(&verts1(k - 1, j - 1, i, 4, 0), &verts1(k - 1, j - 1, i, 5, 0), &verts1(k - 1, j - 1, i, 6, 0), &verts1(k - 1, j - 1, i, 7, 0), p11);
                    // 0,1,2,4交接边点下标
                    double pp0[3];
                    double pp1[3];
                    double pp2[3];
                    double pp4[3];
                    _get_midpoint(&verts1(k - 1, j - 1, i, 4, 0), &verts1(k - 1, j - 1, i, 6, 0), pp0);
                    _get_midpoint(&verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 6, 0), pp2);
                    _get_midpoint(&verts1(k - 1, j - 1, i, 6, 0), &verts1(k - 1, j - 1, i, 7, 0), pp1);
                    _get_midpoint(&verts1(k - 1, j, i, 0, 0), &verts1(k - 1, j, i, 4, 0), pp4);
                    double n2[3];
                    double n3[3];
                    double n5[3];
                    double n8[3];
                    double n11[3];
                    _get_surface_normal(&verts1(k - 1, j - 1, i, 6, 0), pp0, p2, pp4, n2, Axis::XPOSITIVE);
                    _get_surface_normal(&verts1(k - 1, j, i, 4, 0), pp2, p3, pp4, n3, Axis::XPOSITIVE);
                    _get_surface_normal(&verts1(k - 1, j, i, 4, 0), pp4, p5, pp1, n5, Axis::YPOSITIVE);
                    _get_surface_normal(&verts1(k - 1, j, i, 4, 0), pp1, p8, pp2, n8, Axis::ZPOSITIVE);
                    _get_surface_normal(&verts1(k - 1, j - 1, i, 6, 0), pp0, p11, pp1, n11, Axis::ZPOSITIVE);
                    Eigen::Matrix<double, 6, 6> A;
                    A.setZero();
                    A.row(0) << 0, 0, 0, p3[0] - _cx(k - 1, j, i), p3[1] - _cy(k - 1, j, i), p3[2] - _cz(k - 1, j, i);
                    A.row(1) << p2[0] - _cx(k - 1, j - 1, i), p2[1] - _cy(k - 1, j - 1, i), p2[2] - _cz(k - 1, j - 1, i), 0, 0, 0;
                    A.row(2) << _cx(k - 1, j - 1, i) - p5[0], _cy(k - 1, j - 1, i) - p5[1], _cz(k - 1, j - 1, i) - p5[2], p5[0] - _cx(k - 1, j, i), p5[1] - _cy(k - 1, j, i), p5[2] - _cz(k - 1, j, i);
                    Eigen::Matrix3d matK2, matK1;
                    _get_matK(pem1[i2], matK2);
                    _get_matK(pem1[i1], matK1);
                    Eigen::RowVector3d v2, v3, v5, v8, v11;
                    v2 << n2[0], n2[1], n2[2];
                    v3 << n3[0], n3[1], n3[2];
                    v5 << n5[0], n5[1], n5[2];
                    v8 << n8[0], n8[1], n8[2];
                    v11 << n11[0], n11[1], n11[2];
                    A.block(3, 0, 1, 3) = v5 * matK1;
                    A.block(3, 3, 1, 3) = -v5 * matK2;
                    A.block(4, 3, 1, 3) = v8 * matK2;
                    A.block(5, 0, 1, 3) = v11 * matK1;
                    A = A.inverse();
                    Eigen::RowVectorXd r = (v5 - v2) * matK1 * A.topRows(3);
                    Q -= Plow * (r[0] + r[1]);
                    Q -= (-r[1] + r[2]) * p[i1];
                    Q -= (-r[0] - r[2]) * p[i2];
                    r = -(v3 + v5) * matK2 * A.bottomRows(3);
                    Q -= Plow * (r[0] + r[1]);
                    Q -= (-r[1] + r[2]) * p[i1];
                    Q -= (-r[0] - r[2]) * p[i2];
                }
                else
                {
                    int i6 = cur, i5 = i6 - nx, i2 = i6 - nx * ny, i1 = i5 - nx * ny;
                    // 0,1,2,3,4,5,8,11交接面中点下标
                    double p0[3];
                    double p1[3];
                    double p2[3];
                    double p3[3];
                    double p4[3];
                    double p5[3];
                    double p8[3];
                    double p11[3];
                    _get_centroid(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 1, 0), &verts1(k, j, i, 2, 0), &verts1(k, j, i, 3, 0), p8);
                    _get_centroid(&verts1(k, j - 1, i, 0, 0), &verts1(k, j - 1, i, 1, 0), &verts1(k, j - 1, i, 2, 0), &verts1(k, j - 1, i, 3, 0), p11);
                    _get_centroid(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 1, 0), &verts1(k, j, i, 4, 0), &verts1(k, j, i, 5, 0), p4);
                    _get_centroid(&verts1(k - 1, j, i, 0, 0), &verts1(k - 1, j, i, 1, 0), &verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 5, 0), p5);
                    _get_centroid(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 2, 0), &verts1(k, j, i, 4, 0), &verts1(k, j, i, 6, 0), p0);
                    _get_centroid(&verts1(k - 1, j, i, 0, 0), &verts1(k - 1, j, i, 2, 0), &verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 6, 0), p3);
                    _get_centroid(&verts1(k, j - 1, i, 0, 0), &verts1(k, j - 1, i, 2, 0), &verts1(k, j - 1, i, 4, 0), &verts1(k, j - 1, i, 6, 0), p1);
                    _get_centroid(&verts1(k - 1, j - 1, i, 0, 0), &verts1(k - 1, j - 1, i, 2, 0), &verts1(k - 1, j - 1, i, 4, 0), &verts1(k - 1, j - 1, i, 6, 0), p2);
                    // 交接边点下标0,1,2,4,5
                    double pp0[3];
                    double pp1[3];
                    double pp2[3];
                    double pp4[3];
                    double pp5[3];
                    _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 1, 0), pp1);
                    _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 2, 0), pp2);
                    _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 4, 0), pp5);
                    _get_midpoint(&verts1(k - 1, j, i, 0, 0), &verts1(k - 1, j, i, 4, 0), pp4);
                    _get_midpoint(&verts1(k, j - 1, i, 0, 0), &verts1(k, j - 1, i, 2, 0), pp0);
                    double n0[3];
                    double n1[3];
                    double n2[3];
                    double n3[3];
                    double n4[3];
                    double n5[3];
                    double n8[3];
                    double n11[3];
                    _get_surface_normal(&verts1(k, j, i, 0, 0), pp2, p0, pp5, n0, Axis::XPOSITIVE);
                    _get_surface_normal(&verts1(k, j, i, 0, 0), pp0, p1, pp5, n1, Axis::XPOSITIVE);
                    _get_surface_normal(&verts1(k, j, i, 0, 0), pp0, p2, pp4, n2, Axis::XPOSITIVE);
                    _get_surface_normal(&verts1(k, j, i, 0, 0), pp4, p3, pp2, n3, Axis::XPOSITIVE);
                    _get_surface_normal(&verts1(k, j, i, 0, 0), pp1, p4, pp5, n4, Axis::YPOSITIVE);
                    _get_surface_normal(&verts1(k, j, i, 0, 0), pp4, p5, pp1, n5, Axis::YPOSITIVE);
                    _get_surface_normal(&verts1(k, j, i, 0, 0), pp2, p8, pp1, n8, Axis::ZPOSITIVE);
                    _get_surface_normal(&verts1(k, j, i, 0, 0), pp0, p11, pp1, n11, Axis::ZPOSITIVE);
                    Eigen::Matrix<double, 12, 12> A;
                    A.setZero();
                    A.block(0, 9, 1, 3) << p0[0] - _cx(k, j, i), p0[1] - _cy(k, j, i), p0[2] - _cz(k, j, i);
                    A.block(1, 6, 1, 3) << p1[0] - _cx(k, j - 1, i), p1[1] - _cy(k, j - 1, i), p1[2] - _cz(k, j - 1, i);
                    A.block(2, 0, 1, 3) << p2[0] - _cx(k - 1, j - 1, i), p2[1] - _cy(k - 1, j - 1, i), p2[2] - _cz(k - 1, j - 1, i);
                    A.block(3, 3, 1, 3) << p3[0] - _cx(k - 1, j, i), p3[1] - _cy(k - 1, j, i), p3[2] - _cz(k - 1, j, i);
                    Eigen::Matrix3d matK6, matK5, matK2, matK1;
                    _get_matK(pem1[i6], matK6);
                    _get_matK(pem1[i5], matK5);
                    _get_matK(pem1[i2], matK2);
                    _get_matK(pem1[i1], matK1);
                    Eigen::RowVector3d v0, v1, v2, v3, v4, v5, v8, v11;
                    v0 << n0[0], n0[1], n0[2];
                    v1 << n1[0], n1[1], n1[2];
                    v2 << n2[0], n2[1], n2[2];
                    v3 << n3[0], n3[1], n3[2];
                    v4 << n4[0], n4[1], n4[2];
                    v5 << n5[0], n5[1], n5[2];
                    v8 << n8[0], n8[1], n8[2];
                    v11 << n11[0], n11[1], n11[2];
                    A.block(4, 6, 1, 6) << _cx(k, j - 1, i) - p4[0], _cy(k, j - 1, i) - p4[1], _cz(k, j - 1, i) - p4[2], p4[0] - _cx(k, j, i), p4[1] - _cy(k, j, i), p4[2] - _cz(k, j, i);
                    A.block(5, 6, 1, 3) = v4 * matK5;
                    A.block(5, 9, 1, 3) = -v4 * matK6;
                    A.block(6, 0, 1, 6) << _cx(k - 1, j - 1, i) - p5[0], _cy(k - 1, j - 1, i) - p5[1], _cz(k - 1, j - 1, i) - p5[2], p5[0] - _cx(k - 1, j, i), p5[1] - _cy(k - 1, j, i), p5[2] - _cz(k - 1, j, i);
                    A.block(7, 0, 1, 3) = v5 * matK1;
                    A.block(7, 3, 1, 3) = -v5 * matK2;
                    A.block(8, 3, 1, 3) << _cx(k - 1, j, i) - p8[0], _cy(k - 1, j, i) - p8[1], _cz(k - 1, j, i) - p8[2];
                    A.block(8, 9, 1, 3) << p8[0] - _cx(k, j, i), p8[1] - _cy(k, j, i), p8[2] - _cz(k, j, i);
                    A.block(9, 3, 1, 3) = v8 * matK2;
                    A.block(9, 9, 1, 3) = -v8 * matK6;
                    A.block(10, 0, 1, 3) << _cx(k - 1, j - 1, i) - p11[0], _cy(k - 1, j - 1, i) - p11[1], _cz(k - 1, j - 1, i) - p11[2];
                    A.block(10, 6, 1, 3) << p11[0] - _cx(k, j - 1, i), p11[1] - _cy(k, j - 1, i), p11[2] - _cz(k, j - 1, i);
                    A.block(11, 0, 1, 3) = v11 * matK1;
                    A.block(11, 6, 1, 3) = -v11 * matK5;
                    A = A.inverse();
                    Eigen::RowVectorXd r = -v2 * matK1 * A.topRows(3);
                    Q -= Plow * (r[0] + r[1] + r[2] + r[3]);
                    Q -= (-r[2] + r[6] + r[10]) * p[i1];
                    Q -= (-r[3] - r[6] + r[8]) * p[i2];
                    Q -= (-r[1] + r[4] - r[10]) * p[i5];
                    Q -= (-r[0] - r[4] - r[8]) * p[i6];
                    r = -v3 * matK2 * A.block(3, 0, 3, 12);
                    Q -= Plow * (r[0] + r[1] + r[2] + r[3]);
                    Q -= (-r[2] + r[6] + r[10]) * p[i1];
                    Q -= (-r[3] - r[6] + r[8]) * p[i2];
                    Q -= (-r[1] + r[4] - r[10]) * p[i5];
                    Q -= (-r[0] - r[4] - r[8]) * p[i6];
                    r = -v1 * matK5 * A.block(6, 0, 3, 12);
                    Q -= Plow * (r[0] + r[1] + r[2] + r[3]);
                    Q -= (-r[2] + r[6] + r[10]) * p[i1];
                    Q -= (-r[3] - r[6] + r[8]) * p[i2];
                    Q -= (-r[1] + r[4] - r[10]) * p[i5];
                    Q -= (-r[0] - r[4] - r[8]) * p[i6];
                    r = -v0 * matK6 * A.bottomRows(3);
                    Q -= Plow * (r[0] + r[1] + r[2] + r[3]);
                    Q -= (-r[2] + r[6] + r[10]) * p[i1];
                    Q -= (-r[3] - r[6] + r[8]) * p[i2];
                    Q -= (-r[1] + r[4] - r[10]) * p[i5];
                    Q -= (-r[0] - r[4] - r[8]) * p[i6];
                }
            }
    std::cout << "Q= " << Q << std::endl;
    // for (iter = 0; iter < 100; iter++)
    // {
    //     n = 0;
    //     nnz = 0;
    //     Ptr[0] = 0;

    //     /*for (int i = 0; i < 27; i++) {
    //         std::cout << beta[i][0] << " " << beta[i][1] << " " << beta[i][2] << "\n";
    //     }*/

    //     std::cout << "iter=" << iter << "\n";
    //     for (k = 0; k < nz; k++)
    //     {
    //         for (j = 0; j < ny; j++)
    //         {
    //             for (i = 0; i < nx; i++)
    //             {
    //                 // i = 1; j = 1; k = 1;
    //                 // std::cout << "(" << i << "," << j << "," << k << ") : ";
    //                 count = 0;
    //                 double val = 0.0;
    //                 double F = 0.0;
    //                 double F1 = 0.0;
    //                 int i0 = k * nx * ny + j * nx + i;

    //                 _get_F(pem1, p, i, j, k, _verts1, F, divn);
    //                 // std::cout << "F=" << F << std::endl;

    //                 B[n] = -F;
    //                 if (k > 0)
    //                 {
    //                     _get_csr(pem1, p, i, j, k, _verts1, F, i, j, k - 1, Idx, Val, nnz, count, divn);

    //                     if (i > 0)
    //                     {
    //                         _get_csr(pem1, p, i, j, k, _verts1, F, i - 1, j, k - 1, Idx, Val, nnz, count, divn);

    //                         if (j > 0)
    //                         {
    //                             _get_csr(pem1, p, i, j, k, _verts1, F, i - 1, j - 1, k - 1, Idx, Val, nnz, count, divn);
    //                         }
    //                         if (j < ny - 1)
    //                         {
    //                             _get_csr(pem1, p, i, j, k, _verts1, F, i - 1, j + 1, k - 1, Idx, Val, nnz, count, divn);
    //                         }
    //                     }
    //                     if (i < nx - 1)
    //                     {
    //                         _get_csr(pem1, p, i, j, k, _verts1, F, i + 1, j, k - 1, Idx, Val, nnz, count, divn);

    //                         if (j > 0)
    //                         {
    //                             _get_csr(pem1, p, i, j, k, _verts1, F, i + 1, j - 1, k - 1, Idx, Val, nnz, count, divn);
    //                         }
    //                         if (j < ny - 1)
    //                         {
    //                             _get_csr(pem1, p, i, j, k, _verts1, F, i + 1, j + 1, k - 1, Idx, Val, nnz, count, divn);
    //                         }
    //                     }
    //                     if (j > 0)
    //                     {
    //                         _get_csr(pem1, p, i, j, k, _verts1, F, i, j - 1, k - 1, Idx, Val, nnz, count, divn);
    //                     }
    //                     if (j < ny - 1)
    //                     {
    //                         _get_csr(pem1, p, i, j, k, _verts1, F, i, j + 1, k - 1, Idx, Val, nnz, count, divn);
    //                     }
    //                 }

    //                 if (j > 0)
    //                 {
    //                     _get_csr(pem1, p, i, j, k, _verts1, F, i, j - 1, k, Idx, Val, nnz, count, divn);
    //                 }
    //                 if (i > 0)
    //                 {
    //                     _get_csr(pem1, p, i, j, k, _verts1, F, i - 1, j, k, Idx, Val, nnz, count, divn);

    //                     if (j > 0)
    //                     {
    //                         _get_csr(pem1, p, i, j, k, _verts1, F, i - 1, j - 1, k, Idx, Val, nnz, count, divn);
    //                     }
    //                     if (j < ny - 1)
    //                     {
    //                         _get_csr(pem1, p, i, j, k, _verts1, F, i - 1, j + 1, k, Idx, Val, nnz, count, divn);
    //                     }
    //                 }
    //                 _get_csr(pem1, p, i, j, k, _verts1, F, i, j, k, Idx, Val, nnz, count, divn);
    //                 if (i < nx - 1)
    //                 {
    //                     _get_csr(pem1, p, i, j, k, _verts1, F, i + 1, j, k, Idx, Val, nnz, count, divn);

    //                     if (j > 0)
    //                     {
    //                         _get_csr(pem1, p, i, j, k, _verts1, F, i + 1, j - 1, k, Idx, Val, nnz, count, divn);
    //                     }
    //                     if (j < ny - 1)
    //                     {
    //                         _get_csr(pem1, p, i, j, k, _verts1, F, i + 1, j + 1, k, Idx, Val, nnz, count, divn);
    //                     }
    //                 }
    //                 if (j < ny - 1)
    //                 {
    //                     _get_csr(pem1, p, i, j, k, _verts1, F, i, j + 1, k, Idx, Val, nnz, count, divn);
    //                 }

    //                 if (k < nz - 1)
    //                 {
    //                     _get_csr(pem1, p, i, j, k, _verts1, F, i, j, k + 1, Idx, Val, nnz, count, divn);

    //                     if (i > 0)
    //                     {
    //                         _get_csr(pem1, p, i, j, k, _verts1, F, i - 1, j, k + 1, Idx, Val, nnz, count, divn);

    //                         if (j > 0)
    //                         {
    //                             _get_csr(pem1, p, i, j, k, _verts1, F, i - 1, j - 1, k + 1, Idx, Val, nnz, count, divn);
    //                         }
    //                         if (j < ny - 1)
    //                         {
    //                             _get_csr(pem1, p, i, j, k, _verts1, F, i - 1, j + 1, k + 1, Idx, Val, nnz, count, divn);
    //                         }
    //                     }
    //                     if (i < nx - 1)
    //                     {
    //                         _get_csr(pem1, p, i, j, k, _verts1, F, i + 1, j, k + 1, Idx, Val, nnz, count, divn);

    //                         if (j > 0)
    //                         {
    //                             _get_csr(pem1, p, i, j, k, _verts1, F, i + 1, j - 1, k + 1, Idx, Val, nnz, count, divn);
    //                         }
    //                         if (j < ny - 1)
    //                         {
    //                             _get_csr(pem1, p, i, j, k, _verts1, F, i + 1, j + 1, k + 1, Idx, Val, nnz, count, divn);
    //                         }
    //                     }
    //                     if (j > 0)
    //                     {
    //                         _get_csr(pem1, p, i, j, k, _verts1, F, i, j - 1, k + 1, Idx, Val, nnz, count, divn);
    //                     }
    //                     if (j < ny - 1)
    //                     {
    //                         _get_csr(pem1, p, i, j, k, _verts1, F, i, j + 1, k + 1, Idx, Val, nnz, count, divn);
    //                     }
    //                 }
    //                 Ptr[n + 1] = Ptr[n] + count;

    //                 n = n + 1;
    //                 // std::cout << std::endl;
    //             }
    //         }
    //     }

    //     // FILE* fp = std::fopen("Val_MPFA.INC", "w");
    //     // for (int ii = 0; ii < nnz; ii++) {
    //     //     std::fprintf(fp, "%le\n", Val[ii]);
    //     // }
    //     // std::fclose(fp);

    //     pmgmres_ilu_cr(nx * ny * nz, nnz, Ptr, Idx, Val, x, B, 1000, 1000, 1e-5, 1e-5);

    //     double res = 0.0;

    //     for (int i = 0; i < nx * ny * nz; ++i)
    //     {
    //         if (std::abs(x[i]) > res)
    //             res = std::abs(x[i]);
    //     }

    //     for (int i = 0; i < nx * ny * nz; ++i)
    //     {
    //         p[i] += x[i];
    //     }

    //     Q = 0.0;
    //     for (int k = 0; k < nz; k++)
    //     {
    //         for (int j = 0; j < ny; j++)
    //         {
    //             const double *pp[4];
    //             double S = 0.0, D = 0.0;
    //             double pc[3];
    //             int i0 = k * nx * ny + j * nx + nx - 1;

    //             pp[0] = &verts1(k, j, nx - 1, 1, 0);
    //             pp[1] = &verts1(k, j, nx - 1, 3, 0);
    //             pp[2] = &verts1(k, j, nx - 1, 5, 0);
    //             pp[3] = &verts1(k, j, nx - 1, 7, 0);

    //             pc[0] = _cx(k, j, nx - 1);
    //             pc[1] = _cy(k, j, nx - 1);
    //             pc[2] = _cz(k, j, nx - 1);

    //             S = _get_area(pp[0], pp[1], pp[2]) + _get_area(pp[3], pp[1], pp[2]);
    //             D = _get_distance(pc, pp[1], pp[2], pp[3]);
    //             Q = Q + (Phigh - p[i0]) * S / D * pem1[i0].x;
    //         }
    //     }

    //     Iter = iter;
    //     std::cout << "Q=" << Q << "\n";
    //     std::cout << "res=" << res << "\n";
    //     if (res < 1.0e-3)
    //         break;
    // }
    delete[] pem1, p, Ptr, Idx, Val, B, _verts1, cx, cy, cz;
    std::cout << std::endl;
}