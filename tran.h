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

void _get_Qx(const Ktensor *perm, const double *verts, int divn, int &Iter, double &Q)
{
    std::cout << "divn = " << divn << std::endl;
    int nx = nxx * divn;
    int ny = nyy * divn;
    int nz = nzz * divn;
    int n, nnz, count;
    int i, j, k, iter;

    double *_verts1 = new double[nx * ny * nz * 8 * 3];
    std::mdspan<double, std::extents<size_t, -1, -1, -1, 8, 3>> verts1(_verts1, nz, ny, nx);
    Ktensor *pem1 = new Ktensor[nx * ny * nz]();
    double *p = new double[nx * ny * nz]();
    int *Ptr = new int[nx * ny * nz + 1]();
    int *Idx = new int[27 * nx * ny * nz]();
    double *Val = new double[27 * nx * ny * nz]();
    double *B = new double[nx * ny * nz]();
    double *x = new double[nx * ny * nz]();
    _get_divide(verts, _verts1, perm, pem1, divn);

    double *cx = new double[nx * ny * nz];
    double *cy = new double[nx * ny * nz];
    double *cz = new double[nx * ny * nz];

    std::mdspan _cx(cx, nz, ny, nx);
    std::mdspan _cy(cy, nz, ny, nx);
    std::mdspan _cz(cz, nz, ny, nx);

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

    int *actnum = new int[nx * ny * nz];
    std::fill(actnum, actnum + nx * ny * nz, 1);

    //_get_tran(nx, ny, nz, actnum, pem1, _verts1, cx, cy, cz, tranx, trany, tranz);

    double epsilonP = 1.0e-6, residual = 0.0001;

    for (iter = 0; iter < 100; iter++)
    {
        n = 0;
        nnz = 0;
        Ptr[0] = 0;

        /*for (int i = 0; i < 27; i++) {
            std::cout << beta[i][0] << " " << beta[i][1] << " " << beta[i][2] << "\n";
        }*/

        std::cout << "iter=" << iter << "\n";
        for (k = 0; k < nz; k++)
        {
            for (j = 0; j < ny; j++)
            {
                for (i = 0; i < nx; i++)
                {
                    // i = 1; j = 1; k = 1;
                    // std::cout << "(" << i << "," << j << "," << k << ") : ";
                    count = 0;
                    double val = 0.0;
                    double F = 0.0;
                    double F1 = 0.0;
                    int i0 = k * nx * ny + j * nx + i;

                    _get_F(pem1, p, i, j, k, _verts1, F, divn);
                    // std::cout << "F=" << F << std::endl;

                    B[n] = -F;
                    if (k > 0)
                    {
                        _get_csr(pem1, p, i, j, k, _verts1, F, i, j, k - 1, Idx, Val, nnz, count, divn);

                        if (i > 0)
                        {
                            _get_csr(pem1, p, i, j, k, _verts1, F, i - 1, j, k - 1, Idx, Val, nnz, count, divn);

                            if (j > 0)
                            {
                                _get_csr(pem1, p, i, j, k, _verts1, F, i - 1, j - 1, k - 1, Idx, Val, nnz, count, divn);
                            }
                            if (j < ny - 1)
                            {
                                _get_csr(pem1, p, i, j, k, _verts1, F, i - 1, j + 1, k - 1, Idx, Val, nnz, count, divn);
                            }
                        }
                        if (i < nx - 1)
                        {
                            _get_csr(pem1, p, i, j, k, _verts1, F, i + 1, j, k - 1, Idx, Val, nnz, count, divn);

                            if (j > 0)
                            {
                                _get_csr(pem1, p, i, j, k, _verts1, F, i + 1, j - 1, k - 1, Idx, Val, nnz, count, divn);
                            }
                            if (j < ny - 1)
                            {
                                _get_csr(pem1, p, i, j, k, _verts1, F, i + 1, j + 1, k - 1, Idx, Val, nnz, count, divn);
                            }
                        }
                        if (j > 0)
                        {
                            _get_csr(pem1, p, i, j, k, _verts1, F, i, j - 1, k - 1, Idx, Val, nnz, count, divn);
                        }
                        if (j < ny - 1)
                        {
                            _get_csr(pem1, p, i, j, k, _verts1, F, i, j + 1, k - 1, Idx, Val, nnz, count, divn);
                        }
                    }

                    if (j > 0)
                    {
                        _get_csr(pem1, p, i, j, k, _verts1, F, i, j - 1, k, Idx, Val, nnz, count, divn);
                    }
                    if (i > 0)
                    {
                        _get_csr(pem1, p, i, j, k, _verts1, F, i - 1, j, k, Idx, Val, nnz, count, divn);

                        if (j > 0)
                        {
                            _get_csr(pem1, p, i, j, k, _verts1, F, i - 1, j - 1, k, Idx, Val, nnz, count, divn);
                        }
                        if (j < ny - 1)
                        {
                            _get_csr(pem1, p, i, j, k, _verts1, F, i - 1, j + 1, k, Idx, Val, nnz, count, divn);
                        }
                    }
                    _get_csr(pem1, p, i, j, k, _verts1, F, i, j, k, Idx, Val, nnz, count, divn);
                    if (i < nx - 1)
                    {
                        _get_csr(pem1, p, i, j, k, _verts1, F, i + 1, j, k, Idx, Val, nnz, count, divn);

                        if (j > 0)
                        {
                            _get_csr(pem1, p, i, j, k, _verts1, F, i + 1, j - 1, k, Idx, Val, nnz, count, divn);
                        }
                        if (j < ny - 1)
                        {
                            _get_csr(pem1, p, i, j, k, _verts1, F, i + 1, j + 1, k, Idx, Val, nnz, count, divn);
                        }
                    }
                    if (j < ny - 1)
                    {
                        _get_csr(pem1, p, i, j, k, _verts1, F, i, j + 1, k, Idx, Val, nnz, count, divn);
                    }

                    if (k < nz - 1)
                    {
                        _get_csr(pem1, p, i, j, k, _verts1, F, i, j, k + 1, Idx, Val, nnz, count, divn);

                        if (i > 0)
                        {
                            _get_csr(pem1, p, i, j, k, _verts1, F, i - 1, j, k + 1, Idx, Val, nnz, count, divn);

                            if (j > 0)
                            {
                                _get_csr(pem1, p, i, j, k, _verts1, F, i - 1, j - 1, k + 1, Idx, Val, nnz, count, divn);
                            }
                            if (j < ny - 1)
                            {
                                _get_csr(pem1, p, i, j, k, _verts1, F, i - 1, j + 1, k + 1, Idx, Val, nnz, count, divn);
                            }
                        }
                        if (i < nx - 1)
                        {
                            _get_csr(pem1, p, i, j, k, _verts1, F, i + 1, j, k + 1, Idx, Val, nnz, count, divn);

                            if (j > 0)
                            {
                                _get_csr(pem1, p, i, j, k, _verts1, F, i + 1, j - 1, k + 1, Idx, Val, nnz, count, divn);
                            }
                            if (j < ny - 1)
                            {
                                _get_csr(pem1, p, i, j, k, _verts1, F, i + 1, j + 1, k + 1, Idx, Val, nnz, count, divn);
                            }
                        }
                        if (j > 0)
                        {
                            _get_csr(pem1, p, i, j, k, _verts1, F, i, j - 1, k + 1, Idx, Val, nnz, count, divn);
                        }
                        if (j < ny - 1)
                        {
                            _get_csr(pem1, p, i, j, k, _verts1, F, i, j + 1, k + 1, Idx, Val, nnz, count, divn);
                        }
                    }
                    Ptr[n + 1] = Ptr[n] + count;

                    n = n + 1;
                    // std::cout << std::endl;
                }
            }
        }

        // FILE* fp = std::fopen("Val_MPFA.INC", "w");
        // for (int ii = 0; ii < nnz; ii++) {
        //     std::fprintf(fp, "%le\n", Val[ii]);
        // }
        // std::fclose(fp);

        pmgmres_ilu_cr(nx * ny * nz, nnz, Ptr, Idx, Val, x, B, 1000, 1000, 1e-5, 1e-5);

        double res = 0.0;

        for (int i = 0; i < nx * ny * nz; ++i)
        {
            if (std::abs(x[i]) > res)
                res = std::abs(x[i]);
        }

        for (int i = 0; i < nx * ny * nz; ++i)
        {
            p[i] += x[i];
        }

        Q = 0.0;
        for (int k = 0; k < nz; k++)
        {
            for (int j = 0; j < ny; j++)
            {
                const double *pp[4];
                double S = 0.0, D = 0.0;
                double pc[3];
                int i0 = k * nx * ny + j * nx + nx - 1;

                pp[0] = &verts1(k, j, nx - 1, 1, 0);
                pp[1] = &verts1(k, j, nx - 1, 3, 0);
                pp[2] = &verts1(k, j, nx - 1, 5, 0);
                pp[3] = &verts1(k, j, nx - 1, 7, 0);

                pc[0] = _cx(k, j, nx - 1);
                pc[1] = _cy(k, j, nx - 1);
                pc[2] = _cz(k, j, nx - 1);

                S = _get_area(pp[0], pp[1], pp[2]) + _get_area(pp[3], pp[1], pp[2]);
                D = _get_distance(pc, pp[1], pp[2], pp[3]);
                Q = Q + (Phigh - p[i0]) * S / D * pem1[i0].x;
            }
        }

        Iter = iter;
        std::cout << "Q=" << Q << "\n";
        std::cout << "res=" << res << "\n";
        if (res < 1.0e-3)
            break;
    }
    delete[] pem1, p, Ptr, Idx, Val, B, x, _verts1, cx, cy, cz;
    std::cout << std::endl;
}