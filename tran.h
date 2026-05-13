#pragma once
#define _CRT_SECURE_NO_WARNINGS
bool _MPFA_ONLY = false;
#include "Checked.h"
#include "roots.h"
#include "utility.h"
#include "GMRES.h"
#include <Eigen/Dense>
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
            ratio[0] = _get_ratio(n[6], t[0]);
            ratio[1] = _get_ratio(n[9], t[1]);
            ratio[2] = _get_ratio(n[7], t[2]);
            ratio[3] = _get_ratio(n[10], t[3]);
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
            _get_triangle_normal(&verts1(k, j, i, 0, 0), pp2, p0, t[2], Axis::XPOSITIVE);
            ratio[2] = _get_ratio(n0, t[2]);
            _cross_product(n4, n8, n48);
            Eigen::Matrix3d matK;
            _get_matK(pem1[cur], matK);
            Eigen::RowVector3d Dx;
            Eigen::Vector3d Dn;
            Dx << p0[0] - _cx(k, j, i), p0[1] - _cy(k, j, i), p0[2] - _cz(k, j, i);
            Dn << n48[0], n48[1], n48[2];
            double cA = Dx * matK.inverse() * Dn;
            cA = _dot_product(n0, n48) / cA;
            cA *= ratio[2];
            B[cur] += pb * cA;
            C(cur, cur) += cA;
        }
        else if (i == 0 && j == ny && k == 0)
        {
            break;
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
            _get_triangle_normal(&verts1(k - 1, j, i, 4, 0), p3, pp2, t[0], Axis::XPOSITIVE);
            ratio[0] = _get_ratio(n3, t[0]);
            _cross_product(n5, n8, n58);
            Eigen::Matrix3d matK;
            _get_matK(pem1[cur], matK);
            Eigen::RowVector3d Dx;
            Eigen::Vector3d Dn;
            Dx << p3[0] - _cx(k - 1, j, i), p3[1] - _cy(k - 1, j, i), p3[2] - _cz(k - 1, j, i);
            Dn << n58[0], n58[1], n58[2];
            double cA = Dx * matK.inverse() * Dn;
            cA = _dot_product(n3, n58) / cA;
            cA *= ratio[0];
            B[cur] += pb * cA;
            C(cur, cur) += cA;
        }
        else if (i == 0 && j == ny && k == nz)
        {
            break;
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
            _get_triangle_normal(&verts1(k, j, i - 1, 1, 0), pp2, p0, t[2], Axis::XPOSITIVE);
            ratio[2] = _get_ratio(n0, t[2]);
            _cross_product(n7, n9, n79);
            Eigen::Matrix3d matK;
            _get_matK(pem1[cur], matK);
            Eigen::RowVector3d Dx;
            Eigen::Vector3d Dn;
            Dx << p0[0] - _cx(k, j, i - 1), p0[1] - _cy(k, j, i - 1), p0[2] - _cz(k, j, i - 1);
            Dn << n79[0], n79[1], n79[2];
            double cA = Dx * matK.inverse() * Dn;
            cA = _dot_product(n0, n79) / cA;
            cA *= ratio[2];
            B[cur] -= pb * cA;
            C(cur, cur) -= cA;
        }
        else if (i == nx && j == ny && k == 0)
        {
            break;
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
            _get_triangle_normal(&verts1(k - 1, j, i - 1, 5, 0), p3, pp2, t[0], Axis::XPOSITIVE);
            ratio[0] = _get_ratio(n3, t[0]);
            _cross_product(n6, n9, n69);
            Eigen::Matrix3d matK;
            _get_matK(pem1[cur], matK);
            Eigen::RowVector3d Dx;
            Eigen::Vector3d Dn;
            Dx << p3[0] - _cx(k - 1, j, i - 1), p3[1] - _cy(k - 1, j, i - 1), p3[2] - _cz(k - 1, j, i - 1);
            Dn << n69[0], n69[1], n69[2];
            double cA = Dx * matK.inverse() * Dn;
            cA = _dot_product(n3, n69) / cA;
            cA *= ratio[0];
            B[cur] -= pb * cA;
            C(cur, cur) -= cA;
        }
        else if (i == nx && j == ny && k == nz)
        {
            break;
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
            _get_triangle_normal(&verts1(k, j, i, 0, 0), pp2, p0, t[2], Axis::XPOSITIVE);
            ratio[2] = _get_ratio(n0, t[2]);
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
            cA *= ratio[2];
            C(i6, i6) += cA;
            C(i7, i7) += cA;
            C(i6, i7) -= cA;
            C(i7, i6) -= cA;
        }
        else if (j == ny && k == 0)
        {
            break;
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
            _get_triangle_normal(&verts1(k - 1, j, i, 4, 0), p3, pp2, t[0], Axis::XPOSITIVE);
            ratio[0] = _get_ratio(n3, t[0]);
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
            cA *= ratio[0];
            C(i2, i2) += cA;
            C(i3, i3) += cA;
            C(i2, i3) -= cA;
            C(i3, i2) -= cA;
        }
        else if (j == ny && k == nz)
        {
            break;
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
            _get_triangle_normal(&verts1(k, j, i, 0, 0), pp2, p0, t[2], Axis::XPOSITIVE);
            ratio[2] = _get_ratio(n0, t[2]);
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
            Eigen::RowVectorXd r = -(ratio[2] * v0) * matK6 * A.bottomRows(3);
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
            _get_triangle_normal(&verts1(k, j, i - 1, 1, 0), pp2, p0, t[2], Axis::XPOSITIVE);
            ratio[2] = _get_ratio(n0, t[2]);
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
            Eigen::RowVectorXd r = (ratio[2] * v0) * matK7 * A.bottomRows(3);
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
            _get_triangle_normal(&verts1(k - 1, j, i, 4, 0), pp2, p3, t[0], Axis::XPOSITIVE);
            ratio[0] = _get_ratio(n3, t[0]);
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
            Eigen::RowVectorXd r = -(ratio[0] * v3) * matK2 * A.bottomRows(3);
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
            _get_triangle_normal(&verts1(k - 1, j, i - 1, 5, 0), pp2, p3, t[0], Axis::XPOSITIVE);
            ratio[0] = _get_ratio(n3, t[0]);
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
            Eigen::RowVectorXd r = (ratio[0] * v3) * matK3 * A.bottomRows(3);
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
            _get_triangle_normal(&verts1(k, j, i, 0, 0), pp2, p0, t[2], Axis::XPOSITIVE);
            _get_triangle_normal(&verts1(k - 1, j, i, 4, 0), pp2, p3, t[0], Axis::XPOSITIVE);
            _get_triangle_normal(&verts1(k - 1, j, i, 4, 0), p8, pp2, t[1], Axis::ZPOSITIVE);
            ratio[0] = _get_ratio(n3, t[0]);
            ratio[1] = _get_ratio(n8, t[1]);
            ratio[2] = _get_ratio(n0, t[2]);
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
            Eigen::RowVectorXd r = (ratio[1] * v8 - ratio[0] * v3) * matK2 * A.topRows(3);
            B[i2] -= pb * (r[0] + r[1]);
            C(i2, i2) += -r[1] + r[2];
            C(i2, i6) += -r[0] - r[2];
            r = -(ratio[2] * v0 + ratio[1] * v8) * matK6 * A.bottomRows(3);
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
            _get_triangle_normal(&verts1(k, j, i - 1, 1, 0), pp2, p0, t[2], Axis::XPOSITIVE);
            _get_triangle_normal(&verts1(k - 1, j, i - 1, 5, 0), p3, pp2, t[0], Axis::XPOSITIVE);
            _get_triangle_normal(&verts1(k - 1, j, i - 1, 5, 0), pp2, p9, t[3], Axis::ZPOSITIVE);
            ratio[0] = _get_ratio(n3, t[0]);
            ratio[2] = _get_ratio(n0, t[2]);
            ratio[3] = _get_ratio(n9, t[3]);
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
            Eigen::RowVectorXd r = (ratio[3] * v9 + ratio[0] * v3) * matK3 * A.topRows(3);
            B[i3] -= pb * (r[0] + r[1]);
            C(i3, i3) += -r[1] + r[2];
            C(i3, i7) += -r[0] - r[2];
            r = (ratio[2] * v0 - ratio[3] * v9) * matK7 * A.bottomRows(3);
            B[i7] -= pb * (r[0] + r[1]);
            C(i7, i3) += -r[1] + r[2];
            C(i7, i7) += -r[0] - r[2];
            //  To be implemented
        }
        else if (i == 0 && j == ny)
        {
            break;
        }
        else if (i == nx && j == ny)
        {
            break;
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
            _get_triangle_normal(&verts1(k, j, i, 0, 0), pp2, p0, t[2], Axis::XPOSITIVE);
            _get_triangle_normal(&verts1(k, j, i, 0, 0), p3, pp2, t[0], Axis::XPOSITIVE);
            _get_triangle_normal(&verts1(k, j, i, 0, 0), pp2, p8, t[1], Axis::ZPOSITIVE);
            ratio[0] = _get_ratio(n3, t[0]);
            ratio[1] = _get_ratio(n8, t[1]);
            ratio[2] = _get_ratio(n0, t[2]);
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
            Eigen::RowVectorXd r = (-ratio[0] * v3 + ratio[1] * v8) * matK2 * A.block(3, 0, 3, 12);
            B[i2] -= pb * (r[0] + r[1] + r[2] + r[3]);
            C(i2, i1) += -r[2] + r[6] + r[10];
            C(i2, i2) += -r[3] - r[6] + r[8];
            C(i2, i5) += -r[1] + r[4] - r[10];
            C(i2, i6) += -r[0] - r[4] - r[8];
            r = -(ratio[2] * v0 + ratio[1] * v8) * matK6 * A.bottomRows(3);
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
            _get_triangle_normal(&verts1(k, j, i - 1, 1, 0), pp2, p0, t[2], Axis::XPOSITIVE);
            _get_triangle_normal(&verts1(k, j, i - 1, 1, 0), p3, pp2, t[0], Axis::XPOSITIVE);
            _get_triangle_normal(&verts1(k, j, i - 1, 1, 0), pp2, p9, t[3], Axis::ZPOSITIVE);
            ratio[0] = _get_ratio(n3, t[0]);
            ratio[2] = _get_ratio(n0, t[2]);
            ratio[3] = _get_ratio(n9, t[3]);
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
            Eigen::RowVectorXd r = (ratio[0] * v3 + ratio[3] * v9) * matK3 * A.block(3, 0, 3, 12);
            B[i3] -= pb * (r[0] + r[1] + r[2] + r[3]);
            C(i3, i0) += -r[2] + r[4] + r[10];
            C(i3, i3) += -r[3] - r[4] + r[8];
            C(i3, i4) += -r[1] + r[6] - r[10];
            C(i3, i7) += -r[0] - r[6] - r[8];
            r = (ratio[2] * v0 - ratio[3] * v9) * matK7 * A.bottomRows(3);
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
            _get_triangle_normal(&verts1(k, j, i, 0, 0), pp2, p0, t[2], Axis::XPOSITIVE);
            _get_triangle_normal(&verts1(k, j, i, 0, 0), pp2, p3, t[0], Axis::XPOSITIVE);
            _get_triangle_normal(&verts1(k, j, i, 0, 0), pp2, p8, t[1], Axis::ZPOSITIVE);
            _get_triangle_normal(&verts1(k, j, i, 0, 0), pp2, p9, t[3], Axis::ZPOSITIVE);
            ratio[0] = _get_ratio(n3, t[0]);
            ratio[1] = _get_ratio(n8, t[1]);
            ratio[2] = _get_ratio(n0, t[2]);
            ratio[3] = _get_ratio(n9, t[3]);
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
            Eigen::RowVectorXd r = (-ratio[0] * v3 + ratio[1] * v8) * matK2 * A.topRows(3);
            C(i2, i2) += r[2] + r[8];
            C(i2, i3) += -r[2] + r[10];
            C(i2, i6) += r[0] - r[8];
            C(i2, i7) += -r[0] - r[10];
            r = (ratio[0] * v3 + ratio[3] * v9) * matK3 * A.block(3, 0, 3, 12);
            C(i3, i2) += r[2] + r[8];
            C(i3, i3) += -r[2] + r[10];
            C(i3, i6) += r[0] - r[8];
            C(i3, i7) += -r[0] - r[10];
            r = (-ratio[2] * v0 - ratio[1] * v8) * matK6 * A.block(6, 0, 3, 12);
            C(i6, i2) += r[2] + r[8];
            C(i6, i3) += -r[2] + r[10];
            C(i6, i6) += r[0] - r[8];
            C(i6, i7) += -r[0] - r[10];
            r = (ratio[2] * v0 - ratio[3] * v9) * matK7 * A.bottomRows(3);
            C(i7, i2) += r[2] + r[8];
            C(i7, i3) += -r[2] + r[10];
            C(i7, i6) += r[0] - r[8];
            C(i7, i7) += -r[0] - r[10];
            // To be implemented
        }
        else if (j == ny)
        {
            break;
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
            _get_triangle_normal(&verts1(k, j, i, 0, 0), pp2, p0, t[2], Axis::XPOSITIVE);
            ratio[2] = _get_ratio(n0, t[2]);
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
            Eigen::RowVectorXd r = (-ratio[2] * v0) * matK6 * A.block(6, 0, 3, 12);
            C(i6, i4) += -r[2] + r[6];
            C(i6, i5) += r[2] + r[4];
            C(i6, i6) += r[0] - r[4];
            C(i6, i7) += -r[0] - r[6];
            r = (ratio[2] * v0) * matK7 * A.bottomRows(3);
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
            _get_triangle_normal(&verts1(k - 1, j, i, 4, 0), p3, pp2, t[0], Axis::XPOSITIVE);
            ratio[0] = _get_ratio(n3, t[0]);
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
            Eigen::RowVectorXd r = (-ratio[0] * v3) * matK2 * A.block(6, 0, 3, 12);
            C(i2, i0) += r[0] + r[6];
            C(i2, i1) += -r[0] + r[4];
            C(i2, i2) += -r[2] - r[4];
            C(i2, i3) += r[2] - r[6];
            r = (ratio[0] * v3) * matK3 * A.bottomRows(3);
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
            _get_triangle_normal(&verts1(k, j, i, 0, 0), pp[2], p[0], t[2], Axis::XPOSITIVE);
            _get_triangle_normal(&verts1(k, j, i, 0, 0), pp[2], p[3], t[0], Axis::XPOSITIVE);
            _get_triangle_normal(&verts1(k, j, i, 0, 0), pp[2], p[8], t[1], Axis::ZPOSITIVE);
            _get_triangle_normal(&verts1(k, j, i, 0, 0), pp[2], p[9], t[3], Axis::ZPOSITIVE);
            ratio[0] = _get_ratio(n[3], t[0]);
            ratio[1] = _get_ratio(n[8], t[1]);
            ratio[2] = _get_ratio(n[0], t[2]);
            ratio[3] = _get_ratio(n[9], t[3]);
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
            Eigen::RowVectorXd r = (-ratio[0] * v[3] + ratio[1] * v[8]) * matK[2] * A.block(6, 0, 3, 24);
            C(idx[2], idx[0]) += r[4] + r[12] + r[20];
            C(idx[2], idx[1]) += -r[4] + r[10] + r[22];
            C(idx[2], idx[2]) += r[6] - r[10] + r[16];
            C(idx[2], idx[3]) += -r[6] - r[12] + r[18];
            C(idx[2], idx[4]) += r[2] + r[14] - r[20];
            C(idx[2], idx[5]) += -r[2] + r[8] - r[22];
            C(idx[2], idx[6]) += r[0] - r[8] - r[16];
            C(idx[2], idx[7]) += -r[0] - r[14] - r[18];
            r = (ratio[0] * v[3] + ratio[3] * v[9]) * matK[3] * A.block(9, 0, 3, 24);
            C(idx[3], idx[0]) += r[4] + r[12] + r[20];
            C(idx[3], idx[1]) += -r[4] + r[10] + r[22];
            C(idx[3], idx[2]) += r[6] - r[10] + r[16];
            C(idx[3], idx[3]) += -r[6] - r[12] + r[18];
            C(idx[3], idx[4]) += r[2] + r[14] - r[20];
            C(idx[3], idx[5]) += -r[2] + r[8] - r[22];
            C(idx[3], idx[6]) += r[0] - r[8] - r[16];
            C(idx[3], idx[7]) += -r[0] - r[14] - r[18];
            r = (-ratio[2] * v[0] - ratio[1] * v[8]) * matK[6] * A.block(18, 0, 3, 24);
            C(idx[6], idx[0]) += r[4] + r[12] + r[20];
            C(idx[6], idx[1]) += -r[4] + r[10] + r[22];
            C(idx[6], idx[2]) += r[6] - r[10] + r[16];
            C(idx[6], idx[3]) += -r[6] - r[12] + r[18];
            C(idx[6], idx[4]) += r[2] + r[14] - r[20];
            C(idx[6], idx[5]) += -r[2] + r[8] - r[22];
            C(idx[6], idx[6]) += r[0] - r[8] - r[16];
            C(idx[6], idx[7]) += -r[0] - r[14] - r[18];
            r = (ratio[2] * v[0] - ratio[3] * v[9]) * matK[7] * A.block(21, 0, 3, 24);
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
    case Axis::YNEGATIVE:
        assert(j > 0);
        // 3个面的点，2个未知数,8个点
        if (i == 0 && j == 0 && k == 0)
        {
            break;
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
            _get_triangle_normal(&verts1(k, j - 1, i, 2, 0), pp0, p1, t[2], Axis::XPOSITIVE);
            ratio[2] = _get_ratio(n1, t[2]);
            _cross_product(n4, n11, n411);
            Eigen::Matrix3d matK;
            _get_matK(pem1[cur], matK);
            Eigen::RowVector3d Dx;
            Eigen::Vector3d Dn;
            Dx << p1[0] - _cx(k, j - 1, i), p1[1] - _cy(k, j - 1, i), p1[2] - _cz(k, j - 1, i);
            Dn << n411[0], n411[1], n411[2];
            double cA = Dx * matK.inverse() * Dn;
            cA = _dot_product(n1, n411) / cA;
            cA *= ratio[2];
            B[cur] += pb * cA;
            C(cur, cur) += cA;
        }
        else if (i == 0 && j == 0 && k == nz)
        {
            break;
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
            _get_triangle_normal(&verts1(k - 1, j - 1, i, 6, 0), pp0, p2, t[0], Axis::XPOSITIVE);
            ratio[0] = _get_ratio(n2, t[0]);
            _cross_product(n5, n11, n511);
            Eigen::Matrix3d matK;
            _get_matK(pem1[cur], matK);
            Eigen::RowVector3d Dx;
            Eigen::Vector3d Dn;
            Dx << p2[0] - _cx(k - 1, j - 1, i), p2[1] - _cy(k - 1, j - 1, i), p2[2] - _cz(k - 1, j - 1, i);
            Dn << n511[0], n511[1], n511[2];
            double cA = Dx * matK.inverse() * Dn;
            cA = _dot_product(n2, n511) / cA;
            cA *= ratio[0];
            B[cur] += pb * cA;
            C(cur, cur) += cA;
        }
        else if (i == nx && j == 0 && k == 0)
        {
            break;
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
            _get_triangle_normal(&verts1(k, j - 1, i - 1, 3, 0), pp0, p1, t[2], Axis::XPOSITIVE);
            ratio[2] = _get_ratio(n1, t[2]);
            _cross_product(n7, n10, n710);
            Eigen::Matrix3d matK;
            _get_matK(pem1[cur], matK);
            Eigen::RowVector3d Dx;
            Eigen::Vector3d Dn;
            Dx << p1[0] - _cx(k, j - 1, i - 1), p1[1] - _cy(k, j - 1, i - 1), p1[2] - _cz(k, j - 1, i - 1);
            Dn << n710[0], n710[1], n710[2];
            double cA = Dx * matK.inverse() * Dn;
            cA = _dot_product(n1, n710) / cA;
            cA *= ratio[2];
            B[cur] -= pb * cA;
            C(cur, cur) -= cA;
        }
        else if (i == nx && j == 0 && k == nz)
        {
            break;
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
            _get_triangle_normal(&verts1(k - 1, j - 1, i - 1, 7, 0), pp0, p2, t[0], Axis::XPOSITIVE);
            ratio[0] = _get_ratio(n2, t[0]);
            _cross_product(n6, n10, n610);
            Eigen::Matrix3d matK;
            _get_matK(pem1[cur], matK);
            Eigen::RowVector3d Dx;
            Eigen::Vector3d Dn;
            Dx << p2[0] - _cx(k - 1, j - 1, i - 1), p2[1] - _cy(k - 1, j - 1, i - 1), p2[2] - _cz(k - 1, j - 1, i - 1);
            Dn << n610[0], n610[1], n610[2];
            double cA = Dx * matK.inverse() * Dn;
            cA = _dot_product(n2, n610) / cA;
            cA *= ratio[0];
            B[cur] -= pb * cA;
            C(cur, cur) -= cA;
        }
        // 12条边上其余点的处理
        else if (j == 0 && k == 0)
        {
            break;
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
            _get_triangle_normal(&verts1(k, j - 1, i - 1, 3, 0), pp0, p1, t[2], Axis::XPOSITIVE);
            ratio[2] = _get_ratio(n1, t[2]);
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
            cA *= ratio[2];
            C(i5, i5) += cA;
            C(i4, i4) += cA;
            C(i5, i4) -= cA;
            C(i4, i5) -= cA;
        }
        else if (j == 0 && k == nz)
        {
            break;
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
            _get_triangle_normal(&verts1(k - 1, j - 1, i, 6, 0), pp0, p2, t[0], Axis::XPOSITIVE);
            ratio[0] = _get_ratio(n2, t[0]);
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
            cA *= ratio[0];
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
            _get_triangle_normal(&verts1(k, j - 1, i, 2, 0), pp0, p1, t[2], Axis::XPOSITIVE);
            ratio[2] = _get_ratio(n1, t[2]);
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
            Eigen::RowVectorXd r = (-ratio[2] * v1) * matK5 * A.topRows(3);
            B[i5] -= pb * (r[0] + r[1]);
            C(i5, i5) += -r[1] + r[2];
            C(i5, i6) += -r[0] - r[2];
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
            _get_triangle_normal(&verts1(k, j - 1, i - 1, 3, 0), pp0, p1, t[2], Axis::XPOSITIVE);
            ratio[2] = _get_ratio(n1, t[2]);
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
            Eigen::RowVectorXd r = (ratio[2] * v1) * matK4 * A.topRows(3);
            B[i4] -= pb * (r[0] + r[1]);
            C(i4, i4) += -r[1] + r[2];
            C(i4, i7) += -r[0] - r[2];
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
            _get_triangle_normal(&verts1(k - 1, j - 1, i, 6, 0), pp0, p2, t[0], Axis::XPOSITIVE);
            ratio[0] = _get_ratio(n2, t[0]);
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
            Eigen::RowVectorXd r = (-ratio[0] * v2) * matK1 * A.topRows(3);
            B[i1] -= pb * (r[0] + r[1]);
            C(i1, i1) += -r[1] + r[2];
            C(i1, i2) += -r[0] - r[2];
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
            _get_triangle_normal(&verts1(k - 1, j - 1, i - 1, 7, 0), pp0, p2, t[0], Axis::XPOSITIVE);
            ratio[0] = _get_ratio(n2, t[0]);
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
            Eigen::RowVectorXd r = (ratio[0] * v2) * matK0 * A.topRows(3);
            B[i0] -= pb * (r[0] + r[1]);
            C(i0, i0) += -r[1] + r[2];
            C(i0, i3) += -r[0] - r[2];
        }
        else if (i == 0 && j == 0)
        {
            break;
        }
        else if (i == nx && j == 0)
        {
            break;
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
            _get_triangle_normal(&verts1(k, j - 1, i, 2, 0), pp0, p1, t[2], Axis::XPOSITIVE);
            _get_triangle_normal(&verts1(k - 1, j - 1, i, 6, 0), p2, pp0, t[0], Axis::XPOSITIVE);
            _get_triangle_normal(&verts1(k - 1, j - 1, i, 6, 0), pp0, p11, t[1], Axis::ZPOSITIVE);
            ratio[0] = _get_ratio(n2, t[0]);
            ratio[1] = _get_ratio(n11, t[1]);
            ratio[2] = _get_ratio(n1, t[2]);
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
            Eigen::RowVectorXd r = (ratio[1] * v11 - ratio[0] * v2) * matK1 * A.topRows(3);
            B[i1] -= pb * (r[0] + r[1]);
            C(i1, i1) += -r[1] + r[2];
            C(i1, i5) += -r[0] - r[2];
            r = -(ratio[2] * v1 + ratio[1] * v11) * matK5 * A.bottomRows(3);
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
            _get_triangle_normal(&verts1(k, j - 1, i - 1, 3, 0), pp0, p1, t[2], Axis::XPOSITIVE);
            _get_triangle_normal(&verts1(k - 1, j - 1, i - 1, 7, 0), p2, pp0, t[0], Axis::XPOSITIVE);
            _get_triangle_normal(&verts1(k - 1, j - 1, i - 1, 7, 0), pp0, p10, t[3], Axis::ZPOSITIVE);
            ratio[0] = _get_ratio(n2, t[0]);
            ratio[2] = _get_ratio(n1, t[2]);
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
            Eigen::RowVectorXd r = (ratio[3] * v10 + ratio[0] * v2) * matK0 * A.topRows(3);
            B[i0] -= pb * (r[0] + r[1]);
            C(i0, i0) += -r[1] + r[2];
            C(i0, i4) += -r[0] - r[2];
            r = (ratio[2] * v1 - ratio[3] * v10) * matK4 * A.bottomRows(3);
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
            _get_triangle_normal(&verts1(k, j, i, 0, 0), pp0, p1, t[2], Axis::XPOSITIVE);
            _get_triangle_normal(&verts1(k, j, i, 0, 0), pp0, p2, t[0], Axis::XPOSITIVE);
            _get_triangle_normal(&verts1(k, j, i, 0, 0), pp0, p11, t[1], Axis::ZPOSITIVE);
            ratio[0] = _get_ratio(n2, t[0]);
            ratio[1] = _get_ratio(n11, t[1]);
            ratio[2] = _get_ratio(n1, t[2]);
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
            Eigen::RowVectorXd r = (-ratio[0] * v2 + ratio[1] * v11) * matK1 * A.topRows(3);
            B[i1] -= pb * (r[0] + r[1] + r[2] + r[3]);
            C(i1, i1) += -r[2] + r[6] + r[10];
            C(i1, i2) += -r[3] - r[6] + r[8];
            C(i1, i5) += -r[1] + r[4] - r[10];
            C(i1, i6) += -r[0] - r[4] - r[8];
            r = (-ratio[2] * v1 - ratio[1] * v11) * matK5 * A.block(6, 0, 3, 12);
            B[i5] -= pb * (r[0] + r[1] + r[2] + r[3]);
            C(i5, i1) += -r[2] + r[6] + r[10];
            C(i5, i2) += -r[3] - r[6] + r[8];
            C(i5, i5) += -r[1] + r[4] - r[10];
            C(i5, i6) += -r[0] - r[4] - r[8];
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
            _get_triangle_normal(&verts1(k, j, i - 1, 1, 0), pp0, p1, t[2], Axis::XPOSITIVE);
            _get_triangle_normal(&verts1(k, j, i - 1, 1, 0), pp0, p2, t[0], Axis::XPOSITIVE);
            _get_triangle_normal(&verts1(k, j, i - 1, 1, 0), pp0, p10, t[3], Axis::ZPOSITIVE);
            ratio[0] = _get_ratio(n2, t[0]);
            ratio[2] = _get_ratio(n1, t[2]);
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
            Eigen::RowVectorXd r = (ratio[0] * v2 + ratio[3] * v10) * matK0 * A.topRows(3);
            B[i0] -= pb * (r[0] + r[1] + r[2] + r[3]);
            C(i0, i0) += -r[2] + r[4] + r[10];
            C(i0, i3) += -r[3] - r[4] + r[8];
            C(i0, i4) += -r[1] + r[6] - r[10];
            C(i0, i7) += -r[0] - r[6] - r[8];
            r = (ratio[2] * v1 - ratio[3] * v10) * matK4 * A.block(6, 0, 3, 12);
            B[i4] -= pb * (r[0] + r[1] + r[2] + r[3]);
            C(i4, i0) += -r[2] + r[4] + r[10];
            C(i4, i3) += -r[3] - r[4] + r[8];
            C(i4, i4) += -r[1] + r[6] - r[10];
            C(i4, i7) += -r[0] - r[6] - r[8];
            // To be implemented
        }
        else if (j == 0)
        {
            break;
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
            _get_triangle_normal(&verts1(k, j - 1, i, 2, 0), pp0, p1, t[2], Axis::XPOSITIVE);
            _get_triangle_normal(&verts1(k, j - 1, i, 2, 0), pp0, p2, t[0], Axis::XPOSITIVE);
            _get_triangle_normal(&verts1(k, j - 1, i, 2, 0), pp0, p10, t[3], Axis::ZPOSITIVE);
            _get_triangle_normal(&verts1(k, j - 1, i, 2, 0), pp0, p11, t[1], Axis::ZPOSITIVE);
            ratio[0] = _get_ratio(n2, t[0]);
            ratio[1] = _get_ratio(n11, t[1]);
            ratio[2] = _get_ratio(n1, t[2]);
            ratio[3] = _get_ratio(n10, t[3]);
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
            Eigen::RowVectorXd r = (ratio[0] * v2 + ratio[3] * v10) * matK0 * A.topRows(3);
            C(i0, i0) += r[2] + r[8];
            C(i0, i1) += -r[2] + r[10];
            C(i0, i4) += r[0] - r[8];
            C(i0, i5) += -r[0] - r[10];
            r = (-ratio[0] * v2 + ratio[1] * v11) * matK1 * A.block(3, 0, 3, 12);
            C(i1, i0) += r[2] + r[8];
            C(i1, i1) += -r[2] + r[10];
            C(i1, i4) += r[0] - r[8];
            C(i1, i5) += -r[0] - r[10];
            r = (ratio[2] * v1 - ratio[3] * v10) * matK4 * A.block(6, 0, 3, 12);
            C(i4, i0) += r[2] + r[8];
            C(i4, i1) += -r[2] + r[10];
            C(i4, i4) += r[0] - r[8];
            C(i4, i5) += -r[0] - r[10];
            r = (-ratio[2] * v1 - ratio[1] * v11) * matK5 * A.bottomRows(3);
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
            _get_triangle_normal(&verts1(k, j, i, 0, 0), p1, pp0, t[2], Axis::XPOSITIVE);
            ratio[2] = _get_ratio(n1, t[2]);
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
            Eigen::RowVectorXd r = (ratio[2] * v1) * matK4 * A.topRows(3);
            C(i4, i4) += -r[2] + r[6];
            C(i4, i5) += r[2] + r[4];
            C(i4, i6) += r[0] - r[4];
            C(i4, i7) += -r[0] - r[6];
            r = (-ratio[2] * v1) * matK5 * A.block(3, 0, 3, 12);
            C(i5, i4) += -r[2] + r[6];
            C(i5, i5) += r[2] + r[4];
            C(i5, i6) += r[0] - r[4];
            C(i5, i7) += -r[0] - r[6];
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
            _get_triangle_normal(&verts1(k - 1, j, i, 4, 0), pp0, p2, t[0], Axis::XPOSITIVE);
            ratio[0] = _get_ratio(n2, t[0]);
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
            Eigen::RowVectorXd r = (ratio[0] * v2) * matK0 * A.topRows(3);
            C(i0, i0) += r[0] + r[6];
            C(i0, i1) += -r[0] + r[4];
            C(i0, i2) += -r[2] - r[4];
            C(i0, i3) += r[2] - r[6];
            r = (-ratio[0] * v2) * matK1 * A.block(3, 0, 3, 12);
            C(i1, i0) += r[0] + r[6];
            C(i1, i1) += -r[0] + r[4];
            C(i1, i2) += -r[2] - r[4];
            C(i1, i3) += r[2] - r[6];
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
            _get_triangle_normal(&verts1(k, j, i, 0, 0), pp[0], p[1], t[2], Axis::XPOSITIVE);
            _get_triangle_normal(&verts1(k, j, i, 0, 0), pp[0], p[2], t[0], Axis::XPOSITIVE);
            _get_triangle_normal(&verts1(k, j, i, 0, 0), pp[0], p[10], t[3], Axis::ZPOSITIVE);
            _get_triangle_normal(&verts1(k, j, i, 0, 0), pp[0], p[11], t[1], Axis::ZPOSITIVE);
            ratio[0] = _get_ratio(n[2], t[0]);
            ratio[1] = _get_ratio(n[11], t[1]);
            ratio[2] = _get_ratio(n[1], t[2]);
            ratio[3] = _get_ratio(n[10], t[3]);
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
            Eigen::RowVectorXd r = (ratio[0] * v[2] + ratio[3] * v[10]) * matK[0] * A.topRows(3);
            C(idx[0], idx[0]) += r[4] + r[12] + r[20];
            C(idx[0], idx[1]) += -r[4] + r[10] + r[22];
            C(idx[0], idx[2]) += r[6] - r[10] + r[16];
            C(idx[0], idx[3]) += -r[6] - r[12] + r[18];
            C(idx[0], idx[4]) += r[2] + r[14] - r[20];
            C(idx[0], idx[5]) += -r[2] + r[8] - r[22];
            C(idx[0], idx[6]) += r[0] - r[8] - r[16];
            C(idx[0], idx[7]) += -r[0] - r[14] - r[18];
            r = (-ratio[0] * v[2] + ratio[1] * v[11]) * matK[1] * A.block(3, 0, 3, 24);
            C(idx[1], idx[0]) += r[4] + r[12] + r[20];
            C(idx[1], idx[1]) += -r[4] + r[10] + r[22];
            C(idx[1], idx[2]) += r[6] - r[10] + r[16];
            C(idx[1], idx[3]) += -r[6] - r[12] + r[18];
            C(idx[1], idx[4]) += r[2] + r[14] - r[20];
            C(idx[1], idx[5]) += -r[2] + r[8] - r[22];
            C(idx[1], idx[6]) += r[0] - r[8] - r[16];
            C(idx[1], idx[7]) += -r[0] - r[14] - r[18];
            r = (ratio[2] * v[1] - ratio[3] * v[10]) * matK[4] * A.block(12, 0, 3, 24);
            C(idx[4], idx[0]) += r[4] + r[12] + r[20];
            C(idx[4], idx[1]) += -r[4] + r[10] + r[22];
            C(idx[4], idx[2]) += r[6] - r[10] + r[16];
            C(idx[4], idx[3]) += -r[6] - r[12] + r[18];
            C(idx[4], idx[4]) += r[2] + r[14] - r[20];
            C(idx[4], idx[5]) += -r[2] + r[8] - r[22];
            C(idx[4], idx[6]) += r[0] - r[8] - r[16];
            C(idx[4], idx[7]) += -r[0] - r[14] - r[18];
            r = (-ratio[2] * v[1] - ratio[1] * v[11]) * matK[5] * A.block(15, 0, 3, 24);
            C(idx[5], idx[0]) += r[4] + r[12] + r[20];
            C(idx[5], idx[1]) += -r[4] + r[10] + r[22];
            C(idx[5], idx[2]) += r[6] - r[10] + r[16];
            C(idx[5], idx[3]) += -r[6] - r[12] + r[18];
            C(idx[5], idx[4]) += r[2] + r[14] - r[20];
            C(idx[5], idx[5]) += -r[2] + r[8] - r[22];
            C(idx[5], idx[6]) += r[0] - r[8] - r[16];
            C(idx[5], idx[7]) += -r[0] - r[14] - r[18];
            // 内部角点
        }
        break;
    case Axis::ZPOSITIVE:
        assert(k < nz);
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
            _get_triangle_normal(&verts1(k, j, i, 0, 0), p0, pp5, t[2], Axis::XPOSITIVE);
            ratio[2] = _get_ratio(n0, t[2]);
            _cross_product(n4, n8, n48);
            Eigen::Matrix3d matK;
            _get_matK(pem1[cur], matK);
            Eigen::RowVector3d Dx;
            Eigen::Vector3d Dn;
            Dx << p0[0] - _cx(k, j, i), p0[1] - _cy(k, j, i), p0[2] - _cz(k, j, i);
            Dn << n48[0], n48[1], n48[2];
            double cA = Dx * matK.inverse() * Dn;
            cA = _dot_product(n0, n48) / cA;
            cA *= ratio[2];
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
            _get_triangle_normal(&verts1(k, j - 1, i, 2, 0), p1, pp5, t[0], Axis::XPOSITIVE);
            ratio[0] = _get_ratio(n1, t[0]);
            _cross_product(n4, n11, n411);
            Eigen::Matrix3d matK;
            _get_matK(pem1[cur], matK);
            Eigen::RowVector3d Dx;
            Eigen::Vector3d Dn;
            Dx << p1[0] - _cx(k, j - 1, i), p1[1] - _cy(k, j - 1, i), p1[2] - _cz(k, j - 1, i);
            Dn << n411[0], n411[1], n411[2];
            double cA = Dx * matK.inverse() * Dn;
            cA = _dot_product(n1, n411) / cA;
            cA *= ratio[0];
            B[cur] += pb * cA;
            C(cur, cur) += cA;
        }
        else if (i == 0 && j == 0 && k == nz)
        {
            break;
        }
        else if (i == 0 && j == ny && k == nz)
        {
            break;
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
            _get_triangle_normal(&verts1(k, j, i - 1, 1, 0), p0, pp5, t[2], Axis::XPOSITIVE);
            ratio[2] = _get_ratio(n0, t[2]);
            _cross_product(n7, n9, n79);
            Eigen::Matrix3d matK;
            _get_matK(pem1[cur], matK);
            Eigen::RowVector3d Dx;
            Eigen::Vector3d Dn;
            Dx << p0[0] - _cx(k, j, i - 1), p0[1] - _cy(k, j, i - 1), p0[2] - _cz(k, j, i - 1);
            Dn << n79[0], n79[1], n79[2];
            double cA = Dx * matK.inverse() * Dn;
            cA = _dot_product(n0, n79) / cA;
            cA *= ratio[2];
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
            _get_triangle_normal(&verts1(k, j - 1, i - 1, 3, 0), p1, pp5, t[0], Axis::XPOSITIVE);
            ratio[0] = _get_ratio(n1, t[0]);
            _cross_product(n7, n10, n710);
            Eigen::Matrix3d matK;
            _get_matK(pem1[cur], matK);
            Eigen::RowVector3d Dx;
            Eigen::Vector3d Dn;
            Dx << p1[0] - _cx(k, j - 1, i - 1), p1[1] - _cy(k, j - 1, i - 1), p1[2] - _cz(k, j - 1, i - 1);
            Dn << n710[0], n710[1], n710[2];
            double cA = Dx * matK.inverse() * Dn;
            cA = _dot_product(n1, n710) / cA;
            cA *= ratio[0];
            B[cur] -= pb * cA;
            C(cur, cur) -= cA;
        }
        else if (i == nx && j == 0 && k == nz)
        {
            break;
        }
        else if (i == nx && j == ny && k == nz)
        {
            break;
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
            _get_triangle_normal(&verts1(k, j, i, 0, 0), p0, pp5, t[2], Axis::XPOSITIVE);
            ratio[2] = _get_ratio(n0, t[2]);
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
            cA *= ratio[2];
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
            _get_triangle_normal(&verts1(k, j - 1, i - 1, 3, 0), p1, pp5, t[0], Axis::XPOSITIVE);
            ratio[0] = _get_ratio(n1, t[0]);
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
            cA *= ratio[0];
            C(i5, i5) += cA;
            C(i4, i4) += cA;
            C(i5, i4) -= cA;
            C(i4, i5) -= cA;
        }
        else if (j == 0 && k == nz)
        {
            break;
        }
        else if (j == ny && k == nz)
        {
            break;
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
            _get_triangle_normal(&verts1(k, j, i, 0, 0), p0, pp5, t[2], Axis::XPOSITIVE);
            _get_triangle_normal(&verts1(k, j - 1, i, 2, 0), p1, pp5, t[0], Axis::XPOSITIVE);
            _get_triangle_normal(&verts1(k, j, i, 0, 0), pp5, p4, t[1], Axis::YPOSITIVE);
            ratio[0] = _get_ratio(n1, t[0]);
            ratio[1] = _get_ratio(n4, t[1]);
            ratio[2] = _get_ratio(n0, t[2]);
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
            Eigen::RowVectorXd r = (ratio[1] * v4 - ratio[0] * v1) * matK5 * A.topRows(3);
            B[i5] -= pb * (r[0] + r[1]);
            C(i5, i5) += -r[1] + r[2];
            C(i5, i6) += -r[0] - r[2];
            r = -(ratio[2] * v0 + ratio[1] * v4) * matK6 * A.bottomRows(3);
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
            _get_triangle_normal(&verts1(k, j, i - 1, 1, 0), p0, pp5, t[2], Axis::XPOSITIVE);
            _get_triangle_normal(&verts1(k, j - 1, i - 1, 3, 0), p1, pp5, t[0], Axis::XPOSITIVE);
            _get_triangle_normal(&verts1(k, j, i - 1, 1, 0), pp5, p7, t[3], Axis::YPOSITIVE);
            ratio[0] = _get_ratio(n1, t[0]);
            ratio[2] = _get_ratio(n0, t[2]);
            ratio[3] = _get_ratio(n7, t[3]);
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
            Eigen::RowVectorXd r = (ratio[3] * v7 + ratio[0] * v1) * matK4 * A.topRows(3);
            B[i4] -= pb * (r[0] + r[1]);
            C(i4, i4) += -r[1] + r[2];
            C(i4, i7) += -r[0] - r[2];
            r = (ratio[2] * v0 - ratio[3] * v7) * matK7 * A.bottomRows(3);
            B[i7] -= pb * (r[0] + r[1]);
            C(i7, i4) += -r[1] + r[2];
            C(i7, i7) += -r[0] - r[2];
        }
        else if (i == 0 && k == nz)
        {
            break;
        }
        else if (i == nx && k == nz)
        {
            break;
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
            _get_triangle_normal(&verts1(k, j, i, 0, 0), p0, pp5, t[2], Axis::XPOSITIVE);
            ratio[2] = _get_ratio(n0, t[2]);
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
            Eigen::RowVectorXd r = -(ratio[2] * v0) * matK6 * A.bottomRows(3);
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
            _get_triangle_normal(&verts1(k, j, i - 1, 1, 0), p0, pp5, t[2], Axis::XPOSITIVE);
            ratio[2] = _get_ratio(n0, t[2]);
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
            Eigen::RowVectorXd r = (ratio[2] * v0) * matK7 * A.bottomRows(3);
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
            _get_triangle_normal(&verts1(k, j - 1, i, 2, 0), p1, pp5, t[0], Axis::XPOSITIVE);
            ratio[0] = _get_ratio(n1, t[0]);
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
            Eigen::RowVectorXd r = -(ratio[0] * v1) * matK5 * A.bottomRows(3);
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
            _get_triangle_normal(&verts1(k, j - 1, i - 1, 3, 0), p1, pp5, t[0], Axis::XPOSITIVE);
            ratio[0] = _get_ratio(n1, t[0]);
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
            Eigen::RowVectorXd r = (ratio[0] * v1) * matK4 * A.bottomRows(3);
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
            _get_triangle_normal(&verts1(k, j, i, 0, 0), p0, pp5, t[2], Axis::XPOSITIVE);
            _get_triangle_normal(&verts1(k, j, i, 0, 0), p1, pp5, t[0], Axis::XPOSITIVE);
            _get_triangle_normal(&verts1(k, j, i, 0, 0), pp5, p4, t[1], Axis::YPOSITIVE);
            ratio[0] = _get_ratio(n1, t[0]);
            ratio[1] = _get_ratio(n4, t[1]);
            ratio[2] = _get_ratio(n0, t[2]);
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
            Eigen::RowVectorXd r = (-ratio[0] * v1 + ratio[1] * v4) * matK5 * A.block(6, 0, 3, 12);
            B[i5] -= pb * (r[0] + r[1] + r[2] + r[3]);
            C(i5, i1) += -r[2] + r[6] + r[10];
            C(i5, i2) += -r[3] - r[6] + r[8];
            C(i5, i5) += -r[1] + r[4] - r[10];
            C(i5, i6) += -r[0] - r[4] - r[8];
            r = -(ratio[2] * v0 + ratio[1] * v4) * matK6 * A.bottomRows(3);
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
            _get_triangle_normal(&verts1(k, j, i - 1, 1, 0), p0, pp5, t[2], Axis::XPOSITIVE);
            _get_triangle_normal(&verts1(k, j - 1, i - 1, 3, 0), p1, pp5, t[0], Axis::XPOSITIVE);
            _get_triangle_normal(&verts1(k, j, i - 1, 1, 0), pp5, p7, t[3], Axis::YPOSITIVE);
            ratio[0] = _get_ratio(n1, t[0]);
            ratio[2] = _get_ratio(n0, t[2]);
            ratio[3] = _get_ratio(n7, t[3]);
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
            Eigen::RowVectorXd r = (ratio[0] * v1 + ratio[3] * v7) * matK4 * A.block(6, 0, 3, 12);
            B[i4] -= pb * (r[0] + r[1] + r[2] + r[3]);
            C(i4, i0) += -r[2] + r[4] + r[10];
            C(i4, i3) += -r[3] - r[4] + r[8];
            C(i4, i4) += -r[1] + r[6] - r[10];
            C(i4, i7) += -r[0] - r[6] - r[8];
            r = (ratio[2] * v0 - ratio[3] * v7) * matK7 * A.bottomRows(3);
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
            _get_triangle_normal(&verts1(k, j, i, 0, 0), p0, pp5, t[2], Axis::XPOSITIVE);
            ratio[2] = _get_ratio(n0, t[2]);
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
            Eigen::RowVectorXd r = (-ratio[2] * v0) * matK6 * A.block(6, 0, 3, 12);
            C(i6, i2) += r[2] + r[8];
            C(i6, i3) += -r[2] + r[10];
            C(i6, i6) += r[0] - r[8];
            C(i6, i7) += -r[0] - r[10];
            r = (ratio[2] * v0) * matK7 * A.bottomRows(3);
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
            _get_triangle_normal(&verts1(k, j - 1, i, 2, 0), p1, pp5, t[0], Axis::XPOSITIVE);
            ratio[0] = _get_ratio(n1, t[0]);
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
            Eigen::RowVectorXd r = (ratio[0] * v1) * matK4 * A.block(6, 0, 3, 12);
            C(i4, i0) += r[2] + r[8];
            C(i4, i1) += -r[2] + r[10];
            C(i4, i4) += r[0] - r[8];
            C(i4, i5) += -r[0] - r[10];
            r = (-ratio[0] * v1) * matK5 * A.bottomRows(3);
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
            _get_triangle_normal(&verts1(k, j, i, 0, 0), p0, pp5, t[2], Axis::XPOSITIVE);
            _get_triangle_normal(&verts1(k, j, i, 0, 0), pp5, p1, t[0], Axis::XPOSITIVE);
            _get_triangle_normal(&verts1(k, j, i, 0, 0), p4, pp5, t[1], Axis::YPOSITIVE);
            _get_triangle_normal(&verts1(k, j, i, 0, 0), pp5, p7, t[3], Axis::YPOSITIVE);
            ratio[0] = _get_ratio(n1, t[0]);
            ratio[1] = _get_ratio(n4, t[1]);
            ratio[2] = _get_ratio(n0, t[2]);
            ratio[3] = _get_ratio(n7, t[3]);
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
            Eigen::RowVectorXd r = (ratio[0] * v1 + ratio[3] * v7) * matK4 * A.topRows(3);
            C(i4, i4) += -r[2] + r[6];
            C(i4, i5) += r[2] + r[4];
            C(i4, i6) += r[0] - r[4];
            C(i4, i7) += -r[0] - r[6];
            r = (-ratio[0] * v1 + ratio[1] * v4) * matK5 * A.block(3, 0, 3, 12);
            C(i5, i4) += -r[2] + r[6];
            C(i5, i5) += r[2] + r[4];
            C(i5, i6) += r[0] - r[4];
            C(i5, i7) += -r[0] - r[6];
            r = (-ratio[2] * v0 - ratio[1] * v4) * matK6 * A.block(6, 0, 3, 12);
            C(i6, i4) += -r[2] + r[6];
            C(i6, i5) += r[2] + r[4];
            C(i6, i6) += r[0] - r[4];
            C(i6, i7) += -r[0] - r[6];
            r = (ratio[2] * v0 - ratio[3] * v7) * matK7 * A.bottomRows(3);
            C(i7, i4) += -r[2] + r[6];
            C(i7, i5) += r[2] + r[4];
            C(i7, i6) += r[0] - r[4];
            C(i7, i7) += -r[0] - r[6];
            // To be implemented
        }
        else if (k == nz)
        {
            break;
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
            _get_triangle_normal(&verts1(k, j, i, 0, 0), p[0], pp[5], t[2], Axis::XPOSITIVE);
            _get_triangle_normal(&verts1(k, j, i, 0, 0), pp[5], p[1], t[0], Axis::XPOSITIVE);
            _get_triangle_normal(&verts1(k, j, i, 0, 0), p[4], pp[5], t[1], Axis::YPOSITIVE);
            _get_triangle_normal(&verts1(k, j, i, 0, 0), pp[5], p[7], t[3], Axis::YPOSITIVE);
            ratio[0] = _get_ratio(n[1], t[0]);
            ratio[1] = _get_ratio(n[4], t[1]);
            ratio[2] = _get_ratio(n[0], t[2]);
            ratio[3] = _get_ratio(n[7], t[3]);
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
            Eigen::RowVectorXd r = (ratio[0] * v[1] + ratio[3] * v[7]) * matK[4] * A.block(12, 0, 3, 24);
            C(idx[4], idx[0]) += r[4] + r[12] + r[20];
            C(idx[4], idx[1]) += -r[4] + r[10] + r[22];
            C(idx[4], idx[2]) += r[6] - r[10] + r[16];
            C(idx[4], idx[3]) += -r[6] - r[12] + r[18];
            C(idx[4], idx[4]) += r[2] + r[14] - r[20];
            C(idx[4], idx[5]) += -r[2] + r[8] - r[22];
            C(idx[4], idx[6]) += r[0] - r[8] - r[16];
            C(idx[4], idx[7]) += -r[0] - r[14] - r[18];
            r = (-ratio[0] * v[1] + ratio[1] * v[4]) * matK[5] * A.block(15, 0, 3, 24);
            C(idx[5], idx[0]) += r[4] + r[12] + r[20];
            C(idx[5], idx[1]) += -r[4] + r[10] + r[22];
            C(idx[5], idx[2]) += r[6] - r[10] + r[16];
            C(idx[5], idx[3]) += -r[6] - r[12] + r[18];
            C(idx[5], idx[4]) += r[2] + r[14] - r[20];
            C(idx[5], idx[5]) += -r[2] + r[8] - r[22];
            C(idx[5], idx[6]) += r[0] - r[8] - r[16];
            C(idx[5], idx[7]) += -r[0] - r[14] - r[18];
            r = (-ratio[2] * v[0] - ratio[1] * v[4]) * matK[6] * A.block(18, 0, 3, 24);
            C(idx[6], idx[0]) += r[4] + r[12] + r[20];
            C(idx[6], idx[1]) += -r[4] + r[10] + r[22];
            C(idx[6], idx[2]) += r[6] - r[10] + r[16];
            C(idx[6], idx[3]) += -r[6] - r[12] + r[18];
            C(idx[6], idx[4]) += r[2] + r[14] - r[20];
            C(idx[6], idx[5]) += -r[2] + r[8] - r[22];
            C(idx[6], idx[6]) += r[0] - r[8] - r[16];
            C(idx[6], idx[7]) += -r[0] - r[14] - r[18];
            r = (ratio[2] * v[0] - ratio[3] * v[7]) * matK[7] * A.block(21, 0, 3, 24);
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
    case Axis::ZNEGATIVE:
        assert(k > 0);
        if (i == 0 && j == 0 && k == 0)
        {
            break;
        }
        else if (i == 0 && j == ny && k == 0)
        {
            break;
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
            _get_triangle_normal(&verts1(k - 1, j, i, 4, 0), pp4, p3, t[2], Axis::XPOSITIVE);
            ratio[2] = _get_ratio(n3, t[2]);
            _cross_product(n5, n8, n58);
            Eigen::Matrix3d matK;
            _get_matK(pem1[cur], matK);
            Eigen::RowVector3d Dx;
            Eigen::Vector3d Dn;
            Dx << p3[0] - _cx(k - 1, j, i), p3[1] - _cy(k - 1, j, i), p3[2] - _cz(k - 1, j, i);
            Dn << n58[0], n58[1], n58[2];
            double cA = Dx * matK.inverse() * Dn;
            cA = _dot_product(n3, n58) / cA;
            cA *= ratio[2];
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
            _get_triangle_normal(&verts1(k - 1, j - 1, i, 6, 0), p2, pp4, t[0], Axis::XPOSITIVE);
            ratio[0] = _get_ratio(n2, t[0]);
            _cross_product(n5, n11, n511);
            Eigen::Matrix3d matK;
            _get_matK(pem1[cur], matK);
            Eigen::RowVector3d Dx;
            Eigen::Vector3d Dn;
            Dx << p2[0] - _cx(k - 1, j - 1, i), p2[1] - _cy(k - 1, j - 1, i), p2[2] - _cz(k - 1, j - 1, i);
            Dn << n511[0], n511[1], n511[2];
            double cA = Dx * matK.inverse() * Dn;
            cA = _dot_product(n2, n511) / cA;
            cA *= ratio[0];
            B[cur] += pb * cA;
            C(cur, cur) += cA;
        }
        else if (i == nx && j == 0 && k == 0)
        {
            break;
        }
        else if (i == nx && j == ny && k == 0)
        {
            break;
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
            _get_triangle_normal(&verts1(k - 1, j, i - 1, 5, 0), pp4, p3, t[2], Axis::XPOSITIVE);
            ratio[2] = _get_ratio(n3, t[2]);
            _cross_product(n6, n9, n69);
            Eigen::Matrix3d matK;
            _get_matK(pem1[cur], matK);
            Eigen::RowVector3d Dx;
            Eigen::Vector3d Dn;
            Dx << p3[0] - _cx(k - 1, j, i - 1), p3[1] - _cy(k - 1, j, i - 1), p3[2] - _cz(k - 1, j, i - 1);
            Dn << n69[0], n69[1], n69[2];
            double cA = Dx * matK.inverse() * Dn;
            cA = _dot_product(n3, n69) / cA;
            cA *= ratio[2];
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
            _get_triangle_normal(&verts1(k - 1, j - 1, i - 1, 7, 0), p2, pp4, t[0], Axis::XPOSITIVE);
            ratio[0] = _get_ratio(n2, t[0]);
            _cross_product(n6, n10, n610);
            Eigen::Matrix3d matK;
            _get_matK(pem1[cur], matK);
            Eigen::RowVector3d Dx;
            Eigen::Vector3d Dn;
            Dx << p2[0] - _cx(k - 1, j - 1, i - 1), p2[1] - _cy(k - 1, j - 1, i - 1), p2[2] - _cz(k - 1, j - 1, i - 1);
            Dn << n610[0], n610[1], n610[2];
            double cA = Dx * matK.inverse() * Dn;
            cA = _dot_product(n2, n610) / cA;
            cA *= ratio[0];
            B[cur] -= pb * cA;
            C(cur, cur) -= cA;
        }
        // 12条边上其余点的处理
        else if (j == 0 && k == 0)
        {
            break;
        }
        else if (j == ny && k == 0)
        {
            break;
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
            _get_triangle_normal(&verts1(k - 1, j, i, 4, 0), pp4, p3, t[2], Axis::XPOSITIVE);
            ratio[2] = _get_ratio(n3, t[2]);
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
            cA *= ratio[2];
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
            _get_triangle_normal(&verts1(k - 1, j - 1, i, 6, 0), p2, pp4, t[0], Axis::XPOSITIVE);
            ratio[0] = _get_ratio(n2, t[0]);
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
            cA *= ratio[0];
            C(i1, i1) += cA;
            C(i0, i0) += cA;
            C(i1, i0) -= cA;
            C(i0, i1) -= cA;
        }
        else if (i == 0 && k == 0)
        {
            break;
        }
        else if (i == nx && k == 0)
        {
            break;
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
            _get_triangle_normal(&verts1(k - 1, j - 1, i, 6, 0), p2, pp4, t[0], Axis::XPOSITIVE);
            _get_triangle_normal(&verts1(k - 1, j, i, 4, 0), p3, pp4, t[2], Axis::XPOSITIVE);
            _get_triangle_normal(&verts1(k - 1, j, i, 4, 0), pp4, p5, t[1], Axis::YPOSITIVE);
            ratio[0] = _get_ratio(n2, t[0]);
            ratio[1] = _get_ratio(n5, t[1]);
            ratio[2] = _get_ratio(n3, t[2]);
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
            Eigen::RowVectorXd r = (ratio[1] * v5 - ratio[0] * v2) * matK1 * A.topRows(3);
            B[i1] -= pb * (r[0] + r[1]);
            C(i1, i1) += -r[1] + r[2];
            C(i1, i2) += -r[0] - r[2];
            r = -(ratio[2] * v3 + ratio[1] * v5) * matK2 * A.bottomRows(3);
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
            _get_triangle_normal(&verts1(k - 1, j - 1, i - 1, 7, 0), p2, pp4, t[0], Axis::XPOSITIVE);
            _get_triangle_normal(&verts1(k - 1, j, i - 1, 5, 0), p3, pp4, t[2], Axis::XPOSITIVE);
            _get_triangle_normal(&verts1(k - 1, j, i - 1, 5, 0), pp4, p6, t[3], Axis::YPOSITIVE);
            ratio[0] = _get_ratio(n2, t[0]);
            ratio[2] = _get_ratio(n3, t[2]);
            ratio[3] = _get_ratio(n6, t[3]);
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
            Eigen::RowVectorXd r = (ratio[3] * v6 + ratio[0] * v2) * matK0 * A.topRows(3);
            B[i0] -= pb * (r[0] + r[1]);
            C(i0, i0) += -r[1] + r[2];
            C(i0, i3) += -r[0] - r[2];
            r = (ratio[2] * v3 - ratio[3] * v6) * matK3 * A.bottomRows(3);
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
            _get_triangle_normal(&verts1(k - 1, j, i, 4, 0), p3, pp4, t[2], Axis::XPOSITIVE);
            ratio[2] = _get_ratio(n3, t[2]);
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
            Eigen::RowVectorXd r = (-ratio[2] * v3) * matK2 * A.topRows(3);
            B[i2] -= pb * (r[0] + r[1]);
            C(i2, i2) += -r[1] + r[2];
            C(i2, i6) += -r[0] - r[2];
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
            _get_triangle_normal(&verts1(k - 1, j, i - 1, 5, 0), p3, pp4, t[2], Axis::XPOSITIVE);
            ratio[2] = _get_ratio(n3, t[2]);
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
            Eigen::RowVectorXd r = (ratio[2] * v3) * matK3 * A.topRows(3);
            B[i3] -= pb * (r[0] + r[1]);
            C(i3, i3) += -r[1] + r[2];
            C(i3, i7) += -r[0] - r[2];
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
            _get_triangle_normal(&verts1(k - 1, j - 1, i, 6, 0), pp4, p2, t[0], Axis::XPOSITIVE);
            ratio[0] = _get_ratio(n2, t[0]);
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
            Eigen::RowVectorXd r = (-ratio[0] * v2) * matK1 * A.topRows(3);
            B[i1] -= pb * (r[0] + r[1]);
            C(i1, i1) += -r[1] + r[2];
            C(i1, i5) += -r[0] - r[2];
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
            _get_triangle_normal(&verts1(k - 1, j - 1, i - 1, 7, 0), pp4, p2, t[0], Axis::XPOSITIVE);
            ratio[0] = _get_ratio(n2, t[0]);
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
            Eigen::RowVectorXd r = (ratio[0] * v2) * matK0 * A.topRows(3);
            B[i0] -= pb * (r[0] + r[1]);
            C(i0, i0) += -r[1] + r[2];
            C(i0, i4) += -r[0] - r[2];
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
            _get_triangle_normal(&verts1(k - 1, j - 1, i, 6, 0), p2, pp4, t[0], Axis::XPOSITIVE);
            _get_triangle_normal(&verts1(k - 1, j, i, 4, 0), p3, pp4, t[2], Axis::XPOSITIVE);
            _get_triangle_normal(&verts1(k - 1, j, i, 4, 0), pp4, p5, t[1], Axis::YPOSITIVE);
            ratio[0] = _get_ratio(n2, t[0]);
            ratio[1] = _get_ratio(n5, t[1]);
            ratio[2] = _get_ratio(n3, t[2]);
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
            Eigen::RowVectorXd r = (-ratio[0] * v2 + ratio[1] * v5) * matK1 * A.topRows(3);
            B[i1] -= pb * (r[0] + r[1] + r[2] + r[3]);
            C(i1, i1) += -r[2] + r[6] + r[10];
            C(i1, i2) += -r[3] - r[6] + r[8];
            C(i1, i5) += -r[1] + r[4] - r[10];
            C(i1, i6) += -r[0] - r[4] - r[8];
            r = (-ratio[2] * v3 - ratio[1] * v5) * matK2 * A.block(3, 0, 3, 12);
            B[i2] -= pb * (r[0] + r[1] + r[2] + r[3]);
            C(i2, i1) += -r[2] + r[6] + r[10];
            C(i2, i2) += -r[3] - r[6] + r[8];
            C(i2, i5) += -r[1] + r[4] - r[10];
            C(i2, i6) += -r[0] - r[4] - r[8];
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
            _get_triangle_normal(&verts1(k - 1, j - 1, i - 1, 7, 0), p2, pp4, t[0], Axis::XPOSITIVE);
            _get_triangle_normal(&verts1(k - 1, j, i - 1, 5, 0), p3, pp4, t[2], Axis::XPOSITIVE);
            _get_triangle_normal(&verts1(k - 1, j, i - 1, 5, 0), pp4, p6, t[3], Axis::YPOSITIVE);
            ratio[0] = _get_ratio(n2, t[0]);
            ratio[2] = _get_ratio(n3, t[2]);
            ratio[3] = _get_ratio(n6, t[3]);
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
            Eigen::RowVectorXd r = (ratio[0] * v2 + ratio[3] * v6) * matK0 * A.topRows(3);
            B[i0] -= pb * (r[0] + r[1] + r[2] + r[3]);
            C(i0, i0) += -r[2] + r[4] + r[10];
            C(i0, i3) += -r[3] - r[4] + r[8];
            C(i0, i4) += -r[1] + r[6] - r[10];
            C(i0, i7) += -r[0] - r[6] - r[8];
            r = (ratio[2] * v3 - ratio[3] * v6) * matK3 * A.block(3, 0, 3, 12);
            B[i3] -= pb * (r[0] + r[1] + r[2] + r[3]);
            C(i3, i0) += -r[2] + r[4] + r[10];
            C(i3, i3) += -r[3] - r[4] + r[8];
            C(i3, i4) += -r[1] + r[6] - r[10];
            C(i3, i7) += -r[0] - r[6] - r[8];
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
            _get_triangle_normal(&verts1(k, j, i, 0, 0), p3, pp4, t[2], Axis::XPOSITIVE);
            ratio[2] = _get_ratio(n3, t[2]);
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
            Eigen::RowVectorXd r = (-ratio[2] * v3) * matK2 * A.topRows(3);
            C(i2, i2) += r[2] + r[8];
            C(i2, i3) += -r[2] + r[10];
            C(i2, i6) += r[0] - r[8];
            C(i2, i7) += -r[0] - r[10];
            r = (ratio[2] * v3) * matK3 * A.block(3, 0, 3, 12);
            C(i3, i2) += r[2] + r[8];
            C(i3, i3) += -r[2] + r[10];
            C(i3, i6) += r[0] - r[8];
            C(i3, i7) += -r[0] - r[10];
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
            _get_triangle_normal(&verts1(k, j - 1, i, 2, 0), p2, pp4, t[0], Axis::XPOSITIVE);
            ratio[0] = _get_ratio(n2, t[0]);
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
            Eigen::RowVectorXd r = (ratio[0] * v2) * matK0 * A.topRows(3);
            C(i0, i0) += r[2] + r[8];
            C(i0, i1) += -r[2] + r[10];
            C(i0, i4) += r[0] - r[8];
            C(i0, i5) += -r[0] - r[10];
            r = (-ratio[0] * v2) * matK1 * A.block(3, 0, 3, 12);
            C(i1, i0) += r[2] + r[8];
            C(i1, i1) += -r[2] + r[10];
            C(i1, i4) += r[0] - r[8];
            C(i1, i5) += -r[0] - r[10];
            // To be implemented
        }
        else if (k == 0)
        {
            break;
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
            _get_triangle_normal(&verts1(k - 1, j, i, 4, 0), p2, pp4, t[0], Axis::XPOSITIVE);
            _get_triangle_normal(&verts1(k - 1, j, i, 4, 0), pp4, p3, t[2], Axis::XPOSITIVE);
            _get_triangle_normal(&verts1(k - 1, j, i, 4, 0), p5, pp4, t[1], Axis::YPOSITIVE);
            _get_triangle_normal(&verts1(k - 1, j, i, 4, 0), pp4, p6, t[3], Axis::YPOSITIVE);
            ratio[0] = _get_ratio(n2, t[0]);
            ratio[1] = _get_ratio(n5, t[1]);
            ratio[2] = _get_ratio(n3, t[2]);
            ratio[3] = _get_ratio(n6, t[3]);
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
            Eigen::RowVectorXd r = (ratio[0] * v2 + ratio[3] * v6) * matK0 * A.topRows(3);
            C(i0, i0) += r[0] + r[6];
            C(i0, i1) += -r[0] + r[4];
            C(i0, i2) += -r[2] - r[4];
            C(i0, i3) += r[2] - r[6];
            r = (-ratio[0] * v2 + ratio[1] * v5) * matK1 * A.block(3, 0, 3, 12);
            C(i1, i0) += r[0] + r[6];
            C(i1, i1) += -r[0] + r[4];
            C(i1, i2) += -r[2] - r[4];
            C(i1, i3) += r[2] - r[6];
            r = (-ratio[2] * v3 - ratio[1] * v5) * matK2 * A.block(6, 0, 3, 12);
            C(i2, i0) += r[0] + r[6];
            C(i2, i1) += -r[0] + r[4];
            C(i2, i2) += -r[2] - r[4];
            C(i2, i3) += r[2] - r[6];
            r = (ratio[2] * v3 - ratio[3] * v6) * matK3 * A.bottomRows(3);
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
            _get_triangle_normal(&verts1(k - 1, j, i, 4, 0), p[2], pp[4], t[0], Axis::XPOSITIVE);
            _get_triangle_normal(&verts1(k - 1, j, i, 4, 0), pp[4], p[3], t[2], Axis::XPOSITIVE);
            _get_triangle_normal(&verts1(k - 1, j, i, 4, 0), p[5], pp[4], t[1], Axis::YPOSITIVE);
            _get_triangle_normal(&verts1(k - 1, j, i, 4, 0), pp[4], p[6], t[3], Axis::YPOSITIVE);
            ratio[0] = _get_ratio(n[2], t[0]);
            ratio[1] = _get_ratio(n[5], t[1]);
            ratio[2] = _get_ratio(n[3], t[2]);
            ratio[3] = _get_ratio(n[6], t[3]);
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
            Eigen::RowVectorXd r = (ratio[0] * v[2] + ratio[3] * v[6]) * matK[0] * A.topRows(3);
            C(idx[0], idx[0]) += r[4] + r[12] + r[20];
            C(idx[0], idx[1]) += -r[4] + r[10] + r[22];
            C(idx[0], idx[2]) += r[6] - r[10] + r[16];
            C(idx[0], idx[3]) += -r[6] - r[12] + r[18];
            C(idx[0], idx[4]) += r[2] + r[14] - r[20];
            C(idx[0], idx[5]) += -r[2] + r[8] - r[22];
            C(idx[0], idx[6]) += r[0] - r[8] - r[16];
            C(idx[0], idx[7]) += -r[0] - r[14] - r[18];
            r = (-ratio[0] * v[2] + ratio[1] * v[5]) * matK[1] * A.block(3, 0, 3, 24);
            C(idx[1], idx[0]) += r[4] + r[12] + r[20];
            C(idx[1], idx[1]) += -r[4] + r[10] + r[22];
            C(idx[1], idx[2]) += r[6] - r[10] + r[16];
            C(idx[1], idx[3]) += -r[6] - r[12] + r[18];
            C(idx[1], idx[4]) += r[2] + r[14] - r[20];
            C(idx[1], idx[5]) += -r[2] + r[8] - r[22];
            C(idx[1], idx[6]) += r[0] - r[8] - r[16];
            C(idx[1], idx[7]) += -r[0] - r[14] - r[18];
            r = (-ratio[2] * v[3] - ratio[1] * v[5]) * matK[2] * A.block(6, 0, 3, 24);
            C(idx[2], idx[0]) += r[4] + r[12] + r[20];
            C(idx[2], idx[1]) += -r[4] + r[10] + r[22];
            C(idx[2], idx[2]) += r[6] - r[10] + r[16];
            C(idx[2], idx[3]) += -r[6] - r[12] + r[18];
            C(idx[2], idx[4]) += r[2] + r[14] - r[20];
            C(idx[2], idx[5]) += -r[2] + r[8] - r[22];
            C(idx[2], idx[6]) += r[0] - r[8] - r[16];
            C(idx[2], idx[7]) += -r[0] - r[14] - r[18];
            r = (ratio[2] * v[3] - ratio[3] * v[6]) * matK[3] * A.block(9, 0, 3, 24);
            C(idx[3], idx[0]) += r[4] + r[12] + r[20];
            C(idx[3], idx[1]) += -r[4] + r[10] + r[22];
            C(idx[3], idx[2]) += r[6] - r[10] + r[16];
            C(idx[3], idx[3]) += -r[6] - r[12] + r[18];
            C(idx[3], idx[4]) += r[2] + r[14] - r[20];
            C(idx[3], idx[5]) += -r[2] + r[8] - r[22];
            C(idx[3], idx[6]) += r[0] - r[8] - r[16];
            C(idx[3], idx[7]) += -r[0] - r[14] - r[18];
            // 内部角点
        }
        break;
    }
}

// axis: 0 for x, 1 for y, 2 for z
void TPFA(int i, int j, int k, int axis)
{
    // N1 and N2 are the indices of the two neighboring cells along the specified axis. t and p are the normal vectors and centroids of the  interface between the two cells, respectively.
    int cur = k * nx * ny + j * nx + i, N1 = 0, N2 = 0;
    double t[3]{};
    double p[3]{};
    double lamda = 0, lamda1 = 0, lamda2 = 0;
    Eigen::RowVector3d c1, c2, nT;
    Eigen::Matrix3d matK1, matK2;
    if (axis == 0)
    {
        if (j > 0 & j < ny && k > 0)
        {
            N1 = (k - 1) * nx * ny + (j - 1) * nx + i;
            N2 = (k - 1) * nx * ny + j * nx + i;
            _get_centroid(&verts1(k - 1, j, i, 0, 0), &verts1(k - 1, j, i, 1, 0), &verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 5, 0), p);
            _get_triangle_normal(&verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 5, 0), p, t, Axis::YPOSITIVE);
            _get_matK(pem1[N1], matK1);
            _get_matK(pem1[N2], matK2);
            c1 << p[0] - _cx(k - 1, j - 1, i), p[1] - _cy(k - 1, j - 1, i), p[2] - _cz(k - 1, j - 1, i);
            c2 << _cx(k - 1, j, i) - p[0], _cy(k - 1, j, i) - p[1], _cz(k - 1, j, i) - p[2];
            nT << t[0], t[1], t[2];
            lamda1 = (nT * matK1).dot(c1) / c1.squaredNorm();
            lamda2 = (nT * matK2).dot(c2) / c2.squaredNorm();
            lamda = lamda1 * lamda2 / (lamda1 + lamda2);
            C(N1, N1) -= lamda;
            C(N1, N2) += lamda;
            C(N2, N1) += lamda;
            C(N2, N2) -= lamda;
        }
        if (j > 0 & j < ny && k < nz)
        {
            N1 = k * nx * ny + (j - 1) * nx + i;
            N2 = k * nx * ny + j * nx + i;
            _get_centroid(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 1, 0), &verts1(k, j, i, 4, 0), &verts1(k, j, i, 5, 0), p);
            _get_triangle_normal(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 1, 0), p, t, Axis::YPOSITIVE);
            _get_matK(pem1[N1], matK1);
            _get_matK(pem1[N2], matK2);
            c1 << p[0] - _cx(k, j - 1, i), p[1] - _cy(k, j - 1, i), p[2] - _cz(k, j - 1, i);
            c2 << _cx(k, j, i) - p[0], _cy(k, j, i) - p[1], _cz(k, j, i) - p[2];
            nT << t[0], t[1], t[2];
            lamda1 = (nT * matK1).dot(c1) / c1.squaredNorm();
            lamda2 = (nT * matK2).dot(c2) / c2.squaredNorm();
            lamda = lamda1 * lamda2 / (lamda1 + lamda2);
            C(N1, N1) -= lamda;
            C(N1, N2) += lamda;
            C(N2, N1) += lamda;
            C(N2, N2) -= lamda;
        }
        if (k > 0 && k < nz && j > 0)
        {
            N1 = (k - 1) * nx * ny + (j - 1) * nx + i;
            N2 = k * nx * ny + (j - 1) * nx + i;
            _get_centroid(&verts1(k, j - 1, i, 0, 0), &verts1(k, j - 1, i, 1, 0), &verts1(k, j - 1, i, 2, 0), &verts1(k, j - 1, i, 3, 0), p);
            _get_triangle_normal(&verts1(k, j - 1, i, 2, 0), &verts1(k, j - 1, i, 3, 0), p, t, Axis::ZPOSITIVE);
            _get_matK(pem1[N1], matK1);
            _get_matK(pem1[N2], matK2);
            c1 << p[0] - _cx(k - 1, j - 1, i), p[1] - _cy(k - 1, j - 1, i), p[2] - _cz(k - 1, j - 1, i);
            c2 << _cx(k, j - 1, i) - p[0], _cy(k, j - 1, i) - p[1], _cz(k, j - 1, i) - p[2];
            nT << t[0], t[1], t[2];
            lamda1 = (nT * matK1).dot(c1) / c1.squaredNorm();
            lamda2 = (nT * matK2).dot(c2) / c2.squaredNorm();
            lamda = lamda1 * lamda2 / (lamda1 + lamda2);
            C(N1, N1) -= lamda;
            C(N1, N2) += lamda;
            C(N2, N1) += lamda;
            C(N2, N2) -= lamda;
        }
        if (k > 0 && k < nz && j < ny)
        {
            N1 = (k - 1) * nx * ny + j * nx + i;
            N2 = k * nx * ny + j * nx + i;
            _get_centroid(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 1, 0), &verts1(k, j, i, 2, 0), &verts1(k, j, i, 3, 0), p);
            _get_triangle_normal(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 1, 0), p, t, Axis::ZPOSITIVE);
            _get_matK(pem1[N1], matK1);
            _get_matK(pem1[N2], matK2);
            c1 << p[0] - _cx(k - 1, j, i), p[1] - _cy(k - 1, j, i), p[2] - _cz(k - 1, j, i);
            c2 << _cx(k, j, i) - p[0], _cy(k, j, i) - p[1], _cz(k, j, i) - p[2];
            nT << t[0], t[1], t[2];
            lamda1 = (nT * matK1).dot(c1) / c1.squaredNorm();
            lamda2 = (nT * matK2).dot(c2) / c2.squaredNorm();
            lamda = lamda1 * lamda2 / (lamda1 + lamda2);
            C(N1, N1) -= lamda;
            C(N1, N2) += lamda;
            C(N2, N1) += lamda;
            C(N2, N2) -= lamda;
        }
    }
    else if (axis == 1)
    {
        if (i > 0 && i < nx && k > 0)
        {
            N1 = (k - 1) * nx * ny + j * nx + i - 1;
            N2 = (k - 1) * nx * ny + j * nx + i;
            _get_centroid(&verts1(k - 1, j, i, 0, 0), &verts1(k - 1, j, i, 2, 0), &verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 6, 0), p);
            _get_triangle_normal(&verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 6, 0), p, t, Axis::XPOSITIVE);
            _get_matK(pem1[N1], matK1);
            _get_matK(pem1[N2], matK2);
            c1 << p[0] - _cx(k - 1, j, i - 1), p[1] - _cy(k - 1, j, i - 1), p[2] - _cz(k - 1, j, i - 1);
            c2 << _cx(k - 1, j, i) - p[0], _cy(k - 1, j, i) - p[1], _cz(k - 1, j, i) - p[2];
            nT << t[0], t[1], t[2];
            lamda1 = (nT * matK1).dot(c1) / c1.squaredNorm();
            lamda2 = (nT * matK2).dot(c2) / c2.squaredNorm();
            lamda = lamda1 * lamda2 / (lamda1 + lamda2);
            C(N1, N1) -= lamda;
            C(N1, N2) += lamda;
            C(N2, N1) += lamda;
            C(N2, N2) -= lamda;
        }
        if (i > 0 && i < nx && k < nz)
        {
            N1 = k * nx * ny + j * nx + i - 1;
            N2 = k * nx * ny + j * nx + i;
            _get_centroid(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 2, 0), &verts1(k, j, i, 4, 0), &verts1(k, j, i, 6, 0), p);
            _get_triangle_normal(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 2, 0), p, t, Axis::XPOSITIVE);
            _get_matK(pem1[N1], matK1);
            _get_matK(pem1[N2], matK2);
            c1 << p[0] - _cx(k, j, i - 1), p[1] - _cy(k, j, i - 1), p[2] - _cz(k, j, i - 1);
            c2 << _cx(k, j, i) - p[0], _cy(k, j, i) - p[1], _cz(k, j, i) - p[2];
            nT << t[0], t[1], t[2];
            lamda1 = (nT * matK1).dot(c1) / c1.squaredNorm();
            lamda2 = (nT * matK2).dot(c2) / c2.squaredNorm();
            lamda = lamda1 * lamda2 / (lamda1 + lamda2);
            C(N1, N1) -= lamda;
            C(N1, N2) += lamda;
            C(N2, N1) += lamda;
            C(N2, N2) -= lamda;
        }
        if (k > 0 && k < nz && i > 0)
        {
            N1 = (k - 1) * nx * ny + j * nx + i - 1;
            N2 = k * nx * ny + j * nx + i - 1;
            _get_centroid(&verts1(k, j, i - 1, 0, 0), &verts1(k, j, i - 1, 1, 0), &verts1(k, j, i - 1, 2, 0), &verts1(k, j, i - 1, 3, 0), p);
            _get_triangle_normal(&verts1(k, j, i - 1, 1, 0), &verts1(k, j, i - 1, 3, 0), p, t, Axis::ZPOSITIVE);
            _get_matK(pem1[N1], matK1);
            _get_matK(pem1[N2], matK2);
            c1 << p[0] - _cx(k - 1, j, i - 1), p[1] - _cy(k - 1, j, i - 1), p[2] - _cz(k - 1, j, i - 1);
            c2 << _cx(k, j, i - 1) - p[0], _cy(k, j, i - 1) - p[1], _cz(k, j, i - 1) - p[2];
            nT << t[0], t[1], t[2];
            lamda1 = (nT * matK1).dot(c1) / c1.squaredNorm();
            lamda2 = (nT * matK2).dot(c2) / c2.squaredNorm();
            lamda = lamda1 * lamda2 / (lamda1 + lamda2);
            C(N1, N1) -= lamda;
            C(N1, N2) += lamda;
            C(N2, N1) += lamda;
            C(N2, N2) -= lamda;
        }
        if (k > 0 && k < nz && i < nx)
        {
            N1 = (k - 1) * nx * ny + j * nx + i;
            N2 = k * nx * ny + j * nx + i;
            _get_centroid(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 1, 0), &verts1(k, j, i, 2, 0), &verts1(k, j, i, 3, 0), p);
            _get_triangle_normal(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 2, 0), p, t, Axis::ZPOSITIVE);
            _get_matK(pem1[N1], matK1);
            _get_matK(pem1[N2], matK2);
            c1 << p[0] - _cx(k - 1, j, i), p[1] - _cy(k - 1, j, i), p[2] - _cz(k - 1, j, i);
            c2 << _cx(k, j, i) - p[0], _cy(k, j, i) - p[1], _cz(k, j, i) - p[2];
            nT << t[0], t[1], t[2];
            lamda1 = (nT * matK1).dot(c1) / c1.squaredNorm();
            lamda2 = (nT * matK2).dot(c2) / c2.squaredNorm();
            lamda = lamda1 * lamda2 / (lamda1 + lamda2);
            C(N1, N1) -= lamda;
            C(N1, N2) += lamda;
            C(N2, N1) += lamda;
            C(N2, N2) -= lamda;
        }
    }
    else
    {
        if (i > 0 && i < nx && j > 0)
        {
            N1 = k * nx * ny + (j - 1) * nx + i - 1;
            N2 = k * nx * ny + (j - 1) * nx + i;
            _get_centroid(&verts1(k, j - 1, i, 0, 0), &verts1(k, j - 1, i, 2, 0), &verts1(k, j - 1, i, 4, 0), &verts1(k, j - 1, i, 6, 0), p);
            _get_triangle_normal(&verts1(k, j - 1, i, 2, 0), &verts1(k, j - 1, i, 6, 0), p, t, Axis::XPOSITIVE);
            _get_matK(pem1[N1], matK1);
            _get_matK(pem1[N2], matK2);
            c1 << p[0] - _cx(k, j - 1, i - 1), p[1] - _cy(k, j - 1, i - 1), p[2] - _cz(k, j - 1, i - 1);
            c2 << _cx(k, j - 1, i) - p[0], _cy(k, j - 1, i) - p[1], _cz(k, j - 1, i) - p[2];
            nT << t[0], t[1], t[2];
            lamda1 = (nT * matK1).dot(c1) / c1.squaredNorm();
            lamda2 = (nT * matK2).dot(c2) / c2.squaredNorm();
            lamda = lamda1 * lamda2 / (lamda1 + lamda2);
            C(N1, N1) -= lamda;
            C(N1, N2) += lamda;
            C(N2, N1) += lamda;
            C(N2, N2) -= lamda;
        }
        if (i > 0 && i < nx && j < ny)
        {
            N1 = k * nx * ny + j * nx + i - 1;
            N2 = k * nx * ny + j * nx + i;
            _get_centroid(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 2, 0), &verts1(k, j, i, 4, 0), &verts1(k, j, i, 6, 0), p);
            _get_triangle_normal(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 4, 0), p, t, Axis::XPOSITIVE);
            _get_matK(pem1[N1], matK1);
            _get_matK(pem1[N2], matK2);
            c1 << p[0] - _cx(k, j, i - 1), p[1] - _cy(k, j, i - 1), p[2] - _cz(k, j, i - 1);
            c2 << _cx(k, j, i) - p[0], _cy(k, j, i) - p[1], _cz(k, j, i) - p[2];
            nT << t[0], t[1], t[2];
            lamda1 = (nT * matK1).dot(c1) / c1.squaredNorm();
            lamda2 = (nT * matK2).dot(c2) / c2.squaredNorm();
            lamda = lamda1 * lamda2 / (lamda1 + lamda2);
            C(N1, N1) -= lamda;
            C(N1, N2) += lamda;
            C(N2, N1) += lamda;
            C(N2, N2) -= lamda;
        }
        if (j > 0 && j < ny && i > 0)
        {
            N1 = k * nx * ny + (j - 1) * nx + i - 1;
            N2 = k * nx * ny + j * nx + i - 1;
            _get_centroid(&verts1(k, j, i - 1, 0, 0), &verts1(k, j, i - 1, 1, 0), &verts1(k, j, i - 1, 4, 0), &verts1(k, j, i - 1, 5, 0), p);
            _get_triangle_normal(&verts1(k, j, i - 1, 1, 0), &verts1(k, j, i - 1, 5, 0), p, t, Axis::YPOSITIVE);
            _get_matK(pem1[N1], matK1);
            _get_matK(pem1[N2], matK2);
            c1 << p[0] - _cx(k, j - 1, i - 1), p[1] - _cy(k, j - 1, i - 1), p[2] - _cz(k, j - 1, i - 1);
            c2 << _cx(k, j, i - 1) - p[0], _cy(k, j, i - 1) - p[1], _cz(k, j, i - 1) - p[2];
            nT << t[0], t[1], t[2];
            lamda1 = (nT * matK1).dot(c1) / c1.squaredNorm();
            lamda2 = (nT * matK2).dot(c2) / c2.squaredNorm();
            lamda = lamda1 * lamda2 / (lamda1 + lamda2);
            C(N1, N1) -= lamda;
            C(N1, N2) += lamda;
            C(N2, N1) += lamda;
            C(N2, N2) -= lamda;
        }
        if (j > 0 && j < ny && i < nx)
        {
            N1 = k * nx * ny + (j - 1) * nx + i;
            N2 = k * nx * ny + j * nx + i;
            _get_centroid(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 1, 0), &verts1(k, j, i, 4, 0), &verts1(k, j, i, 5, 0), p);
            _get_triangle_normal(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 4, 0), p, t, Axis::YPOSITIVE);
            _get_matK(pem1[N1], matK1);
            _get_matK(pem1[N2], matK2);
            c1 << p[0] - _cx(k, j - 1, i), p[1] - _cy(k, j - 1, i), p[2] - _cz(k, j - 1, i);
            c2 << _cx(k, j, i) - p[0], _cy(k, j, i) - p[1], _cz(k, j, i) - p[2];
            nT << t[0], t[1], t[2];
            lamda1 = (nT * matK1).dot(c1) / c1.squaredNorm();
            lamda2 = (nT * matK2).dot(c2) / c2.squaredNorm();
            lamda = lamda1 * lamda2 / (lamda1 + lamda2);
            C(N1, N1) -= lamda;
            C(N1, N2) += lamda;
            C(N2, N1) += lamda;
            C(N2, N2) -= lamda;
        }
    }
}
void FAM(int *idx, const double lz, const double *l, const double *r, const double *angle, const Eigen::Matrix2d *K, double &alpha, const double *r_p, const double *angle_p)
{
    std::complex<double> c[4];
    for (int i = 0; i < 4; ++i)
    {
        c[i].real(-K[i](0, 1) / K[i](1, 1));
        c[i].imag(std::sqrt(K[i](0, 0) * K[i](1, 1) - K[i](0, 1) * K[i](1, 0)) / K[i](1, 1));
    }
    std::complex<double> zj[3], zj1[3], zp[4];
    double rj[3], rj1[3], rp[4], theta[4], theta1[3], thetap[4];
    for (int i = 0; i < 3; ++i)
    {
        zj[i] = cos(angle[i]) + c[i] * sin(angle[i]);
        zj1[i] = cos(angle[i]) + c[i + 1] * sin(angle[i]);
        rj[i] = std::abs(zj[i]);
        rj1[i] = std::abs(zj1[i]);
        theta[i] = std::arg(zj[i]);
        theta1[i] = std::arg(zj1[i]);
    }
    for (int i = 0; i < 4; ++i)
    {
        zp[i] = cos(angle_p[i]) + c[i] * sin(angle_p[i]);
        rp[i] = std::abs(zp[i]);
        thetap[i] = std::arg(zp[i]);
    }
    Eigen::Matrix<double, 8, 8, Eigen::RowMajor> M;
    M.setZero();
    auto calculation_M = [&]() -> void
    {
        for (int i = 0; i < 3; ++i)
        {
            M(2 * i, 2 * i) = std::pow(rj[i], 1 - alpha) * cos((1 - alpha) * theta[i]);
            M(2 * i, 2 * i + 1) = -std::pow(rj[i], 1 - alpha) * sin((1 - alpha) * theta[i]);
            M(2 * i, 2 * i + 2) = -std::pow(rj1[i], 1 - alpha) * cos((1 - alpha) * theta1[i]);
            M(2 * i, 2 * i + 3) = std::pow(rj1[i], 1 - alpha) * sin((1 - alpha) * theta1[i]);
            M(2 * i + 1, 2 * i) = std::pow(rj[i], -alpha) * ((-sin(angle[i]) * K[i](0, 0) + cos(angle[i]) * K[i](0, 1)) * cos(alpha * theta[i]) + (-sin(angle[i]) * K[i](1, 0) + cos(angle[i]) * K[i](1, 1)) * (cos(alpha * theta[i]) * c[i].real() + sin(alpha * theta[i]) * c[i].imag()));
            M(2 * i + 1, 2 * i + 1) = std::pow(rj[i], -alpha) * ((-sin(angle[i]) * K[i](0, 0) + cos(angle[i]) * K[i](0, 1)) * sin(alpha * theta[i]) + (-sin(angle[i]) * K[i](1, 0) + cos(angle[i]) * K[i](1, 1)) * (sin(alpha * theta[i]) * c[i].real() - cos(alpha * theta[i]) * c[i].imag()));
            M(2 * i + 1, 2 * i + 2) = -std::pow(rj1[i], -alpha) * ((-sin(angle[i]) * K[i + 1](0, 0) + cos(angle[i]) * K[i + 1](0, 1)) * cos(alpha * theta1[i]) + (-sin(angle[i]) * K[i + 1](1, 0) + cos(angle[i]) * K[i + 1](1, 1)) * (cos(alpha * theta1[i]) * c[i + 1].real() + sin(alpha * theta1[i]) * c[i + 1].imag()));
            M(2 * i + 1, 2 * i + 3) = -std::pow(rj1[i], -alpha) * ((-sin(angle[i]) * K[i + 1](0, 0) + cos(angle[i]) * K[i + 1](0, 1)) * sin(alpha * theta1[i]) + (-sin(angle[i]) * K[i + 1](1, 0) + cos(angle[i]) * K[i + 1](1, 1)) * (sin(alpha * theta1[i]) * c[i + 1].real() - cos(alpha * theta1[i]) * c[i + 1].imag()));
        }
        M(6, 6) = cos((1 - alpha) * Pi);
        M(6, 7) = -sin((1 - alpha) * Pi);
        M(6, 0) = -cos((1 - alpha) * Pi);
        M(6, 1) = -sin((1 - alpha) * Pi);
        M(7, 6) = -K[3](0, 1) * cos(alpha * Pi) - K[3](1, 1) * (cos(alpha * Pi) * c[3].real() + sin(alpha * Pi) * c[3].imag());
        M(7, 7) = -K[3](1, 0) * sin(alpha * Pi) - K[3](1, 1) * (sin(alpha * Pi) * c[3].real() - cos(alpha * Pi) * c[3].imag());
        M(7, 0) = K[0](0, 1) * cos(alpha * Pi) + K[0](1, 1) * (cos(alpha * Pi) * c[0].real() - sin(alpha * Pi) * c[0].imag());
        M(7, 1) = -K[0](1, 0) * sin(alpha * Pi) - K[0](1, 1) * (sin(alpha * Pi) * c[0].real() + cos(alpha * Pi) * c[0].imag());
    };
    const auto _nullspace = [](const double *mat, double *x, int m) -> void
    {
        const int n = m * m;
        assert(n > 0);

        double *bak = new double[n];
        dlacpy("N", &m, &m, mat, &m, bak, &m);

        const int t = m * 5;
        const int l = 1;

        double *s = new double[m];
        double *u = new double[n];
        double *v = new double[n];
        double *w = new double[t];

        const double r = dnrm2(&n, bak, &l);
        const double e = m * (std::nextafter(r, std::numeric_limits<double>::max()) - r);
        assert(e > 0.0);

        int err = 0;
        dgesvd("S", "S", &m, &m, bak, &m, s, u, &m, v, &m, w, &t, &err);
        dcopy(&m, v + m - 1, &m, x, &l);
        assert(err == 0);

        delete[] s;
        delete[] u;
        delete[] v;
        delete[] w;
    };
    double delta = 1.0e-4;
    double left = 0.0, right = 0.0, mid = 0.0;
    double lefta = 0.0, righta = 0.0, mida = 0.0;
    for (double x = 0.0; x <= 1.0 - delta; x += delta)
    {
        alpha = x;
        lefta = alpha;
        calculation_M();
        left = M.determinant();
        if (abs(left) < 1.0e-6)
            break;
        righta = alpha = x + delta;
        calculation_M();
        right = M.determinant();
        if (abs(right) < 1.0e-6)
            break;
        if (left * right < 0.0)
        {
            for (int iter = 0; iter < 60; iter++)
            {
                mida = alpha = (lefta + righta) / 2.0;
                calculation_M();
                mid = M.determinant();
                if (abs(mid) < 1.0e-6)
                    break;
                else if (mid * left < 0)
                {
                    righta = mida;
                    right = mid;
                }
                else
                {
                    lefta = mida;
                    left = mid;
                }
            }
            break;
        }
    }
    if (alpha > 0.0 && alpha <= 0.99)
    {
        M.transposeInPlace();
        double null[8], lamdap[4];
        _nullspace(M.data(), null, 8);
        M.transposeInPlace();
        for (int q = 0; q < 4; ++q)
            lamdap[q] = 2.0 * std::pow(r_p[q] * rp[q], 1.0 - alpha) * (null[2 * q] * cos((1.0 - alpha) * thetap[q]) - null[2 * q + 1] * sin((1.0 - alpha) * thetap[q]));
        auto C = [&](int n, int m) -> double &
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
        for (int q = 0; q < 4; ++q)
        { // q+1->q的流量
            double lamda = null[2 * q] * M(2 * q + 1, 2 * q) + null[2 * q + 1] * M(2 * q + 1, 2 * q + 1);
            lamda *= 2.0 * lz * std::pow(r[q], 1.0 - alpha) / (2.0 - alpha);
            lamda /= (lamdap[(q + 1) % 4] - lamdap[q]);
            C(idx[q], idx[q]) -= lamda;
            C(idx[q], idx[(q + 1) % 4]) += lamda;
            C(idx[(q + 1) % 4], idx[(q + 1) % 4]) -= lamda;
            C(idx[(q + 1) % 4], idx[q]) += lamda;
        }
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
        for (int j = 0; j < ny; j++)
            for (int i = 0; i < nx; i++)
            {
                int cur = k * nx * ny + j * nx + i;
                count = 0;
                if (k > 0)
                {
                    Idx[nnz] = cur - nx * ny;
                    nnz++;
                    count++;
                }
                if (j > 0)
                {
                    Idx[nnz] = cur - nx;
                    nnz++;
                    count++;
                }
                if (i > 0)
                {
                    Idx[nnz] = cur - 1;
                    nnz++;
                    count++;
                }
                Idx[nnz] = cur;
                nnz++;
                count++;
                if (i < nx - 1)
                {
                    Idx[nnz] = cur + 1;
                    nnz++;
                    count++;
                }
                if (j < ny - 1)
                {
                    Idx[nnz] = cur + nx;
                    nnz++;
                    count++;
                }
                if (k < nz - 1)
                {
                    Idx[nnz] = cur + nx * ny;
                    nnz++;
                    count++;
                }
                Ptr[n + 1] = Ptr[n] + count;
                ++n;
            }

    for (int k = 0; k <= nz; k++)
        for (int j = 0; j < ny; j++)
            for (int i = 0; i <= nx; i++)
                if (i == 0 || i == nx || k == 0 || k == nz || _MPFA_ONLY)
                {
                    // MPFA(i, j, k, Axis::YPOSITIVE);
                    // MPFA(i, j + 1, k, Axis::YNEGATIVE);
                    TPFA(i, j, k, 1);
                }
                else
                {
                    int idx[4];
                    idx[0] = (k - 1) * nx * ny + j * nx + i - 1;
                    idx[1] = (k - 1) * nx * ny + j * nx + i;
                    idx[2] = k * nx * ny + j * nx + i;
                    idx[3] = k * nx * ny + j * nx + i - 1;
                    // 柱坐标，l,r,angle
                    double new_o[3], new_x[3], new_y[3], new_z[3];
                    _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 2, 0), new_o);
                    _get_normalized(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 2, 0), new_y);
                    double pc[4][3];
                    _get_centroid(&verts1(k - 1, j, i, 0, 0), &verts1(k - 1, j, i, 2, 0), &verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 6, 0), pc[0]);
                    _get_centroid(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 1, 0), &verts1(k, j, i, 2, 0), &verts1(k, j, i, 3, 0), pc[1]);
                    _get_centroid(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 2, 0), &verts1(k, j, i, 4, 0), &verts1(k, j, i, 6, 0), pc[2]);
                    _get_centroid(&verts1(k, j, i - 1, 0, 0), &verts1(k, j, i - 1, 1, 0), &verts1(k, j, i - 1, 2, 0), &verts1(k, j, i - 1, 3, 0), pc[3]);
                    // 求l->求r
                    double l[4];
                    double r[4], r_p[4];
                    double temp[3];
                    for (int q = 0; q < 4; ++q)
                    {
                        _vectorize(new_o, pc[q], temp);
                        l[q] = _dot_product(temp, new_y);
                        r[q] = std::sqrt(_l2norm(temp) * _l2norm(temp) - l[q] * l[q]);
                    }
                    // 先求新坐标轴，再求angle
                    _vectorize(new_o, pc[3], new_x);
                    _cross_product(new_y, new_x, new_z);
                    _get_normalized(new_z);
                    _cross_product(new_y, new_z, new_x);
                    double angle[4], angle_p[4];
                    angle[3] = Pi;
                    for (int q = 0; q < 3; ++q)
                    {
                        _vectorize(new_o, pc[q], temp);
                        _get_angle(new_x, new_z, temp, r[q], angle[q]);
                    }
                    // 渗透率张量变换，变换矩阵为正交矩阵
                    Eigen::Matrix3d J, matK[4];
                    J << new_x[0], new_x[1], new_x[2],
                        new_y[0], new_y[1], new_y[2],
                        new_z[0], new_z[1], new_z[2];
                    Eigen::RowVector3d p_o;
                    p_o << _cx(k - 1, j, i - 1) - new_o[0], _cy(k - 1, j, i - 1) - new_o[1], _cz(k - 1, j, i - 1) - new_o[2];
                    p_o *= J.transpose();
                    r_p[0] = std::sqrt(p_o(0) * p_o(0) + p_o(2) * p_o(2));
                    _get_angle(p_o[0], p_o[2], r_p[0], angle_p[0]);
                    p_o << _cx(k - 1, j, i) - new_o[0], _cy(k - 1, j, i) - new_o[1], _cz(k - 1, j, i) - new_o[2];
                    p_o *= J.transpose();
                    r_p[1] = std::sqrt(p_o(0) * p_o(0) + p_o(2) * p_o(2));
                    _get_angle(p_o[0], p_o[2], r_p[1], angle_p[1]);
                    p_o << _cx(k, j, i) - new_o[0], _cy(k, j, i) - new_o[1], _cz(k, j, i) - new_o[2];
                    p_o *= J.transpose();
                    r_p[2] = std::sqrt(p_o(0) * p_o(0) + p_o(2) * p_o(2));
                    _get_angle(p_o[0], p_o[2], r_p[2], angle_p[2]);
                    p_o << _cx(k, j, i - 1) - new_o[0], _cy(k, j, i - 1) - new_o[1], _cz(k, j, i - 1) - new_o[2];
                    p_o *= J.transpose();
                    r_p[3] = std::sqrt(p_o(0) * p_o(0) + p_o(2) * p_o(2));
                    _get_angle(p_o[0], p_o[2], r_p[3], angle_p[3]);
                    for (int q = 0; q < 4; ++q)
                    {
                        _get_matK(pem1[idx[q]], matK[q]);
                        matK[q] = J * matK[q] * J.transpose();
                    }
                    double alpha = 0.0;
                    Eigen::Matrix2d K[4];
                    for (int q = 0; q < 4; ++q)
                    {
                        K[q](0, 0) = matK[q](0, 0);
                        K[q](0, 1) = matK[q](0, 2);
                        K[q](1, 0) = matK[q](2, 0);
                        K[q](1, 1) = matK[q](2, 2);
                    }

                    FAM(idx, _get_distance(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 2, 0)), l, r, angle, K, alpha, r_p, angle_p);
                    if (alpha == 0.0 || alpha > 0.99)
                    {
                        // MPFA(i, j, k, Axis::YPOSITIVE);
                        // MPFA(i, j + 1, k, Axis::YNEGATIVE);
                        TPFA(i, j, k, 1);
                    }
                }
    for (int k = 0; k <= nz; k++)
        for (int j = 0; j <= ny; j++)
            for (int i = 0; i < nx; i++)
                if (j == 0 || j == ny || k == 0 || k == nz || _MPFA_ONLY)
                {
                    // MPFA(i, j, k, Axis::XPOSITIVE);
                    // MPFA(i + 1, j, k, Axis::XNEGATIVE);
                    TPFA(i, j, k, 0);
                }
                else
                {
                    int idx[4];
                    idx[0] = (k - 1) * nx * ny + (j - 1) * nx + i;
                    idx[1] = (k - 1) * nx * ny + j * nx + i;
                    idx[2] = k * nx * ny + j * nx + i;
                    idx[3] = k * nx * ny + (j - 1) * nx + i;
                    // 柱坐标，l,r,angle
                    double new_o[3], new_x[3], new_y[3], new_z[3];
                    _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 1, 0), new_o);
                    _get_normalized(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 1, 0), new_x);
                    double pc[4][3];
                    _get_centroid(&verts1(k - 1, j, i, 0, 0), &verts1(k - 1, j, i, 1, 0), &verts1(k - 1, j, i, 4, 0), &verts1(k - 1, j, i, 5, 0), pc[0]);
                    _get_centroid(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 1, 0), &verts1(k, j, i, 2, 0), &verts1(k, j, i, 3, 0), pc[1]);
                    _get_centroid(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 1, 0), &verts1(k, j, i, 4, 0), &verts1(k, j, i, 5, 0), pc[2]);
                    _get_centroid(&verts1(k, j - 1, i, 0, 0), &verts1(k, j - 1, i, 1, 0), &verts1(k, j - 1, i, 2, 0), &verts1(k, j - 1, i, 3, 0), pc[3]);
                    // 求l->求r
                    double l[4];
                    double r[4], r_p[4];
                    double temp[3];
                    for (int q = 0; q < 4; ++q)
                    {
                        _vectorize(new_o, pc[q], temp);
                        l[q] = _dot_product(temp, new_x);
                        r[q] = std::sqrt(_l2norm(temp) * _l2norm(temp) - l[q] * l[q]);
                    }
                    // 先求新坐标轴，再求angle
                    _vectorize(new_o, pc[3], new_y);
                    _cross_product(new_y, new_x, new_z);
                    _get_normalized(new_z);
                    _cross_product(new_z, new_x, new_y);
                    double angle[4], angle_p[4];
                    angle[3] = Pi;
                    for (int q = 0; q < 3; ++q)
                    {
                        _vectorize(new_o, pc[q], temp);
                        _get_angle(new_y, new_z, temp, r[q], angle[q]);
                    }
                    // 渗透率张量变换，变换矩阵为正交矩阵
                    Eigen::Matrix3d J, matK[4];
                    J << new_x[0], new_x[1], new_x[2],
                        new_y[0], new_y[1], new_y[2],
                        new_z[0], new_z[1], new_z[2];
                    Eigen::RowVector3d p_o;
                    p_o << _cx(k - 1, j - 1, i) - new_o[0], _cy(k - 1, j - 1, i) - new_o[1], _cz(k - 1, j - 1, i) - new_o[2];
                    p_o *= J.transpose();
                    r_p[0] = std::sqrt(p_o(1) * p_o(1) + p_o(2) * p_o(2));
                    _get_angle(p_o[1], p_o[2], r_p[0], angle_p[0]);
                    p_o << _cx(k - 1, j, i) - new_o[0], _cy(k - 1, j, i) - new_o[1], _cz(k - 1, j, i) - new_o[2];
                    p_o *= J.transpose();
                    r_p[1] = std::sqrt(p_o(1) * p_o(1) + p_o(2) * p_o(2));
                    _get_angle(p_o[1], p_o[2], r_p[1], angle_p[1]);
                    p_o << _cx(k, j, i) - new_o[0], _cy(k, j, i) - new_o[1], _cz(k, j, i) - new_o[2];
                    p_o *= J.transpose();
                    r_p[2] = std::sqrt(p_o(1) * p_o(1) + p_o(2) * p_o(2));
                    _get_angle(p_o[1], p_o[2], r_p[2], angle_p[2]);
                    p_o << _cx(k, j - 1, i) - new_o[0], _cy(k, j - 1, i) - new_o[1], _cz(k, j - 1, i) - new_o[2];
                    p_o *= J.transpose();
                    r_p[3] = std::sqrt(p_o(1) * p_o(1) + p_o(2) * p_o(2));
                    _get_angle(p_o[1], p_o[2], r_p[3], angle_p[3]);
                    for (int q = 0; q < 4; ++q)
                    {
                        _get_matK(pem1[idx[q]], matK[q]);
                        matK[q] = J * matK[q] * J.transpose();
                    }
                    double alpha = 0.0;
                    Eigen::Matrix2d K[4];
                    for (int q = 0; q < 4; ++q)
                        K[q] = matK[q].block(1, 1, 2, 2);
                    FAM(idx, _get_distance(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 1, 0)), l, r, angle, K, alpha, r_p, angle_p);
                    if (alpha == 0.0 || alpha > 0.99)
                    {
                        // MPFA(i, j, k, Axis::XPOSITIVE);
                        // MPFA(i + 1, j, k, Axis::XNEGATIVE);
                        TPFA(i, j, k, 0);
                    }
                }
    for (int k = 0; k < nz; k++)
        for (int j = 0; j <= ny; j++)
            for (int i = 0; i <= nx; i++)
                if (j == 0 || j == ny || i == 0 || i == nx || _MPFA_ONLY)
                {
                    // MPFA(i, j, k, Axis::ZPOSITIVE);
                    // MPFA(i, j, k + 1, Axis::ZNEGATIVE);
                    TPFA(i, j, k, 2);
                }
                else
                {
                    int idx[4];
                    idx[0] = k * nx * ny + (j - 1) * nx + i - 1;
                    idx[1] = k * nx * ny + (j - 1) * nx + i;
                    idx[2] = k * nx * ny + j * nx + i;
                    idx[3] = k * nx * ny + j * nx + i - 1;
                    // 柱坐标，l,r,angle
                    double new_o[3], new_x[3], new_y[3], new_z[3];
                    _get_midpoint(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 4, 0), new_o);
                    _get_normalized(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 4, 0), new_z);
                    double pc[4][3];
                    _get_centroid(&verts1(k, j - 1, i, 0, 0), &verts1(k, j - 1, i, 2, 0), &verts1(k, j - 1, i, 4, 0), &verts1(k, j - 1, i, 6, 0), pc[0]);
                    _get_centroid(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 1, 0), &verts1(k, j, i, 4, 0), &verts1(k, j, i, 5, 0), pc[1]);
                    _get_centroid(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 2, 0), &verts1(k, j, i, 4, 0), &verts1(k, j, i, 6, 0), pc[2]);
                    _get_centroid(&verts1(k, j, i - 1, 0, 0), &verts1(k, j, i - 1, 1, 0), &verts1(k, j, i - 1, 4, 0), &verts1(k, j, i - 1, 5, 0), pc[3]);
                    // 求l->求r
                    double l[4];
                    double r[4], r_p[4];
                    double temp[3];
                    for (int q = 0; q < 4; ++q)
                    {
                        _vectorize(new_o, pc[q], temp);
                        l[q] = _dot_product(temp, new_z);
                        r[q] = std::sqrt(_l2norm(temp) * _l2norm(temp) - l[q] * l[q]);
                    }
                    // 先求新坐标轴，再求angle
                    _vectorize(new_o, pc[3], new_x);
                    _cross_product(new_x, new_z, new_y);
                    _get_normalized(new_y);
                    _cross_product(new_y, new_z, new_x);
                    double angle[4], angle_p[4];
                    angle[3] = Pi;
                    for (int q = 0; q < 3; ++q)
                    {
                        _vectorize(new_o, pc[q], temp);
                        _get_angle(new_x, new_y, temp, r[q], angle[q]);
                    }
                    // 渗透率张量变换，变换矩阵为正交矩阵
                    Eigen::Matrix3d J, matK[4];
                    J << new_x[0], new_x[1], new_x[2],
                        new_y[0], new_y[1], new_y[2],
                        new_z[0], new_z[1], new_z[2];
                    Eigen::RowVector3d p_o;
                    p_o << _cx(k, j - 1, i - 1) - new_o[0], _cy(k, j - 1, i - 1) - new_o[1], _cz(k, j - 1, i - 1) - new_o[2];
                    p_o *= J.transpose();
                    r_p[0] = std::sqrt(p_o(0) * p_o(0) + p_o(1) * p_o(1));
                    _get_angle(p_o[0], p_o[1], r_p[0], angle_p[0]);
                    p_o << _cx(k, j - 1, i) - new_o[0], _cy(k, j - 1, i) - new_o[1], _cz(k, j - 1, i) - new_o[2];
                    p_o *= J.transpose();
                    r_p[1] = std::sqrt(p_o(0) * p_o(0) + p_o(1) * p_o(1));
                    _get_angle(p_o[0], p_o[1], r_p[1], angle_p[1]);
                    p_o << _cx(k, j, i) - new_o[0], _cy(k, j, i) - new_o[1], _cz(k, j, i) - new_o[2];
                    p_o *= J.transpose();
                    r_p[2] = std::sqrt(p_o(0) * p_o(0) + p_o(1) * p_o(1));
                    _get_angle(p_o[0], p_o[1], r_p[2], angle_p[2]);
                    p_o << _cx(k, j, i - 1) - new_o[0], _cy(k, j, i - 1) - new_o[1], _cz(k, j, i - 1) - new_o[2];
                    p_o *= J.transpose();
                    r_p[3] = std::sqrt(p_o(0) * p_o(0) + p_o(1) * p_o(1));
                    _get_angle(p_o[0], p_o[1], r_p[3], angle_p[3]);
                    for (int q = 0; q < 4; ++q)
                    {
                        _get_matK(pem1[idx[q]], matK[q]);
                        matK[q] = J * matK[q] * J.transpose();
                    }
                    double alpha = 0.0;
                    Eigen::Matrix2d K[4];
                    for (int q = 0; q < 4; ++q)
                        K[q] = matK[q].block(0, 0, 2, 2);
                    FAM(idx, _get_distance(&verts1(k, j, i, 0, 0), &verts1(k, j, i, 4, 0)), l, r, angle, K, alpha, r_p, angle_p);
                    if (alpha == 0.0 || alpha > 0.99)
                    {
                        // MPFA(i, j, k, Axis::ZPOSITIVE);
                        // MPFA(i, j, k + 1, Axis::ZNEGATIVE);
                        TPFA(i, j, k, 2);
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
    delete[] pem1, p, Ptr, Idx, Val, B, _verts1, cx, cy, cz;
    std::cout << std::endl;
}