#pragma once

#include <cmath>
#include <math.h>

#include <iomanip>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <mdspan>
#include <numbers>
#include <numeric>
#include <set>
#include <vector>
#include <Eigen/Dense>
#include <mkl.h>

namespace std
{
    using std::experimental::extents;
    using std::experimental::mdspan;
}

constexpr double Pi = 3.141592653589793;
constexpr double Plow = 0.0, Phigh = 1.0;
constexpr double epsilonP = 1.0e-6, residual = 0.0001;
// constexpr double pem0 = 300.0;
// constexpr double tolpem = pem0 + 1.0;
constexpr int nxx = 3, nyy = 3, nzz = 3; // ϸ��ǰ
// constexpr int nxx = 10, nyy = 9, nzz = 8;                          //ϸ��ǰ
int nx, ny, nz;
struct Ktensor;
Ktensor *pem1;
double *_verts1;
int *Ptr;
int *Idx;
double *Val;
double *B;
double *cx;
double *cy;
double *cz;
std::mdspan<double, std::extents<size_t, -1, -1, -1, 8, 3>> verts1;
std::mdspan<double, std::extents<size_t, -1, -1, -1>> _cx;
std::mdspan<double, std::extents<size_t, -1, -1, -1>> _cy;
std::mdspan<double, std::extents<size_t, -1, -1, -1>> _cz;
enum class Axis // Determine the positive or negative of a certain coordinate value of the
                // direction vector
{
    XPOSITIVE,
    YPOSITIVE,
    ZPOSITIVE,
    XNEGATIVE,
    YNEGATIVE,
    ZNEGATIVE
};

enum class Line
{
    LINE31,
    LINE23,
    LINE15,
    LINE26,
    LINE37,
    LINE45,
    LINE75,
    LINE67,
    LINE46
};

struct Ktensor
{
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    double xy = 0.0;
    double yz = 0.0;
    double xz = 0.0;

    Ktensor &operator=(const Ktensor &v)
    {
        if (this != &v)
        {
            x = v.x;
            y = v.y;
            z = v.z;
            xy = v.xy;
            yz = v.yz;
            xz = v.xz;
        }
        return *this;
    }
};

bool _approx(double x, double tol = std::sqrt(std::numeric_limits<double>::epsilon()))
{
    return std::abs(x) < tol;
}

double _l2norm(const double *v)
{
    return std::hypot(v[0], v[1], v[2]);
}

void _vectorize(const double *p1, const double *p2, double *v)
{
    std::transform(p2, p2 + 3, p1, v, std::minus<>{});
}

void _cross_product(const double *v1, const double *v2, double *n)
{
    n[0] = v1[1] * v2[2] - v1[2] * v2[1];
    n[1] = v1[2] * v2[0] - v1[0] * v2[2];
    n[2] = v1[0] * v2[1] - v1[1] * v2[0];
}

double _dot_product(const double *v1, const double *v2)
{
    return std::inner_product(v1, v1 + 3, v2, 0.0);
}

double _get_distance(const double *c, const double *p1, const double *p2, const double *p3)
{
    double v1[3], v2[3], n[3];
    _vectorize(p1, p2, v1);
    _vectorize(p1, p3, v2);
    _cross_product(v1, v2, n);
    return std::abs(_dot_product(c, n) - _dot_product(p1, n)) / _l2norm(n);
}

void _get_midpoint(const double *p1, const double *p2, double *m)
{
    for (int i = 0; i < 3; i++)
    {
        m[i] = (p1[i] + p2[i]) / 2.0;
    }
}

void _get_centroid(const double *p1, const double *p2, const double *p3, const double *p4,
                   const double *p5, const double *p6, const double *p7, const double *p8,
                   double *c)
{
    c[0] = 0.0;
    c[1] = 0.0;
    c[2] = 0.0;

    for (int i = 0; i < 3; ++i)
    {
        c[i] += p1[i];
        c[i] += p2[i];
        c[i] += p3[i];
        c[i] += p4[i];
        c[i] += p5[i];
        c[i] += p6[i];
        c[i] += p7[i];
        c[i] += p8[i];
    }

    c[0] /= 8;
    c[1] /= 8;
    c[2] /= 8;
}

double _harmonic(double x, double y)
{
    return 2.0 * x * y / (x + y);
}

void _normalize(double *v)
{
    using std::placeholders::_1;
    std::transform(v, v + 3, v, std::bind(std::divides<>{}, _1, _l2norm(v)));
}

double _get_distance(const double *p1, const double *p2)
{
    return std::hypot(p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]);
}

double _get_distance(const double *c, const double *p1, const double *p2)
{
    double v1[3];
    _vectorize(p1, p2, v1);
    _normalize(v1);

    double v2[3];
    _vectorize(p1, c, v2);

    double l = _l2norm(v2);
    return std::sin(std::acos(_dot_product(v1, v2) / l)) * l;
}

void _multcross(const double *x, const double *y, double *z) // ���
{
    z[0] = x[1] * y[2] - x[2] * y[1];
    z[1] = x[2] * y[0] - x[0] * y[2];
    z[2] = x[0] * y[1] - x[1] * y[0];
}

double _get_area(const double *p1, const double *p2, const double *p3)
{
    double v1[3], v2[3], n[3];
    for (int i = 0; i < 3; ++i)
    {
        v1[i] = p2[i] - p1[i];
        v2[i] = p3[i] - p1[i];
    }
    _multcross(v1, v2, n);
    return std::hypot(n[0], n[1], n[2]) / 2;
}

double _get_area(const double *p1, const double *p2, const double *p3, const double *p4)
{
    return _get_area(p1, p2, p3) + _get_area(p2, p3, p4);
}

void _get_normalized(const double *p1, const double *p2, double *v)
{
    double d = _get_distance(p1, p2);
    for (int i = 0; i < 3; ++i)
        v[i] = (p2[i] - p1[i]) / d;
}

void _get_normalized(double *v)
{
    double d = std::hypot(v[0], v[1], v[2]);
    for (int i = 0; i < 3; ++i)
        v[i] = v[i] / d;
}

void _get_plane(const double *p1, const double *p2, double *o, double *n)
{
    _get_midpoint(p1, p2, o);
    _get_normalized(p1, p2, n);
}

void _get_reverse(double *v)
{
    for (int i = 0; i < 3; ++i)
        v[i] = -v[i];
}

void _get_centroid(const double *p1, const double *p2, double *c)
{
    c[0] = 0.0;
    c[1] = 0.0;
    c[2] = 0.0;

    for (int i = 0; i < 3; ++i)
    {
        c[i] += p1[i];
        c[i] += p2[i];
    }

    c[0] /= 2;
    c[1] /= 2;
    c[2] /= 2;
}

void _get_centroid(const double *p1, const double *p2, const double *p3, const double *p4,
                   double *c)
{
    c[0] = 0.0;
    c[1] = 0.0;
    c[2] = 0.0;

    for (int i = 0; i < 3; ++i)
    {
        c[i] += p1[i];
        c[i] += p2[i];
        c[i] += p3[i];
        c[i] += p4[i];
    }

    c[0] /= 4;
    c[1] /= 4;
    c[2] /= 4;
}

double _get_plane_angle(const double *p1, const double *p2, const double *s1, const double *s2)
{
    double v1[3], v2[3], v3[3], v4[3];
    _vectorize(s1, p1, v1);
    _vectorize(s1, p2, v2);
    _vectorize(s2, p1, v3);
    _vectorize(s2, p2, v4);

    double n1[3], n2[3];
    _cross_product(v1, v2, n1);
    _cross_product(v3, v4, n2);

    double a = 1.0;
    a *= std::hypot(n1[0], n1[1], n1[2]);
    a *= std::hypot(n2[0], n2[1], n2[2]);
    return std::acos(_dot_product(n1, n2) / a);
}

void _transform(Ktensor *K, int n)
{
    for (int i = 0; i < n; i++)
    {
        K[i].x *= 1.0;
        K[i].y *= 1.0;
        K[i].z *= 0.1;
        K[i].xy = 0.0;
        K[i].yz = 0.0;
        K[i].xz = 0.0;
    }
}

double _assit_divied(const double *a, int i, int j, int k, int divn)
{
    double id = 1.0 * i, jd = 1.0 * j, kd = 1.0 * k, divnd = 1.0 * divn; // ת����double
    double verts1 = a[0] + (id / divnd) * (a[1] - a[0]);
    double verts2 = a[2] + (id / divnd) * (a[3] - a[2]);
    double verts12 = verts1 + (jd / divnd) * (verts2 - verts1);
    double verts3 = a[4] + (id / divnd) * (a[5] - a[4]);
    double verts4 = a[6] + (id / divnd) * (a[7] - a[6]);
    double verts34 = verts3 + (jd / divnd) * (verts4 - verts3);
    return (verts12 + (kd / divn) * (verts34 - verts12));
}

void _get_divide(const double *verts0, double *verts, const Ktensor *perm0, Ktensor *perm, int divn)
{
    std::mdspan<const double, std::extents<size_t, -1, -1, -1, 8, 3>> _verts0(verts0, nzz, nyy,
                                                                              nxx);
    std::mdspan<double, std::extents<size_t, -1, -1, -1, 8, 3>> _verts(verts, nzz * divn, nyy * divn, nxx * divn);
    for (int k = 0; k < nzz; k++)
    {
        for (int j = 0; j < nyy; j++)
        {
            for (int i = 0; i < nxx; i++)
            {
                int n0 = k * nxx * nyy + j * nxx + i;
                // ϸ���µ�ǰ������Ϣ
                int k0 = divn * k;
                int j0 = divn * j;
                int i0 = divn * i;

                for (int kk = 0; kk < divn; kk++)
                {
                    for (int jj = 0; jj < divn; jj++)
                    {
                        for (int ii = 0; ii < divn; ii++)
                        {
                            int n = (k0 + kk) * nxx * divn * nyy * divn + (j0 + jj) * nxx * divn + i0 + ii;
                            perm[n] = perm0[n0];
                        }
                    }
                }

                for (int xyz = 0; xyz < 3; xyz++)
                {
                    double a[8];
                    for (int count = 0; count < 8; count++)
                    {
                        a[count] = _verts0(k, j, i, count, xyz);
                        // std::cout << " " << a[count] << " " ;
                    }

                    for (int kk = 0; kk < divn; kk++)
                    {
                        for (int jj = 0; jj < divn; jj++)
                        {
                            for (int ii = 0; ii < divn; ii++)
                            {

                                _verts(k0 + kk, j0 + jj, i0 + ii, 0, xyz) =
                                    _assit_divied(a, ii, jj, kk, divn);
                                _verts(k0 + kk, j0 + jj, i0 + ii, 1, xyz) =
                                    _assit_divied(a, ii + 1, jj, kk, divn);
                                _verts(k0 + kk, j0 + jj, i0 + ii, 2, xyz) =
                                    _assit_divied(a, ii, jj + 1, kk, divn);
                                _verts(k0 + kk, j0 + jj, i0 + ii, 3, xyz) =
                                    _assit_divied(a, ii + 1, jj + 1, kk, divn);
                                _verts(k0 + kk, j0 + jj, i0 + ii, 4, xyz) =
                                    _assit_divied(a, ii, jj, kk + 1, divn);
                                _verts(k0 + kk, j0 + jj, i0 + ii, 5, xyz) =
                                    _assit_divied(a, ii + 1, jj, kk + 1, divn);
                                _verts(k0 + kk, j0 + jj, i0 + ii, 6, xyz) =
                                    _assit_divied(a, ii, jj + 1, kk + 1, divn);
                                _verts(k0 + kk, j0 + jj, i0 + ii, 7, xyz) =
                                    _assit_divied(a, ii + 1, jj + 1, kk + 1, divn);
                            }
                        }
                    }
                }
            }
        }
    }
}

// 3x3x3 ���̸�
void _get_3x3(double *_verts, Ktensor *pem, double pemx,
              double pemy, double pemz, double pemxy, double pemyz, double pemxz, int Lamda)
{
    std::mdspan<double, std::extents<size_t, -1, -1, -1, 8, 3>> verts(_verts, nzz, nyy, nxx);
    // for (int i = 0; i < nxx; i++)
    // {
    //     verts(0, 0, i, 0, 0) = i * 1.0;
    //     verts(0, 0, i, 0, 1) = i * 0.15;
    //     verts(0, 0, i, 0, 2) = 0.0;
    //     verts(0, 0, i, 1, 0) = i * 1.0 + 1.0;
    //     verts(0, 0, i, 1, 1) = i * 0.15 + 0.15;
    //     verts(0, 0, i, 1, 2) = 0.0;
    //     verts(0, 0, i, 2, 0) = i * 1.0 - 0.1;
    //     verts(0, 0, i, 2, 1) = i * 0.2 + 1.0;
    //     verts(0, 0, i, 2, 2) = 0.0;
    //     verts(0, 0, i, 3, 0) = i * 1.0 + 1.0 - 0.1;
    //     verts(0, 0, i, 3, 1) = i * 0.2 + 1.2;
    //     verts(0, 0, i, 3, 2) = 0.0;

    //     verts(0, 1, i, 0, 0) = verts(0, 0, i, 2, 0);
    //     verts(0, 1, i, 0, 1) = verts(0, 0, i, 2, 1);
    //     verts(0, 1, i, 0, 2) = verts(0, 0, i, 2, 2);
    //     verts(0, 1, i, 1, 0) = verts(0, 0, i, 3, 0);
    //     verts(0, 1, i, 1, 1) = verts(0, 0, i, 3, 1);
    //     verts(0, 1, i, 1, 2) = verts(0, 0, i, 3, 2);
    //     verts(0, 1, i, 2, 0) = i * 1.0;
    //     verts(0, 1, i, 2, 1) = i * 0.2 + 2.25;
    //     verts(0, 1, i, 2, 2) = 0.0;
    //     verts(0, 1, i, 3, 0) = i * 1.0 + 1.0;
    //     verts(0, 1, i, 3, 1) = i * 0.2 + 2.45;
    //     verts(0, 1, i, 3, 2) = 0.0;

    //     verts(0, 2, i, 0, 0) = verts(0, 1, i, 2, 0);
    //     verts(0, 2, i, 0, 1) = verts(0, 1, i, 2, 1);
    //     verts(0, 2, i, 0, 2) = verts(0, 1, i, 2, 2);
    //     verts(0, 2, i, 1, 0) = verts(0, 1, i, 3, 0);
    //     verts(0, 2, i, 1, 1) = verts(0, 1, i, 3, 1);
    //     verts(0, 2, i, 1, 2) = verts(0, 1, i, 3, 2);
    //     verts(0, 2, i, 2, 0) = i * 1.0 - 0.1;
    //     verts(0, 2, i, 2, 1) = i * 0.2 + 3.25;
    //     verts(0, 2, i, 2, 2) = 0.0;
    //     verts(0, 2, i, 3, 0) = i * 1.0 + 1.0 - 0.1;
    //     verts(0, 2, i, 3, 1) = i * 0.2 + 3.45;
    //     verts(0, 2, i, 3, 2) = 0.0;
    // }
    // for (int i = 0; i < nxx; i++)
    // {
    //     for (int j = 0; j < nyy; j++)
    //     {
    //         for (int l = 0; l < 4; l++)
    //         {
    //             verts(0, j, i, l + 4, 0) = verts(0, j, i, l, 0) + 0.1;
    //             verts(0, j, i, l + 4, 1) = verts(0, j, i, l, 1) + 0.1;
    //             verts(0, j, i, l + 4, 2) = verts(0, j, i, l, 2) + 0.9;
    //         }
    //     }
    // }
    // for (int i = 0; i < nxx; i++)
    // {
    //     for (int j = 0; j < nyy; j++)
    //     {
    //         for (int l = 0; l < 4; l++)
    //         {
    //             verts(1, j, i, l, 0) = verts(0, j, i, l + 4, 0);
    //             verts(1, j, i, l, 1) = verts(0, j, i, l + 4, 1);
    //             verts(1, j, i, l, 2) = verts(0, j, i, l + 4, 2);
    //             verts(1, j, i, l + 4, 0) = verts(1, j, i, l, 0) - 0.1;
    //             verts(1, j, i, l + 4, 1) = verts(1, j, i, l, 1) - 0.1;
    //         }
    //     }
    // }
    // for (int i = 0; i < nxx; i++)
    // {
    //     for (int j = 0; j < nyy; j++)
    //     {
    //         verts(1, j, i, 4, 2) = verts(1, j, i, 0, 2) + 1.0;
    //         verts(1, j, i, 5, 2) = verts(1, j, i, 1, 2) + 1.0;
    //         verts(1, j, i, 6, 2) = verts(1, j, i, 2, 2) + 1.1;
    //         verts(1, j, i, 7, 2) = verts(1, j, i, 3, 2) + 1.1;
    //     }
    // }
    // for (int i = 0; i < nxx; i++)
    // {
    //     for (int j = 0; j < nyy; j++)
    //     {
    //         for (int l = 0; l < 4; l++)
    //         {
    //             verts(2, j, i, l, 0) = verts(1, j, i, l + 4, 0);
    //             verts(2, j, i, l, 1) = verts(1, j, i, l + 4, 1);
    //             verts(2, j, i, l, 2) = verts(1, j, i, l + 4, 2);
    //             verts(2, j, i, l + 4, 0) = verts(2, j, i, l, 0) - 0.1;
    //             verts(2, j, i, l + 4, 1) = verts(2, j, i, l, 1) - 0.1;
    //         }
    //     }
    // }
    // for (int i = 0; i < nxx; i++)
    // {
    //     for (int j = 0; j < nyy; j++)
    //     {
    //         verts(2, j, i, 4, 2) = verts(2, j, i, 0, 2) + 0.9;
    //         verts(2, j, i, 5, 2) = verts(2, j, i, 1, 2) + 0.9;
    //         verts(2, j, i, 6, 2) = verts(2, j, i, 2, 2) + 0.8;
    //         verts(2, j, i, 7, 2) = verts(2, j, i, 3, 2) + 0.8;
    //     }
    // }

    double DX = 1.0;
    double DY = 1.0;
    double DZ = 1.0;
    for (int k = 0; k < nzz; ++k) // �ѿ�������
    {
        for (int j = 0; j < nyy; ++j)
        {
            for (int i = 0; i < nxx; ++i)
            {
                verts(k, j, i, 0, 0) = i * DX;
                verts(k, j, i, 0, 1) = j * DY;
                verts(k, j, i, 0, 2) = k * DZ;
                verts(k, j, i, 1, 0) = i * DX + DX;
                verts(k, j, i, 1, 1) = j * DY;
                verts(k, j, i, 1, 2) = k * DZ;
                verts(k, j, i, 2, 0) = i * DX;
                verts(k, j, i, 2, 1) = j * DY + DY;
                verts(k, j, i, 2, 2) = k * DZ;
                verts(k, j, i, 3, 0) = i * DX + DX;
                verts(k, j, i, 3, 1) = j * DY + DY;
                verts(k, j, i, 3, 2) = k * DZ;
                verts(k, j, i, 4, 0) = i * DX;
                verts(k, j, i, 4, 1) = j * DY;
                verts(k, j, i, 4, 2) = k * DZ + DZ;
                verts(k, j, i, 5, 0) = i * DX + DX;
                verts(k, j, i, 5, 1) = j * DY;
                verts(k, j, i, 5, 2) = k * DZ + DZ;
                verts(k, j, i, 6, 0) = i * DX;
                verts(k, j, i, 6, 1) = j * DY + DY;
                verts(k, j, i, 6, 2) = k * DZ + DZ;
                verts(k, j, i, 7, 0) = i * DX + DX;
                verts(k, j, i, 7, 1) = j * DY + DY;
                verts(k, j, i, 7, 2) = k * DZ + DZ;
            }
        }
    }

    for (int k = 0; k < nzz; k++)
    {
        for (int j = 0; j < nyy; j++)
        {
            for (int i = 0; i < nxx; i++)
            {
                int idx = i + j + k;
                int i0 = 9 * k + 3 * j + i;
                if (idx % 2 == 0)
                {
                    pem[i0].x = pemx;
                    pem[i0].y = pemy;
                    pem[i0].z = pemz;
                    pem[i0].xy = pemxy;
                    pem[i0].yz = pemyz;
                    pem[i0].xz = pemxz;
                }
                else
                {
                    pem[i0].x = pemx * Lamda;
                    pem[i0].y = pemy * Lamda;
                    pem[i0].z = pemz * Lamda;
                    pem[i0].xy = pemxy * Lamda;
                    pem[i0].yz = pemyz * Lamda;
                    pem[i0].xz = pemxz * Lamda;
                }
            }
        }
    }
}

void _get_pcood(int nx, int ny, int nz, int i, int j, int k, const double *verts, double *p)
{
    std::mdspan<const double, std::extents<size_t, -1, -1, -1, 8, 3>> _verts(verts, nz, ny, nx);
    const double *pp[8];
    pp[0] = &_verts(k, j, i, 0, 0);
    pp[1] = &_verts(k, j, i, 1, 0);
    pp[2] = &_verts(k, j, i, 2, 0);
    pp[3] = &_verts(k, j, i, 3, 0);
    pp[4] = &_verts(k, j, i, 4, 0);
    pp[5] = &_verts(k, j, i, 5, 0);
    pp[6] = &_verts(k, j, i, 6, 0);
    pp[7] = &_verts(k, j, i, 7, 0);
    double sum[3] = {0.0};
    for (int i = 0; i < 8; i++)
    {
        sum[0] += pp[i][0];
        sum[1] += pp[i][1];
        sum[2] += pp[i][2];
    }
    for (int i = 0; i < 3; ++i)
        p[i] = sum[i] / 8;
}
void _Outputvtk(const double *verts, const double *perm)
{
    constexpr int nx = nxx;
    constexpr int ny = nyy;
    constexpr int nz = nzz;
    std::mdspan<const double, std::extents<size_t, -1, -1, -1, 8, 3>> _verts(verts, nz, ny, nx);

    int SIZE_X = nx;
    int SIZE_Y = ny;
    int SIZE_Z = nz;
    int SIZE_VERTICES = 8;
    int SIZE_DIMENSIONS = 3;

    std::ofstream outputFile("model_data.vtk");

    if (!outputFile)
    {
        std::cout << "Failed to create output file!" << std::endl;
        return;
    }

    // Write header information
    outputFile << "# vtk DataFile Version 3.0" << std::endl;
    outputFile << "3D Model Data" << std::endl;
    outputFile << "ASCII" << std::endl;
    outputFile << "DATASET UNSTRUCTURED_GRID" << std::endl;

    // Calculate total number of cells and vertices
    int numCells = SIZE_X * SIZE_Y * SIZE_Z;
    int numVertices = SIZE_X * SIZE_Y * SIZE_Z * SIZE_VERTICES;

    // Write vertices coordinates
    outputFile << "POINTS " << numVertices << " double" << std::endl;
    for (int k = 0; k < SIZE_Z; k++)
    {
        for (int j = 0; j < SIZE_Y; j++)
        {
            for (int i = 0; i < SIZE_X; i++)
            {
                outputFile << _verts(k, j, i, 0, 0) << " " << _verts(k, j, i, 0, 1) << " " << _verts(k, j, i, 0, 2) << std::endl;
                outputFile << _verts(k, j, i, 1, 0) << " " << _verts(k, j, i, 1, 1) << " " << _verts(k, j, i, 1, 2) << std::endl;
                outputFile << _verts(k, j, i, 3, 0) << " " << _verts(k, j, i, 3, 1) << " " << _verts(k, j, i, 3, 2) << std::endl;
                outputFile << _verts(k, j, i, 2, 0) << " " << _verts(k, j, i, 2, 1) << " " << _verts(k, j, i, 2, 2) << std::endl;

                // Top face vertices
                outputFile << _verts(k, j, i, 4, 0) << " " << _verts(k, j, i, 4, 1) << " " << _verts(k, j, i, 4, 2) << std::endl;
                outputFile << _verts(k, j, i, 5, 0) << " " << _verts(k, j, i, 5, 1) << " " << _verts(k, j, i, 5, 2) << std::endl;
                outputFile << _verts(k, j, i, 7, 0) << " " << _verts(k, j, i, 7, 1) << " " << _verts(k, j, i, 7, 2) << std::endl;
                outputFile << _verts(k, j, i, 6, 0) << " " << _verts(k, j, i, 6, 1) << " " << _verts(k, j, i, 6, 2) << std::endl;
            }
        }
    }

    // Write cell connectivity
    outputFile << "CELLS " << numCells << " " << numCells * (SIZE_VERTICES + 1) << std::endl;
    for (int i = 0; i < numCells; i++)
    {
        int offset = i * SIZE_VERTICES;
        outputFile << SIZE_VERTICES << " ";
        for (int v = 0; v < SIZE_VERTICES; v++)
        {
            int vertexIndex = offset + v;
            outputFile << vertexIndex << " ";
        }
        outputFile << std::endl;
    }

    // Write cell types
    outputFile << "CELL_TYPES " << numCells << std::endl;
    for (int i = 0; i < numCells; i++)
    {
        outputFile << "12" << std::endl; // VTK_HEXAHEDRON type
    }

    // Write parameter values
    outputFile << "CELL_DATA " << numCells << std::endl;
    outputFile << "SCALARS perm int 1" << std::endl;
    outputFile << "LOOKUP_TABLE default" << std::endl;
    for (int k = 0; k < SIZE_Z; k++)
    {
        for (int j = 0; j < SIZE_Y; j++)
        {
            for (int i = 0; i < SIZE_X; i++)
            {
                int n = k * SIZE_X * SIZE_Y + j * SIZE_X + i;
                outputFile << perm[n] << std::endl;
            }
        }
    }

    outputFile.close();
}
// 求取4个点构成曲面的法向量 1->2->3->4构成环
void _get_surface_normal(const double *p1, const double *p2, const double *p3, const double *p4, double *n, const enum Axis axi)
{
    double v1[3], v2[3], v3[3];
    _vectorize(p1, p2, v1);
    _vectorize(p1, p3, v2);
    _vectorize(p1, p4, v3);

    double n1[3], n2[3];
    _cross_product(v1, v2, n1);
    _cross_product(v2, v3, n2);

    for (int i = 0; i < 3; ++i)
    {
        n[i] = (n1[i] + n2[i]) / 2.0;
    }
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
}
void _get_triangle_normal(const double *p1, const double *p2, const double *p3, double *n, const enum Axis axi)
{
    double v1[3], v2[3];
    _vectorize(p1, p2, v1);
    _vectorize(p1, p3, v2);
    _cross_product(v1, v2, n);
    for (int i = 0; i < 3; ++i)
    {
        n[i] /= 2.0;
    }
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
}
double _get_ratio(const double *n, const double *t)
{
    return _dot_product(n, t) / _dot_product(n, n);
}