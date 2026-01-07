#pragma once

#include <cmath>
#include <math.h>
#include <complex>
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
// 定义Point3D结构体
struct Point3D
{
    double x;
    double y;
    double z;

    Point3D(double x_ = 0.0, double y_ = 0.0, double z_ = 0.0)
        : x(x_), y(y_), z(z_) {}
};

// 定义Vector3D结构体
struct Vector3D
{
    double dx; // X方向分量
    double dy; // Y方向分量
    double dz; // Z方向分量

    Vector3D(double dx_ = 0.0, double dy_ = 0.0, double dz_ = 0.0)
        : dx(dx_), dy(dy_), dz(dz_) {}

    // 计算向量长度
    double length() const
    {
        return sqrt(dx * dx + dy * dy + dz * dz);
    }

    // 归一化向量
    Vector3D normalize() const
    {
        double len = length();
        if (len > 0.0)
        {
            return Vector3D(dx / len, dy / len, dz / len);
        }
        return Vector3D(0.0, 0.0, 0.0);
    }

    // 向量加法
    Vector3D operator+(const Vector3D &other) const
    {
        return Vector3D(dx + other.dx, dy + other.dy, dz + other.dz);
    }

    // 向量减法
    Vector3D operator-(const Vector3D &other) const
    {
        return Vector3D(dx - other.dx, dy - other.dy, dz - other.dz);
    }

    // 标量乘法
    Vector3D operator*(double scalar) const
    {
        return Vector3D(dx * scalar, dy * scalar, dz * scalar);
    }

    // 标量除法
    Vector3D operator/(double scalar) const
    {
        if (scalar != 0.0)
        {
            return Vector3D(dx / scalar, dy / scalar, dz / scalar);
        }
        return Vector3D(0.0, 0.0, 0.0);
    }
};

// 点加向量得到新点
Point3D operator+(const Point3D &point, const Vector3D &vec)
{
    return Point3D(point.x + vec.dx, point.y + vec.dy, point.z + vec.dz);
}

// 点减点得到向量
Vector3D operator-(const Point3D &p1, const Point3D &p2)
{
    return Vector3D(p1.x - p2.x, p1.y - p2.y, p1.z - p2.z);
}

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
    using namespace std;
    std::mdspan<double, std::extents<size_t, -1, -1, -1, 8, 3>> verts(_verts, nzz, nyy, nxx);

    // 首先，定义网格尺寸
    double M_PI = Pi;

    {

        // 定义网格节点坐标 (尺寸为(nzz+1) × (nyy+1) × (nxx+1))
        vector<vector<vector<Point3D>>> grid_nodes(nzz + 1,
                                                   vector<vector<Point3D>>(nyy + 1,
                                                                           vector<Point3D>(nxx + 1)));

        // 定义顶层和底层控制点
        vector<vector<Point3D>> top_control_points(nyy + 1, vector<Point3D>(nxx + 1));
        vector<vector<Point3D>> bottom_control_points(nyy + 1, vector<Point3D>(nxx + 1));

        // 初始化顶层控制点（增加网格尺寸变化和波浪形）
        for (int j = 0; j <= nyy; ++j)
        {
            for (int i = 0; i <= nxx; ++i)
            {
                // 增加网格尺寸的不规则性：使用变化的基础间距
                double x_base = i * (80.0 + 40.0 * sin(i * 0.5) * cos(j * 0.3)); // 80-120之间变化
                double y_base = j * (75.0 + 30.0 * sin(i * 0.4) * cos(j * 0.5)); // 75-105之间变化

                // 增加波浪形顶层，适度提高振幅
                double amplitude_x = 25.0;
                double amplitude_y = 20.0;
                double amplitude_z = 35.0;

                top_control_points[j][i].x = x_base +
                                             amplitude_x * sin(i * 0.8) * cos(j * 0.5);
                top_control_points[j][i].y = y_base +
                                             amplitude_y * sin(i * 0.6) * cos(j * 0.7);
                top_control_points[j][i].z = 70.0 +
                                             amplitude_z * sin(i * 0.7) * sin(j * 0.6);
            }
        }
        top_control_points[1][3].y -= 30.0;

        // 初始化底层控制点（显著增加倾斜和扭曲）
        for (int j = 0; j <= nyy; ++j)
        {
            for (int i = 0; i <= nxx; ++i)
            {
                // 底层网格尺寸变化与顶层对应但有所偏移
                double x_base = i * (85.0 + 35.0 * sin(i * 0.6) * cos(j * 0.4));
                double y_base = j * (80.0 + 25.0 * sin(i * 0.5) * cos(j * 0.6));

                // 增加底层变形程度
                double amplitude_x = 35.0;
                double amplitude_y = 30.0;

                bottom_control_points[j][i].x = x_base +
                                                amplitude_x * sin(i * 1.0) * cos(j * 0.8);
                bottom_control_points[j][i].y = y_base +
                                                amplitude_y * sin(i * 0.8) * cos(j * 1.0);

                // 显著增加底层深度变化和倾斜
                double base_z = 400.0;

                // 显著增加倾斜程度
                double tilt_x = 50.0 * sin(i * 0.6);                                         // x方向倾斜增加到50
                double tilt_y = 45.0 * (1.0 - j / (double)nyy) * (1.0 + 0.5 * sin(i * 0.4)); // y方向倾斜增加到45，且变化
                double local_variation = 30.0 * sin(i * 0.7) * cos(j * 0.8);                 // 增加局部变化

                // 添加轻微的断层模拟（非真实断层）
                double offset_effect = 0.0;
                if (i > nxx / 2)
                {
                    offset_effect = -20.0 * sin((i - nxx / 2) * 0.5); // 轻微偏移
                }

                bottom_control_points[j][i].z = base_z + tilt_x + tilt_y + local_variation + offset_effect;
            }
        }
        bottom_control_points[1][3].y += 30.0;
        bottom_control_points[0][3].y += 20.0;
        bottom_control_points[3][3].y -= 20.0;
        bottom_control_points[1][2].z += 25.0;
        bottom_control_points[1][1].z += 25.0;
        bottom_control_points[1][1].y += 25.0;
        bottom_control_points[3][0].y += 25.0;
        bottom_control_points[3][1].y += 25.0;
        // 定义非线性层深度 (nzz+1层，从0到1)
        vector<double> layer_depths(nzz + 1);

        // 创建不均匀的层厚度分布，增加层间差异
        for (int k = 0; k <= nzz; ++k)
        {
            double normalized = (double)k / nzz;

            // 创建有明显厚度变化的层
            if (k == 0)
            {
                layer_depths[k] = 0.0;
            }
            else if (k == nzz)
            {
                layer_depths[k] = 1.0;
            }
            else
            {
                // 使用函数创建不均匀的层厚度
                double exp_factor = exp(1.5 * normalized - 0.75) / exp(0.75);
                double sin_factor = 0.2 * sin(normalized * 4.0 * M_PI);

                // 调整权重，增加非线性
                layer_depths[k] = normalized * 0.4 + exp_factor * 0.5 + sin_factor * 0.1;
            }
        }

        // 对层深度进行归一化
        double min_depth = *min_element(layer_depths.begin(), layer_depths.end());
        double max_depth = *max_element(layer_depths.begin(), layer_depths.end());
        for (int k = 0; k <= nzz; ++k)
        {
            layer_depths[k] = (layer_depths[k] - min_depth) / (max_depth - min_depth);
        }

        // 添加随机扰动到层深度
        srand(12345);
        for (int k = 1; k < nzz; ++k)
        {
            double random_perturbation = 0.12 * ((double)rand() / RAND_MAX - 0.5);
            layer_depths[k] += random_perturbation;
        }

        // 重新排序以确保单调递增
        sort(layer_depths.begin(), layer_depths.end());

        // 计算所有网格节点的坐标，添加全局扰动
        for (int k = 0; k <= nzz; ++k)
        {
            for (int j = 0; j <= nyy; ++j)
            {
                for (int i = 0; i <= nxx; ++i)
                {
                    // 根据角点网格原理，每个节点都在连接顶层和底层对应点的直线上
                    Point3D top_point = top_control_points[j][i];
                    Point3D bottom_point = bottom_control_points[j][i];

                    // 使用层深度比例进行插值
                    double ratio = layer_depths[k];

                    // 基本插值
                    grid_nodes[k][j][i].x = top_point.x + (bottom_point.x - top_point.x) * ratio;
                    grid_nodes[k][j][i].y = top_point.y + (bottom_point.y - top_point.y) * ratio;
                    grid_nodes[k][j][i].z = top_point.z + (bottom_point.z - top_point.z) * ratio;

                    // 添加全局扰动，增加不规则性但不破坏连续性
                    double freq_x = 0.015;
                    double freq_y = 0.013;
                    double freq_z = 0.008;

                    // 使用多层扰动叠加，但幅度适中
                    double perturbation1 = 12.0 * sin(grid_nodes[k][j][i].x * freq_x) *
                                           cos(grid_nodes[k][j][i].y * freq_y) *
                                           sin(grid_nodes[k][j][i].z * freq_z);

                    double perturbation2 = 8.0 * sin(grid_nodes[k][j][i].x * freq_x * 2.0) *
                                           cos(grid_nodes[k][j][i].y * freq_y * 1.8) *
                                           sin(grid_nodes[k][j][i].z * freq_z * 1.2);

                    double total_perturbation = perturbation1 * 0.6 + perturbation2 * 0.4;

                    // 为不同方向添加扰动
                    grid_nodes[k][j][i].x += total_perturbation * 0.25;
                    grid_nodes[k][j][i].y += total_perturbation * 0.18;
                    grid_nodes[k][j][i].z += total_perturbation * 0.12;

                    // 添加基于层高的扰动，使不同层有不同的扭曲
                    double layer_perturbation = 10.0 * sin(ratio * 6.0) * cos(i * 0.4 + j * 0.3 + k * 0.5);
                    grid_nodes[k][j][i].x += layer_perturbation * 0.15;
                    grid_nodes[k][j][i].y += layer_perturbation * 0.1;
                    grid_nodes[k][j][i].z += layer_perturbation * 0.05;
                }
            }
        }

        // 现在，基于这些节点计算每个网格的八个角点
        for (int k = 0; k < nzz; ++k)
        {
            for (int j = 0; j < nyy; ++j)
            {
                for (int i = 0; i < nxx; ++i)
                {
                    // 获取当前网格的八个节点
                    Point3D node000 = grid_nodes[k][j][i];         // 底部左下
                    Point3D node100 = grid_nodes[k][j][i + 1];     // 底部右下
                    Point3D node010 = grid_nodes[k][j + 1][i];     // 底部左上
                    Point3D node110 = grid_nodes[k][j + 1][i + 1]; // 底部右上

                    Point3D node001 = grid_nodes[k + 1][j][i];         // 顶部左下
                    Point3D node101 = grid_nodes[k + 1][j][i + 1];     // 顶部右下
                    Point3D node011 = grid_nodes[k + 1][j + 1][i];     // 顶部左上
                    Point3D node111 = grid_nodes[k + 1][j + 1][i + 1]; // 顶部右上

                    // 将节点坐标赋值给网格顶点
                    // 底部四个角点
                    verts(k, j, i, 0, 0) = node000.x;
                    verts(k, j, i, 0, 1) = node000.y;
                    verts(k, j, i, 0, 2) = node000.z;

                    verts(k, j, i, 1, 0) = node100.x;
                    verts(k, j, i, 1, 1) = node100.y;
                    verts(k, j, i, 1, 2) = node100.z;

                    verts(k, j, i, 2, 0) = node010.x;
                    verts(k, j, i, 2, 1) = node010.y;
                    verts(k, j, i, 2, 2) = node010.z;

                    verts(k, j, i, 3, 0) = node110.x;
                    verts(k, j, i, 3, 1) = node110.y;
                    verts(k, j, i, 3, 2) = node110.z;

                    // 顶部四个角点
                    verts(k, j, i, 4, 0) = node001.x;
                    verts(k, j, i, 4, 1) = node001.y;
                    verts(k, j, i, 4, 2) = node001.z;

                    verts(k, j, i, 5, 0) = node101.x;
                    verts(k, j, i, 5, 1) = node101.y;
                    verts(k, j, i, 5, 2) = node101.z;

                    verts(k, j, i, 6, 0) = node011.x;
                    verts(k, j, i, 6, 1) = node011.y;
                    verts(k, j, i, 6, 2) = node011.z;

                    verts(k, j, i, 7, 0) = node111.x;
                    verts(k, j, i, 7, 1) = node111.y;
                    verts(k, j, i, 7, 2) = node111.z;
                }
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

    std::ofstream outputFile("model_data_deepseek.vtk");

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
                // outputFile << perm[n] << std::endl;
                outputFile << (k + j + i) % 2 << std::endl;
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
void _get_angle(const double *Axis_x, const double *Axis_y, const double *v, const double r, double &ang)
{
    ang = std::acos(_dot_product(Axis_x, v) / r);
    if (_dot_product(Axis_y, v) < 0.0)
        ang = -ang;
}
void _get_angle(const double x, const double y, const double r, double &ang)
{
    ang = std::acos(x / r);
    if (y < 0.0)
        ang = -ang;
}
