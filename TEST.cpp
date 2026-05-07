#define _CRT_SECURE_NO_WARNINGS
#include "roots.h"
#include "tran.h"
#include <mdarray>
#include <mdspan>
#include <random>
#include <chrono>
namespace std
{
    using std::experimental::extents;
    using std::experimental::mdarray;
    using std::experimental::mdspan;
}

int main(int argc, char *argv[])
{
    auto t1 = std::chrono::high_resolution_clock::now();
    using namespace resim;

    double *_verts = new double[nxx * nyy * nzz * 8 * 3];
    _get_verts(_verts);

    Ktensor *pem = new Ktensor[nxx * nyy * nzz]();
    double *kc = new double[nxx * nyy * nzz]{95.80, 65.66, 6.43, 9.12,
                                             25.76, 90.42, 34.83, 49.96,
                                             61.60, 97.87, 62.39, 1.38,
                                             23.73, 84.04, 6.59, 64.21,
                                             79.16, 21.12, 72.82, 37.08,
                                             34.73, 47.48, 52.43, 11.16,
                                             8.35, 69.01, 85.74, 82.65,
                                             78.66, 28.72, 93.24, 72.08,
                                             71.62, 52.34, 38.13, 35.94};
    double *kd = new double[nxx * nyy * nzz]{18.13, 0.46, 0.18, 0.30,
                                             3.07, 16.76, 4.04, 2.68,
                                             1.94, 14.43, 10.61, 0.23,
                                             0.25, 2.91, 0.51, 12.80,
                                             14.87, 3.27, 1.00, 8.86,
                                             8.68, 2.39, 1.68, 1.74,
                                             0.63, 16.38, 2.22, 10.97,
                                             5.18, 2.55, 18.33, 12.78,
                                             14.53, 2.17, 8.33, 7.02};
    _get_pem(pem, kc, kd);

    double Q = 0.0;
    int iter = 0;

    std::ofstream output("output_uniform.txt");

    for (int i = 1; i <= 8; ++i)
    {
        _get_Qx(pem, _verts, i, iter, Q);
        output << Q << ',' << std::flush;
    }

    _get_Qx(pem, _verts, 16, iter, Q);
    output << Q << ',' << std::flush;
    _get_Qx(pem, _verts, 32, iter, Q);
    output << Q << std::endl;

    output.close();
    delete[] pem, _verts;
    delete[] kc, kd;
    // 结束计时
    auto t2 = std::chrono::high_resolution_clock::now();

    // 计算毫秒
    auto ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);

    // 计算秒
    std::chrono::duration<double> sec = t2 - t1;

    std::cout << "运行时间: " << ms_int.count() << "ms\n";
    std::cout << "运行时间: " << sec.count() << "s\n";
    return 0;
}