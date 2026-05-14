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
    double *kc = new double[nxx * nyy * nzz]{};
    double *kd = new double[nxx * nyy * nzz]{};
    _get_pem(pem, kc, kd);

    double Q = 0.0;
    int iter = 0;

    std::ofstream output("output_uniform.txt");

    for (int i = 3; i <= 3; ++i)
    {
        _get_Qx(pem, _verts, i, iter, Q);
        output << Q << ',' << std::flush;
    }
    return 0;
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