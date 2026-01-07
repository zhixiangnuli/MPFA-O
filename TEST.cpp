#define _CRT_SECURE_NO_WARNINGS
#include "roots.h"
#include "tran.h"
#include <mdarray>
#include <mdspan>
#include <random>

namespace std
{
    using std::experimental::extents;
    using std::experimental::mdarray;
    using std::experimental::mdspan;
}

int main(int argc, char *argv[])
{
    using namespace resim;

    double *_verts = new double[nxx * nyy * nzz * 8 * 3];
    std::mdspan<double, std::extents<size_t, -1, -1, -1, 8, 3>> verts(_verts, nzz, nyy, nxx);
    Ktensor *pem = new Ktensor[nxx * nyy * nzz]();

    double Q = 0.0;
    int iter = 0;
    int Lamda = 10;
    double pemx = 10.0, pemy = 20.0, pemz = 20.0;
    double pemxy = 2.0, pemyz = 4.0, pemxz = 5.0;
    _get_3x3(_verts, pem, pemx, pemy, pemz, pemxy, pemyz, pemxz, Lamda);
    _Outputvtk(_verts, (double *)pem);
    for (int divn = 1; divn <= 8; divn++)
        _get_Qx(pem, _verts, divn, iter, Q);
    _get_Qx(pem, _verts, 16, iter, Q);
    _get_Qx(pem, _verts, 32, iter, Q);
    delete[] pem, _verts;
    return 0;
}