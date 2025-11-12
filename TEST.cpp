#define _CRT_SECURE_NO_WARNINGS
#include "roots.h"
#include "tran.h" //尚不兼容无效网格
// #include "F.cpp"

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

    FILE *fp = std::fopen("Qx_MPFA.INC", "w");
    /*_get_2x2(pem);*/
    for (int divn = 1; divn < 8; divn++)
    {
        std::fprintf(fp, "1.1\n");
        // pemx = 10.0;
        // pemy = 20.0;
        // pemz = 20.0;
        // pemxy = 2.0;
        // pemyz = 2.0;
        // pemxz = 2.0;
        Lamda = 10;
        _get_3x3(_verts, pem, pemx, pemy, pemz, pemxy, pemyz, pemxz, Lamda);
        _Outputvtk(_verts, (double *)pem);
        std::fprintf(fp, "divn = %ld\n", divn);
        _get_Qx(pem, _verts, divn, iter, Q);
        std::fprintf(fp, "iter = %ld\n", iter);
        std::fprintf(fp, "Q = %le\n", Q);
        std::fprintf(fp, "/\n");

        /*std::fprintf(fp, "1.2\n");
        pemx = 20.0; pemy = 20.0; pemz = 20.0;
        pemxy = 0.0; pemyz = 0.0; pemxz = 0.0;
        Lamda = 100;
        _get_3x3(_verts, pem, pemx, pemy, pemz, pemxy, pemyz, pemxz, Lamda);
        std::fprintf(fp, "divn = %ld\n", divn);
        _get_Qx(pem, _verts, divn, iter, Q);
        std::fprintf(fp, "iter = %ld\n", iter);
        std::fprintf(fp, "Q = %le\n", Q);
        std::fprintf(fp, "/\n");

        std::fprintf(fp, "2.1\n");
        pemx = 1.0; pemy = 1.0; pemz = 1.0;
        pemxy = 0.125; pemyz = 0.125; pemxz = 0.125;
        Lamda = 20;
        _get_3x3(_verts, pem, pemx, pemy, pemz, pemxy, pemyz, pemxz, Lamda);
        std::fprintf(fp, "divn = %ld\n", divn);
        _get_Qx(pem, _verts, divn, iter, Q);
        std::fprintf(fp, "iter = %ld\n", iter);
        std::fprintf(fp, "Q = %le\n", Q);
        std::fprintf(fp, "/\n");

        std::fprintf(fp, "2.2\n");
        pemx = 1.0; pemy = 1.0; pemz = 1.0;
        pemxy = 0.25; pemyz = 0.25; pemxz = 0.25;
        Lamda = 20;
        _get_3x3(_verts, pem, pemx, pemy, pemz, pemxy, pemyz, pemxz, Lamda);
        std::fprintf(fp, "divn = %ld\n", divn);
        _get_Qx(pem, _verts, divn, iter, Q);
        std::fprintf(fp, "iter = %ld\n", iter);
        std::fprintf(fp, "Q = %le\n", Q);
        std::fprintf(fp, "/\n");

        std::fprintf(fp, "2.3\n");
        pemx = 1.0; pemy = 1.0; pemz = 1.0;
        pemxy = 0.5; pemyz = 0.5; pemxz = 0.5;
        Lamda = 20;
        _get_3x3(_verts, pem, pemx, pemy, pemz, pemxy, pemyz, pemxz, Lamda);
        std::fprintf(fp, "divn = %ld\n", divn);
        _get_Qx(pem, _verts, divn, iter, Q);
        std::fprintf(fp, "iter = %ld\n", iter);
        std::fprintf(fp, "Q = %le\n", Q);
        std::fprintf(fp, "/\n");*/
    }
    // divn = 16
    std::fprintf(fp, "1.1\n");
    pemx = 20.0;
    pemy = 20.0;
    pemz = 20.0;
    // pemxy = 2.0; pemyz = 1.0; pemxz = 2.0;
    Lamda = 10;
    _get_3x3(_verts, pem, pemx, pemy, pemz, pemxy, pemyz, pemxz, Lamda);
    std::fprintf(fp, "divn = %ld\n", 16);
    _get_Qx(pem, _verts, 16, iter, Q);
    std::fprintf(fp, "iter = %ld\n", iter);
    std::fprintf(fp, "Q = %le\n", Q);
    std::fprintf(fp, "/\n");

    std::fprintf(fp, "1.2\n");
    pemx = 20.0;
    pemy = 20.0;
    pemz = 20.0;
    // pemxy = 2.0; pemyz = 1.0; pemxz = 2.0;
    Lamda = 100;
    _get_3x3(_verts, pem, pemx, pemy, pemz, pemxy, pemyz, pemxz, Lamda);
    std::fprintf(fp, "divn = %ld\n", 16);
    _get_Qx(pem, _verts, 16, iter, Q);
    std::fprintf(fp, "iter = %ld\n", iter);
    std::fprintf(fp, "Q = %le\n", Q);
    std::fprintf(fp, "/\n");

    /*std::fprintf(fp, "2.1\n");
    pemx = 1.0; pemy = 1.0; pemz = 1.0;
    pemxy = 0.125; pemyz = 0.125; pemxz = 0.125;
    Lamda = 20;
    _get_3x3(_verts, pem, pemx, pemy, pemz, pemxy, pemyz, pemxz, Lamda);
    std::fprintf(fp, "divn = %ld\n", 16);
    _get_Qx(pem, _verts, 16, iter, Q);
    std::fprintf(fp, "iter = %ld\n", iter);
    std::fprintf(fp, "Q = %le\n", Q);
    std::fprintf(fp, "/\n");

    std::fprintf(fp, "2.2\n");
    pemx = 1.0; pemy = 1.0; pemz = 1.0;
    pemxy = 0.25; pemyz = 0.25; pemxz = 0.25;
    Lamda = 20;
    _get_3x3(_verts, pem, pemx, pemy, pemz, pemxy, pemyz, pemxz, Lamda);
    std::fprintf(fp, "divn = %ld\n", 16);
    _get_Qx(pem, _verts, 16, iter, Q);
    std::fprintf(fp, "iter = %ld\n", iter);
    std::fprintf(fp, "Q = %le\n", Q);
    std::fprintf(fp, "/\n");

    std::fprintf(fp, "2.3\n");
    pemx = 1.0; pemy = 1.0; pemz = 1.0;
    pemxy = 0.5; pemyz = 0.5; pemxz = 0.5;
    Lamda = 20;
    _get_3x3(_verts, pem, pemx, pemy, pemz, pemxy, pemyz, pemxz, Lamda);
    std::fprintf(fp, "divn = %ld\n", 16);
    _get_Qx(pem, _verts, 16, iter, Q);
    std::fprintf(fp, "iter = %ld\n", iter);
    std::fprintf(fp, "Q = %le\n", Q);
    std::fprintf(fp, "/\n");

    //divn = 32
    std::fprintf(fp, "1.1\n");
    pemx = 10.0; pemy = 20.0; pemz = 20.0;
    pemxy = 2.0; pemyz = 4.0; pemxz = 5.0;
    Lamda = 10;
    _get_3x3(_verts, pem, pemx, pemy, pemz, pemxy, pemyz, pemxz, Lamda);
    std::fprintf(fp, "divn = %ld\n", 32);
    _get_Qx(pem, _verts, 32, iter, Q);
    std::fprintf(fp, "iter = %ld\n", iter);
    std::fprintf(fp, "Q = %le\n", Q);
    std::fprintf(fp, "/\n");

    std::fprintf(fp, "1.2\n");
    pemx = 10.0; pemy = 20.0; pemz = 20.0;
    pemxy = 2.0; pemyz = 4.0; pemxz = 5.0;
    Lamda = 100;
    _get_3x3(_verts, pem, pemx, pemy, pemz, pemxy, pemyz, pemxz, Lamda);
    std::fprintf(fp, "divn = %ld\n", 32);
    _get_Qx(pem, _verts, 32, iter, Q);
    std::fprintf(fp, "iter = %ld\n", iter);
    std::fprintf(fp, "Q = %le\n", Q);
    std::fprintf(fp, "/\n");

    std::fprintf(fp, "2.1\n");
    pemx = 1.0; pemy = 1.0; pemz = 1.0;
    pemxy = 0.125; pemyz = 0.125; pemxz = 0.125;
    Lamda = 20;
    _get_3x3(_verts, pem, pemx, pemy, pemz, pemxy, pemyz, pemxz, Lamda);
    std::fprintf(fp, "divn = %ld\n", 32);
    _get_Qx(pem, _verts, 32, iter, Q);
    std::fprintf(fp, "iter = %ld\n", iter);
    std::fprintf(fp, "Q = %le\n", Q);
    std::fprintf(fp, "/\n");

    std::fprintf(fp, "2.2\n");
    pemx = 1.0; pemy = 1.0; pemz = 1.0;
    pemxy = 0.25; pemyz = 0.25; pemxz = 0.25;
    Lamda = 20;
    _get_3x3(_verts, pem, pemx, pemy, pemz, pemxy, pemyz, pemxz, Lamda);
    std::fprintf(fp, "divn = %ld\n", 32);
    _get_Qx(pem, _verts, 32, iter, Q);
    std::fprintf(fp, "iter = %ld\n", iter);
    std::fprintf(fp, "Q = %le\n", Q);
    std::fprintf(fp, "/\n");

    std::fprintf(fp, "2.3\n");
    pemx = 1.0; pemy = 1.0; pemz = 1.0;
    pemxy = 0.5; pemyz = 0.5; pemxz = 0.5;
    Lamda = 20;
    _get_3x3(_verts, pem, pemx, pemy, pemz, pemxy, pemyz, pemxz, Lamda);
    std::fprintf(fp, "divn = %ld\n", 32);
    _get_Qx(pem, _verts, 32, iter, Q);
    std::fprintf(fp, "iter = %ld\n", iter);
    std::fprintf(fp, "Q = %le\n", Q);
    std::fprintf(fp, "/\n");
    */
    std::fclose(fp);

    delete[] pem, _verts;

    return 0;
}