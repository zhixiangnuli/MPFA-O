#include "roots.h"
#include "tran.h"   //尚不兼容无效网格
//#include "F.cpp"

#include <mdarray>
#include <mdspan>
#include <random>


namespace std
{
    using std::experimental::mdarray;
    using std::experimental::mdspan;
    using std::experimental::extents;
}

int main(int argc, char* argv[])
{
    using namespace resim;

    double* _verts = new double[nxx * nyy * nzz * 8 * 3];
    double* _verts1 = new double[nx * ny * nz * 8 * 3];

    std::mdspan<double, std::extents<size_t, -1, -1, -1, 8, 3>> verts(_verts, nzz, nyy, nxx);
    std::mdspan<double, std::extents<size_t, -1, -1, -1, 8, 3>> verts1(_verts1, nz, ny, nx);

    Ktensor* pem = new Ktensor[nxx * nyy * nzz]();
    Ktensor* pem1 = new Ktensor[nx * ny * nz]();

    double* p = new double[nx * ny * nz]();
    double* Qx = new double[nx * ny * nz]();
    double* Qy = new double[nx * ny * nz]();
    double* Qz = new double[nx * ny * nz]();

    std::fill(Qx, Qx + nx * ny * nz, std::numeric_limits<double>::quiet_NaN());
    std::fill(Qy, Qy + nx * ny * nz, std::numeric_limits<double>::quiet_NaN());
    std::fill(Qz, Qz + nx * ny * nz, std::numeric_limits<double>::quiet_NaN()); 

    int n, nnz, count;
    int i, j, k, iter;

    int* Ptr = new int[nx * ny * nz + 1]();
    int* Idx = new int[27 * nx * ny * nz]();
    double* Val = new double[27 * nx * ny * nz]();
    double* B = new double[nx * ny * nz]();
    double* x = new double[nx * ny * nz]();
    /*int Ptr[nx * ny * nz + 1];
    int Idx[7 * nx * ny * nz];
    double Val[7 * nx * ny * nz];
    double B[nx * ny * nz];
    double x[nx * ny * nz];*/

    double DX = 1.0;
    double DY = 1.0;
    double DZ = 1.0;

    //_get_3x3(_verts, pem);
    _get_2x2(pem);
    for (int k = 0; k < nzz; ++k) //笛卡尔网格
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

    _get_divide(_verts, _verts1, pem, pem1);

    double* cx = new double[nx * ny * nz];
    double* cy = new double[nx * ny * nz];
    double* cz = new double[nx * ny * nz];

    std::mdspan _cx(cx, nz, ny, nx);
    std::mdspan _cy(cy, nz, ny, nx);
    std::mdspan _cz(cz, nz, ny, nx);

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

    int* actnum = new int[nx * ny * nz];
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
                    //std::cout << "(" << i << "," << j << "," << k << ") : ";
                    count = 0;
                    double val = 0.0;
                    double F = 0.0;
                    double F1 = 0.0;
                    int    i0 = k * nx *ny + j * nx + i;

                    _get_F(pem1, p, i, j, k, _verts1, F);
                    //std::cout << "F=" << F << std::endl;

                    B[n] = -F;
                    if (k > 0)
                    {
                        int i1 = i0 - nx * ny;
                        p[i1] += epsilonP;
                        _get_F(pem1, p, i, j, k, _verts1, F1);
                        val = (F1 - F) / epsilonP;
                        p[i1] -= epsilonP;

                        if (val != 0.0) {
                            Val[nnz] = val;
                            count = count + 1;
                            Idx[nnz] = i1;
                            nnz = nnz + 1;
                        }

                        if (i > 0) {
                            int i1 = i0 - nx * ny - 1;
                            p[i1] += epsilonP;
                            _get_F(pem1, p, i, j, k, _verts1, F1);
                            val = (F1 - F) / epsilonP;
                            p[i1] -= epsilonP;

                            if (val != 0.0) {
                                Val[nnz] = val;
                                count = count + 1;
                                Idx[nnz] = i1;
                                nnz = nnz + 1;
                            }

                            if (j > 0) {
                                int i1 = i0 - nx * ny - nx - 1;
                                p[i1] += epsilonP;
                                _get_F(pem1, p, i, j, k, _verts1, F1);
                                val = (F1 - F) / epsilonP;
                                p[i1] -= epsilonP;

                                if (val != 0.0) {
                                    Val[nnz] = val;
                                    count = count + 1;
                                    Idx[nnz] = i1;
                                    nnz = nnz + 1;
                                }
                            }
                            if (j < ny - 1) {
                                int i1 = i0 - nx * ny + nx - 1;
                                p[i1] += epsilonP;
                                _get_F(pem1, p, i, j, k, _verts1, F1);
                                val = (F1 - F) / epsilonP;
                                p[i1] -= epsilonP;

                                if (val != 0.0) {
                                    Val[nnz] = val;
                                    count = count + 1;
                                    Idx[nnz] = i1;
                                    nnz = nnz + 1;
                                }
                            }
                        }
                        if (i < nx - 1) {
                            int i1 = i0 - nx * ny + 1;
                            p[i1] += epsilonP;
                            _get_F(pem1, p, i, j, k, _verts1, F1);
                            val = (F1 - F) / epsilonP;
                            p[i1] -= epsilonP;

                            if (val != 0.0) {
                                Val[nnz] = val;
                                count = count + 1;
                                Idx[nnz] = i1;
                                nnz = nnz + 1;
                            }

                            if (j > 0) {
                                int i1 = i0 - nx * ny - nx + 1;
                                p[i1] += epsilonP;
                                _get_F(pem1, p, i, j, k, _verts1, F1);
                                val = (F1 - F) / epsilonP;
                                p[i1] -= epsilonP;

                                if (val != 0.0) {
                                    Val[nnz] = val;
                                    count = count + 1;
                                    Idx[nnz] = i1;
                                    nnz = nnz + 1;
                                }
                            }
                            if (j < ny - 1) {
                                Val[nnz] = val;
                                int i1 = i0 - nx * ny + nx + 1;
                                p[i1] += epsilonP;
                                _get_F(pem1, p, i, j, k, _verts1, F1);
                                val = (F1 - F) / epsilonP;
                                p[i1] -= epsilonP;

                                if (val != 0.0) {
                                    Val[nnz] = val;
                                    count = count + 1;
                                    Idx[nnz] = i1;
                                    nnz = nnz + 1;
                                }
                            }
                        }
                        if (j > 0) {
                            int i1 = i0 - nx * ny - nx;
                            p[i1] += epsilonP;
                            _get_F(pem1, p, i, j, k, _verts1, F1);
                            val = (F1 - F) / epsilonP;
                            p[i1] -= epsilonP;

                            if (val != 0.0) {
                                Val[nnz] = val;
                                count = count + 1;
                                Idx[nnz] = i1;
                                nnz = nnz + 1;
                            }
                        }
                        if (j < ny - 1) {
                            int i1 = i0 - nx * ny + nx;
                            p[i1] += epsilonP;
                            _get_F(pem1, p, i, j, k, _verts1, F1);
                            val = (F1 - F) / epsilonP;
                            p[i1] -= epsilonP;

                            if (val != 0.0) {
                                Val[nnz] = val;
                                count = count + 1;
                                Idx[nnz] = i1;
                                nnz = nnz + 1;
                            }
                        }
                    }

                    p[i0] += epsilonP;
                    _get_F(pem1, p, i, j, k, _verts1, F1);
                    val = (F1 - F) / epsilonP;
                    p[i0] -= epsilonP;

                    if (val != 0.0) {
                        Val[nnz] = val;
                        count = count + 1;
                        Idx[nnz] = i0;
                        nnz = nnz + 1;
                    }

                    if (i > 0) {
                        int i1 = i0 - 1;
                        p[i1] += epsilonP;
                        _get_F(pem1, p, i, j, k, _verts1, F1);
                        val = (F1 - F) / epsilonP;
                        p[i1] -= epsilonP;

                        if (val != 0.0) {
                            Val[nnz] = val;
                            count = count + 1;
                            Idx[nnz] = i1;
                            nnz = nnz + 1;
                        }

                        if (j > 0) {
                            int i1 = i0 - nx - 1;
                            p[i1] += epsilonP;
                            _get_F(pem1, p, i, j, k, _verts1, F1);
                            val = (F1 - F) / epsilonP;
                            p[i1] -= epsilonP;

                            if (val != 0.0) {
                                Val[nnz] = val;
                                count = count + 1;
                                Idx[nnz] = i1;
                                nnz = nnz + 1;
                            }
                        }
                        if (j < ny - 1) {
                            int i1 = i0 + nx - 1;
                            p[i1] += epsilonP;
                            _get_F(pem1, p, i, j, k, _verts1, F1);
                            val = (F1 - F) / epsilonP;
                            p[i1] -= epsilonP;

                            if (val != 0.0) {
                                Val[nnz] = val;
                                count = count + 1;
                                Idx[nnz] = i1;
                                nnz = nnz + 1;
                            }
                        }
                    }
                    if (i < nx - 1) {
                        int i1 = i0 + 1;
                        p[i1] += epsilonP;
                        _get_F(pem1, p, i, j, k, _verts1, F1);
                        val = (F1 - F) / epsilonP;
                        p[i1] -= epsilonP;

                        if (val != 0.0) {
                            Val[nnz] = val;
                            count = count + 1;
                            Idx[nnz] = i1;
                            nnz = nnz + 1;
                        }

                        if (j > 0) {
                            int i1 = i0 - nx + 1;
                            p[i1] += epsilonP;
                            _get_F(pem1, p, i, j, k, _verts1, F1);
                            val = (F1 - F) / epsilonP;
                            p[i1] -= epsilonP;

                            if (val != 0.0) {
                                Val[nnz] = val;
                                count = count + 1;
                                Idx[nnz] = i1;
                                nnz = nnz + 1;
                            }
                        }
                        if (j < ny - 1) {
                            int i1 = i0 + nx + 1;
                            p[i1] += epsilonP;
                            _get_F(pem1, p, i, j, k, _verts1, F1);
                            val = (F1 - F) / epsilonP;
                            p[i1] -= epsilonP;

                            if (val != 0.0) {
                                Val[nnz] = val;
                                count = count + 1;
                                Idx[nnz] = i1;
                                nnz = nnz + 1;
                            }
                        }
                    }
                    if (j > 0) {
                        int i1 = i0 - nx;
                        p[i1] += epsilonP;
                        _get_F(pem1, p, i, j, k, _verts1, F1);
                        val = (F1 - F) / epsilonP;
                        p[i1] -= epsilonP;

                        if (val != 0.0) {
                            Val[nnz] = val;
                            count = count + 1;
                            Idx[nnz] = i1;
                            nnz = nnz + 1;
                        }
                    }
                    if (j < ny - 1) {
                        int i1 = i0 + nx;
                        p[i1] += epsilonP;
                        _get_F(pem1, p, i, j, k, _verts1, F1);
                        val = (F1 - F) / epsilonP;
                        p[i1] -= epsilonP;

                        if (val != 0.0) {
                            Val[nnz] = val;
                            count = count + 1;
                            Idx[nnz] = i1;
                            nnz = nnz + 1;
                        }
                    }

                    if (k < nz - 1)
                    {
                        int i1 = i0 + nx * ny;
                        p[i1] += epsilonP;
                        _get_F(pem1, p, i, j, k, _verts1, F1);
                        val = (F1 - F) / epsilonP;
                        p[i1] -= epsilonP;

                        if (val != 0.0) {
                            Val[nnz] = val;
                            count = count + 1;
                            Idx[nnz] = i1;
                            nnz = nnz + 1;
                        }

                        if (i > 0) {
                            int i1 = i0 + nx * ny - 1;
                            p[i1] += epsilonP;
                            _get_F(pem1, p, i, j, k, _verts1, F1);
                            val = (F1 - F) / epsilonP;
                            p[i1] -= epsilonP;

                            if (val != 0.0) {
                                Val[nnz] = val;
                                count = count + 1;
                                Idx[nnz] = i1;
                                nnz = nnz + 1;
                            }

                            if (j > 0) {
                                int i1 = i0 + nx * ny - nx - 1;
                                p[i1] += epsilonP;
                                _get_F(pem1, p, i, j, k, _verts1, F1);
                                val = (F1 - F) / epsilonP;
                                p[i1] -= epsilonP;

                                if (val != 0.0) {
                                    Val[nnz] = val;
                                    count = count + 1;
                                    Idx[nnz] = i1;
                                    nnz = nnz + 1;
                                }
                            }
                            if (j < ny - 1) {
                                int i1 = i0 + nx * ny + nx - 1;
                                p[i1] += epsilonP;
                                _get_F(pem1, p, i, j, k, _verts1, F1);
                                val = (F1 - F) / epsilonP;
                                p[i1] -= epsilonP;

                                if (val != 0.0) {
                                    Val[nnz] = val;
                                    count = count + 1;
                                    Idx[nnz] = i1;
                                    nnz = nnz + 1;
                                }
                            }
                        }
                        if (i < nx - 1) {
                            int i1 = i0 + nx * ny + 1;
                            p[i1] += epsilonP;
                            _get_F(pem1, p, i, j, k, _verts1, F1);
                            val = (F1 - F) / epsilonP;
                            p[i1] -= epsilonP;

                            if (val != 0.0) {
                                Val[nnz] = val;
                                count = count + 1;
                                Idx[nnz] = i1;
                                nnz = nnz + 1;
                            }

                            if (j > 0) {
                                int i1 = i0 + nx * ny - nx + 1;
                                p[i1] += epsilonP;
                                _get_F(pem1, p, i, j, k, _verts1, F1);
                                val = (F1 - F) / epsilonP;
                                p[i1] -= epsilonP;

                                if (val != 0.0) {
                                    Val[nnz] = val;
                                    count = count + 1;
                                    Idx[nnz] = i1;
                                    nnz = nnz + 1;
                                }
                            }
                            if (j < ny - 1) {
                                int i1 = i0 + nx * ny + nx + 1;
                                p[i1] += epsilonP;
                                _get_F(pem1, p, i, j, k, _verts1, F1);
                                val = (F1 - F) / epsilonP;
                                p[i1] -= epsilonP;

                                if (val != 0.0) {
                                    Val[nnz] = val;
                                    count = count + 1;
                                    Idx[nnz] = i1;
                                    nnz = nnz + 1;
                                }
                            }
                        }
                        if (j > 0) {
                            int i1 = i0 + nx * ny - nx;
                            p[i1] += epsilonP;
                            _get_F(pem1, p, i, j, k, _verts1, F1);
                            val = (F1 - F) / epsilonP;
                            p[i1] -= epsilonP;

                            if (val != 0.0) {
                                Val[nnz] = val;
                                count = count + 1;
                                Idx[nnz] = i1;
                                nnz = nnz + 1;
                            }
                        }
                        if (j < ny - 1) {
                            int i1 = i0 + nx * ny + nx;
                            p[i1] += epsilonP;
                            _get_F(pem1, p, i, j, k, _verts1, F1);
                            val = (F1 - F) / epsilonP;
                            p[i1] -= epsilonP;

                            if (val != 0.0) {
                                Val[nnz] = val;
                                count = count + 1;
                                Idx[nnz] = i1;
                                nnz = nnz + 1;
                            }
                        }
                    }

                    Ptr[n + 1] = Ptr[n] + count;

                    n = n + 1;
                    //std::cout << std::endl;
                }
            }
        }
        
        //int mr = std::min(nx * ny * nz / 100, 40);
        //pmgmres_ilu_cr(nx* ny* nz, nnz, Ptr, Idx, Val, x, B, 1000, 1000, 1e-8, 1e-8);

        pmgmres_ilu_cr(nx* ny* nz, nnz, Ptr, Idx, Val, x, B, 1000, 1000, 1e-5, 1e-5);

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

        double Q = 0.0;
        for (int k = 0; k < nz; k++)
        {
            for (int j = 0; j < ny; j++)
            {
                const double* pp[4];
                double        S = 0.0, D = 0.0;
                double        pc[3];
                int           i0 = k * nx * ny + j * nx + nx - 1;

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

        std::cout << "Q=" << Q << "\n";
        std::cout << "res=" << res << "\n";
        if (res < 1.0e-8) break;
    }
    
    delete[] pem, pem1, Qx, Qx, Qx, p, Ptr, Idx, Val, B, x, _verts, _verts1, cx, cy, cz;

    return 0;
}