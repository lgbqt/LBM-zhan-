#include <iostream>
#include <cmath>
#include <iomanip>
#include <algorithm>

using namespace std;


// ============================================================
// 1. 基本参数
// ============================================================

const double Lx = 2.0;
const double Ly = 2.0;

const double W  = 0.04;
const double dx = 0.01;

const double c  = 100.0;
const double dt = dx / c;

const int NX = static_cast<int>(Lx / dx);
const int NY = static_cast<int>(Ly / dx);

const int Q5 = 5;
const int Q9 = 9;

const double cs2 = c * c / 3.0;


// ============================================================
// 2. 物理参数
// ============================================================

const double nu    = 0.1;
const double sigma = 0.001;

const double rhoA = 1000.0;
const double rhoB = 1.0;

const double phiA =  1.0;
const double phiB = -1.0;

const double Rfinal = 0.5;

const double PI =
    3.14159265358979323846;


// ============================================================
// 3. 自由能参数
//
// beta  = 12 sigma / [W (phiA-phiB)^4]
//
// kappa = 3 sigma W / [2 (phiA-phiB)^2]
// ============================================================

const double betaPF =
    12.0 * sigma
    /
    (
        W
        *
        pow(phiA - phiB, 4.0)
    );

const double kappaPF =
    1.5
    *
    sigma
    *
    W
    /
    (
        (phiA - phiB)
        *
        (phiA - phiB)
    );


// ============================================================
// 4. 初始方形边长
// ============================================================

const double squareLength =
    sqrt(
        PI * Rfinal * Rfinal
        +
        (4.0 - PI) * W * W / 4.0
    );


// ============================================================
// 5. 流场松弛参数
// ============================================================

const double su1 =
    1.0
    /
    (
        nu / (cs2 * dt)
        +
        0.5
    );

const double Su[Q9] =
{
    1.0,
    su1,
    su1,
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
    1.0
};


// ============================================================
// 6. D2Q5
// ============================================================

const int cx5[Q5] =
{
     0,  1,  0, -1,  0
};

const int cy5[Q5] =
{
     0,  0,  1,  0, -1
};

const double w5[Q5] =
{
    1.0 / 3.0,
    1.0 / 6.0,
    1.0 / 6.0,
    1.0 / 6.0,
    1.0 / 6.0
};


// ============================================================
// 7. D2Q9
// ============================================================

const int cx9[Q9] =
{
     0,  1,  0, -1,  0,  1, -1, -1,  1
};

const int cy9[Q9] =
{
     0,  0,  1,  0, -1,  1,  1, -1, -1
};

const double w9[Q9] =
{
    4.0 / 9.0,

    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,

    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0
};


// ============================================================
// 8. Zhan Appendix A: D2Q9 transformation matrix
// ============================================================

const double M9[Q9][Q9] =
{
    {1, 1, 1, 1, 1, 1, 1, 1, 1},
    {0, 1, 0,-1, 0, 1,-1,-1, 1},
    {0, 0, 1, 0,-1, 1, 1,-1,-1},
    {0, 1, 0, 1, 0, 1, 1, 1, 1},
    {0, 0, 0, 0, 0, 1,-1, 1,-1},
    {0, 0, 1, 0, 1, 1, 1, 1, 1},
    {0, 0, 0, 0, 0, 1, 1,-1,-1},
    {0, 0, 0, 0, 0, 1,-1,-1, 1},
    {0, 0, 0, 0, 0, 1, 1, 1, 1}
};


const double Minv9[Q9][Q9] =
{
    {1, 0, 0,-1, 0,-1, 0, 0, 1},

    {0, 0.5, 0, 0.5, 0, 0, 0,-0.5,-0.5},

    {0, 0, 0.5, 0, 0, 0.5,-0.5, 0,-0.5},

    {0,-0.5, 0, 0.5, 0, 0, 0, 0.5,-0.5},

    {0, 0,-0.5, 0, 0, 0.5, 0.5, 0,-0.5},

    {0, 0, 0, 0, 0.25, 0, 0.25, 0.25, 0.25},

    {0, 0, 0, 0,-0.25, 0, 0.25,-0.25, 0.25},

    {0, 0, 0, 0, 0.25, 0,-0.25,-0.25, 0.25},

    {0, 0, 0, 0,-0.25, 0,-0.25, 0.25, 0.25}
};


// ============================================================
// 9. 测试步数
// ============================================================

const int nSteps = 10;


// ============================================================
// 10. 全局数组
// ============================================================

// -------- phase / density --------

double phi[NY][NX];
double rho[NY][NX];

double gradPhiX[NY][NX];
double gradPhiY[NY][NX];

double chemicalPotential[NY][NX];

double gradRhoX[NY][NX];
double gradRhoY[NY][NX];


// -------- flow macros --------

double ux[NY][NX];
double uy[NY][NX];

double pressure[NY][NX];

double dudx[NY][NX];
double dudy[NY][NX];

double dvdx[NY][NX];
double dvdy[NY][NX];


// -------- force --------

double forceX[NY][NX];
double forceY[NY][NX];


// -------- D2Q9 distributions --------

double g1[NY][NX][Q9];
double g2[NY][NX][Q9];

double g1Post[NY][NX][Q9];
double g2Post[NY][NX][Q9];

double g1Next[NY][NX][Q9];
double g2Next[NY][NX][Q9];


// -------- previous Q tensor --------

double oldQxx[NY][NX];
double oldQxy[NY][NX];
double oldQyy[NY][NX];


// -------- initial total momentum --------

double initialMomentumX = 0.0;
double initialMomentumY = 0.0;


// ============================================================
// 11. 函数声明
// ============================================================

void checkM9();

void initializePhi();

void computeDensity();

void computeGradPhi();

void computeChemicalPotential();

void computeSurfaceForce();

void computeGradRho();

void initializeFlowFields();

void initializeFlowDistribution();

void recoverVelocity();

void computeVelocityGradients();

void recoverPressure();

void initializeOldFlowTerms();

void collideFlow();

void updateOldFlowTerms();

void streamFlow();

bool checkNaN(int step);

void printDiagnostics(int step);


// ============================================================
// 12. main
// ============================================================

int main()
{
    cout
        << scientific
        << setprecision(12);

    cout
        << "============================================\n";

    cout
        << "Zhan 2023 MDF-LBM\n";

    cout
        << "Square droplet - FLOW WITH-FORCE test\n";

    cout
        << "============================================\n\n";


    cout << "NX                 = " << NX << "\n";
    cout << "NY                 = " << NY << "\n";
    cout << "dx                 = " << dx << "\n";
    cout << "dt                 = " << dt << "\n";
    cout << "c                  = " << c << "\n";
    cout << "cs^2               = " << cs2 << "\n";
    cout << "nu                 = " << nu << "\n";
    cout << "su1                = " << su1 << "\n";
    cout << "betaPF             = " << betaPF << "\n";
    cout << "kappaPF            = " << kappaPF << "\n\n";


    // --------------------------------------------------------
    // Matrix check
    // --------------------------------------------------------

    checkM9();


    // --------------------------------------------------------
    // Frozen phase/density initialization
    // --------------------------------------------------------

    initializePhi();

    computeDensity();

    computeGradPhi();

    computeChemicalPotential();


    // --------------------------------------------------------
    // Initial flow = zero
    // --------------------------------------------------------

    initializeFlowFields();


    // --------------------------------------------------------
    // Surface force F = mu grad(phi)
    // --------------------------------------------------------

    computeSurfaceForce();


    // --------------------------------------------------------
    // Frozen density gradient
    // --------------------------------------------------------

    computeGradRho();


    // --------------------------------------------------------
    // Initialize transformed g
    //
    // mg = meq - dt*mF/2
    // --------------------------------------------------------

    initializeFlowDistribution();


    // --------------------------------------------------------
    // Recover Step-0 macros
    // --------------------------------------------------------

    recoverVelocity();

    computeVelocityGradients();

    recoverPressure();


    // --------------------------------------------------------
    // Initialize old Q from Step-0 state
    // --------------------------------------------------------

    initializeOldFlowTerms();


    // --------------------------------------------------------
    // Initial momentum
    // --------------------------------------------------------

    initialMomentumX = 0.0;
    initialMomentumY = 0.0;

    for (int y = 0; y < NY; y++)
    {
        for (int x = 0; x < NX; x++)
        {
            initialMomentumX +=
                rho[y][x] * ux[y][x];

            initialMomentumY +=
                rho[y][x] * uy[y][x];
        }
    }


    printDiagnostics(0);


    // ========================================================
    // 10-step flow-with-force test
    // ========================================================

    for (int step = 1; step <= nSteps; step++)
    {
        // current macros
        recoverVelocity();

        computeGradRho();

        computeVelocityGradients();

        recoverPressure();


        // collision at time n
        collideFlow();


        // save Q^n
        updateOldFlowTerms();


        // stream to n+1
        streamFlow();


        // recover n+1 for diagnostics
        recoverVelocity();

        computeGradRho();

        computeVelocityGradients();

        recoverPressure();


        if (checkNaN(step))
        {
            cout
                << "\nSimulation stopped.\n";

            return 1;
        }


        printDiagnostics(step);
    }


    cout
        << "\n============================================\n";

    cout
        << "Flow with-force 10-step test finished.\n";

    cout
        << "============================================\n";


    return 0;
}


// ============================================================
// 13. Check M9 * Minv9
// ============================================================

void checkM9()
{
    double maxError = 0.0;

    for (int i = 0; i < Q9; i++)
    {
        for (int j = 0; j < Q9; j++)
        {
            double value = 0.0;

            for (int k = 0; k < Q9; k++)
            {
                value +=
                    M9[i][k]
                    *
                    Minv9[k][j];
            }

            const double exact =
                (i == j)
                ?
                1.0
                :
                0.0;

            maxError =
                max(
                    maxError,
                    fabs(value - exact)
                );
        }
    }

    cout
        << "Max |M9*Minv9-I|  = "
        << maxError
        << "\n\n";
}


// ============================================================
// 14. Initialize square droplet
// ============================================================

void initializePhi()
{
    const double a =
        0.5
        *
        (
            W
            -
            squareLength
        );

    const double b =
        0.5
        *
        (
            squareLength
            -
            W
        );


    for (int j = 0; j < NY; j++)
    {
        for (int i = 0; i < NX; i++)
        {
            const double x =
                -1.0
                +
                (i + 0.5)
                *
                dx;

            const double y =
                -1.0
                +
                (j + 0.5)
                *
                dx;


            double d = 0.0;


            // bottom-left corner
            if (x <= a && y <= a)
            {
                const double r =
                    sqrt(
                        (x-a)*(x-a)
                        +
                        (y-a)*(y-a)
                    );

                d =
                    0.5 * W
                    -
                    r;
            }

            // bottom-right corner
            else if (x >= b && y <= a)
            {
                const double r =
                    sqrt(
                        (x-b)*(x-b)
                        +
                        (y-a)*(y-a)
                    );

                d =
                    0.5 * W
                    -
                    r;
            }

            // top-left corner
            else if (x <= a && y >= b)
            {
                const double r =
                    sqrt(
                        (x-a)*(x-a)
                        +
                        (y-b)*(y-b)
                    );

                d =
                    0.5 * W
                    -
                    r;
            }

            // top-right corner
            else if (x >= b && y >= b)
            {
                const double r =
                    sqrt(
                        (x-b)*(x-b)
                        +
                        (y-b)*(y-b)
                    );

                d =
                    0.5 * W
                    -
                    r;
            }

            // bottom edge
            else if (x > a && x < b && y <= a)
            {
                d =
                    0.5 * squareLength
                    +
                    y;
            }

            // top edge
            else if (x > a && x < b && y >= b)
            {
                d =
                    0.5 * squareLength
                    -
                    y;
            }

            // left edge
            else if (x <= a && y > a && y < b)
            {
                d =
                    0.5 * squareLength
                    +
                    x;
            }

            // right edge
            else if (x >= b && y > a && y < b)
            {
                d =
                    0.5 * squareLength
                    -
                    x;
            }

            // center region
            else
            {
                phi[j][i] =
                    phiA;

                continue;
            }


            phi[j][i]
                =
            0.5
            *
            (
                phiA
                +
                phiB
            )

            +

            0.5
            *
            (
                phiA
                -
                phiB
            )
            *
            tanh(
                2.0 * d / W
            );
        }
    }
}


// ============================================================
// 15. rho(phi)
// ============================================================

void computeDensity()
{
    for (int y = 0; y < NY; y++)
    {
        for (int x = 0; x < NX; x++)
        {
            rho[y][x]
                =
            (
                rhoA
                -
                rhoB
            )
            /
            (
                phiA
                -
                phiB
            )
            *
            (
                phi[y][x]
                -
                phiB
            )

            +

            rhoB;
        }
    }
}


// ============================================================
// 16. grad(phi)
// ============================================================

void computeGradPhi()
{
    for (int y = 0; y < NY; y++)
    {
        for (int x = 0; x < NX; x++)
        {
            double gx = 0.0;
            double gy = 0.0;

            for (int q = 1; q < Q5; q++)
            {
                const int xn =
                    (
                        x
                        +
                        cx5[q]
                        +
                        NX
                    )
                    %
                    NX;

                const int yn =
                    (
                        y
                        +
                        cy5[q]
                        +
                        NY
                    )
                    %
                    NY;

                const double cix =
                    c * cx5[q];

                const double ciy =
                    c * cy5[q];


                gx +=
                    w5[q]
                    *
                    cix
                    *
                    phi[yn][xn];

                gy +=
                    w5[q]
                    *
                    ciy
                    *
                    phi[yn][xn];
            }


            gradPhiX[y][x]
                =
            gx
            /
            (
                cs2
                *
                dt
            );


            gradPhiY[y][x]
                =
            gy
            /
            (
                cs2
                *
                dt
            );
        }
    }
}


// ============================================================
// 17. chemical potential
// ============================================================

void computeChemicalPotential()
{
    for (int y = 0; y < NY; y++)
    {
        for (int x = 0; x < NX; x++)
        {
            double lapPhi = 0.0;


            for (int q = 1; q < Q5; q++)
            {
                const int xn =
                    (
                        x
                        +
                        cx5[q]
                        +
                        NX
                    )
                    %
                    NX;

                const int yn =
                    (
                        y
                        +
                        cy5[q]
                        +
                        NY
                    )
                    %
                    NY;


                lapPhi +=
                    2.0
                    *
                    w5[q]
                    *
                    (
                        phi[yn][xn]
                        -
                        phi[y][x]
                    );
            }


            lapPhi /=
                (
                    cs2
                    *
                    dt
                    *
                    dt
                );


            const double ph =
                phi[y][x];


            chemicalPotential[y][x]
                =
            4.0
            *
            betaPF
            *
            (
                ph
                -
                phiA
            )
            *
            (
                ph
                -
                phiB
            )
            *
            (
                ph
                -
                0.5
                *
                (
                    phiA
                    +
                    phiB
                )
            )

            -

            kappaPF
            *
            lapPhi;
        }
    }
}


// ============================================================
// 18. F = mu grad(phi)
// ============================================================

void computeSurfaceForce()
{
    for (int y = 0; y < NY; y++)
    {
        for (int x = 0; x < NX; x++)
        {
            forceX[y][x]
                =
            chemicalPotential[y][x]
            *
            gradPhiX[y][x];


            forceY[y][x]
                =
            chemicalPotential[y][x]
            *
            gradPhiY[y][x];
        }
    }
}


// ============================================================
// 19. grad(rho)
// ============================================================

void computeGradRho()
{
    for (int y = 0; y < NY; y++)
    {
        for (int x = 0; x < NX; x++)
        {
            double gx = 0.0;
            double gy = 0.0;


            for (int q = 1; q < Q9; q++)
            {
                const int xn =
                    (
                        x
                        +
                        cx9[q]
                        +
                        NX
                    )
                    %
                    NX;

                const int yn =
                    (
                        y
                        +
                        cy9[q]
                        +
                        NY
                    )
                    %
                    NY;

                const double cix =
                    c * cx9[q];

                const double ciy =
                    c * cy9[q];


                gx +=
                    w9[q]
                    *
                    cix
                    *
                    rho[yn][xn];


                gy +=
                    w9[q]
                    *
                    ciy
                    *
                    rho[yn][xn];
            }


            gradRhoX[y][x]
                =
            gx
            /
            (
                cs2
                *
                dt
            );


            gradRhoY[y][x]
                =
            gy
            /
            (
                cs2
                *
                dt
            );
        }
    }
}


// ============================================================
// 20. initialize u,v,p
// ============================================================

void initializeFlowFields()
{
    for (int y = 0; y < NY; y++)
    {
        for (int x = 0; x < NX; x++)
        {
            ux[y][x]       = 0.0;
            uy[y][x]       = 0.0;
            pressure[y][x] = 0.0;
        }
    }
}


// ============================================================
// 21. initialize transformed g
//
// mg = meq - dt*mF/2
// ============================================================

void initializeFlowDistribution()
{
    for (int y = 0; y < NY; y++)
    {
        for (int x = 0; x < NX; x++)
        {
            const double r =
                rho[y][x];

            const double u =
                ux[y][x];

            const double v =
                uy[y][x];

            const double p =
                pressure[y][x];

            const double Fx =
                forceX[y][x];

            const double Fy =
                forceY[y][x];


            double meq1[Q9];
            double meq2[Q9];

            double mF1[Q9];
            double mF2[Q9];

            double mg1[Q9];
            double mg2[Q9];


            // x momentum equilibrium

            meq1[0] = r * u;
            meq1[1] = r * u * u + p;
            meq1[2] = r * u * v;
            meq1[3] = r * u / 3.0;
            meq1[4] = 0.0;
            meq1[5] = r * u / 3.0;
            meq1[6] = r * u * v / (3.0 * c);
            meq1[7] = (r * u * u + p) / (3.0 * c);
            meq1[8] = r * u / 9.0;


            // y momentum equilibrium

            meq2[0] = r * v;
            meq2[1] = r * u * v;
            meq2[2] = r * v * v + p;
            meq2[3] = r * v / 3.0;
            meq2[4] = 0.0;
            meq2[5] = r * v / 3.0;
            meq2[6] = (r * v * v + p) / (3.0 * c);
            meq2[7] = r * u * v / (3.0 * c);
            meq2[8] = r * v / 9.0;


            // force moments

            mF1[0] = Fx;
            mF1[1] = 0.0;
            mF1[2] = 0.0;
            mF1[3] = Fx / 3.0;
            mF1[4] = 0.0;
            mF1[5] = Fx / 3.0;
            mF1[6] = 0.0;
            mF1[7] = 0.0;
            mF1[8] = Fx / 9.0;


            mF2[0] = Fy;
            mF2[1] = 0.0;
            mF2[2] = 0.0;
            mF2[3] = Fy / 3.0;
            mF2[4] = 0.0;
            mF2[5] = Fy / 3.0;
            mF2[6] = 0.0;
            mF2[7] = 0.0;
            mF2[8] = Fy / 9.0;


            // transformed moments

            for (int k = 0; k < Q9; k++)
            {
                mg1[k]
                    =
                meq1[k]
                -
                0.5
                *
                dt
                *
                mF1[k];


                mg2[k]
                    =
                meq2[k]
                -
                0.5
                *
                dt
                *
                mF2[k];
            }


            // g = Minv * mg

            for (int q = 0; q < Q9; q++)
            {
                g1[y][x][q] = 0.0;
                g2[y][x][q] = 0.0;

                for (int k = 0; k < Q9; k++)
                {
                    g1[y][x][q]
                        +=
                    Minv9[q][k]
                    *
                    mg1[k];


                    g2[y][x][q]
                        +=
                    Minv9[q][k]
                    *
                    mg2[k];
                }
            }
        }
    }
}


// ============================================================
// 22. recover velocity
// ============================================================

void recoverVelocity()
{
    for (int y = 0; y < NY; y++)
    {
        for (int x = 0; x < NX; x++)
        {
            double sumG1 = 0.0;
            double sumG2 = 0.0;


            for (int q = 0; q < Q9; q++)
            {
                sumG1 +=
                    g1[y][x][q];

                sumG2 +=
                    g2[y][x][q];
            }


            ux[y][x]
                =
            (
                sumG1
                +
                0.5
                *
                dt
                *
                forceX[y][x]
            )
            /
            rho[y][x];


            uy[y][x]
                =
            (
                sumG2
                +
                0.5
                *
                dt
                *
                forceY[y][x]
            )
            /
            rho[y][x];
        }
    }
}


// ============================================================
// 23. velocity gradients
// ============================================================

void computeVelocityGradients()
{
    for (int y = 0; y < NY; y++)
    {
        for (int x = 0; x < NX; x++)
        {
            double ux_x = 0.0;
            double ux_y = 0.0;

            double uy_x = 0.0;
            double uy_y = 0.0;


            for (int q = 1; q < Q9; q++)
            {
                const int xn =
                    (
                        x
                        +
                        cx9[q]
                        +
                        NX
                    )
                    %
                    NX;

                const int yn =
                    (
                        y
                        +
                        cy9[q]
                        +
                        NY
                    )
                    %
                    NY;

                const double cix =
                    c * cx9[q];

                const double ciy =
                    c * cy9[q];


                ux_x +=
                    w9[q]
                    *
                    cix
                    *
                    ux[yn][xn];

                ux_y +=
                    w9[q]
                    *
                    ciy
                    *
                    ux[yn][xn];

                uy_x +=
                    w9[q]
                    *
                    cix
                    *
                    uy[yn][xn];

                uy_y +=
                    w9[q]
                    *
                    ciy
                    *
                    uy[yn][xn];
            }


            dudx[y][x]
                =
            ux_x / (cs2 * dt);


            dudy[y][x]
                =
            ux_y / (cs2 * dt);


            dvdx[y][x]
                =
            uy_x / (cs2 * dt);


            dvdy[y][x]
                =
            uy_y / (cs2 * dt);
        }
    }
}


// ============================================================
// 24. pressure recovery
//
// simplified Eq.32
// ============================================================

void recoverPressure()
{
    for (int y = 0; y < NY; y++)
    {
        for (int x = 0; x < NX; x++)
        {
            double firstMomentX = 0.0;
            double firstMomentY = 0.0;


            for (int q = 0; q < Q9; q++)
            {
                const double cix =
                    c * cx9[q];

                const double ciy =
                    c * cy9[q];


                firstMomentX +=
                    cix
                    *
                    g1[y][x][q];


                firstMomentY +=
                    ciy
                    *
                    g2[y][x][q];
            }


            pressure[y][x]
                =
            0.5
            *
            (
                firstMomentX
                +
                firstMomentY

                -

                rho[y][x]
                *
                (
                    ux[y][x] * ux[y][x]
                    +
                    uy[y][x] * uy[y][x]
                )

                +

                0.5
                *
                dt
                *
                cs2
                *
                (
                    ux[y][x]
                    *
                    gradRhoX[y][x]

                    +

                    uy[y][x]
                    *
                    gradRhoY[y][x]
                )
            );
        }
    }
}


// ============================================================
// 25. initialize old Q
// ============================================================

void initializeOldFlowTerms()
{
    for (int y = 0; y < NY; y++)
    {
        for (int x = 0; x < NX; x++)
        {
            oldQxx[y][x]
                =
            rho[y][x]
            *
            ux[y][x]
            *
            ux[y][x]
            +
            pressure[y][x];


            oldQxy[y][x]
                =
            rho[y][x]
            *
            ux[y][x]
            *
            uy[y][x];


            oldQyy[y][x]
                =
            rho[y][x]
            *
            uy[y][x]
            *
            uy[y][x]
            +
            pressure[y][x];
        }
    }
}


// ============================================================
// 26. collide flow
// ============================================================

void collideFlow()
{
    for (int y = 0; y < NY; y++)
    {
        for (int x = 0; x < NX; x++)
        {
            const double r =
                rho[y][x];

            const double u =
                ux[y][x];

            const double v =
                uy[y][x];

            const double p =
                pressure[y][x];

            const double Fx =
                forceX[y][x];

            const double Fy =
                forceY[y][x];


            const double Qxx =
                r * u * u
                +
                p;

            const double Qxy =
                r * u * v;

            const double Qyy =
                r * v * v
                +
                p;


            const double dQxxdt =
                (
                    Qxx
                    -
                    oldQxx[y][x]
                )
                /
                dt;


            const double dQxydt =
                (
                    Qxy
                    -
                    oldQxy[y][x]
                )
                /
                dt;


            const double dQyydt =
                (
                    Qyy
                    -
                    oldQyy[y][x]
                )
                /
                dt;


            double m1[Q9];
            double m2[Q9];

            double meq1[Q9];
            double meq2[Q9];

            double mF1[Q9];
            double mF2[Q9];

            double mG1[Q9];
            double mG2[Q9];

            double mPost1[Q9];
            double mPost2[Q9];


            // ------------------------------------------------
            // m = M g
            // ------------------------------------------------

            for (int k = 0; k < Q9; k++)
            {
                m1[k] = 0.0;
                m2[k] = 0.0;

                for (int q = 0; q < Q9; q++)
                {
                    m1[k]
                        +=
                    M9[k][q]
                    *
                    g1[y][x][q];


                    m2[k]
                        +=
                    M9[k][q]
                    *
                    g2[y][x][q];
                }
            }


            // ------------------------------------------------
            // equilibrium: x momentum
            // ------------------------------------------------

            meq1[0] = r * u;
            meq1[1] = Qxx;
            meq1[2] = Qxy;
            meq1[3] = r * u / 3.0;
            meq1[4] = 0.0;
            meq1[5] = r * u / 3.0;
            meq1[6] = Qxy / (3.0 * c);
            meq1[7] = Qxx / (3.0 * c);
            meq1[8] = r * u / 9.0;


            // ------------------------------------------------
            // equilibrium: y momentum
            // ------------------------------------------------

            meq2[0] = r * v;
            meq2[1] = Qxy;
            meq2[2] = Qyy;
            meq2[3] = r * v / 3.0;
            meq2[4] = 0.0;
            meq2[5] = r * v / 3.0;
            meq2[6] = Qyy / (3.0 * c);
            meq2[7] = Qxy / (3.0 * c);
            meq2[8] = r * v / 9.0;


            // ------------------------------------------------
            // force moments
            // ------------------------------------------------

            mF1[0] = Fx;
            mF1[1] = 0.0;
            mF1[2] = 0.0;
            mF1[3] = Fx / 3.0;
            mF1[4] = 0.0;
            mF1[5] = Fx / 3.0;
            mF1[6] = 0.0;
            mF1[7] = 0.0;
            mF1[8] = Fx / 9.0;


            mF2[0] = Fy;
            mF2[1] = 0.0;
            mF2[2] = 0.0;
            mF2[3] = Fy / 3.0;
            mF2[4] = 0.0;
            mF2[5] = Fy / 3.0;
            mF2[6] = 0.0;
            mF2[7] = 0.0;
            mF2[8] = Fy / 9.0;


            // ------------------------------------------------
            // initialize mG = 0
            // ------------------------------------------------

            for (int k = 0; k < Q9; k++)
            {
                mG1[k] = 0.0;
                mG2[k] = 0.0;
            }

            
            // ------------------------------------------------
            // mG for x momentum
            // ------------------------------------------------

            mG1[1]
                =
            dQxxdt / c

            +

            c / 3.0
            *
            (
                u
                *
                gradRhoX[y][x]

                -

                r
                *
                dudx[y][x]
            );


            mG1[2]
                =
            dQxydt / c

            +

            c / 3.0
            *
            (
                u
                *
                gradRhoY[y][x]

                -

                r
                *
                dvdx[y][x]
            );


            mG1[6] =
                mG1[2] / 3.0;

            mG1[7] =
                mG1[1] / 3.0;


            // ------------------------------------------------
            // mG for y momentum
            // ------------------------------------------------

            mG2[1]
                =
            dQxydt / c

            +

            c / 3.0
            *
            (
                v
                *
                gradRhoX[y][x]

                -

                r
                *
                dudy[y][x]
            );


            mG2[2]
                =
            dQyydt / c

            +

            c / 3.0
            *
            (
                v
                *
                gradRhoY[y][x]

                -

                r
                *
                dvdy[y][x]
            );


            mG2[6] =
                mG2[2] / 3.0;

            mG2[7] =
                mG2[1] / 3.0;

            
            
            
            // ------------------------------------------------
            // MRT collision
            //
            // mPost =
            // (I-S)m
            // + S meq
            // + dt (I-S/2)(mF+mG)
            // ------------------------------------------------

            for (int k = 0; k < Q9; k++)
            {
                mPost1[k]
                    =
                (1.0 - Su[k])
                *
                m1[k]

                +

                Su[k]
                *
                meq1[k]

                +

                dt
                *
                (
                    1.0
                    -
                    0.5
                    *
                    Su[k]
                )
                *
                (
                    mF1[k]
                    +
                    mG1[k]
                );


                mPost2[k]
                    =
                (1.0 - Su[k])
                *
                m2[k]

                +

                Su[k]
                *
                meq2[k]

                +

                dt
                *
                (
                    1.0
                    -
                    0.5
                    *
                    Su[k]
                )
                *
                (
                    mF2[k]
                    +
                    mG2[k]
                );
            }


            // ------------------------------------------------
            // gPost = Minv mPost
            // ------------------------------------------------

            for (int q = 0; q < Q9; q++)
            {
                g1Post[y][x][q] = 0.0;
                g2Post[y][x][q] = 0.0;

                for (int k = 0; k < Q9; k++)
                {
                    g1Post[y][x][q]
                        +=
                    Minv9[q][k]
                    *
                    mPost1[k];


                    g2Post[y][x][q]
                        +=
                    Minv9[q][k]
                    *
                    mPost2[k];
                }
            }
        }
    }
}


// ============================================================
// 27. save Q^n
// ============================================================

void updateOldFlowTerms()
{
    for (int y = 0; y < NY; y++)
    {
        for (int x = 0; x < NX; x++)
        {
            oldQxx[y][x]
                =
            rho[y][x]
            *
            ux[y][x]
            *
            ux[y][x]
            +
            pressure[y][x];


            oldQxy[y][x]
                =
            rho[y][x]
            *
            ux[y][x]
            *
            uy[y][x];


            oldQyy[y][x]
                =
            rho[y][x]
            *
            uy[y][x]
            *
            uy[y][x]
            +
            pressure[y][x];
        }
    }
}


// ============================================================
// 28. periodic streaming
// ============================================================

void streamFlow()
{
    // clear

    for (int y = 0; y < NY; y++)
    {
        for (int x = 0; x < NX; x++)
        {
            for (int q = 0; q < Q9; q++)
            {
                g1Next[y][x][q] = 0.0;
                g2Next[y][x][q] = 0.0;
            }
        }
    }


    // push streaming

    for (int y = 0; y < NY; y++)
    {
        for (int x = 0; x < NX; x++)
        {
            for (int q = 0; q < Q9; q++)
            {
                const int xn =
                    (
                        x
                        +
                        cx9[q]
                        +
                        NX
                    )
                    %
                    NX;

                const int yn =
                    (
                        y
                        +
                        cy9[q]
                        +
                        NY
                    )
                    %
                    NY;


                g1Next[yn][xn][q]
                    =
                g1Post[y][x][q];


                g2Next[yn][xn][q]
                    =
                g2Post[y][x][q];
            }
        }
    }


    // copy back

    for (int y = 0; y < NY; y++)
    {
        for (int x = 0; x < NX; x++)
        {
            for (int q = 0; q < Q9; q++)
            {
                g1[y][x][q]
                    =
                g1Next[y][x][q];


                g2[y][x][q]
                    =
                g2Next[y][x][q];
            }
        }
    }
}


// ============================================================
// 29. NaN / Inf check
// ============================================================

bool checkNaN(int step)
{
    for (int y = 0; y < NY; y++)
    {
        for (int x = 0; x < NX; x++)
        {
            if (
                !isfinite(ux[y][x])
                ||
                !isfinite(uy[y][x])
                ||
                !isfinite(pressure[y][x])
            )
            {
                cout
                    << "\nNaN / Inf detected!\n";

                cout
                    << "step = "
                    << step
                    << "\n";

                cout
                    << "x = "
                    << x
                    << "\n";

                cout
                    << "y = "
                    << y
                    << "\n";

                cout
                    << "u = "
                    << ux[y][x]
                    << "\n";

                cout
                    << "v = "
                    << uy[y][x]
                    << "\n";

                cout
                    << "p = "
                    << pressure[y][x]
                    << "\n";

                return true;
            }
        }
    }

    return false;
}


// ============================================================
// 30. diagnostics
// ============================================================

void printDiagnostics(int step)
{
    double maxUx       = 0.0;
    double maxUy       = 0.0;
    double maxSpeed    = 0.0;
    double maxPressure = 0.0;
    double maxForce    = 0.0;
    double maxMu       = 0.0;

    double momentumX = 0.0;
    double momentumY = 0.0;

    double divergenceSquareSum = 0.0;


    for (int y = 0; y < NY; y++)
    {
        for (int x = 0; x < NX; x++)
        {
            maxUx =
                max(
                    maxUx,
                    fabs(ux[y][x])
                );


            maxUy =
                max(
                    maxUy,
                    fabs(uy[y][x])
                );


            const double speed =
                sqrt(
                    ux[y][x] * ux[y][x]
                    +
                    uy[y][x] * uy[y][x]
                );


            maxSpeed =
                max(
                    maxSpeed,
                    speed
                );


            maxPressure =
                max(
                    maxPressure,
                    fabs(pressure[y][x])
                );


            const double forceMagnitude =
                sqrt(
                    forceX[y][x] * forceX[y][x]
                    +
                    forceY[y][x] * forceY[y][x]
                );


            maxForce =
                max(
                    maxForce,
                    forceMagnitude
                );


            maxMu =
                max(
                    maxMu,
                    fabs(
                        chemicalPotential[y][x]
                    )
                );


            momentumX +=
                rho[y][x]
                *
                ux[y][x];


            momentumY +=
                rho[y][x]
                *
                uy[y][x];


            const double divU =
                dudx[y][x]
                +
                dvdy[y][x];


            divergenceSquareSum +=
                divU
                *
                divU;
        }
    }


    const double meanMomentumChangeX =
        (
            momentumX
            -
            initialMomentumX
        )
        /
        (
            NX
            *
            NY
        );


    const double meanMomentumChangeY =
        (
            momentumY
            -
            initialMomentumY
        )
        /
        (
            NX
            *
            NY
        );


    const double Drms =
        sqrt(
            divergenceSquareSum
            /
            (
                NX
                *
                NY
            )
        );


    cout
        << "\n--------------------------------------------\n";

    cout
        << "Step = "
        << step
        << "\n";

    cout
        << "max|mu|              = "
        << maxMu
        << "\n";

    cout
        << "max|F|               = "
        << maxForce
        << "\n";

    cout
        << "max|u|               = "
        << maxUx
        << "\n";

    cout
        << "max|v|               = "
        << maxUy
        << "\n";

    cout
        << "max speed            = "
        << maxSpeed
        << "\n";

    cout
        << "max|p|               = "
        << maxPressure
        << "\n";

    cout
        << "<Delta rho*u>        = "
        << meanMomentumChangeX
        << "\n";

    cout
        << "<Delta rho*v>        = "
        << meanMomentumChangeY
        << "\n";

    cout
        << "D_rms                = "
        << Drms
        << "\n";
}