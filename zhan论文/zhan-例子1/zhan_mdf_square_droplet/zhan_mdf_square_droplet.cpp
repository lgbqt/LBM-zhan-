#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <algorithm>

using namespace std;


// ============================================================
// 1. 基本参数
// ============================================================

// ------------------------------------------------------------
// Paper fact:
// 计算区域 [-1,1] x [-1,1]
// W = 0.04 = 4 dx
// c = dx/dt = 100
// ------------------------------------------------------------

const double Lx = 2.0;
const double Ly = 2.0;

const double W  = 0.04;

const double dx = 0.01;

const double c  = 100.0;

const double dt = dx / c;


// ------------------------------------------------------------
// Derived from paper:
// dx = 0.01
// dt = 1.0e-4
// NX = NY = 200
// ------------------------------------------------------------

const int NX =
    static_cast<int>(Lx / dx);

const int NY =
    static_cast<int>(Ly / dx);


// ============================================================
// 2. Lattice parameters
// ============================================================

const int Q5 = 5;
const int Q9 = 9;

const double cs2Phase =
    c * c / 3.0;

const double cs2Flow =
    c * c / 3.0;


// ============================================================
// 3. Physical parameters
// ============================================================

// Paper fact
const double mobility =
    0.1;

const double nuA =
    0.1;

const double nuB =
    0.1;

const double sigma =
    0.001;

const double radiusFinal =
    0.5;


// ------------------------------------------------------------
// Paper gives density ratio rhoA/rhoB = 1000 : 1.
//
// Absolute normalization below is implementation choice.
// ------------------------------------------------------------

const double rhoA =
    1000.0;

const double rhoB =
    1.0;


// ------------------------------------------------------------
// Order parameter convention.
//
// implementation choice consistent with previous CHE tests
// ------------------------------------------------------------

const double phiA =
    1.0;

const double phiB =
   -1.0;


// adjustable parameter
const double eta =
    1.0;


// ============================================================
// 4. Phase-field free-energy parameters
// ============================================================

// beta = 12 sigma / [W (phiA-phiB)^4]

const double betaPF =
    12.0
    *
    sigma
    /
    (
        W
        *
        pow(
            phiA - phiB,
            4.0
        )
    );


// kappa = (3/2) sigma W / (phiA-phiB)^2

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
// 5. Initial square-droplet length
// ============================================================

const double pi =
    3.14159265358979323846;


// Paper:
//
// D = sqrt[ pi R^2 + (4-pi) W^2 / 4 ]

const double squareLength =
    sqrt(
        pi
        *
        radiusFinal
        *
        radiusFinal
        +
        (4.0 - pi)
        *
        W
        *
        W
        /
        4.0
    );


// ============================================================
// 6. D2Q5 phase-field lattice
// ============================================================

const int cx5[Q5] =
{
     0,
     1,
     0,
    -1,
     0
};

const int cy5[Q5] =
{
     0,
     0,
     1,
     0,
    -1
};

const double w5[Q5] =
{
    1.0 / 3.0,
    1.0 / 6.0,
    1.0 / 6.0,
    1.0 / 6.0,
    1.0 / 6.0
};


// Zhan Appendix A, Eq. (A1)

const double M5[Q5][Q5] =
{
    {1, 1, 1, 1, 1},

    {0, 1, 0,-1, 0},

    {0, 0, 1, 0,-1},

    {0, 1, 0, 1, 0},

    {0, 0, 1, 0, 1}
};


// inverse of M5

const double Minv5[Q5][Q5] =
{
    {1, 0, 0,-1,-1},

    {0, 0.5, 0, 0.5, 0},

    {0, 0, 0.5, 0, 0.5},

    {0,-0.5, 0, 0.5, 0},

    {0, 0,-0.5, 0, 0.5}
};


// ============================================================
// 7. D2Q9 momentum lattice
// ============================================================

const int cx9[Q9] =
{
     0,
     1,
     0,
    -1,
     0,
     1,
    -1,
    -1,
     1
};

const int cy9[Q9] =
{
     0,
     0,
     1,
     0,
    -1,
     1,
     1,
    -1,
    -1
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
// 8. Zhan D2Q9 transformation matrix
//    Appendix A, Eq. (A4)
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


// inverse of Zhan M9

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
// 9. Relaxation parameters
// ============================================================

// ------------------------------------------------------------
// Phase-field mobility-related relaxation
//
// s1^0 = 1 / [ Mphi/(eta cs^2 dt) + 0.5 ]
// ------------------------------------------------------------

const double sPhase1 =
    1.0
    /
    (
        mobility
        /
        (
            eta
            *
            cs2Phase
            *
            dt
        )
        +
        0.5
    );

const double SPhase[Q5] =
{
    1.0,
    sPhase1,
    sPhase1,
    1.0,
    1.0
};


// ------------------------------------------------------------
// Flow viscosity-related relaxation
//
// su1 = 1 / [ nu/(cs^2 dt) + 0.5 ]
//
// Here nuA = nuB = 0.1,
// so su1 is constant.
// ------------------------------------------------------------

const double su1 =
    1.0
    /
    (
        nuA
        /
        (
            cs2Flow
            *
            dt
        )
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
// 10. Macroscopic arrays
// ============================================================

double phi[NY][NX];

double rho[NY][NX];

double ux[NY][NX];
double uy[NY][NX];

double pressure[NY][NX];

double chemicalPotential[NY][NX];


// ============================================================
// 11. Gradient arrays
// ============================================================

double gradPhiX[NY][NX];
double gradPhiY[NY][NX];

double gradRhoX[NY][NX];
double gradRhoY[NY][NX];

double dudx[NY][NX];
double dudy[NY][NX];

double dvdx[NY][NX];
double dvdy[NY][NX];


// ============================================================
// 12. Force arrays
// ============================================================

double forceX[NY][NX];
double forceY[NY][NX];


// ============================================================
// 13. Distribution functions
// ============================================================

// phase field D2Q5

double g0[NY][NX][Q5];

double g0Post[NY][NX][Q5];
double g0Next[NY][NX][Q5];


// x-momentum D2Q9

double g1[NY][NX][Q9];

double g1Post[NY][NX][Q9];
double g1Next[NY][NX][Q9];


// y-momentum D2Q9

double g2[NY][NX][Q9];

double g2Post[NY][NX][Q9];
double g2Next[NY][NX][Q9];


// ============================================================
// 14. Old-time arrays
// ============================================================

double oldPhiUx[NY][NX];
double oldPhiUy[NY][NX];

double oldQxx[NY][NX];
double oldQxy[NY][NX];
double oldQyy[NY][NX];


// ============================================================
// 15. Initial statistical data
// ============================================================

double rhoInitial[NY][NX];

double rhoUxInitial[NY][NX];
double rhoUyInitial[NY][NX];

double kineticInitial[NY][NX];


// ============================================================
// Function declarations
// ============================================================

void checkM5();
void checkM9();

void initializePhi();

void computeDensity();

void computeGradPhi();

void computeChemicalPotential();

void computeSurfaceForce();

void initializeMacroscopicFields();

void initializePhaseDistribution();

void initializeFlowDistribution();

void initializeOldTimeTerms();

void printStep0Diagnostics();

void savePhiCSV(
    const char* filename
);


// ============================================================
// main
// ============================================================

int main()
{
    cout << setprecision(12);
    cout << scientific;

    cout
        << "============================================\n";

    cout
        << "Zhan 2023 MDF-LBM square droplet\n";

    cout
        << "Initialization / Step-0 test\n";

    cout
        << "============================================\n\n";


    // --------------------------------------------------------
    // Print parameters
    // --------------------------------------------------------

    cout
        << "NX                 = "
        << NX
        << "\n";

    cout
        << "NY                 = "
        << NY
        << "\n";

    cout
        << "dx                 = "
        << dx
        << "\n";

    cout
        << "dt                 = "
        << dt
        << "\n";

    cout
        << "c                  = "
        << c
        << "\n";

    cout
        << "W                  = "
        << W
        << "\n";

    cout
        << "square length D    = "
        << squareLength
        << "\n";

    cout
        << "rhoA               = "
        << rhoA
        << "\n";

    cout
        << "rhoB               = "
        << rhoB
        << "\n";

    cout
        << "density ratio      = "
        << rhoA / rhoB
        << "\n";

    cout
        << "mobility           = "
        << mobility
        << "\n";

    cout
        << "nu                 = "
        << nuA
        << "\n";

    cout
        << "sigma              = "
        << sigma
        << "\n";

    cout
        << "betaPF             = "
        << betaPF
        << "\n";

    cout
        << "kappaPF            = "
        << kappaPF
        << "\n";

    cout
        << "sPhase1            = "
        << sPhase1
        << "\n";

    cout
        << "su1                = "
        << su1
        << "\n\n";


    // --------------------------------------------------------
    // Matrix checks
    // --------------------------------------------------------

    checkM5();

    checkM9();


    // --------------------------------------------------------
    // Initialization
    // --------------------------------------------------------

    initializePhi();

    initializeMacroscopicFields();

    computeDensity();

    computeGradPhi();

    computeChemicalPotential();

    computeSurfaceForce();

    initializePhaseDistribution();

    initializeFlowDistribution();

    initializeOldTimeTerms();


    // --------------------------------------------------------
    // Step 0 diagnostics
    // --------------------------------------------------------

    printStep0Diagnostics();

    savePhiCSV(
        "square_droplet_phi_initial.csv"
    );


    cout
        << "\nInitialization finished.\n";

    return 0;
}


// ============================================================
// Matrix check: D2Q5
// ============================================================

void checkM5()
{
    double maxError =
        0.0;

    for (int i = 0; i < Q5; i++)
    {
        for (int j = 0; j < Q5; j++)
        {
            double value =
                0.0;

            for (int k = 0; k < Q5; k++)
            {
                value +=
                    M5[i][k]
                    *
                    Minv5[k][j];
            }

            double exact =
                (i == j)
                ?
                1.0
                :
                0.0;

            maxError =
                max(
                    maxError,
                    fabs(
                        value - exact
                    )
                );
        }
    }

    cout
        << "Max |M5*Minv5-I|  = "
        << maxError
        << "\n";
}


// ============================================================
// Matrix check: D2Q9
// ============================================================

void checkM9()
{
    double maxError =
        0.0;

    for (int i = 0; i < Q9; i++)
    {
        for (int j = 0; j < Q9; j++)
        {
            double value =
                0.0;

            for (int k = 0; k < Q9; k++)
            {
                value +=
                    M9[i][k]
                    *
                    Minv9[k][j];
            }

            double exact =
                (i == j)
                ?
                1.0
                :
                0.0;

            maxError =
                max(
                    maxError,
                    fabs(
                        value - exact
                    )
                );
        }
    }

    cout
        << "Max |M9*Minv9-I|  = "
        << maxError
        << "\n\n";
}


// ============================================================
// Initialize phi
//
// Direct implementation of the piecewise structure in
// Zhan Eq. (46).
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
            // ------------------------------------------------
            // physical coordinates
            //
            // domain [-1,1] x [-1,1]
            // cell-centered coordinates
            // ------------------------------------------------

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


            double argument =
                1.0e30;


            // ------------------------------------------------
            // bottom-left rounded corner
            // ------------------------------------------------

            if (
                x <= a
                &&
                y <= a
            )
            {
                double r =
                    sqrt(
                        (x - a)
                        *
                        (x - a)
                        +
                        (y - a)
                        *
                        (y - a)
                    );

                argument =
                    0.5 * W
                    -
                    r;
            }

            // ------------------------------------------------
            // bottom-right rounded corner
            // ------------------------------------------------

            else if (
                x >= b
                &&
                y <= a
            )
            {
                double r =
                    sqrt(
                        (x - b)
                        *
                        (x - b)
                        +
                        (y - a)
                        *
                        (y - a)
                    );

                argument =
                    0.5 * W
                    -
                    r;
            }

            // ------------------------------------------------
            // top-left rounded corner
            // ------------------------------------------------

            else if (
                x <= a
                &&
                y >= b
            )
            {
                double r =
                    sqrt(
                        (x - a)
                        *
                        (x - a)
                        +
                        (y - b)
                        *
                        (y - b)
                    );

                argument =
                    0.5 * W
                    -
                    r;
            }

            // ------------------------------------------------
            // top-right rounded corner
            // ------------------------------------------------

            else if (
                x >= b
                &&
                y >= b
            )
            {
                double r =
                    sqrt(
                        (x - b)
                        *
                        (x - b)
                        +
                        (y - b)
                        *
                        (y - b)
                    );

                argument =
                    0.5 * W
                    -
                    r;
            }

            // ------------------------------------------------
            // bottom edge
            // ------------------------------------------------

            else if (
                x > a
                &&
                x < b
                &&
                y <= a
            )
            {
                argument =
                    0.5
                    *
                    squareLength
                    +
                    y;
            }

            // ------------------------------------------------
            // top edge
            // ------------------------------------------------

            else if (
                x > a
                &&
                x < b
                &&
                y >= b
            )
            {
                argument =
                    0.5
                    *
                    squareLength
                    -
                    y;
            }

            // ------------------------------------------------
            // left edge
            // ------------------------------------------------

            else if (
                x <= a
                &&
                y > a
                &&
                y < b
            )
            {
                argument =
                    0.5
                    *
                    squareLength
                    +
                    x;
            }

            // ------------------------------------------------
            // right edge
            // ------------------------------------------------

            else if (
                x >= b
                &&
                y > a
                &&
                y < b
            )
            {
                argument =
                    0.5
                    *
                    squareLength
                    -
                    x;
            }

            // ------------------------------------------------
            // deep inside droplet
            // ------------------------------------------------

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
                2.0
                *
                argument
                /
                W
            );
        }
    }
}


// ============================================================
// Initialize macroscopic flow fields
// ============================================================

void initializeMacroscopicFields()
{
    for (int j = 0; j < NY; j++)
    {
        for (int i = 0; i < NX; i++)
        {
            ux[j][i] =
                0.0;

            uy[j][i] =
                0.0;

            pressure[j][i] =
                0.0;
        }
    }
}


// ============================================================
// rho(phi)
// ============================================================

void computeDensity()
{
    for (int j = 0; j < NY; j++)
    {
        for (int i = 0; i < NX; i++)
        {
            rho[j][i]
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
                phi[j][i]
                -
                phiB
            )
            +
            rhoB;
        }
    }
}


// ============================================================
// grad(phi)
//
// second-order isotropic lattice formula
// periodic boundaries
// ============================================================

void computeGradPhi()
{
    for (int j = 0; j < NY; j++)
    {
        for (int i = 0; i < NX; i++)
        {
            double gx =
                0.0;

            double gy =
                0.0;

            for (int q = 1; q < Q5; q++)
            {
                int in =
                    (i + cx5[q] + NX)
                    %
                    NX;

                int jn =
                    (j + cy5[q] + NY)
                    %
                    NY;


                const double cix =
                    c
                    *
                    cx5[q];

                const double ciy =
                    c
                    *
                    cy5[q];


                gx +=
                    w5[q]
                    *
                    cix
                    *
                    phi[jn][in];

                gy +=
                    w5[q]
                    *
                    ciy
                    *
                    phi[jn][in];
            }

            gradPhiX[j][i]
                =
            gx
            /
            (
                cs2Phase
                *
                dt
            );

            gradPhiY[j][i]
                =
            gy
            /
            (
                cs2Phase
                *
                dt
            );
        }
    }
}


// ============================================================
// chemical potential
//
// mu = 4 beta (phi-phiA)(phi-phiB)
//      [phi-(phiA+phiB)/2]
//      - kappa lap(phi)
// ============================================================

void computeChemicalPotential()
{
    for (int j = 0; j < NY; j++)
    {
        for (int i = 0; i < NX; i++)
        {
            double lapPhi =
                0.0;

            for (int q = 1; q < Q5; q++)
            {
                int in =
                    (i + cx5[q] + NX)
                    %
                    NX;

                int jn =
                    (j + cy5[q] + NY)
                    %
                    NY;


                lapPhi
                    +=
                2.0
                *
                w5[q]
                *
                (
                    phi[jn][in]
                    -
                    phi[j][i]
                );
            }

            lapPhi
                /=
            (
                cs2Phase
                *
                dt
                *
                dt
            );


            const double ph =
                phi[j][i];


            chemicalPotential[j][i]
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
// surface force:
//
// Fs = mu * grad(phi)
// ============================================================

void computeSurfaceForce()
{
    for (int j = 0; j < NY; j++)
    {
        for (int i = 0; i < NX; i++)
        {
            forceX[j][i]
                =
            chemicalPotential[j][i]
            *
            gradPhiX[j][i];

            forceY[j][i]
                =
            chemicalPotential[j][i]
            *
            gradPhiY[j][i];
        }
    }
}


// ============================================================
// Initialize phase distribution g0
//
// Appendix A, Eq. (A2):
//
// m_eq =
// [ phi,
//   phi*u/c,
//   phi*v/c,
//   eta*mu/3,
//   eta*mu/3 ]
//
// S = 0, therefore g0 = f0 at initialization.
// ============================================================

void initializePhaseDistribution()
{
    for (int j = 0; j < NY; j++)
    {
        for (int i = 0; i < NX; i++)
        {
            double meq[Q5];


            meq[0] =
                phi[j][i];

            meq[1] =
                phi[j][i]
                *
                ux[j][i]
                /
                c;

            meq[2] =
                phi[j][i]
                *
                uy[j][i]
                /
                c;

            meq[3] =
                eta
                *
                chemicalPotential[j][i]
                /
                3.0;

            meq[4] =
                eta
                *
                chemicalPotential[j][i]
                /
                3.0;


            for (int q = 0; q < Q5; q++)
            {
                g0[j][i][q] =
                    0.0;

                for (int k = 0; k < Q5; k++)
                {
                    g0[j][i][q]
                        +=
                    Minv5[q][k]
                    *
                    meq[k];
                }
            }
        }
    }
}


// ============================================================
// Initialize flow distributions g1 / g2
//
// Important:
//
// g = f - dt F / 2
//
// We initialize f at equilibrium,
// therefore in moment space:
//
// m_g = m_eq - dt/2 * m_F
//
// This guarantees the macroscopic recovery
//
// u_alpha = [sum_i g_i,alpha + dt F_alpha/2]/rho
//
// gives u=v=0 initially.
// ============================================================

void initializeFlowDistribution()
{
    for (int j = 0; j < NY; j++)
    {
        for (int i = 0; i < NX; i++)
        {
            const double r =
                rho[j][i];

            const double u =
                ux[j][i];

            const double v =
                uy[j][i];

            const double p =
                pressure[j][i];

            const double Fx =
                forceX[j][i];

            const double Fy =
                forceY[j][i];


            double meq1[Q9];
            double meq2[Q9];

            double mF1[Q9];
            double mF2[Q9];

            double mg1[Q9];
            double mg2[Q9];


            // ------------------------------------------------
            // x momentum equilibrium moments
            // ------------------------------------------------

            meq1[0] =
                r * u;

            meq1[1] =
                r * u * u
                +
                p;

            meq1[2] =
                r * u * v;

            meq1[3] =
                r * u
                /
                3.0;

            meq1[4] =
                0.0;

            meq1[5] =
                r * u
                /
                3.0;

            meq1[6] =
                r * u * v
                /
                (3.0 * c);

            meq1[7] =
                (
                    r * u * u
                    +
                    p
                )
                /
                (3.0 * c);

            meq1[8] =
                r * u
                /
                9.0;


            // ------------------------------------------------
            // y momentum equilibrium moments
            // ------------------------------------------------

            meq2[0] =
                r * v;

            meq2[1] =
                r * u * v;

            meq2[2] =
                r * v * v
                +
                p;

            meq2[3] =
                r * v
                /
                3.0;

            meq2[4] =
                0.0;

            meq2[5] =
                r * v
                /
                3.0;

            meq2[6] =
                (
                    r * v * v
                    +
                    p
                )
                /
                (3.0 * c);

            meq2[7] =
                r * u * v
                /
                (3.0 * c);

            meq2[8] =
                r * v
                /
                9.0;


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
            // transformed variable:
            //
            // g = f - dt F / 2
            // ------------------------------------------------

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


            // ------------------------------------------------
            // back to velocity space
            // ------------------------------------------------

            for (int q = 0; q < Q9; q++)
            {
                g1[j][i][q] =
                    0.0;

                g2[j][i][q] =
                    0.0;


                for (int k = 0; k < Q9; k++)
                {
                    g1[j][i][q]
                        +=
                    Minv9[q][k]
                    *
                    mg1[k];

                    g2[j][i][q]
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
// Initialize old-time quantities
// ============================================================

void initializeOldTimeTerms()
{
    for (int j = 0; j < NY; j++)
    {
        for (int i = 0; i < NX; i++)
        {
            oldPhiUx[j][i]
                =
            phi[j][i]
            *
            ux[j][i];

            oldPhiUy[j][i]
                =
            phi[j][i]
            *
            uy[j][i];


            oldQxx[j][i]
                =
            rho[j][i]
            *
            ux[j][i]
            *
            ux[j][i]
            +
            pressure[j][i];

            oldQxy[j][i]
                =
            rho[j][i]
            *
            ux[j][i]
            *
            uy[j][i];

            oldQyy[j][i]
                =
            rho[j][i]
            *
            uy[j][i]
            *
            uy[j][i]
            +
            pressure[j][i];


            rhoInitial[j][i]
                =
            rho[j][i];

            rhoUxInitial[j][i]
                =
            rho[j][i]
            *
            ux[j][i];

            rhoUyInitial[j][i]
                =
            rho[j][i]
            *
            uy[j][i];

            kineticInitial[j][i]
                =
            0.5
            *
            (
                ux[j][i]
                *
                ux[j][i]
                +
                uy[j][i]
                *
                uy[j][i]
            );
        }
    }
}


// ============================================================
// Step-0 diagnostics
// ============================================================

void printStep0Diagnostics()
{
    double minPhi =
        phi[0][0];

    double maxPhi =
        phi[0][0];

    double minRho =
        rho[0][0];

    double maxRho =
        rho[0][0];

    double maxMu =
        0.0;

    double maxForce =
        0.0;

    double maxVelocityRecovered =
        0.0;


    for (int j = 0; j < NY; j++)
    {
        for (int i = 0; i < NX; i++)
        {
            minPhi =
                min(
                    minPhi,
                    phi[j][i]
                );

            maxPhi =
                max(
                    maxPhi,
                    phi[j][i]
                );

            minRho =
                min(
                    minRho,
                    rho[j][i]
                );

            maxRho =
                max(
                    maxRho,
                    rho[j][i]
                );

            maxMu =
                max(
                    maxMu,
                    fabs(
                        chemicalPotential[j][i]
                    )
                );


            double forceMag =
                sqrt(
                    forceX[j][i]
                    *
                    forceX[j][i]
                    +
                    forceY[j][i]
                    *
                    forceY[j][i]
                );

            maxForce =
                max(
                    maxForce,
                    forceMag
                );


            // ------------------------------------------------
            // Recover velocity from g1/g2
            //
            // u = [sum g + dt F/2] / rho
            // ------------------------------------------------

            double sumG1 =
                0.0;

            double sumG2 =
                0.0;

            for (int q = 0; q < Q9; q++)
            {
                sumG1 +=
                    g1[j][i][q];

                sumG2 +=
                    g2[j][i][q];
            }


            double uRecovered =
                (
                    sumG1
                    +
                    0.5
                    *
                    dt
                    *
                    forceX[j][i]
                )
                /
                rho[j][i];

            double vRecovered =
                (
                    sumG2
                    +
                    0.5
                    *
                    dt
                    *
                    forceY[j][i]
                )
                /
                rho[j][i];


            double speed =
                sqrt(
                    uRecovered
                    *
                    uRecovered
                    +
                    vRecovered
                    *
                    vRecovered
                );

            maxVelocityRecovered =
                max(
                    maxVelocityRecovered,
                    speed
                );
        }
    }


    cout
        << "\n============================================\n";

    cout
        << "Step 0 diagnostics\n";

    cout
        << "============================================\n";

    cout
        << "min(phi)           = "
        << minPhi
        << "\n";

    cout
        << "max(phi)           = "
        << maxPhi
        << "\n";

    cout
        << "min(rho)           = "
        << minRho
        << "\n";

    cout
        << "max(rho)           = "
        << maxRho
        << "\n";

    cout
        << "max|mu|            = "
        << maxMu
        << "\n";

    cout
        << "max|F_surface|     = "
        << maxForce
        << "\n";

    cout
        << "max recovered |u|  = "
        << maxVelocityRecovered
        << "\n";

    cout
        << "============================================\n";
}


// ============================================================
// Save phi
// ============================================================

void savePhiCSV(
    const char* filename
)
{
    ofstream fout(
        filename
    );

    for (int j = 0; j < NY; j++)
    {
        for (int i = 0; i < NX; i++)
        {
            fout
                << phi[j][i];

            if (i < NX - 1)
            {
                fout
                    << ",";
            }
        }

        fout
            << "\n";
    }

    fout.close();
}
