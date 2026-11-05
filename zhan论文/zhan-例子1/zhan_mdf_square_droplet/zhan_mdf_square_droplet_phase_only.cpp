#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <algorithm>

using namespace std;


// ============================================================
// 1. Basic parameters
// ============================================================
//
// Paper:
//
// Domain:
// [-1,1] x [-1,1]
//
// W = 0.04 = 4 dx
//
// c = dx / dt = 100
//
// M_phi = 0.1
//
// sigma = 0.001
//
// final radius R = 0.5
//
// This program is only the PHASE-ONLY TEST.
// Fluid velocity is fixed:
//
// u = v = 0
//
// ============================================================

const double Lx = 2.0;
const double Ly = 2.0;

const double W =
    0.04;

const double dx =
    0.01;

const double c =
    100.0;

const double dt =
    dx / c;


// ------------------------------------------------------------
// Derived
// ------------------------------------------------------------

const int NX =
    static_cast<int>(
        Lx / dx
    );

const int NY =
    static_cast<int>(
        Ly / dx
    );


// ============================================================
// 2. Phase-field lattice
// ============================================================

const int Q5 =
    5;

const double cs2 =
    c * c / 3.0;


// ============================================================
// 3. Physical parameters
// ============================================================

const double mobility =
    0.1;

const double sigma =
    0.001;

const double radiusFinal =
    0.5;


// ------------------------------------------------------------
// Implementation choice:
//
// same CHE convention used in our previous tests
// ------------------------------------------------------------

const double phiA =
    1.0;

const double phiB =
   -1.0;


// adjustable parameter eta

const double eta =
    1.0;


// ============================================================
// 4. Density parameters
// ============================================================
//
// Paper gives:
//
// rhoA : rhoB = 1000 : 1
//
// Here rhoB = 1 and rhoA = 1000
// are used as normalization.
// ============================================================

const double rhoA =
    1000.0;

const double rhoB =
    1.0;


// ============================================================
// 5. Free-energy parameters
// ============================================================
//
// beta = 12 sigma /
//        [ W (phiA-phiB)^4 ]
//
// kappa = 3 sigma W /
//         [ 2 (phiA-phiB)^2 ]
// ============================================================

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
// 6. Initial square length
// ============================================================
//
// D = sqrt[
//     pi R^2
//     + (4-pi)W^2/4
// ]
// ============================================================

const double PI =
    3.14159265358979323846;


const double squareLength =
    sqrt(
        PI
        *
        radiusFinal
        *
        radiusFinal

        +

        (4.0 - PI)
        *
        W
        *
        W
        /
        4.0
    );


// ============================================================
// 7. D2Q5 discrete velocities
// ============================================================

const int cx[Q5] =
{
     0,
     1,
     0,
    -1,
     0
};


const int cy[Q5] =
{
     0,
     0,
     1,
     0,
    -1
};


const double w[Q5] =
{
    1.0 / 3.0,

    1.0 / 6.0,
    1.0 / 6.0,
    1.0 / 6.0,
    1.0 / 6.0
};


// ============================================================
// 8. D2Q5 transformation matrix
//
// Zhan Appendix A, Eq. (A1)
// ============================================================

const double M5[Q5][Q5] =
{
    {1, 1, 1, 1, 1},

    {0, 1, 0,-1, 0},

    {0, 0, 1, 0,-1},

    {0, 1, 0, 1, 0},

    {0, 0, 1, 0, 1}
};


// inverse matrix

const double Minv5[Q5][Q5] =
{
    {1, 0, 0,-1,-1},

    {0, 0.5, 0, 0.5, 0},

    {0, 0, 0.5, 0, 0.5},

    {0,-0.5, 0, 0.5, 0},

    {0, 0,-0.5, 0, 0.5}
};


// ============================================================
// 9. Phase relaxation parameters
// ============================================================
//
// s1^0 =
//
// 1 / [
//     M_phi / (eta cs^2 dt)
//     + 0.5
// ]
//
// Appendix:
//
// S0 = diag(
//     s0,
//     s1,
//     s1,
//     s2,
//     s2
// )
//
// Other relaxation parameters = 1
// ============================================================

const double sPhase1 =
    1.0
    /
    (
        mobility
        /
        (
            eta
            *
            cs2
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


// ============================================================
// 10. Number of steps
// ============================================================

const int nSteps =
    10;


// ============================================================
// 11. Macroscopic arrays
// ============================================================

double phi[NY][NX];

double rho[NY][NX];

double chemicalPotential[NY][NX];

double gradPhiX[NY][NX];
double gradPhiY[NY][NX];


// ============================================================
// 12. Phase distribution arrays
// ============================================================

double g0[NY][NX][Q5];

double g0Post[NY][NX][Q5];

double g0Next[NY][NX][Q5];


// ============================================================
// 13. Previous-time arrays
//
// Even though u=v=0 in this test,
// we keep them because the full model later needs:
//
// d_t(phi u)
// ============================================================

double oldPhiUx[NY][NX];

double oldPhiUy[NY][NX];


// ============================================================
// 14. Initial diagnostics
// ============================================================

double initialPhiSum =
    0.0;

double initialRhoSum =
    0.0;


// ============================================================
// Function declarations
// ============================================================

void checkM5();

void initializePhi();

void computeDensity();

void computeGradPhi();

void computeChemicalPotential();

void initializePhaseDistribution();

void initializeOldPhaseTerms();

void recoverPhi();

void collidePhase();

void streamPhase();

void updateOldPhaseTerms();

bool checkNaN(
    int step
);

void printDiagnostics(
    int step
);

void savePhiCSV(
    const char* filename
);


// ============================================================
// main
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
        << "Square droplet - PHASE ONLY test\n";

    cout
        << "============================================\n\n";


    // --------------------------------------------------------
    // Basic parameters
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
        << "cs^2               = "
        << cs2
        << "\n";

    cout
        << "W                  = "
        << W
        << "\n";

    cout
        << "D                  = "
        << squareLength
        << "\n";

    cout
        << "mobility           = "
        << mobility
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
        << "\n\n";


    // --------------------------------------------------------
    // Matrix test
    // --------------------------------------------------------

    checkM5();


    // ========================================================
    // INITIALIZATION
    // ========================================================

    initializePhi();

    computeDensity();

    computeGradPhi();

    computeChemicalPotential();

    initializePhaseDistribution();

    initializeOldPhaseTerms();


    // --------------------------------------------------------
    // Store initial mass
    // --------------------------------------------------------

    initialPhiSum =
        0.0;

    initialRhoSum =
        0.0;


    for (int y = 0; y < NY; y++)
    {
        for (int x = 0; x < NX; x++)
        {
            initialPhiSum
                +=
            phi[y][x];

            initialRhoSum
                +=
            rho[y][x];
        }
    }


    // --------------------------------------------------------
    // Step 0
    // --------------------------------------------------------

    printDiagnostics(
        0
    );


    savePhiCSV(
        "phi_step_000000.csv"
    );


    // ========================================================
    // PHASE-ONLY TIME LOOP
    // ========================================================
    //
    // Important:
    //
    // u = v = 0 during this diagnostic test.
    //
    // Therefore:
    //
    // mG = 0 analytically,
    //
    // but the oldPhiUx / oldPhiUy structure is retained.
    //
    // ========================================================

    for (int step = 1;
         step <= nSteps;
         step++)
    {
        // ----------------------------------------------------
        // At the beginning of a time step:
        //
        // g0^n -> phi^n
        // ----------------------------------------------------

        recoverPhi();


        // ----------------------------------------------------
        // phi -> rho
        // ----------------------------------------------------

        computeDensity();


        // ----------------------------------------------------
        // spatial derivatives and chemical potential
        // ----------------------------------------------------

        computeGradPhi();

        computeChemicalPotential();


        // ----------------------------------------------------
        // collision:
        //
        // g0^n -> g0Post
        // ----------------------------------------------------

        collidePhase();


        // ----------------------------------------------------
        // Save current time-level phi*u
        //
        // Here u=v=0,
        // but this timing is the same structure
        // we will use in the complete coupled model.
        // ----------------------------------------------------

        updateOldPhaseTerms();


        // ----------------------------------------------------
        // streaming:
        //
        // g0Post -> g0Next -> g0
        // ----------------------------------------------------

        streamPhase();


        // ----------------------------------------------------
        // For diagnostics we recover phi^(n+1)
        // from the streamed distribution.
        // ----------------------------------------------------

        recoverPhi();

        computeDensity();

        computeGradPhi();

        computeChemicalPotential();


        // ----------------------------------------------------
        // NaN / Inf detection
        // ----------------------------------------------------

        if (
            checkNaN(
                step
            )
        )
        {
            cout
                << "\nSimulation stopped.\n";

            return 1;
        }


        // ----------------------------------------------------
        // Print every step in this 10-step test
        // ----------------------------------------------------

        printDiagnostics(
            step
        );
    }


    // --------------------------------------------------------
    // Save final phase field
    // --------------------------------------------------------

    savePhiCSV(
        "phi_step_000010.csv"
    );


    cout
        << "\n============================================\n";

    cout
        << "Phase-only 10-step test finished.\n";

    cout
        << "============================================\n";


    return 0;
}


// ============================================================
// Check M5 * Minv5 = I
// ============================================================

void checkM5()
{
    double maxError =
        0.0;


    for (int i = 0;
         i < Q5;
         i++)
    {
        for (int j = 0;
             j < Q5;
             j++)
        {
            double value =
                0.0;


            for (int k = 0;
                 k < Q5;
                 k++)
            {
                value
                    +=
                M5[i][k]
                *
                Minv5[k][j];
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

                    fabs(
                        value
                        -
                        exact
                    )
                );
        }
    }


    cout
        << "Max |M5*Minv5-I|  = "
        << maxError
        << "\n\n";
}


// ============================================================
// Initialize square droplet
//
// Based on Zhan Eq. (46)
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


    for (int yIndex = 0;
         yIndex < NY;
         yIndex++)
    {
        for (int xIndex = 0;
             xIndex < NX;
             xIndex++)
        {
            // ------------------------------------------------
            // cell-centered physical coordinates
            // ------------------------------------------------

            const double x =
                -1.0
                +
                (
                    xIndex
                    +
                    0.5
                )
                *
                dx;


            const double y =
                -1.0
                +
                (
                    yIndex
                    +
                    0.5
                )
                *
                dx;


            double distanceFunction =
                0.0;


            // ------------------------------------------------
            // Four rounded corners
            // ------------------------------------------------

            if (
                x <= a
                &&
                y <= a
            )
            {
                const double r =
                    sqrt(
                        (x - a)
                        *
                        (x - a)

                        +

                        (y - a)
                        *
                        (y - a)
                    );


                distanceFunction =
                    0.5 * W
                    -
                    r;
            }

            else if (
                x >= b
                &&
                y <= a
            )
            {
                const double r =
                    sqrt(
                        (x - b)
                        *
                        (x - b)

                        +

                        (y - a)
                        *
                        (y - a)
                    );


                distanceFunction =
                    0.5 * W
                    -
                    r;
            }

            else if (
                x <= a
                &&
                y >= b
            )
            {
                const double r =
                    sqrt(
                        (x - a)
                        *
                        (x - a)

                        +

                        (y - b)
                        *
                        (y - b)
                    );


                distanceFunction =
                    0.5 * W
                    -
                    r;
            }

            else if (
                x >= b
                &&
                y >= b
            )
            {
                const double r =
                    sqrt(
                        (x - b)
                        *
                        (x - b)

                        +

                        (y - b)
                        *
                        (y - b)
                    );


                distanceFunction =
                    0.5 * W
                    -
                    r;
            }


            // ------------------------------------------------
            // Bottom edge
            // ------------------------------------------------

            else if (
                x > a
                &&
                x < b
                &&
                y <= a
            )
            {
                distanceFunction =
                    0.5
                    *
                    squareLength
                    +
                    y;
            }


            // ------------------------------------------------
            // Top edge
            // ------------------------------------------------

            else if (
                x > a
                &&
                x < b
                &&
                y >= b
            )
            {
                distanceFunction =
                    0.5
                    *
                    squareLength
                    -
                    y;
            }


            // ------------------------------------------------
            // Left edge
            // ------------------------------------------------

            else if (
                x <= a
                &&
                y > a
                &&
                y < b
            )
            {
                distanceFunction =
                    0.5
                    *
                    squareLength
                    +
                    x;
            }


            // ------------------------------------------------
            // Right edge
            // ------------------------------------------------

            else if (
                x >= b
                &&
                y > a
                &&
                y < b
            )
            {
                distanceFunction =
                    0.5
                    *
                    squareLength
                    -
                    x;
            }


            // ------------------------------------------------
            // Deep inside droplet
            // ------------------------------------------------

            else
            {
                phi[yIndex][xIndex] =
                    phiA;

                continue;
            }


            // ------------------------------------------------
            // Diffuse-interface profile
            // ------------------------------------------------

            phi[yIndex][xIndex]
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
                distanceFunction
                /
                W
            );
        }
    }
}


// ============================================================
// rho(phi)
// ============================================================

void computeDensity()
{
    for (int y = 0;
         y < NY;
         y++)
    {
        for (int x = 0;
             x < NX;
             x++)
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
// Isotropic gradient of phi
//
// grad Phi =
//
// sum_i w_i c_i Phi(x+c_i dt)
//
// / (cs^2 dt)
//
// Periodic boundary
// ============================================================

void computeGradPhi()
{
    for (int y = 0;
         y < NY;
         y++)
    {
        for (int x = 0;
             x < NX;
             x++)
        {
            double gx =
                0.0;

            double gy =
                0.0;


            for (int q = 1;
                 q < Q5;
                 q++)
            {
                const int xn =
                    (
                        x
                        +
                        cx[q]
                        +
                        NX
                    )
                    %
                    NX;


                const int yn =
                    (
                        y
                        +
                        cy[q]
                        +
                        NY
                    )
                    %
                    NY;


                const double cix =
                    c
                    *
                    cx[q];


                const double ciy =
                    c
                    *
                    cy[q];


                gx
                    +=
                w[q]
                *
                cix
                *
                phi[yn][xn];


                gy
                    +=
                w[q]
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
// Chemical potential
//
// mu =
// 4 beta
// (phi-phiA)
// (phi-phiB)
// [phi-(phiA+phiB)/2]
//
// - kappa lap(phi)
// ============================================================

void computeChemicalPotential()
{
    for (int y = 0;
         y < NY;
         y++)
    {
        for (int x = 0;
             x < NX;
             x++)
        {
            double lapPhi =
                0.0;


            for (int q = 1;
                 q < Q5;
                 q++)
            {
                const int xn =
                    (
                        x
                        +
                        cx[q]
                        +
                        NX
                    )
                    %
                    NX;


                const int yn =
                    (
                        y
                        +
                        cy[q]
                        +
                        NY
                    )
                    %
                    NY;


                lapPhi
                    +=
                2.0
                *
                w[q]
                *
                (
                    phi[yn][xn]
                    -
                    phi[y][x]
                );
            }


            lapPhi
                /=
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
// Initialize phase distribution
//
// Appendix A:
//
// m_eq =
//
// [ phi,
//   phi*u/c,
//   phi*v/c,
//   eta*mu/3,
//   eta*mu/3 ]
//
// In this test:
//
// u = 0
// v = 0
//
// S = 0
//
// Therefore g0 = f_eq initially.
// ============================================================

void initializePhaseDistribution()
{
    for (int y = 0;
         y < NY;
         y++)
    {
        for (int x = 0;
             x < NX;
             x++)
        {
            double meq[Q5];


            meq[0] =
                phi[y][x];


            meq[1] =
                0.0;


            meq[2] =
                0.0;


            meq[3] =
                eta
                *
                chemicalPotential[y][x]
                /
                3.0;


            meq[4] =
                eta
                *
                chemicalPotential[y][x]
                /
                3.0;


            for (int q = 0;
                 q < Q5;
                 q++)
            {
                g0[y][x][q] =
                    0.0;


                for (int k = 0;
                     k < Q5;
                     k++)
                {
                    g0[y][x][q]
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
// Initialize old phi*u
//
// u=v=0
//
// so both arrays are zero.
// ============================================================

void initializeOldPhaseTerms()
{
    for (int y = 0;
         y < NY;
         y++)
    {
        for (int x = 0;
             x < NX;
             x++)
        {
            oldPhiUx[y][x] =
                0.0;


            oldPhiUy[y][x] =
                0.0;
        }
    }
}


// ============================================================
// Recover phi
//
// Eq. (35a):
//
// phi = sum_i g_i,0 + dt*S/2
//
// S = 0
//
// therefore:
//
// phi = sum_i g_i,0
// ============================================================

void recoverPhi()
{
    for (int y = 0;
         y < NY;
         y++)
    {
        for (int x = 0;
             x < NX;
             x++)
        {
            double sum =
                0.0;


            for (int q = 0;
                 q < Q5;
                 q++)
            {
                sum
                    +=
                g0[y][x][q];
            }


            phi[y][x] =
                sum;
        }
    }
}


// ============================================================
// Phase collision
//
// Zhan Eq. (36):
//
// mPost =
// (I-S)m
// + S meq
// + dt(I-S/2)(mF+mG)
//
// Immiscible:
// S(source) = 0
//
// therefore:
//
// mF = 0
//
// CHE Appendix A:
//
// meq =
// [phi,
//  phi*u/c,
//  phi*v/c,
//  eta*mu/3,
//  eta*mu/3]
//
// mG =
// [0,
//  d_t(phi*u)/c,
//  d_t(phi*v)/c,
//  0,
//  0]
//
// Phase-only:
//
// u=v=0
//
// ============================================================

void collidePhase()
{
    for (int y = 0;
         y < NY;
         y++)
    {
        for (int x = 0;
             x < NX;
             x++)
        {
            double m[Q5];

            double meq[Q5];

            double mG[Q5];

            double mPost[Q5];


            // ------------------------------------------------
            // m = M g
            // ------------------------------------------------

            for (int k = 0;
                 k < Q5;
                 k++)
            {
                m[k] =
                    0.0;


                for (int q = 0;
                     q < Q5;
                     q++)
                {
                    m[k]
                        +=
                    M5[k][q]
                    *
                    g0[y][x][q];
                }
            }


            // ------------------------------------------------
            // Equilibrium moments
            // ------------------------------------------------

            meq[0] =
                phi[y][x];


            meq[1] =
                0.0;


            meq[2] =
                0.0;


            meq[3] =
                eta
                *
                chemicalPotential[y][x]
                /
                3.0;


            meq[4] =
                eta
                *
                chemicalPotential[y][x]
                /
                3.0;


            // ------------------------------------------------
            // Time derivative:
            //
            // d_t(phi*u)
            //
            // Current u=v=0.
            // ------------------------------------------------

            const double currentPhiUx =
                0.0;


            const double currentPhiUy =
                0.0;


            const double dPhiUxDt =
                (
                    currentPhiUx
                    -
                    oldPhiUx[y][x]
                )
                /
                dt;


            const double dPhiUyDt =
                (
                    currentPhiUy
                    -
                    oldPhiUy[y][x]
                )
                /
                dt;


            mG[0] =
                0.0;


            mG[1] =
                dPhiUxDt
                /
                c;


            mG[2] =
                dPhiUyDt
                /
                c;


            mG[3] =
                0.0;


            mG[4] =
                0.0;


            // ------------------------------------------------
            // MRT collision
            // ------------------------------------------------

            for (int k = 0;
                 k < Q5;
                 k++)
            {
                mPost[k]
                    =
                (
                    1.0
                    -
                    SPhase[k]
                )
                *
                m[k]

                +

                SPhase[k]
                *
                meq[k]

                +

                dt
                *
                (
                    1.0
                    -
                    0.5
                    *
                    SPhase[k]
                )
                *
                mG[k];
            }


            // ------------------------------------------------
            // gPost = M^-1 mPost
            // ------------------------------------------------

            for (int q = 0;
                 q < Q5;
                 q++)
            {
                g0Post[y][x][q] =
                    0.0;


                for (int k = 0;
                     k < Q5;
                     k++)
                {
                    g0Post[y][x][q]
                        +=
                    Minv5[q][k]
                    *
                    mPost[k];
                }
            }
        }
    }
}


// ============================================================
// Update previous-time phi*u
//
// Current time level n.
//
// u=v=0 in this test.
// ============================================================

void updateOldPhaseTerms()
{
    for (int y = 0;
         y < NY;
         y++)
    {
        for (int x = 0;
             x < NX;
             x++)
        {
            oldPhiUx[y][x] =
                0.0;


            oldPhiUy[y][x] =
                0.0;
        }
    }
}


// ============================================================
// Periodic streaming
// ============================================================

void streamPhase()
{
    // --------------------------------------------------------
    // Clear next array
    // --------------------------------------------------------

    for (int y = 0;
         y < NY;
         y++)
    {
        for (int x = 0;
             x < NX;
             x++)
        {
            for (int q = 0;
                 q < Q5;
                 q++)
            {
                g0Next[y][x][q] =
                    0.0;
            }
        }
    }


    // --------------------------------------------------------
    // Push streaming
    // --------------------------------------------------------

    for (int y = 0;
         y < NY;
         y++)
    {
        for (int x = 0;
             x < NX;
             x++)
        {
            for (int q = 0;
                 q < Q5;
                 q++)
            {
                const int xn =
                    (
                        x
                        +
                        cx[q]
                        +
                        NX
                    )
                    %
                    NX;


                const int yn =
                    (
                        y
                        +
                        cy[q]
                        +
                        NY
                    )
                    %
                    NY;


                g0Next[yn][xn][q]
                    =
                g0Post[y][x][q];
            }
        }
    }


    // --------------------------------------------------------
    // next -> current
    // --------------------------------------------------------

    for (int y = 0;
         y < NY;
         y++)
    {
        for (int x = 0;
             x < NX;
             x++)
        {
            for (int q = 0;
                 q < Q5;
                 q++)
            {
                g0[y][x][q]
                    =
                g0Next[y][x][q];
            }
        }
    }
}


// ============================================================
// NaN / Inf check
// ============================================================

bool checkNaN(
    int step
)
{
    for (int y = 0;
         y < NY;
         y++)
    {
        for (int x = 0;
             x < NX;
             x++)
        {
            if (
                !isfinite(
                    phi[y][x]
                )
                ||
                !isfinite(
                    rho[y][x]
                )
                ||
                !isfinite(
                    chemicalPotential[y][x]
                )
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
                    << "phi = "
                    << phi[y][x]
                    << "\n";

                cout
                    << "rho = "
                    << rho[y][x]
                    << "\n";

                cout
                    << "mu = "
                    << chemicalPotential[y][x]
                    << "\n";


                return true;
            }
        }
    }


    return false;
}


// ============================================================
// Diagnostics
// ============================================================

void printDiagnostics(
    int step
)
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

    double phiSum =
        0.0;

    double rhoSum =
        0.0;


    for (int y = 0;
         y < NY;
         y++)
    {
        for (int x = 0;
             x < NX;
             x++)
        {
            minPhi =
                min(
                    minPhi,
                    phi[y][x]
                );


            maxPhi =
                max(
                    maxPhi,
                    phi[y][x]
                );


            minRho =
                min(
                    minRho,
                    rho[y][x]
                );


            maxRho =
                max(
                    maxRho,
                    rho[y][x]
                );


            maxMu =
                max(
                    maxMu,
                    fabs(
                        chemicalPotential[y][x]
                    )
                );


            phiSum
                +=
            phi[y][x];


            rhoSum
                +=
            rho[y][x];
        }
    }


    double relativePhiMassChange =
        0.0;


    double relativeRhoMassChange =
        0.0;


    if (
        fabs(
            initialPhiSum
        )
        >
        1.0e-30
    )
    {
        relativePhiMassChange =
            (
                phiSum
                -
                initialPhiSum
            )
            /
            initialPhiSum;
    }


    if (
        fabs(
            initialRhoSum
        )
        >
        1.0e-30
    )
    {
        relativeRhoMassChange =
            (
                rhoSum
                -
                initialRhoSum
            )
            /
            initialRhoSum;
    }


    cout
        << "\n--------------------------------------------\n";

    cout
        << "Step = "
        << step
        << "\n";

    cout
        << "min(phi)              = "
        << minPhi
        << "\n";

    cout
        << "max(phi)              = "
        << maxPhi
        << "\n";

    cout
        << "min(rho)              = "
        << minRho
        << "\n";

    cout
        << "max(rho)              = "
        << maxRho
        << "\n";

    cout
        << "max|mu|               = "
        << maxMu
        << "\n";

    cout
        << "relative phi mass     = "
        << relativePhiMassChange
        << "\n";

    cout
        << "relative rho mass     = "
        << relativeRhoMassChange
        << "\n";
}


// ============================================================
// Save phi to CSV
// ============================================================

void savePhiCSV(
    const char* filename
)
{
    ofstream fout(
        filename
    );


    fout
        << scientific
        << setprecision(12);


    for (int y = 0;
         y < NY;
         y++)
    {
        for (int x = 0;
             x < NX;
             x++)
        {
            fout
                << phi[y][x];


            if (
                x
                <
                NX - 1
            )
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