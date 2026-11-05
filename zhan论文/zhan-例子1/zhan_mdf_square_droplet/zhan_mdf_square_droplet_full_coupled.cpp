#include <iostream>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <cstdio>
using namespace std;


// ============================================================
// 1. Basic parameters
// ============================================================

const double Lx = 2.0;
const double Ly = 2.0;

const double dx = 0.01;
const double W  = 0.04;

const double c  = 100.0;
const double dt = dx / c;

const int NX = static_cast<int>(Lx / dx);
const int NY = static_cast<int>(Ly / dx);

const int Q5 = 5;
const int Q9 = 9;

const double cs2 =
    c * c / 3.0;


// ============================================================
// 2. Physical parameters
// ============================================================

const double rhoA = 1000.0;
const double rhoB = 1.0;

const double phiA =  1.0;
const double phiB = -1.0;

const double nuA = 0.1;
const double nuB = 0.1;

const double mobility =
    0.1;

const double eta =
    1.0;

const double sigma =
    0.001;

const double Rfinal =
    0.5;

const double PI =
    3.14159265358979323846;


// ============================================================
// 3. Free-energy parameters
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
// 4. Initial square size
// ============================================================

const double squareLength =
    sqrt(
        PI
        *
        Rfinal
        *
        Rfinal

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
// 5. Relaxation parameters
// ============================================================

// ------------------------------------------------------------
// CHE:
// M_phi = eta * cs^2 * (1/s1 - 1/2) * dt
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
            cs2
            *
            dt
        )

        +

        0.5
    );


const double Sphase[Q5] =
{
    1.0,
    sPhase1,
    sPhase1,
    1.0,
    1.0
};


// ------------------------------------------------------------
// Flow:
// nu = cs^2 * (1/su1 - 1/2) * dt
//
// Here nuA = nuB = 0.1.
// ------------------------------------------------------------

const double nu =
    0.1;


const double su1 =
    1.0
    /
    (
        nu
        /
        (
            cs2
            *
            dt
        )

        +

        0.5
    );


const double Sflow[Q9] =
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


// ============================================================
// 7. D2Q9
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
// 8. D2Q5 transformation matrix
// ============================================================

const double M5[Q5][Q5] =
{
    {1,1,1,1,1},

    {0,1,0,-1,0},

    {0,0,1,0,-1},

    {0,1,0,1,0},

    {0,0,1,0,1}
};


const double Minv5[Q5][Q5] =
{
    {1,0,0,-1,-1},

    {0,0.5,0,0.5,0},

    {0,0,0.5,0,0.5},

    {0,-0.5,0,0.5,0},

    {0,0,-0.5,0,0.5}
};


// ============================================================
// 9. D2Q9 transformation matrix
// ============================================================

const double M9[Q9][Q9] =
{
    {1,1,1,1,1,1,1,1,1},

    {0,1,0,-1,0,1,-1,-1,1},

    {0,0,1,0,-1,1,1,-1,-1},

    {0,1,0,1,0,1,1,1,1},

    {0,0,0,0,0,1,-1,1,-1},

    {0,0,1,0,1,1,1,1,1},

    {0,0,0,0,0,1,1,-1,-1},

    {0,0,0,0,0,1,-1,-1,1},

    {0,0,0,0,0,1,1,1,1}
};


const double Minv9[Q9][Q9] =
{
    {1,0,0,-1,0,-1,0,0,1},

    {0,0.5,0,0.5,0,0,0,-0.5,-0.5},

    {0,0,0.5,0,0,0.5,-0.5,0,-0.5},

    {0,-0.5,0,0.5,0,0,0,0.5,-0.5},

    {0,0,-0.5,0,0,0.5,0.5,0,-0.5},

    {0,0,0,0,0.25,0,0.25,0.25,0.25},

    {0,0,0,0,-0.25,0,0.25,-0.25,0.25},

    {0,0,0,0,0.25,0,-0.25,-0.25,0.25},

    {0,0,0,0,-0.25,0,-0.25,0.25,0.25}
};


// ============================================================
// 10. Number of steps
// ============================================================

const int nSteps =
    500000;


// ============================================================
// 11. Macroscopic fields
// ============================================================

double phi[NY][NX];

double rho[NY][NX];

double mu[NY][NX];

double gradPhiX[NY][NX];
double gradPhiY[NY][NX];

double gradRhoX[NY][NX];
double gradRhoY[NY][NX];

double forceX[NY][NX];
double forceY[NY][NX];

double ux[NY][NX];
double uy[NY][NX];

double pressure[NY][NX];

double dudx[NY][NX];
double dudy[NY][NX];

double dvdx[NY][NX];
double dvdy[NY][NX];


// ============================================================
// 12. Old-time fields
// ============================================================

// ------------------------------------------------------------
// Phase G:
//
// d_t(phi u)
// ------------------------------------------------------------

double oldPhiUx[NY][NX];
double oldPhiUy[NY][NX];


// ------------------------------------------------------------
// Flow G:
//
// d_t(Q)
// ------------------------------------------------------------

double oldQxx[NY][NX];
double oldQxy[NY][NX];
double oldQyy[NY][NX];


// ============================================================
// 13. Phase distributions
// ============================================================

double gPhase[NY][NX][Q5];

double gPhasePost[NY][NX][Q5];

double gPhaseNext[NY][NX][Q5];


// ============================================================
// 14. Flow distributions
// ============================================================

double g1[NY][NX][Q9];
double g2[NY][NX][Q9];

double g1Post[NY][NX][Q9];
double g2Post[NY][NX][Q9];

double g1Next[NY][NX][Q9];
double g2Next[NY][NX][Q9];


// ============================================================
// 15. Initial diagnostics
// ============================================================

double initialPhaseMass =
    0.0;

double initialDensityMass =
    0.0;

double initialMomentumX =
    0.0;

double initialMomentumY =
    0.0;


// ============================================================
// Function declarations
// ============================================================

void checkMatrices();

void initializePhi();

void recoverPhi();

void computeDensity();

void computePhaseDerivatives();

void computeSurfaceForce();

void computeGradRho();

void initializeVelocityPressure();

void buildPhaseEquilibriumMoments(
    int y,
    int x,
    double meq[Q5]
);

void buildPhaseAuxiliaryMoments(
    int y,
    int x,
    double mG[Q5]
);

void buildFlowEquilibriumMoments(
    int y,
    int x,
    double meq1[Q9],
    double meq2[Q9]
);

void buildFlowForceMoments(
    int y,
    int x,
    double mF1[Q9],
    double mF2[Q9]
);

void buildFlowAuxiliaryMoments(
    int y,
    int x,
    double mG1[Q9],
    double mG2[Q9]
);

void initializePhaseDistribution();

void initializeFlowDistribution();

void recoverVelocity();

void computeVelocityGradients();

void recoverPressure();

void initializeOldFields();

void saveCurrentOldFields();

void collidePhase();

void collideFlow();

void streamPhase();

void streamFlow();

bool checkNaN(
    int step
);

void printDiagnostics(
    int step
);
void savePhi(
    int step
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
        << "====================================================\n";

    cout
        << "Zhan 2023 square droplet\n";

    cout
        << "FULL COUPLED CHE + MDF FLOW TEST\n";

    cout
        << "====================================================\n\n";


    cout
        << "NX             = "
        << NX
        << "\n";

    cout
        << "NY             = "
        << NY
        << "\n";

    cout
        << "dx             = "
        << dx
        << "\n";

    cout
        << "dt             = "
        << dt
        << "\n";

    cout
        << "c              = "
        << c
        << "\n";

    cout
        << "cs2            = "
        << cs2
        << "\n";

    cout
        << "W              = "
        << W
        << "\n";

    cout
        << "mobility       = "
        << mobility
        << "\n";

    cout
        << "nu             = "
        << nu
        << "\n";

    cout
        << "sigma          = "
        << sigma
        << "\n";

    cout
        << "betaPF         = "
        << betaPF
        << "\n";

    cout
        << "kappaPF        = "
        << kappaPF
        << "\n";

    cout
        << "sPhase1        = "
        << sPhase1
        << "\n";

    cout
        << "su1            = "
        << su1
        << "\n\n";


    // ========================================================
    // Check matrices
    // ========================================================

    checkMatrices();


    // ========================================================
    // Initial phase
    // ========================================================

    initializePhi();

    computeDensity();

    computePhaseDerivatives();

    computeSurfaceForce();

    computeGradRho();


    // ========================================================
    // Initial velocity and pressure
    // ========================================================

    initializeVelocityPressure();


    // ========================================================
    // Initialize phase distribution
    // ========================================================

    initializePhaseDistribution();


    // ========================================================
    // Initialize flow distribution
    // ========================================================

    initializeFlowDistribution();


    // ========================================================
    // Recover initial macros
    // ========================================================

    recoverPhi();

    computeDensity();

    computePhaseDerivatives();

    computeSurfaceForce();

    computeGradRho();

    recoverVelocity();

    computeVelocityGradients();

    recoverPressure();


    // ========================================================
    // Initialize old-time quantities
    // ========================================================

    initializeOldFields();


    // ========================================================
    // Initial conserved quantities
    // ========================================================

    initialPhaseMass =
        0.0;

    initialDensityMass =
        0.0;

    initialMomentumX =
        0.0;

    initialMomentumY =
        0.0;


    for (int y = 0;
         y < NY;
         y++)
    {
        for (int x = 0;
             x < NX;
             x++)
        {
            initialPhaseMass
                +=
            phi[y][x];


            initialDensityMass
                +=
            rho[y][x];


            initialMomentumX
                +=
            rho[y][x]
            *
            ux[y][x];


            initialMomentumY
                +=
            rho[y][x]
            *
            uy[y][x];
        }
    }


    printDiagnostics(
        0
    );

    savePhi(
        0
    );
    // ========================================================
    // Time loop
    // ========================================================

    for (int step = 1;
         step <= nSteps;
         step++)
    {
        // ====================================================
        // 1. Recover phi^n
        // ====================================================

        recoverPhi();


        // ====================================================
        // 2. rho(phi^n)
        // ====================================================

        computeDensity();


        // ====================================================
        // 3. grad(phi), lap(phi), mu
        // ====================================================

        computePhaseDerivatives();


        // ====================================================
        // 4. Surface force
        // ====================================================

        computeSurfaceForce();


        // ====================================================
        // 5. grad(rho)
        // ====================================================

        computeGradRho();


        // ====================================================
        // 6. Recover velocity
        // ====================================================

        recoverVelocity();


        // ====================================================
        // 7. velocity gradients
        // ====================================================

        computeVelocityGradients();


        // ====================================================
        // 8. pressure
        // ====================================================

        recoverPressure();


        // ====================================================
        // 9. Phase collision
        // ====================================================

        collidePhase();


        // ====================================================
        // 10. Flow collision
        // ====================================================

        collideFlow();


        // ====================================================
        // 11. Save n-level old quantities
        // ====================================================

        saveCurrentOldFields();


        // ====================================================
        // 12. Streaming
        // ====================================================

        streamPhase();

        streamFlow();


        // ====================================================
        // 13. Recover n+1 for diagnostics
        // ====================================================

        recoverPhi();

        computeDensity();

        computePhaseDerivatives();

        computeSurfaceForce();

        computeGradRho();

        recoverVelocity();

        computeVelocityGradients();

        recoverPressure();


        // ====================================================
        // 14. Safety check
        // ====================================================

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


        if (step % 100000 == 0)
        {
          printDiagnostics(step);
          savePhi(step);
        }
        
    }


    cout
        << "\n====================================================\n";

    cout
        << "FULL COUPLED 10-step test finished.\n";

    cout
        << "====================================================\n";


    return 0;
}


// ============================================================
// Check M5 and M9
// ============================================================

void checkMatrices()
{
    double maxError5 =
        0.0;

    double maxError9 =
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
                (
                    i == j
                )
                ?
                1.0
                :
                0.0;


            maxError5 =
                max(
                    maxError5,
                    fabs(
                        value
                        -
                        exact
                    )
                );
        }
    }


    for (int i = 0;
         i < Q9;
         i++)
    {
        for (int j = 0;
             j < Q9;
             j++)
        {
            double value =
                0.0;


            for (int k = 0;
                 k < Q9;
                 k++)
            {
                value
                    +=
                M9[i][k]
                *
                Minv9[k][j];
            }


            const double exact =
                (
                    i == j
                )
                ?
                1.0
                :
                0.0;


            maxError9 =
                max(
                    maxError9,
                    fabs(
                        value
                        -
                        exact
                    )
                );
        }
    }


    cout
        << "Max |M5*Minv5-I| = "
        << maxError5
        << "\n";

    cout
        << "Max |M9*Minv9-I| = "
        << maxError9
        << "\n\n";
}


// ============================================================
// Initial square droplet
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


    for (int j = 0;
         j < NY;
         j++)
    {
        for (int i = 0;
             i < NX;
             i++)
        {
            const double x =
                -1.0
                +
                (
                    i
                    +
                    0.5
                )
                *
                dx;


            const double y =
                -1.0
                +
                (
                    j
                    +
                    0.5
                )
                *
                dx;


            double d =
                0.0;


            if (
                x <= a
                &&
                y <= a
            )
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

            else if (
                x >= b
                &&
                y <= a
            )
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

            else if (
                x <= a
                &&
                y >= b
            )
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

            else if (
                x >= b
                &&
                y >= b
            )
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

            else if (
                x > a
                &&
                x < b
                &&
                y <= a
            )
            {
                d =
                    0.5
                    *
                    squareLength
                    +
                    y;
            }

            else if (
                x > a
                &&
                x < b
                &&
                y >= b
            )
            {
                d =
                    0.5
                    *
                    squareLength
                    -
                    y;
            }

            else if (
                x <= a
                &&
                y > a
                &&
                y < b
            )
            {
                d =
                    0.5
                    *
                    squareLength
                    +
                    x;
            }

            else if (
                x >= b
                &&
                y > a
                &&
                y < b
            )
            {
                d =
                    0.5
                    *
                    squareLength
                    -
                    x;
            }

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
                d
                /
                W
            );
        }
    }
}


// ============================================================
// Recover phi = sum_i g_i0
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
            double value =
                0.0;


            for (int q = 0;
                 q < Q5;
                 q++)
            {
                value
                    +=
                gPhase[y][x][q];
            }


            phi[y][x] =
                value;
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
// grad(phi), lap(phi), chemical potential
// ============================================================

void computePhaseDerivatives()
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

            double lap =
                0.0;


            for (int q = 1;
                 q < Q5;
                 q++)
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
                    c
                    *
                    cx5[q];


                const double ciy =
                    c
                    *
                    cy5[q];


                gx
                    +=
                w5[q]
                *
                cix
                *
                phi[yn][xn];


                gy
                    +=
                w5[q]
                *
                ciy
                *
                phi[yn][xn];


                lap
                    +=
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


            lap
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


            mu[y][x]
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
            lap;
        }
    }
}


// ============================================================
// Surface force F_s = mu grad(phi)
// ============================================================

void computeSurfaceForce()
{
    for (int y = 0;
         y < NY;
         y++)
    {
        for (int x = 0;
             x < NX;
             x++)
        {
            forceX[y][x]
                =
            mu[y][x]
            *
            gradPhiX[y][x];


            forceY[y][x]
                =
            mu[y][x]
            *
            gradPhiY[y][x];
        }
    }
}


// ============================================================
// grad(rho)
// ============================================================

void computeGradRho()
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
                 q < Q9;
                 q++)
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
                    c
                    *
                    cx9[q];


                const double ciy =
                    c
                    *
                    cy9[q];


                gx
                    +=
                w9[q]
                *
                cix
                *
                rho[yn][xn];


                gy
                    +=
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
// Initial velocity / pressure
// ============================================================

void initializeVelocityPressure()
{
    for (int y = 0;
         y < NY;
         y++)
    {
        for (int x = 0;
             x < NX;
             x++)
        {
            ux[y][x] =
                0.0;

            uy[y][x] =
                0.0;

            pressure[y][x] =
                0.0;

            dudx[y][x] =
                0.0;

            dudy[y][x] =
                0.0;

            dvdx[y][x] =
                0.0;

            dvdy[y][x] =
                0.0;
        }
    }
}


// ============================================================
// Phase equilibrium
//
// Eq.(9a) -> velocity space -> M5
// ============================================================

void buildPhaseEquilibriumMoments(
    int y,
    int x,
    double meq[Q5]
)
{
    double feq[Q5];


    const double ph =
        phi[y][x];

    const double chem =
        mu[y][x];

    const double u =
        ux[y][x];

    const double v =
        uy[y][x];


    for (int q = 0;
         q < Q5;
         q++)
    {
        if (
            q == 0
        )
        {
            feq[q]
                =
            ph

            +

            (
                w5[q]
                -
                1.0
            )
            *
            eta
            *
            chem;
        }
        else
        {
            const double cix =
                c
                *
                cx5[q];


            const double ciy =
                c
                *
                cy5[q];


            feq[q]
                =
            w5[q]
            *
            eta
            *
            chem

            +

            w5[q]
            *
            (
                cix
                *
                ph
                *
                u

                +

                ciy
                *
                ph
                *
                v
            )
            /
            cs2;
        }
    }


    for (int k = 0;
         k < Q5;
         k++)
    {
        meq[k] =
            0.0;


        for (int q = 0;
             q < Q5;
             q++)
        {
            meq[k]
                +=
            M5[k][q]
            *
            feq[q];
        }
    }
}


// ============================================================
// Phase G
//
// Eq.(9a):
//
// Gi = wi ci_beta d_t(phi u_beta) / cs^2
//
// Then mG = M5 * G
// ============================================================

void buildPhaseAuxiliaryMoments(
    int y,
    int x,
    double mG[Q5]
)
{
    const double currentPhiUx =
        phi[y][x]
        *
        ux[y][x];


    const double currentPhiUy =
        phi[y][x]
        *
        uy[y][x];


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


    double Gi[Q5];


    for (int q = 0;
         q < Q5;
         q++)
    {
        const double cix =
            c
            *
            cx5[q];


        const double ciy =
            c
            *
            cy5[q];


        Gi[q]
            =
        w5[q]
        *
        (
            cix
            *
            dPhiUxDt

            +

            ciy
            *
            dPhiUyDt
        )
        /
        cs2;
    }


    for (int k = 0;
         k < Q5;
         k++)
    {
        mG[k] =
            0.0;


        for (int q = 0;
             q < Q5;
             q++)
        {
            mG[k]
                +=
            M5[k][q]
            *
            Gi[q];
        }
    }
}


// ============================================================
// Flow equilibrium
//
// Eq.(10), C = 0
//
// velocity space -> M9
// ============================================================

void buildFlowEquilibriumMoments(
    int y,
    int x,
    double meq1[Q9],
    double meq2[Q9]
)
{
    const double r =
        rho[y][x];

    const double u =
        ux[y][x];

    const double v =
        uy[y][x];

    const double p =
        pressure[y][x];


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


    double feq1[Q9];
    double feq2[Q9];


    for (int q = 0;
         q < Q9;
         q++)
    {
        const double cix =
            c
            *
            cx9[q];


        const double ciy =
            c
            *
            cy9[q];


        feq1[q]
            =
        w9[q]
        *
        (
            r * u

            +

            (
                cix * Qxx
                +
                ciy * Qxy
            )
            /
            cs2
        );


        feq2[q]
            =
        w9[q]
        *
        (
            r * v

            +

            (
                cix * Qxy
                +
                ciy * Qyy
            )
            /
            cs2
        );
    }


    for (int k = 0;
         k < Q9;
         k++)
    {
        meq1[k] =
            0.0;

        meq2[k] =
            0.0;


        for (int q = 0;
             q < Q9;
             q++)
        {
            meq1[k]
                +=
            M9[k][q]
            *
            feq1[q];


            meq2[k]
                +=
            M9[k][q]
            *
            feq2[q];
        }
    }
}


// ============================================================
// Flow force
//
// Eq.(11):
//
// Fi,alpha = wi F_alpha
//
// velocity space -> M9
// ============================================================

void buildFlowForceMoments(
    int y,
    int x,
    double mF1[Q9],
    double mF2[Q9]
)
{
    double Fi1[Q9];
    double Fi2[Q9];


    for (int q = 0;
         q < Q9;
         q++)
    {
        Fi1[q]
            =
        w9[q]
        *
        forceX[y][x];


        Fi2[q]
            =
        w9[q]
        *
        forceY[y][x];
    }


    for (int k = 0;
         k < Q9;
         k++)
    {
        mF1[k] =
            0.0;

        mF2[k] =
            0.0;


        for (int q = 0;
             q < Q9;
             q++)
        {
            mF1[k]
                +=
            M9[k][q]
            *
            Fi1[q];


            mF2[k]
                +=
            M9[k][q]
            *
            Fi2[q];
        }
    }
}


// ============================================================
// Flow auxiliary G
//
// Eq.(37), C = 0
//
// velocity space -> M9
// ============================================================

void buildFlowAuxiliaryMoments(
    int y,
    int x,
    double mG1[Q9],
    double mG2[Q9]
)
{
    const double r =
        rho[y][x];

    const double u =
        ux[y][x];

    const double v =
        uy[y][x];

    const double p =
        pressure[y][x];


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


    const double dQxxDt =
        (
            Qxx
            -
            oldQxx[y][x]
        )
        /
        dt;


    const double dQxyDt =
        (
            Qxy
            -
            oldQxy[y][x]
        )
        /
        dt;


    const double dQyyDt =
        (
            Qyy
            -
            oldQyy[y][x]
        )
        /
        dt;


    // --------------------------------------------------------
    // alpha = x
    // --------------------------------------------------------

    const double A1x =
        u
        *
        gradRhoX[y][x]

        -

        r
        *
        dudx[y][x];


    const double A1y =
        u
        *
        gradRhoY[y][x]

        -

        r
        *
        dvdx[y][x];


    // --------------------------------------------------------
    // alpha = y
    // --------------------------------------------------------

    const double A2x =
        v
        *
        gradRhoX[y][x]

        -

        r
        *
        dudy[y][x];


    const double A2y =
        v
        *
        gradRhoY[y][x]

        -

        r
        *
        dvdy[y][x];


    double Gi1[Q9];
    double Gi2[Q9];


    for (int q = 0;
         q < Q9;
         q++)
    {
        const double cix =
            c
            *
            cx9[q];


        const double ciy =
            c
            *
            cy9[q];


        Gi1[q]
            =
        w9[q]
        *
        (
            (
                cix
                *
                dQxxDt

                +

                ciy
                *
                dQxyDt
            )
            /
            cs2

            +

            cix
            *
            A1x

            +

            ciy
            *
            A1y
        );


        Gi2[q]
            =
        w9[q]
        *
        (
            (
                cix
                *
                dQxyDt

                +

                ciy
                *
                dQyyDt
            )
            /
            cs2

            +

            cix
            *
            A2x

            +

            ciy
            *
            A2y
        );
    }


    for (int k = 0;
         k < Q9;
         k++)
    {
        mG1[k] =
            0.0;

        mG2[k] =
            0.0;


        for (int q = 0;
             q < Q9;
             q++)
        {
            mG1[k]
                +=
            M9[k][q]
            *
            Gi1[q];


            mG2[k]
                +=
            M9[k][q]
            *
            Gi2[q];
        }
    }
}


// ============================================================
// Initialize phase distribution
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


            buildPhaseEquilibriumMoments(
                y,
                x,
                meq
            );


            // ------------------------------------------------
            // S = 0
            //
            // therefore transformed phase g = f
            // ------------------------------------------------

            for (int q = 0;
                 q < Q5;
                 q++)
            {
                gPhase[y][x][q] =
                    0.0;


                for (int k = 0;
                     k < Q5;
                     k++)
                {
                    gPhase[y][x][q]
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
// Initialize flow distribution
//
// g = f - dt F / 2
// ============================================================

void initializeFlowDistribution()
{
    for (int y = 0;
         y < NY;
         y++)
    {
        for (int x = 0;
             x < NX;
             x++)
        {
            double meq1[Q9];
            double meq2[Q9];

            double mF1[Q9];
            double mF2[Q9];

            double mg1[Q9];
            double mg2[Q9];


            buildFlowEquilibriumMoments(
                y,
                x,
                meq1,
                meq2
            );


            buildFlowForceMoments(
                y,
                x,
                mF1,
                mF2
            );


            for (int k = 0;
                 k < Q9;
                 k++)
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


            for (int q = 0;
                 q < Q9;
                 q++)
            {
                g1[y][x][q] =
                    0.0;

                g2[y][x][q] =
                    0.0;


                for (int k = 0;
                     k < Q9;
                     k++)
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
// Recover velocity
// ============================================================

void recoverVelocity()
{
    for (int y = 0;
         y < NY;
         y++)
    {
        for (int x = 0;
             x < NX;
             x++)
        {
            double sum1 =
                0.0;

            double sum2 =
                0.0;


            for (int q = 0;
                 q < Q9;
                 q++)
            {
                sum1
                    +=
                g1[y][x][q];


                sum2
                    +=
                g2[y][x][q];
            }


            ux[y][x]
                =
            (
                sum1

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
                sum2

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
// Velocity gradients
// ============================================================

void computeVelocityGradients()
{
    for (int y = 0;
         y < NY;
         y++)
    {
        for (int x = 0;
             x < NX;
             x++)
        {
            double ux_x =
                0.0;

            double ux_y =
                0.0;

            double uy_x =
                0.0;

            double uy_y =
                0.0;


            for (int q = 1;
                 q < Q9;
                 q++)
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
                    c
                    *
                    cx9[q];


                const double ciy =
                    c
                    *
                    cy9[q];


                ux_x
                    +=
                w9[q]
                *
                cix
                *
                ux[yn][xn];


                ux_y
                    +=
                w9[q]
                *
                ciy
                *
                ux[yn][xn];


                uy_x
                    +=
                w9[q]
                *
                cix
                *
                uy[yn][xn];


                uy_y
                    +=
                w9[q]
                *
                ciy
                *
                uy[yn][xn];
            }


            dudx[y][x]
                =
            ux_x
            /
            (
                cs2
                *
                dt
            );


            dudy[y][x]
                =
            ux_y
            /
            (
                cs2
                *
                dt
            );


            dvdx[y][x]
                =
            uy_x
            /
            (
                cs2
                *
                dt
            );


            dvdy[y][x]
                =
            uy_y
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
// Pressure recovery
//
// simplified Eq.(32)
//
// C = 0
//
// last temporal term omitted
// ============================================================

void recoverPressure()
{
    for (int y = 0;
         y < NY;
         y++)
    {
        for (int x = 0;
             x < NX;
             x++)
        {
            double firstMomentX =
                0.0;

            double firstMomentY =
                0.0;


            for (int q = 0;
                 q < Q9;
                 q++)
            {
                const double cix =
                    c
                    *
                    cx9[q];


                const double ciy =
                    c
                    *
                    cy9[q];


                firstMomentX
                    +=
                cix
                *
                g1[y][x][q];


                firstMomentY
                    +=
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
                    ux[y][x]
                    *
                    ux[y][x]

                    +

                    uy[y][x]
                    *
                    uy[y][x]
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
// Initialize old fields
//
// First temporal derivatives = 0
// ============================================================

void initializeOldFields()
{
    for (int y = 0;
         y < NY;
         y++)
    {
        for (int x = 0;
             x < NX;
             x++)
        {
            oldPhiUx[y][x]
                =
            phi[y][x]
            *
            ux[y][x];


            oldPhiUy[y][x]
                =
            phi[y][x]
            *
            uy[y][x];


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
// Save current n-level fields
// ============================================================

void saveCurrentOldFields()
{
    for (int y = 0;
         y < NY;
         y++)
    {
        for (int x = 0;
             x < NX;
             x++)
        {
            oldPhiUx[y][x]
                =
            phi[y][x]
            *
            ux[y][x];


            oldPhiUy[y][x]
                =
            phi[y][x]
            *
            uy[y][x];


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
// Phase collision
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
                    gPhase[y][x][q];
                }
            }


            buildPhaseEquilibriumMoments(
                y,
                x,
                meq
            );


            buildPhaseAuxiliaryMoments(
                y,
                x,
                mG
            );


            // ------------------------------------------------
            // S = 0
            //
            // mPost =
            // (I-S)m
            // + S meq
            // + dt(I-S/2)mG
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
                    Sphase[k]
                )
                *
                m[k]

                +

                Sphase[k]
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
                    Sphase[k]
                )
                *
                mG[k];
            }


            for (int q = 0;
                 q < Q5;
                 q++)
            {
                gPhasePost[y][x][q] =
                    0.0;


                for (int k = 0;
                     k < Q5;
                     k++)
                {
                    gPhasePost[y][x][q]
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
// Flow collision
// ============================================================

void collideFlow()
{
    for (int y = 0;
         y < NY;
         y++)
    {
        for (int x = 0;
             x < NX;
             x++)
        {
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


            for (int k = 0;
                 k < Q9;
                 k++)
            {
                m1[k] =
                    0.0;

                m2[k] =
                    0.0;


                for (int q = 0;
                     q < Q9;
                     q++)
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


            buildFlowEquilibriumMoments(
                y,
                x,
                meq1,
                meq2
            );


            buildFlowForceMoments(
                y,
                x,
                mF1,
                mF2
            );


            buildFlowAuxiliaryMoments(
                y,
                x,
                mG1,
                mG2
            );


            for (int k = 0;
                 k < Q9;
                 k++)
            {
                mPost1[k]
                    =
                (
                    1.0
                    -
                    Sflow[k]
                )
                *
                m1[k]

                +

                Sflow[k]
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
                    Sflow[k]
                )
                *
                (
                    mF1[k]
                    +
                    mG1[k]
                );


                mPost2[k]
                    =
                (
                    1.0
                    -
                    Sflow[k]
                )
                *
                m2[k]

                +

                Sflow[k]
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
                    Sflow[k]
                )
                *
                (
                    mF2[k]
                    +
                    mG2[k]
                );
            }


            for (int q = 0;
                 q < Q9;
                 q++)
            {
                g1Post[y][x][q] =
                    0.0;

                g2Post[y][x][q] =
                    0.0;


                for (int k = 0;
                     k < Q9;
                     k++)
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
// Phase streaming
// ============================================================

void streamPhase()
{
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
                gPhaseNext[y][x][q] =
                    0.0;
            }
        }
    }


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


                gPhaseNext[yn][xn][q]
                    =
                gPhasePost[y][x][q];
            }
        }
    }


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
                gPhase[y][x][q]
                    =
                gPhaseNext[y][x][q];
            }
        }
    }
}


// ============================================================
// Flow streaming
// ============================================================

void streamFlow()
{
    for (int y = 0;
         y < NY;
         y++)
    {
        for (int x = 0;
             x < NX;
             x++)
        {
            for (int q = 0;
                 q < Q9;
                 q++)
            {
                g1Next[y][x][q] =
                    0.0;

                g2Next[y][x][q] =
                    0.0;
            }
        }
    }


    for (int y = 0;
         y < NY;
         y++)
    {
        for (int x = 0;
             x < NX;
             x++)
        {
            for (int q = 0;
                 q < Q9;
                 q++)
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


    for (int y = 0;
         y < NY;
         y++)
    {
        for (int x = 0;
             x < NX;
             x++)
        {
            for (int q = 0;
                 q < Q9;
                 q++)
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
                    ux[y][x]
                )
                ||
                !isfinite(
                    uy[y][x]
                )
                ||
                !isfinite(
                    pressure[y][x]
                )
                ||
                !isfinite(
                    mu[y][x]
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

    double maxForce =
        0.0;

    double maxUx =
        0.0;

    double maxUy =
        0.0;

    double maxSpeed =
        0.0;

    double maxP =
        0.0;

    double phaseMass =
        0.0;

    double densityMass =
        0.0;

    double momentumX =
        0.0;

    double momentumY =
        0.0;

    double divSquareSum =
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
                        mu[y][x]
                    )
                );


            const double Fmag =
                sqrt(
                    forceX[y][x]
                    *
                    forceX[y][x]

                    +

                    forceY[y][x]
                    *
                    forceY[y][x]
                );


            maxForce =
                max(
                    maxForce,
                    Fmag
                );


            maxUx =
                max(
                    maxUx,
                    fabs(
                        ux[y][x]
                    )
                );


            maxUy =
                max(
                    maxUy,
                    fabs(
                        uy[y][x]
                    )
                );


            const double speed =
                sqrt(
                    ux[y][x]
                    *
                    ux[y][x]

                    +

                    uy[y][x]
                    *
                    uy[y][x]
                );


            maxSpeed =
                max(
                    maxSpeed,
                    speed
                );


            maxP =
                max(
                    maxP,
                    fabs(
                        pressure[y][x]
                    )
                );


            phaseMass
                +=
            phi[y][x];


            densityMass
                +=
            rho[y][x];


            momentumX
                +=
            rho[y][x]
            *
            ux[y][x];


            momentumY
                +=
            rho[y][x]
            *
            uy[y][x];


            const double divU =
                dudx[y][x]
                +
                dvdy[y][x];


            divSquareSum
                +=
            divU
            *
            divU;
        }
    }


    const double relativePhaseMass =
        (
            phaseMass
            -
            initialPhaseMass
        )
        /
        (
            fabs(
                initialPhaseMass
            )
            +
            1.0e-30
        );


    const double relativeDensityMass =
        (
            densityMass
            -
            initialDensityMass
        )
        /
        (
            fabs(
                initialDensityMass
            )
            +
            1.0e-30
        );


    const double deltaMomentumX =
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


    const double deltaMomentumY =
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
            divSquareSum
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
        << "phi range              = ["
        << minPhi
        << ", "
        << maxPhi
        << "]\n";


    cout
        << "rho range              = ["
        << minRho
        << ", "
        << maxRho
        << "]\n";


    cout
        << "max|mu|                = "
        << maxMu
        << "\n";


    cout
        << "max|F|                 = "
        << maxForce
        << "\n";


    cout
        << "max|u|                 = "
        << maxUx
        << "\n";


    cout
        << "max|v|                 = "
        << maxUy
        << "\n";


    cout
        << "max speed              = "
        << maxSpeed
        << "\n";


    cout
        << "max|p|                 = "
        << maxP
        << "\n";


    cout
        << "relative phase mass    = "
        << relativePhaseMass
        << "\n";


    cout
        << "relative density mass  = "
        << relativeDensityMass
        << "\n";


    cout
        << "<Delta rho*u>          = "
        << deltaMomentumX
        << "\n";


    cout
        << "<Delta rho*v>          = "
        << deltaMomentumY
        << "\n";


    cout
        << "D_rms                  = "
        << Drms
        << "\n";
}
// ============================================================
// Save phi field
// ============================================================

void savePhi(
    int step
)
{
    char fileName[100];

    sprintf(
        fileName,
        "phi_%05d.dat",
        step
    );


    FILE *fp =
        fopen(
            fileName,
            "w"
        );


    if (
        fp == nullptr
    )
    {
        cout
            << "Cannot open "
            << fileName
            << "\n";

        return;
    }


    for (int y = 0;
         y < NY;
         y++)
    {
        for (int x = 0;
             x < NX;
             x++)
        {
            fprintf(
                fp,
                "%.15e",
                phi[y][x]
            );


            if (
                x < NX - 1
            )
            {
                fprintf(
                    fp,
                    " "
                );
            }
        }


        fprintf(
            fp,
            "\n"
        );
    }


    fclose(
        fp
    );
}