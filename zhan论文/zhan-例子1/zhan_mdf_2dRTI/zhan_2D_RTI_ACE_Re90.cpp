#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <cstdio>
#include <cstdlib>

using namespace std;

// ============================================================
// Zhan 2023 - 2D Rayleigh-Taylor instability
// ACE + MDF-LBM
//
// Present case:
//     Re = 90
//     At = 0.5
//
// Phase:
//     D2Q5 conservative Allen-Cahn
//
// Flow:
//     D2Q9 MDF-LBM
//
// Velocity gradient in Eq.(37):
//     isotropic Eq.(41)
//
// ACE interface normal:
//     local Eq.(43)
//     IMPORTANT:
//     |grad(phi)| = (-|C| - B) / A
//
// Boundary:
//     x : periodic
//     y : no-flux phase + no-slip flow
// ============================================================


// ============================================================
// 1. Basic parameters
// ============================================================

const int NX = 150;
const int NY = 600;

const int Q5 = 5;
const int Q9 = 9;

const double dx = 1.0;
const double dt = 1.0;

const double c = dx / dt;
const double cs2 = c * c / 3.0;

const double PI =
    3.1415926535897932384626433832795;


// ============================================================
// 2. RTI geometry
//
// domain:
//
// [-d/2,d/2] x [-2d,2d]
// ============================================================

const double d = 150.0;

const double L = d;

const double H = 4.0 * d;


// ============================================================
// 3. Phase field
// ============================================================

const double phiA = 1.0;
const double phiB = 0.0;

const double W = 4.0;


// ============================================================
// 4. Dimensionless parameters
// ============================================================

const double Atwood = 0.5;

const double Reynolds = 3000.0;

const double Peclet = 1000.0;


// ============================================================
// Reference velocity
//
// Ug = sqrt(g*d)
//
// Zu 2020 RTI reference:
// lattice characteristic velocity ~ 0.0707
// ============================================================

// ============================================================
// Gravity and RTI characteristic velocity
//
// Zhan Sec.4.5 does not explicitly give the lattice value of g.
//
// Use a smaller lattice gravity for the MDF-LBM.
// Re and Pe are still preserved because nu and M are
// recomputed from sqrt(g*d).
// ============================================================

const double gravity = 1.0e-5;

const double Ug =
    sqrt(
        gravity
        *
        d
    );


// ============================================================
// RTI dimensionless time
//
// T = t / sqrt[d/(g*At)]
// ============================================================

const double rtiTimeScale =
    sqrt(
        d
        /
        (
            gravity
            *
            Atwood
        )
    );


// ============================================================
// 5. Density
//
// At = (rhoA-rhoB)/(rhoA+rhoB) = 0.5
//
// rhoA = 3
// rhoB = 1
// ============================================================

const double rhoA = 3.0;

const double rhoB = 1.0;

const double rhoReference =
    0.5
    *
    (
        rhoA
        +
        rhoB
    );


// ============================================================
// 6. Viscosity
//
// Re = d*sqrt(gd)/nuA
//
// equal dynamic viscosity:
//
// muA = muB
// ============================================================

const double nuA =
    d
    *
    Ug
    /
    Reynolds;


const double muA =
    rhoA
    *
    nuA;


const double muB =
    muA;


const double nuB =
    muB
    /
    rhoB;


// ============================================================
// 7. Weak surface tension
//
// Zhan Sec.4.5 does not explicitly state sigma.
//
// Following Zu 2020:
//
// Ca = muA*Ug/sigma = 100
// ============================================================

const double Capillary = 100.0;


const double sigma =
    muA
    *
    Ug
    /
    Capillary;


// ============================================================
// 8. Free-energy parameters
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
// 9. ACE mobility
//
// Pe = d*sqrt(gd)/M
// ============================================================

const double mobility =
    d
    *
    Ug
    /
    Peclet;


const double eta = 1.0;

const double normalEpsilon = 1.0e-14;


// ============================================================
// 10. Phase relaxation
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


// ============================================================
// Phase MRT relaxation
//
// Zhan Sec.4.5:
//
//     s_2^0 = 1.7
//
// D2Q5:
// S^0 = diag(
//     s0^0,
//     s1^0,
//     s1^0,
//     s2^0,
//     s2^0
// )
//
// s1^0 is determined by mobility.
// ============================================================

const double sPhase2 = 1.7;


const double Sphase[Q5] =
{
    1.0,

    sPhase1,
    sPhase1,

    sPhase2,
    sPhase2
};


// ============================================================
// Flow MRT relaxation
//
// Zhan states that unspecified relaxation parameters
// are set to 1.
//
// s1^u is still calculated locally from viscosity.
//
// IMPORTANT:
// s2^0 = 1.7 belongs to PHASE FIELD,
// not to FLOW FIELD.
// ============================================================

const double su0 = 1.0;

const double su2 = 1.0;

const double su3 = 1.0;

const double su4 = 1.0;

// ============================================================
// 12. Paper output times
// ============================================================

const int stepT1 =
    static_cast<int>(
        lround(
            1.0
            *
            rtiTimeScale
            /
            dt
        )
    );


const int stepT175 =
    static_cast<int>(
        lround(
            1.75
            *
            rtiTimeScale
            /
            dt
        )
    );


const int stepT25 =
    static_cast<int>(
        lround(
            2.5
            *
            rtiTimeScale
            /
            dt
        )
    );


const int stepT3 =
    static_cast<int>(
        lround(
            3.0
            *
            rtiTimeScale
            /
            dt
        )
    );


// ============================================================
// 13. D2Q5
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


const int opp5[Q5] =
{
    0,
    3,
    4,
    1,
    2
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
// 14. D2Q9
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


const int opp9[Q9] =
{
    0,
    3,
    4,
    1,
    2,
    7,
    8,
    5,
    6
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
// 15. D2Q5 MRT matrix
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
// 16. D2Q9 MRT matrix
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
// 17. Macroscopic fields
// ============================================================

static double phi[NY][NX];

static double rho[NY][NX];

static double nu[NY][NX];

static double muDyn[NY][NX];

static double chemical[NY][NX];


static double gradPhiX[NY][NX];

static double gradPhiY[NY][NX];


static double gradRhoX[NY][NX];

static double gradRhoY[NY][NX];


static double forceX[NY][NX];

static double forceY[NY][NX];


static double ux[NY][NX];

static double uy[NY][NX];


static double pressure[NY][NX];


static double dudx[NY][NX];

static double dudy[NY][NX];

static double dvdx[NY][NX];

static double dvdy[NY][NX];


// ============================================================
// 18. Previous-time fields
// ============================================================

static double oldPhiUx[NY][NX];

static double oldPhiUy[NY][NX];


static double oldQxx[NY][NX];

static double oldQxy[NY][NX];

static double oldQyx[NY][NX];

static double oldQyy[NY][NX];


// ============================================================
// 19. Phase distributions
// ============================================================

static double gPhase[NY][NX][Q5];

static double gPhasePost[NY][NX][Q5];

static double gPhaseNext[NY][NX][Q5];


// ============================================================
// 20. Flow distributions
// ============================================================

static double g1[NY][NX][Q9];

static double g2[NY][NX][Q9];


static double g1Post[NY][NX][Q9];

static double g2Post[NY][NX][Q9];


static double g1Next[NY][NX][Q9];

static double g2Next[NY][NX][Q9];


// ============================================================
// Diagnostic
// ============================================================

double initialPhaseMass = 0.0;


// ============================================================
// Function declarations
// ============================================================

void checkMatrices();

void initializeMacroscopicFields();

void computeDensityViscosity();

void computePhaseDerivatives();

void computeGradRho();

void computeSurfaceForce();

void initializePhaseDistribution();

void initializeFlowDistribution();

void initializeOldFields();

void saveCurrentOldFields();

void recoverPhi();

void recoverVelocity();

void recoverPressure();

void computeVelocityGradientsIsotropic();

void collidePhase();

void collideFlow();

void streamPhase();

void streamFlow();

bool checkNaN(int step);

void printDiagnostics(int step);

void saveSnapshot(int step);


// ============================================================
// main
// ============================================================

int main(
    int argc,
    char* argv[]
)
{
    int maxSteps =
        stepT3;


    if(argc >= 2)
    {
        maxSteps =
            atoi(
                argv[1]
            );
    }


    cout
        << scientific
        << setprecision(12);


    cout
        << "=============================================\n";

    cout
        << "Zhan 2023 2D Rayleigh-Taylor instability\n";

    cout
        << "ACE + MDF-LBM\n";

    cout
        << "Velocity gradient in Eq.(37): isotropic Eq.(41)\n";

    cout
        << "ACE interface normal: local Eq.(43), corrected sign\n";

    cout
        << "Re = "
        << Reynolds
        << "\n";

    cout
        << "At = "
        << Atwood
        << "\n";

    cout
        << "Pe = "
        << Peclet
        << "\n";

    cout
        << "=============================================\n";


    cout << "NX              = " << NX << "\n";
    cout << "NY              = " << NY << "\n";
    cout << "d               = " << d << "\n";
    cout << "dx              = " << dx << "\n";
    cout << "dt              = " << dt << "\n";
    cout << "cs2             = " << cs2 << "\n";
    cout << "W               = " << W << "\n";

    cout
        << "Ug=sqrt(gd)     = "
        << Ug
        << "\n";

    cout
        << "gravity g       = "
        << gravity
        << "\n";

    cout << "rhoA            = " << rhoA << "\n";
    cout << "rhoB            = " << rhoB << "\n";

    cout << "muA             = " << muA << "\n";
    cout << "muB             = " << muB << "\n";

    cout << "nuA             = " << nuA << "\n";
    cout << "nuB             = " << nuB << "\n";

    cout
        << "Ca (assumed)    = "
        << Capillary
        << "\n";

    cout << "sigma           = " << sigma << "\n";
    cout << "betaPF          = " << betaPF << "\n";
    cout << "kappaPF         = " << kappaPF << "\n";

    cout << "Mphi            = " << mobility << "\n";

    cout << "sPhase1         = " << sPhase1 << "\n";

    cout << "sPhase2         = " << sPhase2 << "\n";
    
    cout << "flow other rates = 1.0\n";
    cout
        << "RTI time scale  = "
        << rtiTimeScale
        << "\n";

    cout
        << "step(T=1)       = "
        << stepT1
        << "\n";

    cout
        << "step(T=1.75)    = "
        << stepT175
        << "\n";

    cout
        << "step(T=2.5)     = "
        << stepT25
        << "\n";

    cout
        << "step(T=3)       = "
        << stepT3
        << "\n";

    cout
        << "maxSteps        = "
        << maxSteps
        << "\n";

    cout
        << "=============================================\n\n";


    checkMatrices();


    initializeMacroscopicFields();

    computeDensityViscosity();

    computePhaseDerivatives();

    computeGradRho();

    computeSurfaceForce();


    initializePhaseDistribution();

    initializeFlowDistribution();


    recoverPhi();

    computeDensityViscosity();

    computePhaseDerivatives();

    computeGradRho();

    computeSurfaceForce();

    recoverVelocity();

    recoverPressure();

    computeVelocityGradientsIsotropic();


    initializeOldFields();


    initialPhaseMass =
        0.0;


    for(int y = 0;
        y < NY;
        y++)
    {
        for(int x = 0;
            x < NX;
            x++)
        {
            initialPhaseMass
                +=
            phi[y][x];
        }
    }


    printDiagnostics(0);

    saveSnapshot(0);


    for(int step = 1;
        step <= maxSteps;
        step++)
    {
        collidePhase();

        collideFlow();


        saveCurrentOldFields();


        streamPhase();

        streamFlow();


        recoverPhi();

        computeDensityViscosity();

        computePhaseDerivatives();

        computeGradRho();

        computeSurfaceForce();

        recoverVelocity();

        recoverPressure();

        computeVelocityGradientsIsotropic();


        if(checkNaN(step))
        {
            cout
                << "\nSimulation stopped because NaN / Inf was detected.\n";

            return 1;
        }


        if(
            step % 100 == 0
            ||
            step == maxSteps
        )
        {
            printDiagnostics(step);
        }


        if(
            step == stepT1
            ||
            step == stepT175
            ||
            step == stepT25
            ||
            step == stepT3
            ||
            step == maxSteps
        )
        {
            saveSnapshot(step);
        }
    }


    cout
        << "\n=============================================\n";

    cout
        << "Simulation finished.\n";

    cout
        << "=============================================\n";


    return 0;
}


// ============================================================
// Check MRT matrices
// ============================================================

void checkMatrices()
{
    double error5 = 0.0;

    double error9 = 0.0;


    for(int i = 0;
        i < Q5;
        i++)
    {
        for(int j = 0;
            j < Q5;
            j++)
        {
            double value = 0.0;


            for(int k = 0;
                k < Q5;
                k++)
            {
                value
                    +=
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


            error5 =
                max(
                    error5,
                    fabs(
                        value - exact
                    )
                );
        }
    }


    for(int i = 0;
        i < Q9;
        i++)
    {
        for(int j = 0;
            j < Q9;
            j++)
        {
            double value = 0.0;


            for(int k = 0;
                k < Q9;
                k++)
            {
                value
                    +=
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


            error9 =
                max(
                    error9,
                    fabs(
                        value - exact
                    )
                );
        }
    }


    cout
        << "max |M5*Minv5-I| = "
        << error5
        << "\n";


    cout
        << "max |M9*Minv9-I| = "
        << error9
        << "\n\n";
}


// ============================================================
// Initial field
//
// Zhan Eq.(57)
//
// h = -0.1 d cos(2 pi x / d)
//
// Heavy fluid A above
// Light fluid B below
// ============================================================

void initializeMacroscopicFields()
{
    for(int y = 0;
        y < NY;
        y++)
    {
        const double yy =
            -2.0
            *
            d

            +

            (y + 0.5)
            *
            dx;


        for(int x = 0;
            x < NX;
            x++)
        {
            const double xx =
                -0.5
                *
                d

                +

                (x + 0.5)
                *
                dx;


            const double hInterface =
                -0.1
                *
                d

                *

                cos(
                    2.0
                    *
                    PI
                    *
                    xx
                    /
                    d
                );


            phi[y][x] =
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
                    (
                        yy
                        -
                        hInterface
                    )
                    /
                    W
                );


            ux[y][x] =
                0.0;


            uy[y][x] =
                0.0;


            pressure[y][x] =
                0.0;
        }
    }
}


// ============================================================
// Density and viscosity
//
// NOTE:
// phi itself is NOT clipped.
//
// Only property interpolation uses phi limited to [0,1].
// ============================================================

void computeDensityViscosity()
{
    for(int y = 0;
        y < NY;
        y++)
    {
        for(int x = 0;
            x < NX;
            x++)
        {
            double phiProperty =
                phi[y][x];


            if(phiProperty > phiA)
            {
                phiProperty =
                    phiA;
            }


            if(phiProperty < phiB)
            {
                phiProperty =
                    phiB;
            }


            rho[y][x] =
                (rhoA - rhoB)
                /
                (phiA - phiB)

                *

                (
                    phiProperty
                    -
                    phiB
                )

                +

                rhoB;


            // equal dynamic viscosity
            muDyn[y][x] =
                muA;


            nu[y][x] =
                muDyn[y][x]
                /
                rho[y][x];
        }
    }
}


// ============================================================
// grad(phi), laplacian(phi), chemical potential
//
// x periodic
// y no-flux
// ============================================================

void computePhaseDerivatives()
{
    for(int y = 0;
        y < NY;
        y++)
    {
        for(int x = 0;
            x < NX;
            x++)
        {
            double gx = 0.0;

            double gy = 0.0;

            double lap = 0.0;


            for(int q = 1;
                q < Q9;
                q++)
            {
                int xn =
                    x
                    +
                    cx9[q];


                if(xn < 0)
                {
                    xn += NX;
                }


                if(xn >= NX)
                {
                    xn -= NX;
                }


                int yn =
                    y
                    +
                    cy9[q];


                if(yn < 0)
                {
                    yn =
                        0;
                }


                if(yn >= NY)
                {
                    yn =
                        NY - 1;
                }


                double cix =
                    c
                    *
                    cx9[q];


                double ciy =
                    c
                    *
                    cy9[q];


                double phin =
                    phi[yn][xn];


                gx
                    +=
                w9[q]
                *
                cix
                *
                phin;


                gy
                    +=
                w9[q]
                *
                ciy
                *
                phin;


                lap
                    +=
                2.0
                *
                w9[q]
                *
                (
                    phin
                    -
                    phi[y][x]
                );
            }


            gradPhiX[y][x] =
                gx
                /
                (
                    cs2
                    *
                    dt
                );


            gradPhiY[y][x] =
                gy
                /
                (
                    cs2
                    *
                    dt
                );


            lap =
                lap
                /
                (
                    cs2
                    *
                    dt
                    *
                    dt
                );


            chemical[y][x] =
                4.0
                *
                betaPF
                *
                (
                    phi[y][x]
                    -
                    phiA
                )
                *
                (
                    phi[y][x]
                    -
                    phiB
                )
                *
                (
                    phi[y][x]

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
// grad(rho)
// ============================================================

void computeGradRho()
{
    for(int y = 0;
        y < NY;
        y++)
    {
        for(int x = 0;
            x < NX;
            x++)
        {
            double gx = 0.0;

            double gy = 0.0;


            for(int q = 1;
                q < Q9;
                q++)
            {
                int xn =
                    x
                    +
                    cx9[q];


                if(xn < 0)
                {
                    xn += NX;
                }


                if(xn >= NX)
                {
                    xn -= NX;
                }


                int yn =
                    y
                    +
                    cy9[q];


                if(yn < 0)
                {
                    yn = 0;
                }


                if(yn >= NY)
                {
                    yn =
                        NY - 1;
                }


                double cix =
                    c
                    *
                    cx9[q];


                double ciy =
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


            gradRhoX[y][x] =
                gx
                /
                (
                    cs2
                    *
                    dt
                );


            gradRhoY[y][x] =
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
// Surface force + gravity
//
// G = -[0, rho-rhoRef] g
// ============================================================

void computeSurfaceForce()
{
    for(int y = 0;
        y < NY;
        y++)
    {
        for(int x = 0;
            x < NX;
            x++)
        {
            const double FsX =
                chemical[y][x]
                *
                gradPhiX[y][x];


            const double FsY =
                chemical[y][x]
                *
                gradPhiY[y][x];


            const double gravityY =
                -(
                    rho[y][x]
                    -
                    rhoReference
                )
                *
                gravity;


            forceX[y][x] =
                FsX;


            forceY[y][x] =
                FsY
                +
                gravityY;
        }
    }
}


// ============================================================
// Initialize phase distribution
// ============================================================

void initializePhaseDistribution()
{
    for(int y = 0;
        y < NY;
        y++)
    {
        for(int x = 0;
            x < NX;
            x++)
        {
            for(int q = 0;
                q < Q5;
                q++)
            {
                double cix =
                    c
                    *
                    cx5[q];


                double ciy =
                    c
                    *
                    cy5[q];


                if(q == 0)
                {
                    gPhase[y][x][q] =
                        phi[y][x]

                        +

                        (
                            w5[q]
                            -
                            1.0
                        )
                        *
                        eta
                        *
                        phi[y][x];
                }
                else
                {
                    gPhase[y][x][q] =
                        w5[q]
                        *
                        eta
                        *
                        phi[y][x]

                        +

                        w5[q]
                        *
                        (
                            cix
                            *
                            phi[y][x]
                            *
                            ux[y][x]

                            +

                            ciy
                            *
                            phi[y][x]
                            *
                            uy[y][x]
                        )
                        /
                        cs2;
                }
            }
        }
    }
}


// ============================================================
// Initialize flow distributions
// ============================================================

void initializeFlowDistribution()
{
    for(int y = 0;
        y < NY;
        y++)
    {
        for(int x = 0;
            x < NX;
            x++)
        {
            double r =
                rho[y][x];

            double u =
                ux[y][x];

            double v =
                uy[y][x];

            double p =
                pressure[y][x];


            for(int q = 0;
                q < Q9;
                q++)
            {
                double cix =
                    c
                    *
                    cx9[q];


                double ciy =
                    c
                    *
                    cy9[q];


                double feq1 =
                    w9[q]
                    *
                    (
                        r
                        *
                        u

                        +

                        (
                            cix
                            *
                            (
                                r
                                *
                                u
                                *
                                u

                                +

                                p
                            )

                            +

                            ciy
                            *
                            (
                                r
                                *
                                u
                                *
                                v
                            )
                        )
                        /
                        cs2
                    );


                double feq2 =
                    w9[q]
                    *
                    (
                        r
                        *
                        v

                        +

                        (
                            cix
                            *
                            (
                                r
                                *
                                v
                                *
                                u
                            )

                            +

                            ciy
                            *
                            (
                                r
                                *
                                v
                                *
                                v

                                +

                                p
                            )
                        )
                        /
                        cs2
                    );


                double Fi1 =
                    w9[q]
                    *
                    forceX[y][x];


                double Fi2 =
                    w9[q]
                    *
                    forceY[y][x];


                g1[y][x][q] =
                    feq1
                    -
                    0.5
                    *
                    dt
                    *
                    Fi1;


                g2[y][x][q] =
                    feq2
                    -
                    0.5
                    *
                    dt
                    *
                    Fi2;
            }
        }
    }
}


// ============================================================
// Save previous-time fields
// ============================================================

void initializeOldFields()
{
    saveCurrentOldFields();
}


void saveCurrentOldFields()
{
    for(int y = 0;
        y < NY;
        y++)
    {
        for(int x = 0;
            x < NX;
            x++)
        {
            oldPhiUx[y][x] =
                phi[y][x]
                *
                ux[y][x];


            oldPhiUy[y][x] =
                phi[y][x]
                *
                uy[y][x];


            oldQxx[y][x] =
                rho[y][x]
                *
                ux[y][x]
                *
                ux[y][x]

                +

                pressure[y][x];


            oldQxy[y][x] =
                rho[y][x]
                *
                ux[y][x]
                *
                uy[y][x];


            oldQyx[y][x] =
                rho[y][x]
                *
                uy[y][x]
                *
                ux[y][x];


            oldQyy[y][x] =
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
// Recover phi
// ============================================================

void recoverPhi()
{
    for(int y = 0;
        y < NY;
        y++)
    {
        for(int x = 0;
            x < NX;
            x++)
        {
            double value =
                0.0;


            for(int q = 0;
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
// Recover velocity
// ============================================================

void recoverVelocity()
{
    for(int y = 0;
        y < NY;
        y++)
    {
        for(int x = 0;
            x < NX;
            x++)
        {
            double sumX =
                0.0;


            double sumY =
                0.0;


            for(int q = 0;
                q < Q9;
                q++)
            {
                sumX
                    +=
                g1[y][x][q];


                sumY
                    +=
                g2[y][x][q];
            }


            ux[y][x] =
                (
                    sumX

                    +

                    0.5
                    *
                    dt
                    *
                    forceX[y][x]
                )
                /
                rho[y][x];


            uy[y][x] =
                (
                    sumY

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
// Recover pressure
// ============================================================

void recoverPressure()
{
    for(int y = 0;
        y < NY;
        y++)
    {
        for(int x = 0;
            x < NX;
            x++)
        {
            double momentX =
                0.0;


            double momentY =
                0.0;


            for(int q = 0;
                q < Q9;
                q++)
            {
                momentX
                    +=
                c
                *
                cx9[q]
                *
                g1[y][x][q];


                momentY
                    +=
                c
                *
                cy9[q]
                *
                g2[y][x][q];
            }


            pressure[y][x] =
                0.5
                *
                (
                    momentX

                    +

                    momentY

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
// Velocity gradient
//
// Zhan Eq.(41)
// ============================================================

void computeVelocityGradientsIsotropic()
{
    for(int y = 0;
        y < NY;
        y++)
    {
        for(int x = 0;
            x < NX;
            x++)
        {
            double ux_x = 0.0;

            double ux_y = 0.0;

            double uy_x = 0.0;

            double uy_y = 0.0;


            for(int q = 1;
                q < Q9;
                q++)
            {
                int xn =
                    x
                    +
                    cx9[q];


                if(xn < 0)
                {
                    xn += NX;
                }


                if(xn >= NX)
                {
                    xn -= NX;
                }


                int yn =
                    y
                    +
                    cy9[q];


                if(yn < 0)
                {
                    yn = 0;
                }


                if(yn >= NY)
                {
                    yn = NY - 1;
                }


                double cix =
                    c
                    *
                    cx9[q];


                double ciy =
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


            dudx[y][x] =
                ux_x
                /
                (
                    cs2
                    *
                    dt
                );


            dudy[y][x] =
                ux_y
                /
                (
                    cs2
                    *
                    dt
                );


            dvdx[y][x] =
                uy_x
                /
                (
                    cs2
                    *
                    dt
                );


            dvdy[y][x] =
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
// ACE collision
//
// local normal:
// Zhan Eq.(43)
// ============================================================

void collidePhase()
{
    for(int y = 0;
        y < NY;
        y++)
    {
        for(int x = 0;
            x < NX;
            x++)
        {
            double m[Q5] =
                {0};


            double meq[Q5] =
                {0};


            double mG[Q5] =
                {0};


            double mPost[Q5] =
                {0};


            double feq[Q5];

            double Gi[Q5];


            for(int k = 0;
                k < Q5;
                k++)
            {
                for(int q = 0;
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


            double dPhiUx =
                (
                    phi[y][x]
                    *
                    ux[y][x]

                    -

                    oldPhiUx[y][x]
                )
                /
                dt;


            double dPhiUy =
                (
                    phi[y][x]
                    *
                    uy[y][x]

                    -

                    oldPhiUy[y][x]
                )
                /
                dt;


            double lambdaACE =
                sqrt(
                    2.0
                    *
                    betaPF
                    /
                    kappaPF
                )

                *

                (
                    phiA
                    -
                    phi[y][x]
                )

                *

                (
                    phi[y][x]
                    -
                    phiB
                );


            // ------------------------------------------------
            // Eq.(43)
            // ------------------------------------------------

            double Cx = 0.0;

            double Cy = 0.0;


            for(int q = 0;
                q < Q5;
                q++)
            {
                double cix =
                    c
                    *
                    cx5[q];


                double ciy =
                    c
                    *
                    cy5[q];


                Cx
                    +=
                cix
                *
                gPhase[y][x][q];


                Cy
                    +=
                ciy
                *
                gPhase[y][x][q];
            }


            Cx =
                Cx

                -

                phi[y][x]
                *
                ux[y][x]

                +

                0.5
                *
                dt
                *
                dPhiUx;


            Cy =
                Cy

                -

                phi[y][x]
                *
                uy[y][x]

                +

                0.5
                *
                dt
                *
                dPhiUy;


            double A =
                -mobility

                -

                0.5
                *
                eta
                *
                cs2
                *
                dt;


            double B =
                mobility
                *
                lambdaACE
                /
                eta;


            double Cmag =
                sqrt(
                    Cx
                    *
                    Cx

                    +

                    Cy
                    *
                    Cy
                );


            // CORRECT SIGN
            double gradPhiMagnitudeLocal =
                (
                    -Cmag
                    -
                    B
                )
                /
                A;


            double nx = 0.0;

            double ny = 0.0;


            if(
                gradPhiMagnitudeLocal
                >
                normalEpsilon
            )
            {
                double denominator =
                    A

                    +

                    B
                    /
                    gradPhiMagnitudeLocal;


                if(
                    fabs(denominator)
                    >
                    normalEpsilon
                )
                {
                    double gradPhiXLocal =
                        Cx
                        /
                        denominator;


                    double gradPhiYLocal =
                        Cy
                        /
                        denominator;


                    double localNorm =
                        sqrt(
                            gradPhiXLocal
                            *
                            gradPhiXLocal

                            +

                            gradPhiYLocal
                            *
                            gradPhiYLocal
                        );


                    if(
                        localNorm
                        >
                        normalEpsilon
                    )
                    {
                        nx =
                            gradPhiXLocal
                            /
                            localNorm;


                        ny =
                            gradPhiYLocal
                            /
                            localNorm;
                    }
                }
            }


            for(int q = 0;
                q < Q5;
                q++)
            {
                double cix =
                    c
                    *
                    cx5[q];


                double ciy =
                    c
                    *
                    cy5[q];


                if(q == 0)
                {
                    feq[q] =
                        phi[y][x]

                        +

                        (
                            w5[q]
                            -
                            1.0
                        )
                        *
                        eta
                        *
                        phi[y][x];
                }
                else
                {
                    feq[q] =
                        w5[q]
                        *
                        eta
                        *
                        phi[y][x]

                        +

                        w5[q]
                        *
                        (
                            cix
                            *
                            phi[y][x]
                            *
                            ux[y][x]

                            +

                            ciy
                            *
                            phi[y][x]
                            *
                            uy[y][x]
                        )
                        /
                        cs2;
                }


                Gi[q] =
                    w5[q]
                    *
                    (
                        (
                            cix
                            *
                            dPhiUx

                            +

                            ciy
                            *
                            dPhiUy
                        )
                        /
                        cs2

                        +

                        eta
                        *
                        lambdaACE
                        *
                        (
                            cix
                            *
                            nx

                            +

                            ciy
                            *
                            ny
                        )
                    );
            }


            for(int k = 0;
                k < Q5;
                k++)
            {
                for(int q = 0;
                    q < Q5;
                    q++)
                {
                    meq[k]
                        +=
                    M5[k][q]
                    *
                    feq[q];


                    mG[k]
                        +=
                    M5[k][q]
                    *
                    Gi[q];
                }


                double s =
                    Sphase[k];


                mPost[k] =
                    (1.0 - s)
                    *
                    m[k]

                    +

                    s
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
                        s
                    )
                    *
                    mG[k];
            }


            for(int q = 0;
                q < Q5;
                q++)
            {
                double value =
                    0.0;


                for(int k = 0;
                    k < Q5;
                    k++)
                {
                    value
                        +=
                    Minv5[q][k]
                    *
                    mPost[k];
                }


                gPhasePost[y][x][q] =
                    value;
            }
        }
    }
}


// ============================================================
// Flow collision
// ============================================================

void collideFlow()
{
    for(int y = 0;
        y < NY;
        y++)
    {
        for(int x = 0;
            x < NX;
            x++)
        {
            double r =
                rho[y][x];


            double u =
                ux[y][x];


            double v =
                uy[y][x];


            double p =
                pressure[y][x];


            double Qxx =
                r
                *
                u
                *
                u

                +

                p;


            double Qxy =
                r
                *
                u
                *
                v;


            double Qyx =
                r
                *
                v
                *
                u;


            double Qyy =
                r
                *
                v
                *
                v

                +

                p;


            double dQxx =
                (
                    Qxx
                    -
                    oldQxx[y][x]
                )
                /
                dt;


            double dQxy =
                (
                    Qxy
                    -
                    oldQxy[y][x]
                )
                /
                dt;


            double dQyx =
                (
                    Qyx
                    -
                    oldQyx[y][x]
                )
                /
                dt;


            double dQyy =
                (
                    Qyy
                    -
                    oldQyy[y][x]
                )
                /
                dt;


            double feq1[Q9];

            double feq2[Q9];


            double Fi1[Q9];

            double Fi2[Q9];


            double Gi1[Q9];

            double Gi2[Q9];


            for(int q = 0;
                q < Q9;
                q++)
            {
                double cix =
                    c
                    *
                    cx9[q];


                double ciy =
                    c
                    *
                    cy9[q];


                feq1[q] =
                    w9[q]
                    *
                    (
                        r
                        *
                        u

                        +

                        (
                            cix
                            *
                            (
                                r
                                *
                                u
                                *
                                u

                                +

                                p
                            )

                            +

                            ciy
                            *
                            (
                                r
                                *
                                u
                                *
                                v
                            )
                        )
                        /
                        cs2
                    );


                feq2[q] =
                    w9[q]
                    *
                    (
                        r
                        *
                        v

                        +

                        (
                            cix
                            *
                            (
                                r
                                *
                                v
                                *
                                u
                            )

                            +

                            ciy
                            *
                            (
                                r
                                *
                                v
                                *
                                v

                                +

                                p
                            )
                        )
                        /
                        cs2
                    );


                Fi1[q] =
                    w9[q]
                    *
                    forceX[y][x];


                Fi2[q] =
                    w9[q]
                    *
                    forceY[y][x];


                Gi1[q] =
                    w9[q]
                    *
                    (
                        (
                            cix
                            *
                            dQxx

                            +

                            ciy
                            *
                            dQxy
                        )
                        /
                        cs2

                        +

                        u
                        *
                        (
                            cix
                            *
                            gradRhoX[y][x]

                            +

                            ciy
                            *
                            gradRhoY[y][x]
                        )

                        -

                        r
                        *
                        (
                            cix
                            *
                            dudx[y][x]

                            +

                            ciy
                            *
                            dvdx[y][x]
                        )
                    );


                Gi2[q] =
                    w9[q]
                    *
                    (
                        (
                            cix
                            *
                            dQyx

                            +

                            ciy
                            *
                            dQyy
                        )
                        /
                        cs2

                        +

                        v
                        *
                        (
                            cix
                            *
                            gradRhoX[y][x]

                            +

                            ciy
                            *
                            gradRhoY[y][x]
                        )

                        -

                        r
                        *
                        (
                            cix
                            *
                            dudy[y][x]

                            +

                            ciy
                            *
                            dvdy[y][x]
                        )
                    );
            }


            double m1[Q9] = {0};

            double m2[Q9] = {0};

            double meq1[Q9] = {0};

            double meq2[Q9] = {0};

            double mF1[Q9] = {0};

            double mF2[Q9] = {0};

            double mG1[Q9] = {0};

            double mG2[Q9] = {0};

            double mPost1[Q9] = {0};

            double mPost2[Q9] = {0};


            for(int k = 0;
                k < Q9;
                k++)
            {
                for(int q = 0;
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


            double localS1 =
                1.0
                /
                (
                    nu[y][x]
                    /
                    (
                        cs2
                        *
                        dt
                    )

                    +

                    0.5
                );


            double Sflow[Q9] =
            {
                su0,

                localS1,
                localS1,

                su2,
                su2,
                su2,

                su3,
                su3,

                su4
            };


            for(int k = 0;
                k < Q9;
                k++)
            {
                double s =
                    Sflow[k];


                mPost1[k] =
                    (1.0 - s)
                    *
                    m1[k]

                    +

                    s
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
                        s
                    )
                    *
                    (
                        mF1[k]
                        +
                        mG1[k]
                    );


                mPost2[k] =
                    (1.0 - s)
                    *
                    m2[k]

                    +

                    s
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
                        s
                    )
                    *
                    (
                        mF2[k]
                        +
                        mG2[k]
                    );
            }


            for(int q = 0;
                q < Q9;
                q++)
            {
                double value1 =
                    0.0;


                double value2 =
                    0.0;


                for(int k = 0;
                    k < Q9;
                    k++)
                {
                    value1
                        +=
                    Minv9[q][k]
                    *
                    mPost1[k];


                    value2
                        +=
                    Minv9[q][k]
                    *
                    mPost2[k];
                }


                g1Post[y][x][q] =
                    value1;


                g2Post[y][x][q] =
                    value2;
            }
        }
    }
}


// ============================================================
// Phase streaming
//
// x periodic
// y no-flux
// ============================================================

void streamPhase()
{
    for(int y = 0;
        y < NY;
        y++)
    {
        for(int x = 0;
            x < NX;
            x++)
        {
            for(int q = 0;
                q < Q5;
                q++)
            {
                gPhaseNext[y][x][q] =
                    0.0;
            }
        }
    }


    for(int y = 0;
        y < NY;
        y++)
    {
        for(int x = 0;
            x < NX;
            x++)
        {
            for(int q = 0;
                q < Q5;
                q++)
            {
                int xn =
                    x
                    +
                    cx5[q];


                if(xn < 0)
                {
                    xn += NX;
                }


                if(xn >= NX)
                {
                    xn -= NX;
                }


                int yn =
                    y
                    +
                    cy5[q];


                if(
                    yn >= 0
                    &&
                    yn < NY
                )
                {
                    gPhaseNext[yn][xn][q] =
                        gPhasePost[y][x][q];
                }
                else
                {
                    gPhaseNext[y][x][opp5[q]] =
                        gPhasePost[y][x][q];
                }
            }
        }
    }


    for(int y = 0;
        y < NY;
        y++)
    {
        for(int x = 0;
            x < NX;
            x++)
        {
            for(int q = 0;
                q < Q5;
                q++)
            {
                gPhase[y][x][q] =
                    gPhaseNext[y][x][q];
            }
        }
    }
}


// ============================================================
// Flow streaming
//
// x periodic
// top/bottom no-slip
// ============================================================

void streamFlow()
{
    for(int y = 0;
        y < NY;
        y++)
    {
        for(int x = 0;
            x < NX;
            x++)
        {
            for(int q = 0;
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


    for(int y = 0;
        y < NY;
        y++)
    {
        for(int x = 0;
            x < NX;
            x++)
        {
            for(int q = 0;
                q < Q9;
                q++)
            {
                int xn =
                    x
                    +
                    cx9[q];


                if(xn < 0)
                {
                    xn += NX;
                }


                if(xn >= NX)
                {
                    xn -= NX;
                }


                const int yn =
                    y
                    +
                    cy9[q];


                if(
                    yn >= 0
                    &&
                    yn < NY
                )
                {
                    g1Next[yn][xn][q] =
                        g1Post[y][x][q];


                    g2Next[yn][xn][q] =
                        g2Post[y][x][q];
                }
                else
                {
                    const int qo =
                        opp9[q];


                    g1Next[y][x][qo] =
                        -g1Post[y][x][q]

                        -

                        dt
                        *
                        w9[q]
                        *
                        forceX[y][x];


                    g2Next[y][x][qo] =
                        -g2Post[y][x][q]

                        -

                        dt
                        *
                        w9[q]
                        *
                        forceY[y][x];
                }
            }
        }
    }


    for(int y = 0;
        y < NY;
        y++)
    {
        for(int x = 0;
            x < NX;
            x++)
        {
            for(int q = 0;
                q < Q9;
                q++)
            {
                g1[y][x][q] =
                    g1Next[y][x][q];


                g2[y][x][q] =
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
    for(int y = 0;
        y < NY;
        y++)
    {
        for(int x = 0;
            x < NX;
            x++)
        {
            if(
                !isfinite(phi[y][x])
                ||
                !isfinite(rho[y][x])
                ||
                !isfinite(nu[y][x])
                ||
                !isfinite(chemical[y][x])
                ||
                !isfinite(ux[y][x])
                ||
                !isfinite(uy[y][x])
                ||
                !isfinite(pressure[y][x])
            )
            {
                const double xx =
                    -0.5
                    *
                    d

                    +

                    (x + 0.5)
                    *
                    dx;


                const double yy =
                    -2.0
                    *
                    d

                    +

                    (y + 0.5)
                    *
                    dx;


                cout
                    << "\n=============================================\n";

                cout
                    << "NaN / Inf detected\n";

                cout
                    << "step = "
                    << step
                    << "\n";

                cout
                    << "T = "
                    << (
                        step
                        *
                        dt
                        /
                        rtiTimeScale
                    )
                    << "\n";

                cout
                    << "x index = "
                    << x
                    << "\n";

                cout
                    << "y index = "
                    << y
                    << "\n";

                cout
                    << "physical x = "
                    << xx
                    << "\n";

                cout
                    << "physical y = "
                    << yy
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
                    << "nu = "
                    << nu[y][x]
                    << "\n";

                cout
                    << "chemical = "
                    << chemical[y][x]
                    << "\n";

                cout
                    << "ux = "
                    << ux[y][x]
                    << "\n";

                cout
                    << "uy = "
                    << uy[y][x]
                    << "\n";

                cout
                    << "pressure = "
                    << pressure[y][x]
                    << "\n";

                cout
                    << "=============================================\n";


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
    double phiMin =
        phi[0][0];

    double phiMax =
        phi[0][0];

    double rhoMin =
        rho[0][0];

    double rhoMax =
        rho[0][0];

    double nuMin =
        nu[0][0];

    double nuMax =
        nu[0][0];


    double maxSpeed = 0.0;

    double maxChemical = 0.0;

    double maxGradPhi = 0.0;

    double maxLambdaACE = 0.0;

    double maxForce = 0.0;

    double phaseMass = 0.0;


    for(int y = 0;
        y < NY;
        y++)
    {
        for(int x = 0;
            x < NX;
            x++)
        {
            phiMin =
                min(
                    phiMin,
                    phi[y][x]
                );


            phiMax =
                max(
                    phiMax,
                    phi[y][x]
                );


            rhoMin =
                min(
                    rhoMin,
                    rho[y][x]
                );


            rhoMax =
                max(
                    rhoMax,
                    rho[y][x]
                );


            nuMin =
                min(
                    nuMin,
                    nu[y][x]
                );


            nuMax =
                max(
                    nuMax,
                    nu[y][x]
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


            maxChemical =
                max(
                    maxChemical,
                    fabs(
                        chemical[y][x]
                    )
                );


            const double gradPhiMagnitude =
                sqrt(
                    gradPhiX[y][x]
                    *
                    gradPhiX[y][x]

                    +

                    gradPhiY[y][x]
                    *
                    gradPhiY[y][x]
                );


            maxGradPhi =
                max(
                    maxGradPhi,
                    gradPhiMagnitude
                );


            const double lambdaACE =
                sqrt(
                    2.0
                    *
                    betaPF
                    /
                    kappaPF
                )

                *

                (
                    phiA
                    -
                    phi[y][x]
                )

                *

                (
                    phi[y][x]
                    -
                    phiB
                );


            maxLambdaACE =
                max(
                    maxLambdaACE,
                    fabs(
                        lambdaACE
                    )
                );


            const double ff =
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
                    ff
                );


            phaseMass
                +=
            phi[y][x];
        }
    }


    const double relativeMassChange =
        (
            phaseMass
            -
            initialPhaseMass
        )
        /
        initialPhaseMass;


    const double T =
        step
        *
        dt
        /
        rtiTimeScale;


    cout
        << "\n---------------------------------------------\n";

    cout
        << "step              = "
        << step
        << "\n";

    cout
        << "T = t/t_RTI       = "
        << T
        << "\n";

    cout
        << "phi range         = ["
        << phiMin
        << ", "
        << phiMax
        << "]\n";

    cout
        << "rho range         = ["
        << rhoMin
        << ", "
        << rhoMax
        << "]\n";

    cout
        << "nu range          = ["
        << nuMin
        << ", "
        << nuMax
        << "]\n";

    cout
        << "max |chemical|    = "
        << maxChemical
        << "\n";

    cout
        << "max |grad phi|    = "
        << maxGradPhi
        << "\n";

    cout
        << "max |lambda|      = "
        << maxLambdaACE
        << "\n";

    cout
        << "max |force|       = "
        << maxForce
        << "\n";

    cout
        << "max |u|           = "
        << maxSpeed
        << "\n";

    cout
        << "relative phi mass = "
        << relativeMassChange
        << "\n";

    cout
        << "---------------------------------------------\n";
}


// ============================================================
// Save phase snapshots
// ============================================================

void saveSnapshot(
    int step
)
{
    char filename[160];


    sprintf(
        filename,
        "RTI_phi_ACE_Re90_step_%05d.dat",
        step
    );


    ofstream fout(
        filename
    );


    fout
        << scientific
        << setprecision(10);


    for(int y = 0;
        y < NY;
        y++)
    {
        for(int x = 0;
            x < NX;
            x++)
        {
            fout
                << phi[y][x];


            if(x < NX - 1)
            {
                fout
                    << " ";
            }
        }


        fout
            << "\n";
    }


    fout.close();


    cout
        << "Saved: "
        << filename
        << "\n";
}