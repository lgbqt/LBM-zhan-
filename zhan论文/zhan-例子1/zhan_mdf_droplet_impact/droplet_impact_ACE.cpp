#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <cstdio>
#include <cstdlib>

using namespace std;


// ============================================================
// Zhan 2023
//
// Droplet impact on a thin liquid film
//
// ACE + MDF-LBM
//
// rhoA : rhoB = 1000 : 1
// Re = 500
// We = 8000
//
// IMPORTANT:
//
// 1. Phase:
//    D2Q5 conservative Allen-Cahn equation (Zhan Eq.(9b))
//
// 2. Flow:
//    D2Q9 MDF-LBM
//
// 3. Velocity gradient used in Eq.(37):
//    Zhan Remark 4, Eq.(41)
//    second-order isotropic scheme
//
// 4. ACE interface normal n:
//    Zhan Eq.(43) local phase-gradient scheme
//
// 5. Mphi = 0.1:
//    implementation assumption
//
// 6. Top flow boundary:
//    open boundary, implemented by zero-gradient copy
// ============================================================


// ============================================================
// 1. Grid
// ============================================================

const int NX = 1500;
const int NY = 500;

const int Q5 = 5;
const int Q9 = 9;


const double dx = 1.0;
const double dt = 1.0;


const double c =
    dx / dt;


const double cs2 =
    c * c / 3.0;


// ============================================================
// 2. Physical domain
// ============================================================

const double L = 1500.0;
const double H = 500.0;


// Film height

const double h =
    0.1 * H;


// ============================================================
// 3. Droplet
// ============================================================

const double D = 200.0;

const double R =
    0.5 * D;


const double xc = 0.0;

const double yc =
    h + 0.5 * D;


// ============================================================
// 4. Phase field
// ============================================================

// liquid

const double phiA = 1.0;


// gas

const double phiB = 0.0;


// interface width

const double W = 5.0;


// ============================================================
// 5. Dimensionless parameters
// ============================================================

const double Reynolds = 500.0;

const double Weber = 8000.0;


// Impact velocity

const double V0 = 0.05;


// ============================================================
// 6. Density
//
// rhoA : rhoB = 1000 : 1
//
// normalization used here:
// rhoA = 1000
// rhoB = 1
// ============================================================

const double rhoA = 1000.0;

const double rhoB = 1.0;


// ============================================================
// 7. Viscosity
// ============================================================

// Paper gives muB = 0.2

const double muB = 0.2;


// Re = rhoA D V0 / muA

const double muA =

    rhoA
    *
    D
    *
    V0

    /

    Reynolds;


// kinematic viscosities

const double nuA =
    muA / rhoA;


const double nuB =
    muB / rhoB;


// ============================================================
// 8. Surface tension
//
// We = rhoA D V0^2 / sigma
// ============================================================

const double sigma =

    rhoA
    *
    D
    *
    V0
    *
    V0

    /

    Weber;


// ============================================================
// 9. Free-energy parameters
//
// Zhan Eq.(4)
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
// 10. ACE mobility
//
// IMPLEMENTATION ASSUMPTION
//
// Zhan Sec.4.4 does not explicitly give Mphi.
// ============================================================

const double mobility = 0.10;

const double eta = 1.0;


// Avoid division by zero in n = grad(phi) / |grad(phi)|.
// At points with a vanishing phase gradient, n is set to zero.

const double normalEpsilon = 1.0e-14;


// ============================================================
// 11. Phase relaxation
//
// Mphi = eta cs^2 (1/s1 - 1/2) dt
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


const double Sphase[Q5] =
{
    1.0,

    sPhase1,
    sPhase1,

    1.0,
    1.0
};


// ============================================================
// 12. Flow relaxation parameters
//
// su1 depends on local nu(phi)
//
// Zhan Sec.4.4:
//
// su2 = su3 = 1.2
// ============================================================

const double su0 = 1.0;

const double su2 = 1.2;

const double su3 = 1.2;

const double su4 = 1.0;


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
// 15. D2Q5 transformation matrix
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
// 16. D2Q9 transformation matrix
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
// 18. Old-time quantities
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
// 21. Diagnostics
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


bool checkNaN(
    int step
);


void printDiagnostics(
    int step
);


void saveSnapshot(
    int step
);


// ============================================================
// main
// ============================================================

int main(
    int argc,
    char* argv[]
)
{
    int maxSteps = 6400;


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
        << "Zhan 2023 droplet impact\n";

    cout
        << "ACE + MDF-LBM (Zhan 2023 Eq.(9b))\n";

    cout
        << "Velocity gradient in Eq.(37): isotropic Eq.(41)\n";

    cout
        << "ACE interface normal: local Eq.(43)\n";

    cout
        << "rhoA : rhoB = "
        << rhoA / rhoB
        << " : 1\n";

    cout
        << "Re = "
        << Reynolds
        << "\n";

    cout
        << "We = 8000\n";

    cout
        << "=============================================\n";


    cout
        << "NX              = "
        << NX
        << "\n";


    cout
        << "NY              = "
        << NY
        << "\n";


    cout
        << "dx              = "
        << dx
        << "\n";


    cout
        << "dt              = "
        << dt
        << "\n";


    cout
        << "c               = "
        << c
        << "\n";


    cout
        << "cs2             = "
        << cs2
        << "\n";


    cout
        << "h               = "
        << h
        << "\n";


    cout
        << "D               = "
        << D
        << "\n";


    cout
        << "W               = "
        << W
        << "\n";


    cout
        << "V0              = "
        << V0
        << "\n";


    cout
        << "rhoA            = "
        << rhoA
        << "\n";


    cout
        << "rhoB            = "
        << rhoB
        << "\n";


    cout
        << "muA             = "
        << muA
        << "\n";


    cout
        << "muB             = "
        << muB
        << "\n";


    cout
        << "nuA             = "
        << nuA
        << "\n";


    cout
        << "nuB             = "
        << nuB
        << "\n";


    cout
        << "sigma           = "
        << sigma
        << "\n";


    cout
        << "betaPF          = "
        << betaPF
        << "\n";


    cout
        << "kappaPF         = "
        << kappaPF
        << "\n";


    cout
        << "Mphi (assumed)  = "
        << mobility
        << "\n";


    cout
        << "sPhase1         = "
        << sPhase1
        << "\n";


    cout
        << "su2 = su3       = "
        << su2
        << "\n";


    cout
        << "maxSteps        = "
        << maxSteps
        << "\n";


    cout
        << "=============================================\n\n";


    // ========================================================
    // Check matrices
    // ========================================================

    checkMatrices();


    // ========================================================
    // Initialize macroscopic fields
    // ========================================================

    initializeMacroscopicFields();


    computeDensityViscosity();


    computePhaseDerivatives();


    computeGradRho();


    computeSurfaceForce();


    // ========================================================
    // Initialize distributions
    // ========================================================

    initializePhaseDistribution();


    initializeFlowDistribution();


    // ========================================================
    // Recover initial state from distributions
    // ========================================================

    recoverPhi();


    computeDensityViscosity();


    computePhaseDerivatives();


    computeGradRho();


    computeSurfaceForce();


    recoverVelocity();


    recoverPressure();


    computeVelocityGradientsIsotropic();


    // ========================================================
    // Initialize explicit old-time quantities
    // ========================================================

    initializeOldFields();


    // ========================================================
    // Initial phase mass
    // ========================================================

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


    printDiagnostics(
        0
    );


    saveSnapshot(
        0
    );


    // ========================================================
    // Main loop
    // ========================================================

    for(int step = 1;
        step <= maxSteps;
        step++)
    {
        // ----------------------------------------------------
        // Current macroscopic quantities are already known.
        //
        // Collision uses current n-level fields.
        // ----------------------------------------------------

        collidePhase();


        collideFlow();


        // ----------------------------------------------------
        // Save n-level values for explicit time derivatives
        // ----------------------------------------------------

        saveCurrentOldFields();


        // ----------------------------------------------------
        // Streaming
        // ----------------------------------------------------

        streamPhase();


        streamFlow();


        // ====================================================
        // Recover n+1 state
        // ====================================================

        recoverPhi();


        computeDensityViscosity();


        computePhaseDerivatives();


        computeGradRho();


        computeSurfaceForce();


        recoverVelocity();


        recoverPressure();


        // ----------------------------------------------------
        // IMPORTANT:
        //
        // Eq.(37) velocity gradient:
        // isotropic Eq.(41)
        // ----------------------------------------------------

        computeVelocityGradientsIsotropic();


        // ====================================================
        // NaN check EVERY STEP
        // ====================================================

        if(
            checkNaN(
                step
            )
        )
        {
            cout
                << "\nSimulation stopped because "
                << "NaN / Inf was detected.\n";


            return 1;
        }


        // ====================================================
        // Diagnostics
        // ====================================================

        if(
            step % 10 == 0
            ||
            step == maxSteps
        )
        {
            printDiagnostics(
                step
            );
        }


        // ====================================================
        // Save paper comparison times
        // ====================================================

        if(
            step == 400
            ||
            step == 800
            ||
            step == 1600
            ||
            step == 3200
            ||
            step == 6400
            ||
            step == maxSteps
        )
        {
            saveSnapshot(
                step
            );
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
// Check matrices
// ============================================================

void checkMatrices()
{
    double error5 =
        0.0;


    double error9 =
        0.0;


    for(int i = 0;
        i < Q5;
        i++)
    {
        for(int j = 0;
            j < Q5;
            j++)
        {
            double value =
                0.0;


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
                        value
                        -
                        exact
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
            double value =
                0.0;


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
                        value
                        -
                        exact
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
// Initial macroscopic fields
//
// Zhan Eq.(53)-(55)
// ============================================================

void initializeMacroscopicFields()
{
    for(int y = 0;
        y < NY;
        y++)
    {
        double yy =

            (y + 0.5)
            *
            dx;


        for(int x = 0;
            x < NX;
            x++)
        {
            double xx =

                -0.5
                *
                L

                +

                (x + 0.5)
                *
                dx;


            // =================================================
            // 1. Liquid-film profile
            //
            // Zhan Eq.(53)
            //
            // IMPORTANT:
            //
            // Do not cut this profile at y = h.
            // Let the diffuse film interface continue smoothly
            // into the gas region.
            // =================================================

            double phiFilm =

                0.5
                *
                (phiA + phiB)

                +

                0.5
                *
                (phiA - phiB)

                *
                tanh
                (
                    2.0
                    *
                    (h - yy)

                    /

                    W
                );


            // =================================================
            // 2. Droplet profile
            //
            // Zhan Eq.(54)
            // =================================================

            double rr =

                sqrt
                (
                    (xx - xc)
                    *
                    (xx - xc)

                    +

                    (yy - yc)
                    *
                    (yy - yc)
                );


            double phiDrop =

                0.5
                *
                (phiA + phiB)

                +

                0.5
                *
                (phiA - phiB)

                *
                tanh
                (
                    2.0
                    *
                    (
                        0.5 * D
                        -
                        rr
                    )

                    /

                    W
                );


            // =================================================
            // 3. Union of film and droplet
            //
            // IMPLEMENTATION CHOICE:
            //
            // phi = max(phiFilm, phiDrop)
            //
            // This preserves both diffuse interfaces and avoids
            // the artificial discontinuity at y = h.
            // =================================================

            phi[y][x] =

                max(
                    phiFilm,
                    phiDrop
                );


            // =================================================
            // Initial velocity
            //
            // Zhan Eq.(55)
            // =================================================

            ux[y][x] =
                0.0;


            if(yy > h)
            {
                // --------------------------------------------
                // Keep the paper's velocity initialization:
                //
                // v = -phi V0
                // --------------------------------------------

                uy[y][x] =

                    -phiDrop
                    *
                    V0;
            }
            else
            {
                uy[y][x] =
                    0.0;
            }


            pressure[y][x] =
                0.0;
        }
    }
}

// ============================================================
// rho(phi), nu(phi)
// ============================================================

// ============================================================
// Density and viscosity from phase field
//
// IMPORTANT:
//
// phi[y][x] itself is NOT clipped.
//
// Only the value used to calculate physical properties
// is limited to [phiB, phiA].
//
// This is a diagnostic test.
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
            // =================================================
            // Real phase field
            // =================================================

            double phiReal =
                phi[y][x];


            // =================================================
            // Diagnostic bounded phase field
            //
            // phiProperty is used ONLY for rho and nu.
            //
            // The real phi field is still used by:
            //
            // chemical potential
            // phase-field collision
            // surface force
            // gradients
            //
            // =================================================

            double phiProperty =
                phiReal;


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


            // =================================================
            // Density
            //
            // rho =
            //
            // (rhoA-rhoB)/(phiA-phiB)
            // *
            // (phi-phiB)
            // +
            // rhoB
            // =================================================

            rho[y][x] =

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
                    phiProperty
                    -
                    phiB
                )

                +

                rhoB;


            // =================================================
            // Kinematic viscosity
            //
            // nu =
            //
            // (nuA-nuB)/(phiA-phiB)
            // *
            // (phi-phiB)
            // +
            // nuB
            // =================================================

            nu[y][x] =

                (
                    nuA
                    -
                    nuB
                )

                /

                (
                    phiA
                    -
                    phiB
                )

                *

                (
                    phiProperty
                    -
                    phiB
                )

                +

                nuB;
        }
    }
}

// ============================================================
// grad(phi), lap(phi), chemical potential
//
// Eq.(41)-(42)
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
            double gx =
                0.0;


            double gy =
                0.0;


            double lap =
                0.0;


            for(int q = 1;
                q < Q9;
                q++)
            {
                int xn =
                    x + cx9[q];


                if(xn < 0)
                {
                    xn +=
                        NX;
                }


                if(xn >= NX)
                {
                    xn -=
                        NX;
                }


                int yn =
                    y + cy9[q];


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
                    (phiA + phiB)
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
            double gx =
                0.0;


            double gy =
                0.0;


            for(int q = 1;
                q < Q9;
                q++)
            {
                int xn =
                    x + cx9[q];


                if(xn < 0)
                {
                    xn +=
                        NX;
                }


                if(xn >= NX)
                {
                    xn -=
                        NX;
                }


                int yn =
                    y + cy9[q];


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
// Surface force
//
// Fs = mu_phi grad(phi)
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
            forceX[y][x] =

                chemical[y][x]
                *
                gradPhiX[y][x];


            forceY[y][x] =

                chemical[y][x]
                *
                gradPhiY[y][x];
        }
    }
}


// ============================================================
// Initialize ACE phase distribution (Zhan Eq.(9b))
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
//
// g = f - dt F / 2
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
// Initialize old fields
// ============================================================

void initializeOldFields()
{
    saveCurrentOldFields();
}


// ============================================================
// Save old fields
// ============================================================

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
//
// Zhan Eq.(35b)
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
//
// Zhan Eq.(32)
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
// Zhan 2023 Remark 4, Eq.(41)
//
// Used ONLY for velocity gradients in Eq.(37).
//
// Second-order isotropic scheme.
//
// x direction: periodic
// y direction: nearest-value ghost treatment
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
                int xn = x + cx9[q];

                if(xn < 0)
                {
                    xn += NX;
                }

                if(xn >= NX)
                {
                    xn -= NX;
                }

                int yn = y + cy9[q];

                if(yn < 0)
                {
                    yn = 0;
                }

                if(yn >= NY)
                {
                    yn = NY - 1;
                }

                double cix = c * cx9[q];
                double ciy = c * cy9[q];

                ux_x += w9[q] * cix * ux[yn][xn];
                ux_y += w9[q] * ciy * ux[yn][xn];

                uy_x += w9[q] * cix * uy[yn][xn];
                uy_y += w9[q] * ciy * uy[yn][xn];
            }

            dudx[y][x] = ux_x / (cs2 * dt);
            dudy[y][x] = ux_y / (cs2 * dt);

            dvdx[y][x] = uy_x / (cs2 * dt);
            dvdy[y][x] = uy_y / (cs2 * dt);
        }
    }
}



// ============================================================
// ACE phase collision (Zhan Eq.(9b))
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


            // =================================================
            // m = M g
            // =================================================

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


            // =================================================
            // Explicit d_t(phi u)
            // =================================================

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


            // =================================================
            // Conservative Allen-Cahn interface-restoring term
            //
            // lambda = sqrt(2 beta / kappa)
            //          (phiA - phi)(phi - phiB)
            //
            // IMPORTANT:
            // For the ACE normal vector n, use the local
            // phase-gradient scheme of Zhan Eq.(43).
            //
            // Chemical potential and surface force still use
            // the Eq.(41)-(42) finite-difference derivatives
            // computed in computePhaseDerivatives().
            // =================================================

            double lambdaACE =

                sqrt
                (
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


            // =================================================
            // Zhan Eq.(43)
            //
            // C_beta = sum_i c_i,beta g_i,0
            //          - phi u_beta
            //          + dt/2 * d_t(phi u_beta)
            //
            // Here g_i,0 is gPhase.
            // =================================================

            double Cx = 0.0;
            double Cy = 0.0;

            for(int q = 0;
                q < Q5;
                q++)
            {
                double cix = c * cx5[q];
                double ciy = c * cy5[q];

                Cx += cix * gPhase[y][x][q];
                Cy += ciy * gPhase[y][x][q];
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


            // =================================================
            // Eq.(43a,b)
            //
            // A = -M_phi - 0.5 eta cs^2 dt
            // B =  M_phi lambda / eta
            //
            // |grad phi| = -( |C| - B ) / A
            // grad phi   = C / ( A + B/|grad phi| )
            // =================================================

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

                sqrt
                (
                    Cx
                    *
                    Cx

                    +

                    Cy
                    *
                    Cy
                );

            double gradPhiMagnitudeLocal =
                (
                    -Cmag
                    -
                    B
                )

                /

                A;


            // =================================================
            // Local phase gradient and unit normal
            // =================================================

            double nx = 0.0;
            double ny = 0.0;

            if(gradPhiMagnitudeLocal > normalEpsilon)
            {
                double denominator =

                    A

                    +

                    B
                    /
                    gradPhiMagnitudeLocal;

                if(fabs(denominator) > normalEpsilon)
                {
                    double gradPhiXLocal = Cx / denominator;
                    double gradPhiYLocal = Cy / denominator;

                    double localNorm =

                        sqrt
                        (
                            gradPhiXLocal
                            *
                            gradPhiXLocal

                            +

                            gradPhiYLocal
                            *
                            gradPhiYLocal
                        );

                    if(localNorm > normalEpsilon)
                    {
                        nx = gradPhiXLocal / localNorm;
                        ny = gradPhiYLocal / localNorm;
                    }
                }
            }


            // =================================================
            // Equilibrium and auxiliary term
            // =================================================

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


            // =================================================
            // Moment transform
            // =================================================

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
                        0.5 * s
                    )

                    *
                    mG[k];
            }


            // =================================================
            // inverse transform
            // =================================================

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
//
// C = 0
//
// Eq.(37)
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


            // =================================================
            // Q tensors
            // =================================================

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


            // =================================================
            // Explicit temporal derivatives
            // =================================================

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


            // =================================================
            // Velocity-space construction
            // =================================================

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


                // ---------------------------------------------
                // alpha = x
                // ---------------------------------------------

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


                // ---------------------------------------------
                // alpha = y
                // ---------------------------------------------

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


                // =================================================
                // Eq.(37)
                //
                // Velocity gradients here are from Eq.(41)
                // =================================================

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


            // =================================================
            // Moment arrays
            // =================================================

            double m1[Q9] =
                {0};


            double m2[Q9] =
                {0};


            double meq1[Q9] =
                {0};


            double meq2[Q9] =
                {0};


            double mF1[Q9] =
                {0};


            double mF2[Q9] =
                {0};


            double mG1[Q9] =
                {0};


            double mG2[Q9] =
                {0};


            double mPost1[Q9] =
                {0};


            double mPost2[Q9] =
                {0};


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


            // =================================================
            // Local viscosity relaxation
            // =================================================

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


            // =================================================
            // Collision
            // =================================================

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
                        0.5 * s
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
                        0.5 * s
                    )

                    *
                    (
                        mF2[k]
                        +
                        mG2[k]
                    );
            }


            // =================================================
            // Inverse transform
            // =================================================

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
//
// top/bottom:
//
// half-way bounce-back
//
// no-flux phase boundary
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
                    x + cx5[q];


                if(xn < 0)
                {
                    xn +=
                        NX;
                }


                if(xn >= NX)
                {
                    xn -=
                        NX;
                }


                int yn =
                    y + cy5[q];


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
//
// bottom:
// no-slip transformed-g anti-bounce-back
//
// top:
// open boundary
//
// IMPLEMENTATION CHOICE:
// zero-gradient copy
// ============================================================

void streamFlow()
{
    // --------------------------------------------------------
    // Clear
    // --------------------------------------------------------

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


    // --------------------------------------------------------
    // Push streaming
    // --------------------------------------------------------

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
                    x + cx9[q];


                if(xn < 0)
                {
                    xn +=
                        NX;
                }


                if(xn >= NX)
                {
                    xn -=
                        NX;
                }


                int yn =
                    y + cy9[q];


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

                else if(yn < 0)
                {
                    // ========================================
                    // Bottom no-slip
                    // ========================================

                    int qo =
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

                else
                {
                    // Top outgoing populations leave domain.
                }
            }
        }
    }


    // ========================================================
    // Top open boundary
    // ========================================================

    int yTop =
        NY - 1;


    int yInner =
        NY - 2;


    for(int x = 0;
        x < NX;
        x++)
    {
        for(int q = 0;
            q < Q9;
            q++)
        {
            // Incoming from top

            if(cy9[q] < 0)
            {
                g1Next[yTop][x][q] =

                    g1Next[yInner][x][q];


                g2Next[yTop][x][q] =

                    g2Next[yInner][x][q];
            }
        }
    }


    // --------------------------------------------------------
    // next -> current
    // --------------------------------------------------------

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
//
// EVERY STEP
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
                double xx =

                    -0.5 * L

                    +

                    (x + 0.5)
                    *
                    dx;


                double yy =

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


    double maxSpeed =
        0.0;


    double maxChemical =
        0.0;


    double maxGradPhi =
        0.0;


    double maxLambdaACE =
        0.0;


    double maxForce =
        0.0;


    double phaseMass =
        0.0;


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


            double speed =

                sqrt
                (
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


            double gradPhiMagnitude =

                sqrt
                (
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


            double lambdaACE =

                sqrt
                (
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
                    fabs(lambdaACE)
                );


            double ff =

                sqrt
                (
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


    double relativeMassChange =

        (
            phaseMass
            -
            initialPhaseMass
        )

        /

        initialPhaseMass;


    double T =

        V0

        *
        (
            step
            *
            dt
        )

        /

        D;


    cout
        << "\n---------------------------------------------\n";


    cout
        << "step              = "
        << step
        << "\n";


    cout
        << "T = V0*t/D        = "
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
// Save phase snapshot
// ============================================================

void saveSnapshot(
    int step
)
{
    char filename[128];


    sprintf(
        filename,
        "droplet_phi_ACE_iso_step_%05d.dat",
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
