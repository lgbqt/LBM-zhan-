#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <fstream>

using namespace std;


// ============================================================
// Zhan 2023
// Layered Poiseuille flow
//
// First case:
// rhoA : rhoB = 100 : 1
//
// Phase field:
//     D2Q5 CHE
//
// Flow field:
//     D2Q9 MDF-LBM
//
// Grid:
//     10 x 100
//
// Physical domain:
//     x: periodic
//     y in [-0.5, 0.5]
//
// Boundary:
//     phase: half-way bounce-back
//     flow : half-way anti-bounce-back
//
// NOTE:
// Flow equilibrium / F / G are first constructed in velocity
// space, then transformed by the numerical M matrix.
// ============================================================


// ============================================================
// 1. Basic parameters
// ============================================================

const int NX = 10;
const int NY = 100;

const int Q5 = 5;
const int Q9 = 9;

const double h  = 0.5;

const double dx = 0.01;

const double c  = 100.0;

const double dt = 1.0e-4;

const double cs2 = c * c / 3.0;


// ------------------------------------------------------------
// Phase parameters
// ------------------------------------------------------------

const double phiA =  1.0;
const double phiB = -1.0;

const double W = 0.04;

const double sigma = 0.001;

const double Mphi = 0.1;

const double eta = 1.0;


// beta = 12 sigma / [W (phiA-phiB)^4]
const double betaPF =
    12.0 * sigma
    /
    (
        W
        *
        pow(phiA - phiB, 4.0)
    );


// kappa = 3/2 sigma W / (phiA-phiB)^2
const double kappa =
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


// ------------------------------------------------------------
// Fluid parameters
// ------------------------------------------------------------

const double rhoA = 100.0;
const double rhoB = 1.0;

const double nuA = 0.1;
const double nuB = 0.1;

const double Uc = 1.0e-4;


// ------------------------------------------------------------
// Relaxation rates
// ------------------------------------------------------------

const double sPhase1 =
    1.0
    /
    (
        Mphi
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


// nuA = nuB = 0.1 here,
// therefore one constant flow relaxation rate is enough.
const double sFlow1 =
    1.0
    /
    (
        nuA
        /
        (
            cs2
            *
            dt
        )
        +
        0.5
    );


// ============================================================
// 2. Lattice directions
// ============================================================


// ------------------------------------------------------------
// D2Q5
//
// 0 rest
// 1 east
// 2 north
// 3 west
// 4 south
// ------------------------------------------------------------

const int ex5[Q5] =
{
     0,
     1,
     0,
    -1,
     0
};

const int ey5[Q5] =
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


// ------------------------------------------------------------
// D2Q9
// ------------------------------------------------------------

const int ex9[Q9] =
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

const int ey9[Q9] =
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
// 3. Transformation matrices
// ============================================================


// ------------------------------------------------------------
// D2Q5 matrix
// ------------------------------------------------------------

const double M5[Q5][Q5] =
{
    { 1,  1,  1,  1,  1 },

    { 0,  1,  0, -1,  0 },

    { 0,  0,  1,  0, -1 },

    { 0,  1,  0,  1,  0 },

    { 0,  0,  1,  0,  1 }
};


const double invM5[Q5][Q5] =
{
    { 1,  0,    0,   -1,   -1 },

    { 0,  0.5,  0,    0.5,  0 },

    { 0,  0,    0.5,  0,    0.5 },

    { 0, -0.5,  0,    0.5,  0 },

    { 0,  0,   -0.5,  0,    0.5 }
};


// ------------------------------------------------------------
// D2Q9 transformation matrix from Zhan Appendix A
// ------------------------------------------------------------

const double M9[Q9][Q9] =
{
    { 1, 1, 1, 1, 1, 1, 1, 1, 1 },

    { 0, 1, 0,-1, 0, 1,-1,-1, 1 },

    { 0, 0, 1, 0,-1, 1, 1,-1,-1 },

    { 0, 1, 0, 1, 0, 1, 1, 1, 1 },

    { 0, 0, 0, 0, 0, 1,-1, 1,-1 },

    { 0, 0, 1, 0, 1, 1, 1, 1, 1 },

    { 0, 0, 0, 0, 0, 1, 1,-1,-1 },

    { 0, 0, 0, 0, 0, 1,-1,-1, 1 },

    { 0, 0, 0, 0, 0, 1, 1, 1, 1 }
};


const double invM9[Q9][Q9] =
{
    { 1, 0,    0,   -1,    0,   -1,    0,    0,    1 },

    { 0, 0.5,  0,    0.5,  0,    0,    0,   -0.5, -0.5 },

    { 0, 0,    0.5,  0,    0,    0.5, -0.5,  0,   -0.5 },

    { 0,-0.5,  0,    0.5,  0,    0,    0,    0.5, -0.5 },

    { 0, 0,   -0.5,  0,    0,    0.5,  0.5,  0,   -0.5 },

    { 0, 0,    0,    0,    0.25, 0,    0.25, 0.25, 0.25 },

    { 0, 0,    0,    0,   -0.25, 0,    0.25,-0.25, 0.25 },

    { 0, 0,    0,    0,    0.25, 0,   -0.25,-0.25, 0.25 },

    { 0, 0,    0,    0,   -0.25, 0,   -0.25, 0.25, 0.25 }
};


// ============================================================
// 4. Relaxation matrices
// ============================================================

const double sPhase[Q5] =
{
    1.0,
    sPhase1,
    sPhase1,
    1.0,
    1.0
};


const double sFlow[Q9] =
{
    1.0,

    sFlow1,
    sFlow1,

    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
    1.0
};


// ============================================================
// 5. Macroscopic arrays
// ============================================================

static double phi[NY][NX];

static double rho[NY][NX];

static double chem[NY][NX];

static double ux[NY][NX];
static double uy[NY][NX];

static double pressure[NY][NX];

static double forceX[NY][NX];
static double forceY[NY][NX];

static double gradRhoX[NY][NX];
static double gradRhoY[NY][NX];


// central-gradient version,
// used in G_i during evolution
static double dudx[NY][NX];
static double dudy[NY][NX];

static double dvdx[NY][NX];
static double dvdy[NY][NX];


// Eq.(40) local gradient,
// used for final Table-4 comparison
static double dudyLocal[NY][NX];


// ============================================================
// 6. Old-time arrays
// ============================================================

static double oldPhiUx[NY][NX];
static double oldPhiUy[NY][NX];

static double oldQxx[NY][NX];
static double oldQxy[NY][NX];

static double oldQyx[NY][NX];
static double oldQyy[NY][NX];


// ============================================================
// 7. Distribution functions
// ============================================================


// ------------------------------------------------------------
// Phase distribution
// ------------------------------------------------------------

static double gPhase[Q5][NY][NX];

static double gPhasePost[Q5][NY][NX];

static double gPhaseNew[Q5][NY][NX];


// ------------------------------------------------------------
// Momentum distributions
//
// gX -> rho*u equation
// gY -> rho*v equation
// ------------------------------------------------------------

static double gX[Q9][NY][NX];
static double gY[Q9][NY][NX];

static double gXPost[Q9][NY][NX];
static double gYPost[Q9][NY][NX];

static double gXNew[Q9][NY][NX];
static double gYNew[Q9][NY][NX];


// ============================================================
// 8. Analytical reference
// ============================================================

static double yCell[NY];

static double uExact[NY];

static double duExact[NY];


// Body force calculated from analytical profile
double GxBody = 0.0;


// ============================================================
// 9. Function declarations
// ============================================================

int wrapX(int x);

double phiAnalytical(double y);

double dynamicMuAnalytical(double y);

double simpsonIntegral(
    double a,
    double b,
    bool multiplyY,
    int nInterval
);

void buildAnalyticalReference();

double scalarNeighborNeumann(
    const double a[NY][NX],
    int x,
    int y,
    int q
);

double velocityNeighborDirichletZero(
    const double a[NY][NX],
    int x,
    int y,
    int q
);

void scalarGradientLaplacian(
    const double a[NY][NX],
    int x,
    int y,
    double &gx,
    double &gy,
    double &lap
);

void velocityGradient(
    const double a[NY][NX],
    int x,
    int y,
    double &gx,
    double &gy
);

void initializeMacroscopicFields();

void computeDensity();

void computeChemicalPotentialAndForce();

void initializeDistributions();

void recoverPhi();

void recoverVelocity();

void computeDensityGradient();

void computeVelocityGradients();

void computePressure();

void updateMacroscopicFields();

void initializeOldFields();

void phaseEquilibrium(
    int x,
    int y,
    double feq[Q5]
);

void flowEquilibrium(
    int alpha,
    int x,
    int y,
    double feq[Q9]
);

void collidePhase();

void collideFlow();

void saveCurrentAsOld();

void streamPhase();

void streamFlow();

void computeLocalDudyEq40();

bool checkFinite();

void printDiagnostics(int step);

void computeErrorsAndSave();


// ============================================================
// 10. Main
// ============================================================

int main(int argc, char *argv[])
{
    int maxSteps = 50000;

    if(argc >= 2)
    {
        maxSteps = atoi(argv[1]);
    }


    cout << scientific << setprecision(12);


    cout << "=============================================\n";
    cout << "Zhan 2023 layered Poiseuille flow\n";
    cout << "CHE + MDF-LBM\n";
    cout << "rhoA : rhoB = 100 : 1\n";
    cout << "=============================================\n";

    cout << "NX              = " << NX << "\n";
    cout << "NY              = " << NY << "\n";

    cout << "dx              = " << dx << "\n";
    cout << "dt              = " << dt << "\n";

    cout << "c               = " << c << "\n";
    cout << "cs2             = " << cs2 << "\n";

    cout << "beta            = " << betaPF << "\n";
    cout << "kappa           = " << kappa << "\n";

    cout << "sPhase1         = " << sPhase1 << "\n";
    cout << "sFlow1          = " << sFlow1 << "\n";

    cout << "maxSteps        = " << maxSteps << "\n";

    cout << "=============================================\n\n";


    // --------------------------------------------------------
    // Analytical reference
    // --------------------------------------------------------

    buildAnalyticalReference();


    // --------------------------------------------------------
    // Initial macroscopic fields
    // --------------------------------------------------------

    initializeMacroscopicFields();


    // --------------------------------------------------------
    // Initial distributions
    // --------------------------------------------------------

    initializeDistributions();


    // --------------------------------------------------------
    // Recover from actual g distributions
    // --------------------------------------------------------

    updateMacroscopicFields();


    // --------------------------------------------------------
    // Old-time quantities
    //
    // Therefore the first explicit time derivative is zero.
    // --------------------------------------------------------

    initializeOldFields();


    printDiagnostics(0);


    // --------------------------------------------------------
    // Convergence monitor
    // --------------------------------------------------------

    static double oldCheckU[NY][NX];

    for(int y = 0; y < NY; y++)
    {
        for(int x = 0; x < NX; x++)
        {
            oldCheckU[y][x] = ux[y][x];
        }
    }


    const int CHECK_INTERVAL = 1000;

    const double convergenceTolerance = 1.0e-8;


    // ========================================================
    // Time loop
    // ========================================================

    for(int step = 1; step <= maxSteps; step++)
    {
        // ----------------------------------------------------
        // Collision at t^n
        // ----------------------------------------------------

        collidePhase();

        collideFlow();


        // ----------------------------------------------------
        // Save old explicit-time quantities
        // ----------------------------------------------------

        saveCurrentAsOld();


        // ----------------------------------------------------
        // Streaming + boundaries
        // ----------------------------------------------------

        streamPhase();

        streamFlow();


        // ----------------------------------------------------
        // Recover t^(n+1)
        // ----------------------------------------------------

        updateMacroscopicFields();


        // ----------------------------------------------------
        // Stability check
        // ----------------------------------------------------

        if(step % CHECK_INTERVAL == 0)
        {
            if(!checkFinite())
            {
                cout << "\nERROR: NaN/Inf detected.\n";
                cout << "Stop at step " << step << "\n";

                return 1;
            }


            double maxU = 0.0;

            double maxDifference = 0.0;


            for(int y = 0; y < NY; y++)
            {
                for(int x = 0; x < NX; x++)
                {
                    double au = fabs(ux[y][x]);

                    if(au > maxU)
                    {
                        maxU = au;
                    }


                    double du =
                        fabs(
                            ux[y][x]
                            -
                            oldCheckU[y][x]
                        );

                    if(du > maxDifference)
                    {
                        maxDifference = du;
                    }
                }
            }


            double relativeChange = 0.0;

            if(maxU > 1.0e-30)
            {
                relativeChange =
                    maxDifference
                    /
                    maxU;
            }


            cout << "\n";
            cout << "Step = " << step << "\n";

            cout << "physical t       = "
                 << step * dt
                 << "\n";

            cout << "max |u|          = "
                 << maxU
                 << "\n";

            cout << "max delta u      = "
                 << maxDifference
                 << "\n";

            cout << "relative change  = "
                 << relativeChange
                 << "\n";


            printDiagnostics(step);


            for(int y = 0; y < NY; y++)
            {
                for(int x = 0; x < NX; x++)
                {
                    oldCheckU[y][x] =
                        ux[y][x];
                }
            }


            if(step >= 5000)
            {
                if(relativeChange < convergenceTolerance)
                {
                    cout << "\n";
                    cout << "Steady criterion reached.\n";
                    cout << "Stop automatically.\n";

                    maxSteps = step;

                    break;
                }
            }
        }
    }


    // --------------------------------------------------------
    // Final local gradient Eq.(40)
    // --------------------------------------------------------

    computeLocalDudyEq40();


    // --------------------------------------------------------
    // Error + output
    // --------------------------------------------------------

    computeErrorsAndSave();


    return 0;
}


// ============================================================
// 11. Simple helpers
// ============================================================

int wrapX(int x)
{
    if(x < 0)
    {
        return x + NX;
    }

    if(x >= NX)
    {
        return x - NX;
    }

    return x;
}


// ============================================================
// 12. Analytical profile
// ============================================================

double phiAnalytical(double y)
{
    return

        0.5
        *
        (phiA + phiB)

        +

        0.5
        *
        (phiA - phiB)
        *
        tanh(
            -2.0
            *
            y
            /
            W
        );
}


double dynamicMuAnalytical(double y)
{
    double phiValue =
        phiAnalytical(y);


    double muA =
        rhoA
        *
        nuA;


    double muB =
        rhoB
        *
        nuB;


    double muValue =

        (muA - muB)
        /
        (phiA - phiB)

        *

        (phiValue - phiB)

        +

        muB;


    return muValue;
}


// ============================================================
// 13. Simpson integral
//
// If multiplyY == false:
//       integral 1/mu dy
//
// If multiplyY == true:
//       integral y/mu dy
// ============================================================

double simpsonIntegral(
    double a,
    double b,
    bool multiplyY,
    int nInterval
)
{
    if(a == b)
    {
        return 0.0;
    }


    if(nInterval % 2 != 0)
    {
        nInterval++;
    }


    double ds =
        (b - a)
        /
        nInterval;


    double sum = 0.0;


    for(int n = 0; n <= nInterval; n++)
    {
        double y =
            a
            +
            n
            *
            ds;


        double value =
            1.0
            /
            dynamicMuAnalytical(y);


        if(multiplyY)
        {
            value *= y;
        }


        double coefficient = 2.0;


        if(n == 0 || n == nInterval)
        {
            coefficient = 1.0;
        }
        else if(n % 2 == 1)
        {
            coefficient = 4.0;
        }


        sum +=
            coefficient
            *
            value;
    }


    return

        ds
        /
        3.0
        *
        sum;
}


// ============================================================
// 14. Analytical Eq.(50)-(51)
// ============================================================

void buildAnalyticalReference()
{
    const int NINT = 100000;


    double A =
        simpsonIntegral(
            -h,
             h,
             true,
             NINT
        );


    double B =
        simpsonIntegral(
            -h,
             h,
             false,
             NINT
        );


    double C =
        simpsonIntegral(
            -h,
             0.0,
             true,
             NINT
        );


    double D =
        simpsonIntegral(
            -h,
             0.0,
             false,
             NINT
        );


    GxBody =
        Uc
        /
        (
            A
            *
            D
            /
            B

            -

            C
        );


    cout << "Analytical reference\n";

    cout << "A               = "
         << A
         << "\n";

    cout << "B               = "
         << B
         << "\n";

    cout << "C               = "
         << C
         << "\n";

    cout << "D               = "
         << D
         << "\n";

    cout << "Gx              = "
         << GxBody
         << "\n\n";


    for(int j = 0; j < NY; j++)
    {
        double y =
            -h
            +
            (j + 0.5)
            *
            dx;


        yCell[j] = y;


        double I1 =
            simpsonIntegral(
                -h,
                 y,
                 true,
                 20000
            );


        double I2 =
            simpsonIntegral(
                -h,
                 y,
                 false,
                 20000
            );


        uExact[j] =

            -GxBody
            *
            I1

            +

            GxBody
            *
            A
            /
            B
            *
            I2;


        double mu =
            dynamicMuAnalytical(y);


        duExact[j] =

            -GxBody
            *
            y
            /
            mu

            +

            GxBody
            *
            A
            /
            B
            /
            mu;
    }


    double maxU = 0.0;
    double yMax = 0.0;


    for(int j = 0; j < NY; j++)
    {
        if(uExact[j] > maxU)
        {
            maxU = uExact[j];
            yMax = yCell[j];
        }
    }


    cout << "Analytical max u = "
         << maxU
         << "\n";

    cout << "y at max u       = "
         << yMax
         << "\n";

    cout << "=============================================\n\n";
}


// ============================================================
// 15. Scalar ghost treatment: no flux
//
// At y outside:
//
// phi_ghost = phi_boundary_cell
//
// This gives zero normal derivative at the half-way wall.
// ============================================================

double scalarNeighborNeumann(
    const double a[NY][NX],
    int x,
    int y,
    int q
)
{
    int xp =
        wrapX(
            x + ex9[q]
        );


    int yp =
        y + ey9[q];


    if(yp < 0 || yp >= NY)
    {
        return a[y][xp];
    }


    return a[yp][xp];
}


// ============================================================
// 16. Velocity ghost treatment: no-slip u_wall = 0
//
// Wall is half-way between fluid node and ghost node.
//
// Linear reflection:
//
// u_ghost = -u_fluid
// ============================================================

double velocityNeighborDirichletZero(
    const double a[NY][NX],
    int x,
    int y,
    int q
)
{
    int xp =
        wrapX(
            x + ex9[q]
        );


    int yp =
        y + ey9[q];


    if(yp < 0 || yp >= NY)
    {
        return
            -a[y][xp];
    }


    return a[yp][xp];
}


// ============================================================
// 17. Isotropic scalar gradient + Laplacian
//
// Zhan Eq.(41)-(42)
// ============================================================

void scalarGradientLaplacian(
    const double a[NY][NX],
    int x,
    int y,
    double &gx,
    double &gy,
    double &lap
)
{
    gx = 0.0;
    gy = 0.0;
    lap = 0.0;


    double center =
        a[y][x];


    for(int q = 1; q < Q9; q++)
    {
        double neighbor =
            scalarNeighborNeumann(
                a,
                x,
                y,
                q
            );


        double cx =
            c
            *
            ex9[q];


        double cy =
            c
            *
            ey9[q];


        gx +=
            w9[q]
            *
            cx
            *
            neighbor;


        gy +=
            w9[q]
            *
            cy
            *
            neighbor;


        lap +=
            2.0
            *
            w9[q]
            *
            (
                neighbor
                -
                center
            );
    }


    gx /=
        cs2
        *
        dt;


    gy /=
        cs2
        *
        dt;


    lap /=
        cs2
        *
        dt
        *
        dt;
}


// ============================================================
// 18. Velocity gradient used inside G_i
//
// This is Eq.(41)-type second-order isotropic gradient.
// Paper allows this when C_{alpha beta gamma} = 0.
//
// Final Table-4 du/dy will be computed separately by Eq.(40).
// ============================================================

void velocityGradient(
    const double a[NY][NX],
    int x,
    int y,
    double &gx,
    double &gy
)
{
    gx = 0.0;
    gy = 0.0;


    for(int q = 1; q < Q9; q++)
    {
        double neighbor =
            velocityNeighborDirichletZero(
                a,
                x,
                y,
                q
            );


        double cx =
            c
            *
            ex9[q];


        double cy =
            c
            *
            ey9[q];


        gx +=
            w9[q]
            *
            cx
            *
            neighbor;


        gy +=
            w9[q]
            *
            cy
            *
            neighbor;
    }


    gx /=
        cs2
        *
        dt;


    gy /=
        cs2
        *
        dt;
}


// ============================================================
// 19. Initial macroscopic fields
// ============================================================

void initializeMacroscopicFields()
{
    for(int y = 0; y < NY; y++)
    {
        double yy =
            -h
            +
            (y + 0.5)
            *
            dx;


        for(int x = 0; x < NX; x++)
        {
            phi[y][x] =
                phiAnalytical(yy);


            ux[y][x] = 0.0;
            uy[y][x] = 0.0;

            pressure[y][x] = 0.0;
        }
    }


    computeDensity();

    computeChemicalPotentialAndForce();

    computeDensityGradient();

    computeVelocityGradients();
}


// ============================================================
// 20. Density
// ============================================================

void computeDensity()
{
    for(int y = 0; y < NY; y++)
    {
        for(int x = 0; x < NX; x++)
        {
            rho[y][x] =

                (rhoA - rhoB)
                /
                (phiA - phiB)

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
// 21. Chemical potential and total force
//
// mu_phi =
// 4 beta (phi-phiA)(phi-phiB)
//          (phi-(phiA+phiB)/2)
// - kappa lap(phi)
//
// F = mu_phi grad(phi) + (Gx, 0)
// ============================================================

void computeChemicalPotentialAndForce()
{
    static double gradPhiX[NY][NX];
    static double gradPhiY[NY][NX];
    static double lapPhi[NY][NX];


    for(int y = 0; y < NY; y++)
    {
        for(int x = 0; x < NX; x++)
        {
            double gx, gy, lap;


            scalarGradientLaplacian(
                phi,
                x,
                y,
                gx,
                gy,
                lap
            );


            gradPhiX[y][x] = gx;
            gradPhiY[y][x] = gy;
            lapPhi[y][x] = lap;
        }
    }


    for(int y = 0; y < NY; y++)
    {
        for(int x = 0; x < NX; x++)
        {
            double ph =
                phi[y][x];


            chem[y][x] =

                4.0
                *
                betaPF
                *
                (ph - phiA)
                *
                (ph - phiB)
                *
                (
                    ph
                    -
                    0.5
                    *
                    (phiA + phiB)
                )

                -

                kappa
                *
                lapPhi[y][x];


            double surfaceFx =
                chem[y][x]
                *
                gradPhiX[y][x];


            double surfaceFy =
                chem[y][x]
                *
                gradPhiY[y][x];


            forceX[y][x] =
                surfaceFx
                +
                GxBody;


            forceY[y][x] =
                surfaceFy;
        }
    }
}


// ============================================================
// 22. Phase equilibrium Eq.(9a)
// ============================================================

void phaseEquilibrium(
    int x,
    int y,
    double feq[Q5]
)
{
    double ph =
        phi[y][x];


    double mu =
        chem[y][x];


    for(int q = 0; q < Q5; q++)
    {
        double cx =
            c
            *
            ex5[q];


        double cy =
            c
            *
            ey5[q];


        if(q == 0)
        {
            feq[q] =

                ph

                +

                (w5[q] - 1.0)
                *
                eta
                *
                mu;
        }
        else
        {
            double cu =

                cx
                *
                ph
                *
                ux[y][x]

                +

                cy
                *
                ph
                *
                uy[y][x];


            feq[q] =

                w5[q]
                *
                eta
                *
                mu

                +

                w5[q]
                *
                cu
                /
                cs2;
        }
    }
}


// ============================================================
// 23. Flow equilibrium Eq.(10), C = 0
//
// alpha = 0 -> rho*u equation
// alpha = 1 -> rho*v equation
// ============================================================

void flowEquilibrium(
    int alpha,
    int x,
    int y,
    double feq[Q9]
)
{
    double r =
        rho[y][x];


    double u =
        ux[y][x];


    double v =
        uy[y][x];


    double p =
        pressure[y][x];


    double uAlpha =
        (alpha == 0)
        ?
        u
        :
        v;


    for(int q = 0; q < Q9; q++)
    {
        double cx =
            c
            *
            ex9[q];


        double cy =
            c
            *
            ey9[q];


        double fluxX;
        double fluxY;


        if(alpha == 0)
        {
            // rho*u equation

            fluxX =
                r
                *
                u
                *
                u
                +
                p;


            fluxY =
                r
                *
                u
                *
                v;
        }
        else
        {
            // rho*v equation

            fluxX =
                r
                *
                v
                *
                u;


            fluxY =
                r
                *
                v
                *
                v
                +
                p;
        }


        feq[q] =

            w9[q]
            *
            (
                r
                *
                uAlpha

                +

                (
                    cx
                    *
                    fluxX

                    +

                    cy
                    *
                    fluxY
                )
                /
                cs2
            );
    }
}


// ============================================================
// 24. Initialize distributions
// ============================================================

void initializeDistributions()
{
    for(int y = 0; y < NY; y++)
    {
        for(int x = 0; x < NX; x++)
        {
            // ------------------------------------------------
            // Phase
            // ------------------------------------------------

            double feqP[Q5];


            phaseEquilibrium(
                x,
                y,
                feqP
            );


            for(int q = 0; q < Q5; q++)
            {
                // S = 0 -> no F_i0 transform correction
                gPhase[q][y][x] =
                    feqP[q];
            }


            // ------------------------------------------------
            // Flow
            // ------------------------------------------------

            double feqX[Q9];
            double feqY[Q9];


            flowEquilibrium(
                0,
                x,
                y,
                feqX
            );


            flowEquilibrium(
                1,
                x,
                y,
                feqY
            );


            for(int q = 0; q < Q9; q++)
            {
                double FiX =
                    w9[q]
                    *
                    forceX[y][x];


                double FiY =
                    w9[q]
                    *
                    forceY[y][x];


                // g = f - dt F/2

                gX[q][y][x] =

                    feqX[q]

                    -

                    0.5
                    *
                    dt
                    *
                    FiX;


                gY[q][y][x] =

                    feqY[q]

                    -

                    0.5
                    *
                    dt
                    *
                    FiY;
            }
        }
    }
}


// ============================================================
// 25. Recover phi
// ============================================================

void recoverPhi()
{
    for(int y = 0; y < NY; y++)
    {
        for(int x = 0; x < NX; x++)
        {
            double sum = 0.0;


            for(int q = 0; q < Q5; q++)
            {
                sum +=
                    gPhase[q][y][x];
            }


            // S = 0
            phi[y][x] = sum;
        }
    }
}


// ============================================================
// 26. Recover velocity Eq.(35b)
// ============================================================

void recoverVelocity()
{
    for(int y = 0; y < NY; y++)
    {
        for(int x = 0; x < NX; x++)
        {
            double sumX = 0.0;
            double sumY = 0.0;


            for(int q = 0; q < Q9; q++)
            {
                sumX +=
                    gX[q][y][x];


                sumY +=
                    gY[q][y][x];
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
// 27. Density gradient
// ============================================================

void computeDensityGradient()
{
    for(int y = 0; y < NY; y++)
    {
        for(int x = 0; x < NX; x++)
        {
            double gx, gy, lap;


            scalarGradientLaplacian(
                rho,
                x,
                y,
                gx,
                gy,
                lap
            );


            gradRhoX[y][x] = gx;
            gradRhoY[y][x] = gy;
        }
    }
}


// ============================================================
// 28. Velocity gradients used in G_i
// ============================================================

void computeVelocityGradients()
{
    for(int y = 0; y < NY; y++)
    {
        for(int x = 0; x < NX; x++)
        {
            velocityGradient(
                ux,
                x,
                y,
                dudx[y][x],
                dudy[y][x]
            );


            velocityGradient(
                uy,
                x,
                y,
                dvdx[y][x],
                dvdy[y][x]
            );
        }
    }
}


// ============================================================
// 29. Pressure Eq.(32)
//
// d = 2
// C = 0
//
// Last time derivative term omitted,
// following paper low-Mach implementation.
// ============================================================

void computePressure()
{
    for(int y = 0; y < NY; y++)
    {
        for(int x = 0; x < NX; x++)
        {
            double firstMoments = 0.0;


            for(int q = 0; q < Q9; q++)
            {
                double cx =
                    c
                    *
                    ex9[q];


                double cy =
                    c
                    *
                    ey9[q];


                firstMoments +=

                    cx
                    *
                    gX[q][y][x]

                    +

                    cy
                    *
                    gY[q][y][x];
            }


            double u =
                ux[y][x];


            double v =
                uy[y][x];


            double densityTerm =

                0.5
                *
                dt
                *
                cs2
                *
                (
                    u
                    *
                    gradRhoX[y][x]

                    +

                    v
                    *
                    gradRhoY[y][x]
                );


            pressure[y][x] =

                0.5
                *
                (
                    firstMoments

                    -

                    rho[y][x]
                    *
                    (
                        u*u
                        +
                        v*v
                    )

                    +

                    densityTerm
                );
        }
    }
}


// ============================================================
// 30. Full macroscopic recovery
// ============================================================

void updateMacroscopicFields()
{
    recoverPhi();

    computeDensity();

    computeChemicalPotentialAndForce();

    recoverVelocity();

    computeDensityGradient();

    computeVelocityGradients();

    computePressure();
}


// ============================================================
// 31. Initialize old explicit-time quantities
// ============================================================

void initializeOldFields()
{
    for(int y = 0; y < NY; y++)
    {
        for(int x = 0; x < NX; x++)
        {
            double r =
                rho[y][x];


            double u =
                ux[y][x];


            double v =
                uy[y][x];


            double p =
                pressure[y][x];


            oldPhiUx[y][x] =
                phi[y][x]
                *
                u;


            oldPhiUy[y][x] =
                phi[y][x]
                *
                v;


            oldQxx[y][x] =
                r
                *
                u
                *
                u
                +
                p;


            oldQxy[y][x] =
                r
                *
                u
                *
                v;


            oldQyx[y][x] =
                r
                *
                v
                *
                u;


            oldQyy[y][x] =
                r
                *
                v
                *
                v
                +
                p;
        }
    }
}


// ============================================================
// 32. Phase collision
// ============================================================

void collidePhase()
{
    for(int y = 0; y < NY; y++)
    {
        for(int x = 0; x < NX; x++)
        {
            double currentPhiUx =
                phi[y][x]
                *
                ux[y][x];


            double currentPhiUy =
                phi[y][x]
                *
                uy[y][x];


            double dPhiUx =

                (
                    currentPhiUx

                    -

                    oldPhiUx[y][x]
                )

                /
                dt;


            double dPhiUy =

                (
                    currentPhiUy

                    -

                    oldPhiUy[y][x]
                )

                /
                dt;


            double feq[Q5];
            double Gi[Q5];


            phaseEquilibrium(
                x,
                y,
                feq
            );


            for(int q = 0; q < Q5; q++)
            {
                double cx =
                    c
                    *
                    ex5[q];


                double cy =
                    c
                    *
                    ey5[q];


                Gi[q] =

                    w5[q]

                    *

                    (
                        cx
                        *
                        dPhiUx

                        +

                        cy
                        *
                        dPhiUy
                    )

                    /
                    cs2;
            }


            double m[Q5] = {0};
            double meq[Q5] = {0};
            double mG[Q5] = {0};

            double mPost[Q5] = {0};


            for(int a = 0; a < Q5; a++)
            {
                for(int q = 0; q < Q5; q++)
                {
                    m[a] +=
                        M5[a][q]
                        *
                        gPhase[q][y][x];


                    meq[a] +=
                        M5[a][q]
                        *
                        feq[q];


                    mG[a] +=
                        M5[a][q]
                        *
                        Gi[q];
                }
            }


            for(int a = 0; a < Q5; a++)
            {
                double s =
                    sPhase[a];


                mPost[a] =

                    (1.0 - s)
                    *
                    m[a]

                    +

                    s
                    *
                    meq[a]

                    +

                    dt
                    *
                    (1.0 - 0.5*s)
                    *
                    mG[a];
            }


            for(int q = 0; q < Q5; q++)
            {
                double value = 0.0;


                for(int a = 0; a < Q5; a++)
                {
                    value +=
                        invM5[q][a]
                        *
                        mPost[a];
                }


                gPhasePost[q][y][x] =
                    value;
            }
        }
    }
}


// ============================================================
// 33. Flow collision
//
// Eq.(37), C = 0
// ============================================================

void collideFlow()
{
    for(int y = 0; y < NY; y++)
    {
        for(int x = 0; x < NX; x++)
        {
            double r =
                rho[y][x];


            double u =
                ux[y][x];


            double v =
                uy[y][x];


            double p =
                pressure[y][x];


            // ------------------------------------------------
            // Current flux tensors
            // ------------------------------------------------

            double Qxx =
                r*u*u
                +
                p;


            double Qxy =
                r*u*v;


            double Qyx =
                r*v*u;


            double Qyy =
                r*v*v
                +
                p;


            // ------------------------------------------------
            // Explicit time derivatives
            // ------------------------------------------------

            double dQxx =
                (Qxx - oldQxx[y][x])
                /
                dt;


            double dQxy =
                (Qxy - oldQxy[y][x])
                /
                dt;


            double dQyx =
                (Qyx - oldQyx[y][x])
                /
                dt;


            double dQyy =
                (Qyy - oldQyy[y][x])
                /
                dt;


            // ------------------------------------------------
            // Velocity-space equilibrium / F / G
            // ------------------------------------------------

            double feqX[Q9];
            double feqY[Q9];

            double FiX[Q9];
            double FiY[Q9];

            double GiX[Q9];
            double GiY[Q9];


            flowEquilibrium(
                0,
                x,
                y,
                feqX
            );


            flowEquilibrium(
                1,
                x,
                y,
                feqY
            );


            for(int q = 0; q < Q9; q++)
            {
                double cx =
                    c
                    *
                    ex9[q];


                double cy =
                    c
                    *
                    ey9[q];


                FiX[q] =
                    w9[q]
                    *
                    forceX[y][x];


                FiY[q] =
                    w9[q]
                    *
                    forceY[y][x];


                // --------------------------------------------
                // Eq.(37), alpha = x
                //
                // -rho c_beta d_x u_beta
                //
                // = -rho(
                //      cx * du/dx
                //    + cy * dv/dx )
                // --------------------------------------------

                GiX[q] =

                    w9[q]

                    *

                    (

                        (
                            cx
                            *
                            dQxx

                            +

                            cy
                            *
                            dQxy
                        )
                        /
                        cs2

                        +

                        u
                        *
                        (
                            cx
                            *
                            gradRhoX[y][x]

                            +

                            cy
                            *
                            gradRhoY[y][x]
                        )

                        -

                        r
                        *
                        (
                            cx
                            *
                            dudx[y][x]

                            +

                            cy
                            *
                            dvdx[y][x]
                        )
                    );


                // --------------------------------------------
                // Eq.(37), alpha = y
                //
                // -rho c_beta d_y u_beta
                //
                // = -rho(
                //      cx * du/dy
                //    + cy * dv/dy )
                // --------------------------------------------

                GiY[q] =

                    w9[q]

                    *

                    (

                        (
                            cx
                            *
                            dQyx

                            +

                            cy
                            *
                            dQyy
                        )
                        /
                        cs2

                        +

                        v
                        *
                        (
                            cx
                            *
                            gradRhoX[y][x]

                            +

                            cy
                            *
                            gradRhoY[y][x]
                        )

                        -

                        r
                        *
                        (
                            cx
                            *
                            dudy[y][x]

                            +

                            cy
                            *
                            dvdy[y][x]
                        )
                    );
            }


            // =================================================
            // alpha = x
            // =================================================

            double mX[Q9] = {0};
            double meqX[Q9] = {0};

            double mFX[Q9] = {0};
            double mGX[Q9] = {0};

            double mPostX[Q9] = {0};


            for(int a = 0; a < Q9; a++)
            {
                for(int q = 0; q < Q9; q++)
                {
                    mX[a] +=
                        M9[a][q]
                        *
                        gX[q][y][x];


                    meqX[a] +=
                        M9[a][q]
                        *
                        feqX[q];


                    mFX[a] +=
                        M9[a][q]
                        *
                        FiX[q];


                    mGX[a] +=
                        M9[a][q]
                        *
                        GiX[q];
                }
            }


            for(int a = 0; a < Q9; a++)
            {
                double s =
                    sFlow[a];


                mPostX[a] =

                    (1.0 - s)
                    *
                    mX[a]

                    +

                    s
                    *
                    meqX[a]

                    +

                    dt
                    *
                    (1.0 - 0.5*s)
                    *
                    (
                        mFX[a]
                        +
                        mGX[a]
                    );
            }


            for(int q = 0; q < Q9; q++)
            {
                double value = 0.0;


                for(int a = 0; a < Q9; a++)
                {
                    value +=
                        invM9[q][a]
                        *
                        mPostX[a];
                }


                gXPost[q][y][x] =
                    value;
            }


            // =================================================
            // alpha = y
            // =================================================

            double mY[Q9] = {0};
            double meqY[Q9] = {0};

            double mFY[Q9] = {0};
            double mGY[Q9] = {0};

            double mPostY[Q9] = {0};


            for(int a = 0; a < Q9; a++)
            {
                for(int q = 0; q < Q9; q++)
                {
                    mY[a] +=
                        M9[a][q]
                        *
                        gY[q][y][x];


                    meqY[a] +=
                        M9[a][q]
                        *
                        feqY[q];


                    mFY[a] +=
                        M9[a][q]
                        *
                        FiY[q];


                    mGY[a] +=
                        M9[a][q]
                        *
                        GiY[q];
                }
            }


            for(int a = 0; a < Q9; a++)
            {
                double s =
                    sFlow[a];


                mPostY[a] =

                    (1.0 - s)
                    *
                    mY[a]

                    +

                    s
                    *
                    meqY[a]

                    +

                    dt
                    *
                    (1.0 - 0.5*s)
                    *
                    (
                        mFY[a]
                        +
                        mGY[a]
                    );
            }


            for(int q = 0; q < Q9; q++)
            {
                double value = 0.0;


                for(int a = 0; a < Q9; a++)
                {
                    value +=
                        invM9[q][a]
                        *
                        mPostY[a];
                }


                gYPost[q][y][x] =
                    value;
            }
        }
    }
}


// ============================================================
// 34. Save current state for next explicit derivative
// ============================================================

void saveCurrentAsOld()
{
    for(int y = 0; y < NY; y++)
    {
        for(int x = 0; x < NX; x++)
        {
            double r =
                rho[y][x];


            double u =
                ux[y][x];


            double v =
                uy[y][x];


            double p =
                pressure[y][x];


            oldPhiUx[y][x] =
                phi[y][x]
                *
                u;


            oldPhiUy[y][x] =
                phi[y][x]
                *
                v;


            oldQxx[y][x] =
                r*u*u
                +
                p;


            oldQxy[y][x] =
                r*u*v;


            oldQyx[y][x] =
                r*v*u;


            oldQyy[y][x] =
                r*v*v
                +
                p;
        }
    }
}


// ============================================================
// 35. Phase streaming
//
// x: periodic
// y: half-way bounce-back
//
// No-flux phase boundary.
// ============================================================

void streamPhase()
{
    memset(
        gPhaseNew,
        0,
        sizeof(gPhaseNew)
    );


    for(int q = 0; q < Q5; q++)
    {
        for(int y = 0; y < NY; y++)
        {
            for(int x = 0; x < NX; x++)
            {
                int xn =
                    wrapX(
                        x
                        +
                        ex5[q]
                    );


                int yn =
                    y
                    +
                    ey5[q];


                if(yn >= 0 && yn < NY)
                {
                    gPhaseNew[q][yn][xn] =
                        gPhasePost[q][y][x];
                }
                else
                {
                    // half-way bounce-back

                    int qb =
                        opp5[q];


                    gPhaseNew[qb][y][x] =
                        gPhasePost[q][y][x];
                }
            }
        }
    }


    memcpy(
        gPhase,
        gPhaseNew,
        sizeof(gPhase)
    );
}


// ============================================================
// 36. Flow streaming
//
// x: periodic
//
// y:
// half-way anti-bounce-back
//
// Wall velocity = 0.
//
// ------------------------------------------------------------
// IMPORTANT
//
// The paper explicitly specifies half-way anti-bounce-back,
// but does not print the transformed-g boundary formula.
//
// Here we use the transformed-variable correction:
//
// g_in
// = -g_out^*
//   + 2 w_i (rho u_wall)
//   - dt w_i F_alpha
//
// Since u_wall = 0:
//
// g_in
// = -g_out^*
//   - dt w_i F_alpha
//
// This is an implementation choice derived from
//
// g = f - dt F/2.
//
// If Table-4 errors do not agree later,
// this boundary formula is one of the first items to recheck.
// ============================================================

void streamFlow()
{
    memset(
        gXNew,
        0,
        sizeof(gXNew)
    );


    memset(
        gYNew,
        0,
        sizeof(gYNew)
    );


    for(int q = 0; q < Q9; q++)
    {
        for(int y = 0; y < NY; y++)
        {
            for(int x = 0; x < NX; x++)
            {
                int xn =
                    wrapX(
                        x
                        +
                        ex9[q]
                    );


                int yn =
                    y
                    +
                    ey9[q];


                if(yn >= 0 && yn < NY)
                {
                    gXNew[q][yn][xn] =
                        gXPost[q][y][x];


                    gYNew[q][yn][xn] =
                        gYPost[q][y][x];
                }
                else
                {
                    int qb =
                        opp9[q];


                    // wall momentum rho*u_wall = 0
                    // wall momentum rho*v_wall = 0

                    gXNew[qb][y][x] =

                        -gXPost[q][y][x];


                    gYNew[qb][y][x] =

                        -gYPost[q][y][x];
                }
            }
        }
    }


    memcpy(
        gX,
        gXNew,
        sizeof(gX)
    );


    memcpy(
        gY,
        gYNew,
        sizeof(gY)
    );
}


// ============================================================
// 37. Eq.(40) local du/dy
//
// C_{alpha beta gamma} = 0
//
// Unknown vector:
//
// [ dv/dx ]
// [ du/dy ]
//
// Matrix:
//
// [ 1/s - 1/2     1/s       ]
// [ 1/s           1/s - 1/2 ]
// ============================================================

void computeLocalDudyEq40()
{
    double a =
        1.0
        /
        sFlow1

        -

        0.5;


    double b =
        1.0
        /
        sFlow1;


    double determinant =
        a*a
        -
        b*b;


    for(int y = 0; y < NY; y++)
    {
        for(int x = 0; x < NX; x++)
        {
            double sumCyGX = 0.0;

            double sumCxGY = 0.0;


            for(int q = 0; q < Q9; q++)
            {
                double cx =
                    c
                    *
                    ex9[q];


                double cy =
                    c
                    *
                    ey9[q];


                sumCyGX +=
                    cy
                    *
                    gX[q][y][x];


                sumCxGY +=
                    cx
                    *
                    gY[q][y][x];
            }


            double r =
                rho[y][x];


            double u =
                ux[y][x];


            double v =
                uy[y][x];


            double R1 =

                sumCyGX

                -

                r
                *
                u
                *
                v

                +

                0.5
                *
                dt
                *
                cs2
                *
                u
                *
                gradRhoY[y][x];


            double R2 =

                sumCxGY

                -

                r
                *
                u
                *
                v

                +

                0.5
                *
                dt
                *
                cs2
                *
                v
                *
                gradRhoX[y][x];


            double rhs1 =

                -R1
                /
                (
                    r
                    *
                    cs2
                    *
                    dt
                );


            double rhs2 =

                -R2
                /
                (
                    r
                    *
                    cs2
                    *
                    dt
                );


            // second component:
            // du/dy

            dudyLocal[y][x] =

                (
                    -b
                    *
                    rhs1

                    +

                    a
                    *
                    rhs2
                )

                /
                determinant;
        }
    }
}


// ============================================================
// 38. Stability check
// ============================================================

bool checkFinite()
{
    for(int y = 0; y < NY; y++)
    {
        for(int x = 0; x < NX; x++)
        {
            if(
                !isfinite(phi[y][x])
                ||
                !isfinite(rho[y][x])
                ||
                !isfinite(ux[y][x])
                ||
                !isfinite(uy[y][x])
                ||
                !isfinite(pressure[y][x])
            )
            {
                return false;
            }
        }
    }


    return true;
}


// ============================================================
// 39. Diagnostics
// ============================================================

void printDiagnostics(int step)
{
    double minPhi = 1.0e100;
    double maxPhi = -1.0e100;

    double minRho = 1.0e100;
    double maxRho = -1.0e100;

    double maxU = 0.0;
    double maxV = 0.0;

    double maxChem = 0.0;

    double maxSurfaceForce = 0.0;


    for(int y = 0; y < NY; y++)
    {
        for(int x = 0; x < NX; x++)
        {
            if(phi[y][x] < minPhi)
            {
                minPhi = phi[y][x];
            }

            if(phi[y][x] > maxPhi)
            {
                maxPhi = phi[y][x];
            }


            if(rho[y][x] < minRho)
            {
                minRho = rho[y][x];
            }

            if(rho[y][x] > maxRho)
            {
                maxRho = rho[y][x];
            }


            if(fabs(ux[y][x]) > maxU)
            {
                maxU =
                    fabs(ux[y][x]);
            }


            if(fabs(uy[y][x]) > maxV)
            {
                maxV =
                    fabs(uy[y][x]);
            }


            if(fabs(chem[y][x]) > maxChem)
            {
                maxChem =
                    fabs(chem[y][x]);
            }


            double surfaceFx =
                forceX[y][x]
                -
                GxBody;


            double surfaceFy =
                forceY[y][x];


            double surfaceMagnitude =
                sqrt(
                    surfaceFx*surfaceFx
                    +
                    surfaceFy*surfaceFy
                );


            if(surfaceMagnitude > maxSurfaceForce)
            {
                maxSurfaceForce =
                    surfaceMagnitude;
            }
        }
    }


    cout << "Diagnostics step " << step << "\n";

    cout << "phi range        = ["
         << minPhi
         << ", "
         << maxPhi
         << "]\n";

    cout << "rho range        = ["
         << minRho
         << ", "
         << maxRho
         << "]\n";

    cout << "max |chem|       = "
         << maxChem
         << "\n";

    cout << "max |Fs|         = "
         << maxSurfaceForce
         << "\n";

    cout << "max |u|          = "
         << maxU
         << "\n";

    cout << "max |v|          = "
         << maxV
         << "\n";
}


// ============================================================
// 40. Final errors + save profile
//
// Eq.(52):
//
// Err_chi =
// sum_j |chi_num - chi_exact|
// --------------------------------
// sum_j |chi_exact|
//
// Numerical profile is averaged over x.
// ============================================================

void computeErrorsAndSave()
{
    double numeratorU = 0.0;
    double denominatorU = 0.0;

    double numeratorDu = 0.0;
    double denominatorDu = 0.0;


    ofstream fout(
        "poiseuille_numerical_100to1_B.dat"
    );


    fout << scientific
         << setprecision(15);


    for(int y = 0; y < NY; y++)
    {
        double uAverage = 0.0;

        double duAverage = 0.0;

        double phiAverage = 0.0;

        double rhoAverage = 0.0;


        for(int x = 0; x < NX; x++)
        {
            uAverage +=
                ux[y][x];


            duAverage +=
                dudyLocal[y][x];


            phiAverage +=
                phi[y][x];


            rhoAverage +=
                rho[y][x];
        }


        uAverage /= NX;

        duAverage /= NX;

        phiAverage /= NX;

        rhoAverage /= NX;


        numeratorU +=
            fabs(
                uAverage
                -
                uExact[y]
            );


        denominatorU +=
            fabs(
                uExact[y]
            );


        numeratorDu +=
            fabs(
                duAverage
                -
                duExact[y]
            );


        denominatorDu +=
            fabs(
                duExact[y]
            );


        // columns:
        //
        // y
        // phi
        // rho
        // u_num
        // u_exact
        // du_num_local
        // du_exact

        fout
            << yCell[y] << " "
            << phiAverage << " "
            << rhoAverage << " "
            << uAverage << " "
            << uExact[y] << " "
            << duAverage << " "
            << duExact[y] << "\n";
    }


    fout.close();


    double errU =
        numeratorU
        /
        denominatorU;


    double errDu =
        numeratorDu
        /
        denominatorDu;


    cout << "\n";
    cout << "=============================================\n";
    cout << "Final Poiseuille errors\n";
    cout << "=============================================\n";

    cout << "Err_u numerical        = "
         << errU
         << "\n";

    cout << "Err_u paper CHE        = "
         << 3.8011e-2
         << "\n";


    cout << "\n";


    cout << "Err_du/dy numerical    = "
         << errDu
         << "\n";

    cout << "Err_du/dy paper CHE    = "
         << 2.9178e-2
         << "\n";


    cout << "=============================================\n";

    cout << "Saved:\n";
    cout << "poiseuille_numerical_100to1_B.dat\n";
}