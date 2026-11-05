#include <iostream>
#include <cmath>
#include <iomanip>
#include <algorithm>

using namespace std;


// ============================================================
// 1. Basic parameters
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

const double cs2 =
    c * c / 3.0;


// ============================================================
// 2. Physical parameters
// ============================================================

const double nu =
    0.1;

const double sigma =
    0.001;

const double rhoA =
    1000.0;

const double rhoB =
    1.0;

const double phiA =
    1.0;

const double phiB =
   -1.0;

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
// 4. Initial square length
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
// 5. Flow relaxation
// ============================================================

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
     0, 1, 0, -1, 0
};

const int cy5[Q5] =
{
     0, 0, 1, 0, -1
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
// 8. Zhan Appendix A D2Q9 matrix
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
// 9. Test steps
// ============================================================

const int nSteps =
    10;


// ============================================================
// 10. Global arrays
// ============================================================

double phi[NY][NX];
double rho[NY][NX];

double gradPhiX[NY][NX];
double gradPhiY[NY][NX];

double chemicalPotential[NY][NX];

double forceX[NY][NX];
double forceY[NY][NX];

double ux[NY][NX];
double uy[NY][NX];

double dudx[NY][NX];
double dudy[NY][NX];

double dvdx[NY][NX];
double dvdy[NY][NX];

double g1[NY][NX][Q9];
double g2[NY][NX][Q9];

double g1Post[NY][NX][Q9];
double g2Post[NY][NX][Q9];

double g1Next[NY][NX][Q9];
double g2Next[NY][NX][Q9];

double initialMomentumX =
    0.0;

double initialMomentumY =
    0.0;


// ============================================================
// Function declarations
// ============================================================

void checkM9();

void initializePhi();

void computeDensity();

void computeGradPhi();

void computeChemicalPotential();

void computeSurfaceForce();

void initializeFlowFields();

void initializeFlowDistribution();

void recoverVelocity();

void computeVelocityGradients();

void collideFlow();

void streamFlow();

bool checkNaN(
    int step
);

void printDiagnostics(
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
        << "============================================\n";

    cout
        << "Zhan 2023 MDF-LBM diagnostic\n";

    cout
        << "Force = ON\n";

    cout
        << "G     = OFF\n";

    cout
        << "p     = FIXED ZERO\n";

    cout
        << "============================================\n\n";


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
        << "nu                 = "
        << nu
        << "\n";

    cout
        << "su1                = "
        << su1
        << "\n";

    cout
        << "betaPF             = "
        << betaPF
        << "\n";

    cout
        << "kappaPF            = "
        << kappaPF
        << "\n\n";


    // ========================================================
    // Matrix check
    // ========================================================

    checkM9();


    // ========================================================
    // Frozen phase / density
    // ========================================================

    initializePhi();

    computeDensity();

    computeGradPhi();

    computeChemicalPotential();

    computeSurfaceForce();


    // ========================================================
    // Initial flow
    // ========================================================

    initializeFlowFields();


    // ========================================================
    // transformed g initialization
    //
    // g = feq - dt F/2
    //
    // p = 0
    // ========================================================

    initializeFlowDistribution();


    // ========================================================
    // Recover Step 0
    // ========================================================

    recoverVelocity();

    computeVelocityGradients();


    // ========================================================
    // Initial momentum
    // ========================================================

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


    // ========================================================
    // Time loop
    // ========================================================

    for (int step = 1;
         step <= nSteps;
         step++)
    {
        recoverVelocity();

        computeVelocityGradients();


        // ----------------------------------------------------
        // Collision:
        //
        // F ON
        // G OFF
        // p = 0
        // ----------------------------------------------------

        collideFlow();


        // ----------------------------------------------------
        // streaming
        // ----------------------------------------------------

        streamFlow();


        // ----------------------------------------------------
        // recover new state
        // ----------------------------------------------------

        recoverVelocity();

        computeVelocityGradients();


        // ----------------------------------------------------
        // check
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


        printDiagnostics(
            step
        );
    }


    cout
        << "\n============================================\n";

    cout
        << "Diagnostic test finished.\n";

    cout
        << "============================================\n";


    return 0;
}


// ============================================================
// Matrix check
// ============================================================

void checkM9()
{
    double maxError =
        0.0;


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
        << "Max |M9*Minv9-I|  = "
        << maxError
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
// grad(phi)
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


                gx
                    +=
                w5[q]
                *
                (
                    c
                    *
                    cx5[q]
                )
                *
                phi[yn][xn];


                gy
                    +=
                w5[q]
                *
                (
                    c
                    *
                    cy5[q]
                )
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


                lapPhi
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
// Surface force
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
// Initial flow
//
// p is NOT stored.
// In this diagnostic:
//
// p = 0 always.
// ============================================================

void initializeFlowFields()
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
        }
    }
}


// ============================================================
// Initialize transformed distributions
//
// p = 0
//
// mg = meq - dt*mF/2
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
            const double r =
                rho[y][x];

            const double u =
                ux[y][x];

            const double v =
                uy[y][x];

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


            // ------------------------------------------------
            // p = 0
            // ------------------------------------------------

            const double Qxx =
                r * u * u;

            const double Qxy =
                r * u * v;

            const double Qyy =
                r * v * v;


            // x momentum

            meq1[0] =
                r * u;

            meq1[1] =
                Qxx;

            meq1[2] =
                Qxy;

            meq1[3] =
                r * u / 3.0;

            meq1[4] =
                0.0;

            meq1[5] =
                r * u / 3.0;

            meq1[6] =
                Qxy
                /
                (3.0 * c);

            meq1[7] =
                Qxx
                /
                (3.0 * c);

            meq1[8] =
                r * u / 9.0;


            // y momentum

            meq2[0] =
                r * v;

            meq2[1] =
                Qxy;

            meq2[2] =
                Qyy;

            meq2[3] =
                r * v / 3.0;

            meq2[4] =
                0.0;

            meq2[5] =
                r * v / 3.0;

            meq2[6] =
                Qyy
                /
                (3.0 * c);

            meq2[7] =
                Qxy
                /
                (3.0 * c);

            meq2[8] =
                r * v / 9.0;


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


            // velocity space

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
            double sumG1 =
                0.0;

            double sumG2 =
                0.0;


            for (int q = 0;
                 q < Q9;
                 q++)
            {
                sumG1
                    +=
                g1[y][x][q];

                sumG2
                    +=
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
// Velocity gradients
//
// only used for diagnostics
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
// Collision
//
// Force = ON
// G     = OFF
// p     = 0
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
            const double r =
                rho[y][x];

            const double u =
                ux[y][x];

            const double v =
                uy[y][x];

            const double Fx =
                forceX[y][x];

            const double Fy =
                forceY[y][x];


            // ------------------------------------------------
            // PRESSURE OFF
            // ------------------------------------------------

            const double Qxx =
                r
                *
                u
                *
                u;

            const double Qxy =
                r
                *
                u
                *
                v;

            const double Qyy =
                r
                *
                v
                *
                v;


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


            // =================================================
            // m = M g
            // =================================================

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


            // =================================================
            // equilibrium x momentum
            //
            // p = 0
            // =================================================

            meq1[0] =
                r * u;

            meq1[1] =
                Qxx;

            meq1[2] =
                Qxy;

            meq1[3] =
                r * u / 3.0;

            meq1[4] =
                0.0;

            meq1[5] =
                r * u / 3.0;

            meq1[6] =
                Qxy
                /
                (3.0 * c);

            meq1[7] =
                Qxx
                /
                (3.0 * c);

            meq1[8] =
                r * u / 9.0;


            // =================================================
            // equilibrium y momentum
            //
            // p = 0
            // =================================================

            meq2[0] =
                r * v;

            meq2[1] =
                Qxy;

            meq2[2] =
                Qyy;

            meq2[3] =
                r * v / 3.0;

            meq2[4] =
                0.0;

            meq2[5] =
                r * v / 3.0;

            meq2[6] =
                Qyy
                /
                (3.0 * c);

            meq2[7] =
                Qxy
                /
                (3.0 * c);

            meq2[8] =
                r * v / 9.0;


            // =================================================
            // force moments
            // =================================================

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


            // =================================================
            // G OFF
            // =================================================

            for (int k = 0;
                 k < Q9;
                 k++)
            {
                mG1[k] =
                    0.0;

                mG2[k] =
                    0.0;
            }


            // =================================================
            // MRT collision
            //
            // mPost =
            // (I-S)m
            // + S meq
            // + dt(I-S/2)mF
            //
            // because mG = 0
            // =================================================

            for (int k = 0;
                 k < Q9;
                 k++)
            {
                mPost1[k]
                    =
                (
                    1.0
                    -
                    Su[k]
                )
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
                mF1[k];


                mPost2[k]
                    =
                (
                    1.0
                    -
                    Su[k]
                )
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
                mF2[k];
            }


            // =================================================
            // back to velocity space
            // =================================================

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
// Streaming
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
// NaN check
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
                    ux[y][x]
                )
                ||
                !isfinite(
                    uy[y][x]
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
                    << "u = "
                    << ux[y][x]
                    << "\n";

                cout
                    << "v = "
                    << uy[y][x]
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
    double maxUx =
        0.0;

    double maxUy =
        0.0;

    double maxSpeed =
        0.0;

    double maxForce =
        0.0;

    double maxMu =
        0.0;

    double momentumX =
        0.0;

    double momentumY =
        0.0;

    double divergenceSquareSum =
        0.0;


    for (int y = 0;
         y < NY;
         y++)
    {
        for (int x = 0;
             x < NX;
             x++)
        {
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


            const double forceMagnitude =
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
                    forceMagnitude
                );


            maxMu =
                max(
                    maxMu,
                    fabs(
                        chemicalPotential[y][x]
                    )
                );


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


            divergenceSquareSum
                +=
            divU
            *
            divU;
        }
    }


    const double dMomentumX =
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


    const double dMomentumY =
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
        << "G                    = OFF\n";

    cout
        << "p                    = 0\n";

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
        << "<Delta rho*u>        = "
        << dMomentumX
        << "\n";

    cout
        << "<Delta rho*v>        = "
        << dMomentumY
        << "\n";

    cout
        << "D_rms                = "
        << Drms
        << "\n";
}