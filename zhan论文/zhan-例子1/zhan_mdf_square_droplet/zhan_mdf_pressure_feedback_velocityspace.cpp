#include <iostream>
#include <cmath>
#include <iomanip>
#include <algorithm>

using namespace std;


// ============================================================
// 1. Parameters
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
// 2. Free-energy parameters
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


const double squareLength =
    sqrt(
        PI * Rfinal * Rfinal
        +
        (4.0 - PI) * W * W / 4.0
    );


// ============================================================
// 3. Relaxation
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
// 4. D2Q5
// ============================================================

const int cx5[Q5] =
{
     0, 1, 0,-1, 0
};

const int cy5[Q5] =
{
     0, 0, 1, 0,-1
};

const double w5[Q5] =
{
    1.0/3.0,
    1.0/6.0,
    1.0/6.0,
    1.0/6.0,
    1.0/6.0
};


// ============================================================
// 5. D2Q9
// ============================================================

const int cx9[Q9] =
{
     0, 1, 0,-1, 0, 1,-1,-1, 1
};

const int cy9[Q9] =
{
     0, 0, 1, 0,-1, 1, 1,-1,-1
};

const double w9[Q9] =
{
    4.0/9.0,

    1.0/9.0,
    1.0/9.0,
    1.0/9.0,
    1.0/9.0,

    1.0/36.0,
    1.0/36.0,
    1.0/36.0,
    1.0/36.0
};


// ============================================================
// 6. Zhan Appendix A M9
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


const int nSteps = 10;


// ============================================================
// 7. Global arrays
// ============================================================

double phi[NY][NX];
double rho[NY][NX];

double gradPhiX[NY][NX];
double gradPhiY[NY][NX];

double chemicalPotential[NY][NX];

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

double g1[NY][NX][Q9];
double g2[NY][NX][Q9];

double g1Post[NY][NX][Q9];
double g2Post[NY][NX][Q9];

double g1Next[NY][NX][Q9];
double g2Next[NY][NX][Q9];

double initialMomentumX = 0.0;
double initialMomentumY = 0.0;


// ============================================================
// 8. Function declarations
// ============================================================

void checkM9();

void initializePhi();

void computeDensity();

void computeGradPhi();

void computeChemicalPotential();

void computeSurfaceForce();

void computeGradRho();

void initializeFlowFields();

void buildEquilibriumMoments(
    int y,
    int x,
    double meq1[Q9],
    double meq2[Q9]
);

void buildForceMoments(
    int y,
    int x,
    double mF1[Q9],
    double mF2[Q9]
);

void initializeFlowDistribution();

void recoverVelocity();

void computeVelocityGradients();

void recoverPressure();

void collideFlow();

void streamFlow();

bool checkNaN(int step);

void printDiagnostics(int step);


// ============================================================
// main
// ============================================================

int main()
{
    cout
        << scientific
        << setprecision(12);


    cout
        << "================================================\n";

    cout
        << "Zhan 2023 MDF diagnostic\n";

    cout
        << "Force             = ON\n";

    cout
        << "G                 = OFF\n";

    cout
        << "Pressure recover  = ON\n";

    cout
        << "Pressure feedback = ON\n";

    cout
        << "meq construction  = Eq.(10) -> M * feq\n";

    cout
        << "================================================\n\n";


    cout << "NX       = " << NX << "\n";
    cout << "NY       = " << NY << "\n";
    cout << "dx       = " << dx << "\n";
    cout << "dt       = " << dt << "\n";
    cout << "c        = " << c << "\n";
    cout << "cs2      = " << cs2 << "\n";
    cout << "nu       = " << nu << "\n";
    cout << "su1      = " << su1 << "\n";
    cout << "betaPF   = " << betaPF << "\n";
    cout << "kappaPF  = " << kappaPF << "\n\n";


    checkM9();


    initializePhi();

    computeDensity();

    computeGradPhi();

    computeChemicalPotential();

    computeSurfaceForce();

    computeGradRho();

    initializeFlowFields();

    initializeFlowDistribution();


    recoverVelocity();

    computeVelocityGradients();

    recoverPressure();


    initialMomentumX = 0.0;
    initialMomentumY = 0.0;


    for (int y = 0; y < NY; y++)
    {
        for (int x = 0; x < NX; x++)
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


    printDiagnostics(0);


    for (int step = 1;
         step <= nSteps;
         step++)
    {
        recoverVelocity();

        computeGradRho();

        computeVelocityGradients();

        recoverPressure();


        collideFlow();

        streamFlow();


        recoverVelocity();

        computeGradRho();

        computeVelocityGradients();

        recoverPressure();


        if (
            checkNaN(step)
        )
        {
            cout
                << "\nSimulation stopped.\n";

            return 1;
        }


        printDiagnostics(step);
    }


    cout
        << "\n================================================\n";

    cout
        << "Diagnostic finished.\n";

    cout
        << "================================================\n";


    return 0;
}


// ============================================================
// Check M9 Minv9
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
                    fabs(
                        value
                        -
                        exact
                    )
                );
        }
    }


    cout
        << "Max |M9*Minv9-I| = "
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
// Density
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
// grad(phi)
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


                gx +=
                    w5[q]
                    *
                    (
                        c
                        *
                        cx5[q]
                    )
                    *
                    phi[yn][xn];


                gy +=
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
// grad(rho)
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
                    c
                    *
                    cx9[q];


                const double ciy =
                    c
                    *
                    cy9[q];


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
// Initial flow
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
// Eq.(10) -> feq -> meq = M feq
//
// This is the key diagnostic.
// ============================================================

void buildEquilibriumMoments(
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


    for (int q = 0; q < Q9; q++)
    {
        const double cix =
            c
            *
            cx9[q];

        const double ciy =
            c
            *
            cy9[q];


        // ----------------------------------------------------
        // Zhan Eq.(10), C = 0
        // ----------------------------------------------------

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


    for (int k = 0; k < Q9; k++)
    {
        meq1[k] = 0.0;
        meq2[k] = 0.0;


        for (int q = 0; q < Q9; q++)
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
// Eq.(11) -> Fi -> mF = M Fi
// ============================================================

void buildForceMoments(
    int y,
    int x,
    double mF1[Q9],
    double mF2[Q9]
)
{
    double Fi1[Q9];
    double Fi2[Q9];


    for (int q = 0; q < Q9; q++)
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


    for (int k = 0; k < Q9; k++)
    {
        mF1[k] = 0.0;
        mF2[k] = 0.0;


        for (int q = 0; q < Q9; q++)
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
// Initialize g
// ============================================================

void initializeFlowDistribution()
{
    for (int y = 0; y < NY; y++)
    {
        for (int x = 0; x < NX; x++)
        {
            double meq1[Q9];
            double meq2[Q9];

            double mF1[Q9];
            double mF2[Q9];

            double mg1[Q9];
            double mg2[Q9];


            buildEquilibriumMoments(
                y,
                x,
                meq1,
                meq2
            );


            buildForceMoments(
                y,
                x,
                mF1,
                mF2
            );


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
// Velocity recovery
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
// Velocity gradients
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


            dudx[y][x] =
                ux_x
                /
                (
                    cs2 * dt
                );


            dudy[y][x] =
                ux_y
                /
                (
                    cs2 * dt
                );


            dvdx[y][x] =
                uy_x
                /
                (
                    cs2 * dt
                );


            dvdy[y][x] =
                uy_y
                /
                (
                    cs2 * dt
                );
        }
    }
}


// ============================================================
// Pressure recovery
//
// Eq.(32), g replacing f.
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
// Collision
//
// G OFF
// ============================================================

void collideFlow()
{
    for (int y = 0; y < NY; y++)
    {
        for (int x = 0; x < NX; x++)
        {
            double m1[Q9];
            double m2[Q9];

            double meq1[Q9];
            double meq2[Q9];

            double mF1[Q9];
            double mF2[Q9];

            double mPost1[Q9];
            double mPost2[Q9];


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


            buildEquilibriumMoments(
                y,
                x,
                meq1,
                meq2
            );


            buildForceMoments(
                y,
                x,
                mF1,
                mF2
            );


            for (int k = 0; k < Q9; k++)
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
// Streaming
// ============================================================

void streamFlow()
{
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


    for (int y = 0; y < NY; y++)
    {
        for (int x = 0; x < NX; x++)
        {
            for (int q = 0; q < Q9; q++)
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
// NaN check
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
                    << "\nNaN/Inf at step "
                    << step
                    << ", x="
                    << x
                    << ", y="
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

void printDiagnostics(int step)
{
    double maxU = 0.0;
    double maxV = 0.0;
    double maxSpeed = 0.0;

    double maxP = 0.0;
    double maxF = 0.0;

    double maxMeq = 0.0;
    double maxGdist = 0.0;

    double momentumX = 0.0;
    double momentumY = 0.0;

    double div2 = 0.0;


    for (int y = 0; y < NY; y++)
    {
        for (int x = 0; x < NX; x++)
        {
            maxU =
                max(
                    maxU,
                    fabs(
                        ux[y][x]
                    )
                );


            maxV =
                max(
                    maxV,
                    fabs(
                        uy[y][x]
                    )
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


            maxP =
                max(
                    maxP,
                    fabs(
                        pressure[y][x]
                    )
                );


            const double Fmag =
                sqrt(
                    forceX[y][x] * forceX[y][x]
                    +
                    forceY[y][x] * forceY[y][x]
                );


            maxF =
                max(
                    maxF,
                    Fmag
                );


            double meq1[Q9];
            double meq2[Q9];


            buildEquilibriumMoments(
                y,
                x,
                meq1,
                meq2
            );


            for (int k = 0; k < Q9; k++)
            {
                maxMeq =
                    max(
                        maxMeq,
                        fabs(
                            meq1[k]
                        )
                    );


                maxMeq =
                    max(
                        maxMeq,
                        fabs(
                            meq2[k]
                        )
                    );
            }


            for (int q = 0; q < Q9; q++)
            {
                maxGdist =
                    max(
                        maxGdist,
                        fabs(
                            g1[y][x][q]
                        )
                    );


                maxGdist =
                    max(
                        maxGdist,
                        fabs(
                            g2[y][x][q]
                        )
                    );
            }


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


            div2 +=
                divU
                *
                divU;
        }
    }


    const double dPx =
        (
            momentumX
            -
            initialMomentumX
        )
        /
        (
            NX * NY
        );


    const double dPy =
        (
            momentumY
            -
            initialMomentumY
        )
        /
        (
            NX * NY
        );


    const double Drms =
        sqrt(
            div2
            /
            (
                NX * NY
            )
        );


    cout
        << "\n--------------------------------------------\n";

    cout
        << "Step = "
        << step
        << "\n";

    cout
        << "max|F|               = "
        << maxF
        << "\n";

    cout
        << "max|u|               = "
        << maxU
        << "\n";

    cout
        << "max|v|               = "
        << maxV
        << "\n";

    cout
        << "max speed            = "
        << maxSpeed
        << "\n";

    cout
        << "max|p|               = "
        << maxP
        << "\n";

    cout
        << "max|meq=M*feq|       = "
        << maxMeq
        << "\n";

    cout
        << "max|g|               = "
        << maxGdist
        << "\n";

    cout
        << "<Delta rho*u>        = "
        << dPx
        << "\n";

    cout
        << "<Delta rho*v>        = "
        << dPy
        << "\n";

    cout
        << "D_rms                = "
        << Drms
        << "\n";
}