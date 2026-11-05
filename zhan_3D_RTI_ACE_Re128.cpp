#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include <string>

using namespace std;


// ============================================================
// Zhan et al. 2023
// 3D Rayleigh-Taylor instability
//
// ACE + MDF-LBM
//
// Phase : D3Q7
// Flow  : D3Q15 full MRT
//
// ------------------------------------------------------------
// VERIFIED STABLE VERSION
//
// Surface force         : ON
// dQ/dt term in Eq.(37) : OFF
//
// ------------------------------------------------------------
// NEW:
//
// checkpoint / restart support
//
// Fresh run:
//
//   .\rti3d.exe 1600
//
// Restart:
//
//   .\rti3d.exe 3200 RTI3D_checkpoint_step_01600.bin
//
// ============================================================


// ============================================================
// 1. GRID
// ============================================================

const int NX = 64;
const int NY = 64;
const int NZ = 256;

const int Q7  = 7;
const int Q15 = 15;

const double dx = 1.0;
const double dt = 1.0;

const double PI =
    3.14159265358979323846;


// ============================================================
// MODEL SWITCHES
// ============================================================

// Surface tension force
//
// Fs = chemical * grad(phi)

const bool USE_SURFACE_FORCE =
    true;


// IMPORTANT:
//
// We already diagnosed that explicitly inserting
//
// d/dt(rho*u_alpha*u_beta + p delta_ab)
//
// into Eq.(37) generates numerical positive feedback
// in the present 3D implementation.
//
// Therefore it remains OFF in this stable version.

const bool USE_DQDT_IN_G =
    false;


// ============================================================
// 2. LATTICE SOUND SPEED
// ============================================================

const double cs2Phase =
    1.0 / 4.0;

const double cs2Flow =
    1.0 / 3.0;


// ============================================================
// 3. PHYSICAL PARAMETERS
// ============================================================

const double d =
    64.0;

const double W =
    4.0;

const double phiA =
    1.0;

const double phiB =
    0.0;

const double Atwood =
    0.5;

const double Reynolds =
    128.0;

const double Peclet =
    1024.0;


// ============================================================
// Characteristic lattice velocity
//
// Ug = sqrt(g*d)
// ============================================================

const double Ug =
    0.04;

const double gravity =
    Ug
    *
    Ug
    /
    d;


// ============================================================
// Density
//
// rhoA = 3
// rhoB = 1
//
// At = (3-1)/(3+1) = 0.5
// ============================================================

const double rhoA =
    3.0;

const double rhoB =
    1.0;

const double rhoRef =
    0.5
    *
    (
        rhoA
        +
        rhoB
    );


// ============================================================
// Equal kinematic viscosity
//
// Re = d*Ug/nu
// ============================================================

const double nu =
    d
    *
    Ug
    /
    Reynolds;


const double muA =
    rhoA
    *
    nu;


const double muB =
    rhoB
    *
    nu;


// ============================================================
// Surface tension parameters
//
// Ca = 960 is our reproduction choice.
// ============================================================

const double Capillary =
    960.0;


const double sigma =
    muA
    *
    Ug
    /
    Capillary;


// ============================================================
// Free energy parameters
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
// ACE mobility
//
// Pe = d*Ug/M
// ============================================================

const double mobility =
    d
    *
    Ug
    /
    Peclet;


const double eta =
    1.0;


const double epsNormal =
    1.0e-14;


// ============================================================
// 4. RELAXATION PARAMETERS
// ============================================================

// Phase D3Q7

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


// Flow D3Q15

const double su1 =
    1.0
    /
    (
        nu
        /
        (
            cs2Flow
            *
            dt
        )
        +
        0.5
    );


const double su0 =
    1.0;

const double su2 =
    1.0;

const double su3 =
    1.0;

const double su4 =
    1.0;


// ============================================================
// 5. DIMENSIONLESS TIME
//
// T = t / sqrt(d/g)
//
// sqrt(d/g) = 1600
//
// T=1 -> 1600
// T=2 -> 3200
// T=3 -> 4800
// T=4 -> 6400
// ============================================================

const double timeScale =
    sqrt(
        d
        /
        gravity
    );


const int stepT1 =
    static_cast<int>(
        lround(
            timeScale
        )
    );


const int stepT2 =
    static_cast<int>(
        lround(
            2.0
            *
            timeScale
        )
    );


const int stepT3 =
    static_cast<int>(
        lround(
            3.0
            *
            timeScale
        )
    );


const int stepT4 =
    static_cast<int>(
        lround(
            4.0
            *
            timeScale
        )
    );


// ============================================================
// NUMBER OF NODES
// ============================================================

const size_t N =
    static_cast<size_t>(NX)
    *
    NY
    *
    NZ;


// ============================================================
// 6. D3Q7
// ============================================================

const int cx7[Q7] =
{
     0,
     1,-1,
     0, 0,
     0, 0
};


const int cy7[Q7] =
{
     0,
     0, 0,
     1,-1,
     0, 0
};


const int cz7[Q7] =
{
     0,
     0, 0,
     0, 0,
     1,-1
};


const int opp7[Q7] =
{
    0,
    2,1,
    4,3,
    6,5
};


const double w7[Q7] =
{
    1.0 / 4.0,

    1.0 / 8.0,
    1.0 / 8.0,

    1.0 / 8.0,
    1.0 / 8.0,

    1.0 / 8.0,
    1.0 / 8.0
};


// ============================================================
// 7. D3Q15
// ============================================================

const int cx15[Q15] =
{
     0,

     1,-1,
     0, 0,
     0, 0,

     1,-1, 1,-1,
     1,-1, 1,-1
};


const int cy15[Q15] =
{
     0,

     0, 0,
     1,-1,
     0, 0,

     1, 1,-1,-1,
     1, 1,-1,-1
};


const int cz15[Q15] =
{
     0,

     0, 0,
     0, 0,
     1,-1,

     1, 1, 1, 1,
    -1,-1,-1,-1
};


const int opp15[Q15] =
{
     0,

     2,1,
     4,3,
     6,5,

    14,13,12,11,
    10, 9, 8, 7
};


const double w15[Q15] =
{
    2.0 / 9.0,

    1.0 / 9.0,
    1.0 / 9.0,

    1.0 / 9.0,
    1.0 / 9.0,

    1.0 / 9.0,
    1.0 / 9.0,

    1.0 / 72.0,
    1.0 / 72.0,
    1.0 / 72.0,
    1.0 / 72.0,
    1.0 / 72.0,
    1.0 / 72.0,
    1.0 / 72.0,
    1.0 / 72.0
};


// ============================================================
// 8. TRANSFORMATION MATRICES
// ============================================================

const double M7[Q7][Q7] =
{
    {1,1,1,1,1,1,1},

    {0,1,-1,0,0,0,0},

    {0,0,0,1,-1,0,0},

    {0,0,0,0,0,1,-1},

    {0,1,1,0,0,0,0},

    {0,0,0,1,1,0,0},

    {0,0,0,0,0,1,1}
};


const double M15[Q15][Q15] =
{
    {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1},

    {0,1,-1,0,0,0,0,1,-1,1,-1,1,-1,1,-1},

    {0,0,0,1,-1,0,0,1,1,-1,-1,1,1,-1,-1},

    {0,0,0,0,0,1,-1,1,1,1,1,-1,-1,-1,-1},

    {0,1,1,0,0,0,0,1,1,1,1,1,1,1,1},

    {0,0,0,1,1,0,0,1,1,1,1,1,1,1,1},

    {0,0,0,0,0,1,1,1,1,1,1,1,1,1,1},

    {0,0,0,0,0,0,0,1,-1,-1,1,1,-1,-1,1},

    {0,0,0,0,0,0,0,1,1,-1,-1,-1,-1,1,1},

    {0,0,0,0,0,0,0,1,-1,1,-1,-1,1,-1,1},

    {0,0,0,0,0,0,0,1,1,-1,-1,1,1,-1,-1},

    {0,0,0,0,0,0,0,1,1,1,1,-1,-1,-1,-1},

    {0,0,0,0,0,0,0,1,-1,1,-1,1,-1,1,-1},

    {0,0,0,0,0,0,0,1,-1,-1,1,-1,1,1,-1},

    {0,0,0,0,0,0,0,1,1,1,1,1,1,1,1}
};


double Minv7[Q7][Q7];

double Minv15[Q15][Q15];


// ============================================================
// 9. INDEX
// ============================================================

inline size_t id(
    int x,
    int y,
    int z
)
{
    return
        (
            (
                static_cast<size_t>(z)
                *
                NY
                +
                y
            )
            *
            NX
            +
            x
        );
}


inline size_t phaseID(
    size_t n,
    int q
)
{
    return
        n
        *
        Q7
        +
        q;
}


inline size_t flowID(
    size_t n,
    int alpha,
    int q
)
{
    return
        (
            n
            *
            3
            +
            alpha
        )
        *
        Q15
        +
        q;
}


// ============================================================
// 10. FIELDS
// ============================================================

vector<double> phi(N);

vector<double> rho(N);

vector<double> chemical(N);

vector<double> gradPhiX(N);
vector<double> gradPhiY(N);
vector<double> gradPhiZ(N);

vector<double> gradRhoX(N);
vector<double> gradRhoY(N);
vector<double> gradRhoZ(N);

vector<double> forceX(N);
vector<double> forceY(N);
vector<double> forceZ(N);

vector<double> ux(N);
vector<double> uy(N);
vector<double> uz(N);

vector<double> pressure(N);


// gradU:
//
// gradU[(component*3+direction)*N+n]
//
// = partial_direction u_component

vector<double> gradU(
    9
    *
    N
);


// old phi*u

vector<double> oldPhiUx(N);
vector<double> oldPhiUy(N);
vector<double> oldPhiUz(N);


// old Q tensor:
//
// 0 = xx
// 1 = xy
// 2 = xz
// 3 = yy
// 4 = yz
// 5 = zz

vector<double> oldQ(
    6
    *
    N
);


// Distribution functions

vector<double> phaseDist(
    N
    *
    Q7
);

vector<double> phaseNext(
    N
    *
    Q7
);


vector<double> flowDist(
    N
    *
    3
    *
    Q15
);

vector<double> flowNext(
    N
    *
    3
    *
    Q15
);


double initialMass =
    0.0;


// ============================================================
// 11. MATRIX INVERSE
// ============================================================

template<int NN>
void invertMatrix(
    const double A[NN][NN],
    double invA[NN][NN]
)
{
    double aug[NN][2 * NN];


    for(int i = 0;
        i < NN;
        i++)
    {
        for(int j = 0;
            j < NN;
            j++)
        {
            aug[i][j] =
                A[i][j];


            aug[i][j + NN] =
                (i == j)
                ?
                1.0
                :
                0.0;
        }
    }


    for(int col = 0;
        col < NN;
        col++)
    {
        int pivot =
            col;


        double maxValue =
            fabs(
                aug[col][col]
            );


        for(int row = col + 1;
            row < NN;
            row++)
        {
            double value =
                fabs(
                    aug[row][col]
                );


            if(value > maxValue)
            {
                maxValue =
                    value;

                pivot =
                    row;
            }
        }


        if(maxValue < 1.0e-14)
        {
            cerr
                << "Matrix inversion failed.\n";

            exit(1);
        }


        if(pivot != col)
        {
            for(int j = 0;
                j < 2 * NN;
                j++)
            {
                swap(
                    aug[pivot][j],
                    aug[col][j]
                );
            }
        }


        double divisor =
            aug[col][col];


        for(int j = 0;
            j < 2 * NN;
            j++)
        {
            aug[col][j]
                /=
            divisor;
        }


        for(int row = 0;
            row < NN;
            row++)
        {
            if(row == col)
            {
                continue;
            }


            double factor =
                aug[row][col];


            for(int j = 0;
                j < 2 * NN;
                j++)
            {
                aug[row][j]
                    -=
                factor
                *
                aug[col][j];
            }
        }
    }


    for(int i = 0;
        i < NN;
        i++)
    {
        for(int j = 0;
            j < NN;
            j++)
        {
            invA[i][j] =
                aug[i][j + NN];
        }
    }
}


// ============================================================
// 12. MATRIX CHECK
// ============================================================

void checkMatrices()
{
    double err7 =
        0.0;

    double err15 =
        0.0;


    for(int i = 0;
        i < Q7;
        i++)
    {
        for(int j = 0;
            j < Q7;
            j++)
        {
            double value =
                0.0;


            for(int k = 0;
                k < Q7;
                k++)
            {
                value
                    +=
                M7[i][k]
                *
                Minv7[k][j];
            }


            double exact =
                (i == j)
                ?
                1.0
                :
                0.0;


            err7 =
                max(
                    err7,
                    fabs(
                        value
                        -
                        exact
                    )
                );
        }
    }


    for(int i = 0;
        i < Q15;
        i++)
    {
        for(int j = 0;
            j < Q15;
            j++)
        {
            double value =
                0.0;


            for(int k = 0;
                k < Q15;
                k++)
            {
                value
                    +=
                M15[i][k]
                *
                Minv15[k][j];
            }


            double exact =
                (i == j)
                ?
                1.0
                :
                0.0;


            err15 =
                max(
                    err15,
                    fabs(
                        value
                        -
                        exact
                    )
                );
        }
    }


    cout
        << "M7 inverse error  = "
        << err7
        << "\n";


    cout
        << "M15 inverse error = "
        << err15
        << "\n";
}


// ============================================================
// 13. BOUNDARY HELPERS
// ============================================================

inline int periodicX(
    int x
)
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


inline int periodicY(
    int y
)
{
    if(y < 0)
    {
        return y + NY;
    }


    if(y >= NY)
    {
        return y - NY;
    }


    return y;
}


inline int clampZ(
    int z
)
{
    if(z < 0)
    {
        return 0;
    }


    if(z >= NZ)
    {
        return NZ - 1;
    }


    return z;
}


// ============================================================
// 14. INITIAL FIELD
// ============================================================

void initializeFields()
{
    #pragma omp parallel for collapse(3)

    for(int z = 0;
        z < NZ;
        z++)
    {
        for(int y = 0;
            y < NY;
            y++)
        {
            for(int x = 0;
                x < NX;
                x++)
            {
                size_t n =
                    id(
                        x,
                        y,
                        z
                    );


                double xx =
                    -0.5
                    *
                    d

                    +

                    (
                        x
                        +
                        0.5
                    )
                    *
                    dx;


                double yy =
                    -0.5
                    *
                    d

                    +

                    (
                        y
                        +
                        0.5
                    )
                    *
                    dx;


                double zz =
                    -2.0
                    *
                    d

                    +

                    (
                        z
                        +
                        0.5
                    )
                    *
                    dx;


                double h =
                    0.05
                    *
                    d
                    *
                    (
                        cos(
                            2.0
                            *
                            PI
                            *
                            xx
                            /
                            d
                        )

                        +

                        cos(
                            2.0
                            *
                            PI
                            *
                            yy
                            /
                            d
                        )
                    );


                phi[n] =
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
                            zz
                            -
                            h
                        )
                        /
                        W
                    );


                ux[n] =
                    0.0;


                uy[n] =
                    0.0;


                uz[n] =
                    0.0;


                pressure[n] =
                    0.0;
            }
        }
    }
}


// ============================================================
// 15. DENSITY
// ============================================================

void updateDensity()
{
    #pragma omp parallel for

    for(size_t n = 0;
        n < N;
        n++)
    {
        double value =
            phi[n];


        // Bound only property interpolation.

        if(value < phiB)
        {
            value =
                phiB;
        }


        if(value > phiA)
        {
            value =
                phiA;
        }


        rho[n] =
            rhoB

            +

            (
                rhoA
                -
                rhoB
            )
            *
            (
                value
                -
                phiB
            )
            /
            (
                phiA
                -
                phiB
            );
    }
}


// ============================================================
// 16. DERIVATIVES
// ============================================================

void computeDerivatives()
{
    #pragma omp parallel for collapse(3)

    for(int z = 0;
        z < NZ;
        z++)
    {
        for(int y = 0;
            y < NY;
            y++)
        {
            for(int x = 0;
                x < NX;
                x++)
            {
                size_t n =
                    id(
                        x,
                        y,
                        z
                    );


                double gpx =
                    0.0;

                double gpy =
                    0.0;

                double gpz =
                    0.0;


                double grx =
                    0.0;

                double gry =
                    0.0;

                double grz =
                    0.0;


                double lapPhi =
                    0.0;


                for(int q = 1;
                    q < Q15;
                    q++)
                {
                    int xn =
                        periodicX(
                            x
                            +
                            cx15[q]
                        );


                    int yn =
                        periodicY(
                            y
                            +
                            cy15[q]
                        );


                    int zn =
                        clampZ(
                            z
                            +
                            cz15[q]
                        );


                    size_t nn =
                        id(
                            xn,
                            yn,
                            zn
                        );


                    double cix =
                        static_cast<double>(
                            cx15[q]
                        );


                    double ciy =
                        static_cast<double>(
                            cy15[q]
                        );


                    double ciz =
                        static_cast<double>(
                            cz15[q]
                        );


                    gpx
                        +=
                    w15[q]
                    *
                    cix
                    *
                    phi[nn];


                    gpy
                        +=
                    w15[q]
                    *
                    ciy
                    *
                    phi[nn];


                    gpz
                        +=
                    w15[q]
                    *
                    ciz
                    *
                    phi[nn];


                    grx
                        +=
                    w15[q]
                    *
                    cix
                    *
                    rho[nn];


                    gry
                        +=
                    w15[q]
                    *
                    ciy
                    *
                    rho[nn];


                    grz
                        +=
                    w15[q]
                    *
                    ciz
                    *
                    rho[nn];


                    lapPhi
                        +=
                    2.0
                    *
                    w15[q]
                    *
                    (
                        phi[nn]
                        -
                        phi[n]
                    );
                }


                gradPhiX[n] =
                    gpx
                    /
                    (
                        cs2Flow
                        *
                        dt
                    );


                gradPhiY[n] =
                    gpy
                    /
                    (
                        cs2Flow
                        *
                        dt
                    );


                gradPhiZ[n] =
                    gpz
                    /
                    (
                        cs2Flow
                        *
                        dt
                    );


                gradRhoX[n] =
                    grx
                    /
                    (
                        cs2Flow
                        *
                        dt
                    );


                gradRhoY[n] =
                    gry
                    /
                    (
                        cs2Flow
                        *
                        dt
                    );


                gradRhoZ[n] =
                    grz
                    /
                    (
                        cs2Flow
                        *
                        dt
                    );


                lapPhi
                    /=
                (
                    cs2Flow
                    *
                    dt
                    *
                    dt
                );


                chemical[n] =
                    4.0
                    *
                    betaPF
                    *
                    (
                        phi[n]
                        -
                        phiA
                    )
                    *
                    (
                        phi[n]
                        -
                        phiB
                    )
                    *
                    (
                        phi[n]

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
}


// ============================================================
// 17. FORCE
// ============================================================

void computeForce()
{
    #pragma omp parallel for

    for(size_t n = 0;
        n < N;
        n++)
    {
        if(USE_SURFACE_FORCE)
        {
            forceX[n] =
                chemical[n]
                *
                gradPhiX[n];


            forceY[n] =
                chemical[n]
                *
                gradPhiY[n];


            forceZ[n] =
                chemical[n]
                *
                gradPhiZ[n]

                -

                (
                    rho[n]
                    -
                    rhoRef
                )
                *
                gravity;
        }
        else
        {
            forceX[n] =
                0.0;


            forceY[n] =
                0.0;


            forceZ[n] =
                -
                (
                    rho[n]
                    -
                    rhoRef
                )
                *
                gravity;
        }
    }
}


// ============================================================
// 18. HYDROSTATIC INITIAL PRESSURE
// ============================================================

void initializeHydrostaticPressure()
{
    #pragma omp parallel for collapse(2)

    for(int y = 0;
        y < NY;
        y++)
    {
        for(int x = 0;
            x < NX;
            x++)
        {
            size_t n0 =
                id(
                    x,
                    y,
                    0
                );


            pressure[n0] =
                0.0;


            for(int z = 1;
                z < NZ;
                z++)
            {
                size_t nBelow =
                    id(
                        x,
                        y,
                        z - 1
                    );


                size_t n =
                    id(
                        x,
                        y,
                        z
                    );


                pressure[n] =
                    pressure[nBelow]

                    +

                    0.5
                    *
                    (
                        forceZ[nBelow]
                        +
                        forceZ[n]
                    )
                    *
                    dx;
            }
        }
    }
}


// ============================================================
// 19. PHASE EQUILIBRIUM
// ============================================================

inline double phaseEquilibrium(
    size_t n,
    int q
)
{
    double cu =
        cx7[q]
        *
        ux[n]

        +

        cy7[q]
        *
        uy[n]

        +

        cz7[q]
        *
        uz[n];


    if(q == 0)
    {
        return
            phi[n]

            +

            (
                w7[q]
                -
                1.0
            )
            *
            eta
            *
            phi[n];
    }


    return
        w7[q]
        *
        eta
        *
        phi[n]

        +

        w7[q]
        *
        phi[n]
        *
        cu
        /
        cs2Phase;
}


// ============================================================
// 20. FLOW EQUILIBRIUM
// ============================================================

inline double flowEquilibrium(
    size_t n,
    int alpha,
    int q
)
{
    double u[3] =
    {
        ux[n],
        uy[n],
        uz[n]
    };


    double ci[3] =
    {
        static_cast<double>(
            cx15[q]
        ),

        static_cast<double>(
            cy15[q]
        ),

        static_cast<double>(
            cz15[q]
        )
    };


    double value =
        rho[n]
        *
        u[alpha];


    for(int beta = 0;
        beta < 3;
        beta++)
    {
        double Qab =
            rho[n]
            *
            u[alpha]
            *
            u[beta];


        if(alpha == beta)
        {
            Qab
                +=
            pressure[n];
        }


        value
            +=
        ci[beta]
        *
        Qab
        /
        cs2Flow;
    }


    return
        w15[q]
        *
        value;
}


// ============================================================
// 21. INITIAL DISTRIBUTIONS
// ============================================================

void initializeDistributions()
{
    #pragma omp parallel for

    for(size_t n = 0;
        n < N;
        n++)
    {
        for(int q = 0;
            q < Q7;
            q++)
        {
            phaseDist[
                phaseID(
                    n,
                    q
                )
            ]
                =
            phaseEquilibrium(
                n,
                q
            );
        }


        double F[3] =
        {
            forceX[n],
            forceY[n],
            forceZ[n]
        };


        for(int alpha = 0;
            alpha < 3;
            alpha++)
        {
            for(int q = 0;
                q < Q15;
                q++)
            {
                double Fi =
                    w15[q]
                    *
                    F[alpha];


                flowDist[
                    flowID(
                        n,
                        alpha,
                        q
                    )
                ]
                    =
                flowEquilibrium(
                    n,
                    alpha,
                    q
                )

                -

                0.5
                *
                dt
                *
                Fi;
            }
        }
    }
}


// ============================================================
// 22. RECOVER PHI
// ============================================================

void recoverPhi()
{
    #pragma omp parallel for

    for(size_t n = 0;
        n < N;
        n++)
    {
        double value =
            0.0;


        for(int q = 0;
            q < Q7;
            q++)
        {
            value
                +=
            phaseDist[
                phaseID(
                    n,
                    q
                )
            ];
        }


        phi[n] =
            value;
    }
}


// ============================================================
// 23. RECOVER VELOCITY
// ============================================================

void recoverVelocity()
{
    #pragma omp parallel for

    for(size_t n = 0;
        n < N;
        n++)
    {
        double sum[3] =
        {
            0.0,
            0.0,
            0.0
        };


        for(int alpha = 0;
            alpha < 3;
            alpha++)
        {
            for(int q = 0;
                q < Q15;
                q++)
            {
                sum[alpha]
                    +=
                flowDist[
                    flowID(
                        n,
                        alpha,
                        q
                    )
                ];
            }
        }


        ux[n] =
            (
                sum[0]

                +

                0.5
                *
                dt
                *
                forceX[n]
            )
            /
            rho[n];


        uy[n] =
            (
                sum[1]

                +

                0.5
                *
                dt
                *
                forceY[n]
            )
            /
            rho[n];


        uz[n] =
            (
                sum[2]

                +

                0.5
                *
                dt
                *
                forceZ[n]
            )
            /
            rho[n];
    }
}


// ============================================================
// 24. RECOVER PRESSURE
// ============================================================

void recoverPressure()
{
    #pragma omp parallel for

    for(size_t n = 0;
        n < N;
        n++)
    {
        double traceMoment =
            0.0;


        for(int q = 0;
            q < Q15;
            q++)
        {
            traceMoment
                +=
            cx15[q]
            *
            flowDist[
                flowID(
                    n,
                    0,
                    q
                )
            ];


            traceMoment
                +=
            cy15[q]
            *
            flowDist[
                flowID(
                    n,
                    1,
                    q
                )
            ];


            traceMoment
                +=
            cz15[q]
            *
            flowDist[
                flowID(
                    n,
                    2,
                    q
                )
            ];
        }


        double u2 =
            ux[n]
            *
            ux[n]

            +

            uy[n]
            *
            uy[n]

            +

            uz[n]
            *
            uz[n];


        double uGradRho =
            ux[n]
            *
            gradRhoX[n]

            +

            uy[n]
            *
            gradRhoY[n]

            +

            uz[n]
            *
            gradRhoZ[n];


        pressure[n] =
            (
                traceMoment

                -

                rho[n]
                *
                u2

                +

                0.5
                *
                dt
                *
                cs2Flow
                *
                uGradRho
            )
            /
            3.0;
    }
}


// ============================================================
// 25. VELOCITY GRADIENT
// ============================================================

void computeVelocityGradient()
{
    #pragma omp parallel for collapse(3)

    for(int z = 0;
        z < NZ;
        z++)
    {
        for(int y = 0;
            y < NY;
            y++)
        {
            for(int x = 0;
                x < NX;
                x++)
            {
                size_t n =
                    id(
                        x,
                        y,
                        z
                    );


                double grad[3][3] =
                {
                    {0,0,0},
                    {0,0,0},
                    {0,0,0}
                };


                for(int q = 1;
                    q < Q15;
                    q++)
                {
                    int xn =
                        periodicX(
                            x
                            +
                            cx15[q]
                        );


                    int yn =
                        periodicY(
                            y
                            +
                            cy15[q]
                        );


                    int zn =
                        clampZ(
                            z
                            +
                            cz15[q]
                        );


                    size_t nn =
                        id(
                            xn,
                            yn,
                            zn
                        );


                    double uq[3] =
                    {
                        ux[nn],
                        uy[nn],
                        uz[nn]
                    };


                    double ci[3] =
                    {
                        static_cast<double>(
                            cx15[q]
                        ),

                        static_cast<double>(
                            cy15[q]
                        ),

                        static_cast<double>(
                            cz15[q]
                        )
                    };


                    for(int alpha = 0;
                        alpha < 3;
                        alpha++)
                    {
                        for(int beta = 0;
                            beta < 3;
                            beta++)
                        {
                            grad[alpha][beta]
                                +=
                            w15[q]
                            *
                            ci[beta]
                            *
                            uq[alpha];
                        }
                    }
                }


                for(int alpha = 0;
                    alpha < 3;
                    alpha++)
                {
                    for(int beta = 0;
                        beta < 3;
                        beta++)
                    {
                        gradU[
                            (
                                alpha
                                *
                                3
                                +
                                beta
                            )
                            *
                            N
                            +
                            n
                        ]
                            =
                        grad[alpha][beta]
                        /
                        (
                            cs2Flow
                            *
                            dt
                        );
                    }
                }
            }
        }
    }
}


// ============================================================
// 26. SAVE OLD FIELDS
// ============================================================

void saveOldFields()
{
    #pragma omp parallel for

    for(size_t n = 0;
        n < N;
        n++)
    {
        oldPhiUx[n] =
            phi[n]
            *
            ux[n];


        oldPhiUy[n] =
            phi[n]
            *
            uy[n];


        oldPhiUz[n] =
            phi[n]
            *
            uz[n];


        oldQ[
            0 * N + n
        ]
            =
        rho[n]
        *
        ux[n]
        *
        ux[n]

        +

        pressure[n];


        oldQ[
            1 * N + n
        ]
            =
        rho[n]
        *
        ux[n]
        *
        uy[n];


        oldQ[
            2 * N + n
        ]
            =
        rho[n]
        *
        ux[n]
        *
        uz[n];


        oldQ[
            3 * N + n
        ]
            =
        rho[n]
        *
        uy[n]
        *
        uy[n]

        +

        pressure[n];


        oldQ[
            4 * N + n
        ]
            =
        rho[n]
        *
        uy[n]
        *
        uz[n];


        oldQ[
            5 * N + n
        ]
            =
        rho[n]
        *
        uz[n]
        *
        uz[n]

        +

        pressure[n];
    }
}


// ============================================================
// 27. PHASE COLLISION + STREAMING
// ============================================================

void collideStreamPhase()
{
    #pragma omp parallel for collapse(3)

    for(int z = 0;
        z < NZ;
        z++)
    {
        for(int y = 0;
            y < NY;
            y++)
        {
            for(int x = 0;
                x < NX;
                x++)
            {
                size_t n =
                    id(
                        x,
                        y,
                        z
                    );


                double dPhiU[3] =
                {
                    (
                        phi[n]
                        *
                        ux[n]

                        -

                        oldPhiUx[n]
                    )
                    /
                    dt,


                    (
                        phi[n]
                        *
                        uy[n]

                        -

                        oldPhiUy[n]
                    )
                    /
                    dt,


                    (
                        phi[n]
                        *
                        uz[n]

                        -

                        oldPhiUz[n]
                    )
                    /
                    dt
                };


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
                        phi[n]
                    )

                    *
                    (
                        phi[n]
                        -
                        phiB
                    );


                // ==========================================
                // Eq.(43)
                // ==========================================

                double C[3] =
                {
                    0.0,
                    0.0,
                    0.0
                };


                for(int q = 0;
                    q < Q7;
                    q++)
                {
                    double f =
                        phaseDist[
                            phaseID(
                                n,
                                q
                            )
                        ];


                    C[0]
                        +=
                    cx7[q]
                    *
                    f;


                    C[1]
                        +=
                    cy7[q]
                    *
                    f;


                    C[2]
                        +=
                    cz7[q]
                    *
                    f;
                }


                C[0] =
                    C[0]

                    -

                    phi[n]
                    *
                    ux[n]

                    +

                    0.5
                    *
                    dt
                    *
                    dPhiU[0];


                C[1] =
                    C[1]

                    -

                    phi[n]
                    *
                    uy[n]

                    +

                    0.5
                    *
                    dt
                    *
                    dPhiU[1];


                C[2] =
                    C[2]

                    -

                    phi[n]
                    *
                    uz[n]

                    +

                    0.5
                    *
                    dt
                    *
                    dPhiU[2];


                double A =
                    -mobility

                    -

                    0.5
                    *
                    eta
                    *
                    cs2Phase
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
                        C[0]
                        *
                        C[0]

                        +

                        C[1]
                        *
                        C[1]

                        +

                        C[2]
                        *
                        C[2]
                    );


                double gradMag =
                    (
                        -Cmag
                        -
                        B
                    )
                    /
                    A;


                double normal[3] =
                {
                    0.0,
                    0.0,
                    0.0
                };


                if(gradMag > epsNormal)
                {
                    double denominator =
                        A
                        +
                        B
                        /
                        gradMag;


                    if(
                        fabs(
                            denominator
                        )
                        >
                        epsNormal
                    )
                    {
                        double gx =
                            C[0]
                            /
                            denominator;


                        double gy =
                            C[1]
                            /
                            denominator;


                        double gz =
                            C[2]
                            /
                            denominator;


                        double norm =
                            sqrt(
                                gx * gx
                                +
                                gy * gy
                                +
                                gz * gz
                            );


                        if(norm > epsNormal)
                        {
                            normal[0] =
                                gx
                                /
                                norm;


                            normal[1] =
                                gy
                                /
                                norm;


                            normal[2] =
                                gz
                                /
                                norm;
                        }
                    }
                }


                double mCurrent[3] =
                {
                    0.0,
                    0.0,
                    0.0
                };


                for(int q = 0;
                    q < Q7;
                    q++)
                {
                    double f =
                        phaseDist[
                            phaseID(
                                n,
                                q
                            )
                        ];


                    mCurrent[0]
                        +=
                    cx7[q]
                    *
                    f;


                    mCurrent[1]
                        +=
                    cy7[q]
                    *
                    f;


                    mCurrent[2]
                        +=
                    cz7[q]
                    *
                    f;
                }


                double mEq[3] =
                {
                    phi[n]
                    *
                    ux[n],

                    phi[n]
                    *
                    uy[n],

                    phi[n]
                    *
                    uz[n]
                };


                double mG[3] =
                {
                    dPhiU[0]

                    +

                    eta
                    *
                    lambdaACE
                    *
                    cs2Phase
                    *
                    normal[0],


                    dPhiU[1]

                    +

                    eta
                    *
                    lambdaACE
                    *
                    cs2Phase
                    *
                    normal[1],


                    dPhiU[2]

                    +

                    eta
                    *
                    lambdaACE
                    *
                    cs2Phase
                    *
                    normal[2]
                };


                for(int q = 0;
                    q < Q7;
                    q++)
                {
                    double feq =
                        phaseEquilibrium(
                            n,
                            q
                        );


                    double ciDotDPhiU =
                        cx7[q]
                        *
                        dPhiU[0]

                        +

                        cy7[q]
                        *
                        dPhiU[1]

                        +

                        cz7[q]
                        *
                        dPhiU[2];


                    double ciDotNormal =
                        cx7[q]
                        *
                        normal[0]

                        +

                        cy7[q]
                        *
                        normal[1]

                        +

                        cz7[q]
                        *
                        normal[2];


                    double Gi =
                        w7[q]
                        *
                        (
                            ciDotDPhiU
                            /
                            cs2Phase

                            +

                            eta
                            *
                            lambdaACE
                            *
                            ciDotNormal
                        );


                    double correctionNoneq =
                        0.0;


                    double correctionSource =
                        0.0;


                    for(int k = 0;
                        k < 3;
                        k++)
                    {
                        correctionNoneq
                            +=
                        Minv7[q][k + 1]
                        *
                        (
                            mCurrent[k]
                            -
                            mEq[k]
                        );


                        correctionSource
                            +=
                        Minv7[q][k + 1]
                        *
                        mG[k];
                    }


                    double post =
                        feq

                        +

                        (
                            1.0
                            -
                            sPhase1
                        )
                        *
                        correctionNoneq

                        +

                        dt
                        *
                        (
                            0.5
                            *
                            Gi

                            +

                            0.5
                            *
                            (
                                1.0
                                -
                                sPhase1
                            )
                            *
                            correctionSource
                        );


                    int xn =
                        periodicX(
                            x
                            +
                            cx7[q]
                        );


                    int yn =
                        periodicY(
                            y
                            +
                            cy7[q]
                        );


                    int zn =
                        z
                        +
                        cz7[q];


                    if(
                        zn >= 0
                        &&
                        zn < NZ
                    )
                    {
                        size_t nn =
                            id(
                                xn,
                                yn,
                                zn
                            );


                        phaseNext[
                            phaseID(
                                nn,
                                q
                            )
                        ]
                            =
                        post;
                    }
                    else
                    {
                        phaseNext[
                            phaseID(
                                n,
                                opp7[q]
                            )
                        ]
                            =
                        post;
                    }
                }
            }
        }
    }


    phaseDist.swap(
        phaseNext
    );
}


// ============================================================
// 28. FLOW COLLISION + STREAMING
//
// Full D3Q15 MRT
//
// Stable Eq.(37):
//
// Gi = wi [
//     u_alpha c_i beta grad_beta rho
//     - rho c_i beta grad_alpha u_beta
// ]
//
// dQ/dt is NOT inserted.
// ============================================================

void collideStreamFlow()
{
    #pragma omp parallel for collapse(3)

    for(int z = 0;
        z < NZ;
        z++)
    {
        for(int y = 0;
            y < NY;
            y++)
        {
            for(int x = 0;
                x < NX;
                x++)
            {
                size_t n =
                    id(
                        x,
                        y,
                        z
                    );


                double u[3] =
                {
                    ux[n],
                    uy[n],
                    uz[n]
                };


                double gradRho[3] =
                {
                    gradRhoX[n],
                    gradRhoY[n],
                    gradRhoZ[n]
                };


                double F[3] =
                {
                    forceX[n],
                    forceY[n],
                    forceZ[n]
                };


                // ------------------------------------------
                // Q tensor
                // ------------------------------------------

                double Qxx =
                    rho[n]
                    *
                    u[0]
                    *
                    u[0]

                    +

                    pressure[n];


                double Qxy =
                    rho[n]
                    *
                    u[0]
                    *
                    u[1];


                double Qxz =
                    rho[n]
                    *
                    u[0]
                    *
                    u[2];


                double Qyy =
                    rho[n]
                    *
                    u[1]
                    *
                    u[1]

                    +

                    pressure[n];


                double Qyz =
                    rho[n]
                    *
                    u[1]
                    *
                    u[2];


                double Qzz =
                    rho[n]
                    *
                    u[2]
                    *
                    u[2]

                    +

                    pressure[n];


                double dQ[6] =
                {
                    (
                        Qxx
                        -
                        oldQ[
                            0 * N + n
                        ]
                    )
                    /
                    dt,


                    (
                        Qxy
                        -
                        oldQ[
                            1 * N + n
                        ]
                    )
                    /
                    dt,


                    (
                        Qxz
                        -
                        oldQ[
                            2 * N + n
                        ]
                    )
                    /
                    dt,


                    (
                        Qyy
                        -
                        oldQ[
                            3 * N + n
                        ]
                    )
                    /
                    dt,


                    (
                        Qyz
                        -
                        oldQ[
                            4 * N + n
                        ]
                    )
                    /
                    dt,


                    (
                        Qzz
                        -
                        oldQ[
                            5 * N + n
                        ]
                    )
                    /
                    dt
                };


                double Sflow[Q15] =
                {
                    su0,

                    su1,
                    su1,
                    su1,

                    su2,
                    su2,
                    su2,
                    su2,
                    su2,
                    su2,

                    su3,
                    su3,
                    su3,
                    su3,

                    su4
                };


                for(int alpha = 0;
                    alpha < 3;
                    alpha++)
                {
                    double dQrow[3];


                    if(alpha == 0)
                    {
                        dQrow[0] =
                            dQ[0];

                        dQrow[1] =
                            dQ[1];

                        dQrow[2] =
                            dQ[2];
                    }
                    else if(alpha == 1)
                    {
                        dQrow[0] =
                            dQ[1];

                        dQrow[1] =
                            dQ[3];

                        dQrow[2] =
                            dQ[4];
                    }
                    else
                    {
                        dQrow[0] =
                            dQ[2];

                        dQrow[1] =
                            dQ[4];

                        dQrow[2] =
                            dQ[5];
                    }


                    double gNow[Q15];

                    double feq[Q15];

                    double Fi[Q15];

                    double Gi[Q15];


                    for(int q = 0;
                        q < Q15;
                        q++)
                    {
                        gNow[q] =
                            flowDist[
                                flowID(
                                    n,
                                    alpha,
                                    q
                                )
                            ];


                        feq[q] =
                            flowEquilibrium(
                                n,
                                alpha,
                                q
                            );


                        Fi[q] =
                            w15[q]
                            *
                            F[alpha];


                        double ci[3] =
                        {
                            static_cast<double>(
                                cx15[q]
                            ),

                            static_cast<double>(
                                cy15[q]
                            ),

                            static_cast<double>(
                                cz15[q]
                            )
                        };


                        // ----------------------------------
                        // A: dQ/dt
                        //
                        // calculated but currently NOT used
                        // ----------------------------------

                        double ciDotDQ =
                            ci[0]
                            *
                            dQrow[0]

                            +

                            ci[1]
                            *
                            dQrow[1]

                            +

                            ci[2]
                            *
                            dQrow[2];


                        double termA =
                            ciDotDQ
                            /
                            cs2Flow;


                        // ----------------------------------
                        // B
                        // ----------------------------------

                        double ciDotGradRho =
                            ci[0]
                            *
                            gradRho[0]

                            +

                            ci[1]
                            *
                            gradRho[1]

                            +

                            ci[2]
                            *
                            gradRho[2];


                        double termB =
                            u[alpha]
                            *
                            ciDotGradRho;


                        // ----------------------------------
                        // C
                        // ----------------------------------

                        double ciDotGradU =
                            0.0;


                        for(int beta = 0;
                            beta < 3;
                            beta++)
                        {
                            double partialAlphaUBeta =
                                gradU[
                                    (
                                        beta
                                        *
                                        3
                                        +
                                        alpha
                                    )
                                    *
                                    N
                                    +
                                    n
                                ];


                            ciDotGradU
                                +=
                            ci[beta]
                            *
                            partialAlphaUBeta;
                        }


                        double termC =
                            -
                            rho[n]
                            *
                            ciDotGradU;


                        double totalG =
                            termB
                            +
                            termC;


                        if(USE_DQDT_IN_G)
                        {
                            totalG
                                +=
                            termA;
                        }


                        Gi[q] =
                            w15[q]
                            *
                            totalG;
                    }


                    // --------------------------------------
                    // Moment transform
                    // --------------------------------------

                    double m[Q15] =
                    {
                        0.0
                    };


                    double meq[Q15] =
                    {
                        0.0
                    };


                    double mF[Q15] =
                    {
                        0.0
                    };


                    double mG[Q15] =
                    {
                        0.0
                    };


                    for(int k = 0;
                        k < Q15;
                        k++)
                    {
                        for(int q = 0;
                            q < Q15;
                            q++)
                        {
                            m[k]
                                +=
                            M15[k][q]
                            *
                            gNow[q];


                            meq[k]
                                +=
                            M15[k][q]
                            *
                            feq[q];


                            mF[k]
                                +=
                            M15[k][q]
                            *
                            Fi[q];


                            mG[k]
                                +=
                            M15[k][q]
                            *
                            Gi[q];
                        }
                    }


                    // --------------------------------------
                    // MRT collision
                    // --------------------------------------

                    double mPost[Q15];


                    for(int k = 0;
                        k < Q15;
                        k++)
                    {
                        double s =
                            Sflow[k];


                        mPost[k] =
                            (
                                1.0
                                -
                                s
                            )
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
                            (
                                mF[k]
                                +
                                mG[k]
                            );
                    }


                    // --------------------------------------
                    // Inverse transform
                    // --------------------------------------

                    double post[Q15];


                    for(int q = 0;
                        q < Q15;
                        q++)
                    {
                        post[q] =
                            0.0;


                        for(int k = 0;
                            k < Q15;
                            k++)
                        {
                            post[q]
                                +=
                            Minv15[q][k]
                            *
                            mPost[k];
                        }
                    }


                    // --------------------------------------
                    // Streaming
                    // --------------------------------------

                    for(int q = 0;
                        q < Q15;
                        q++)
                    {
                        int xn =
                            periodicX(
                                x
                                +
                                cx15[q]
                            );


                        int yn =
                            periodicY(
                                y
                                +
                                cy15[q]
                            );


                        int zn =
                            z
                            +
                            cz15[q];


                        if(
                            zn >= 0
                            &&
                            zn < NZ
                        )
                        {
                            size_t nn =
                                id(
                                    xn,
                                    yn,
                                    zn
                                );


                            flowNext[
                                flowID(
                                    nn,
                                    alpha,
                                    q
                                )
                            ]
                                =
                            post[q];
                        }
                        else
                        {
                            flowNext[
                                flowID(
                                    n,
                                    alpha,
                                    opp15[q]
                                )
                            ]
                                =
                            -post[q]

                            -

                            dt
                            *
                            w15[q]
                            *
                            F[alpha];
                        }
                    }
                }
            }
        }
    }


    flowDist.swap(
        flowNext
    );
}


// ============================================================
// 29. NaN / INF CHECK
// ============================================================

bool checkNaN(
    int step
)
{
    for(size_t n = 0;
        n < N;
        n++)
    {
        if(
            !isfinite(
                phi[n]
            )
            ||
            !isfinite(
                rho[n]
            )
            ||
            !isfinite(
                ux[n]
            )
            ||
            !isfinite(
                uy[n]
            )
            ||
            !isfinite(
                uz[n]
            )
            ||
            !isfinite(
                pressure[n]
            )
        )
        {
            size_t xy =
                static_cast<size_t>(
                    NX
                    *
                    NY
                );


            int z =
                static_cast<int>(
                    n
                    /
                    xy
                );


            size_t remain =
                n
                -
                static_cast<size_t>(
                    z
                )
                *
                xy;


            int y =
                static_cast<int>(
                    remain
                    /
                    NX
                );


            int x =
                static_cast<int>(
                    remain
                    %
                    NX
                );


            cout
                << "\n========================================\n";


            cout
                << "NaN / Inf detected\n";


            cout
                << "step = "
                << step
                << "\n";


            cout
                << "T = "
                << step
                /
                timeScale
                << "\n";


            cout
                << "x,y,z = "
                << x
                << ", "
                << y
                << ", "
                << z
                << "\n";


            cout
                << "phi = "
                << phi[n]
                << "\n";


            cout
                << "rho = "
                << rho[n]
                << "\n";


            cout
                << "u = ("
                << ux[n]
                << ", "
                << uy[n]
                << ", "
                << uz[n]
                << ")\n";


            cout
                << "p = "
                << pressure[n]
                << "\n";


            cout
                << "========================================\n";


            return true;
        }
    }


    return false;
}


// ============================================================
// 30. DIAGNOSTICS
// ============================================================

void diagnostics(
    int step
)
{
    double phiMin =
        phi[0];

    double phiMax =
        phi[0];


    double pMin =
        pressure[0];

    double pMax =
        pressure[0];


    double maxSpeed =
        0.0;


    size_t maxSpeedID =
        0;


    double maxChemical =
        0.0;


    double maxForce =
        0.0;


    double mass =
        0.0;


    for(size_t n = 0;
        n < N;
        n++)
    {
        phiMin =
            min(
                phiMin,
                phi[n]
            );


        phiMax =
            max(
                phiMax,
                phi[n]
            );


        pMin =
            min(
                pMin,
                pressure[n]
            );


        pMax =
            max(
                pMax,
                pressure[n]
            );


        double speed =
            sqrt(
                ux[n]
                *
                ux[n]

                +

                uy[n]
                *
                uy[n]

                +

                uz[n]
                *
                uz[n]
            );


        if(speed > maxSpeed)
        {
            maxSpeed =
                speed;

            maxSpeedID =
                n;
        }


        maxChemical =
            max(
                maxChemical,
                fabs(
                    chemical[n]
                )
            );


        double forceMag =
            sqrt(
                forceX[n]
                *
                forceX[n]

                +

                forceY[n]
                *
                forceY[n]

                +

                forceZ[n]
                *
                forceZ[n]
            );


        maxForce =
            max(
                maxForce,
                forceMag
            );


        mass
            +=
        phi[n];
    }


    double massError =
        (
            mass
            -
            initialMass
        )
        /
        initialMass;


    size_t xy =
        static_cast<size_t>(
            NX
            *
            NY
        );


    int maxUZ =
        static_cast<int>(
            maxSpeedID
            /
            xy
        );


    size_t remain =
        maxSpeedID

        -

        static_cast<size_t>(
            maxUZ
        )
        *
        xy;


    int maxUY =
        static_cast<int>(
            remain
            /
            NX
        );


    int maxUX =
        static_cast<int>(
            remain
            %
            NX
        );


    cout
        << "\n---------------------------------------------\n";


    cout
        << "step              = "
        << step
        << "\n";


    cout
        << "T                 = "
        << step
        /
        timeScale
        << "\n";


    cout
        << "phi range         = ["
        << phiMin
        << ", "
        << phiMax
        << "]\n";


    cout
        << "pressure range    = ["
        << pMin
        << ", "
        << pMax
        << "]\n";


    cout
        << "max |chemical|    = "
        << maxChemical
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
        << "max |u| location  = ("
        << maxUX
        << ", "
        << maxUY
        << ", "
        << maxUZ
        << ")\n";


    cout
        << "relative phi mass = "
        << massError
        << "\n";


    cout
        << "---------------------------------------------\n";
}


// ============================================================
// 31. SAVE PHI ONLY
//
// Used for MATLAB visualization.
// ============================================================

void savePhi(
    int step
)
{
    char filename[256];


    sprintf(
        filename,
        "RTI3D_ACE_Re128_stable_step_%05d.bin",
        step
    );


    ofstream fout(
        filename,
        ios::binary
    );


    if(!fout)
    {
        cerr
            << "Cannot save phi file: "
            << filename
            << "\n";

        return;
    }


    fout.write(
        reinterpret_cast<const char*>(
            phi.data()
        ),

        static_cast<streamsize>(
            phi.size()
            *
            sizeof(double)
        )
    );


    fout.close();


    cout
        << "Saved phi: "
        << filename
        << "\n";
}


// ============================================================
// 32. SAVE CHECKPOINT
//
// This file stores everything required to continue
// the time integration.
//
// IMPORTANT:
//
// This is NOT the same as savePhi().
//
// savePhi:
//   visualization only.
//
// checkpoint:
//   full restart data.
// ============================================================

void saveCheckpoint(
    int step
)
{
    char filename[256];


    sprintf(
        filename,
        "RTI3D_checkpoint_step_%05d.bin",
        step
    );


    ofstream fout(
        filename,
        ios::binary
    );


    if(!fout)
    {
        cerr
            << "Cannot create checkpoint: "
            << filename
            << "\n";

        return;
    }


    // ========================================================
    // HEADER
    // ========================================================

    const char magic[8] =
    {
        'R',
        'T',
        'I',
        '3',
        'D',
        'C',
        'P',
        '1'
    };


    fout.write(
        magic,
        8
    );


    int nx =
        NX;

    int ny =
        NY;

    int nz =
        NZ;


    fout.write(
        reinterpret_cast<const char*>(
            &nx
        ),
        sizeof(int)
    );


    fout.write(
        reinterpret_cast<const char*>(
            &ny
        ),
        sizeof(int)
    );


    fout.write(
        reinterpret_cast<const char*>(
            &nz
        ),
        sizeof(int)
    );


    fout.write(
        reinterpret_cast<const char*>(
            &step
        ),
        sizeof(int)
    );


    fout.write(
        reinterpret_cast<const char*>(
            &initialMass
        ),
        sizeof(double)
    );


    // ========================================================
    // Helper
    // ========================================================

    auto writeVector =
        [&fout](
            const vector<double>& v
        )
        {
            fout.write(
                reinterpret_cast<const char*>(
                    v.data()
                ),

                static_cast<streamsize>(
                    v.size()
                    *
                    sizeof(double)
                )
            );
        };


    // ========================================================
    // DISTRIBUTIONS
    // ========================================================

    writeVector(
        phaseDist
    );


    writeVector(
        flowDist
    );


    // ========================================================
    // CURRENT MACROSCOPIC STATE
    // ========================================================

    writeVector(
        phi
    );


    writeVector(
        ux
    );


    writeVector(
        uy
    );


    writeVector(
        uz
    );


    writeVector(
        pressure
    );


    // ========================================================
    // PREVIOUS-TIME STATE
    //
    // Required by d(phi*u)/dt.
    //
    // oldQ is also retained so that checkpoint format
    // remains compatible if the Q diagnostic is used later.
    // ========================================================

    writeVector(
        oldPhiUx
    );


    writeVector(
        oldPhiUy
    );


    writeVector(
        oldPhiUz
    );


    writeVector(
        oldQ
    );


    if(!fout)
    {
        cerr
            << "Error while writing checkpoint.\n";

        fout.close();

        return;
    }


    fout.close();


    cout
        << "\n=============================================\n";

    cout
        << "Checkpoint saved\n";

    cout
        << "step = "
        << step
        << "\n";

    cout
        << "T    = "
        << step
        /
        timeScale
        << "\n";

    cout
        << "file = "
        << filename
        << "\n";

    cout
        << "=============================================\n";
}


// ============================================================
// 33. LOAD CHECKPOINT
//
// Returns the saved step number.
// ============================================================

int loadCheckpoint(
    const string& filename
)
{
    ifstream fin(
        filename,
        ios::binary
    );


    if(!fin)
    {
        cerr
            << "Cannot open checkpoint:\n"
            << filename
            << "\n";

        exit(1);
    }


    // ========================================================
    // CHECK MAGIC
    // ========================================================

    char magic[8];


    fin.read(
        magic,
        8
    );


    const char correctMagic[8] =
    {
        'R',
        'T',
        'I',
        '3',
        'D',
        'C',
        'P',
        '1'
    };


    for(int i = 0;
        i < 8;
        i++)
    {
        if(
            magic[i]
            !=
            correctMagic[i]
        )
        {
            cerr
                << "Invalid checkpoint format.\n";

            exit(1);
        }
    }


    // ========================================================
    // HEADER
    // ========================================================

    int nx;
    int ny;
    int nz;

    int savedStep;


    fin.read(
        reinterpret_cast<char*>(
            &nx
        ),
        sizeof(int)
    );


    fin.read(
        reinterpret_cast<char*>(
            &ny
        ),
        sizeof(int)
    );


    fin.read(
        reinterpret_cast<char*>(
            &nz
        ),
        sizeof(int)
    );


    fin.read(
        reinterpret_cast<char*>(
            &savedStep
        ),
        sizeof(int)
    );


    fin.read(
        reinterpret_cast<char*>(
            &initialMass
        ),
        sizeof(double)
    );


    // ========================================================
    // GRID CHECK
    // ========================================================

    if(
        nx != NX
        ||
        ny != NY
        ||
        nz != NZ
    )
    {
        cerr
            << "Checkpoint grid mismatch.\n";


        cerr
            << "Checkpoint grid = "
            << nx
            << " x "
            << ny
            << " x "
            << nz
            << "\n";


        cerr
            << "Current code grid = "
            << NX
            << " x "
            << NY
            << " x "
            << NZ
            << "\n";


        exit(1);
    }


    // ========================================================
    // Helper
    // ========================================================

    auto readVector =
        [&fin](
            vector<double>& v
        )
        {
            fin.read(
                reinterpret_cast<char*>(
                    v.data()
                ),

                static_cast<streamsize>(
                    v.size()
                    *
                    sizeof(double)
                )
            );
        };


    // ========================================================
    // DISTRIBUTIONS
    // ========================================================

    readVector(
        phaseDist
    );


    readVector(
        flowDist
    );


    // ========================================================
    // CURRENT MACROSCOPIC STATE
    // ========================================================

    readVector(
        phi
    );


    readVector(
        ux
    );


    readVector(
        uy
    );


    readVector(
        uz
    );


    readVector(
        pressure
    );


    // ========================================================
    // PREVIOUS-TIME STATE
    // ========================================================

    readVector(
        oldPhiUx
    );


    readVector(
        oldPhiUy
    );


    readVector(
        oldPhiUz
    );


    readVector(
        oldQ
    );


    if(!fin)
    {
        cerr
            << "Checkpoint is incomplete or corrupted.\n";

        exit(1);
    }


    fin.close();


    // ========================================================
    // RECONSTRUCT DERIVED FIELDS
    //
    // These arrays do not need to be stored:
    //
    // rho
    // chemical
    // grad(phi)
    // grad(rho)
    // force
    // grad(u)
    // ========================================================

    updateDensity();


    computeDerivatives();


    computeForce();


    computeVelocityGradient();


    // Temporary arrays can be cleared safely.

    fill(
        phaseNext.begin(),
        phaseNext.end(),
        0.0
    );


    fill(
        flowNext.begin(),
        flowNext.end(),
        0.0
    );


    cout
        << "\n=============================================\n";

    cout
        << "CHECKPOINT LOADED SUCCESSFULLY\n";

    cout
        << "file        = "
        << filename
        << "\n";

    cout
        << "restart step= "
        << savedStep
        << "\n";

    cout
        << "restart T   = "
        << savedStep
        /
        timeScale
        << "\n";

    cout
        << "=============================================\n";


    return savedStep;
}


// ============================================================
// MAIN
// ============================================================

int main(
    int argc,
    char* argv[]
)
{
    // ========================================================
    // COMMAND LINE
    //
    // Fresh:
    //
    // .\rti3d.exe 1600
    //
    // Restart:
    //
    // .\rti3d.exe 3200 RTI3D_checkpoint_step_01600.bin
    // ========================================================

    int maxSteps =
        stepT4;


    if(argc >= 2)
    {
        maxSteps =
            atoi(
                argv[1]
            );
    }


    bool restartMode =
        false;


    string restartFile;


    if(argc >= 3)
    {
        restartMode =
            true;


        restartFile =
            argv[2];
    }


    cout
        << scientific
        << setprecision(12);


    cout
        << "=============================================\n";

    cout
        << "Zhan 2023 - 3D RTI\n";

    cout
        << "ACE + MDF-LBM\n";

    cout
        << "D3Q7 phase + D3Q15 flow\n";

    cout
        << "Full MRT\n";

    cout
        << "Checkpoint / Restart version\n";

    cout
        << "=============================================\n";


    cout
        << "Grid            = "
        << NX
        << " x "
        << NY
        << " x "
        << NZ
        << "\n";


    cout
        << "d               = "
        << d
        << "\n";


    cout
        << "W               = "
        << W
        << "\n";


    cout
        << "At              = "
        << Atwood
        << "\n";


    cout
        << "Re              = "
        << Reynolds
        << "\n";


    cout
        << "Pe              = "
        << Peclet
        << "\n";


    cout
        << "Ug              = "
        << Ug
        << "\n";


    cout
        << "g               = "
        << gravity
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
        << "nuA=nuB         = "
        << nu
        << "\n";


    cout
        << "sigma           = "
        << sigma
        << "\n";


    cout
        << "Mphi            = "
        << mobility
        << "\n";


    cout
        << "sPhase1         = "
        << sPhase1
        << "\n";


    cout
        << "su1             = "
        << su1
        << "\n";


    cout
        << "su2             = "
        << su2
        << "\n";


    cout
        << "surface force   = "
        << (
            USE_SURFACE_FORCE
            ?
            "ON"
            :
            "OFF"
        )
        << "\n";


    cout
        << "dQ/dt in G      = "
        << (
            USE_DQDT_IN_G
            ?
            "ON"
            :
            "OFF"
        )
        << "\n";


    cout
        << "timeScale       = "
        << timeScale
        << "\n";


    cout
        << "T=1 step        = "
        << stepT1
        << "\n";


    cout
        << "T=2 step        = "
        << stepT2
        << "\n";


    cout
        << "T=3 step        = "
        << stepT3
        << "\n";


    cout
        << "T=4 step        = "
        << stepT4
        << "\n";


    cout
        << "target step     = "
        << maxSteps
        << "\n";


    cout
        << "restart mode    = "
        << (
            restartMode
            ?
            "YES"
            :
            "NO"
        )
        << "\n";


    cout
        << "=============================================\n\n";


    // ========================================================
    // MATRICES
    // ========================================================

    invertMatrix<Q7>(
        M7,
        Minv7
    );


    invertMatrix<Q15>(
        M15,
        Minv15
    );


    checkMatrices();


    // ========================================================
    // INITIALIZATION OR RESTART
    // ========================================================

    int startStep =
        0;


    if(restartMode)
    {
        // ----------------------------------------------------
        // Restart
        // ----------------------------------------------------

        startStep =
            loadCheckpoint(
                restartFile
            );


        if(maxSteps < startStep)
        {
            cerr
                << "\nTarget step "
                << maxSteps
                << " is smaller than checkpoint step "
                << startStep
                << ".\n";

            return 1;
        }


        diagnostics(
            startStep
        );
    }
    else
    {
        // ----------------------------------------------------
        // Fresh simulation
        // ----------------------------------------------------

        initializeFields();


        updateDensity();


        computeDerivatives();


        computeForce();


        initializeHydrostaticPressure();


        initializeDistributions();


        recoverPhi();


        updateDensity();


        computeDerivatives();


        computeForce();


        recoverVelocity();


        recoverPressure();


        computeVelocityGradient();


        saveOldFields();


        // ----------------------------------------------------
        // Initial mass
        // ----------------------------------------------------

        initialMass =
            0.0;


        for(size_t n = 0;
            n < N;
            n++)
        {
            initialMass
                +=
            phi[n];
        }


        diagnostics(
            0
        );


        savePhi(
            0
        );
    }


    // ========================================================
    // TIME LOOP
    //
    // Fresh:
    //
    // startStep = 0
    //
    // Restart from 1600:
    //
    // startStep = 1600
    //
    // therefore next step is 1601.
    // ========================================================

    for(int step = startStep + 1;
        step <= maxSteps;
        step++)
    {
        // ----------------------------------------------------
        // Collision + streaming
        // ----------------------------------------------------

        collideStreamPhase();


        collideStreamFlow();


        // ----------------------------------------------------
        // Save current t^n fields
        //
        // These become "old" after recovering t^(n+1).
        // ----------------------------------------------------

        saveOldFields();


        // ----------------------------------------------------
        // Recover new state
        // ----------------------------------------------------

        recoverPhi();


        updateDensity();


        computeDerivatives();


        computeForce();


        recoverVelocity();


        recoverPressure();


        computeVelocityGradient();


        // ----------------------------------------------------
        // Stability check
        // ----------------------------------------------------

        if(
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
        // Diagnostics every 100 steps
        // ----------------------------------------------------

        if(
            step % 100 == 0
            ||
            step == maxSteps
        )
        {
            diagnostics(
                step
            );
        }


        // ----------------------------------------------------
        // PAPER SNAPSHOTS
        // ----------------------------------------------------

        if(
            step == stepT1
            ||
            step == stepT2
            ||
            step == stepT3
            ||
            step == stepT4
        )
        {
            savePhi(
                step
            );
        }


        // ----------------------------------------------------
        // CHECKPOINT
        //
        // Save at paper times:
        //
        // T=1,2,3,4
        //
        // Also save at the requested final step.
        // ----------------------------------------------------

        if(
            step == stepT1
            ||
            step == stepT2
            ||
            step == stepT3
            ||
            step == stepT4
            ||
            step == maxSteps
        )
        {
            saveCheckpoint(
                step
            );
        }
    }


    cout
        << "\n=============================================\n";

    cout
        << "Simulation finished.\n";

    cout
        << "Final step = "
        << maxSteps
        << "\n";

    cout
        << "Final T    = "
        << maxSteps
        /
        timeScale
        << "\n";

    cout
        << "=============================================\n";


    return 0;
}