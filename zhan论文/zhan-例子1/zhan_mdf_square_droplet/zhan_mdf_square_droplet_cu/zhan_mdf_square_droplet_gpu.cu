#include <cuda_runtime.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <utility>

using namespace std;

// ============================================================
// Zhan 2023 square droplet: CUDA GPU version
// Full coupled CHE + MDF flow
//
// Target:
// GTX 1060, sm_61
// CUDA 12.9
// MSVC v143
//
// Precision:
// double
// ============================================================


// ============================================================
// 1. Basic parameters
// ============================================================

constexpr double Lx = 2.0;
constexpr double Ly = 2.0;

constexpr double dx = 0.01;
constexpr double W  = 0.04;

constexpr double c  = 100.0;
constexpr double dt = dx / c;

constexpr int NX = 200;
constexpr int NY = 200;

constexpr int N = NX * NY;

constexpr int Q5 = 5;
constexpr int Q9 = 9;

constexpr double cs2 =
    c * c / 3.0;


// ============================================================
// 2. Physical parameters
// ============================================================

constexpr double rhoA = 1000.0;
constexpr double rhoB = 1.0;

constexpr double phiA =  1.0;
constexpr double phiB = -1.0;

constexpr double mobility = 0.1;

constexpr double eta = 1.0;

constexpr double sigma = 0.001;

constexpr double Rfinal = 0.5;

constexpr double PI =
    3.1415926535897932384626433832795;


// ============================================================
// 3. Free-energy parameters
// ============================================================

constexpr double betaPF =
    0.01875;

constexpr double kappaPF =
    1.5e-5;


// ============================================================
// 4. Relaxation parameters
// ============================================================

constexpr double nu =
    0.1;

constexpr double sPhase1 =
    1.25;

constexpr double su1 =
    1.25;


// ============================================================
// 5. Output interval
// ============================================================

// Save phi and print diagnostics every 100000 steps
constexpr int OUTPUT_INTERVAL =
    100000;

// NaN check every 10000 steps
constexpr int CHECK_INTERVAL =
    10000;


// ============================================================
// 6. Host D2Q5 constants
// ============================================================

static const int h_cx5[Q5] =
{
     0,
     1,
     0,
    -1,
     0
};


static const int h_cy5[Q5] =
{
     0,
     0,
     1,
     0,
    -1
};


static const double h_w5[Q5] =
{
    1.0 / 3.0,

    1.0 / 6.0,
    1.0 / 6.0,
    1.0 / 6.0,
    1.0 / 6.0
};


// ============================================================
// 7. Host D2Q9 constants
// ============================================================

static const int h_cx9[Q9] =
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


static const int h_cy9[Q9] =
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


static const double h_w9[Q9] =
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
// 8. Relaxation matrices
// ============================================================

static const double h_Sphase[Q5] =
{
    1.0,
    sPhase1,
    sPhase1,
    1.0,
    1.0
};


static const double h_Sflow[Q9] =
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
// 9. D2Q5 transformation matrices
// ============================================================

static const double h_M5[Q5 * Q5] =
{
    1, 1, 1, 1, 1,

    0, 1, 0,-1, 0,

    0, 0, 1, 0,-1,

    0, 1, 0, 1, 0,

    0, 0, 1, 0, 1
};


static const double h_Minv5[Q5 * Q5] =
{
    1, 0, 0,-1,-1,

    0, 0.5, 0, 0.5, 0,

    0, 0, 0.5, 0, 0.5,

    0,-0.5, 0, 0.5, 0,

    0, 0,-0.5, 0, 0.5
};


// ============================================================
// 10. D2Q9 transformation matrices
// ============================================================

static const double h_M9[Q9 * Q9] =
{
    1,1,1,1,1,1,1,1,1,

    0,1,0,-1,0,1,-1,-1,1,

    0,0,1,0,-1,1,1,-1,-1,

    0,1,0,1,0,1,1,1,1,

    0,0,0,0,0,1,-1,1,-1,

    0,0,1,0,1,1,1,1,1,

    0,0,0,0,0,1,1,-1,-1,

    0,0,0,0,0,1,-1,-1,1,

    0,0,0,0,0,1,1,1,1
};


static const double h_Minv9[Q9 * Q9] =
{
    1,0,0,-1,0,-1,0,0,1,

    0,0.5,0,0.5,0,0,0,-0.5,-0.5,

    0,0,0.5,0,0,0.5,-0.5,0,-0.5,

    0,-0.5,0,0.5,0,0,0,0.5,-0.5,

    0,0,-0.5,0,0,0.5,0.5,0,-0.5,

    0,0,0,0,0.25,0,0.25,0.25,0.25,

    0,0,0,0,-0.25,0,0.25,-0.25,0.25,

    0,0,0,0,0.25,0,-0.25,-0.25,0.25,

    0,0,0,0,-0.25,0,-0.25,0.25,0.25
};


// ============================================================
// 11. GPU constant memory
// ============================================================

__constant__ int d_cx5[Q5];
__constant__ int d_cy5[Q5];

__constant__ double d_w5[Q5];

__constant__ double d_Sphase[Q5];

__constant__ double d_M5[Q5 * Q5];

__constant__ double d_Minv5[Q5 * Q5];


__constant__ int d_cx9[Q9];
__constant__ int d_cy9[Q9];

__constant__ double d_w9[Q9];

__constant__ double d_Sflow[Q9];

__constant__ double d_M9[Q9 * Q9];

__constant__ double d_Minv9[Q9 * Q9];


// ============================================================
// 12. CUDA error checking
// ============================================================

#define CUDA_CHECK(call)                                      \
do                                                            \
{                                                             \
    cudaError_t err__ = (call);                               \
                                                              \
    if (err__ != cudaSuccess)                                 \
    {                                                         \
        cerr                                                  \
            << "CUDA error at "                               \
            << __FILE__                                       \
            << ":"                                            \
            << __LINE__                                       \
            << " -> "                                         \
            << cudaGetErrorString(err__)                      \
            << "\n";                                          \
                                                              \
        exit(EXIT_FAILURE);                                   \
    }                                                         \
}                                                             \
while (0)


inline void checkKernel(
    const char* name
)
{
    cudaError_t err =
        cudaGetLastError();

    if (
        err
        !=
        cudaSuccess
    )
    {
        cerr
            << "Kernel launch failed ["
            << name
            << "]: "
            << cudaGetErrorString(err)
            << "\n";

        exit(EXIT_FAILURE);
    }
}


// ============================================================
// 13. GPU index helpers
// ============================================================

__device__
__forceinline__
int idx2D(
    int x,
    int y
)
{
    return
        y * NX
        +
        x;
}


// ------------------------------------------------------------
// SoA storage:
//
// distribution[q][cell]
//
// This gives more coalesced memory access on GPU.
// ------------------------------------------------------------

__device__
__forceinline__
int idxQ5(
    int q,
    int id
)
{
    return
        q * N
        +
        id;
}


__device__
__forceinline__
int idxQ9(
    int q,
    int id
)
{
    return
        q * N
        +
        id;
}


// ============================================================
// 14. Periodic BC
// ============================================================

__device__
__forceinline__
int periodicX(
    int x
)
{
    if (
        x < 0
    )
    {
        return
            x + NX;
    }

    if (
        x >= NX
    )
    {
        return
            x - NX;
    }

    return
        x;
}


__device__
__forceinline__
int periodicY(
    int y
)
{
    if (
        y < 0
    )
    {
        return
            y + NY;
    }

    if (
        y >= NY
    )
    {
        return
            y - NY;
    }

    return
        y;
}


// ============================================================
// 15. Initialize phi
// ============================================================

__global__
void kernelInitializePhi(
    double* phi,
    double squareLength
)
{
    int xIndex =
        blockIdx.x
        *
        blockDim.x
        +
        threadIdx.x;

    int yIndex =
        blockIdx.y
        *
        blockDim.y
        +
        threadIdx.y;


    if (
        xIndex >= NX
        ||
        yIndex >= NY
    )
    {
        return;
    }


    int id =
        idx2D(
            xIndex,
            yIndex
        );


    double x =
        -1.0
        +
        (
            xIndex
            +
            0.5
        )
        *
        dx;


    double y =
        -1.0
        +
        (
            yIndex
            +
            0.5
        )
        *
        dx;


    double a =
        0.5
        *
        (
            W
            -
            squareLength
        );


    double b =
        0.5
        *
        (
            squareLength
            -
            W
        );


    double dist =
        0.0;


    bool interfaceRegion =
        true;


    if (
        x <= a
        &&
        y <= a
    )
    {
        double rx =
            x - a;

        double ry =
            y - a;

        dist =
            0.5
            *
            W

            -

            sqrt(
                rx * rx
                +
                ry * ry
            );
    }

    else if (
        x >= b
        &&
        y <= a
    )
    {
        double rx =
            x - b;

        double ry =
            y - a;

        dist =
            0.5
            *
            W

            -

            sqrt(
                rx * rx
                +
                ry * ry
            );
    }

    else if (
        x <= a
        &&
        y >= b
    )
    {
        double rx =
            x - a;

        double ry =
            y - b;

        dist =
            0.5
            *
            W

            -

            sqrt(
                rx * rx
                +
                ry * ry
            );
    }

    else if (
        x >= b
        &&
        y >= b
    )
    {
        double rx =
            x - b;

        double ry =
            y - b;

        dist =
            0.5
            *
            W

            -

            sqrt(
                rx * rx
                +
                ry * ry
            );
    }

    else if (
        x > a
        &&
        x < b
        &&
        y <= a
    )
    {
        dist =
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
        dist =
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
        dist =
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
        dist =
            0.5
            *
            squareLength

            -

            x;
    }

    else
    {
        interfaceRegion =
            false;
    }


    if (
        interfaceRegion
    )
    {
        phi[id]
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
            dist
            /
            W
        );
    }
    else
    {
        phi[id] =
            phiA;
    }
}


// ============================================================
// 16. Initialize velocity / pressure = 0
// ============================================================

__global__
void kernelZeroVelocityPressure(
    double* ux,
    double* uy,
    double* pressure,
    double* dudx,
    double* dudy,
    double* dvdx,
    double* dvdy
)
{
    int id =
        blockIdx.x
        *
        blockDim.x
        +
        threadIdx.x;


    if (
        id >= N
    )
    {
        return;
    }


    ux[id] =
        0.0;

    uy[id] =
        0.0;

    pressure[id] =
        0.0;

    dudx[id] =
        0.0;

    dudy[id] =
        0.0;

    dvdx[id] =
        0.0;

    dvdy[id] =
        0.0;
}


// ============================================================
// 17. Density from phi
// ============================================================

__global__
void kernelDensityFromPhi(
    const double* phi,
    double* rho
)
{
    int id =
        blockIdx.x
        *
        blockDim.x
        +
        threadIdx.x;


    if (
        id >= N
    )
    {
        return;
    }


    rho[id]
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
        phi[id]
        -
        phiB
    )

    +

    rhoB;
}


// ============================================================
// 18. Recover phi + rho
// ============================================================

__global__
void kernelRecoverPhiDensity(
    const double* gPhase,
    double* phi,
    double* rho
)
{
    int id =
        blockIdx.x
        *
        blockDim.x
        +
        threadIdx.x;


    if (
        id >= N
    )
    {
        return;
    }


    double ph =
        0.0;


    #pragma unroll

    for (
        int q = 0;
        q < Q5;
        q++
    )
    {
        ph
            +=
        gPhase[
            idxQ5(
                q,
                id
            )
        ];
    }


    phi[id] =
        ph;


    rho[id]
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
        ph
        -
        phiB
    )

    +

    rhoB;
}


// ============================================================
// 19. Phase derivatives + chemical potential + force
// ============================================================

__global__
void kernelPhaseDerivativesAndForce(
    const double* phi,
    double* mu,
    double* gradPhiX,
    double* gradPhiY,
    double* forceX,
    double* forceY
)
{
    int x =
        blockIdx.x
        *
        blockDim.x
        +
        threadIdx.x;

    int y =
        blockIdx.y
        *
        blockDim.y
        +
        threadIdx.y;


    if (
        x >= NX
        ||
        y >= NY
    )
    {
        return;
    }


    int id =
        idx2D(
            x,
            y
        );


    double gx =
        0.0;

    double gy =
        0.0;

    double lap =
        0.0;


    double ph0 =
        phi[id];


    #pragma unroll

    for (
        int q = 1;
        q < Q5;
        q++
    )
    {
        int xn =
            periodicX(
                x
                +
                d_cx5[q]
            );

        int yn =
            periodicY(
                y
                +
                d_cy5[q]
            );


        int nid =
            idx2D(
                xn,
                yn
            );


        double cix =
            c
            *
            d_cx5[q];

        double ciy =
            c
            *
            d_cy5[q];


        double phn =
            phi[nid];


        gx
            +=
        d_w5[q]
        *
        cix
        *
        phn;


        gy
            +=
        d_w5[q]
        *
        ciy
        *
        phn;


        lap
            +=
        2.0
        *
        d_w5[q]
        *
        (
            phn
            -
            ph0
        );
    }


    gx
        /=
    (
        cs2
        *
        dt
    );


    gy
        /=
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


    double chem =
        4.0
        *
        betaPF
        *
        (
            ph0
            -
            phiA
        )
        *
        (
            ph0
            -
            phiB
        )
        *
        (
            ph0
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


    gradPhiX[id] =
        gx;

    gradPhiY[id] =
        gy;

    mu[id] =
        chem;


    forceX[id]
        =
    chem
    *
    gx;


    forceY[id]
        =
    chem
    *
    gy;
}


// ============================================================
// 20. grad rho
// ============================================================

__global__
void kernelGradRho(
    const double* rho,
    double* gradRhoX,
    double* gradRhoY
)
{
    int x =
        blockIdx.x
        *
        blockDim.x
        +
        threadIdx.x;

    int y =
        blockIdx.y
        *
        blockDim.y
        +
        threadIdx.y;


    if (
        x >= NX
        ||
        y >= NY
    )
    {
        return;
    }


    int id =
        idx2D(
            x,
            y
        );


    double gx =
        0.0;

    double gy =
        0.0;


    #pragma unroll

    for (
        int q = 1;
        q < Q9;
        q++
    )
    {
        int xn =
            periodicX(
                x
                +
                d_cx9[q]
            );

        int yn =
            periodicY(
                y
                +
                d_cy9[q]
            );


        int nid =
            idx2D(
                xn,
                yn
            );


        double cix =
            c
            *
            d_cx9[q];

        double ciy =
            c
            *
            d_cy9[q];


        gx
            +=
        d_w9[q]
        *
        cix
        *
        rho[nid];


        gy
            +=
        d_w9[q]
        *
        ciy
        *
        rho[nid];
    }


    gradRhoX[id]
        =
    gx
    /
    (
        cs2
        *
        dt
    );


    gradRhoY[id]
        =
    gy
    /
    (
        cs2
        *
        dt
    );
}


// ============================================================
// 21. Recover velocity
// ============================================================

__global__
void kernelRecoverVelocity(
    const double* g1,
    const double* g2,
    const double* rho,
    const double* forceX,
    const double* forceY,
    double* ux,
    double* uy
)
{
    int id =
        blockIdx.x
        *
        blockDim.x
        +
        threadIdx.x;


    if (
        id >= N
    )
    {
        return;
    }


    double sum1 =
        0.0;

    double sum2 =
        0.0;


    #pragma unroll

    for (
        int q = 0;
        q < Q9;
        q++
    )
    {
        sum1
            +=
        g1[
            idxQ9(
                q,
                id
            )
        ];


        sum2
            +=
        g2[
            idxQ9(
                q,
                id
            )
        ];
    }


    double r =
        rho[id];


    ux[id]
        =
    (
        sum1

        +

        0.5
        *
        dt
        *
        forceX[id]
    )
    /
    r;


    uy[id]
        =
    (
        sum2

        +

        0.5
        *
        dt
        *
        forceY[id]
    )
    /
    r;
}


// ============================================================
// 22. Velocity gradients
// ============================================================

__global__
void kernelVelocityGradients(
    const double* ux,
    const double* uy,
    double* dudx,
    double* dudy,
    double* dvdx,
    double* dvdy
)
{
    int x =
        blockIdx.x
        *
        blockDim.x
        +
        threadIdx.x;

    int y =
        blockIdx.y
        *
        blockDim.y
        +
        threadIdx.y;


    if (
        x >= NX
        ||
        y >= NY
    )
    {
        return;
    }


    int id =
        idx2D(
            x,
            y
        );


    double ux_x =
        0.0;

    double ux_y =
        0.0;

    double uy_x =
        0.0;

    double uy_y =
        0.0;


    #pragma unroll

    for (
        int q = 1;
        q < Q9;
        q++
    )
    {
        int xn =
            periodicX(
                x
                +
                d_cx9[q]
            );

        int yn =
            periodicY(
                y
                +
                d_cy9[q]
            );


        int nid =
            idx2D(
                xn,
                yn
            );


        double cix =
            c
            *
            d_cx9[q];

        double ciy =
            c
            *
            d_cy9[q];


        ux_x
            +=
        d_w9[q]
        *
        cix
        *
        ux[nid];


        ux_y
            +=
        d_w9[q]
        *
        ciy
        *
        ux[nid];


        uy_x
            +=
        d_w9[q]
        *
        cix
        *
        uy[nid];


        uy_y
            +=
        d_w9[q]
        *
        ciy
        *
        uy[nid];
    }


    double denom =
        cs2
        *
        dt;


    dudx[id] =
        ux_x
        /
        denom;


    dudy[id] =
        ux_y
        /
        denom;


    dvdx[id] =
        uy_x
        /
        denom;


    dvdy[id] =
        uy_y
        /
        denom;
}


// ============================================================
// 23. Recover pressure
// ============================================================

__global__
void kernelRecoverPressure(
    const double* g1,
    const double* g2,
    const double* rho,
    const double* ux,
    const double* uy,
    const double* gradRhoX,
    const double* gradRhoY,
    double* pressure
)
{
    int id =
        blockIdx.x
        *
        blockDim.x
        +
        threadIdx.x;


    if (
        id >= N
    )
    {
        return;
    }


    double firstMomentX =
        0.0;

    double firstMomentY =
        0.0;


    #pragma unroll

    for (
        int q = 0;
        q < Q9;
        q++
    )
    {
        double cix =
            c
            *
            d_cx9[q];

        double ciy =
            c
            *
            d_cy9[q];


        firstMomentX
            +=
        cix
        *
        g1[
            idxQ9(
                q,
                id
            )
        ];


        firstMomentY
            +=
        ciy
        *
        g2[
            idxQ9(
                q,
                id
            )
        ];
    }


    double u =
        ux[id];

    double v =
        uy[id];


    pressure[id]
        =
    0.5
    *
    (
        firstMomentX

        +

        firstMomentY

        -

        rho[id]
        *
        (
            u * u
            +
            v * v
        )

        +

        0.5
        *
        dt
        *
        cs2
        *
        (
            u
            *
            gradRhoX[id]

            +

            v
            *
            gradRhoY[id]
        )
    );
}


// ============================================================
// 24. Initialize phase distribution
// ============================================================

__global__
void kernelInitializePhaseDistribution(
    const double* phi,
    const double* mu,
    const double* ux,
    const double* uy,
    double* gPhase
)
{
    int id =
        blockIdx.x
        *
        blockDim.x
        +
        threadIdx.x;


    if (
        id >= N
    )
    {
        return;
    }


    double feq[Q5];

    double meq[Q5];


    double ph =
        phi[id];

    double chem =
        mu[id];

    double u =
        ux[id];

    double v =
        uy[id];


    #pragma unroll

    for (
        int q = 0;
        q < Q5;
        q++
    )
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
                d_w5[q]
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
            double cix =
                c
                *
                d_cx5[q];

            double ciy =
                c
                *
                d_cy5[q];


            feq[q]
                =
            d_w5[q]
            *
            eta
            *
            chem

            +

            d_w5[q]
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


    #pragma unroll

    for (
        int k = 0;
        k < Q5;
        k++
    )
    {
        double sum =
            0.0;


        #pragma unroll

        for (
            int q = 0;
            q < Q5;
            q++
        )
        {
            sum
                +=
            d_M5[
                k * Q5
                +
                q
            ]
            *
            feq[q];
        }


        meq[k] =
            sum;
    }


    #pragma unroll

    for (
        int q = 0;
        q < Q5;
        q++
    )
    {
        double sum =
            0.0;


        #pragma unroll

        for (
            int k = 0;
            k < Q5;
            k++
        )
        {
            sum
                +=
            d_Minv5[
                q * Q5
                +
                k
            ]
            *
            meq[k];
        }


        gPhase[
            idxQ5(
                q,
                id
            )
        ]
            =
        sum;
    }
}


// ============================================================
// 25. Initialize flow distribution
// ============================================================

__global__
void kernelInitializeFlowDistribution(
    const double* rho,
    const double* ux,
    const double* uy,
    const double* pressure,
    const double* forceX,
    const double* forceY,
    double* g1,
    double* g2
)
{
    int id =
        blockIdx.x
        *
        blockDim.x
        +
        threadIdx.x;


    if (
        id >= N
    )
    {
        return;
    }


    double r =
        rho[id];

    double u =
        ux[id];

    double v =
        uy[id];

    double p =
        pressure[id];


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


    double Qyy =
        r
        *
        v
        *
        v

        +

        p;


    double feq1[Q9];

    double feq2[Q9];

    double Fi1[Q9];

    double Fi2[Q9];


    double meq1[Q9];

    double meq2[Q9];

    double mF1[Q9];

    double mF2[Q9];


    double mg1[Q9];

    double mg2[Q9];


    #pragma unroll

    for (
        int q = 0;
        q < Q9;
        q++
    )
    {
        double cix =
            c
            *
            d_cx9[q];

        double ciy =
            c
            *
            d_cy9[q];


        feq1[q]
            =
        d_w9[q]
        *
        (
            r
            *
            u

            +

            (
                cix
                *
                Qxx

                +

                ciy
                *
                Qxy
            )
            /
            cs2
        );


        feq2[q]
            =
        d_w9[q]
        *
        (
            r
            *
            v

            +

            (
                cix
                *
                Qxy

                +

                ciy
                *
                Qyy
            )
            /
            cs2
        );


        Fi1[q]
            =
        d_w9[q]
        *
        forceX[id];


        Fi2[q]
            =
        d_w9[q]
        *
        forceY[id];
    }


    #pragma unroll

    for (
        int k = 0;
        k < Q9;
        k++
    )
    {
        double a1 =
            0.0;

        double a2 =
            0.0;

        double b1 =
            0.0;

        double b2 =
            0.0;


        #pragma unroll

        for (
            int q = 0;
            q < Q9;
            q++
        )
        {
            double M =
                d_M9[
                    k * Q9
                    +
                    q
                ];


            a1
                +=
            M
            *
            feq1[q];


            a2
                +=
            M
            *
            feq2[q];


            b1
                +=
            M
            *
            Fi1[q];


            b2
                +=
            M
            *
            Fi2[q];
        }


        meq1[k] =
            a1;

        meq2[k] =
            a2;

        mF1[k] =
            b1;

        mF2[k] =
            b2;


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


    #pragma unroll

    for (
        int q = 0;
        q < Q9;
        q++
    )
    {
        double a1 =
            0.0;

        double a2 =
            0.0;


        #pragma unroll

        for (
            int k = 0;
            k < Q9;
            k++
        )
        {
            double Mi =
                d_Minv9[
                    q * Q9
                    +
                    k
                ];


            a1
                +=
            Mi
            *
            mg1[k];


            a2
                +=
            Mi
            *
            mg2[k];
        }


        g1[
            idxQ9(
                q,
                id
            )
        ]
            =
        a1;


        g2[
            idxQ9(
                q,
                id
            )
        ]
            =
        a2;
    }
}


// ============================================================
// 26. Initialize old-time fields
// ============================================================

__global__
void kernelInitializeOldFields(
    const double* phi,
    const double* rho,
    const double* ux,
    const double* uy,
    const double* pressure,
    double* oldPhiUx,
    double* oldPhiUy,
    double* oldQxx,
    double* oldQxy,
    double* oldQyy
)
{
    int id =
        blockIdx.x
        *
        blockDim.x
        +
        threadIdx.x;


    if (
        id >= N
    )
    {
        return;
    }


    double r =
        rho[id];

    double u =
        ux[id];

    double v =
        uy[id];

    double p =
        pressure[id];


    oldPhiUx[id]
        =
    phi[id]
    *
    u;


    oldPhiUy[id]
        =
    phi[id]
    *
    v;


    oldQxx[id]
        =
    r
    *
    u
    *
    u

    +

    p;


    oldQxy[id]
        =
    r
    *
    u
    *
    v;


    oldQyy[id]
        =
    r
    *
    v
    *
    v

    +

    p;
}


// ============================================================
// 27. Phase collision
// ============================================================

__global__
void kernelCollidePhase(
    const double* gPhase,
    double* gPhasePost,
    const double* phi,
    const double* mu,
    const double* ux,
    const double* uy,
    const double* oldPhiUx,
    const double* oldPhiUy
)
{
    int id =
        blockIdx.x
        *
        blockDim.x
        +
        threadIdx.x;


    if (
        id >= N
    )
    {
        return;
    }


    double m[Q5];

    double feq[Q5];

    double meq[Q5];

    double Gi[Q5];

    double mG[Q5];

    double mPost[Q5];


    // --------------------------------------------------------
    // m = M g
    // --------------------------------------------------------

    #pragma unroll

    for (
        int k = 0;
        k < Q5;
        k++
    )
    {
        double sum =
            0.0;


        #pragma unroll

        for (
            int q = 0;
            q < Q5;
            q++
        )
        {
            sum
                +=
            d_M5[
                k * Q5
                +
                q
            ]
            *
            gPhase[
                idxQ5(
                    q,
                    id
                )
            ];
        }


        m[k] =
            sum;
    }


    // --------------------------------------------------------
    // equilibrium
    // --------------------------------------------------------

    double ph =
        phi[id];

    double chem =
        mu[id];

    double u =
        ux[id];

    double v =
        uy[id];


    #pragma unroll

    for (
        int q = 0;
        q < Q5;
        q++
    )
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
                d_w5[q]
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
            double cix =
                c
                *
                d_cx5[q];

            double ciy =
                c
                *
                d_cy5[q];


            feq[q]
                =
            d_w5[q]
            *
            eta
            *
            chem

            +

            d_w5[q]
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


    #pragma unroll

    for (
        int k = 0;
        k < Q5;
        k++
    )
    {
        double sum =
            0.0;


        #pragma unroll

        for (
            int q = 0;
            q < Q5;
            q++
        )
        {
            sum
                +=
            d_M5[
                k * Q5
                +
                q
            ]
            *
            feq[q];
        }


        meq[k] =
            sum;
    }


    // --------------------------------------------------------
    // Phase G
    // --------------------------------------------------------

    double dPhiUxDt
        =
    (
        ph
        *
        u

        -

        oldPhiUx[id]
    )
    /
    dt;


    double dPhiUyDt
        =
    (
        ph
        *
        v

        -

        oldPhiUy[id]
    )
    /
    dt;


    #pragma unroll

    for (
        int q = 0;
        q < Q5;
        q++
    )
    {
        double cix =
            c
            *
            d_cx5[q];

        double ciy =
            c
            *
            d_cy5[q];


        Gi[q]
            =
        d_w5[q]
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


    #pragma unroll

    for (
        int k = 0;
        k < Q5;
        k++
    )
    {
        double sum =
            0.0;


        #pragma unroll

        for (
            int q = 0;
            q < Q5;
            q++
        )
        {
            sum
                +=
            d_M5[
                k * Q5
                +
                q
            ]
            *
            Gi[q];
        }


        mG[k] =
            sum;


        double sk =
            d_Sphase[k];


        mPost[k]
            =
        (
            1.0
            -
            sk
        )
        *
        m[k]

        +

        sk
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
            sk
        )
        *
        mG[k];
    }


    #pragma unroll

    for (
        int q = 0;
        q < Q5;
        q++
    )
    {
        double sum =
            0.0;


        #pragma unroll

        for (
            int k = 0;
            k < Q5;
            k++
        )
        {
            sum
                +=
            d_Minv5[
                q * Q5
                +
                k
            ]
            *
            mPost[k];
        }


        gPhasePost[
            idxQ5(
                q,
                id
            )
        ]
            =
        sum;
    }
}


// ============================================================
// 28. Flow collision + save old fields
// ============================================================

__global__
void kernelCollideFlowAndSaveOld(
    const double* g1,
    const double* g2,
    double* g1Post,
    double* g2Post,

    const double* phi,
    const double* rho,
    const double* ux,
    const double* uy,
    const double* pressure,

    const double* forceX,
    const double* forceY,

    const double* gradRhoX,
    const double* gradRhoY,

    const double* dudx,
    const double* dudy,
    const double* dvdx,
    const double* dvdy,

    double* oldPhiUx,
    double* oldPhiUy,

    double* oldQxx,
    double* oldQxy,
    double* oldQyy
)
{
    int id =
        blockIdx.x
        *
        blockDim.x
        +
        threadIdx.x;


    if (
        id >= N
    )
    {
        return;
    }


    double m1[Q9];

    double m2[Q9];


    double feq1[Q9];

    double feq2[Q9];


    double Fi1[Q9];

    double Fi2[Q9];


    double Gi1[Q9];

    double Gi2[Q9];


    double meq1[Q9];

    double meq2[Q9];


    double mF1[Q9];

    double mF2[Q9];


    double mG1[Q9];

    double mG2[Q9];


    double mPost1[Q9];

    double mPost2[Q9];


    // --------------------------------------------------------
    // m = M g
    // --------------------------------------------------------

    #pragma unroll

    for (
        int k = 0;
        k < Q9;
        k++
    )
    {
        double a1 =
            0.0;

        double a2 =
            0.0;


        #pragma unroll

        for (
            int q = 0;
            q < Q9;
            q++
        )
        {
            double M =
                d_M9[
                    k * Q9
                    +
                    q
                ];


            a1
                +=
            M
            *
            g1[
                idxQ9(
                    q,
                    id
                )
            ];


            a2
                +=
            M
            *
            g2[
                idxQ9(
                    q,
                    id
                )
            ];
        }


        m1[k] =
            a1;

        m2[k] =
            a2;
    }


    // --------------------------------------------------------
    // Macroscopic variables
    // --------------------------------------------------------

    double r =
        rho[id];

    double u =
        ux[id];

    double v =
        uy[id];

    double p =
        pressure[id];


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


    double Qyy =
        r
        *
        v
        *
        v

        +

        p;


    // --------------------------------------------------------
    // d_t Q
    // --------------------------------------------------------

    double dQxxDt
        =
    (
        Qxx
        -
        oldQxx[id]
    )
    /
    dt;


    double dQxyDt
        =
    (
        Qxy
        -
        oldQxy[id]
    )
    /
    dt;


    double dQyyDt
        =
    (
        Qyy
        -
        oldQyy[id]
    )
    /
    dt;


    // --------------------------------------------------------
    // Eq.(37)
    // --------------------------------------------------------

    double A1x
        =
    u
    *
    gradRhoX[id]

    -

    r
    *
    dudx[id];


    double A1y
        =
    u
    *
    gradRhoY[id]

    -

    r
    *
    dvdx[id];


    double A2x
        =
    v
    *
    gradRhoX[id]

    -

    r
    *
    dudy[id];


    double A2y
        =
    v
    *
    gradRhoY[id]

    -

    r
    *
    dvdy[id];


    // --------------------------------------------------------
    // velocity-space equilibrium / F / G
    // --------------------------------------------------------

    #pragma unroll

    for (
        int q = 0;
        q < Q9;
        q++
    )
    {
        double cix =
            c
            *
            d_cx9[q];

        double ciy =
            c
            *
            d_cy9[q];


        feq1[q]
            =
        d_w9[q]
        *
        (
            r
            *
            u

            +

            (
                cix
                *
                Qxx

                +

                ciy
                *
                Qxy
            )
            /
            cs2
        );


        feq2[q]
            =
        d_w9[q]
        *
        (
            r
            *
            v

            +

            (
                cix
                *
                Qxy

                +

                ciy
                *
                Qyy
            )
            /
            cs2
        );


        Fi1[q]
            =
        d_w9[q]
        *
        forceX[id];


        Fi2[q]
            =
        d_w9[q]
        *
        forceY[id];


        Gi1[q]
            =
        d_w9[q]
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
        d_w9[q]
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


    // --------------------------------------------------------
    // M-transform + collision
    // --------------------------------------------------------

    #pragma unroll

    for (
        int k = 0;
        k < Q9;
        k++
    )
    {
        double e1 =
            0.0;

        double e2 =
            0.0;

        double f1 =
            0.0;

        double f2 =
            0.0;

        double gg1 =
            0.0;

        double gg2 =
            0.0;


        #pragma unroll

        for (
            int q = 0;
            q < Q9;
            q++
        )
        {
            double M =
                d_M9[
                    k * Q9
                    +
                    q
                ];


            e1
                +=
            M
            *
            feq1[q];


            e2
                +=
            M
            *
            feq2[q];


            f1
                +=
            M
            *
            Fi1[q];


            f2
                +=
            M
            *
            Fi2[q];


            gg1
                +=
            M
            *
            Gi1[q];


            gg2
                +=
            M
            *
            Gi2[q];
        }


        meq1[k] =
            e1;

        meq2[k] =
            e2;


        mF1[k] =
            f1;

        mF2[k] =
            f2;


        mG1[k] =
            gg1;

        mG2[k] =
            gg2;


        double sk =
            d_Sflow[k];


        mPost1[k]
            =
        (
            1.0
            -
            sk
        )
        *
        m1[k]

        +

        sk
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
            sk
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
            sk
        )
        *
        m2[k]

        +

        sk
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
            sk
        )
        *
        (
            mF2[k]
            +
            mG2[k]
        );
    }


    // --------------------------------------------------------
    // inverse transform
    // --------------------------------------------------------

    #pragma unroll

    for (
        int q = 0;
        q < Q9;
        q++
    )
    {
        double a1 =
            0.0;

        double a2 =
            0.0;


        #pragma unroll

        for (
            int k = 0;
            k < Q9;
            k++
        )
        {
            double Mi =
                d_Minv9[
                    q * Q9
                    +
                    k
                ];


            a1
                +=
            Mi
            *
            mPost1[k];


            a2
                +=
            Mi
            *
            mPost2[k];
        }


        g1Post[
            idxQ9(
                q,
                id
            )
        ]
            =
        a1;


        g2Post[
            idxQ9(
                q,
                id
            )
        ]
            =
        a2;
    }


    // --------------------------------------------------------
    // Save current n-level values
    // Same position as CPU saveCurrentOldFields()
    // --------------------------------------------------------

    oldPhiUx[id]
        =
    phi[id]
    *
    u;


    oldPhiUy[id]
        =
    phi[id]
    *
    v;


    oldQxx[id] =
        Qxx;

    oldQxy[id] =
        Qxy;

    oldQyy[id] =
        Qyy;
}


// ============================================================
// 29. Pull streaming
//
// GPU uses pull streaming:
//
// current cell reads from source neighbour.
//
// This avoids write conflicts.
// ============================================================

__global__
void kernelStreamAll(
    const double* gPhasePost,
    const double* g1Post,
    const double* g2Post,

    double* gPhaseNext,
    double* g1Next,
    double* g2Next
)
{
    int x =
        blockIdx.x
        *
        blockDim.x
        +
        threadIdx.x;

    int y =
        blockIdx.y
        *
        blockDim.y
        +
        threadIdx.y;


    if (
        x >= NX
        ||
        y >= NY
    )
    {
        return;
    }


    int id =
        idx2D(
            x,
            y
        );


    #pragma unroll

    for (
        int q = 0;
        q < Q5;
        q++
    )
    {
        int xs =
            periodicX(
                x
                -
                d_cx5[q]
            );

        int ys =
            periodicY(
                y
                -
                d_cy5[q]
            );


        int sid =
            idx2D(
                xs,
                ys
            );


        gPhaseNext[
            idxQ5(
                q,
                id
            )
        ]
            =
        gPhasePost[
            idxQ5(
                q,
                sid
            )
        ];
    }


    #pragma unroll

    for (
        int q = 0;
        q < Q9;
        q++
    )
    {
        int xs =
            periodicX(
                x
                -
                d_cx9[q]
            );

        int ys =
            periodicY(
                y
                -
                d_cy9[q]
            );


        int sid =
            idx2D(
                xs,
                ys
            );


        g1Next[
            idxQ9(
                q,
                id
            )
        ]
            =
        g1Post[
            idxQ9(
                q,
                sid
            )
        ];


        g2Next[
            idxQ9(
                q,
                id
            )
        ]
            =
        g2Post[
            idxQ9(
                q,
                sid
            )
        ];
    }
}


// ============================================================
// 30. NaN / Inf check
// ============================================================

__global__
void kernelCheckFinite(
    const double* phi,
    const double* rho,
    const double* ux,
    const double* uy,
    const double* pressure,
    const double* mu,

    int* badFlag,
    int* badIndex
)
{
    int id =
        blockIdx.x
        *
        blockDim.x
        +
        threadIdx.x;


    if (
        id >= N
    )
    {
        return;
    }


    bool bad =
        !isfinite(
            phi[id]
        )

        ||

        !isfinite(
            rho[id]
        )

        ||

        !isfinite(
            ux[id]
        )

        ||

        !isfinite(
            uy[id]
        )

        ||

        !isfinite(
            pressure[id]
        )

        ||

        !isfinite(
            mu[id]
        );


    if (
        bad
    )
    {
        if (
            atomicCAS(
                badFlag,
                0,
                1
            )
            ==
            0
        )
        {
            *badIndex =
                id;
        }
    }
}


// ============================================================
// 31. Copy constants to GPU
// ============================================================

void copyConstantsToGPU()
{
    CUDA_CHECK(
        cudaMemcpyToSymbol(
            d_cx5,
            h_cx5,
            sizeof(h_cx5)
        )
    );


    CUDA_CHECK(
        cudaMemcpyToSymbol(
            d_cy5,
            h_cy5,
            sizeof(h_cy5)
        )
    );


    CUDA_CHECK(
        cudaMemcpyToSymbol(
            d_w5,
            h_w5,
            sizeof(h_w5)
        )
    );


    CUDA_CHECK(
        cudaMemcpyToSymbol(
            d_Sphase,
            h_Sphase,
            sizeof(h_Sphase)
        )
    );


    CUDA_CHECK(
        cudaMemcpyToSymbol(
            d_M5,
            h_M5,
            sizeof(h_M5)
        )
    );


    CUDA_CHECK(
        cudaMemcpyToSymbol(
            d_Minv5,
            h_Minv5,
            sizeof(h_Minv5)
        )
    );


    CUDA_CHECK(
        cudaMemcpyToSymbol(
            d_cx9,
            h_cx9,
            sizeof(h_cx9)
        )
    );


    CUDA_CHECK(
        cudaMemcpyToSymbol(
            d_cy9,
            h_cy9,
            sizeof(h_cy9)
        )
    );


    CUDA_CHECK(
        cudaMemcpyToSymbol(
            d_w9,
            h_w9,
            sizeof(h_w9)
        )
    );


    CUDA_CHECK(
        cudaMemcpyToSymbol(
            d_Sflow,
            h_Sflow,
            sizeof(h_Sflow)
        )
    );


    CUDA_CHECK(
        cudaMemcpyToSymbol(
            d_M9,
            h_M9,
            sizeof(h_M9)
        )
    );


    CUDA_CHECK(
        cudaMemcpyToSymbol(
            d_Minv9,
            h_Minv9,
            sizeof(h_Minv9)
        )
    );
}


// ============================================================
// 32. Check M * Minv
// ============================================================

void checkMatricesHost()
{
    double maxError5 =
        0.0;

    double maxError9 =
        0.0;


    for (
        int i = 0;
        i < Q5;
        i++
    )
    {
        for (
            int j = 0;
            j < Q5;
            j++
        )
        {
            double value =
                0.0;


            for (
                int k = 0;
                k < Q5;
                k++
            )
            {
                value
                    +=
                h_M5[
                    i * Q5
                    +
                    k
                ]
                *
                h_Minv5[
                    k * Q5
                    +
                    j
                ];
            }


            double exact =
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


    for (
        int i = 0;
        i < Q9;
        i++
    )
    {
        for (
            int j = 0;
            j < Q9;
            j++
        )
        {
            double value =
                0.0;


            for (
                int k = 0;
                k < Q9;
                k++
            )
            {
                value
                    +=
                h_M9[
                    i * Q9
                    +
                    k
                ]
                *
                h_Minv9[
                    k * Q9
                    +
                    j
                ];
            }


            double exact =
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
// 33. Save phi
// ============================================================

void savePhi(
    const double* h_phi,
    int step
)
{
    char fileName[100];


    sprintf(
        fileName,
        "phi_%05d.dat",
        step
    );


    FILE* fp =
        fopen(
            fileName,
            "w"
        );


    if (
        fp == nullptr
    )
    {
        cerr
            << "Cannot open "
            << fileName
            << "\n";

        return;
    }


    for (
        int y = 0;
        y < NY;
        y++
    )
    {
        for (
            int x = 0;
            x < NX;
            x++
        )
        {
            fprintf(
                fp,
                "%.15e",
                h_phi[
                    y * NX
                    +
                    x
                ]
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


// ============================================================
// 34. Initial conserved quantities
// ============================================================

struct InitialConserved
{
    double phaseMass =
        0.0;

    double densityMass =
        0.0;

    double momentumX =
        0.0;

    double momentumY =
        0.0;
};


// ============================================================
// 35. Host diagnostics
// ============================================================

void printDiagnosticsHost(
    int step,

    const double* phi,
    const double* rho,
    const double* mu,

    const double* forceX,
    const double* forceY,

    const double* ux,
    const double* uy,

    const double* pressure,

    const double* dudx,
    const double* dvdy,

    const InitialConserved& init
)
{
    double minPhi =
        phi[0];

    double maxPhi =
        phi[0];

    double minRho =
        rho[0];

    double maxRho =
        rho[0];


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


    for (
        int id = 0;
        id < N;
        id++
    )
    {
        minPhi =
            min(
                minPhi,
                phi[id]
            );


        maxPhi =
            max(
                maxPhi,
                phi[id]
            );


        minRho =
            min(
                minRho,
                rho[id]
            );


        maxRho =
            max(
                maxRho,
                rho[id]
            );


        maxMu =
            max(
                maxMu,
                fabs(
                    mu[id]
                )
            );


        double Fmag =
            sqrt(
                forceX[id]
                *
                forceX[id]

                +

                forceY[id]
                *
                forceY[id]
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
                    ux[id]
                )
            );


        maxUy =
            max(
                maxUy,
                fabs(
                    uy[id]
                )
            );


        double speed =
            sqrt(
                ux[id]
                *
                ux[id]

                +

                uy[id]
                *
                uy[id]
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
                    pressure[id]
                )
            );


        phaseMass
            +=
        phi[id];


        densityMass
            +=
        rho[id];


        momentumX
            +=
        rho[id]
        *
        ux[id];


        momentumY
            +=
        rho[id]
        *
        uy[id];


        double divU =
            dudx[id]
            +
            dvdy[id];


        divSquareSum
            +=
        divU
        *
        divU;
    }


    double relPhase
        =
    (
        phaseMass
        -
        init.phaseMass
    )
    /
    (
        fabs(
            init.phaseMass
        )
        +
        1.0e-30
    );


    double relDensity
        =
    (
        densityMass
        -
        init.densityMass
    )
    /
    (
        fabs(
            init.densityMass
        )
        +
        1.0e-30
    );


    double dMomX
        =
    (
        momentumX
        -
        init.momentumX
    )
    /
    N;


    double dMomY
        =
    (
        momentumY
        -
        init.momentumY
    )
    /
    N;


    double Drms =
        sqrt(
            divSquareSum
            /
            N
        );


    cout
        << "\n--------------------------------------------\n";


    cout
        << "Step = "
        << step
        << "   physical t = "
        << step * dt
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
        << relPhase
        << "\n";


    cout
        << "relative density mass  = "
        << relDensity
        << "\n";


    cout
        << "<Delta rho*u>          = "
        << dMomX
        << "\n";


    cout
        << "<Delta rho*v>          = "
        << dMomY
        << "\n";


    cout
        << "D_rms                  = "
        << Drms
        << "\n";
}


// ============================================================
// 36. Copy diagnostics GPU -> CPU
// ============================================================

void copyDiagnosticsToHost(
    double* h_phi,
    double* h_rho,
    double* h_mu,

    double* h_forceX,
    double* h_forceY,

    double* h_ux,
    double* h_uy,

    double* h_pressure,

    double* h_dudx,
    double* h_dvdy,

    const double* d_phi,
    const double* d_rho,
    const double* d_mu,

    const double* d_forceX,
    const double* d_forceY,

    const double* d_ux,
    const double* d_uy,

    const double* d_pressure,

    const double* d_dudx,
    const double* d_dvdy
)
{
    size_t bytes =
        N
        *
        sizeof(double);


    CUDA_CHECK(
        cudaMemcpy(
            h_phi,
            d_phi,
            bytes,
            cudaMemcpyDeviceToHost
        )
    );


    CUDA_CHECK(
        cudaMemcpy(
            h_rho,
            d_rho,
            bytes,
            cudaMemcpyDeviceToHost
        )
    );


    CUDA_CHECK(
        cudaMemcpy(
            h_mu,
            d_mu,
            bytes,
            cudaMemcpyDeviceToHost
        )
    );


    CUDA_CHECK(
        cudaMemcpy(
            h_forceX,
            d_forceX,
            bytes,
            cudaMemcpyDeviceToHost
        )
    );


    CUDA_CHECK(
        cudaMemcpy(
            h_forceY,
            d_forceY,
            bytes,
            cudaMemcpyDeviceToHost
        )
    );


    CUDA_CHECK(
        cudaMemcpy(
            h_ux,
            d_ux,
            bytes,
            cudaMemcpyDeviceToHost
        )
    );


    CUDA_CHECK(
        cudaMemcpy(
            h_uy,
            d_uy,
            bytes,
            cudaMemcpyDeviceToHost
        )
    );


    CUDA_CHECK(
        cudaMemcpy(
            h_pressure,
            d_pressure,
            bytes,
            cudaMemcpyDeviceToHost
        )
    );


    CUDA_CHECK(
        cudaMemcpy(
            h_dudx,
            d_dudx,
            bytes,
            cudaMemcpyDeviceToHost
        )
    );


    CUDA_CHECK(
        cudaMemcpy(
            h_dvdy,
            d_dvdy,
            bytes,
            cudaMemcpyDeviceToHost
        )
    );
}


// ============================================================
// 37. main
// ============================================================

int main(
    int argc,
    char** argv
)
{
    cout
        << scientific
        << setprecision(12);


    // --------------------------------------------------------
    // Default:
    // 1000 steps for first GPU validation.
    //
    // Examples:
    //
    // zhan_gpu.exe 1000
    //
    // zhan_gpu.exe 500000
    //
    // zhan_gpu.exe 2000000
    // --------------------------------------------------------

    int nSteps =
        1000;


    if (
        argc >= 2
    )
    {
        nSteps =
            atoi(
                argv[1]
            );


        if (
            nSteps <= 0
        )
        {
            cerr
                << "nSteps must be positive.\n";

            return
                1;
        }
    }


    cout
        << "====================================================\n";

    cout
        << "Zhan 2023 square droplet - CUDA GPU\n";

    cout
        << "FULL COUPLED CHE + MDF FLOW\n";

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
        << "\n";

    cout
        << "nSteps         = "
        << nSteps
        << "\n\n";


    // ========================================================
    // GPU information
    // ========================================================

    int device =
        0;


    cudaDeviceProp prop{};


    CUDA_CHECK(
        cudaGetDeviceProperties(
            &prop,
            device
        )
    );


    cout
        << "GPU            = "
        << prop.name
        << "\n";


    cout
        << "Compute cap.   = "
        << prop.major
        << "."
        << prop.minor
        << "\n";


    cout
        << "Global memory  = "
        << (
            prop.totalGlobalMem
            /
            (
                1024.0
                *
                1024.0
            )
        )
        << " MiB\n\n";


    // ========================================================
    // Matrix check + copy constants
    // ========================================================

    checkMatricesHost();


    copyConstantsToGPU();


    // ========================================================
    // Square length
    // ========================================================

    const double squareLength =
        sqrt(
            PI
            *
            Rfinal
            *
            Rfinal

            +

            (
                4.0
                -
                PI
            )
            *
            W
            *
            W
            /
            4.0
        );


    // ========================================================
    // Memory sizes
    // ========================================================

    size_t bytes =
        N
        *
        sizeof(double);


    size_t bytesQ5 =
        static_cast<size_t>(
            Q5
        )
        *
        N
        *
        sizeof(double);


    size_t bytesQ9 =
        static_cast<size_t>(
            Q9
        )
        *
        N
        *
        sizeof(double);


    // ========================================================
    // GPU macroscopic arrays
    // ========================================================

    double* d_phi;

    double* d_rho;

    double* d_mu;


    double* d_gradPhiX;

    double* d_gradPhiY;


    double* d_gradRhoX;

    double* d_gradRhoY;


    double* d_forceX;

    double* d_forceY;


    double* d_ux;

    double* d_uy;


    double* d_pressure;


    double* d_dudx;

    double* d_dudy;

    double* d_dvdx;

    double* d_dvdy;


    // ========================================================
    // Old fields
    // ========================================================

    double* d_oldPhiUx;

    double* d_oldPhiUy;


    double* d_oldQxx;

    double* d_oldQxy;

    double* d_oldQyy;


    // ========================================================
    // Distribution functions
    // ========================================================

    double* d_gPhase;

    double* d_gPhasePost;

    double* d_gPhaseNext;


    double* d_g1;

    double* d_g2;


    double* d_g1Post;

    double* d_g2Post;


    double* d_g1Next;

    double* d_g2Next;


    // ========================================================
    // Allocate simple fields
    // ========================================================

    auto allocField =
        [&](
            double** p
        )
        {
            CUDA_CHECK(
                cudaMalloc(
                    reinterpret_cast<void**>(
                        p
                    ),
                    bytes
                )
            );
        };


    allocField(
        &d_phi
    );

    allocField(
        &d_rho
    );

    allocField(
        &d_mu
    );


    allocField(
        &d_gradPhiX
    );

    allocField(
        &d_gradPhiY
    );


    allocField(
        &d_gradRhoX
    );

    allocField(
        &d_gradRhoY
    );


    allocField(
        &d_forceX
    );

    allocField(
        &d_forceY
    );


    allocField(
        &d_ux
    );

    allocField(
        &d_uy
    );


    allocField(
        &d_pressure
    );


    allocField(
        &d_dudx
    );

    allocField(
        &d_dudy
    );

    allocField(
        &d_dvdx
    );

    allocField(
        &d_dvdy
    );


    allocField(
        &d_oldPhiUx
    );

    allocField(
        &d_oldPhiUy
    );


    allocField(
        &d_oldQxx
    );

    allocField(
        &d_oldQxy
    );

    allocField(
        &d_oldQyy
    );


    // ========================================================
    // Allocate distribution arrays
    // ========================================================

    CUDA_CHECK(
        cudaMalloc(
            reinterpret_cast<void**>(
                &d_gPhase
            ),
            bytesQ5
        )
    );


    CUDA_CHECK(
        cudaMalloc(
            reinterpret_cast<void**>(
                &d_gPhasePost
            ),
            bytesQ5
        )
    );


    CUDA_CHECK(
        cudaMalloc(
            reinterpret_cast<void**>(
                &d_gPhaseNext
            ),
            bytesQ5
        )
    );


    CUDA_CHECK(
        cudaMalloc(
            reinterpret_cast<void**>(
                &d_g1
            ),
            bytesQ9
        )
    );


    CUDA_CHECK(
        cudaMalloc(
            reinterpret_cast<void**>(
                &d_g2
            ),
            bytesQ9
        )
    );


    CUDA_CHECK(
        cudaMalloc(
            reinterpret_cast<void**>(
                &d_g1Post
            ),
            bytesQ9
        )
    );


    CUDA_CHECK(
        cudaMalloc(
            reinterpret_cast<void**>(
                &d_g2Post
            ),
            bytesQ9
        )
    );


    CUDA_CHECK(
        cudaMalloc(
            reinterpret_cast<void**>(
                &d_g1Next
            ),
            bytesQ9
        )
    );


    CUDA_CHECK(
        cudaMalloc(
            reinterpret_cast<void**>(
                &d_g2Next
            ),
            bytesQ9
        )
    );


    // ========================================================
    // NaN flag
    // ========================================================

    int* d_badFlag;

    int* d_badIndex;


    CUDA_CHECK(
        cudaMalloc(
            reinterpret_cast<void**>(
                &d_badFlag
            ),
            sizeof(int)
        )
    );


    CUDA_CHECK(
        cudaMalloc(
            reinterpret_cast<void**>(
                &d_badIndex
            ),
            sizeof(int)
        )
    );


    // ========================================================
    // Host diagnostic buffers
    // ========================================================

    double* h_phi =
        new double[N];

    double* h_rho =
        new double[N];

    double* h_mu =
        new double[N];

    double* h_forceX =
        new double[N];

    double* h_forceY =
        new double[N];

    double* h_ux =
        new double[N];

    double* h_uy =
        new double[N];

    double* h_pressure =
        new double[N];

    double* h_dudx =
        new double[N];

    double* h_dvdy =
        new double[N];


    // ========================================================
    // CUDA blocks
    // ========================================================

    dim3 block2D(
        16,
        16
    );


    dim3 grid2D(
        (
            NX
            +
            block2D.x
            -
            1
        )
        /
        block2D.x,

        (
            NY
            +
            block2D.y
            -
            1
        )
        /
        block2D.y
    );


    int threads =
        256;


    int blocks =
        (
            N
            +
            threads
            -
            1
        )
        /
        threads;


    // ========================================================
    // Helper:
    // recover all macroscopic values
    // ========================================================

    auto recoverAllCurrent =
        [&]()
        {
            kernelRecoverPhiDensity
                <<<blocks, threads>>>
                (
                    d_gPhase,
                    d_phi,
                    d_rho
                );


            checkKernel(
                "kernelRecoverPhiDensity"
            );


            kernelPhaseDerivativesAndForce
                <<<grid2D, block2D>>>
                (
                    d_phi,
                    d_mu,

                    d_gradPhiX,
                    d_gradPhiY,

                    d_forceX,
                    d_forceY
                );


            checkKernel(
                "kernelPhaseDerivativesAndForce"
            );


            kernelGradRho
                <<<grid2D, block2D>>>
                (
                    d_rho,

                    d_gradRhoX,
                    d_gradRhoY
                );


            checkKernel(
                "kernelGradRho"
            );


            kernelRecoverVelocity
                <<<blocks, threads>>>
                (
                    d_g1,
                    d_g2,

                    d_rho,

                    d_forceX,
                    d_forceY,

                    d_ux,
                    d_uy
                );


            checkKernel(
                "kernelRecoverVelocity"
            );


            kernelVelocityGradients
                <<<grid2D, block2D>>>
                (
                    d_ux,
                    d_uy,

                    d_dudx,
                    d_dudy,

                    d_dvdx,
                    d_dvdy
                );


            checkKernel(
                "kernelVelocityGradients"
            );


            kernelRecoverPressure
                <<<blocks, threads>>>
                (
                    d_g1,
                    d_g2,

                    d_rho,

                    d_ux,
                    d_uy,

                    d_gradRhoX,
                    d_gradRhoY,

                    d_pressure
                );


            checkKernel(
                "kernelRecoverPressure"
            );
        };


    // ========================================================
    // Initialization
    // ========================================================

    kernelInitializePhi
        <<<grid2D, block2D>>>
        (
            d_phi,
            squareLength
        );


    checkKernel(
        "kernelInitializePhi"
    );


    kernelDensityFromPhi
        <<<blocks, threads>>>
        (
            d_phi,
            d_rho
        );


    checkKernel(
        "kernelDensityFromPhi"
    );


    kernelPhaseDerivativesAndForce
        <<<grid2D, block2D>>>
        (
            d_phi,
            d_mu,

            d_gradPhiX,
            d_gradPhiY,

            d_forceX,
            d_forceY
        );


    checkKernel(
        "kernelPhaseDerivativesAndForce(init)"
    );


    kernelGradRho
        <<<grid2D, block2D>>>
        (
            d_rho,

            d_gradRhoX,
            d_gradRhoY
        );


    checkKernel(
        "kernelGradRho(init)"
    );


    kernelZeroVelocityPressure
        <<<blocks, threads>>>
        (
            d_ux,
            d_uy,

            d_pressure,

            d_dudx,
            d_dudy,

            d_dvdx,
            d_dvdy
        );


    checkKernel(
        "kernelZeroVelocityPressure"
    );


    kernelInitializePhaseDistribution
        <<<blocks, threads>>>
        (
            d_phi,
            d_mu,

            d_ux,
            d_uy,

            d_gPhase
        );


    checkKernel(
        "kernelInitializePhaseDistribution"
    );


    kernelInitializeFlowDistribution
        <<<blocks, threads>>>
        (
            d_rho,

            d_ux,
            d_uy,

            d_pressure,

            d_forceX,
            d_forceY,

            d_g1,
            d_g2
        );


    checkKernel(
        "kernelInitializeFlowDistribution"
    );


    // ========================================================
    // Recover initial macros
    // ========================================================

    recoverAllCurrent();


    // ========================================================
    // Initialize old values
    // ========================================================

    kernelInitializeOldFields
        <<<blocks, threads>>>
        (
            d_phi,
            d_rho,

            d_ux,
            d_uy,

            d_pressure,

            d_oldPhiUx,
            d_oldPhiUy,

            d_oldQxx,
            d_oldQxy,
            d_oldQyy
        );


    checkKernel(
        "kernelInitializeOldFields"
    );


    CUDA_CHECK(
        cudaDeviceSynchronize()
    );


    // ========================================================
    // Initial diagnostics
    // ========================================================

    copyDiagnosticsToHost(
        h_phi,
        h_rho,
        h_mu,

        h_forceX,
        h_forceY,

        h_ux,
        h_uy,

        h_pressure,

        h_dudx,
        h_dvdy,

        d_phi,
        d_rho,
        d_mu,

        d_forceX,
        d_forceY,

        d_ux,
        d_uy,

        d_pressure,

        d_dudx,
        d_dvdy
    );


    InitialConserved init;


    for (
        int id = 0;
        id < N;
        id++
    )
    {
        init.phaseMass
            +=
        h_phi[id];


        init.densityMass
            +=
        h_rho[id];


        init.momentumX
            +=
        h_rho[id]
        *
        h_ux[id];


        init.momentumY
            +=
        h_rho[id]
        *
        h_uy[id];
    }


    printDiagnosticsHost(
        0,

        h_phi,
        h_rho,
        h_mu,

        h_forceX,
        h_forceY,

        h_ux,
        h_uy,

        h_pressure,

        h_dudx,
        h_dvdy,

        init
    );


    savePhi(
        h_phi,
        0
    );


    // ========================================================
    // GPU timing
    // ========================================================

    cudaEvent_t evStart;

    cudaEvent_t evStop;


    CUDA_CHECK(
        cudaEventCreate(
            &evStart
        )
    );


    CUDA_CHECK(
        cudaEventCreate(
            &evStop
        )
    );


    CUDA_CHECK(
        cudaEventRecord(
            evStart
        )
    );


    // ========================================================
    // Time loop
    // ========================================================

    for (
        int step = 1;
        step <= nSteps;
        step++
    )
    {
        // ====================================================
        // 1. Phase collision
        // ====================================================

        kernelCollidePhase
            <<<blocks, threads>>>
            (
                d_gPhase,
                d_gPhasePost,

                d_phi,
                d_mu,

                d_ux,
                d_uy,

                d_oldPhiUx,
                d_oldPhiUy
            );


        checkKernel(
            "kernelCollidePhase"
        );


        // ====================================================
        // 2. Flow collision + save old fields
        // ====================================================

        kernelCollideFlowAndSaveOld
            <<<blocks, threads>>>
            (
                d_g1,
                d_g2,

                d_g1Post,
                d_g2Post,

                d_phi,
                d_rho,

                d_ux,
                d_uy,

                d_pressure,

                d_forceX,
                d_forceY,

                d_gradRhoX,
                d_gradRhoY,

                d_dudx,
                d_dudy,

                d_dvdx,
                d_dvdy,

                d_oldPhiUx,
                d_oldPhiUy,

                d_oldQxx,
                d_oldQxy,
                d_oldQyy
            );


        checkKernel(
            "kernelCollideFlowAndSaveOld"
        );


        // ====================================================
        // 3. Pull streaming
        // ====================================================

        kernelStreamAll
            <<<grid2D, block2D>>>
            (
                d_gPhasePost,

                d_g1Post,
                d_g2Post,

                d_gPhaseNext,

                d_g1Next,
                d_g2Next
            );


        checkKernel(
            "kernelStreamAll"
        );


        // ====================================================
        // 4. Swap distribution pointers
        // ====================================================

        std::swap(
            d_gPhase,
            d_gPhaseNext
        );


        std::swap(
            d_g1,
            d_g1Next
        );


        std::swap(
            d_g2,
            d_g2Next
        );


        // ====================================================
        // 5. Recover n+1 macroscopic fields
        // ====================================================

        recoverAllCurrent();


        // ====================================================
        // 6. Safety check
        // ====================================================

        if (
            step
            %
            CHECK_INTERVAL
            ==
            0
        )
        {
            int zero =
                0;


            CUDA_CHECK(
                cudaMemcpy(
                    d_badFlag,

                    &zero,

                    sizeof(int),

                    cudaMemcpyHostToDevice
                )
            );


            kernelCheckFinite
                <<<blocks, threads>>>
                (
                    d_phi,
                    d_rho,

                    d_ux,
                    d_uy,

                    d_pressure,

                    d_mu,

                    d_badFlag,
                    d_badIndex
                );


            checkKernel(
                "kernelCheckFinite"
            );


            int bad =
                0;


            CUDA_CHECK(
                cudaMemcpy(
                    &bad,

                    d_badFlag,

                    sizeof(int),

                    cudaMemcpyDeviceToHost
                )
            );


            if (
                bad
            )
            {
                int badId =
                    -1;


                CUDA_CHECK(
                    cudaMemcpy(
                        &badId,

                        d_badIndex,

                        sizeof(int),

                        cudaMemcpyDeviceToHost
                    )
                );


                cerr
                    << "\nNaN/Inf detected at step "
                    << step

                    << ", x="
                    << (
                        badId
                        %
                        NX
                    )

                    << ", y="
                    << (
                        badId
                        /
                        NX
                    )

                    << "\n";


                return
                    2;
            }
        }


        // ====================================================
        // 7. Output
        // ====================================================

        if (
            step
            %
            OUTPUT_INTERVAL
            ==
            0

            ||

            step
            ==
            nSteps
        )
        {
            CUDA_CHECK(
                cudaDeviceSynchronize()
            );


            copyDiagnosticsToHost(
                h_phi,
                h_rho,
                h_mu,

                h_forceX,
                h_forceY,

                h_ux,
                h_uy,

                h_pressure,

                h_dudx,
                h_dvdy,

                d_phi,
                d_rho,
                d_mu,

                d_forceX,
                d_forceY,

                d_ux,
                d_uy,

                d_pressure,

                d_dudx,
                d_dvdy
            );


            printDiagnosticsHost(
                step,

                h_phi,
                h_rho,
                h_mu,

                h_forceX,
                h_forceY,

                h_ux,
                h_uy,

                h_pressure,

                h_dudx,
                h_dvdy,

                init
            );


            savePhi(
                h_phi,
                step
            );
        }
    }


    // ========================================================
    // Stop timer
    // ========================================================

    CUDA_CHECK(
        cudaEventRecord(
            evStop
        )
    );


    CUDA_CHECK(
        cudaEventSynchronize(
            evStop
        )
    );


    float elapsedMs =
        0.0f;


    CUDA_CHECK(
        cudaEventElapsedTime(
            &elapsedMs,
            evStart,
            evStop
        )
    );


    // ========================================================
    // Timing report
    // ========================================================

    cout
        << "\n====================================================\n";


    cout
        << "CUDA simulation finished.\n";


    cout
        << "steps          = "
        << nSteps
        << "\n";


    cout
        << "GPU loop time  = "
        << elapsedMs
        /
        1000.0
        << " s\n";


    cout
        << "steps / second = "
        << nSteps
        /
        (
            elapsedMs
            /
            1000.0
        )
        << "\n";


    cout
        << "====================================================\n";


    // ========================================================
    // Destroy events
    // ========================================================

    CUDA_CHECK(
        cudaEventDestroy(
            evStart
        )
    );


    CUDA_CHECK(
        cudaEventDestroy(
            evStop
        )
    );


    // ========================================================
    // Free GPU memory
    // ========================================================

    cudaFree(
        d_phi
    );

    cudaFree(
        d_rho
    );

    cudaFree(
        d_mu
    );


    cudaFree(
        d_gradPhiX
    );

    cudaFree(
        d_gradPhiY
    );


    cudaFree(
        d_gradRhoX
    );

    cudaFree(
        d_gradRhoY
    );


    cudaFree(
        d_forceX
    );

    cudaFree(
        d_forceY
    );


    cudaFree(
        d_ux
    );

    cudaFree(
        d_uy
    );


    cudaFree(
        d_pressure
    );


    cudaFree(
        d_dudx
    );

    cudaFree(
        d_dudy
    );

    cudaFree(
        d_dvdx
    );

    cudaFree(
        d_dvdy
    );


    cudaFree(
        d_oldPhiUx
    );

    cudaFree(
        d_oldPhiUy
    );


    cudaFree(
        d_oldQxx
    );

    cudaFree(
        d_oldQxy
    );

    cudaFree(
        d_oldQyy
    );


    cudaFree(
        d_gPhase
    );

    cudaFree(
        d_gPhasePost
    );

    cudaFree(
        d_gPhaseNext
    );


    cudaFree(
        d_g1
    );

    cudaFree(
        d_g2
    );


    cudaFree(
        d_g1Post
    );

    cudaFree(
        d_g2Post
    );


    cudaFree(
        d_g1Next
    );

    cudaFree(
        d_g2Next
    );


    cudaFree(
        d_badFlag
    );

    cudaFree(
        d_badIndex
    );


    // ========================================================
    // Free host memory
    // ========================================================

    delete[]
        h_phi;

    delete[]
        h_rho;

    delete[]
        h_mu;


    delete[]
        h_forceX;

    delete[]
        h_forceY;


    delete[]
        h_ux;

    delete[]
        h_uy;


    delete[]
        h_pressure;


    delete[]
        h_dudx;

    delete[]
        h_dvdy;


    return
        0;
}
