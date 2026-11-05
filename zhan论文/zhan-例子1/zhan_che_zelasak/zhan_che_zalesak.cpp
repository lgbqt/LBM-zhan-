#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

// ============================================================
// Zhan 2023 - D2Q5 CHE diagonal translation
//
// Paper-based:
//   D2Q5 transformation matrix: Appendix A
//   CHE moments: Appendix A
//   Evolution in moment space: Eq. (36)
//   W = 3 dx
//   dx = dt = 1
//   L0 = 200
//   U0 = 0.02
//   periodic boundaries
//   Pe = U0 L0 / M_phi
//   s2^0 = s1^0 for interface-capturing tests
//
// Implementation choices in this first version:
//   eta = 1
//   phiA = +1, phiB = -1
//   sigma = 0.04
//   circular tanh initial profile
// ============================================================


// ============================================================
// 1. Basic parameters
// ============================================================

const int NX = 200;
const int NY = 200;
const int Q  = 5;

const double dx = 1.0;
const double dt = 1.0;
const double c  = dx / dt;

const double cs2 = c * c / 3.0;

const double L0 = 200.0;

const double radius =80.0;
    
const double slotWidth = 15.0;
const double W = 3.0;

const double U0 = 0.02;




// ------------------------------------------------------------
// Paper target:
// Table 1 includes Pe = 125, 500, 2000
// First reproduce Pe = 125
// ------------------------------------------------------------

#ifndef PE_VALUE
#define PE_VALUE 80.0
#endif

const double Pe = PE_VALUE;


// ------------------------------------------------------------
// From paper definition:
// Pe = U0 * L0 / M_phi
// ------------------------------------------------------------

const double mobility =
    U0 * L0 / Pe;


// ------------------------------------------------------------
// Implementation choice
// ------------------------------------------------------------

const double eta = 1.0;


// ------------------------------------------------------------
// Paper relation:
//
// s1^0 = 1 / [ M_phi / (eta cs^2 dt) + 0.5 ]
//
// and in interface-capturing tests:
//
// s2^0 = s1^0
// ------------------------------------------------------------

const double s1 =
    1.0
    /
    (
        mobility
        /
        (eta * cs2 * dt)
        +
        0.5
    );

const double s2 = s1;


// ------------------------------------------------------------
// s0 is irrelevant for conserved moment if m0 = meq0.
// We use 1.0.
// ------------------------------------------------------------

const double s0 = 1.0;


// ============================================================
// 2. Phase-field parameters
// ============================================================

const double phiA = 1.0;
const double phiB = -1.0;


// ------------------------------------------------------------
// Implementation choice
// ------------------------------------------------------------

const double sigma = 0.04;


// ------------------------------------------------------------
// Zhan Eq. (4):
//
// beta = 12 sigma / [ W (phiA-phiB)^4 ]
//
// kappa = (3/2) sigma W / (phiA-phiB)^2
// ------------------------------------------------------------

const double deltaPhi =
    phiA - phiB;

const double betaPF =
    12.0
    *
    sigma
    /
    (
        W
        *
        deltaPhi
        *
        deltaPhi
        *
        deltaPhi
        *
        deltaPhi
    );

const double kappaPF =
    1.5
    *
    sigma
    *
    W
    /
    (
        deltaPhi
        *
        deltaPhi
    );


// ============================================================
// 3. D2Q5 discrete velocities
// ============================================================

const int cx[Q] =
{
     0,
     1,
     0,
    -1,
     0
};

const int cy[Q] =
{
     0,
     0,
     1,
     0,
    -1
};


// ============================================================
// 4. Zhan Appendix A D2Q5 matrix
// ============================================================

const double Mmat[Q][Q] =
{
    {1.0, 1.0, 1.0,  1.0,  1.0},
    {0.0, 1.0, 0.0, -1.0,  0.0},
    {0.0, 0.0, 1.0,  0.0, -1.0},
    {0.0, 1.0, 0.0,  1.0,  0.0},
    {0.0, 0.0, 1.0,  0.0,  1.0}
};


// ============================================================
// 5. Inverse of Zhan D2Q5 matrix
//
// From:
// m0 = g0+g1+g2+g3+g4
// m1 = g1-g3
// m2 = g2-g4
// m3 = g1+g3
// m4 = g2+g4
//
// therefore:
//
// g0 = m0-m3-m4
// g1 = (m1+m3)/2
// g2 = (m2+m4)/2
// g3 = (m3-m1)/2
// g4 = (m4-m2)/2
// ============================================================

const double Minv[Q][Q] =
{
    {1.0,  0.0,  0.0, -1.0, -1.0},
    {0.0,  0.5,  0.0,  0.5,  0.0},
    {0.0,  0.0,  0.5,  0.0,  0.5},
    {0.0, -0.5,  0.0,  0.5,  0.0},
    {0.0,  0.0, -0.5,  0.0,  0.5}
};


// ============================================================
// 6. Global arrays
// ============================================================

double ux[NY][NX];
double uy[NY][NX];

double phi[NY][NX];
double phiInitial[NY][NX];

double lapPhi[NY][NX];
double mu[NY][NX];

double oldPhiUx[NY][NX];
double oldPhiUy[NY][NX];

double g[Q][NY][NX];
double gPost[Q][NY][NX];
double gNext[Q][NY][NX];


// ============================================================
// 7. Function declarations
// ============================================================

void initializePhi();

void computeLaplacian();

void computeChemicalPotential();

void initializeDistribution();

void recoverPhi();

void collide();

void streamPeriodic();

double computeErrorL2();

double computeMass();

bool checkFinitePhi();

double checkMatrix();

void initializeVelocity();

void printDiagnostics(
    int step,
    double initialMass
);

void savePhiCSV(
    const char* fileName,
    double field[NY][NX]
);


// ============================================================
// 8. Main
// ============================================================

int main()
{
    std::cout
        << std::scientific
        << std::setprecision(12);

    std::cout
        << "============================================\n";

    std::cout
        << "Zhan 2023 D2Q5-CHE diagonal translation\n";

    std::cout
        << "============================================\n\n";


    std::cout
        << "NX                 = "
        << NX
        << "\n";

    std::cout
        << "NY                 = "
        << NY
        << "\n";

    std::cout
        << "L0                 = "
        << L0
        << "\n";

    std::cout
        << "Radius             = "
        << radius
        << "\n";

    std::cout
        << "W                  = "
        << W
        << "\n";

    std::cout
        << "U0                 = "
        << U0
        << "\n";

    std::cout
        << "Pe                 = "
        << Pe
        << "\n";

    std::cout
        << "Mobility           = "
        << mobility
        << "\n";

    std::cout
        << "eta                = "
        << eta
        << "  [implementation choice]\n";

    std::cout
        << "sigma              = "
        << sigma
        << "  [implementation choice]\n";

    std::cout
        << "betaPF             = "
        << betaPF
        << "\n";

    std::cout
        << "kappaPF            = "
        << kappaPF
        << "\n";

    std::cout
        << "s0                 = "
        << s0
        << "\n";

    std::cout
        << "s1                 = "
        << s1
        << "\n";

    std::cout
        << "s2                 = "
        << s2
        << "\n\n";


    const double matrixError =
        checkMatrix();

    std::cout
        << "Max |M*Minv-I|     = "
        << matrixError
        << "\n";

    if (
        matrixError
        <
        1.0e-12
    )
    {
        std::cout
            << "Matrix check       = PASS\n\n";
    }
    else
    {
        std::cout
            << "Matrix check       = FAIL\n";

        return 1;
    }


    // --------------------------------------------------------
    // Initialization
    // --------------------------------------------------------

    initializeVelocity();
    initializePhi();

    computeLaplacian();

    computeChemicalPotential();

    initializeDistribution();


    const double initialMass =
        computeMass();


    savePhiCSV(
        "zhan_che_phi_initial.csv",
        phi
    );


    printDiagnostics(
        0,
        initialMass
    );


    // --------------------------------------------------------
    // One period:
    //
    // T = L0 / U0 = 10000
    // --------------------------------------------------------

    const int numberOfSteps =
        static_cast<int>(
            L0 / U0
        );


    std::cout
        << "\nNumber of steps    = "
        << numberOfSteps
        << "\n";

    std::cout
        << "Start iteration...\n\n";


    // ========================================================
    // Time loop
    // ========================================================

    for (
        int step = 1;
        step <= numberOfSteps;
        step++
    )
    {
        // ----------------------------------------------------
        // phi^n -> mu^n
        // ----------------------------------------------------

        computeLaplacian();

        computeChemicalPotential();


        // ----------------------------------------------------
        // collision
        // ----------------------------------------------------

        collide();


        // ----------------------------------------------------
        // streaming
        // ----------------------------------------------------

        streamPeriodic();


        // ----------------------------------------------------
        // g <- gNext
        // ----------------------------------------------------

        for (
            int i = 0;
            i < Q;
            i++
        )
        {
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
                    g[i][y][x]
                        =
                    gNext[i][y][x];
                }
            }
        }


        // ----------------------------------------------------
        // phi^{n+1} = sum_i g_i
        //
        // Here S = 0, so Eq. (35a) becomes simply sum g_i.
        // ----------------------------------------------------

        recoverPhi();


        if (!checkFinitePhi())
        {
            std::cout
                << "ERROR: non-finite phi at step = "
                << step
                << "\n";

            return 2;
        }


        if (
            step == 10
            ||
            step == 100
            ||
            step == 500
            ||
            step == 1000
            ||
            step == 2500
            ||
            step == 5000
            ||
            step == 7500
            ||
            step == 10000
        )
        {
            printDiagnostics(
                step,
                initialMass
            );
        }
    }


    // ========================================================
    // Final result
    // ========================================================

    const double errorL2 =
        computeErrorL2();


    std::cout
        << "\n============================================\n";

    std::cout
        << "Final result\n";

    std::cout
        << "============================================\n";

    std::cout
        << "Relative L2 error = "
        << errorL2
        << "\n";

    double paperTarget = 0.0;

    if (std::fabs(Pe - 80.0)<1.0e-12)
    {
        paperTarget =
            0.1608;
    }
    else if (std::fabs(Pe - 400.0)<1.0e-12)    
    {
        paperTarget =
            0.1728;
    }
    else if (std::fabs(Pe - 800.0)<1.0e-12) 
    {
        paperTarget =
            0.1744;
    }
    else if (std::fabs(Pe - 4000.0)<1.0e-12)
    {
        paperTarget =
            0.1808;
    }
    std::cout
        << "Paper CHE target = "
        << paperTarget
        << "\n";


    savePhiCSV(
        "zhan_che_phi_final.csv",
        phi
    );


    std::cout
        << "============================================\n";

    return 0;
}


// ============================================================
// 9. Initialize circular interface
//
// Implementation choice:
// use a standard diffuse circular tanh profile
// ============================================================

void initializePhi()
{
    const double xc =
        0.5 * NX;

    const double yc =
        0.5 * NY;

    const double halfSlot =
        0.5 * slotWidth;


    for (int y = 0; y < NY; y++)
    {
        for (int x = 0; x < NX; x++)
        {
            const double xp =
                x + 0.5;

            const double yp =
                y + 0.5;


            // --------------------------------------------
            // 圆盘
            // --------------------------------------------

            const double rx =
                xp - xc;

            const double ry =
                yp - yc;

            const double r =
                std::sqrt(
                    rx * rx
                    +
                    ry * ry
                );


            const bool insideCircle =
                (
                    r <= radius
                );


            // --------------------------------------------
            // 槽
            //
            // 暂时仍保持：
            // 槽宽 = 15
            // 槽底 = 圆心
            // 槽向上延伸至圆周
            // --------------------------------------------

            const bool insideSlot =
                (
                    std::fabs(
                        xp - xc
                    )
                    <=
                    halfSlot
                    &&
                    yp >= yc
                );


            // --------------------------------------------
            // sharp Zalesak disk
            //
            // 圆内且不在槽内 -> phiA
            // 其他 -> phiB
            // --------------------------------------------

            if (
                insideCircle
                &&
                !insideSlot
            )
            {
                phi[y][x] =
                    phiA;
            }
            else
            {
                phi[y][x] =
                    phiB;
            }


            phiInitial[y][x] =
                phi[y][x];


            oldPhiUx[y][x] =
                phi[y][x]
                *
                ux[y][x];

            oldPhiUy[y][x] =
                phi[y][x]
                *
                uy[y][x];
        }
    }
}

// ============================================================
// 10. Laplacian
//
// D2Q5 five-point Laplacian
// ============================================================

void computeLaplacian()
{
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
            const int xp =
                (x + 1) % NX;

            const int xm =
                (x - 1 + NX) % NX;

            const int yp =
                (y + 1) % NY;

            const int ym =
                (y - 1 + NY) % NY;


            lapPhi[y][x]
                =
            phi[y][xp]
            +
            phi[y][xm]
            +
            phi[yp][x]
            +
            phi[ym][x]
            -
            4.0
            *
            phi[y][x];
        }
    }
}


// ============================================================
// 11. Chemical potential
//
// Zhan Eq. (3)
//
// mu =
// 4 beta (phi-phiA)(phi-phiB)
// [phi-(phiA+phiB)/2]
// - kappa lap(phi)
// ============================================================

void computeChemicalPotential()
{
    const double midPhi =
        0.5
        *
        (
            phiA
            +
            phiB
        );


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
            const double p =
                phi[y][x];


            mu[y][x]
                =
            4.0
            *
            betaPF
            *
            (
                p - phiA
            )
            *
            (
                p - phiB
            )
            *
            (
                p - midPhi
            )
            -
            kappaPF
            *
            lapPhi[y][x];
        }
    }
}


// ============================================================
// 12. Initialize distributions from equilibrium moments
// ============================================================

void initializeDistribution()
{
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
            double meq[Q];


            // Zhan Appendix A CHE moments

            meq[0]
                =
            phi[y][x];

            meq[1]
                =
            phi[y][x]
            *
            ux[y][x]
            /
            c;

            meq[2]
                =
            phi[y][x]
            *
            uy[y][x]
            /
            c;

            meq[3]
                =
            eta
            *
            mu[y][x]
            /
            3.0;

            meq[4]
                =
            eta
            *
            mu[y][x]
            /
            3.0;


            for (
                int i = 0;
                i < Q;
                i++
            )
            {
                double value = 0.0;


                for (
                    int a = 0;
                    a < Q;
                    a++
                )
                {
                    value +=
                        Minv[i][a]
                        *
                        meq[a];
                }


                g[i][y][x]
                    =
                value;

                gPost[i][y][x]
                    =
                value;

                gNext[i][y][x]
                    =
                value;
            }
        }
    }
}


// ============================================================
// 13. Recover phi
// ============================================================

void recoverPhi()
{
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
            double sum = 0.0;


            for (
                int i = 0;
                i < Q;
                i++
            )
            {
                sum +=
                    g[i][y][x];
            }


            phi[y][x]
                =
            sum;
        }
    }
}


// ============================================================
// 14. Collision
//
// Zhan Eq. (36), for CHE and S=0:
//
// m* =
// (I-S)m
// + S meq
// + dt (I-S/2) mG
// ============================================================

void collide()
{
    const double relax[Q] =
    {
        s0,
        s1,
        s1,
        s2,
        s2
    };


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
            // ------------------------------------------------
            // Explicit temporal derivative:
            //
            // d_t(phi*u)
            // ------------------------------------------------

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


            // ------------------------------------------------
            // Current moments m = M g
            // ------------------------------------------------

            double m[Q] =
            {
                0.0,
                0.0,
                0.0,
                0.0,
                0.0
            };


            for (
                int a = 0;
                a < Q;
                a++
            )
            {
                for (
                    int i = 0;
                    i < Q;
                    i++
                )
                {
                    m[a]
                        +=
                    Mmat[a][i]
                    *
                    g[i][y][x];
                }
            }


            // ------------------------------------------------
            // Zhan equilibrium moments
            // ------------------------------------------------

            double meq[Q];


            meq[0]
                =
            phi[y][x];

            meq[1]
                =
            phi[y][x]
            *
            ux[y][x]
            /
            c;

            meq[2]
                =
            phi[y][x]
            *
            uy[y][x]
            /
            c;

            meq[3]
                =
            eta
            *
            mu[y][x]
            /
            3.0;

            meq[4]
                =
            eta
            *
            mu[y][x]
            /
            3.0;


            // ------------------------------------------------
            // Zhan Appendix A auxiliary moments
            //
            // mG =
            // [0,
            //  d_t(phi u)/c,
            //  d_t(phi v)/c,
            //  0,
            //  0]^T
            // ------------------------------------------------

            double mG[Q];


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
            // Collision in moment space
            // ------------------------------------------------

            double mPost[Q];


            for (
                int a = 0;
                a < Q;
                a++
            )
            {
                mPost[a]
                    =
                (
                    1.0
                    -
                    relax[a]
                )
                *
                m[a]
                +
                relax[a]
                *
                meq[a]
                +
                dt
                *
                (
                    1.0
                    -
                    0.5
                    *
                    relax[a]
                )
                *
                mG[a];
            }


            // ------------------------------------------------
            // Back to distribution space
            // ------------------------------------------------

            for (
                int i = 0;
                i < Q;
                i++
            )
            {
                double value = 0.0;


                for (
                    int a = 0;
                    a < Q;
                    a++
                )
                {
                    value +=
                        Minv[i][a]
                        *
                        mPost[a];
                }


                gPost[i][y][x]
                    =
                value;
            }


            // ------------------------------------------------
            // Save current phi*u
            // ------------------------------------------------

            oldPhiUx[y][x]
                =
            currentPhiUx;

            oldPhiUy[y][x]
                =
            currentPhiUy;
        }
    }
}


// ============================================================
// 15. Periodic streaming
// ============================================================

void streamPeriodic()
{
    for (
        int i = 0;
        i < Q;
        i++
    )
    {
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
                int xn =
                    x + cx[i];

                int yn =
                    y + cy[i];


                if (xn < 0)
                {
                    xn += NX;
                }

                if (xn >= NX)
                {
                    xn -= NX;
                }

                if (yn < 0)
                {
                    yn += NY;
                }

                if (yn >= NY)
                {
                    yn -= NY;
                }


                gNext[i][yn][xn]
                    =
                gPost[i][y][x];
            }
        }
    }
}


// ============================================================
// 16. Relative L2 error
//
// Zhan Eq. (44)
// ============================================================

double computeErrorL2()
{
    double numerator = 0.0;
    double denominator = 0.0;


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
            const double diff =
                phi[y][x]
                -
                phiInitial[y][x];


            numerator +=
                diff
                *
                diff;


            denominator +=
                phiInitial[y][x]
                *
                phiInitial[y][x];
        }
    }


    return
        std::sqrt(
            numerator
            /
            denominator
        );
}


// ============================================================
// 17. Mass
// ============================================================

double computeMass()
{
    double mass = 0.0;


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
            mass +=
                phi[y][x];
        }
    }


    return mass;
}


// ============================================================
// 18. Finite check
// ============================================================

bool checkFinitePhi()
{
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
            if (
                !std::isfinite(
                    phi[y][x]
                )
            )
            {
                return false;
            }
        }
    }


    return true;
}


// ============================================================
// 19. Matrix check
// ============================================================

double checkMatrix()
{
    double maxError = 0.0;


    for (
        int i = 0;
        i < Q;
        i++
    )
    {
        for (
            int j = 0;
            j < Q;
            j++
        )
        {
            double value = 0.0;


            for (
                int k = 0;
                k < Q;
                k++
            )
            {
                value +=
                    Mmat[i][k]
                    *
                    Minv[k][j];
            }


            double exact = 0.0;


            if (i == j)
            {
                exact = 1.0;
            }


            const double error =
                std::fabs(
                    value
                    -
                    exact
                );


            if (
                error
                >
                maxError
            )
            {
                maxError =
                    error;
            }
        }
    }


    return maxError;
}


// ============================================================
// 20. Diagnostics
// ============================================================

void printDiagnostics(
    int step,
    double initialMass
)
{
    double minPhi =
        phi[0][0];

    double maxPhi =
        phi[0][0];

    double maxAbsMu = 0.0;


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
            if (
                phi[y][x]
                <
                minPhi
            )
            {
                minPhi =
                    phi[y][x];
            }


            if (
                phi[y][x]
                >
                maxPhi
            )
            {
                maxPhi =
                    phi[y][x];
            }


            const double absMu =
                std::fabs(
                    mu[y][x]
                );


            if (
                absMu
                >
                maxAbsMu
            )
            {
                maxAbsMu =
                    absMu;
            }
        }
    }


    const double mass =
        computeMass();


    const double relativeMassChange =
        (
            mass
            -
            initialMass
        )
        /
        std::fabs(
            initialMass
        );


    std::cout
        << "Step = "
        << step
        << "\n";

    std::cout
        << "Relative mass change = "
        << relativeMassChange
        << "\n";

    std::cout
        << "min(phi) = "
        << minPhi
        << "\n";

    std::cout
        << "max(phi) = "
        << maxPhi
        << "\n";

    std::cout
        << "max|mu|  = "
        << maxAbsMu
        << "\n";
}


// ============================================================
// 21. Save CSV
// ============================================================

void savePhiCSV(
    const char* fileName,
    double field[NY][NX]
)
{
    std::ofstream file(
        fileName
    );


    file
        << std::scientific
        << std::setprecision(12);


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
            file
                << field[y][x];


            if (
                x
                <
                NX - 1
            )
            {
                file
                    << ",";
            }
        }


        file
            << "\n";
    }


    file.close();
}
void initializeVelocity()
{
    const double pi =
        3.14159265358979323846;

    for (int y = 0; y < NY; y++)
    {
        for (int x = 0; x < NX; x++)
        {
            const double xp = x + 0.5;

            const double yp = y + 0.5;

            ux[y][x] = -U0 * 2.0 * pi * ( yp / L0 - 0.5);

            uy[y][x] = U0 * 2.0 * pi * ( xp / L0 - 0.5 );
        }
    }
}