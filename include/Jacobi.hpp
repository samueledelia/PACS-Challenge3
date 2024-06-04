#ifndef JACOBI_HPP
#define JACOBI_HPP

#include "Mesh.hpp"
#include "utils.hpp"
#include <muParser.h>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <mpi.h>

class Jacobi {
    using solution = std::vector<std::vector<double>>;

private:
    Mesh mesh; 
    solution sol;

    mu::Parser f; 
    mu::Parser fBound; 
    double x; 
    double y; 

    unsigned int nMax; 
    double tol; 

    size_t nx; 
    size_t ny; 
    double hx;
    double hy;

public:
    Jacobi(const Mesh & mesh_, const std::string & function, unsigned int nMax_, double tol_, const std::string & bCond = "0");

    std::vector<std::vector<double>> solve();
    inline solution getSolution() { return sol;}

private:
    void initSolution();
    inline void setDim(solution & s) { s.resize(ny, std::vector<double>(nx, 0)); };
    void setParser(mu::Parser & parser, const std::string & expr);
};



Jacobi::Jacobi(const Mesh & mesh_, const std::string & function, unsigned int nMax_, double tol_, const std::string & bCond) 
    : mesh{mesh_}, nMax{nMax_}, tol{tol_} {
    nx = mesh.getnx();
    hx = mesh.gethx();
    ny = mesh.getny();
    hy = mesh.gethy();

    setParser(this->f, function);
    setParser(this->fBound, bCond);

    initSolution();
}

void Jacobi::setParser(mu::Parser & parser, const std::string & expr) {
    try {
        parser.DefineVar("x", &x); 
        parser.DefineVar("y", &y);
        parser.SetExpr(expr);
    } catch (mu::Parser::exception_type &e) {
        std::cerr << "Error: " << e.GetMsg() << std::endl;
    }
}

void Jacobi::initSolution() {
    setDim(sol);
    for (size_t j = 0; j < nx; ++j) {
        x = mesh(0, j).getX();
        y = mesh(0, j).getY();
        sol[0][j] = fBound.Eval();
        x = mesh(ny - 1, j).getX();
        y = mesh(ny - 1, j).getY();
        sol[ny - 1][j] = fBound.Eval();
    }
    for (size_t i = 0; i < ny; ++i) {
        x = mesh(i, 0).getX();
        y = mesh(i, 0).getY();
        sol[i][0] = fBound.Eval();
        x = mesh(i, nx - 1).getX();
        y = mesh(i, nx - 1).getY();
        sol[i][nx - 1] = fBound.Eval();
    }
}

std::vector<std::vector<double>> Jacobi::solve() {
    solution newSol;
    setDim(newSol);
    double e = tol;
    for (size_t k = 0; k < nMax && e >= tol; ++k) {
        for (size_t i = 1; i < ny - 1; ++i) {
            for (size_t j = 1; j < nx - 1; ++j) {
                x = mesh(i, j).getX();
                y = mesh(i, j).getY();
                newSol[i][j] = 0.25 * (sol[i + 1][j] + sol[i - 1][j] + sol[i][j + 1] + sol[i][j - 1] + hx * hy * f.Eval());
            }
        }
        std::swap(sol, newSol);
        e = norm(newSol, sol, std::max(hx, hy));
    }
    std::string extra = "_solution";
    generateVTK(sol, mesh.getx0(), mesh.gety0(), nx, ny, hx, hy, extra);
    return sol;
}

std::vector<double> parallelSolve(const Mesh & mesh, std::string & function, unsigned int nMax, double tol) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    size_t nx = mesh.getnx();
    double hx = mesh.gethx();
    size_t ny = mesh.getny();
    double hy = mesh.gethy();

    if (ny / size < 2)
        throw std::invalid_argument("Too much processes wrt to rows, at least two rows at each process");
    
    mu::Parser f;
    double x;
    double y;
    try {
        f.DefineVar("x", &x); 
        f.DefineVar("y", &y);
        f.SetExpr(function);
    } catch (mu::Parser::exception_type &e) {
        std::cerr << "Error: " << e.GetMsg() << std::endl;
    }

    int allConverged = 0;
    int nBaseRows = ny / size;
    int nExtraRows = ny % size;
    int firstRow = rank * nBaseRows + std::min(rank, nExtraRows);
    int lastRow = firstRow + ((nExtraRows > rank) ? nBaseRows + 1 : nBaseRows) - 1;
    int localRows = lastRow - firstRow + 1;

    std::vector<std::vector<double>> localSol(localRows, std::vector<double>(nx, 0));
    std::vector<std::vector<double>> localNewSol(localRows, std::vector<double>(nx, 0));
    std::vector<double> belowRow(nx, 0);
    std::vector<double> aboveRow(nx, 0);

    for (size_t k = 0; k < nMax && allConverged != size; ++k) {
        allConverged = 0;
        if (rank != size - 1) {
            MPI_Sendrecv(localSol[localRows - 1].data(), nx, MPI_DOUBLE, rank + 1, 0, aboveRow.data(), nx, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        if (rank != 0) {
            MPI_Sendrecv(localSol[0].data(), nx, MPI_DOUBLE, rank - 1, 0, belowRow.data(), nx, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        if (rank != size - 1) {
            int i = localRows - 1;
            for (size_t j = 1; j < nx - 1; ++j) {
                x = mesh(lastRow, j).getX();
                y = mesh(lastRow, j).getY();
                localNewSol[i][j] = 0.25 * (aboveRow[j] + localSol[i - 1][j] + localSol[i][j + 1] + localSol[i][j - 1] + hx * hy * f.Eval());
            }
        }

        if (rank != 0) {
            int i = 0;
            for (size_t j = 1; j < nx - 1; ++j) {
                x = mesh(firstRow, j).getX();
                y = mesh(firstRow, j).getY();
                localNewSol[i][j] = 0.25 * (localSol[i + 1][j] + belowRow[j] + localSol[i][j + 1] + localSol[i][j - 1] + hx * hy * f.Eval());
            }
        }

        for (size_t i = 1; i < localRows - 1; ++i) {
            for (size_t j = 1; j < nx - 1; ++j) {
                x = mesh(firstRow + i, j).getX();
                y = mesh(firstRow + i, j).getY();
                localNewSol[i][j] = 0.25 * (localSol[i + 1][j] + localSol[i - 1][j] + localSol[i][j + 1] + localSol[i][j - 1] + hx * hy * f.Eval());
            }
        }

        std::swap(localSol, localNewSol);

        if (norm(localNewSol, localSol, std::max(hx, hy)) < tol) {
            allConverged = 1;
        }
        MPI_Allreduce(MPI_IN_PLACE, &allConverged, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    }

    std::vector<double> result;
    if (rank == 0) {
        result.reserve(ny * nx);
    }

    for (int i = 0; i < size; ++i) {
        if (rank == i) {
            for (int j = 0; j < localRows; ++j) {
                result.insert(result.end(), localSol[j].begin(), localSol[j].end());
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    return result;
}

#endif // JACOBI_HPP