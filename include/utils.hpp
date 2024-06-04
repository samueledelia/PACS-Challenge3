#ifndef UTILS_HPP
#define UTILS_HPP

#include "Mesh.hpp"
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <muParser.h>
#include <algorithm>
#include <fstream>

// Function to generate VTK file for flat vector of values
void generateVTK(const std::vector<double>& values, double x0, double y0, size_t nx, size_t ny, double hx, double hy, const std::string& extra) {
    std::string filename = "../VTK/solution" + extra + ".vtk"; // name of the solution file
    std::ofstream svtkFile(filename); // open file

    // check if the files were opened
    if (!svtkFile.is_open()) {
        std::cerr << "Error: could not open file " << filename << std::endl;
        return;
    }

    // write VTK header 
    svtkFile << "# vtk DataFile Version 3.0\n";
    svtkFile << "Discrete solution\n";
    svtkFile << "ASCII\n";                                

    // write grid data 
    svtkFile << "DATASET STRUCTURED_GRID\n";
    svtkFile << "DIMENSIONS " << nx << " " << ny << " 1\n";
    svtkFile << "POINTS " << nx * ny << " double\n";

    // Write the points
    for (size_t j = 0; j < ny; j++) {
        for (size_t i = 0; i < nx; i++) {
            svtkFile << x0 + i * hx << " " << y0 + j * hy << " 0\n"; // Z-coordinate should be zero
        }
    }

    // include scalar field for coloring solution
    svtkFile << "POINT_DATA " << nx * ny << "\n";
    svtkFile << "SCALARS solution double 1\n";
    svtkFile << "LOOKUP_TABLE default\n";
    for (size_t j = 0; j < ny; j++) {
        for (size_t i = 0; i < nx; i++) {
            svtkFile << values[j * nx + i] << "\n"; // Correct indexing
        }
    }

    std::cout << "VTK produced: " << filename << std::endl;
    svtkFile.close();
}

// Function to generate VTK file for vector of vectors of values
void generateVTK(const std::vector<std::vector<double>>& values, double x0, double y0, size_t nx, size_t ny, double hx, double hy, const std::string& extra) {
    // Flatten the vector of vectors and rely on writing of vector as values type
    std::vector<double> vecValues;
    for (const auto& row : values) 
        vecValues.insert(vecValues.end(), row.begin(), row.end());
    generateVTK(vecValues, x0, y0, nx, ny, hx, hy, extra);
}
// Function to compute the L2 norm difference for vector of vectors
double norm(const std::vector<std::vector<double>>& m1, const std::vector<std::vector<double>>& m2, double h) {
    size_t nRows = m1.size();

    if (nRows != m2.size()) // Check row size compatibility
        throw std::invalid_argument("norm size mismatch - row");

    double ss = 0.0; // Initialize sum of squares

    for (size_t i = 0; i < m1.size(); ++i) {
        size_t nCols = m1[i].size();

        if (nCols != m2[i].size())  // Check column size compatibility
            throw std::invalid_argument("norm size mismatch - col");

        for (size_t j = 0; j < nCols; ++j) {
            ss += (m1[i][j] - m2[i][j]) * (m1[i][j] - m2[i][j]); // Sum of squares
        }
    }

    return std::sqrt(h * ss);
}

// Function to compute the L2 norm difference for flat vectors
double norm(const std::vector<double>& v1, const std::vector<double>& v2, double h, size_t nR, size_t nC) {
    // Transform flat vectors into vectors of vectors
    std::vector<std::vector<double>> m1(nR, std::vector<double>(nC, 0));
    std::vector<std::vector<double>> m2(nR, std::vector<double>(nC, 0));

    for (size_t i = 0; i < nR; ++i) {
        for (size_t j = 0; j < nC; ++j) {
            m1[i][j] = v1[i * nC + j];
            m2[i][j] = v2[i * nC + j];
        }
    }

    // Call the overloaded norm function for vector of vectors
    return norm(m1, m2, h);
}

// compare the given solution vectorial form to the exact one wrt l2 norm
void compareSolution(const Mesh & mesh, const std::vector<double> & solution, const std::string & exact, int nProcess) {
    mu::Parser f; 

    // parser variable
    double x;
    double y;

    // set the parser
    try {
        // set parser references variable
        f.DefineVar("x", &x); 
        f.DefineVar("y", &y);
        f.SetExpr(exact); ///< set parser expression
    } catch (mu::Parser::exception_type &e) {
        std::cerr << "Error: " << e.GetMsg() << std::endl;
    }

    // generate discrete solution using the exact result
    size_t nx = mesh.getnx();
    size_t ny = mesh.getny();

    std::vector<double> exacSol(nx * ny, 0);
    for (size_t i = 0 ; i < ny ; ++i) {
        for (size_t j = 0 ; i < nx ; ++j) {
            x = mesh(i, j).getX(); // coordinate x value on the mesh, offset wrt to first row of the current process
            y = mesh(i, j).getY(); // coordinate y value on the mesh, offset wrt to first row of the current process
            exacSol[i * nx + j] = f.Eval(); // exact value on point (i,j)
        }
    }

    double error = norm(solution, exacSol, std::max(mesh.gethx(), mesh.gethy()), ny, nx); // error between exact and computed solution
    std::cout << "Error using " << nx << " points along X, " << ny << " points along Y and " << nProcess << " processes: " << error << std::endl;
}

// compare the given solution vectorial form to the exact one wrt l2 norm
void compareSolution(const Mesh & mesh, const std::vector<std::vector<double>> & solution, const std::string & exact, int nProcess) {
    // flatten local solution into a vector for gathering
    std::vector<double> vecSol;
    for (const auto& row : solution) 
        vecSol.insert(vecSol.end(), row.begin(), row.end());
    compareSolution(mesh, vecSol, exact, nProcess);       
}

#endif



/*#ifndef UTILS_HPP
#define UTILS_HPP

#include "Mesh.hpp"
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <muParser.h>
#include <algorithm>
#include <fstream>


// generate VTK file of the solution
    void generateVTK(const std::vector<double> & values, double x0, double y0, size_t nx, size_t ny, double hx, double hy,  const std::string & extra){

        std::string filename = "../VTK/solution"+extra+".vtk"; // name of the solution file

        std::ofstream svtkFile(filename); // open file

        // check if the files were opened
        if (!svtkFile.is_open()) {
            std::cerr << "Error: could not open file " << filename << std::endl;
            return;
        }

        // write VTK header 
        svtkFile << "# vtk DataFile Version 3.0\n";
        svtkFile << "Discrete solution\n";
        svtkFile << "ASCII\n";                                

        // write grid data 
        svtkFile << "DATASET STRUCTURED_GRID\n";
        svtkFile << "DIMENSIONS " << nx << " " << ny << " 1\n";
        svtkFile << "POINTS " << nx * ny << " double\n";


        // Write the points and scalar values
        for (size_t j = 0; j < ny; j++) {
            for (size_t i = 0; i < nx; i++) {
                svtkFile << x0 + i * hx << " " << y0 + j * hy << " " << values[i*nx+j] << "\n";
            }
        }

        // include scalar field for coloring solution
        svtkFile << "POINT_DATA " << nx * ny << "\n";
        svtkFile << "SCALARS solution double 1\n";
        svtkFile << "LOOKUP_TABLE default\n";
        for (size_t j = 0; j < ny; j++) {
            for (size_t i = 0; i < nx; i++) {
                svtkFile << values[i*nx+j] << "\n";
            }
        }

        std::cout << "VTK produced: " << filename << std::endl;
        svtkFile.close();

    };


    void generateVTK(const std::vector<std::vector<double>> & values, double x0, double y0, size_t nx, size_t ny, double hx, double hy,  const std::string & extra){

        // flatten the vector of vector and rely on writing of vector as values type
        std::vector<double> vecValues;
        for (const auto& row : values) 
            vecValues.insert(vecValues.end(), row.begin(), row.end());
        generateVTK(vecValues, x0, y0, nx, ny, hx, hy, extra);

    }

    // l2 matrix norm of the difference
    double norm(const std::vector<double> & v1 , const std::vector<double> & v2, double h, size_t nR, size_t nC) {

        // transform vector in vector of vector and use defintion for this type
        std::vector<std::vector<double>> m1;
        m1.resize(nR,std::vector<double>(nC,0));
        std::vector<std::vector<double>> m2;
        m2.resize(nR,std::vector<double>(nC,0));

        for (size_t i = 0; i < nR; ++i) {
            for(size_t j = 0; j < nC; ++j) {

                m1[i][j] = v1[i * nC + j];
                m2[i][j] = v2[i * nC + j];

            }
        }

        return norm(m1, m2, h); 

    };

    // l2 matrix norm of the difference
    double norm(const std::vector<std::vector<double>> & m1 , const std::vector<std::vector<double>> & m2, double h) {

        size_t nRows = m1.size();

        if(nRows != m2.size()) // check row size compatibility
            throw std::invalid_argument("norm size mismatch - row");

        double ss{0.}; // init SS

        for(size_t i = 0 ; i < m1.size() ; ++i){

            size_t nCols = m1[i].size(); 

            if(nCols != m2[i].size())  // check col size compatibility
                throw std::invalid_argument("norm size mismatch - col");

            for(size_t j = 0 ; j < nCols ; ++j){
               ss += (m1[i][j]-m2[i][j])*(m1[i][j]-m2[i][j]); // sum of SS 
            }
        }

        return sqrt(h*ss);

    };

    // compare the given solution vectorial form to the exact one wrt l2 norm
    void compareSolution(const Mesh & mesh, const std::vector<double> & solution, const std::string & exact, int nProcess) {
       
       mu::Parser f; 

        // parser variable
        double x;
        double y;

        // set the parser
        try {
            
            // set parser references variable
            f.DefineVar("x", &x); 
            f.DefineVar("y", &y);

            f.SetExpr(exact); ///< set parser expression

        }
        catch (mu::Parser::exception_type &e) {

            std::cerr << "Error: " << e.GetMsg() << std::endl;

        }

        // generate discrete solution using the exact result
        size_t nx = mesh.getnx();
        size_t ny = mesh.getny();

        std::vector<double> exacSol(nx*ny,0);
        for(size_t i = 0 ; i < ny ; ++i){
            for(size_t j = 0 ; j < nx ; ++j){

                x = mesh(i, j).getX(); // coordinate x value on the mesh, offset wrt to first row of the current process
                y = mesh(i, j).getY(); // coordinate y value on the mesh, offset wrt to first row of the current process

                exacSol[i*nx + j] = f.Eval(); //  exact value on point (i,j)

            }
        }

        double error = norm(solution, exacSol, std::max(mesh.gethx(), mesh.gethy()), ny, nx); // error between exact and computed solution
        std::cout << "Error using " << nx << " point along X, " << ny << " points along Y and " << nProcess << " processes: " << error << std::endl;

    };

    // compare the given solution vectorial form to the exact one wrt l2 norm
    void compareSolution(const Mesh & mesh, const std::vector<std::vector<double>> & solution, const std::string & exact, int nProcess) {
       
        // flatten local solution into a vector for gathering
        std::vector<double> vecSol;
        for (const auto& row : solution) 
            vecSol.insert(vecSol.end(), row.begin(), row.end());

        compareSolution(mesh, vecSol, exact, nProcess);       

    };

#endif*/