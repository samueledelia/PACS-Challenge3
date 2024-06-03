#ifndef JACOBI_HPP
#define JACOBI_HPP

#include "Mesh.hpp"

#include <muParser.h>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <mpi.h>


    // Jacobi Iteration methoed for solving Laplace equation
    class Jacobi {
        
        // Discrete solution matrix
        using solution = std::vector<std::vector<double>>;

    private:

        Mesh mesh; 
        solution sol; // Discrete solution 

        mu::Parser f; // function for which solving
        mu::Parser fBound; // function used for setting boundary condition
        double x; // x variable for parser
        double y; // y variable for parser

        unsigned int nMax; // maximum number of iteration
        double tol; // tolerance of the method

        size_t nx; 
        size_t ny; 
        double hx;
        double hy;

    public:
        // Jocobi method constructor
        Jacobi(const Mesh & mesh_, const std::string & function, unsigned int nMax_, double tol_ , const std::string & bCond = "0");

        // Solver method
        std::vector<std::vector<double>> solve();

        // Solution getter
        inline solution getSolution() { return sol;}

    private:

        // Initialize solution to zero and set Dirichlet boundary condition
        void initSolution();

        // Set solution size according to mesh dimension
        inline void setDim(solution & s) { s.resize(ny,std::vector<double>(nx,0)); };

        // Tie x,y variable to the parser and set its expression
        void setParser(mu::Parser & parser, const std::string & expr);
        
};
    

Jacobi::Jacobi(const Mesh & mesh_, const std::string & function, unsigned int nMax_, double tol_, const std::string & bCond) : mesh{mesh_},nMax{nMax_},tol{tol_} {
        
        /// set number of points and spacing
        nx = mesh.getnx();
        hx = mesh.gethx();
        ny = mesh.getny();
        hy = mesh.gethy();

        /// define parsers expression
        setParser(this->f,function);
        setParser(this->fBound,bCond);

        initSolution(); // initialize the solution matrix
    };

    void Jacobi::setParser(mu::Parser & parser, const std::string & expr){
        try {
            
            /// set parser references variable
            parser.DefineVar("x", &x); 
            parser.DefineVar("y", &y);

            parser.SetExpr(expr); // set parser expression

        }
        catch (mu::Parser::exception_type &e) {

            std::cerr << "Error: " << e.GetMsg() << std::endl;

        }
};

    void Jacobi::initSolution(){

        setDim(sol); // adjust dimension

        // set boundary condition
        for(size_t j = 0 ; j < nx ; ++j){
            x = mesh(0,j).getX();
            y = mesh(0,j).getY();
            sol[0][j] = fBound.Eval(); // first row
            x = mesh(ny-1,j).getX();
            y = mesh(ny-1,j).getY();
            sol[ny-1][j] = fBound.Eval();; // last row
        }

        for(size_t i = 0 ; i < ny ; ++i){
            x = mesh(i,0).getX();
            y = mesh(i,0).getY();
            sol[i][0] = fBound.Eval(); // first column
            x = mesh(i,nx-1).getX();
            y = mesh(i,nx-1).getY();
            sol[i][nx-1] = fBound.Eval(); // last column
        }    

    };

    std::vector<std::vector<double>> Jacobi::solve(){
        
        solution newSol;
        setDim(newSol);
        double e = tol;

        // Jacobi iteration
        for(size_t k = 0 ; k < nMax && e >= tol; ++k){

            for(size_t i = 1 ; i < ny-1 ; ++i){
                for(size_t j = 1 ; j < nx-1; ++j){
                    x = mesh(i,j).getX(); ///< coordinate x value on the mesh
                    y = mesh(i,j).getY(); ///< coordinate y value on the mesh
                    newSol[i][j] = .25*(sol[i+1][j]+sol[i-1][j]+sol[i][j+1]+sol[i][j-1]+hx*hy*f.Eval()); // four point stencil
                }
            }
            std::swap(sol,newSol); // set current solution to the one calculated
            e = norm(newSol,sol,std::max(hx,hy)); // get the norm of the increment of updated solution, use h as th maximum between the spacings
        }

        generateVTK(sol,mesh.getx0(),mesh.gety0(),nx,ny,hx,hy); // generate VTK file

        return sol;

    };









    // parallel solver
    std::vector<double> parallelSolve(const Mesh & mesh, std::string & function, unsigned int nMax, double tol){

        // setup MPI variables
        int rank, size;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);

        // set number of points and spacing
        size_t nx = mesh.getnx();
        double hx = mesh.gethx();
        size_t ny = mesh.getny();
        double hy = mesh.gethy();

        // if not at least two row per process
        if(ny / size < 2)
            throw std::invalid_argument("Too much processes wrt to rows, at least two rows at each process");
        
        mu::Parser f; // parser of the function

        // parser variable
        double x;
        double y;

        // set the parser
        try {
            
            // set parser references variable
            f.DefineVar("x", &x); 
            f.DefineVar("y", &y);

            f.SetExpr(function); // set parser expression

        }
        catch (mu::Parser::exception_type &e) {

            std::cerr << "Error: " << e.GetMsg() << std::endl;

        }

        int allConverged = 0; // to check convergence of all processes

        int nBaseRows = ny / size; // common number of rows for process, even spread part
        int nExtraRows = ny % size; // extra number of rows to be accounted

        int firstRow = rank * nBaseRows + std::min(rank, nExtraRows); // set first row, consider the number of extra rows already assigned
        int lastRow = firstRow + ((nExtraRows > rank) ? nBaseRows + 1 : nBaseRows) - 1; // last row considring if it has extra rows to be considered

        int localRows = lastRow - firstRow + 1; // number of local rows

        std::vector<std::vector<double>> localSol; // local solution
        localSol.resize(localRows, std::vector<double>(nx,0)); // size solution

        std::vector<std::vector<double>> localNewSol; // local new solution for update
        localNewSol.resize(localRows, std::vector<double>(nx,0)); // size new solution

        std::vector<double> belowRow(nx,0); // for receiving row below evaluated by previous rank process
        std::vector<double> aboveRow(nx,0); // for receiving row above evaluated by succesive rank process

        // Jacobi iteration
        for(size_t k = 0; k < nMax && allConverged != size ; ++k){

            allConverged = 0; // reset convergence flag

            if(rank != size - 1){ // if not last process send solution of last of your row and receive the first row of the following process
                MPI_Sendrecv(localSol[localRows-1].data(), nx, MPI_DOUBLE, rank+1, 0, aboveRow.data(), nx, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            if(rank != 0){ // if not first process send solution of first of your row and receive the last row of the previous process
                MPI_Sendrecv(localSol[0].data(), nx, MPI_DOUBLE, rank-1, 0, belowRow.data(), nx, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }

            // last row to be considered for non last rank process (boundary)
            if(rank != size - 1){

                int i = localRows - 1;

                for(size_t j = 1 ; j < nx-1; ++j){

                    x = mesh(lastRow,j).getX(); // coordinate x value last row mesh correspondence
                    y = mesh(lastRow,j).getY(); // coordinate y value last row mesh correspondence

                    localNewSol[i][j] = .25*(aboveRow[j]+localSol[i-1][j]+localSol[i][j+1]+localSol[i][j-1]+hx*hy*f.Eval()); ///<  four point stencil
                }
            }
            
            // first row to be considered for non first process (boundary)
            if(rank != 0){

                int i = 0;

                for(size_t j = 1 ; j < nx-1; ++j){

                    x = mesh(firstRow,j).getX(); // coordinate x value first row mesh correspondence
                    y = mesh(firstRow,j).getY(); // coordinate y value first row mesh correspondence

                    //std::cout << "Process: " << rank << "- new sol below" << std::endl;
                    localNewSol[i][j] = .25*(localSol[i+1][j]+belowRow[j]+localSol[i][j+1]+localSol[i][j-1]+hx*hy*f.Eval()); //  four point stencil
                }
            }

            // middle rows
            for(size_t i = 1 ; i < localRows - 1 ; ++i){

                for(size_t j = 1 ; j < nx-1; ++j){

                    x = mesh(firstRow + i,j).getX(); // coordinate x value on the mesh, offset wrt to first row of the current process
                    y = mesh(firstRow + i,j).getY(); // coordinate y value on the mesh, offset wrt to first row of the current process

                    localNewSol[i][j] = .25*(localSol[i+1][j]+localSol[i-1][j]+localSol[i][j+1]+localSol[i][j-1]+hx*hy*f.Eval()); //  four point stencil
                }
            }

            std::swap(localSol,localNewSol); // set current solution to the one calculated

            if(norm(localNewSol,localSol,std::max(hx,hy)) < tol) // check local tolerance convergence
                allConverged = 1;

            MPI_Barrier(MPI_COMM_WORLD); 
            
            MPI_Allreduce(MPI_IN_PLACE, &allConverged, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); // sum to consider how much processes converged
            
        }

        MPI_Barrier(MPI_COMM_WORLD);

        // assemble the global solution
        std::vector<int> nReceive(size,0); // vector containing number of rows to be received from rank [.]
        std::vector<int> vDisplacement(size,0); // vector containing displacement of where to start writing recived row from rank [.]
        
        int nElem = localRows * nx; // total number of element to send
        int displacement = firstRow * nx; // evaluate relative displacement
        
        MPI_Allgather(&nElem, 1, MPI_INT, nReceive.data(), 1, MPI_INT, MPI_COMM_WORLD); // gather the number of rows to receive for each rank
        MPI_Allgather(&displacement, 1, MPI_INT, vDisplacement.data(), 1, MPI_INT, MPI_COMM_WORLD); // gather the number of rows to receive for each rank

        // global solution that collect local ones
        std::vector<double>  globalSol;
        globalSol.resize(nx*ny);

        // flatten local solution into a vector for gathering
        std::vector<double> vecLocalS;
        for (const auto& row : localSol) 
            vecLocalS.insert(vecLocalS.end(), row.begin(), row.end());
  
        MPI_Allgatherv(vecLocalS.data(), nElem, MPI_DOUBLE, globalSol.data(), nReceive.data(), vDisplacement.data(), MPI_DOUBLE, MPI_COMM_WORLD);
        
        if(rank == 0) ///< first process writes the vtk file
            generateVTK(globalSol, mesh.getMinX(), mesh.getMinY(), nx, ny, hx, hy, "_"+std::to_string(size)+"parallel_n_"+std::to_string(nx)+"_"+std::to_string(ny));
        
        return globalSol;
            
    };


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

        return norm(m1, m2, h); ///< multiply by h and square the SS according to instruction in ../doc

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


#endif