#include "Point.hpp"
#include "Mesh.hpp"
#include "Jacobi.hpp"
#include "json.hpp"
#include <iostream>
#include <mpi.h>

// Parameters for mesh, function and method option
struct Parameters{

  /// mesh  
  double x0; // min x
  double y0; // min y
  double xn; // max x
  double yn; // max y
  size_t nx; // number of points along x
  size_t ny; // number of points along y
  double hx; // x spacing
  double hy; // y spacing

  /// iteration option
  unsigned int maxIter; // number of maximum iteration
  double tol; // tolerance

  /// functions
  std::string f; // function for solving
  std::string fBound; // function for boundary condition
  std::string exactSol; // exact solution

};

// reading parameters from json
Parameters readParameters(const std::string & parFileName){

  Parameters parameters; ///< Parameters object to hold the read parameters

  // Reading JSON file with data
  std::ifstream f(parFileName);
  nlohmann::json parFile = nlohmann::json::parse(f);   

  /// read and set values
  parameters.x0 = parFile["mesh"].value("minX", 0);
  parameters.y0 = parFile["mesh"].value("minY", 0);
  parameters.xn = parFile["mesh"].value("maxX", 1);
  parameters.yn = parFile["mesh"].value("maxY", 1);
  parameters.nx = parFile["mesh"].value("nx", 2);
  parameters.hx = parFile["mesh"].value("hx", 0.5);
  parameters.ny = parFile["mesh"].value("ny", 2);
  parameters.hy = parFile["mesh"].value("hy", 0.5);
  parameters.maxIter = parFile["option"].value("iteration", 1);
  parameters.tol = parFile["option"].value("tolerance", 1e-3);
  parameters.f = parFile["function"]["f"];
  parameters.fBound = parFile["function"]["fbound"];
  parameters.exactSol = parFile["function"]["exact"];

  return parameters; // Return the read parameters
}






int main(int argc, char* argv[]){

  // init mpi
  MPI_Init(&argc,&argv);

  // setup MPI variables
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // reading parameters
  Parameters p = readParameters("../data/param.json");

  // initialize mesh
  Mesh mesh = Mesh::createWithPoints(p.x0,p.xn,p.y0,p.yn,p.nx,p.ny);

  // Start time measurement
  double startTime = MPI_Wtime();

  // parallel solve
  std::vector<double> parallelSol = parallelSolve(mesh,p.f,p.maxIter,p.tol);

  // End time measurement
  double endTime = MPI_Wtime();

  // comapre solutions
  if(rank == 0){

    compareSolution(mesh,parallelSol,p.exactSol,size);

    std::cout << "Elapsed time: " << endTime - startTime << " seconds" << std::endl << std::endl;

  }

  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Finalize();

  return 0;
}

