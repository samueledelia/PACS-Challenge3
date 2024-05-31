# PACS Course Challenge 3: A Matrix-Free Parallel Solver for the Laplace Equation

## Project Objective

The goal of this project is to develop a parallel solver for the Laplace equation using the discrete Jacobi iterative method.

## Problem Statement

Solve the Laplace equation:

$$ -\Delta u = f(x), \quad \text{in} \; \Omega = (0, 1)^2, $$

with the following boundary conditions:

$$ u = 0, \quad \text{on} \; x = 0, $$
$$ u = 0, \quad \text{on} \; x = 1, $$
$$ u = 0, \quad \text{on} \; y = 0, $$
$$ u = 0, \quad \text{on} \; y = 1, $$

where the solution $ u $ is represented as a dense matrix $ U $ of size $ n \times n $.

The function $ f $ is initially defined as:

$$ f(x) = 8\pi^2 \sin(2\pi x) \sin(2\pi y) $$

However, you are free to modify this function as needed.
