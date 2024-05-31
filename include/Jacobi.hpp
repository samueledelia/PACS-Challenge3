#ifndef HH_JACOBI_HH
#define HH_JACOBI_HH
#include <map>
#include <array>
#include <functional>
namespace JacobiMethod{

template <class T>
    using ElemType = std::map<std::array<std::size_t,2>,T>;

template <class T>
class Jacobi {
    private: 
        std::function<T(T,T)> m_fun;
        int m_nx, m_ny;
        T m_h;
        T m_tol;
        int m_maxIter;
        ElemType<T> m_U;
    public:
        Jacobi(const std::function<T(T,T)>& fun, int nx, int ny,T h,T tol, int maxIter);
        void solve(ElemType<T>& sol, unsigned int num_threads = 1);
};

template <class T>
Jacobi<T>::Jacobi(const std::function<T(T,T)>& fun, int nx, int ny, T h,T tol, int maxIter): 
    m_fun(fun), 
    m_nx(nx),
    m_ny(ny), 
    m_h(h),
    m_tol(tol),
    m_maxIter(maxIter){
    //initialize the solution with zeros
    //account for the boundary conditions of homogeous Dirichlet type
    for (size_t i = 0; i < unsigned(nx+1); i++) {
        for (size_t j = 0; j < unsigned(ny+1); j++) {
            m_U[{i,j}] = 0.0;
        }
    }
}

template <class T>
void Jacobi<T>::solve(ElemType<T> &sol, unsigned int num_threads){
    int iter = 0;  //set the iteration counter to zero
    T error = 1.0; //initialize the error to a value greater than the tolerance
    ElemType<T> U_old = m_U; //initialize the old solution to the initial guess
    while (iter < m_maxIter && error > m_tol){
        error = 0.0; //reset the error
        #pragma omp parallel for num_threads(num_threads) shared(m_nx,m_ny, iter, U_old)   \
        reduction(+:error)
        for (size_t i = 1; i < unsigned(m_nx); i++) {
            for (size_t j = 1; j < unsigned(m_ny); j++) {
                //update the solution with the Jacobi iteration
                m_U[{i,j}] = 0.25*(U_old[{i-1,j}] + U_old[{i+1,j}] + U_old[{i,j-1}] + U_old[{i,j+1}] + m_h*m_h*m_fun(i*m_h,j*m_h));
                //compute the error
                error +=(m_U[{i,j}]-U_old[{i,j}]) * (m_U[{i,j}]-U_old[{i,j}]);
                }
        }
        error= std::sqrt(m_h*error);//compute the norm of the error
        U_old = m_U;    //update the old solution
        iter++;//increment the iteration counter
        
    }
    sol = m_U;//return the solution
}
}//namespace JacobiMethod
#endif//HH_JACOBI_HH