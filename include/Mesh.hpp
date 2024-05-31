#ifndef MESH_HPP
#define MESH_HPP

#include "Point.hpp"

#include <iostream>
#include <iomanip>


    /**
     * @brief Grid mesh class on rectangular domain
    */

    class Mesh {

    private:
        double x0{0};
        double xn{0};
        double y0{0};
        double yn{0};

        std::vector<Point> grid; // coordinates of the Points stored left to right, bottom to top 
        size_t nx; // number of points in x direction
        size_t ny; // number of points in y direction
        double hx; // spacing along x
        double hy; // spacing along y
    public:

        // A series of different Constructor
        Mesh (double x0_, double xn_, double y0_, double yn_, size_t nx_, size_t ny_);
        Mesh (double x0_, double xn_, double y0_, double yn_, double hx_, double hy_);
        // Constructors given top right and bottom left points and others parameters
        Mesh2D (const Point & tr, const Point & bl, size_t nx_, size_t ny_);
        Mesh2D (const Point & tr, const Point & bl, double hx_, double hy_);


        // getters
        inline double getx0() const { return x0; }; 
        inline double gety0() const { return y0; }; 
        inline double getxn() const { return xn; }; 
        inline double getyn() const { return yn; }; 

        inline size_t getnx() const { return nx; }
        inline size_t getny() const { return ny; }
        inline double gethx() const { return hx; }
        inline double gethy() const { return hy; }

        // Print domain four vertices 
        friend std::ostream& operator<<(std::ostream& os, const Domain2D & domain);




        /**
         * @brief Create the grid given ranges and number of points
        */
        static Mesh createWithPoints(double x0_, double xn_, double y0_, double yn_, size_t nx, size_t ny);

        /**
         * @brief factory given ranges and spacing
        */
        static Mesh createWithSpacing(double x0_, double xn_, double y0_, double yn_, double hx, double hy);

        // ***** factory methods to select explicitly the strategy of the creation of the mesh (vertices) ***** //

        /**
         * @brief factory given vertices and number of points
        */
        static Mesh createWithPoints(const Point & tr, const Point & bl, size_t nx, size_t ny);

        /**
         * @brief factory given vertices and spacing
        */
        static Mesh createWithSpacing(const Point & tr, const Point & bl, double hx, double hy);

        /**
         * @brief point by its indexes in the mesh
        */
        const Point & operator()(size_t i, size_t j) const;


    private:
        // Check that the entries are valid
        inline bool check_entries(double x0_, double xn_, double y0_, double yn_) const { return x0_ < xn_ && y0_ < yn_; }; 

    };


Mesh::Mesh(double x0_, double xn_, double y0_, double yn_, size_t nx_, size_t ny_) : nx{nx_},ny{ny_} { 

        if(isValid(x0_, xn_, y0_, yn_)){ 
            x0 = x0_;
            xn = xn_;
            y0 = y0_;
            yn = yn_;
        }
        else{
            throw std::invalid_argument("Not a proper rectangular domain!");
        } 

        hx = (xn-x0)/(nx-1); 
        hy = (yn-y0)/(ny-1);

        coordinate.reserve(nx*ny); // Preallocate space to get a faster code!

        // populate points grid
        for(double y = y0 ; y <= yn ; y+=hy){
            for(double x = x0 ; x <= xn ; x+=hx){
                coordinate.emplace_back(x,y);
            }
        }

    };

    Mesh::Mesh(double x0_, double xn_, double y0_, double yn_, double hx_, double hy_) : 
                Mesh2D(x0_, xn_, y0_, yn_, static_cast<size_t>(std::ceil((xn - x0) / hx_)+1), static_cast<size_t>(std::ceil((yn - y0) / hy_)+1)) {
        
        if (std::fmod((xn - x0),hx_) >= std::numeric_limits<double>::epsilon() ){
            std::cout << "Cannot evenly divide with given x spacing, correction occurred on spacing from " << hx_ << "---> " << hx << std::endl;
        }

        if (std::fmod((yn - y0),hy_) >= std::numeric_limits<double>::epsilon() ){ ///< if no possible evenly divide along y with given spacing
            std::cout << "Cannot evenly divide with given y spacing, corrrection occured on spacing from " << hy_ << "---> " << hy << std::endl;
        }

    };

    Mesh::Mesh(const Point & tr, const Point & bl, size_t nx_, size_t ny_) : 
                Mesh(bl.getX(), tr.getX(), bl.getY(), tr.getY(), nx_, ny) { 
    };

    Mesh::Mesh(const Point & tr, const Point & bl, double hx_, double hy_) : 
                Mesh(bl.getX(), tr.getX(), bl.getY(), tr.getY(), hx_, hy_) {
    };

std::ostream& operator<<(std::ostream& os, const Mesh & mesh) {

        os << std::fixed <<std::setprecision(2) << static_cast<const Mesh &>(mesh) 
            << "Number of points in x direction: " << mesh.getnx() << " , Spacing: " << mesh.gethx() << std::endl
            << "Number of points in y direction: " << mesh.getny() << " , Spacing: " << mesh.gethy() << std::endl << std::endl;

        std::cout << std::setw(3);

        for(int i = mesh.ny-1 ; i >= 0 ; --i){
            for(int j = 0 ; j < mesh.nx ; ++j){
                std::cout << "(" << mesh(i,j).getX() << "," <<  mesh(i,j).getY() << ") " << std::setw(3);
            }
            std::cout << std::endl;
        }

        std::cout << std::endl;

        return os;
        
    };




    Mesh Mesh::createWithPoints(double x0_, double xn_, double y0_, double yn_, size_t nx, size_t ny) {
        return Mesh(x0_, xn_, y0_, yn_, nx, ny);
    };

    Mesh Mesh::createWithSpacing(double x0_, double xn_, double y0_, double yn_, double hx, double hy) {
        return Mesh(x0_, xn_, y0_, yn_, hx, hy);
    };

    Mesh Mesh::createWithPoints(const Point & tr, const Point & bl, size_t nx, size_t ny) {
        return Mesh(tr, bl, nx, ny);
    };

    Mesh Mesh::createWithSpacing(const Point & tr, const Point & bl, double hx, double hy) {
        return Mesh(tr, bl, hx, hy);
    };



#endif