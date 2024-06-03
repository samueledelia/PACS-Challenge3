#ifndef POINT_HPP
#define POINT_HPP

class Point{

    private:

        double x;
        double y;

    public:

        // constructor
        explicit Point(double x, double y) : x{x},y{y} {};

        // getters
        inline double getX() const { return x; }
        inline double getY() const { return y; }

};
    
#endif