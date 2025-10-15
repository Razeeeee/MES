#pragma once

#include <functional>
#include <vector>
#include <stdexcept>

class GaussianIntegration {
public:
    // Function types for 1D and 2D integration
    using Function1D = std::function<double(double)>;
    using Function2D = std::function<double(double, double)>;

    // Constructor
    GaussianIntegration();

    // 1D Integration methods
    double integrate1D(const Function1D& func, double a, double b, int numPoints = 2);
    
    // 2D Integration methods  
    double integrate2D(const Function2D& func, double x1, double x2, double y1, double y2, int numPoints = 2);

    // Generic integration method that auto-detects 1D vs 2D
    double integrate(const Function1D& func, const std::vector<double>& bounds, int numPoints = 2);
    double integrate(const Function2D& func, const std::vector<double>& bounds, int numPoints = 2);

private:
    // Gaussian quadrature weights and points for 2-point and 3-point rules
    struct GaussianPoints {
        std::vector<double> points;
        std::vector<double> weights;
    };

    GaussianPoints gauss2Point;
    GaussianPoints gauss3Point;

    // Initialize Gaussian quadrature points and weights
    void initializeGaussianData();

    // Helper method to get appropriate Gaussian data
    const GaussianPoints& getGaussianData(int numPoints) const;

    // Transform coordinates from [-1, 1] to [a, b]
    double transformCoordinate(double xi, double a, double b) const;
};
