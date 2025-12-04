#pragma once

#include <functional>
#include <vector>
#include <stdexcept>

/**
 * @brief Gaussian quadrature integration for FEM
 * 
 * Implements numerical integration using Gaussian quadrature with 2, 3, or 4 points.
 * Supports both 1D (edge) and 2D (volume) integration needed for FEM calculations.
 * 
 * Theory:
 * Gaussian quadrature approximates integrals as weighted sums:
 * ∫f(x)dx ≈ Σ(w_i * f(x_i)) where w_i are weights and x_i are integration points
 * 
 * For n-point quadrature:
 * - 2-point: Exact for polynomials up to degree 3
 * - 3-point: Exact for polynomials up to degree 5
 * - 4-point: Exact for polynomials up to degree 7
 */
class GaussianIntegration {
public:
    // Function types for 1D and 2D integration
    using Function1D = std::function<double(double)>;
    using Function2D = std::function<double(double, double)>;

    // Constructor
    GaussianIntegration();

    // 1D Integration methods (for boundary/edge integrals)
    double integrate1D(const Function1D& func, double a, double b, int numPoints = 2);
    
    // 2D Integration methods (for volume integrals)
    double integrate2D(const Function2D& func, double x1, double x2, double y1, double y2, int numPoints = 2);

    // Generic integration method that auto-detects 1D vs 2D
    double integrate(const Function1D& func, const std::vector<double>& bounds, int numPoints = 2);
    double integrate(const Function2D& func, const std::vector<double>& bounds, int numPoints = 2);

private:
    // Gaussian quadrature weights and points for 2, 3, and 4-point rules
    struct GaussianPoints {
        std::vector<double> points;   // Integration points in [-1, 1]
        std::vector<double> weights;  // Corresponding weights
    };

    GaussianPoints gauss2Point;
    GaussianPoints gauss3Point;
    GaussianPoints gauss4Point;

    // Initialize Gaussian quadrature points and weights
    void initializeGaussianData();

    // Helper method to get appropriate Gaussian data
    const GaussianPoints& getGaussianData(int numPoints) const;

    // Transform coordinates from [-1, 1] to [a, b]
    double transformCoordinate(double xi, double a, double b) const;
};
