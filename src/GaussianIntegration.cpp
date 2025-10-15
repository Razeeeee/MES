#include "GaussianIntegration.h"
#include <cmath>

GaussianIntegration::GaussianIntegration() {
    initializeGaussianData();
}

void GaussianIntegration::initializeGaussianData() {
    // 2-point Gaussian quadrature
    gauss2Point.points = {-1.0/std::sqrt(3.0), 1.0/std::sqrt(3.0)};
    gauss2Point.weights = {1.0, 1.0};
    
    // 3-point Gaussian quadrature
    gauss3Point.points = {-std::sqrt(3.0/5.0), 0.0, std::sqrt(3.0/5.0)};
    gauss3Point.weights = {5.0/9.0, 8.0/9.0, 5.0/9.0};
}

const GaussianIntegration::GaussianPoints& GaussianIntegration::getGaussianData(int numPoints) const {
    if (numPoints == 2) {
        return gauss2Point;
    } else if (numPoints == 3) {
        return gauss3Point;
    } else {
        throw std::invalid_argument("Only 2-point and 3-point Gaussian quadrature are supported");
    }
}

double GaussianIntegration::transformCoordinate(double xi, double a, double b) const {
    return 0.5 * (b - a) * xi + 0.5 * (a + b);
}

double GaussianIntegration::integrate1D(const Function1D& func, double a, double b, int numPoints) {
    const auto& gaussData = getGaussianData(numPoints);
    double result = 0.0;
    double jacobian = (b - a) / 2.0;  // Simple transformation jacobian
    
    for (size_t i = 0; i < gaussData.points.size(); ++i) {
        double x = transformCoordinate(gaussData.points[i], a, b);
        result += gaussData.weights[i] * func(x);
    }
    
    return result * jacobian;
}

double GaussianIntegration::integrate2D(const Function2D& func, double x1, double x2, double y1, double y2, int numPoints) {
    const auto& gaussData = getGaussianData(numPoints);
    double result = 0.0;
    double jacobian = ((x2 - x1) / 2.0) * ((y2 - y1) / 2.0);  // 2D transformation jacobian
    
    for (size_t i = 0; i < gaussData.points.size(); ++i) {
        for (size_t j = 0; j < gaussData.points.size(); ++j) {
            double x = transformCoordinate(gaussData.points[i], x1, x2);
            double y = transformCoordinate(gaussData.points[j], y1, y2);
            result += gaussData.weights[i] * gaussData.weights[j] * func(x, y);
        }
    }
    
    return result * jacobian;
}

double GaussianIntegration::integrate(const Function1D& func, const std::vector<double>& bounds, int numPoints) {
    if (bounds.size() != 2) {
        throw std::invalid_argument("1D integration requires exactly 2 boundaries [a, b]");
    }
    return integrate1D(func, bounds[0], bounds[1], numPoints);
}

double GaussianIntegration::integrate(const Function2D& func, const std::vector<double>& bounds, int numPoints) {
    if (bounds.size() != 4) {
        throw std::invalid_argument("2D integration requires exactly 4 boundaries [x1, x2, y1, y2]");
    }
    return integrate2D(func, bounds[0], bounds[1], bounds[2], bounds[3], numPoints);
}
