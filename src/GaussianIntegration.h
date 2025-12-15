#pragma once

#include <functional>
#include <vector>
#include <stdexcept>

class GaussianIntegration {
public:
    using Function1D = std::function<double(double)>;
    using Function2D = std::function<double(double, double)>;

    GaussianIntegration();

    double integrate1D(const Function1D& func, double a, double b, int numPoints = 2);
    double integrate2D(const Function2D& func, double x1, double x2, double y1, double y2, int numPoints = 2);
    double integrate(const Function1D& func, const std::vector<double>& bounds, int numPoints = 2);
    double integrate(const Function2D& func, const std::vector<double>& bounds, int numPoints = 2);

private:
    struct GaussianPoints {
        std::vector<double> points;
        std::vector<double> weights;
    };

    GaussianPoints gauss2Point;
    GaussianPoints gauss3Point;
    GaussianPoints gauss4Point;

    void initializeGaussianData();
    const GaussianPoints& getGaussianData(int numPoints) const;
    double transformCoordinate(double xi, double a, double b) const;
};
