#include "GaussianIntegration.h"
#include <iostream>
#include <iomanip>

int main() {
    GaussianIntegration gauss;
    
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Gaussian Integration Tests\n";
    std::cout << "==========================\n\n";
    
    // Define 1D function: f(x) = 5x^2 + 3x + 6
    auto f1D = [](double x) {
        return 5.0 * x * x + 3.0 * x + 6.0;
    };
    
    // Define 2D function: f(x,y) = 5x^2*y^2 + 3xy + 6
    auto f2D = [](double x, double y) {
        return 5.0 * x * x * y * y + 3.0 * x * y + 6.0;
    };
    
    // 1D Integration tests (boundaries: -1 to 1)
    std::cout << "1D Function: f(x) = 5x^2 + 3x + 6, bounds: [-1, 1]\n";
    double result1D_2pt = gauss.integrate1D(f1D, -1.0, 1.0, 2);
    double result1D_3pt = gauss.integrate1D(f1D, -1.0, 1.0, 3);
    
    std::cout << "2-point Gaussian: " << result1D_2pt << std::endl;
    std::cout << "3-point Gaussian: " << result1D_3pt << std::endl;
    std::cout << std::endl;
    
    // 2D Integration tests (boundaries: x: -1 to 1, y: -1 to 1)
    std::cout << "2D Function: f(x,y) = 5x^2*y^2 + 3xy + 6, bounds: x[-1,1], y[-1,1]\n";
    double result2D_2pt = gauss.integrate2D(f2D, -1.0, 1.0, -1.0, 1.0, 2);
    double result2D_3pt = gauss.integrate2D(f2D, -1.0, 1.0, -1.0, 1.0, 3);
    
    std::cout << "2-point Gaussian: " << result2D_2pt << std::endl;
    std::cout << "3-point Gaussian: " << result2D_3pt << std::endl;
    
    return 0;
}
