#include "Element.h"
#include <iostream>
#include <iomanip>
#include <cmath>

int main() {
    std::cout << std::fixed << std::setprecision(6);
    
    std::cout << "=== Shape Function Derivative Matrices for Jacobian ===\n\n";
    
    // 2x2 Gauss points (4 points total)
    std::cout << "2x2 Gauss Points - dN/dXi Matrix (4x4):\n";
    std::cout << "Integration points: ±1/sqrt(3) = ±" << 1.0/std::sqrt(3.0) << "\n";
    
    std::vector<std::string> point2Names = {
        "(-1/sqrt(3), -1/sqrt(3))",
        "( 1/sqrt(3), -1/sqrt(3))",
        "(-1/sqrt(3),  1/sqrt(3))",
        "( 1/sqrt(3),  1/sqrt(3))"
    };
    
    auto dN_dXi_2pt = Element::calculateShapeFunctionDerivativeMatrix_Xi(2);
    for (size_t i = 0; i < dN_dXi_2pt.size(); ++i) {
        std::cout << "Point " << i+1 << " " << point2Names[i] << ": [";
        for (size_t j = 0; j < dN_dXi_2pt[i].size(); ++j) {
            std::cout << std::setw(8) << dN_dXi_2pt[i][j];
            if (j < dN_dXi_2pt[i].size() - 1) std::cout << ", ";
        }
        std::cout << " ]\n";
    }
    
    std::cout << "\n2x2 Gauss Points - dN/dEta Matrix (4x4):\n";
    auto dN_dEta_2pt = Element::calculateShapeFunctionDerivativeMatrix_Eta(2);
    for (size_t i = 0; i < dN_dEta_2pt.size(); ++i) {
        std::cout << "Point " << i+1 << " " << point2Names[i] << ": [";
        for (size_t j = 0; j < dN_dEta_2pt[i].size(); ++j) {
            std::cout << std::setw(8) << dN_dEta_2pt[i][j];
            if (j < dN_dEta_2pt[i].size() - 1) std::cout << ", ";
        }
        std::cout << " ]\n";
    }
    
    // 3x3 Gauss points (9 points total)
    std::cout << "\n3x3 Gauss Points - dN/dXi Matrix (9x4):\n";
    std::cout << "Integration points: -sqrt(3/5) = " << -std::sqrt(3.0/5.0) << ", 0, +sqrt(3/5) = " << std::sqrt(3.0/5.0) << "\n";
    
    std::vector<std::string> point3Names = {
        "(-sqrt(3/5), -sqrt(3/5))",
        "( 0,         -sqrt(3/5))",
        "(+sqrt(3/5), -sqrt(3/5))",
        "(-sqrt(3/5),  0        )",
        "( 0,          0        )",
        "(+sqrt(3/5),  0        )",
        "(-sqrt(3/5), +sqrt(3/5))",
        "( 0,         +sqrt(3/5))",
        "(+sqrt(3/5), +sqrt(3/5))"
    };
    
    auto dN_dXi_3pt = Element::calculateShapeFunctionDerivativeMatrix_Xi(3);
    for (size_t i = 0; i < dN_dXi_3pt.size(); ++i) {
        std::cout << "Point " << std::setw(2) << i+1 << " " << point3Names[i] << ": [";
        for (size_t j = 0; j < dN_dXi_3pt[i].size(); ++j) {
            std::cout << std::setw(8) << dN_dXi_3pt[i][j];
            if (j < dN_dXi_3pt[i].size() - 1) std::cout << ", ";
        }
        std::cout << " ]\n";
    }
    
    std::cout << "\n3x3 Gauss Points - dN/dEta Matrix (9x4):\n";
    auto dN_dEta_3pt = Element::calculateShapeFunctionDerivativeMatrix_Eta(3);
    for (size_t i = 0; i < dN_dEta_3pt.size(); ++i) {
        std::cout << "Point " << std::setw(2) << i+1 << " " << point3Names[i] << ": [";
        for (size_t j = 0; j < dN_dEta_3pt[i].size(); ++j) {
            std::cout << std::setw(8) << dN_dEta_3pt[i][j];
            if (j < dN_dEta_3pt[i].size() - 1) std::cout << ", ";
        }
        std::cout << " ]\n";
    }
    
    return 0;
}
