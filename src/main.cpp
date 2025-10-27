#include "Element.h"
#include "Grid.h"
#include <iostream>
#include <iomanip>
#include <cmath>

void printMatrix(const std::vector<std::vector<double>>& matrix, const std::string& name) {
    std::cout << name << ":\n";
    for (const auto& row : matrix) {
        std::cout << "  [";
        for (size_t j = 0; j < row.size(); ++j) {
            std::cout << std::setw(12) << std::fixed << std::setprecision(6) << row[j];
            if (j < row.size() - 1) std::cout << ", ";
        }
        std::cout << " ]\n";
    }
}

int main() {
    std::cout << std::fixed << std::setprecision(6);
    
    // Load the grid
    Grid grid;
    if (!grid.loadFromFile("grids/Test2_4_4_MixGrid.txt")) {
        std::cerr << "Failed to load grid file!\n";
        return 1;
    }
    
    // Get global data including conductivity
    const auto& globalData = grid.getGlobalData();
    double conductivity = globalData.getConductivity();
    
    std::cout << "=== H Matrix Calculation ===\n\n";
    std::cout << "Conductivity: " << conductivity << "\n\n";
    std::cout << "=== Shape Function Derivatives with Respect to Physical Coordinates (x, y) ===\n\n";
    
    // First show the natural coordinate derivative matrices
    std::cout << "2x2 Gauss Points Integration Scheme:\n";
    std::cout << "Gauss points: ±1/sqrt(3) = ±" << 1.0/std::sqrt(3.0) << "\n\n";
    
    std::vector<std::string> point2Names = {
        "(-1/sqrt(3), -1/sqrt(3))",
        "( 1/sqrt(3), -1/sqrt(3))",
        "(-1/sqrt(3),  1/sqrt(3))",
        "( 1/sqrt(3),  1/sqrt(3))"
    };
    
    // Show dN/dXi and dN/dEta matrices (natural coordinates)
    std::cout << "Shape function derivatives with respect to natural coordinates:\n\n";
    std::cout << "dN/dXi Matrix (4x4):\n";
    auto dN_dXi_2pt = Element::calculateShapeFunctionDerivativeMatrix_Xi(2);
    for (size_t i = 0; i < dN_dXi_2pt.size(); ++i) {
        std::cout << "Point " << (i+1) << " " << point2Names[i] << ": [";
        for (size_t j = 0; j < dN_dXi_2pt[i].size(); ++j) {
            std::cout << std::setw(8) << dN_dXi_2pt[i][j];
            if (j < dN_dXi_2pt[i].size() - 1) std::cout << ", ";
        }
        std::cout << " ]\n";
    }
    
    std::cout << "\ndN/dEta Matrix (4x4):\n";
    auto dN_dEta_2pt = Element::calculateShapeFunctionDerivativeMatrix_Eta(2);
    for (size_t i = 0; i < dN_dEta_2pt.size(); ++i) {
        std::cout << "Point " << (i+1) << " " << point2Names[i] << ": [";
        for (size_t j = 0; j < dN_dEta_2pt[i].size(); ++j) {
            std::cout << std::setw(8) << dN_dEta_2pt[i][j];
            if (j < dN_dEta_2pt[i].size() - 1) std::cout << ", ";
        }
        std::cout << " ]\n";
    }
    
    std::cout << "\n" << std::string(80, '=') << "\n\n";
    
    const auto& nodes = grid.getNodes();
    const auto& elements = grid.getElements();
    
    std::cout << "Physical derivatives (dN/dx, dN/dy) for elements:\n\n";
    
    // Process first 3 elements to demonstrate physical derivatives
    for (int elemIdx = 0; elemIdx < 3 && elemIdx < (int)elements.size(); ++elemIdx) {
        const auto& element = elements[elemIdx];
        
        std::cout << "=== Element " << element.getId() << " ===\n";
        
        // Get node IDs for this element
        const auto& nodeIds = element.getNodeIds();
        std::cout << "Nodes: ";
        for (size_t i = 0; i < nodeIds.size(); ++i) {
            std::cout << nodeIds[i];
            if (i < nodeIds.size() - 1) std::cout << ", ";
        }
        std::cout << "\n";
        
        // Extract node coordinates
        std::vector<double> nodeX(4), nodeY(4);
        for (size_t i = 0; i < nodeIds.size(); ++i) {
            const Node& node = nodes[nodeIds[i] - 1];
            nodeX[i] = node.getX();
            nodeY[i] = node.getY();
        }
        
        std::cout << "Node coordinates:\n";
        for (size_t i = 0; i < 4; ++i) {
            std::cout << "  Node " << nodeIds[i] << ": (" << std::setw(10) << nodeX[i] 
                      << ", " << std::setw(10) << nodeY[i] << ")\n";
        }
        
        try {
            // Calculate physical derivatives using the new method
            auto physicalDerivatives = element.calculatePhysicalDerivativeMatrices(nodeX, nodeY, 2);
            auto dN_dx_matrix = physicalDerivatives.first;
            auto dN_dy_matrix = physicalDerivatives.second;
            
            std::cout << "\nShape function derivatives with respect to physical coordinates:\n\n";
            
            std::cout << "dN/dx Matrix (4x4) - derivatives with respect to x:\n";
            for (size_t i = 0; i < dN_dx_matrix.size(); ++i) {
                std::cout << "Point " << (i+1) << " " << point2Names[i] << ": [";
                for (size_t j = 0; j < dN_dx_matrix[i].size(); ++j) {
                    std::cout << std::setw(10) << dN_dx_matrix[i][j];
                    if (j < dN_dx_matrix[i].size() - 1) std::cout << ", ";
                }
                std::cout << " ]\n";
            }
            
            std::cout << "\ndN/dy Matrix (4x4) - derivatives with respect to y:\n";
            for (size_t i = 0; i < dN_dy_matrix.size(); ++i) {
                std::cout << "Point " << (i+1) << " " << point2Names[i] << ": [";
                for (size_t j = 0; j < dN_dy_matrix[i].size(); ++j) {
                    std::cout << std::setw(10) << dN_dy_matrix[i][j];
                    if (j < dN_dy_matrix[i].size() - 1) std::cout << ", ";
                }
                std::cout << " ]\n";
            }
            
            // Also show the Jacobian data for verification
            std::cout << "\nJacobian verification for this element:\n";
            auto jacobianData = element.calculateElementJacobians(nodeX, nodeY, 2);
            for (size_t i = 0; i < jacobianData.size(); ++i) {
                std::cout << "\nPoint " << (i+1) << " " << point2Names[i] << ":\n";
                std::cout << "  Determinant |J| = " << std::setw(12) << jacobianData[i].detJ << "\n";
                printMatrix(jacobianData[i].inverseJacobian, "  Inverse Jacobian [J^-1]");
            }
            
            // Calculate and display H matrix
            std::cout << "\nH Matrix Calculation:\n";
            std::cout << "H = conductivity * (dN/dx * dN/dx^T + dN/dy * dN/dy^T) * detJ * weight\n";
            std::cout << "Summed over all integration points\n\n";
            
            auto H_matrix = element.calculateHMatrix(nodeX, nodeY, conductivity, 2);
            printMatrix(H_matrix, "H Matrix [4x4]");
            
        } catch (const std::exception& e) {
            std::cerr << "Error calculating physical derivatives: " << e.what() << "\n";
        }
        
        std::cout << "\n" << std::string(80, '-') << "\n\n";
    }
    
    // Summary showing all elements' physical derivatives at first integration point
    std::cout << "=== Summary: Physical Derivatives at First Integration Point for All Elements ===\n\n";
    std::cout << "Point 1: ξ = -1/√3, η = -1/√3\n\n";
    std::cout << "Element | dN1/dx    | dN2/dx    | dN3/dx    | dN4/dx    | dN1/dy    | dN2/dy    | dN3/dy    | dN4/dy\n";
    std::cout << std::string(100, '-') << "\n";
    
    for (const auto& element : elements) {
        const auto& nodeIds = element.getNodeIds();
        
        // Extract node coordinates
        std::vector<double> nodeX(4), nodeY(4);
        for (size_t i = 0; i < nodeIds.size(); ++i) {
            const Node& node = nodes[nodeIds[i] - 1];
            nodeX[i] = node.getX();
            nodeY[i] = node.getY();
        }
        
        try {
            auto physicalDerivatives = element.calculatePhysicalDerivativeMatrices(nodeX, nodeY, 2);
            auto dN_dx_matrix = physicalDerivatives.first;
            auto dN_dy_matrix = physicalDerivatives.second;
            
            std::cout << std::setw(7) << element.getId() << " |";
            // Show dN/dx for all shape functions at first integration point
            for (int j = 0; j < 4; ++j) {
                std::cout << std::setw(10) << dN_dx_matrix[0][j] << " |";
            }
            // Show dN/dy for all shape functions at first integration point
            for (int j = 0; j < 4; ++j) {
                std::cout << std::setw(10) << dN_dy_matrix[0][j] << " |";
            }
            std::cout << "\n";
            
        } catch (const std::exception& e) {
            std::cout << std::setw(7) << element.getId() << " | Error: " << e.what() << "\n";
        }
    }
    
    // H Matrix Summary for all elements
    std::cout << "\n" << std::string(80, '=') << "\n";
    std::cout << "\n=== H Matrix Summary for All Elements ===\n\n";
    std::cout << "H Matrix diagonal elements (H[0][0], H[1][1], H[2][2], H[3][3]) for each element:\n\n";
    std::cout << "Element | H[0][0]     | H[1][1]     | H[2][2]     | H[3][3]     | Matrix Sum\n";
    std::cout << std::string(80, '-') << "\n";
    
    for (const auto& element : elements) {
        const auto& nodeIds = element.getNodeIds();
        
        // Extract node coordinates
        std::vector<double> nodeX(4), nodeY(4);
        for (size_t i = 0; i < nodeIds.size(); ++i) {
            const Node& node = nodes[nodeIds[i] - 1];
            nodeX[i] = node.getX();
            nodeY[i] = node.getY();
        }
        
        try {
            auto H_matrix = element.calculateHMatrix(nodeX, nodeY, conductivity, 2);
            
            // Calculate sum of all matrix elements
            double matrixSum = 0.0;
            for (const auto& row : H_matrix) {
                for (double val : row) {
                    matrixSum += val;
                }
            }
            
            std::cout << std::setw(7) << element.getId() << " |";
            for (int i = 0; i < 4; ++i) {
                std::cout << std::setw(12) << H_matrix[i][i] << " |";
            }
            std::cout << std::setw(12) << matrixSum << "\n";
            
        } catch (const std::exception& e) {
            std::cout << std::setw(7) << element.getId() << " | Error: " << e.what() << "\n";
        }
    }
    
    return 0;
}
