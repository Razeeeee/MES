#include "Element.h"
#include <iomanip>
#include <cmath>

Element::Element() : id(0), type("DC2D4") {}

Element::Element(int id, const std::vector<int>& nodeIds, const std::string& type) 
    : id(id), nodeIds(nodeIds), type(type) {}

void Element::print() const {
    std::cout << "Element " << std::setw(3) << id << "]: ";
    for (size_t i = 0; i < nodeIds.size(); ++i) {
        std::cout << nodeIds[i];
        if (i < nodeIds.size() - 1) {
            std::cout << ", ";
        }
    }
    std::cout << std::endl;
}

// Shape function derivatives with respect to Xi
// For 4-node quadrilateral element: N1 = 0.25*(1-xi)*(1-eta), N2 = 0.25*(1+xi)*(1-eta)
//                                   N3 = 0.25*(1+xi)*(1+eta), N4 = 0.25*(1-xi)*(1+eta)

double Element::dN1_dXi(double /*xi*/, double eta) {
    return -0.25 * (1.0 - eta);
}

double Element::dN2_dXi(double /*xi*/, double eta) {
    return 0.25 * (1.0 - eta);
}

double Element::dN3_dXi(double /*xi*/, double eta) {
    return 0.25 * (1.0 + eta);
}

double Element::dN4_dXi(double /*xi*/, double eta) {
    return -0.25 * (1.0 + eta);
}

// Shape function derivatives with respect to Eta

double Element::dN1_dEta(double xi, double /*eta*/) {
    return -0.25 * (1.0 - xi);
}

double Element::dN2_dEta(double xi, double /*eta*/) {
    return -0.25 * (1.0 + xi);
}

double Element::dN3_dEta(double xi, double /*eta*/) {
    return 0.25 * (1.0 + xi);
}

double Element::dN4_dEta(double xi, double /*eta*/) {
    return 0.25 * (1.0 - xi);
}

// Calculate matrix of shape function derivatives with respect to Xi
std::vector<std::vector<double>> Element::calculateShapeFunctionDerivativeMatrix_Xi(int numPoints) {
    std::vector<double> gaussPoints;
    
    // Set up Gaussian points based on numPoints
    if (numPoints == 2) {
        gaussPoints = {-1.0/std::sqrt(3.0), 1.0/std::sqrt(3.0)};
    } else if (numPoints == 3) {
        gaussPoints = {-std::sqrt(3.0/5.0), 0.0, std::sqrt(3.0/5.0)};
    } else {
        throw std::invalid_argument("Only 2-point and 3-point Gaussian quadrature are supported");
    }
    
    // Create matrix: rows = integration points, columns = shape functions (4 for DC2D4)
    int totalPoints = numPoints * numPoints;
    std::vector<std::vector<double>> matrix(totalPoints, std::vector<double>(4));
    
    int pointIndex = 0;
    for (int j = 0; j < numPoints; ++j) {        // eta (j) is outer loop
        for (int i = 0; i < numPoints; ++i) {    // xi (i) is inner loop
            double xi = gaussPoints[i];
            double eta = gaussPoints[j];
            
            matrix[pointIndex][0] = dN1_dXi(xi, eta);
            matrix[pointIndex][1] = dN2_dXi(xi, eta);
            matrix[pointIndex][2] = dN3_dXi(xi, eta);
            matrix[pointIndex][3] = dN4_dXi(xi, eta);
            
            pointIndex++;
        }
    }
    
    return matrix;
}

// Calculate matrix of shape function derivatives with respect to Eta
std::vector<std::vector<double>> Element::calculateShapeFunctionDerivativeMatrix_Eta(int numPoints) {
    std::vector<double> gaussPoints;
    
    // Set up Gaussian points based on numPoints
    if (numPoints == 2) {
        gaussPoints = {-1.0/std::sqrt(3.0), 1.0/std::sqrt(3.0)};
    } else if (numPoints == 3) {
        gaussPoints = {-std::sqrt(3.0/5.0), 0.0, std::sqrt(3.0/5.0)};
    } else {
        throw std::invalid_argument("Only 2-point and 3-point Gaussian quadrature are supported");
    }
    
    // Create matrix: rows = integration points, columns = shape functions (4 for DC2D4)
    int totalPoints = numPoints * numPoints;
    std::vector<std::vector<double>> matrix(totalPoints, std::vector<double>(4));
    
    int pointIndex = 0;
    for (int j = 0; j < numPoints; ++j) {        // eta (j) is outer loop
        for (int i = 0; i < numPoints; ++i) {    // xi (i) is inner loop
            double xi = gaussPoints[i];
            double eta = gaussPoints[j];
            
            matrix[pointIndex][0] = dN1_dEta(xi, eta);
            matrix[pointIndex][1] = dN2_dEta(xi, eta);
            matrix[pointIndex][2] = dN3_dEta(xi, eta);
            matrix[pointIndex][3] = dN4_dEta(xi, eta);
            
            pointIndex++;
        }
    }
    
    return matrix;
}

// Calculate Jacobian matrix (2x2) at a specific integration point
std::vector<std::vector<double>> Element::calculateJacobian(const std::vector<double>& nodeX, 
                                                           const std::vector<double>& nodeY, 
                                                           double xi, double eta) const {
    // Jacobian matrix is 2x2: [dx/dxi  dy/dxi ]
    //                        [dx/deta dy/deta]
    std::vector<std::vector<double>> jacobian(2, std::vector<double>(2, 0.0));
    
    // dx/dxi = sum(dNi/dxi * xi)
    jacobian[0][0] = dN1_dXi(xi, eta) * nodeX[0] + dN2_dXi(xi, eta) * nodeX[1] + 
                     dN3_dXi(xi, eta) * nodeX[2] + dN4_dXi(xi, eta) * nodeX[3];
    
    // dy/dxi = sum(dNi/dxi * yi)
    jacobian[0][1] = dN1_dXi(xi, eta) * nodeY[0] + dN2_dXi(xi, eta) * nodeY[1] + 
                     dN3_dXi(xi, eta) * nodeY[2] + dN4_dXi(xi, eta) * nodeY[3];
    
    // dx/deta = sum(dNi/deta * xi)
    jacobian[1][0] = dN1_dEta(xi, eta) * nodeX[0] + dN2_dEta(xi, eta) * nodeX[1] + 
                     dN3_dEta(xi, eta) * nodeX[2] + dN4_dEta(xi, eta) * nodeX[3];
    
    // dy/deta = sum(dNi/deta * yi)
    jacobian[1][1] = dN1_dEta(xi, eta) * nodeY[0] + dN2_dEta(xi, eta) * nodeY[1] + 
                     dN3_dEta(xi, eta) * nodeY[2] + dN4_dEta(xi, eta) * nodeY[3];
    
    return jacobian;
}

// Calculate inverse Jacobian matrix (2x2)
std::vector<std::vector<double>> Element::calculateInverseJacobian(const std::vector<std::vector<double>>& jacobian) const {
    std::vector<std::vector<double>> invJacobian(2, std::vector<double>(2, 0.0));
    
    // Calculate determinant
    double detJ = jacobian[0][0] * jacobian[1][1] - jacobian[0][1] * jacobian[1][0];
    
    if (std::abs(detJ) < 1e-12) {
        throw std::runtime_error("Jacobian determinant is too small (singular matrix)");
    }
    
    // Calculate inverse using the formula for 2x2 matrix inverse
    // J^-1 = (1/detJ) * [ J22  -J12]
    //                   [-J21   J11]
    double invDetJ = 1.0 / detJ;
    
    invJacobian[0][0] =  jacobian[1][1] * invDetJ;  // J22 * (1/detJ)
    invJacobian[0][1] = -jacobian[0][1] * invDetJ;  // -J12 * (1/detJ)
    invJacobian[1][0] = -jacobian[1][0] * invDetJ;  // -J21 * (1/detJ)
    invJacobian[1][1] =  jacobian[0][0] * invDetJ;  // J11 * (1/detJ)
    
    return invJacobian;
}

// Calculate determinant of Jacobian matrix
double Element::calculateDetJ(const std::vector<std::vector<double>>& jacobian) const {
    return jacobian[0][0] * jacobian[1][1] - jacobian[0][1] * jacobian[1][0];
}

// Calculate Jacobians for all integration points for this element
std::vector<Element::JacobianData> Element::calculateElementJacobians(const std::vector<double>& nodeX, 
                                                                      const std::vector<double>& nodeY, 
                                                                      int numGaussPoints) const {
    std::vector<double> gaussPoints;
    
    // Set up Gaussian points based on numGaussPoints
    if (numGaussPoints == 2) {
        gaussPoints = {-1.0/std::sqrt(3.0), 1.0/std::sqrt(3.0)};
    } else if (numGaussPoints == 3) {
        gaussPoints = {-std::sqrt(3.0/5.0), 0.0, std::sqrt(3.0/5.0)};
    } else {
        throw std::invalid_argument("Only 2-point and 3-point Gaussian quadrature are supported");
    }
    
    int totalPoints = numGaussPoints * numGaussPoints;
    std::vector<JacobianData> jacobianData(totalPoints);
    
    int pointIndex = 0;
    for (int j = 0; j < numGaussPoints; ++j) {        // eta (j) is outer loop
        for (int i = 0; i < numGaussPoints; ++i) {    // xi (i) is inner loop
            double xi = gaussPoints[i];
            double eta = gaussPoints[j];
            
            // Calculate Jacobian matrix
            jacobianData[pointIndex].jacobian = calculateJacobian(nodeX, nodeY, xi, eta);
            
            // Calculate determinant
            jacobianData[pointIndex].detJ = calculateDetJ(jacobianData[pointIndex].jacobian);
            
            // Calculate inverse Jacobian
            jacobianData[pointIndex].inverseJacobian = calculateInverseJacobian(jacobianData[pointIndex].jacobian);
            
            pointIndex++;
        }
    }
    
    return jacobianData;
}

// Calculate dN/dx for all shape functions at a specific integration point
std::vector<double> Element::calculateDN_Dx(double xi, double eta, 
                                            const std::vector<std::vector<double>>& inverseJacobian) {
    std::vector<double> dN_dx(4);
    
    // Get derivatives with respect to natural coordinates
    std::vector<double> dN_dxi = {dN1_dXi(xi, eta), dN2_dXi(xi, eta), dN3_dXi(xi, eta), dN4_dXi(xi, eta)};
    std::vector<double> dN_deta = {dN1_dEta(xi, eta), dN2_dEta(xi, eta), dN3_dEta(xi, eta), dN4_dEta(xi, eta)};
    
    // Transform using inverse Jacobian: dN/dx = J^-1[0][0] * dN/dxi + J^-1[0][1] * dN/deta
    for (int i = 0; i < 4; ++i) {
        dN_dx[i] = inverseJacobian[0][0] * dN_dxi[i] + inverseJacobian[0][1] * dN_deta[i];
    }
    
    return dN_dx;
}

// Calculate dN/dy for all shape functions at a specific integration point
std::vector<double> Element::calculateDN_Dy(double xi, double eta, 
                                            const std::vector<std::vector<double>>& inverseJacobian) {
    std::vector<double> dN_dy(4);
    
    // Get derivatives with respect to natural coordinates
    std::vector<double> dN_dxi = {dN1_dXi(xi, eta), dN2_dXi(xi, eta), dN3_dXi(xi, eta), dN4_dXi(xi, eta)};
    std::vector<double> dN_deta = {dN1_dEta(xi, eta), dN2_dEta(xi, eta), dN3_dEta(xi, eta), dN4_dEta(xi, eta)};
    
    // Transform using inverse Jacobian: dN/dy = J^-1[1][0] * dN/dxi + J^-1[1][1] * dN/deta
    for (int i = 0; i < 4; ++i) {
        dN_dy[i] = inverseJacobian[1][0] * dN_dxi[i] + inverseJacobian[1][1] * dN_deta[i];
    }
    
    return dN_dy;
}

// Calculate matrices of shape function derivatives with respect to x and y
// at all integration points for this element
std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> 
Element::calculatePhysicalDerivativeMatrices(const std::vector<double>& nodeX, 
                                             const std::vector<double>& nodeY, 
                                             int numGaussPoints) const {
    std::vector<double> gaussPoints;
    
    // Set up Gaussian points based on numGaussPoints
    if (numGaussPoints == 2) {
        gaussPoints = {-1.0/std::sqrt(3.0), 1.0/std::sqrt(3.0)};
    } else if (numGaussPoints == 3) {
        gaussPoints = {-std::sqrt(3.0/5.0), 0.0, std::sqrt(3.0/5.0)};
    } else {
        throw std::invalid_argument("Only 2-point and 3-point Gaussian quadrature are supported");
    }
    
    int totalPoints = numGaussPoints * numGaussPoints;
    
    // Initialize matrices: rows = integration points, columns = shape functions (4 for DC2D4)
    std::vector<std::vector<double>> dN_dx_matrix(totalPoints, std::vector<double>(4));
    std::vector<std::vector<double>> dN_dy_matrix(totalPoints, std::vector<double>(4));
    
    int pointIndex = 0;
    for (int j = 0; j < numGaussPoints; ++j) {        // eta (j) is outer loop
        for (int i = 0; i < numGaussPoints; ++i) {    // xi (i) is inner loop
            double xi = gaussPoints[i];
            double eta = gaussPoints[j];
            
            // Calculate Jacobian and its inverse for this integration point
            std::vector<std::vector<double>> jacobian = calculateJacobian(nodeX, nodeY, xi, eta);
            std::vector<std::vector<double>> inverseJacobian = calculateInverseJacobian(jacobian);
            
            // Calculate physical derivatives for all shape functions at this point
            std::vector<double> dN_dx = calculateDN_Dx(xi, eta, inverseJacobian);
            std::vector<double> dN_dy = calculateDN_Dy(xi, eta, inverseJacobian);
            
            // Store in matrices
            for (int k = 0; k < 4; ++k) {
                dN_dx_matrix[pointIndex][k] = dN_dx[k];
                dN_dy_matrix[pointIndex][k] = dN_dy[k];
            }
            
            pointIndex++;
        }
    }
    
    return std::make_pair(dN_dx_matrix, dN_dy_matrix);
}

std::vector<std::vector<double>> Element::calculateHMatrix(const std::vector<double>& nodeX, 
                                                          const std::vector<double>& nodeY, 
                                                          double conductivity,
                                                          int numGaussPoints) const {
    // Initialize 4x4 H matrix with zeros
    std::vector<std::vector<double>> H(4, std::vector<double>(4, 0.0));
    
    // Get physical derivative matrices and Jacobian data
    auto physicalDerivatives = calculatePhysicalDerivativeMatrices(nodeX, nodeY, numGaussPoints);
    auto dN_dx_matrix = physicalDerivatives.first;  // [integration_points x 4]
    auto dN_dy_matrix = physicalDerivatives.second; // [integration_points x 4]
    auto jacobianData = calculateElementJacobians(nodeX, nodeY, numGaussPoints);
    
    // Get Gauss weights for integration
    std::vector<double> weights;
    if (numGaussPoints == 2) {
        weights = {1.0, 1.0, 1.0, 1.0}; // 2x2 Gauss quadrature weights
    } else if (numGaussPoints == 3) {
        // 3x3 Gauss quadrature weights
        double w5_9 = 5.0/9.0;
        double w8_9 = 8.0/9.0;
        weights = {w5_9*w5_9, w5_9*w8_9, w5_9*w5_9,
                   w8_9*w5_9, w8_9*w8_9, w8_9*w5_9,
                   w5_9*w5_9, w5_9*w8_9, w5_9*w5_9};
    }
    
    // Calculate H matrix by summing contributions from all integration points
    // H = sum over all integration points of: conductivity * (dN/dx * dN/dx^T + dN/dy * dN/dy^T) * detJ * weight
    for (size_t point = 0; point < dN_dx_matrix.size(); ++point) {
        double detJ = jacobianData[point].detJ;
        double weight = weights[point];
        
        // Get derivatives at this integration point
        const auto& dN_dx = dN_dx_matrix[point]; // [4] - derivatives of shape functions w.r.t. x
        const auto& dN_dy = dN_dy_matrix[point]; // [4] - derivatives of shape functions w.r.t. y
        
        // Calculate contribution: conductivity * (dN/dx * dN/dx^T + dN/dy * dN/dy^T) * detJ * weight
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                double contribution = conductivity * (dN_dx[i] * dN_dx[j] + dN_dy[i] * dN_dy[j]) * detJ * weight;
                H[i][j] += contribution;
            }
        }
    }
    
    return H;
}

std::ostream& operator<<(std::ostream& os, const Element& element) {
    os << "Element " << std::setw(3) << element.id << " [" << element.type << "]: ";
    for (size_t i = 0; i < element.nodeIds.size(); ++i) {
        os << element.nodeIds[i];
        if (i < element.nodeIds.size() - 1) {
            os << ", ";
        }
    }
    return os;
}
