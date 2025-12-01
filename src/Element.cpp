#include "Element.h"
#include "GaussianIntegration.h"
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

// Shape functions
double Element::N1(double xi, double eta) {
    return 0.25 * (1.0 - xi) * (1.0 - eta);
}

double Element::N2(double xi, double eta) {
    return 0.25 * (1.0 + xi) * (1.0 - eta);
}

double Element::N3(double xi, double eta) {
    return 0.25 * (1.0 + xi) * (1.0 + eta);
}

double Element::N4(double xi, double eta) {
    return 0.25 * (1.0 - xi) * (1.0 + eta);
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

std::vector<std::vector<double>> Element::calculateHbcMatrix(const std::vector<double>& nodeX,
                                                            const std::vector<double>& nodeY,
                                                            double alfa,
                                                            const std::vector<bool>& boundaryEdges) const {
    // Initialize 4x4 Hbc matrix with zeros
    std::vector<std::vector<double>> Hbc(4, std::vector<double>(4, 0.0));
    
    // Edge definitions for 4-node quadrilateral:
    // Edge 0 (bottom): nodes 1-2 (local 0-1), eta = -1
    // Edge 1 (right):  nodes 2-3 (local 1-2), xi = 1
    // Edge 2 (top):    nodes 3-4 (local 2-3), eta = 1
    // Edge 3 (left):   nodes 4-1 (local 3-0), xi = -1
    
    // 2-point Gaussian integration for 1D (edge integration)
    double gaussPoint1D = 1.0 / std::sqrt(3.0);
    std::vector<double> gaussPoints1D = {-gaussPoint1D, gaussPoint1D};
    std::vector<double> weights1D = {1.0, 1.0};
    
    // Process each edge
    for (int edge = 0; edge < 4; ++edge) {
        if (!boundaryEdges[edge]) continue;
        
        // Define which local nodes are on this edge and the parametric coordinate
        std::vector<int> edgeNodes(2);
        int fixedCoord; // 0 for xi, 1 for eta
        double fixedValue;
        
        if (edge == 0) { // Bottom edge: nodes 0-1, eta = -1
            edgeNodes = {0, 1};
            fixedCoord = 1; // eta is fixed
            fixedValue = -1.0;
        } else if (edge == 1) { // Right edge: nodes 1-2, xi = 1
            edgeNodes = {1, 2};
            fixedCoord = 0; // xi is fixed
            fixedValue = 1.0;
        } else if (edge == 2) { // Top edge: nodes 2-3, eta = 1
            edgeNodes = {2, 3};
            fixedCoord = 1; // eta is fixed
            fixedValue = 1.0;
        } else { // edge == 3, Left edge: nodes 3-0, xi = -1
            edgeNodes = {3, 0};
            fixedCoord = 0; // xi is fixed
            fixedValue = -1.0;
        }
        
        // Integrate along the edge using 2-point Gaussian quadrature
        for (size_t gp = 0; gp < gaussPoints1D.size(); ++gp) {
            double xi, eta;
            
            if (fixedCoord == 0) { // xi is fixed
                xi = fixedValue;
                eta = gaussPoints1D[gp];
            } else { // eta is fixed
                xi = gaussPoints1D[gp];
                eta = fixedValue;
            }
            
            // Calculate shape functions at this integration point
            std::vector<double> N = {N1(xi, eta), N2(xi, eta), N3(xi, eta), N4(xi, eta)};
            
            // Calculate edge length element (detJ for 1D)
            // For edge integration, we need the length of the edge in physical space
            double detJ_edge;
            
            if (fixedCoord == 0) { // xi is fixed, eta varies
                // d/deta of x and y coordinates
                double dx_deta = dN1_dEta(xi, eta) * nodeX[0] + dN2_dEta(xi, eta) * nodeX[1] +
                                dN3_dEta(xi, eta) * nodeX[2] + dN4_dEta(xi, eta) * nodeX[3];
                double dy_deta = dN1_dEta(xi, eta) * nodeY[0] + dN2_dEta(xi, eta) * nodeY[1] +
                                dN3_dEta(xi, eta) * nodeY[2] + dN4_dEta(xi, eta) * nodeY[3];
                detJ_edge = std::sqrt(dx_deta * dx_deta + dy_deta * dy_deta);
            } else { // eta is fixed, xi varies
                // d/dxi of x and y coordinates
                double dx_dxi = dN1_dXi(xi, eta) * nodeX[0] + dN2_dXi(xi, eta) * nodeX[1] +
                               dN3_dXi(xi, eta) * nodeX[2] + dN4_dXi(xi, eta) * nodeX[3];
                double dy_dxi = dN1_dXi(xi, eta) * nodeY[0] + dN2_dXi(xi, eta) * nodeY[1] +
                               dN3_dXi(xi, eta) * nodeY[2] + dN4_dXi(xi, eta) * nodeY[3];
                detJ_edge = std::sqrt(dx_dxi * dx_dxi + dy_dxi * dy_dxi);
            }
            
            double weight = weights1D[gp];
            
            // Calculate contribution: alfa * N * N^T * detJ_edge * weight
            for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 4; ++j) {
                    Hbc[i][j] += alfa * N[i] * N[j] * detJ_edge * weight;
                }
            }
        }
    }
    
    return Hbc;
}

std::vector<double> Element::calculatePVector(const std::vector<double>& nodeX,
                                               const std::vector<double>& nodeY,
                                               double alfa,
                                               double tot,
                                               const std::vector<bool>& boundaryEdges) const {
    // Initialize 4-element P vector with zeros
    std::vector<double> P(4, 0.0);
    
    // Edge definitions for 4-node quadrilateral:
    // Edge 0 (bottom): nodes 1-2 (local 0-1), eta = -1
    // Edge 1 (right):  nodes 2-3 (local 1-2), xi = 1
    // Edge 2 (top):    nodes 3-4 (local 2-3), eta = 1
    // Edge 3 (left):   nodes 4-1 (local 3-0), xi = -1
    
    // 2-point Gaussian integration for 1D (edge integration)
    double gaussPoint1D = 1.0 / std::sqrt(3.0);
    std::vector<double> gaussPoints1D = {-gaussPoint1D, gaussPoint1D};
    std::vector<double> weights1D = {1.0, 1.0};
    
    // Process each edge
    for (int edge = 0; edge < 4; ++edge) {
        if (!boundaryEdges[edge]) continue;
        
        // Define which local nodes are on this edge and the parametric coordinate
        std::vector<int> edgeNodes(2);
        int fixedCoord; // 0 for xi, 1 for eta
        double fixedValue;
        
        if (edge == 0) { // Bottom edge: nodes 0-1, eta = -1
            edgeNodes = {0, 1};
            fixedCoord = 1; // eta is fixed
            fixedValue = -1.0;
        } else if (edge == 1) { // Right edge: nodes 1-2, xi = 1
            edgeNodes = {1, 2};
            fixedCoord = 0; // xi is fixed
            fixedValue = 1.0;
        } else if (edge == 2) { // Top edge: nodes 2-3, eta = 1
            edgeNodes = {2, 3};
            fixedCoord = 1; // eta is fixed
            fixedValue = 1.0;
        } else { // edge == 3, Left edge: nodes 3-0, xi = -1
            edgeNodes = {3, 0};
            fixedCoord = 0; // xi is fixed
            fixedValue = -1.0;
        }
        
        // Integrate along the edge using 2-point Gaussian quadrature
        for (size_t gp = 0; gp < gaussPoints1D.size(); ++gp) {
            double xi, eta;
            
            if (fixedCoord == 0) { // xi is fixed
                xi = fixedValue;
                eta = gaussPoints1D[gp];
            } else { // eta is fixed
                xi = gaussPoints1D[gp];
                eta = fixedValue;
            }
            
            // Calculate shape functions at this integration point
            std::vector<double> N = {N1(xi, eta), N2(xi, eta), N3(xi, eta), N4(xi, eta)};
            
            // Calculate edge length element (detJ for 1D)
            // For edge integration, we need the length of the edge in physical space
            double detJ_edge;
            
            if (fixedCoord == 0) { // xi is fixed, eta varies
                // d/deta of x and y coordinates
                double dx_deta = dN1_dEta(xi, eta) * nodeX[0] + dN2_dEta(xi, eta) * nodeX[1] +
                                dN3_dEta(xi, eta) * nodeX[2] + dN4_dEta(xi, eta) * nodeX[3];
                double dy_deta = dN1_dEta(xi, eta) * nodeY[0] + dN2_dEta(xi, eta) * nodeY[1] +
                                dN3_dEta(xi, eta) * nodeY[2] + dN4_dEta(xi, eta) * nodeY[3];
                detJ_edge = std::sqrt(dx_deta * dx_deta + dy_deta * dy_deta);
            } else { // eta is fixed, xi varies
                // d/dxi of x and y coordinates
                double dx_dxi = dN1_dXi(xi, eta) * nodeX[0] + dN2_dXi(xi, eta) * nodeX[1] +
                               dN3_dXi(xi, eta) * nodeX[2] + dN4_dXi(xi, eta) * nodeX[3];
                double dy_dxi = dN1_dXi(xi, eta) * nodeY[0] + dN2_dXi(xi, eta) * nodeY[1] +
                               dN3_dXi(xi, eta) * nodeY[2] + dN4_dXi(xi, eta) * nodeY[3];
                detJ_edge = std::sqrt(dx_dxi * dx_dxi + dy_dxi * dy_dxi);
            }
            
            double weight = weights1D[gp];
            
            // Calculate contribution: alfa * tot * N * detJ_edge * weight
            for (int i = 0; i < 4; ++i) {
                P[i] += alfa * tot * N[i] * detJ_edge * weight;
            }
        }
    }
    
    return P;
}

std::vector<std::vector<double>> Element::calculateCMatrix(const std::vector<double>& nodeX,
                                                           const std::vector<double>& nodeY,
                                                           double density,
                                                           double specificHeat,
                                                           int numGaussPoints) const {
    // Initialize 4x4 C matrix with zeros
    std::vector<std::vector<double>> C(4, std::vector<double>(4, 0.0));
    
    // Get Gaussian integration points and weights
    std::vector<double> gaussPoints, weights;
    
    if (numGaussPoints == 2) {
        double gp = 1.0 / std::sqrt(3.0);
        gaussPoints = {-gp, gp};
        weights = {1.0, 1.0};
    } else if (numGaussPoints == 3) {
        gaussPoints = {-std::sqrt(3.0/5.0), 0.0, std::sqrt(3.0/5.0)};
        weights = {5.0/9.0, 8.0/9.0, 5.0/9.0};
    } else {
        throw std::runtime_error("Invalid number of Gauss points. Use 2 or 3.");
    }
    
    // Iterate through all integration points (2D: numGaussPoints x numGaussPoints)
    for (size_t ip_xi = 0; ip_xi < gaussPoints.size(); ++ip_xi) {
        for (size_t ip_eta = 0; ip_eta < gaussPoints.size(); ++ip_eta) {
            double xi = gaussPoints[ip_xi];
            double eta = gaussPoints[ip_eta];
            double weight_xi = weights[ip_xi];
            double weight_eta = weights[ip_eta];
            
            // Calculate shape functions at this integration point
            std::vector<double> N(4);
            N[0] = N1(xi, eta);
            N[1] = N2(xi, eta);
            N[2] = N3(xi, eta);
            N[3] = N4(xi, eta);
            
            // Calculate Jacobian at this integration point
            auto jacobianData = calculateJacobian(nodeX, nodeY, xi, eta);
            double detJ = calculateDetJ(jacobianData);
            
            // Calculate C contribution at this integration point
            // C = density * specificHeat * N * transpose(N) * detJ * weight
            double factor = density * specificHeat * detJ * weight_xi * weight_eta;
            
            // Build N * transpose(N) and add to C matrix
            for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 4; ++j) {
                    C[i][j] += factor * N[i] * N[j];
                }
            }
        }
    }
    
    return C;
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
