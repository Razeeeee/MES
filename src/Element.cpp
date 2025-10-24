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
    invJacobian[0][0] =  jacobian[1][1] / detJ;  // J22/detJ
    invJacobian[0][1] = -jacobian[0][1] / detJ;  // -J12/detJ
    invJacobian[1][0] = -jacobian[1][0] / detJ;  // -J21/detJ
    invJacobian[1][1] =  jacobian[0][0] / detJ;  // J11/detJ
    
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
