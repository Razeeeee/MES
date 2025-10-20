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
