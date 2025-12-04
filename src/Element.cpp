#include "Element.h"
#include <iomanip>
#include <cmath>
#include <stdexcept>

// Constructors

Element::Element() : id(0), type("DC2D4") {}

Element::Element(int id, const std::vector<int>& nodeIds, const std::string& type) 
    : id(id), nodeIds(nodeIds), type(type) {}

void Element::print() const {
    std::cout << "Element " << std::setw(3) << id << "]: ";
    for (size_t i = 0; i < nodeIds.size(); ++i) {
        std::cout << nodeIds[i];
        if (i < nodeIds.size() - 1) std::cout << ", ";
    }
    std::cout << std::endl;
}

// Shape functions for 4-node quadrilateral (Q4)
// Node numbering counterclockwise: 1(-1,-1), 2(1,-1), 3(1,1), 4(-1,1)
// Nᵢ(ξ,η) = ¼(1±ξ)(1±η)

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

// ∂N/∂ξ derivatives

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

// ∂N/∂η derivatives

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

// Gaussian quadrature points and weights (1D)
// ============================================================================

void Element::getGaussianPoints(int numPoints, std::vector<double>& points, std::vector<double>& weights) {
    if (numPoints == 2) {
        // 2-point: exact for polynomials up to degree 3
        double gp = 1.0 / std::sqrt(3.0);
        points = {-gp, gp};
        weights = {1.0, 1.0};
    } else if (numPoints == 3) {
        // 3-point: exact for polynomials up to degree 5
        points = {-std::sqrt(3.0/5.0), 0.0, std::sqrt(3.0/5.0)};
        weights = {5.0/9.0, 8.0/9.0, 5.0/9.0};
    } else if (numPoints == 4) {
        // 4-point: exact for polynomials up to degree 7
        double p1 = std::sqrt(3.0/7.0 - 2.0/7.0 * std::sqrt(6.0/5.0));
        double p2 = std::sqrt(3.0/7.0 + 2.0/7.0 * std::sqrt(6.0/5.0));
        double w1 = (18.0 + std::sqrt(30.0)) / 36.0;
        double w2 = (18.0 - std::sqrt(30.0)) / 36.0;
        points = {-p2, -p1, p1, p2};
        weights = {w2, w1, w1, w2};
    } else {
        throw std::invalid_argument("Only 2, 3, and 4-point Gaussian quadrature are supported");
    }
}

// Jacobian transformation: J = [∂x/∂ξ ∂y/∂ξ; ∂x/∂η ∂y/∂η]

std::vector<std::vector<double>> Element::calculateJacobian(
    const std::vector<double>& nodeX, 
    const std::vector<double>& nodeY, 
    double xi, double eta) const {
    
    std::vector<std::vector<double>> J(2, std::vector<double>(2, 0.0));
    
    J[0][0] = dN1_dXi(xi, eta) * nodeX[0] + dN2_dXi(xi, eta) * nodeX[1] + 
              dN3_dXi(xi, eta) * nodeX[2] + dN4_dXi(xi, eta) * nodeX[3];
    J[0][1] = dN1_dXi(xi, eta) * nodeY[0] + dN2_dXi(xi, eta) * nodeY[1] + 
              dN3_dXi(xi, eta) * nodeY[2] + dN4_dXi(xi, eta) * nodeY[3];
    J[1][0] = dN1_dEta(xi, eta) * nodeX[0] + dN2_dEta(xi, eta) * nodeX[1] + 
              dN3_dEta(xi, eta) * nodeX[2] + dN4_dEta(xi, eta) * nodeX[3];
    J[1][1] = dN1_dEta(xi, eta) * nodeY[0] + dN2_dEta(xi, eta) * nodeY[1] + 
              dN3_dEta(xi, eta) * nodeY[2] + dN4_dEta(xi, eta) * nodeY[3];
    
    return J;
}

double Element::calculateDetJ(const std::vector<std::vector<double>>& J) const {
    return J[0][0] * J[1][1] - J[0][1] * J[1][0];
}

std::vector<std::vector<double>> Element::calculateInverseJacobian(
    const std::vector<std::vector<double>>& J) const {
    
    double detJ = calculateDetJ(J);
    if (std::abs(detJ) < 1e-12) {
        throw std::runtime_error("Jacobian determinant is too small (singular element)");
    }
    
    double invDetJ = 1.0 / detJ;
    std::vector<std::vector<double>> invJ(2, std::vector<double>(2));
    
    invJ[0][0] =  J[1][1] * invDetJ;
    invJ[0][1] = -J[0][1] * invDetJ;
    invJ[1][0] = -J[1][0] * invDetJ;
    invJ[1][1] =  J[0][0] * invDetJ;
    
    return invJ;
}

// Transform ∂N/∂ξ,η to ∂N/∂x,y using inverse Jacobian

std::vector<double> Element::calculateDN_Dx(double xi, double eta, 
                                            const std::vector<std::vector<double>>& invJ) {
    std::vector<double> dN_dxi = {dN1_dXi(xi, eta), dN2_dXi(xi, eta), 
                                   dN3_dXi(xi, eta), dN4_dXi(xi, eta)};
    std::vector<double> dN_deta = {dN1_dEta(xi, eta), dN2_dEta(xi, eta), 
                                    dN3_dEta(xi, eta), dN4_dEta(xi, eta)};
    
    std::vector<double> dN_dx(4);
    for (int i = 0; i < 4; ++i) {
        dN_dx[i] = invJ[0][0] * dN_dxi[i] + invJ[0][1] * dN_deta[i];
    }
    return dN_dx;
}

std::vector<double> Element::calculateDN_Dy(double xi, double eta, 
                                            const std::vector<std::vector<double>>& invJ) {
    std::vector<double> dN_dxi = {dN1_dXi(xi, eta), dN2_dXi(xi, eta), 
                                   dN3_dXi(xi, eta), dN4_dXi(xi, eta)};
    std::vector<double> dN_deta = {dN1_dEta(xi, eta), dN2_dEta(xi, eta), 
                                    dN3_dEta(xi, eta), dN4_dEta(xi, eta)};
    
    std::vector<double> dN_dy(4);
    for (int i = 0; i < 4; ++i) {
        dN_dy[i] = invJ[1][0] * dN_dxi[i] + invJ[1][1] * dN_deta[i];
    }
    return dN_dy;
}

// ============================================================================
// COMBINED H AND C MATRIX CALCULATION (OPTIMIZED)
// ============================================================================
// Calculates both conductivity (H) and capacity (C) matrices in a single loop
// to avoid redundant Jacobian calculations.
//
// FEM Theory:
// -----------
// H Matrix (Conductivity/Stiffness):
//   Hᵢⱼ = ∫∫_Ω k·(∂Nᵢ/∂x·∂Nⱼ/∂x + ∂Nᵢ/∂y·∂Nⱼ/∂y) dΩ
//   where k is thermal conductivity
//   Represents heat conduction within the element
//
// C Matrix (Heat Capacity):
//   Cᵢⱼ = ∫∫_Ω ρ·c·Nᵢ·Nⱼ dΩ
//   where ρ is density, c is specific heat
//   Represents thermal energy storage capacity
//
// Both use Gaussian quadrature: ∫∫f dΩ ≈ Σᵢ Σⱼ wᵢ·wⱼ·f(ξᵢ,ηⱼ)·|J|
// ============================================================================

std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> 
Element::calculateHAndCMatrices(
    const std::vector<double>& nodeX,
    const std::vector<double>& nodeY,
    double conductivity,
    double density,
    double specificHeat,
    int numGaussPoints) const {
    
    // Initialize 4x4 matrices
    std::vector<std::vector<double>> H(4, std::vector<double>(4, 0.0));
    std::vector<std::vector<double>> C(4, std::vector<double>(4, 0.0));
    
    // Get Gaussian integration points and weights
    std::vector<double> gp, gw;
    getGaussianPoints(numGaussPoints, gp, gw);
    
    // Loop over all 2D integration points (numGaussPoints × numGaussPoints)
    for (size_t ip_eta = 0; ip_eta < gp.size(); ++ip_eta) {
        for (size_t ip_xi = 0; ip_xi < gp.size(); ++ip_xi) {
            double xi = gp[ip_xi];
            double eta = gp[ip_eta];
            double weight = gw[ip_xi] * gw[ip_eta];
            
            // Calculate Jacobian and its determinant (shared for H and C)
            auto J = calculateJacobian(nodeX, nodeY, xi, eta);
            double detJ = calculateDetJ(J);
            auto invJ = calculateInverseJacobian(J);
            
            // Get physical derivatives for H matrix
            auto dN_dx = calculateDN_Dx(xi, eta, invJ);
            auto dN_dy = calculateDN_Dy(xi, eta, invJ);
            
            // Get shape functions for C matrix
            std::vector<double> N = {N1(xi, eta), N2(xi, eta), N3(xi, eta), N4(xi, eta)};
            
            // Accumulate H and C matrix contributions at this integration point
            double h_factor = conductivity * detJ * weight;
            double c_factor = density * specificHeat * detJ * weight;
            
            for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 4; ++j) {
                    // H matrix: k · (∂Nᵢ/∂x·∂Nⱼ/∂x + ∂Nᵢ/∂y·∂Nⱼ/∂y) · |J| · w
                    H[i][j] += h_factor * (dN_dx[i] * dN_dx[j] + dN_dy[i] * dN_dy[j]);
                    
                    // C matrix: ρ·c · Nᵢ·Nⱼ · |J| · w
                    C[i][j] += c_factor * N[i] * N[j];
                }
            }
        }
    }
    
    return std::make_pair(H, C);
}

// Calculate Hbc (boundary convection) and P (heat flux) on element edges

std::pair<std::vector<std::vector<double>>, std::vector<double>>
Element::calculateHbcAndPVector(
    const std::vector<double>& nodeX,
    const std::vector<double>& nodeY,
    double alfa,
    double tot,
    const std::vector<bool>& boundaryEdges,
    int numGaussPoints) const {
    
    std::vector<std::vector<double>> Hbc(4, std::vector<double>(4, 0.0));
    std::vector<double> P(4, 0.0);
    
    std::vector<double> gp, gw;
    getGaussianPoints(numGaussPoints, gp, gw);
    
    // Edge order: 0=bottom(η=-1), 1=right(ξ=1), 2=top(η=1), 3=left(ξ=-1)
    for (int edge = 0; edge < 4; ++edge) {
        if (!boundaryEdges[edge]) continue;
        
        int fixedCoord;
        double fixedValue;
        
        if (edge == 0) { fixedCoord = 1; fixedValue = -1.0; }
        else if (edge == 1) { fixedCoord = 0; fixedValue = 1.0; }
        else if (edge == 2) { fixedCoord = 1; fixedValue = 1.0; }
        else { fixedCoord = 0; fixedValue = -1.0; }
        
        for (size_t igp = 0; igp < gp.size(); ++igp) {
            double xi = (fixedCoord == 0) ? fixedValue : gp[igp];
            double eta = (fixedCoord == 1) ? fixedValue : gp[igp];
            
            std::vector<double> N = {N1(xi, eta), N2(xi, eta), N3(xi, eta), N4(xi, eta)};
            
            double detJ_edge;
            if (fixedCoord == 0) {
                double dx_deta = dN1_dEta(xi, eta) * nodeX[0] + dN2_dEta(xi, eta) * nodeX[1] +
                                 dN3_dEta(xi, eta) * nodeX[2] + dN4_dEta(xi, eta) * nodeX[3];
                double dy_deta = dN1_dEta(xi, eta) * nodeY[0] + dN2_dEta(xi, eta) * nodeY[1] +
                                 dN3_dEta(xi, eta) * nodeY[2] + dN4_dEta(xi, eta) * nodeY[3];
                detJ_edge = std::sqrt(dx_deta * dx_deta + dy_deta * dy_deta);
            } else {
                double dx_dxi = dN1_dXi(xi, eta) * nodeX[0] + dN2_dXi(xi, eta) * nodeX[1] +
                                dN3_dXi(xi, eta) * nodeX[2] + dN4_dXi(xi, eta) * nodeX[3];
                double dy_dxi = dN1_dXi(xi, eta) * nodeY[0] + dN2_dXi(xi, eta) * nodeY[1] +
                                dN3_dXi(xi, eta) * nodeY[2] + dN4_dXi(xi, eta) * nodeY[3];
                detJ_edge = std::sqrt(dx_dxi * dx_dxi + dy_dxi * dy_dxi);
            }
            
            double weight = gw[igp];
            double hbc_factor = alfa * detJ_edge * weight;
            double p_factor = alfa * tot * detJ_edge * weight;
            
            for (int i = 0; i < 4; ++i) {
                P[i] += p_factor * N[i];
                for (int j = 0; j < 4; ++j) {
                    Hbc[i][j] += hbc_factor * N[i] * N[j];
                }
            }
        }
    }
    
    return std::make_pair(Hbc, P);
}

std::ostream& operator<<(std::ostream& os, const Element& element) {
    os << "Element " << std::setw(3) << element.id << " [" << element.type << "]: ";
    for (size_t i = 0; i < element.nodeIds.size(); ++i) {
        os << element.nodeIds[i];
        if (i < element.nodeIds.size() - 1) os << ", ";
    }
    return os;
}
