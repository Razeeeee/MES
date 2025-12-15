#include "Element.h"
#include <iomanip>
#include <cmath>
#include <stdexcept>

ShapeFunctions ShapeFunctions::calculate(double xi, double eta) {
    ShapeFunctions sf(4);
    
    sf.N[0] = 0.25 * (1.0 - xi) * (1.0 - eta);
    sf.N[1] = 0.25 * (1.0 + xi) * (1.0 - eta);
    sf.N[2] = 0.25 * (1.0 + xi) * (1.0 + eta);
    sf.N[3] = 0.25 * (1.0 - xi) * (1.0 + eta);
    
    sf.dN_dXi[0] = -0.25 * (1.0 - eta);
    sf.dN_dXi[1] =  0.25 * (1.0 - eta);
    sf.dN_dXi[2] =  0.25 * (1.0 + eta);
    sf.dN_dXi[3] = -0.25 * (1.0 + eta);
    
    sf.dN_dEta[0] = -0.25 * (1.0 - xi);
    sf.dN_dEta[1] = -0.25 * (1.0 + xi);
    sf.dN_dEta[2] =  0.25 * (1.0 + xi);
    sf.dN_dEta[3] =  0.25 * (1.0 - xi);
    
    return sf;
}

void BoundarySurface::getGaussianPoints(int numPoints, std::vector<double>& points, std::vector<double>& weights) {
    if (numPoints == 2) {
        double gp = 1.0 / std::sqrt(3.0);
        points = {-gp, gp};
        weights = {1.0, 1.0};
    } else if (numPoints == 3) {
        points = {-std::sqrt(3.0/5.0), 0.0, std::sqrt(3.0/5.0)};
        weights = {5.0/9.0, 8.0/9.0, 5.0/9.0};
    } else if (numPoints == 4) {
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

BoundarySurface::BoundarySurface(int edge, int numGaussPoints) 
    : edgeIndex(edge) {
    BoundarySurface::getGaussianPoints(numGaussPoints, gaussPoints, gaussWeights);
    
    for (size_t i = 0; i < gaussPoints.size(); ++i) {
        double xi, eta;
        if (edge == 0) { xi = gaussPoints[i]; eta = -1.0; }
        else if (edge == 1) { xi = 1.0; eta = gaussPoints[i]; }
        else if (edge == 2) { xi = gaussPoints[i]; eta = 1.0; }
        else { xi = -1.0; eta = gaussPoints[i]; }
        
        shapeFuncsAtGP.push_back(ShapeFunctions::calculate(xi, eta));
    }
}

Element::Element() : id(0), type("DC2D4"), 
    H_local(4, std::vector<double>(4, 0.0)),
    C_local(4, std::vector<double>(4, 0.0)),
    Hbc_local(4, std::vector<double>(4, 0.0)),
    P_local(4, 0.0) {}

Element::Element(int id, const std::vector<int>& nodeIds, const std::string& type) 
    : id(id), nodeIds(nodeIds), type(type),
    H_local(4, std::vector<double>(4, 0.0)),
    C_local(4, std::vector<double>(4, 0.0)),
    Hbc_local(4, std::vector<double>(4, 0.0)),
    P_local(4, 0.0) {}

void Element::print() const {
    std::cout << "Element " << std::setw(3) << id << "]: ";
    for (size_t i = 0; i < nodeIds.size(); ++i) {
        std::cout << nodeIds[i];
        if (i < nodeIds.size() - 1) std::cout << ", ";
    }
    std::cout << std::endl;
}

JacobianData Element::calculateJacobian(const std::vector<double>& nodeX, 
                                        const std::vector<double>& nodeY, 
                                        double xi, double eta) const {
    JacobianData jac;
    ShapeFunctions sf = ShapeFunctions::calculate(xi, eta);
    
    jac.J[0][0] = sf.dN_dXi[0] * nodeX[0] + sf.dN_dXi[1] * nodeX[1] + 
                  sf.dN_dXi[2] * nodeX[2] + sf.dN_dXi[3] * nodeX[3];
    jac.J[0][1] = sf.dN_dXi[0] * nodeY[0] + sf.dN_dXi[1] * nodeY[1] + 
                  sf.dN_dXi[2] * nodeY[2] + sf.dN_dXi[3] * nodeY[3];
    jac.J[1][0] = sf.dN_dEta[0] * nodeX[0] + sf.dN_dEta[1] * nodeX[1] + 
                  sf.dN_dEta[2] * nodeX[2] + sf.dN_dEta[3] * nodeX[3];
    jac.J[1][1] = sf.dN_dEta[0] * nodeY[0] + sf.dN_dEta[1] * nodeY[1] + 
                  sf.dN_dEta[2] * nodeY[2] + sf.dN_dEta[3] * nodeY[3];
    
    jac.detJ = jac.J[0][0] * jac.J[1][1] - jac.J[0][1] * jac.J[1][0];
    
    if (std::abs(jac.detJ) < 1e-12) {
        throw std::runtime_error("Jacobian determinant is too small");
    }
    
    double invDetJ = 1.0 / jac.detJ;
    jac.invJ[0][0] =  jac.J[1][1] * invDetJ;
    jac.invJ[0][1] = -jac.J[0][1] * invDetJ;
    jac.invJ[1][0] = -jac.J[1][0] * invDetJ;
    jac.invJ[1][1] =  jac.J[0][0] * invDetJ;
    
    return jac;
}

std::vector<double> Element::calculateDN_Dx(const ShapeFunctions& sf, 
                                            const JacobianData& jac) const {
    std::vector<double> dN_dx(4);
    for (int i = 0; i < 4; ++i) {
        dN_dx[i] = jac.invJ[0][0] * sf.dN_dXi[i] + jac.invJ[0][1] * sf.dN_dEta[i];
    }
    return dN_dx;
}

std::vector<double> Element::calculateDN_Dy(const ShapeFunctions& sf, 
                                            const JacobianData& jac) const {
    std::vector<double> dN_dy(4);
    for (int i = 0; i < 4; ++i) {
        dN_dy[i] = jac.invJ[1][0] * sf.dN_dXi[i] + jac.invJ[1][1] * sf.dN_dEta[i];
    }
    return dN_dy;
}

void Element::calculateLocalMatrices(const std::vector<double>& nodeX,
                                     const std::vector<double>& nodeY,
                                     double conductivity,
                                     double density,
                                     double specificHeat,
                                     double alfa,
                                     double tot,
                                     const std::vector<bool>& boundaryEdges,
                                     int numGaussPoints) {
    
    H_local.assign(4, std::vector<double>(4, 0.0));
    C_local.assign(4, std::vector<double>(4, 0.0));
    Hbc_local.assign(4, std::vector<double>(4, 0.0));
    P_local.assign(4, 0.0);
    
    std::vector<double> gp, gw;
    BoundarySurface::getGaussianPoints(numGaussPoints, gp, gw);
    
    for (size_t ip_eta = 0; ip_eta < gp.size(); ++ip_eta) {
        for (size_t ip_xi = 0; ip_xi < gp.size(); ++ip_xi) {
            double xi = gp[ip_xi];
            double eta = gp[ip_eta];
            double weight = gw[ip_xi] * gw[ip_eta];
            
            JacobianData jac = calculateJacobian(nodeX, nodeY, xi, eta);
            ShapeFunctions sf = ShapeFunctions::calculate(xi, eta);
            
            std::vector<double> dN_dx = calculateDN_Dx(sf, jac);
            std::vector<double> dN_dy = calculateDN_Dy(sf, jac);
            
            double h_factor = conductivity * jac.detJ * weight;
            double c_factor = density * specificHeat * jac.detJ * weight;
            
            for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 4; ++j) {
                    H_local[i][j] += h_factor * (dN_dx[i] * dN_dx[j] + dN_dy[i] * dN_dy[j]);
                    C_local[i][j] += c_factor * sf.N[i] * sf.N[j];
                }
            }
        }
    }
    
    for (int edge = 0; edge < 4; ++edge) {
        if (!boundaryEdges[edge]) continue;
        
        BoundarySurface surface(edge, numGaussPoints);
        
        for (size_t igp = 0; igp < surface.gaussPoints.size(); ++igp) {
            ShapeFunctions sf = surface.shapeFuncsAtGP[igp];
            
            double detJ_edge;
            if (edge == 0 || edge == 2) {
                double dx_dxi = sf.dN_dXi[0] * nodeX[0] + sf.dN_dXi[1] * nodeX[1] +
                                sf.dN_dXi[2] * nodeX[2] + sf.dN_dXi[3] * nodeX[3];
                double dy_dxi = sf.dN_dXi[0] * nodeY[0] + sf.dN_dXi[1] * nodeY[1] +
                                sf.dN_dXi[2] * nodeY[2] + sf.dN_dXi[3] * nodeY[3];
                detJ_edge = std::sqrt(dx_dxi * dx_dxi + dy_dxi * dy_dxi);
            } else {
                double dx_deta = sf.dN_dEta[0] * nodeX[0] + sf.dN_dEta[1] * nodeX[1] +
                                 sf.dN_dEta[2] * nodeX[2] + sf.dN_dEta[3] * nodeX[3];
                double dy_deta = sf.dN_dEta[0] * nodeY[0] + sf.dN_dEta[1] * nodeY[1] +
                                 sf.dN_dEta[2] * nodeY[2] + sf.dN_dEta[3] * nodeY[3];
                detJ_edge = std::sqrt(dx_deta * dx_deta + dy_deta * dy_deta);
            }
            
            double weight = surface.gaussWeights[igp];
            double hbc_factor = alfa * detJ_edge * weight;
            double p_factor = alfa * tot * detJ_edge * weight;
            
            for (int i = 0; i < 4; ++i) {
                P_local[i] += p_factor * sf.N[i];
                for (int j = 0; j < 4; ++j) {
                    Hbc_local[i][j] += hbc_factor * sf.N[i] * sf.N[j];
                }
            }
        }
    }
}

std::ostream& operator<<(std::ostream& os, const Element& element) {
    os << "Element " << std::setw(3) << element.id << " [" << element.type << "]: ";
    for (size_t i = 0; i < element.nodeIds.size(); ++i) {
        os << element.nodeIds[i];
        if (i < element.nodeIds.size() - 1) os << ", ";
    }
    return os;
}
