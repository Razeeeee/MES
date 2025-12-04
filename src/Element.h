#pragma once
#include <vector>
#include <iostream>
#include <utility>

// 4-node quadrilateral element for 2D heat transfer (supports 2/3/4-point Gaussian quadrature)
class Element {
private:
    int id;
    std::vector<int> nodeIds;
    std::string type;
    
    static void getGaussianPoints(int numPoints, std::vector<double>& points, std::vector<double>& weights);

public:
    Element();
    Element(int id, const std::vector<int>& nodeIds, const std::string& type = "DC2D4");
    
    // Getters
    int getId() const noexcept { return id; }
    const std::vector<int>& getNodeIds() const noexcept { return nodeIds; }
    const std::string& getType() const noexcept { return type; }
    
    // Setters
    void setId(int id) noexcept { this->id = id; }
    void setNodeIds(const std::vector<int>& nodeIds) { this->nodeIds = nodeIds; }
    void setType(const std::string& type) { this->type = type; }
    
    void print() const;
    
    // Shape functions N(ξ,η) in natural coordinates
    static double N1(double xi, double eta);
    static double N2(double xi, double eta);
    static double N3(double xi, double eta);
    static double N4(double xi, double eta);
    
    // Shape function derivatives ∂N/∂ξ
    static double dN1_dXi(double xi, double eta);
    static double dN2_dXi(double xi, double eta);
    static double dN3_dXi(double xi, double eta);
    static double dN4_dXi(double xi, double eta);
    
    // Shape function derivatives ∂N/∂η
    static double dN1_dEta(double xi, double eta);
    static double dN2_dEta(double xi, double eta);
    static double dN3_dEta(double xi, double eta);
    static double dN4_dEta(double xi, double eta);
    
    // Jacobian transformation
    std::vector<std::vector<double>> calculateJacobian(
        const std::vector<double>& nodeX, 
        const std::vector<double>& nodeY, 
        double xi, double eta) const;
    
    double calculateDetJ(const std::vector<std::vector<double>>& jacobian) const;
    
    std::vector<std::vector<double>> calculateInverseJacobian(
        const std::vector<std::vector<double>>& jacobian) const;
    
    // Physical coordinate derivatives ∂N/∂x and ∂N/∂y
    static std::vector<double> calculateDN_Dx(
        double xi, double eta, 
        const std::vector<std::vector<double>>& inverseJacobian);
    
    static std::vector<double> calculateDN_Dy(
        double xi, double eta, 
        const std::vector<std::vector<double>>& inverseJacobian);
    
    // Calculate H (conductivity) and C (capacity) matrices simultaneously
    std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> 
    calculateHAndCMatrices(
        const std::vector<double>& nodeX,
        const std::vector<double>& nodeY,
        double conductivity,
        double density,
        double specificHeat,
        int numGaussPoints) const;
    
    // Calculate Hbc (boundary convection) and P (heat flux) simultaneously
    std::pair<std::vector<std::vector<double>>, std::vector<double>>
    calculateHbcAndPVector(
        const std::vector<double>& nodeX,
        const std::vector<double>& nodeY,
        double alfa,
        double tot,
        const std::vector<bool>& boundaryEdges,
        int numGaussPoints) const;
    
    // Output stream operator
    friend std::ostream& operator<<(std::ostream& os, const Element& element);
};
