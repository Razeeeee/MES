#pragma once
#include <vector>
#include <iostream>
#include <utility>

struct ShapeFunctions {
    std::vector<double> N;
    std::vector<double> dN_dXi;
    std::vector<double> dN_dEta;
    
    ShapeFunctions(int numNodes = 4) 
        : N(numNodes), dN_dXi(numNodes), dN_dEta(numNodes) {}
    
    static ShapeFunctions calculate(double xi, double eta);
};

struct JacobianData {
    std::vector<std::vector<double>> J;
    std::vector<std::vector<double>> invJ;
    double detJ;
    
    JacobianData() : J(2, std::vector<double>(2)), invJ(2, std::vector<double>(2)), detJ(0.0) {}
};

struct BoundarySurface {
    int edgeIndex;
    std::vector<double> gaussPoints;
    std::vector<double> gaussWeights;
    std::vector<ShapeFunctions> shapeFuncsAtGP;
    
    BoundarySurface(int edge, int numGaussPoints);
    
    static void getGaussianPoints(int numPoints, std::vector<double>& points, std::vector<double>& weights);
};

class Element {
private:
    int id;
    std::vector<int> nodeIds;
    std::string type;
    
    std::vector<std::vector<double>> H_local;
    std::vector<std::vector<double>> C_local;
    std::vector<std::vector<double>> Hbc_local;
    std::vector<double> P_local;
    
    std::vector<JacobianData> jacobiansAtGP;

public:
    Element();
    Element(int id, const std::vector<int>& nodeIds, const std::string& type = "DC2D4");
    
    int getId() const noexcept { return id; }
    const std::vector<int>& getNodeIds() const noexcept { return nodeIds; }
    const std::string& getType() const noexcept { return type; }
    
    const std::vector<std::vector<double>>& getHLocal() const noexcept { return H_local; }
    const std::vector<std::vector<double>>& getCLocal() const noexcept { return C_local; }
    const std::vector<std::vector<double>>& getHbcLocal() const noexcept { return Hbc_local; }
    const std::vector<double>& getPLocal() const noexcept { return P_local; }
    
    void setId(int id) noexcept { this->id = id; }
    void setNodeIds(const std::vector<int>& nodeIds) { this->nodeIds = nodeIds; }
    void setType(const std::string& type) { this->type = type; }
    
    void print() const;
    
    JacobianData calculateJacobian(const std::vector<double>& nodeX, 
                                    const std::vector<double>& nodeY, 
                                    double xi, double eta) const;
    
    std::vector<double> calculateDN_Dx(const ShapeFunctions& sf, 
                                       const JacobianData& jac) const;
    
    std::vector<double> calculateDN_Dy(const ShapeFunctions& sf, 
                                       const JacobianData& jac) const;
    
    void calculateLocalMatrices(const std::vector<double>& nodeX,
                               const std::vector<double>& nodeY,
                               double conductivity,
                               double density,
                               double specificHeat,
                               double alfa,
                               double tot,
                               const std::vector<bool>& boundaryEdges,
                               int numGaussPoints);
    
    friend std::ostream& operator<<(std::ostream& os, const Element& element);
};
