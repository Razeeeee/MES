#pragma once
#include <vector>
#include <iostream>

/**
 * @brief Represents a finite element in the mesh
 * 
 * Contains element ID and node IDs that form the element
 * Supports DC2D4 (4-node quadrilateral) elements
 */
class Element {
private:
    int id;
    std::vector<int> nodeIds;
    std::string type;

public:
    Element();
    Element(int id, const std::vector<int>& nodeIds, const std::string& type = "DC2D4");
    
    // Getters
    int getId() const { return id; }
    const std::vector<int>& getNodeIds() const { return nodeIds; }
    const std::string& getType() const { return type; }
    
    // Setters
    void setId(int id) { this->id = id; }
    void setNodeIds(const std::vector<int>& nodeIds) { this->nodeIds = nodeIds; }
    void setType(const std::string& type) { this->type = type; }
    
    // Display element information
    void print() const;
    
    // Shape function derivatives for Jacobian calculation
    // For 4-node quadrilateral element
    static double dN1_dXi(double xi, double eta);
    static double dN2_dXi(double xi, double eta);
    static double dN3_dXi(double xi, double eta);
    static double dN4_dXi(double xi, double eta);
    
    static double dN1_dEta(double xi, double eta);
    static double dN2_dEta(double xi, double eta);
    static double dN3_dEta(double xi, double eta);
    static double dN4_dEta(double xi, double eta);
    
    // Calculate matrix of shape function derivatives at integration points
    // Returns matrix where rows are integration points and columns are shape functions
    // For 2x2 integration: 4x4 matrix, for 3x3 integration: 9x4 matrix
    static std::vector<std::vector<double>> calculateShapeFunctionDerivativeMatrix_Xi(int numPoints);
    static std::vector<std::vector<double>> calculateShapeFunctionDerivativeMatrix_Eta(int numPoints);
    
    // Jacobian calculation methods
    // Calculate Jacobian matrix (2x2) at a specific integration point
    std::vector<std::vector<double>> calculateJacobian(const std::vector<double>& nodeX, 
                                                       const std::vector<double>& nodeY, 
                                                       double xi, double eta) const;
    
    // Calculate inverse Jacobian matrix (2x2)
    std::vector<std::vector<double>> calculateInverseJacobian(const std::vector<std::vector<double>>& jacobian) const;
    
    // Calculate determinant of Jacobian matrix
    double calculateDetJ(const std::vector<std::vector<double>>& jacobian) const;
    
    // Calculate Jacobians for all integration points for this element
    // Returns vector of {jacobian, inverse_jacobian, detJ} for each integration point
    struct JacobianData {
        std::vector<std::vector<double>> jacobian;
        std::vector<std::vector<double>> inverseJacobian;
        double detJ;
    };
    
    std::vector<JacobianData> calculateElementJacobians(const std::vector<double>& nodeX, 
                                                        const std::vector<double>& nodeY, 
                                                        int numGaussPoints) const;
    
    // Friend operator for output stream
    friend std::ostream& operator<<(std::ostream& os, const Element& element);
};
