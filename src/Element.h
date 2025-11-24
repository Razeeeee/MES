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
    static double dN4_dEta(double xi, double /*eta*/);
    
    // Shape functions for boundary condition matrix
    static double N1(double xi, double eta);
    static double N2(double xi, double eta);
    static double N3(double xi, double eta);
    static double N4(double xi, double eta);
    
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
    
    // Calculate shape function derivatives with respect to physical coordinates (x, y)
    // Using chain rule: [dN/dx] = [J^-1] * [dN/dxi ]
    //                  [dN/dy]           [dN/deta]
    
    // Calculate dN/dx for all shape functions at a specific integration point
    static std::vector<double> calculateDN_Dx(double xi, double eta, 
                                              const std::vector<std::vector<double>>& inverseJacobian);
    
    // Calculate dN/dy for all shape functions at a specific integration point  
    static std::vector<double> calculateDN_Dy(double xi, double eta, 
                                              const std::vector<std::vector<double>>& inverseJacobian);
    
    // Calculate matrices of shape function derivatives with respect to x and y
    // at all integration points for this element
    // Returns pair of matrices: {dN_dx_matrix, dN_dy_matrix}
    // Each matrix has rows = integration points, columns = shape functions (4 for DC2D4)
    std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> 
    calculatePhysicalDerivativeMatrices(const std::vector<double>& nodeX, 
                                       const std::vector<double>& nodeY, 
                                       int numGaussPoints) const;
    
    // Calculate H matrix for heat conduction using physical derivatives
    // H = conductivity * (dN/dx * transpose(dN/dx) + dN/dy * transpose(dN/dy)) * detJ
    // Sums contributions from all integration points
    std::vector<std::vector<double>> calculateHMatrix(const std::vector<double>& nodeX, 
                                                     const std::vector<double>& nodeY, 
                                                     double conductivity,
                                                     int numGaussPoints) const;
    
    // Calculate Hbc matrix for boundary conditions
    // Hbc = alfa * N * transpose(N) * detJ
    // Only calculates for edges with boundary conditions
    // Parameters: nodeX, nodeY - coordinates of element nodes
    //            alfa - heat transfer coefficient
    //            boundaryEdges - which edges have BC (0=bottom, 1=right, 2=top, 3=left)
    std::vector<std::vector<double>> calculateHbcMatrix(const std::vector<double>& nodeX,
                                                        const std::vector<double>& nodeY,
                                                        double alfa,
                                                        const std::vector<bool>& boundaryEdges) const;
    
    // Calculate Hbc matrix with node-specific alfa values
    // nodeAlfa - alfa value for each node (4 values)
    std::vector<std::vector<double>> calculateHbcMatrix(const std::vector<double>& nodeX,
                                                        const std::vector<double>& nodeY,
                                                        const std::vector<double>& nodeAlfa,
                                                        const std::vector<bool>& boundaryEdges) const;
    
    // Calculate P vector for boundary conditions
    // P = alfa * tot * N * detJ
    // Only calculates for edges with boundary conditions
    // Parameters: nodeX, nodeY - coordinates of element nodes
    //            alfa - heat transfer coefficient
    //            tot - ambient temperature
    //            boundaryEdges - which edges have BC (0=bottom, 1=right, 2=top, 3=left)
    std::vector<double> calculatePVector(const std::vector<double>& nodeX,
                                        const std::vector<double>& nodeY,
                                        double alfa,
                                        double tot,
                                        const std::vector<bool>& boundaryEdges) const;
    
    // Calculate P vector with node-specific alfa and tot values
    // nodeAlfa - alfa value for each node (4 values)
    // nodeTot - tot value for each node (4 values)
    std::vector<double> calculatePVector(const std::vector<double>& nodeX,
                                        const std::vector<double>& nodeY,
                                        const std::vector<double>& nodeAlfa,
                                        const std::vector<double>& nodeTot,
                                        const std::vector<bool>& boundaryEdges) const;
    
    // Friend operator for output stream
    friend std::ostream& operator<<(std::ostream& os, const Element& element);
};
