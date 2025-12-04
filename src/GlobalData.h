#pragma once
#include <string>
#include <map>

/**
 * @brief Material properties for heat transfer
 */
struct Material {
    int id;
    std::string name;
    double conductivity;
    double density;
    double specificHeat;
    
    Material() : id(0), name(""), conductivity(0.0), density(0.0), specificHeat(0.0) {}
    Material(int id, const std::string& name, double k, double rho, double c)
        : id(id), name(name), conductivity(k), density(rho), specificHeat(c) {}
};

/**
 * @brief Stores global simulation parameters
 * 
 * Contains all the global data parameters from the grid file
 * such as simulation time, material properties, etc.
 * Supports both single-material (legacy) and multi-material definitions
 */
class GlobalData {
private:
    double simulationTime;
    double simulationStepTime;
    double conductivity;        // Legacy: default material conductivity
    double alfa;
    double tot;
    double initialTemp;
    double density;              // Legacy: default material density
    double specificHeat;         // Legacy: default material specific heat
    int nodesNumber;
    int elementsNumber;
    
    // Multi-material support
    std::map<int, Material> materials;  // Material ID -> Material
    bool hasMultipleMaterials;

public:
    GlobalData();
    
    // Getters
    double getSimulationTime() const noexcept { return simulationTime; }
    double getSimulationStepTime() const noexcept { return simulationStepTime; }
    double getConductivity() const noexcept { return conductivity; }
    double getAlfa() const noexcept { return alfa; }
    double getTot() const noexcept { return tot; }
    double getInitialTemp() const noexcept { return initialTemp; }
    double getDensity() const noexcept { return density; }
    double getSpecificHeat() const noexcept { return specificHeat; }
    int getNodesNumber() const noexcept { return nodesNumber; }
    int getElementsNumber() const noexcept { return elementsNumber; }
    
    // Setters
    void setSimulationTime(double value) noexcept { simulationTime = value; }
    void setSimulationStepTime(double value) noexcept { simulationStepTime = value; }
    void setConductivity(double value) noexcept { conductivity = value; }
    void setAlfa(double value) noexcept { alfa = value; }
    void setTot(double value) noexcept { tot = value; }
    void setInitialTemp(double value) noexcept { initialTemp = value; }
    void setDensity(double value) noexcept { density = value; }
    void setSpecificHeat(double value) noexcept { specificHeat = value; }
    void setNodesNumber(int value) noexcept { nodesNumber = value; }
    void setElementsNumber(int value) noexcept { elementsNumber = value; }
    
    // Multi-material support
    void addMaterial(const Material& material);
    bool hasMaterial(int materialId) const;
    const Material& getMaterial(int materialId) const;
    const std::map<int, Material>& getMaterials() const noexcept { return materials; }
    bool hasMultiMaterialSupport() const noexcept { return hasMultipleMaterials; }
    
    // Get material properties (with fallback to legacy for backward compatibility)
    double getMaterialConductivity(int materialId = 0) const;
    double getMaterialDensity(int materialId = 0) const;
    double getMaterialSpecificHeat(int materialId = 0) const;
    
    // Display global data information
    void print() const;
};
