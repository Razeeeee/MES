#include "GlobalData.h"
#include <iostream>
#include <iomanip>
#include <stdexcept>

GlobalData::GlobalData() 
    : simulationTime(0.0), simulationStepTime(0.0), conductivity(0.0),
      alfa(0.0), tot(0.0), initialTemp(0.0), density(0.0), specificHeat(0.0),
      nodesNumber(0), elementsNumber(0), hasMultipleMaterials(false) {}

void GlobalData::print() const {
    std::cout << "\n=== GLOBAL DATA ===" << std::endl;
    std::cout << std::left;
    std::cout << std::setw(20) << "Simulation Time:" << simulationTime << std::endl;
    std::cout << std::setw(20) << "Step Time:" << simulationStepTime << std::endl;
    std::cout << std::setw(20) << "Conductivity:" << conductivity << std::endl;
    std::cout << std::setw(20) << "Alfa:" << alfa << std::endl;
    std::cout << std::setw(20) << "Tot:" << tot << std::endl;
    std::cout << std::setw(20) << "Initial Temp:" << initialTemp << std::endl;
    std::cout << std::setw(20) << "Density:" << density << std::endl;
    std::cout << std::setw(20) << "Specific Heat:" << specificHeat << std::endl;
    std::cout << std::setw(20) << "Nodes Number:" << nodesNumber << std::endl;
    std::cout << std::setw(20) << "Elements Number:" << elementsNumber << std::endl;
    std::cout << "===================" << std::endl;
    
    // Print materials if multi-material mode
    if (hasMultipleMaterials && !materials.empty()) {
        std::cout << "\n=== MATERIALS ===" << std::endl;
        for (const auto& [id, mat] : materials) {
            std::cout << "Material " << id << " (" << mat.name << "):" << std::endl;
            std::cout << "  Conductivity: " << mat.conductivity << " W/(m·K)" << std::endl;
            std::cout << "  Density: " << mat.density << " kg/m³" << std::endl;
            std::cout << "  Specific Heat: " << mat.specificHeat << " J/(kg·K)" << std::endl;
        }
        std::cout << "==================" << std::endl;
    }
}

void GlobalData::addMaterial(const Material& material) {
    materials[material.id] = material;
    hasMultipleMaterials = true;
}

bool GlobalData::hasMaterial(int materialId) const {
    return materials.find(materialId) != materials.end();
}

const Material& GlobalData::getMaterial(int materialId) const {
    auto it = materials.find(materialId);
    if (it != materials.end()) {
        return it->second;
    }
    throw std::runtime_error("Material ID " + std::to_string(materialId) + " not found");
}

double GlobalData::getMaterialConductivity(int materialId) const {
    if (hasMultipleMaterials && materialId > 0 && hasMaterial(materialId)) {
        return getMaterial(materialId).conductivity;
    }
    return conductivity;  // Fallback to legacy
}

double GlobalData::getMaterialDensity(int materialId) const {
    if (hasMultipleMaterials && materialId > 0 && hasMaterial(materialId)) {
        return getMaterial(materialId).density;
    }
    return density;  // Fallback to legacy
}

double GlobalData::getMaterialSpecificHeat(int materialId) const {
    if (hasMultipleMaterials && materialId > 0 && hasMaterial(materialId)) {
        return getMaterial(materialId).specificHeat;
    }
    return specificHeat;  // Fallback to legacy
}
