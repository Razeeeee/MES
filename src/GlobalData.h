#pragma once
#include <string>

class GlobalData {
private:
    double simulationTime;
    double simulationStepTime;
    double conductivity;
    double alfa;
    double tot;
    double initialTemp;
    double density;
    double specificHeat;
    int nodesNumber;
    int elementsNumber;

public:
    GlobalData();
    
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
    
    void print() const;
};
