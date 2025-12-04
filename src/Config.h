#pragma once
#include <string>
#include <map>

// FEM simulation configuration (Gaussian points and output control)
class Config {
private:
    int gaussPointsH;
    int gaussPointsC;
    int gaussPointsBoundary;
    
    bool printLocalMatrices;
    bool printGlobalMatrices;
    bool printAssemblyMapping;
    bool printSystemInfo;
    bool printBoundaryConditions;
    bool printTimeSteps;
    bool printFinalResults;
    bool printSteadyState;
    
    std::string configFilePath;
    
    void setDefaults();
    std::string trim(const std::string& str);
    void parseLine(const std::string& line);
    bool parseBool(const std::string& value);
    int parseInt(const std::string& value);

public:
    Config();
    explicit Config(const std::string& configFile);
    
    bool loadFromFile(const std::string& configFile);
    
    int getGaussPointsH() const noexcept { return gaussPointsH; }
    int getGaussPointsC() const noexcept { return gaussPointsC; }
    int getGaussPointsBoundary() const noexcept { return gaussPointsBoundary; }
    
    bool shouldPrintLocalMatrices() const noexcept { return printLocalMatrices; }
    bool shouldPrintGlobalMatrices() const noexcept { return printGlobalMatrices; }
    bool shouldPrintAssemblyMapping() const noexcept { return printAssemblyMapping; }
    bool shouldPrintSystemInfo() const noexcept { return printSystemInfo; }
    bool shouldPrintBoundaryConditions() const noexcept { return printBoundaryConditions; }
    bool shouldPrintTimeSteps() const noexcept { return printTimeSteps; }
    bool shouldPrintFinalResults() const noexcept { return printFinalResults; }
    bool shouldPrintSteadyState() const noexcept { return printSteadyState; }
    
    // Display current configuration
    void print() const;
};
