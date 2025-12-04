#include "Config.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <cctype>

Config::Config() {
    setDefaults();
}

Config::Config(const std::string& configFile) {
    setDefaults();
    loadFromFile(configFile);
}

void Config::setDefaults() {
    // Default Gaussian integration points
    gaussPointsH = 2;
    gaussPointsC = 2;
    gaussPointsBoundary = 2;
    
    // Default output flags (all enabled for backward compatibility)
    printLocalMatrices = true;
    printGlobalMatrices = true;
    printAssemblyMapping = true;
    printSystemInfo = true;
    printBoundaryConditions = true;
    printTimeSteps = true;
    printFinalResults = true;
    printSteadyState = true;
}

std::string Config::trim(const std::string& str) {
    size_t first = str.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return "";
    size_t last = str.find_last_not_of(" \t\r\n");
    return str.substr(first, last - first + 1);
}

bool Config::parseBool(const std::string& value) {
    std::string lower = value;
    std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
    return (lower == "true" || lower == "1" || lower == "yes" || lower == "on");
}

int Config::parseInt(const std::string& value) {
    try {
        return std::stoi(value);
    } catch (...) {
        return 2; // Default to 2 points
    }
}

void Config::parseLine(const std::string& line) {
    // Skip empty lines and comments
    std::string trimmedLine = trim(line);
    if (trimmedLine.empty() || trimmedLine[0] == '#') {
        return;
    }
    
    // Find the '=' separator
    size_t equalPos = trimmedLine.find('=');
    if (equalPos == std::string::npos) {
        return;
    }
    
    std::string key = trim(trimmedLine.substr(0, equalPos));
    std::string value = trim(trimmedLine.substr(equalPos + 1));
    
    // Parse Gaussian integration points
    if (key == "gauss_points_h" || key == "gauss_h") {
        gaussPointsH = parseInt(value);
        if (gaussPointsH < 2 || gaussPointsH > 4) gaussPointsH = 2;
    }
    else if (key == "gauss_points_c" || key == "gauss_c") {
        gaussPointsC = parseInt(value);
        if (gaussPointsC < 2 || gaussPointsC > 4) gaussPointsC = 2;
    }
    else if (key == "gauss_points_boundary" || key == "gauss_boundary") {
        gaussPointsBoundary = parseInt(value);
        if (gaussPointsBoundary < 2 || gaussPointsBoundary > 4) gaussPointsBoundary = 2;
    }
    // Parse output flags
    else if (key == "print_local_matrices") {
        printLocalMatrices = parseBool(value);
    }
    else if (key == "print_global_matrices") {
        printGlobalMatrices = parseBool(value);
    }
    else if (key == "print_assembly_mapping") {
        printAssemblyMapping = parseBool(value);
    }
    else if (key == "print_system_info") {
        printSystemInfo = parseBool(value);
    }
    else if (key == "print_boundary_conditions") {
        printBoundaryConditions = parseBool(value);
    }
    else if (key == "print_time_steps") {
        printTimeSteps = parseBool(value);
    }
    else if (key == "print_final_results") {
        printFinalResults = parseBool(value);
    }
    else if (key == "print_steady_state") {
        printSteadyState = parseBool(value);
    }
}

bool Config::loadFromFile(const std::string& configFile) {
    configFilePath = configFile;
    std::ifstream file(configFile);
    
    if (!file.is_open()) {
        std::cerr << "Warning: Cannot open " << configFile << "\n";
        std::cerr << "Using defaults\n";
        return false;
    }
    
    std::string line;
    while (std::getline(file, line)) {
        parseLine(line);
    }
    
    file.close();
    return true;
}

void Config::print() const {
    std::cout << "\nCONFIGURATION\n";
    std::cout << "Gauss points: H=" << gaussPointsH << ", C=" << gaussPointsC 
              << ", BC=" << gaussPointsBoundary << "\n";
    std::cout << "Output: "
              << (printLocalMatrices ? "Local " : "")
              << (printGlobalMatrices ? "Global " : "")
              << (printAssemblyMapping ? "Mapping " : "")
              << (printSystemInfo ? "SysInfo " : "")
              << (printBoundaryConditions ? "BC " : "")
              << (printTimeSteps ? "TimeSteps " : "")
              << (printFinalResults ? "Final " : "")
              << (printSteadyState ? "Steady" : "") << "\n";
    std::cout << std::string(40, '=') << "\n\n";
}
