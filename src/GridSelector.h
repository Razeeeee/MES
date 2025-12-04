#pragma once
#include <string>
#include <vector>
#include <filesystem>

/**
 * @brief Class responsible for handling grid file selection
 * 
 * Manages the display of available grid files, user input,
 * and provides the selected grid filename
 */
class GridSelector {
private:
    std::string gridDirectory;
    std::vector<std::string> availableGridFiles;
    
    // Helper methods
    void scanGridDirectory();
    bool isValidGridFile(const std::filesystem::path& filePath) const;
    
public:
    GridSelector(const std::string& gridDir = "grids/");
    
    // Main functionality
    void displayMenu() const;
    std::string selectGrid();
    bool hasGridFiles() const;
    
    // Getters
    const std::vector<std::string>& getAvailableFiles() const noexcept { return availableGridFiles; }
    const std::string& getGridDirectory() const noexcept { return gridDirectory; }
    
    // Utility methods
    void refreshFileList();
    std::string getGridPath(int choice) const;
};
