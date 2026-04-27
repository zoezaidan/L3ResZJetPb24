#ifndef __input_config_h__
#define __input_config_h__

#include <algorithm>
#include <dirent.h>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <sys/stat.h>
#include <vector>

// Input type enumeration
enum class InputType {
  FILE,      // Single ROOT file
  DIRECTORY, // Directory containing ROOT files
  FILELIST   // Text file with list of ROOT file paths
};

// Input configuration structure
struct InputConfig {
  InputType type;
  std::string path;
  int maxFiles;  // -1 for all files
  int skipFiles; // Number of files to skip (for batching)
  int maxEvents; // -1 for all events
  std::string outputDir;
  std::string outputTag;
  int batchIndex;   // -1 for non-batch mode
  int totalBatches; // Total number of batches
};

// Get list of ROOT files from directory (recursive)
inline void GetFilesFromDirectory(const std::string &dirPath,
                                  std::vector<std::string> &files,
                                  int maxFiles = -1) {
  DIR *dir = opendir(dirPath.c_str());
  if (!dir) {
    std::cerr << "ERROR: Cannot open directory: " << dirPath << std::endl;
    return;
  }

  struct dirent *entry;
  while ((entry = readdir(dir)) != nullptr) {
    std::string name = entry->d_name;
    if (name == "." || name == "..")
      continue;

    std::string fullPath = dirPath + "/" + name;
    struct stat statbuf;
    if (stat(fullPath.c_str(), &statbuf) == 0) {
      if (S_ISDIR(statbuf.st_mode)) {
        // Recursively search subdirectories
        GetFilesFromDirectory(fullPath, files, maxFiles);
      } else if (name.size() > 5 && name.substr(name.size() - 5) == ".root") {
        files.push_back(fullPath);
        if (maxFiles > 0 && (int)files.size() >= maxFiles) {
          closedir(dir);
          return;
        }
      }
    }
  }
  closedir(dir);
}

// Get list of ROOT files from filelist text file
inline void GetFilesFromFilelist(const std::string &listPath,
                                 std::vector<std::string> &files,
                                 int maxFiles = -1) {
  std::ifstream infile(listPath);
  if (!infile.is_open()) {
    std::cerr << "ERROR: Cannot open filelist: " << listPath << std::endl;
    return;
  }

  std::string line;
  while (std::getline(infile, line)) {
    // Skip empty lines and comments
    if (line.empty() || line[0] == '#')
      continue;
    // Trim whitespace
    size_t start = line.find_first_not_of(" \t");
    size_t end = line.find_last_not_of(" \t");
    if (start == std::string::npos)
      continue;
    line = line.substr(start, end - start + 1);

    files.push_back(line);
    if (maxFiles > 0 && (int)files.size() >= maxFiles)
      break;
  }
  infile.close();
}

// Main function to get files based on input type
inline std::vector<std::string> GetInputFiles(const InputConfig &config) {
  std::vector<std::string> files;

  switch (config.type) {
  case InputType::FILE:
    files.push_back(config.path);
    break;
  case InputType::DIRECTORY:
    GetFilesFromDirectory(config.path, files, -1); // Get all first
    break;
  case InputType::FILELIST:
    GetFilesFromFilelist(config.path, files, -1); // Get all first
    break;
  }

  // Sort for reproducibility
  std::sort(files.begin(), files.end());

  // Apply skip for batch mode
  if (config.skipFiles > 0 && config.skipFiles < (int)files.size()) {
    files.erase(files.begin(), files.begin() + config.skipFiles);
  }

  // Apply maxFiles limit
  if (config.maxFiles > 0 && (int)files.size() > config.maxFiles) {
    files.resize(config.maxFiles);
  }

  return files;
}

// Auto-detect input type from path
inline InputType DetectInputType(const std::string &path) {
  struct stat statbuf;
  if (stat(path.c_str(), &statbuf) == 0) {
    if (S_ISDIR(statbuf.st_mode)) {
      return InputType::DIRECTORY;
    }
  }
  // Check file extension
  if (path.size() > 4) {
    std::string ext = path.substr(path.size() - 4);
    if (ext == ".txt" || path.substr(path.size() - 5) == ".list") {
      return InputType::FILELIST;
    }
  }
  return InputType::FILE;
}

#endif
