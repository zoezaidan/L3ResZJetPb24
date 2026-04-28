#pragma once

#include "TMath.h"

#include <cmath>
#include <iostream>
#include <sstream>
#include <string>

#define RESET_COLOR "\033[0m"
#define RED_COLOR "\033[91m"
#define GREEN_COLOR "\033[92m"
#define YELLOW_COLOR "\033[93m"
#define BLUE_COLOR "\033[94m"
#define MAGENTA_COLOR "\033[95m"
#define CYAN_COLOR "\033[96m"
#define WHITE_COLOR "\033[97m"
#define BOLD_COLOR "\033[1m"

inline int g_verbosity = 2;

enum LogLevel {
  LOG_ERROR = 0,
  LOG_WARNING = 1,
  LOG_INFO = 2,
  LOG_DEBUG = 3,
  LOG_TRACE = 4
};

inline void log(LogLevel level, const std::string &message) {
  if (static_cast<int>(level) > g_verbosity)
    return;

  const char *colorCode = RESET_COLOR;
  const char *levelStr = "LOG";

  switch (level) {
  case LOG_ERROR:
    colorCode = RED_COLOR;
    levelStr = "ERROR";
    break;
  case LOG_WARNING:
    colorCode = YELLOW_COLOR;
    levelStr = "WARNING";
    break;
  case LOG_INFO:
    colorCode = GREEN_COLOR;
    levelStr = "INFO";
    break;
  case LOG_DEBUG:
    colorCode = MAGENTA_COLOR;
    levelStr = "DEBUG";
    break;
  case LOG_TRACE:
    colorCode = CYAN_COLOR;
    levelStr = "TRACE";
    break;
  }

  std::cout << colorCode << "[" << levelStr << "] " << message << RESET_COLOR
            << std::endl;
}

inline void log_progress_every(long long currentEvent, long long totalEvents,
                               long long interval = 1000,
                               LogLevel level = LOG_INFO) {
  if (totalEvents <= 0 || currentEvent <= 0 || interval <= 0)
    return;

  if (currentEvent % interval != 0 && currentEvent != totalEvents)
    return;

  const double percent = 100.0 * static_cast<double>(currentEvent) /
                         static_cast<double>(totalEvents);
  std::ostringstream progress;
  progress << "Processed event " << currentEvent << "/" << totalEvents << " ("
           << percent << "%)";
  log(level, progress.str());
}

inline double DPhi(double phi1, double phi2) {
  double dphi = std::fabs(phi1 - phi2);
  return (dphi <= TMath::Pi()) ? dphi : TMath::TwoPi() - dphi;
}
