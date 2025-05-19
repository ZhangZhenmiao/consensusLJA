#ifndef COMMAND_UTILS_H
#define COMMAND_UTILS_H

#include <string>

int execute_command(const std::string& command);

std::string getExecutablePath();

#endif // COMMAND_UTILS_H