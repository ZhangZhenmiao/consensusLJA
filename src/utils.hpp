#ifndef COMMAND_UTILS_H
#define COMMAND_UTILS_H

#include <string>

int execute_command(const std::string& command);

double unialigner_identity(const std::string& seq1, const std::string& seq2, int prefix = 1000000);

std::string getExecutablePath();

#endif // COMMAND_UTILS_H