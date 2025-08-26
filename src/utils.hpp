#ifndef COMMAND_UTILS_H
#define COMMAND_UTILS_H

#include <string>

int execute_command(const std::string& command, bool exit_when_fail = true);

double unialigner_identity(const std::string& seq1, const std::string& seq2, int prefix = 1000000);

std::string getExecutablePath();

std::string replaceNsWithRandomBases(const std::string& seq);

std::string format_with_commas(int value);

std::pair<double, double> calculate_identities_from_cigar(const std::string& cigar, int gap_threshold = 10);

#endif // COMMAND_UTILS_H