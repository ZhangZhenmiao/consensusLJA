#include "utils.hpp"
#include <cstdlib>
#include <fstream>
#include <unistd.h>
#include <sys/wait.h>
#include <fcntl.h>
#include <iostream>
#include <limits.h>

int execute_command(const std::string& command) {
    struct CommandResult {
        bool success;
        int exit_code;
    };
    CommandResult result;
    result.success = false;
    result.exit_code = -1;

    result.exit_code = system(command.c_str());

    // Check if command executed successfully
    if (WIFEXITED(result.exit_code)) {
        result.exit_code = WEXITSTATUS(result.exit_code);
        result.success = (result.exit_code == 0);
    }

    std::string status;

    status += "Command " + command + " " + std::string(result.success ? "succeeded" : "failed");
    status += " with exit code " + std::to_string(result.exit_code) + "\n";
    // std::cout << status << std::endl;

    if (!result.success)
        std::cout << "Failed to execute " + command << std::endl;

    return result.exit_code;
}

std::string getExecutablePath() {
    char buffer[1024];
    ssize_t len = readlink("/proc/self/exe", buffer, sizeof(buffer) - 1);
    if (len != -1) {
        buffer[len] = '\0';
        std::string execPath = std::string(buffer);
        return execPath.substr(0, execPath.find_last_of("/"));
    }
    return "";
}