#pragma once

#include <string>
#include <ctime>
#include <iostream>

void Log(double value);
void Log(std::string const& message);
void Log(int argc, char** argv);

double elapsedTime(std::clock_t t0, std::clock_t t1);