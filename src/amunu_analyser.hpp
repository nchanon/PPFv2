#pragma once 

#include <vector>
#include <array>
#include <string>

std::vector<std::array<std::array<double, 4>, 4> > readVecMatrix(std::string const& s, int n);

void CreateTH1(std::vector<std::array<std::array<double, 4>, 4> > const& vm,
               std::string const& name, int nbin, int binMin, int binMax);

bool isZero(std::array<std::array<double, 4>, 4> const& matrix);

void CreateComparaison(std::vector<std::array<std::array<double, 4>, 4> > const& matrix, std::string const& wilson);