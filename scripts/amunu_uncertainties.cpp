#include <iostream>
#include <fstream>
#include <array>
#include <vector>
#include <cmath>

using Matrix = std::array<std::array<float, 4>, 4>;

void Print(Matrix const& m)
{
    for(size_t j = 0; j < 4; ++j){
        for(size_t k = 0; k < 4; ++k){
                std::cout << m[j][k] << " ";
        }
        std::cout << std::endl;
    }
}

void Lorentz(Matrix & m)
{
    std::array<std::array<float, 4>, 4> foo = m;
    m[0][0] = foo[3][3];
    for(size_t i = 1; i <= 3; ++i){
        m[0][i] = foo[i-1][3];
        m[i][0] = foo[3][i-1];
        for(size_t j = 1; j <= 3; ++j){
            m[i][j] = foo[i-1][j-1];
        }
    }
}

bool nonZero(Matrix const& m)
{
    bool foo = false;
    for(size_t j = 0; j < 4; ++j){
        for(size_t k = 0; k < 4; ++k){
            if(m[j][k] != 0)
                foo = true;
        }
    } 
    return foo;
}

Matrix average(std::vector<Matrix> const& m, float size)
{
    Matrix foo;
    for(size_t i = 0; i < m.size(); ++i){
        for(size_t j = 0; j < 4; ++j){
            for(size_t k = 0; k < 4; ++k){
                foo[j][k] += m[i][j][k];
            }
        }
    }
    for(size_t j = 0; j < 4; ++j){
        for(size_t k = 0; k < 4; ++k){
            foo[j][k] /= size;
        }
    }
    return foo;
}

Matrix stdDev(std::vector<Matrix> const& m, Matrix const& av, float size)
{
    Matrix foo;
    for(size_t i = 0; i < m.size(); ++i){
        for(size_t j = 0; j < 4; ++j){
            for(size_t k = 0; k < 4; ++k){
                if(m[i][j][k] != 0)
                    foo[j][k] += std::pow(m[i][j][k] - av[j][k], 2);
            }
        }
    }
    for(size_t j = 0; j < 4; ++j){
        for(size_t k = 0; k < 4; ++k){
            foo[j][k] /= size;
            foo[j][k] = std::sqrt(foo[j][k]);
            foo[j][k] /= std::sqrt(size);
        }
    }
    return foo;
}

int main(int argc, char* argv[])
{
    int n = 5000000;
    std::string input;

    if(argc == 2){
        std::string option = argv[1];
        if(option == "f"){
            input = "VMAF13TeVCMS.txt";
        }
        else if (option == "gg"){   
            input = "VMAP2g13TeVCMS.txt";
        }  
        else if (option == "qq"){
            input = "VMAPqqbar13TeVCMS.txt";
        }
        else if (option == "test"){
            input = "VMAF100TeVCMS.txt";
            n = 1000000;
        }
        else{
            std::cout << "Wrong argument (\"f\", \"gg\", \"qq\")" << std::endl;
            return 0;
        }
    } 
    else{
        std::cout << "No good number of arguments"<< std::endl;
        return 0;
    }

    std::vector<Matrix> matrix(n);

//--------------------------------//

    std::cout << "Get inputs" << std::endl;

    std::ifstream file("./inputs/pheno/VecMatrix/"+input);
    
    for(int i=0; i<n; ++i){
        for(size_t j = 0; j < 4; ++j){
            for(size_t k = 0; k < 4; ++k){
                file >> matrix[i][j][k];
            }
        }
        if(i%1000000 == 0)
            std::cout << "--> 1 millions events passed" << std::endl;
    }
    std::cout << std::endl;

//--------------------------------//

    std::cout << "Number of event" << std::endl;
    float nevent = 0;
    for(int i=0; i<n; ++i){
        if(nonZero(matrix[i]))
            nevent++;
    }
    std::cout << "--> n = " << nevent << std::endl;
    std::cout << std::endl;

//--------------------------------//

    Matrix matrixAverage = average(matrix, nevent);
    Matrix matrixDeviation = stdDev(matrix, matrixAverage, nevent);

    std::cout << "Average" << std::endl;
    std::cout << std::endl;
    Lorentz(matrixAverage);
    Print(matrixAverage);
    std::cout << std::endl;
    
    std::cout << "Standard deviation" << std::endl;
    std::cout << std::endl;
    Lorentz(matrixDeviation);
    Print(matrixDeviation);
    std::cout << std::endl;

    return 0;
}