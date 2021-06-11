#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <chrono>

#include <CL/cl.hpp>

#include "lodepng.h"

std::vector<unsigned char> sequentialImplementation(const std::vector<unsigned char>& img1, const std::vector<unsigned char>& img2)
{
    std::vector<unsigned char> res;
    // Последовательная реаизация алгоритма
    return res;
}

std::vector<unsigned char> vectorizedImplementation(const std::vector<unsigned char>& img1, const std::vector<unsigned char>& img2)
{
    std::vector<unsigned char> res;
    // Реализация алгоритма с использованием векторизации
    return res;
}

std::vector<unsigned char> openmpImplementation(const std::vector<unsigned char>& img1, const std::vector<unsigned char>& img2)
{
    std::vector<unsigned char> res;
    // Реализация алгоритма с использованием OpenMP
    return res;
}

std::vector<unsigned char> openclImplementation(const std::vector<unsigned char>& img1, const std::vector<unsigned char>& img2,
                                                const cl::Context& context, const cl::Device& device, const cl::Program &program)
{
    std::vector<unsigned char> res(img1.size());
    // Реализация алгоритма с использованием OpenCL
    return res;
}

std::string readCLSources(const std::string& fileName)
{
    std::stringstream buffer;
    std::ifstream f(fileName);
    if(f.is_open())
    {
        buffer <<  f.rdbuf();
        f.close();
    }

    return buffer.str();
}


int main()
{
    std::string file1("img1.png"), file2("img2.png");

    unsigned int width, height;

    std::vector<unsigned char> inImg1, inImg2;

    auto error = lodepng::decode(inImg1, width, height, file1.c_str(), LCT_RGB);
    if(error)
    {
        std::cout << "inImg1 decode error. Error code:" << error << std::endl;
        return 1;
    }

    error = lodepng::decode(inImg2, width, height, file2.c_str(), LCT_RGB);
    if(error)
    {
        std::cout << "inImg2 decode error. Error code:" << error << std::endl;
        return 1;
    }

    const uint64_t amountOfRuns = 1000;

    long seqAvg = 0;
    long vecAvg = 0;
    long openmpAvg = 0;
    long openclAvg = 0;

    std::chrono::high_resolution_clock::time_point t1, t2;

    for(uint64_t i = 0; i < amountOfRuns; ++i)
    {
        t1 = std::chrono::high_resolution_clock::now();
        std::vector<unsigned char> seqImg = sequentialImplementation(inImg1, inImg2);
        t2 = std::chrono::high_resolution_clock::now();
        seqAvg += std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
        lodepng::encode("outseq.png", seqImg, width, height, LCT_RGB);

        t1 = std::chrono::high_resolution_clock::now();
        std::vector<unsigned char> vecImg = vectorizedImplementation(inImg1, inImg2);
        t2 = std::chrono::high_resolution_clock::now();
        vecAvg += std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
        lodepng::encode("outvec.png", vecImg, width, height, LCT_RGB);

        t1 = std::chrono::high_resolution_clock::now();
        std::vector<unsigned char> openmpImg = openmpImplementation(inImg1, inImg2);
        t2 = std::chrono::high_resolution_clock::now();
        openmpAvg += std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
        lodepng::encode("outopenmp.png", vecImg, width, height, LCT_RGB);

        t1 = std::chrono::high_resolution_clock::now();
        std::vector<unsigned char> openclImg = openclImplementation(inImg1, inImg2, context, default_device, program);
        t2 = std::chrono::high_resolution_clock::now();
        openclAvg += std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
        lodepng::encode("outopencl.png", vecImg, width, height, LCT_RGB);
    }

    //Вычисление среднего времени выполненине.

    return 0;
}
