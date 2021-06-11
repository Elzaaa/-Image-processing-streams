#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <chrono>
#include "SimpleTimer.h"
#include <CL/cl.hpp>

#include "lodepng.h"

// opencl
#define VECTOR_SIZE 625 * 625 * 3
#define LOCAL_VECTOR_SIZE 5

SimpleTimer GlobalTimer;

std::vector<unsigned char> sequentialImplementation(const std::vector<unsigned char>& img1, const std::vector<unsigned char>& img2, unsigned int w, unsigned int h)
{
    std::vector<unsigned char> res(h * w * 3);
    // Последовательная реализация алгоритма
    for (unsigned y = 0; y < h; y++)
        for (unsigned x = 0; x < w; x++) {
            // get RGB components
            res[3 * y * w + 3 * x + 0] = (((int)img1[3 * y * w + 3 * x + 0] - (int)img2[3 * y * w + 3 * x + 0]) > 0) ? img1[3 * y * w + 3 * x + 0] - img2[3 * y * w + 3 * x + 0] : 0; //red
            res[3 * y * w + 3 * x + 1] = (((int)img1[3 * y * w + 3 * x + 1] - (int)img2[3 * y * w + 3 * x + 1]) > 0) ? img1[3 * y * w + 3 * x + 1] - img2[3 * y * w + 3 * x + 1] : 0; //green
            res[3 * y * w + 3 * x + 2] = (((int)img1[3 * y * w + 3 * x + 2] - (int)img2[3 * y * w + 3 * x + 2]) > 0) ? img1[3 * y * w + 3 * x + 2] - img2[3 * y * w + 3 * x + 2] : 0; //blue
        }
    return res;
}

std::vector<unsigned char> vectorizedImplementation(const std::vector<unsigned char>& img1, const std::vector<unsigned char>& img2, unsigned int w, unsigned int h)
{
    std::vector<unsigned char> res(h * w * 3);
    // Реализация алгоритма с использованием векторизации

    __m256i _img1, _img2, _res, _sub, _zero, _const256 = _mm256_set1_epi16(0xff);
    _zero = _mm256_set1_epi16(0);
    const __m256i;
    //part == часть
    int part = res.size() / 32 * 32;
    for (int i = 0; i < part; i += 32) {
        _img1 = _mm256_and_si256(_mm256_loadu_epi16(&img1[i]), _const256);
        _img2 = _mm256_and_si256(_mm256_loadu_epi16(&img2[i]), _const256);
        _res = _mm256_subs_epu16(_img1, _img2);


        for (int j = 0; j < 32; j += 2) res[i + j] = _res.m256i_u16[j / 2];

        _img1 = _mm256_and_si256(_mm256_loadu_epi16(&img1[i + 1]), _const256);
        _img2 = _mm256_and_si256(_mm256_loadu_epi16(&img2[i + 1]), _const256);
        _res = _mm256_subs_epu16(_img1, _img2);


        for (int j = 1; j < 32; j += 2) res[i + j] = _res.m256i_u16[j / 2];
    }
    if (res.size() - part) {
        _img1 = _mm256_and_si256(_mm256_loadu_epi16(&img1[part]), _const256);
        _img2 = _mm256_and_si256(_mm256_loadu_epi16(&img2[part]), _const256);
        _res = _mm256_subs_epu16(_img1, _img2);


        for (int j = 0; j < res.size() - part; j += 2) res[part + j] = _res.m256i_u16[j / 2];

        _img1 = _mm256_and_si256(_mm256_loadu_epi16(&img1[part + 1]), _const256);
        _img2 = _mm256_and_si256(_mm256_loadu_epi16(&img2[part + 1]), _const256);
        _res = _mm256_subs_epu16(_img1, _img2);


        for (int j = 1; j < res.size() - part; j += 2) res[part + j] = _res.m256i_u16[j / 2];
    }

    return res;
}

std::vector<unsigned char> openmpImplementation(const std::vector<unsigned char>& img1, const std::vector<unsigned char>& img2, unsigned int w, unsigned int h)
{
    std::vector<unsigned char> res(h * w * 3);
    // Реализация алгоритма с использованием OpenMP

#pragma omp parallel for shared(img1, img2, res)
    for (size_t i = 0; i < res.size(); i++)
        res[i] = (((int)img1[i] - (int)img2[i]) > 0) ? img1[i] - img2[i] : 0;
    return res;
}

std::vector<unsigned char> openclImplementation(const std::vector<unsigned char>& img1, const std::vector<unsigned char>& img2,
    cl_kernel& kernel, cl_command_queue& cq, cl_mem& res_clmem)
{
    std::vector<unsigned char> res(img1.size());
    // Реализация алгоритма с использованием OpenCL
    size_t global_size = VECTOR_SIZE; // Process the entire lists // Обработать все списки
    size_t local_size = LOCAL_VECTOR_SIZE; // Process one item at a time // Обрабатывать по одному элементу за раз
    clEnqueueNDRangeKernel(cq, kernel, 1, NULL, &global_size, &local_size, 0, NULL, NULL);

    // Read the cl memory res_clmem on device to the host variable res
    //Прочтите cl memory res_clmem на устройстве в переменную хоста res
    clEnqueueReadBuffer(cq, res_clmem, CL_TRUE, 0, VECTOR_SIZE * sizeof(unsigned char), &res[0], 0, NULL, NULL);

    //Clean upand wait for all the comands to complete.
    //Очистите и дождитесь завершения всех команд.
    clFlush(cq);
    clFinish(cq);

    return res;
}

std::string readCLSources(const std::string& fileName)
{
    std::stringstream buffer;
    std::ifstream f(fileName);
    if (f.is_open())
    {
        buffer << f.rdbuf();
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
    if (error)
    {
        std::cout << "inImg1 decode error. Error code:" << error << std::endl;
        return 1;
    }

    error = lodepng::decode(inImg2, width, height, file2.c_str(), LCT_RGB);
    if (error)
    {
        std::cout << "inImg2 decode error. Error code:" << error << std::endl;
        return 1;
    }

    const uint64_t amountOfRuns = 1001;

    SimpleTimer seqTimer;
    for (uint64_t i = 0; i < amountOfRuns; i++) {
        std::vector<unsigned char> seqImg = sequentialImplementation(inImg1, inImg2, width, height);
        if (i == 0)
            lodepng::encode("outoseq.png", seqImg, width, height, LCT_RGB);
    }
    seqTimer.StopTime();

    SimpleTimer vecTimer;
    for (uint64_t i = 0; i < amountOfRuns; i++)
    {
        std::vector<unsigned char> vecImg = vectorizedImplementation(inImg1, inImg2, width, height);
        if (i == 0)
            lodepng::encode("outovec.png", vecImg, width, height, LCT_RGB);
    }
    vecTimer.StopTime();

    SimpleTimer openmpTimer;
    for (uint64_t i = 0; i < amountOfRuns; i++)
    {
        std::vector<unsigned char> openmpImg = openmpImplementation(inImg1, inImg2, width, height);
        if (i == 0)
            lodepng::encode("outopenmp.png", openmpImg, width, height, LCT_RGB);

    }
    openmpTimer.StopTime();




    // OPENCL
    //OpenCL kernel which is run for every work item created.
    const char* saxpy_kernel =
        "__kernel                                   \n"
        "void saxpy_kernel(                         \n"
        "           __global unsigned char *img1,       \n"
        "           __global unsigned char *img2,       \n"
        "           __global unsigned char *res)       \n"
        "{                                          \n"
        "    //Get the index of the work-item       \n"
        "    int index = get_global_id(0);          \n"
        "    res[index] = (((int)img1[index] - (int)img2[index]) > 0) ? img1[index] - img2[index] : 0; \n"
        "}                                          \n";

    // Get platform and device information
    cl_platform_id* platforms = NULL;
    cl_uint     num_platforms;
    // Set up the Platform
    cl_int clStatus = clGetPlatformIDs(0, NULL, &num_platforms);
    platforms = (cl_platform_id*)
        malloc(sizeof(cl_platform_id) * num_platforms);
    clStatus = clGetPlatformIDs(num_platforms, platforms, NULL);

    // Get the devices list and choose the device you want to run on
    cl_device_id* device_list = NULL;
    cl_uint           num_devices;

    clStatus = clGetDeviceIDs(platforms[0], CL_DEVICE_TYPE_GPU, 0, NULL, &num_devices);
    device_list = (cl_device_id*)malloc(sizeof(cl_device_id) * num_devices);
    clStatus = clGetDeviceIDs(platforms[0], CL_DEVICE_TYPE_GPU, num_devices, device_list, NULL);

    // Create one OpenCL context for each device in the platform
    cl_context context;
    context = clCreateContext(NULL, num_devices, device_list, NULL, NULL, &clStatus);

    // Create a command queue
    cl_command_queue command_queue = clCreateCommandQueue(context, device_list[0], 0, &clStatus);

    // Create memory buffers on the device for each vector
    cl_mem img1_clmem = clCreateBuffer(context, CL_MEM_READ_ONLY, VECTOR_SIZE * sizeof(unsigned char), NULL, &clStatus);
    cl_mem img2_clmem = clCreateBuffer(context, CL_MEM_READ_ONLY, VECTOR_SIZE * sizeof(unsigned char), NULL, &clStatus);
    cl_mem res_clmem = clCreateBuffer(context, CL_MEM_WRITE_ONLY, VECTOR_SIZE * sizeof(unsigned char), NULL, &clStatus);

    // Copy the Buffer img1 and img2 to the device
    clStatus = clEnqueueWriteBuffer(command_queue, img1_clmem, CL_TRUE, 0, VECTOR_SIZE * sizeof(unsigned char), &inImg1[0], 0, NULL, NULL);
    clStatus = clEnqueueWriteBuffer(command_queue, img2_clmem, CL_TRUE, 0, VECTOR_SIZE * sizeof(unsigned char), &inImg2[0], 0, NULL, NULL);

    // Create a program from the kernel source
    cl_program program = clCreateProgramWithSource(context, 1, (const char**)&saxpy_kernel, NULL, &clStatus);

    // Build the program
    clStatus = clBuildProgram(program, 1, device_list, NULL, NULL, NULL);

    // Create the OpenCL kernel
    cl_kernel kernel = clCreateKernel(program, "saxpy_kernel", &clStatus);

    // Set the arguments of the kernel
    clStatus = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void*)&img1_clmem);
    clStatus = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void*)&img2_clmem);
    clStatus = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void*)&res_clmem);


    SimpleTimer openclTimer;
    for (uint64_t i = 0; i < amountOfRuns; i++)
    {
        std::vector<unsigned char> openclImg = openclImplementation(inImg1, inImg2, kernel, command_queue, res_clmem);
        if (i == 0)
            lodepng::encode("outopencl.png", openclImg, width, height, LCT_RGB);

    }
    openclTimer.StopTime();

    // Finally release all OpenCL allocated objects and host buffers.
    clStatus = clReleaseKernel(kernel);
    clStatus = clReleaseProgram(program);
    clStatus = clReleaseMemObject(img1_clmem);
    clStatus = clReleaseMemObject(img2_clmem);
    clStatus = clReleaseMemObject(res_clmem);
    clStatus = clReleaseCommandQueue(command_queue);
    clStatus = clReleaseContext(context);
    free(platforms);
    free(device_list);


    //Вычисление среднего времени выполнениня.
    std::cout << "\nRuns done:\t " << amountOfRuns << "\n\n";

    std::cout << "sequentialImplementation Average Time:\t " << seqTimer.GetTime() / amountOfRuns << "s" << std::endl;
    std::cout << "vectorizedImplementation Average Time:\t " << vecTimer.GetTime() / amountOfRuns << "s" << std::endl;
    std::cout << "openMPImplementation Average Time:\t " << openmpTimer.GetTime() / amountOfRuns << "s" << std::endl;
    std::cout << "openCLImplementation Average Time:\t " << openclTimer.GetTime() / amountOfRuns << "s" << std::endl;

    std::cout << "\nWhole program:\t " << GlobalTimer.GetTime() << "s" << std::endl;
    return 0;
}
