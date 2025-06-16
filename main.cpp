#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <chrono>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <random>
#include <sstream>
#include <ctime>
#include <limits>
#include <array>
#include <cpuid.h>
#include <zlib.h>
#include <openssl/evp.h>
#include <openssl/err.h>
#include <openssl/rand.h>
#include <cstring>

// Constants and configurations
const std::string CSV_FILENAME = "benchmark_results.csv";
const size_t FLOAT_ARRAY_SIZE = 1000000;
const int ITERATIONS = 10;

// Benchmark types
enum BenchmarkType {
    BENCH_ZLIB,
    BENCH_OPENSSL,
    BENCH_FLOAT,
    BENCH_COUNT
};

const char* benchmarkNames[BENCH_COUNT] = {
    "ZLIB Compression",
    "OpenSSL Encryption",
    "Floating Point Operations"
};

// ZLIB compression methods
enum CompressionMethod {
    ZLIB_DEFAULT,
    ZLIB_BEST_SPEED,
    ZLIB_BEST_COMPRESSION,
    ZLIB_NO_COMPRESSION,
    ZLIB_HUFFMAN_ONLY,
    ZLIB_RLE,
    GZIP_DEFAULT,
    ZLIB_METHOD_COUNT
};

const char* zlibMethodNames[ZLIB_METHOD_COUNT] = {
    "ZLIB Default",
    "ZLIB Best Speed",
    "ZLIB Best Compression",
    "ZLIB No Compression",
    "ZLIB Huffman Only",
    "ZLIB RLE",
    "GZIP Default"
};

// OpenSSL encryption methods
enum EncryptionMethod {
    AES_128_CBC,
    AES_256_CBC,
    AES_128_GCM,
    AES_256_GCM,
    CHACHA20,
    OPENSSL_METHOD_COUNT
};

const char* opensslMethodNames[OPENSSL_METHOD_COUNT] = {
    "AES-128-CBC",
    "AES-256-CBC",
    "AES-128-GCM",
    "AES-256-GCM",
    "ChaCha20"
};

// Floating point operations
enum OperationType {
    OP_ADD,
    OP_SUB,
    OP_MUL,
    OP_DIV,
    OP_FMA,
    OP_SQRT,
    OP_TRIG,
    OP_LOG,
    OP_MIXED,
    OP_COUNT
};

const char* floatOpNames[OP_COUNT] = {
    "Addition",
    "Subtraction",
    "Multiplication",
    "Division",
    "Fused Multiply-Add",
    "Square Root",
    "Trigonometric",
    "Logarithmic",
    "Mixed Operations"
};

// Utility functions
std::string getCurrentTimestamp() {
    auto now = std::chrono::system_clock::now();
    auto in_time_t = std::chrono::system_clock::to_time_t(now);

    std::stringstream ss;
    ss << std::put_time(std::localtime(&in_time_t), "%Y-%m-%d %X");
    return ss.str();
}

std::string getCPUModel() {
    char CPUBrandString[0x40];
    unsigned int CPUInfo[4] = {0, 0, 0, 0};

    __cpuid(0x80000002, CPUInfo[0], CPUInfo[1], CPUInfo[2], CPUInfo[3]);
    std::memcpy(CPUBrandString, CPUInfo, sizeof(CPUInfo));
    __cpuid(0x80000003, CPUInfo[0], CPUInfo[1], CPUInfo[2], CPUInfo[3]);
    std::memcpy(CPUBrandString + 16, CPUInfo, sizeof(CPUInfo));
    __cpuid(0x80000004, CPUInfo[0], CPUInfo[1], CPUInfo[2], CPUInfo[3]);
    std::memcpy(CPUBrandString + 32, CPUInfo, sizeof(CPUInfo));

    return std::string(CPUBrandString);
}

void writeCSVHeader(std::ofstream& csvFile) {
    csvFile << "Timestamp,CPU Model,Benchmark Type,Operation/Method,"
            << "Time (ns),Operations,Size Before (bytes),Size After (bytes),"
            << "Ratio (%),Notes\n";
}

void appendToCSV(std::ofstream& csvFile, 
                const std::string& benchmarkType,
                const std::string& operation,
                uint64_t timeNs,
                size_t operations = 0,
                size_t sizeBefore = 0,
                size_t sizeAfter = 0,
                const std::string& notes = "") {
    csvFile << getCurrentTimestamp() << ","
            << "\"" << getCPUModel() << "\","
            << "\"" << benchmarkType << "\","
            << "\"" << operation << "\","
            << timeNs << ","
            << operations << ","
            << sizeBefore << ","
            << sizeAfter << ","
            << (sizeBefore > 0 ? std::to_string((double)sizeAfter / sizeBefore * 100.0) : "0") << ","
            << "\"" << notes << "\"\n";
}

// ZLIB benchmark functions
std::vector<std::string> readFileByLine(const std::string& filename) {
    std::vector<std::string> lines;
    std::ifstream file(filename);
    std::string line;
    
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        exit(1);
    }
    
    while (std::getline(file, line)) {
        lines.push_back(line);
    }
    
    file.close();
    return lines;
}

std::string compressData(const std::string& data, int method, int level) {
    z_stream zs;
    memset(&zs, 0, sizeof(zs));
    
    int windowBits = 15;
    if (method == GZIP_DEFAULT) {
        windowBits += 16;
    }
    
    if (deflateInit2(&zs, level, Z_DEFLATED, windowBits, 8, Z_DEFAULT_STRATEGY) != Z_OK) {
        return "";
    }
    
    zs.next_in = (Bytef*)data.data();
    zs.avail_in = data.size();
    
    int ret;
    char outbuffer[32768];
    std::string outstring;
    
    do {
        zs.next_out = reinterpret_cast<Bytef*>(outbuffer);
        zs.avail_out = sizeof(outbuffer);
        
        ret = deflate(&zs, Z_FINISH);
        
        if (outstring.size() < zs.total_out) {
            outstring.append(outbuffer, zs.total_out - outstring.size());
        }
    } while (ret == Z_OK);
    
    deflateEnd(&zs);
    
    if (ret != Z_STREAM_END) {
        return "";
    }
    
    return outstring;
}

std::string decompressData(const std::string& data, int method) {
    z_stream zs;
    memset(&zs, 0, sizeof(zs));
    
    int windowBits = 15;
    if (method == GZIP_DEFAULT) {
        windowBits += 16;
    }
    
    if (inflateInit2(&zs, windowBits) != Z_OK) {
        return "";
    }
    
    zs.next_in = (Bytef*)data.data();
    zs.avail_in = data.size();
    
    int ret;
    char outbuffer[32768];
    std::string outstring;
    
    do {
        zs.next_out = reinterpret_cast<Bytef*>(outbuffer);
        zs.avail_out = sizeof(outbuffer);
        
        ret = inflate(&zs, 0);
        
        if (outstring.size() < zs.total_out) {
            outstring.append(outbuffer, zs.total_out - outstring.size());
        }
    } while (ret == Z_OK);
    
    inflateEnd(&zs);
    
    if (ret != Z_STREAM_END) {
        return "";
    }
    
    return outstring;
}

void benchmarkZLIB(std::ofstream& csvFile, const std::string& filename) {
    std::vector<std::string> lines = readFileByLine(filename);
    std::string data;
    for (const auto& line : lines) {
        data += line + "\n";
    }
    
    std::cout << "\n=== ZLIB Compression Benchmark ===\n";
    std::cout << "Original data size: " << data.size() << " bytes\n";
    std::cout << "Number of lines: " << lines.size() << "\n";
    
    for (int i = 0; i < ZLIB_METHOD_COUNT; ++i) {
        int level = Z_DEFAULT_COMPRESSION;
        int strategy = Z_DEFAULT_STRATEGY;
        
        switch (i) {
            case ZLIB_BEST_SPEED: level = Z_BEST_SPEED; break;
            case ZLIB_BEST_COMPRESSION: level = Z_BEST_COMPRESSION; break;
            case ZLIB_NO_COMPRESSION: level = Z_NO_COMPRESSION; break;
            case ZLIB_HUFFMAN_ONLY: strategy = Z_HUFFMAN_ONLY; break;
            case ZLIB_RLE: strategy = Z_RLE; break;
        }
        
        // Compression
        auto start = std::chrono::high_resolution_clock::now();
        std::string compressed = compressData(data, i, level);
        auto end = std::chrono::high_resolution_clock::now();
        auto compressTime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
        
        if (!compressed.empty()) {
            // Decompression
            start = std::chrono::high_resolution_clock::now();
            std::string decompressed = decompressData(compressed, i);
            end = std::chrono::high_resolution_clock::now();
            auto decompressTime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
            
            if (!decompressed.empty() && decompressed == data) {
                appendToCSV(csvFile, "ZLIB", zlibMethodNames[i], compressTime, 
                          1, data.size(), compressed.size(), "Compression");
                appendToCSV(csvFile, "ZLIB", zlibMethodNames[i], decompressTime, 
                          1, compressed.size(), data.size(), "Decompression");
                
                std::cout << zlibMethodNames[i] << ": "
                          << "Compress " << compressTime / 1000000.0 << " ms, "
                          << "Decompress " << decompressTime / 1000000.0 << " ms, "
                          << "Ratio " << std::fixed << std::setprecision(2) 
                          << (double)compressed.size() / data.size() * 100.0 << "%\n";
            }
        }
    }
}

// OpenSSL benchmark functions
EVP_CIPHER_CTX* createCipherContext(int method, bool encrypt, 
                                   const unsigned char* key, 
                                   const unsigned char* iv) {
    const EVP_CIPHER* cipher = nullptr;
    
    switch (method) {
        case AES_128_CBC: cipher = EVP_aes_128_cbc(); break;
        case AES_256_CBC: cipher = EVP_aes_256_cbc(); break;
        case AES_128_GCM: cipher = EVP_aes_128_gcm(); break;
        case AES_256_GCM: cipher = EVP_aes_256_gcm(); break;
        case CHACHA20: cipher = EVP_chacha20(); break;
        default: return nullptr;
    }
    
    EVP_CIPHER_CTX* ctx = EVP_CIPHER_CTX_new();
    if (!ctx) return nullptr;
    
    int ret = encrypt ? 
        EVP_EncryptInit_ex(ctx, cipher, nullptr, key, iv) :
        EVP_DecryptInit_ex(ctx, cipher, nullptr, key, iv);
    
    if (ret != 1) {
        EVP_CIPHER_CTX_free(ctx);
        return nullptr;
    }
    
    return ctx;
}

std::string encryptData(const std::string& data, int method, 
                        const unsigned char* key, const unsigned char* iv) {
    EVP_CIPHER_CTX* ctx = createCipherContext(method, true, key, iv);
    if (!ctx) return "";
    
    int blockSize = EVP_CIPHER_CTX_block_size(ctx);
    std::string output;
    output.resize(data.size() + blockSize);
    
    int outLen = 0;
    if (EVP_EncryptUpdate(ctx, 
                         reinterpret_cast<unsigned char*>(&output[0]), 
                         &outLen,
                         reinterpret_cast<const unsigned char*>(data.data()), 
                         data.size()) != 1) {
        EVP_CIPHER_CTX_free(ctx);
        return "";
    }
    
    int finalLen = 0;
    if (EVP_EncryptFinal_ex(ctx, 
                           reinterpret_cast<unsigned char*>(&output[outLen]), 
                           &finalLen) != 1) {
        EVP_CIPHER_CTX_free(ctx);
        return "";
    }
    
    output.resize(outLen + finalLen);
    EVP_CIPHER_CTX_free(ctx);
    return output;
}

std::string decryptData(const std::string& data, int method, 
                        const unsigned char* key, const unsigned char* iv) {
    EVP_CIPHER_CTX* ctx = createCipherContext(method, false, key, iv);
    if (!ctx) return "";
    
    std::string output;
    output.resize(data.size());
    
    int outLen = 0;
    if (EVP_DecryptUpdate(ctx, 
                         reinterpret_cast<unsigned char*>(&output[0]), 
                         &outLen,
                         reinterpret_cast<const unsigned char*>(data.data()), 
                         data.size()) != 1) {
        EVP_CIPHER_CTX_free(ctx);
        return "";
    }
    
    int finalLen = 0;
    if (EVP_DecryptFinal_ex(ctx, 
                           reinterpret_cast<unsigned char*>(&output[outLen]), 
                           &finalLen) != 1) {
        EVP_CIPHER_CTX_free(ctx);
        return "";
    }
    
    output.resize(outLen + finalLen);
    EVP_CIPHER_CTX_free(ctx);
    return output;
}

void generateKeyAndIV(int method, unsigned char* key, unsigned char* iv) {
    int keyLen = 0;
    int ivLen = 0;
    
    switch (method) {
        case AES_128_CBC: keyLen = 16; ivLen = 16; break;
        case AES_256_CBC: keyLen = 32; ivLen = 16; break;
        case AES_128_GCM: keyLen = 16; ivLen = 12; break;
        case AES_256_GCM: keyLen = 32; ivLen = 12; break;
        case CHACHA20: keyLen = 32; ivLen = 16; break;
    }
    
    if (!RAND_bytes(key, keyLen) || !RAND_bytes(iv, ivLen)) {
        std::cerr << "Failed to generate random key/IV" << std::endl;
        exit(1);
    }
}

void benchmarkOpenSSL(std::ofstream& csvFile, const std::string& filename) {
    std::vector<std::string> lines = readFileByLine(filename);
    std::string data;
    for (const auto& line : lines) {
        data += line + "\n";
    }
    
    std::cout << "\n=== OpenSSL Encryption Benchmark ===\n";
    std::cout << "Original data size: " << data.size() << " bytes\n";
    std::cout << "Number of lines: " << lines.size() << "\n";
    
    OpenSSL_add_all_algorithms();
    ERR_load_crypto_strings();
    
    for (int i = 0; i < OPENSSL_METHOD_COUNT; ++i) {
        unsigned char key[32];
        unsigned char iv[16];
        generateKeyAndIV(i, key, iv);
        
        // Encryption
        auto start = std::chrono::high_resolution_clock::now();
        std::string encrypted = encryptData(data, i, key, iv);
        auto end = std::chrono::high_resolution_clock::now();
        auto encryptTime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
        
        if (!encrypted.empty()) {
            // Decryption
            start = std::chrono::high_resolution_clock::now();
            std::string decrypted = decryptData(encrypted, i, key, iv);
            end = std::chrono::high_resolution_clock::now();
            auto decryptTime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
            
            if (!decrypted.empty() && decrypted == data) {
                appendToCSV(csvFile, "OpenSSL", opensslMethodNames[i], encryptTime, 
                          1, data.size(), encrypted.size(), "Encryption");
                appendToCSV(csvFile, "OpenSSL", opensslMethodNames[i], decryptTime, 
                          1, encrypted.size(), data.size(), "Decryption");
                
                std::cout << opensslMethodNames[i] << ": "
                          << "Encrypt " << encryptTime / 1000000.0 << " ms, "
                          << "Decrypt " << decryptTime / 1000000.0 << " ms, "
                          << "Size " << encrypted.size() << " bytes\n";
            }
        }
    }
    
    EVP_cleanup();
    ERR_free_strings();
}

// Floating point benchmark functions
std::vector<double> generateRandomDoubles(size_t count, double min = 0.1, double max = 100.0) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(min, max);

    std::vector<double> result(count);
    for (auto& val : result) {
        val = dis(gen);
    }
    return result;
}

void benchmarkFloatingPoint(std::ofstream& csvFile) {
    std::cout << "\n=== Floating Point Operations Benchmark ===\n";
    std::cout << "Array size: " << FLOAT_ARRAY_SIZE << "\n";
    std::cout << "Iterations: " << ITERATIONS << "\n";
    
    auto data1 = generateRandomDoubles(FLOAT_ARRAY_SIZE);
    auto data2 = generateRandomDoubles(FLOAT_ARRAY_SIZE);
    auto data3 = generateRandomDoubles(FLOAT_ARRAY_SIZE);
    
    volatile double sink;
    
    for (int op = 0; op < OP_COUNT; ++op) {
        uint64_t totalTime = 0;
        uint64_t minTime = UINT64_MAX;
        uint64_t maxTime = 0;
        
        for (int iter = 0; iter < ITERATIONS; ++iter) {
            auto start = std::chrono::high_resolution_clock::now();
            
            switch (op) {
                case OP_ADD: {
                    for (size_t i = 0; i < FLOAT_ARRAY_SIZE; ++i) {
                        sink = data1[i] + data2[i];
                    }
                    break;
                }
                case OP_SUB: {
                    for (size_t i = 0; i < FLOAT_ARRAY_SIZE; ++i) {
                        sink = data1[i] - data2[i];
                    }
                    break;
                }
                case OP_MUL: {
                    for (size_t i = 0; i < FLOAT_ARRAY_SIZE; ++i) {
                        sink = data1[i] * data2[i];
                    }
                    break;
                }
                case OP_DIV: {
                    for (size_t i = 0; i < FLOAT_ARRAY_SIZE; ++i) {
                        sink = data1[i] / data2[i];
                    }
                    break;
                }
                case OP_FMA: {
                    for (size_t i = 0; i < FLOAT_ARRAY_SIZE; ++i) {
                        sink = std::fma(data1[i], data2[i], data3[i]);
                    }
                    break;
                }
                case OP_SQRT: {
                    for (size_t i = 0; i < FLOAT_ARRAY_SIZE; ++i) {
                        sink = std::sqrt(data1[i]);
                    }
                    break;
                }
                case OP_TRIG: {
                    for (size_t i = 0; i < FLOAT_ARRAY_SIZE; ++i) {
                        sink = std::sin(data1[i]) + std::cos(data2[i]) + std::tan(data3[i]);
                    }
                    break;
                }
                case OP_LOG: {
                    for (size_t i = 0; i < FLOAT_ARRAY_SIZE; ++i) {
                        sink = std::log(data1[i]) + std::log10(data2[i]);
                    }
                    break;
                }
                case OP_MIXED: {
                    for (size_t i = 0; i < FLOAT_ARRAY_SIZE; ++i) {
                        sink = (data1[i] + data2[i]) * data3[i] - 
                               std::sqrt(data1[i]) / std::log(data2[i]) + 
                               std::sin(data3[i]);
                    }
                    break;
                }
            }
            
            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
            
            totalTime += duration;
            minTime = std::min<float>(minTime, duration);
            maxTime = std::max<float>(maxTime, duration);
        }
        
        double avgTime = totalTime / (double)ITERATIONS;
        appendToCSV(csvFile, "Floating Point", floatOpNames[op], avgTime, 
                   FLOAT_ARRAY_SIZE * ITERATIONS);
        
        std::cout << floatOpNames[op] << ": "
                  << "Avg " << avgTime / 1000000.0 << " ms, "
                  << "Min " << minTime / 1000000.0 << " ms, "
                  << "Max " << maxTime / 1000000.0 << " ms, "
                  << FLOAT_ARRAY_SIZE / (avgTime / 1000000000.0) / 1000000.0 << " M ops/s\n";
    }
}

int main(int argc, char* argv[]) {
    std::cout << "=== System Information ===\n";
    std::string cpuModel = getCPUModel();
    std::cout << "CPU: " << cpuModel << "\n";
    std::cout << "Timestamp: " << getCurrentTimestamp() << "\n\n";
    
    // Open CSV file for writing
    std::ofstream csvFile(CSV_FILENAME);
    if (!csvFile.is_open()) {
        std::cerr << "Failed to open CSV file for writing" << std::endl;
        return 1;
    }
    writeCSVHeader(csvFile);
    
    // Run benchmarks
    if (argc > 1) {
        benchmarkZLIB(csvFile, argv[1]);
        benchmarkOpenSSL(csvFile, argv[1]);
    } else {
        std::cout << "No input file provided, skipping ZLIB and OpenSSL benchmarks\n";
    }
    
    benchmarkFloatingPoint(csvFile);
    
    csvFile.close();
    std::cout << "\nResults saved to " << CSV_FILENAME << std::endl;
    
    return 0;
}