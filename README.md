# System Benchmark Tool

A comprehensive benchmarking tool that tests:
- ZLIB compression/decompression performance
- OpenSSL encryption/decryption performance
- Floating-point operation performance

## Features

- CPU model detection
- Timestamped results
- CSV output for easy analysis
- Multiple algorithms/methods tested
- Detailed performance metrics

## Requirements

- Linux system
- C++17 compatible compiler
- CMake 3.10+
- ZLIB development libraries
- OpenSSL development libraries

## Installation

1. Clone the repository:
   ```bash
git clone https://github.com/809898292/system-benchmark.git
cd system-benchmark 

python3 generate.py
mkdir build
mv -f stock_data.csv build/
cd build
cmake .. && make

echo "benchmark_results.csv is the generated test result"
echo "Navigate to the /system-benchmark/build directory"
echo "Run ./SystemBenchmark stock_data.csv to test the CPU's performance for high-frequency trading"
