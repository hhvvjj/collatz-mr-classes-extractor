# mr Classes Extractor

A high-performance parallel C program for discovering and classifying n values that generate specific `mr` values in Collatz sequences using the tuple-based transform with detailed taxonomy analysis.

## Features

This program analyzes Collatz sequences to find values of `n` that produce a specific repeated `m` parameter (`mr`) and classifies them into taxonomic categories based on the position of maximum values relative to repetition points. For each positive integer `n`, it computes the Collatz sequence, applies the transformation `m = (c - p) / 2` where `p = 2` if `c` is even, `p = 1` if `c` is odd, and analyzes the resulting sequence structure.

- **High Performance**: Optimized parallel implementation using OpenMP
- **Taxonomy Classification**: Automatic classification into types A, B, C
- **Detailed Analysis**: Tracks positions, maximum values, and sequence statistics
- **Adaptive Scheduling**: Automatically chooses optimal work distribution strategy
- **Memory Efficient**: Dynamic memory allocation with intelligent growth
- **Large Scale**: Supports analysis up to 2^64 numbers
- **Thread Safe**: Lock-optimized concurrent access to shared data structures
- **Structured Output**: CSV format grouped by taxonomy for easy analysis

## Taxonomy Classifications

- **Type A**: Maximum m occurs BEFORE the first mr repetition
- **Type B**: Maximum m occurs BETWEEN first and second mr repetition  
- **Type C**: Maximum m occurs AT or AFTER the second mr repetition

## Dependencies

### Ubuntu/Debian
```bash
sudo apt-get install gcc libomp-dev
```

### CentOS/RHEL
```bash
sudo yum install gcc libgomp-devel
```

## Installation

```bash
# Clone the repository
git clone https://github.com/hhvvjj/collatz-mr-classes-extractor.git
cd collatz-mr-classes-extractor

# Compile 
gcc -fopenmp -O3 -march=native mr_classes_extractor.c -lgomp -lpthread -o mr_classes_extractor
```

## Usage

```bash
./mr_classes_extractor <exponent> <target_mr>
```

### Examples

```bash
# Find all n < 2^20 that generate mr = 444
./mr_classes_extractor 20 444

# Find all n < 2^25 that generate mr = 0 (trivial cycle)
./mr_classes_extractor 25 0

# Large scale analysis - find n < 2^30 that generate mr = 166
./mr_classes_extractor 30 141
```

## Performance Guide

| Exponent | Range | Numbers | Est. Time* | Memory Usage |      Use Cases     |
|----------|-------|---------|------------|--------------|--------------------|
|    15    |  2^15 |   32K   |    < 1 sec |    < 1 MB    | Testing            |
|    20    |  2^20 |    1M   |     ~5 sec |    ~10 MB    | Quick analysis     |
|    25    |  2^25 |   33M   |   ~2-5 min |    ~50 MB    | Standard research  |
|    30    |  2^30 |    1B   | ~1-3 hours |   ~200 MB    | Intensive analysis |
|    35    |  2^35 |   34B   |  ~1-2 days |     ~1 GB    | Extreme analysis   |

*Times are approximate and depend on hardware (8-core CPU assumed)

## Output

### Real-time Progress
```
Using 2 threads for parallelization
Searching for all n < 2^20 that generate mr = 141
Range: from 1 to 1048575 

Using scheduling strategy: static (small range)
Number of threads: 2

 [*] Found n = 189: mr = 141, taxonomy = C, first_mr_pos = 3, second_mr_pos = 47, max_m = 4615 at pos 73 (total: 1)
 [*] Found n = 284: mr = 141, taxonomy = C, first_mr_pos = 1, second_mr_pos = 45, max_m = 4615 at pos 71 (total: 2)
 [*] Found n = 568: mr = 141, taxonomy = C, first_mr_pos = 2, second_mr_pos = 46, max_m = 4615 at pos 72 (total: 3)
 [*] Found n = 757: mr = 141, taxonomy = C, first_mr_pos = 5, second_mr_pos = 49, max_m = 4615 at pos 75 (total: 4)
 [*] Found n = 1009: mr = 141, taxonomy = C, first_mr_pos = 8, second_mr_pos = 52, max_m = 4615 at pos 78 (total: 5)
 [*] Found n = 1136: mr = 141, taxonomy = C, first_mr_pos = 3, second_mr_pos = 47, max_m = 4615 at pos 73 (total: 6)
```

### Final Results
```
SEARCH RESULTS:

All n values that generate mr = 141:

Showing first 50 results (total found: 2132):

        n      Tax    1st_mr     2nd_mr      Max_m   Max_m_pos    Max_collatz  Max_collatz_pos
 -------- -------- --------- ---------- ---------- ----------- -------------- ----------------
      189        C         3         47       4615          73           9232               73
      284        C         1         45       4615          71           9232               71
      568        C         2         46       4615          72           9232               72
      757        C         5         49       4615          75           9232               75
     1009        C         8         52       4615          78           9232               78
     1136        C         3         47       4615          73           9232               73
     1195        C        16         60       4615          86           9232               86
     1345        C        11         55       4615          81           9232               81
     1514        C         6         50       4615          76           9232               76
 ... and 8852 more values

Taxonomy Distribution:
 Type A: 8821 (99.6%)
 Type B: 0 (0.0%)
 Type C: 32 (0.4%)

Statistics:
 Average first mr position: 64.7
 Average second mr position: 108.7
 Average distance between mr: 44.0
 Smallest maximum m value: 4615
 Largest maximum m value: 260895247
 Smallest maximum Collatz value: 9232
 Largest maximum Collatz value: 521790496

First few n values: 189, 284, 568, 757, 1009, 1136, 1195, 1345, 1514, 1593, ...

Final Summary:

 [*] Target mr value: 141
 [*] Search range: 1 to 2^20 = 1048575
 [*] Total time: 5.74 seconds
 [*] Speed: 182730.2 numbers/second
 [*] Numbers processed: 1048575
 [*] Numbers matching criteria: 8853
 [*] Percentage matching criteria: 0.844289%
```

### CSV Output Format

The program generates a CSV file with results grouped by taxonomy:

```csv
mr,taxonomy,n_values
141,A,"1767,2391,2651,2691,2983,3187,..."
141,C,"189,284,568,757,1009,1136,1195,..."
```

## Research Applications

This tool is designed for researchers studying:

- **Collatz Sequence Structure**: Understanding how different starting values of n lead to finite mr values repetition patterns
- **Taxonomy Distribution**: Analyzing the relative frequency of different sequence types
- **Scaling Patterns**: Investigating how mr values and their associated n values scale
- **Convergence Analysis**: Studying the relationship between sequence maximum values and repetition positions

## Algorithm Details

### Core Algorithm
1. **Sequence Generation**: Standard 3n+1 Collatz sequence
2. **Transformation**: Apply m = (c-p)/2 for each step
3. **Repetition Detection**: Track all m values until first repetition occurs
4. **Taxonomy Classification**: Analyze relative positions of maximum m and mr occurrences
5. **Statistical Collection**: Track detailed metrics for comprehensive analysis

### Parallel Strategy
- **Work Distribution**: Adaptive scheduling (static/guided/dynamic)
- **Thread Safety**: Lock-optimized result collection and progress tracking
- **Memory Management**: Dynamic allocation with intelligent growth
- **Real-time Reporting**: Live discovery updates with detailed statistics

### Memory Optimization
- **Hash Tables**: Fast repetition detection using optimized hash functions
- **Dynamic Arrays**: Intelligent capacity growth to minimize reallocations
- **Lock Minimization**: Careful synchronization design for maximum throughput

## Research Background

This implementation extends the foundational work on mr value discovery by adding taxonomic classification based on sequence structure analysis. By examining the relative positions of maximum values and repetition points, we can categorize sequences into distinct behavioral classes that reveal deeper patterns in Collatz dynamics.

The taxonomy system provides insight into how different starting values approach their first repetition, with implications for understanding convergence patterns and sequence complexity.

**Reference:** The complete theoretical framework is detailed in the research paper: http://dx.doi.org/10.5281/zenodo.15546925

## Contributing

Contributions are welcome! Areas of interest:

- Extended validation testing using large exponents
- Additional statistical analysis features
- Performance optimizations for specific architectures
- Documentation improvements and example analyses

## Academic Citation

For academic use, please cite both the original theoretical work and this implementation.

## License

CC-BY-NC-SA 4 License - see LICENSE file for details.

## Contact

For questions about the algorithm implementation or mathematical details drop me an email.