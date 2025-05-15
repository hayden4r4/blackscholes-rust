# BlackScholes Performance Benchmarking

This document outlines the benchmarking infrastructure for the BlackScholes library, explaining how to run benchmarks, analyze results, and track performance improvements.

## Overview

The benchmarking system is built around Criterion.rs, a statistics-driven benchmarking framework for Rust. We have extended it with custom scripts for visualization, regression detection, and reporting.

## System Dependencies

Before running benchmarks, ensure you have the following system dependencies installed:

- **Fontconfig**: Required for plotting and visualization
  - Ubuntu/Debian: `sudo apt-get install libfontconfig1-dev`
  - Fedora/RHEL: `sudo dnf install fontconfig-devel`
  - Arch Linux: `sudo pacman -S fontconfig`
  - macOS: `brew install fontconfig`

- **Gnuplot**: Required for generating benchmark plots
  - Ubuntu/Debian: `sudo apt-get install gnuplot`
  - Fedora/RHEL: `sudo dnf install gnuplot`
  - Arch Linux: `sudo pacman -S gnuplot`
  - macOS: `brew install gnuplot`

The benchmark scripts will check for these dependencies and provide instructions if they're missing.

## Directory Structure

```
benches/
├── single/            # Benchmarks for single option calculations
│   ├── option_pricing.rs
│   ├── greeks.rs
│   └── implied_volatility.rs
├── batch/             # Benchmarks for batch operations
│   ├── pricing.rs
│   └── greeks.rs
├── throughput/        # Benchmarks focused on high-volume operations
│   ├── option_pricing.rs
│   ├── scaling.rs
│   └── batch_size_study.rs
└── common.rs          # Shared benchmark utilities
```

## Running Benchmarks

### Basic Usage

Run a specific benchmark:
```bash
cargo bench --bench single_option
```

Run all benchmarks (including throughput benchmarks):
```bash
cargo bench --features bench_all
```

### Using the Benchmark Report Script

The `scripts/benchmark_report.sh` script automates running all benchmarks and generates an HTML report:

```bash
./scripts/benchmark_report.sh
```

This will:
1. Run all benchmarks
2. Generate HTML reports in `target/criterion/`
3. Start a local HTTP server to view the reports

### Performance Profiling with Flamegraphs

Generate flamegraphs to identify performance hotspots:

```bash
# Install cargo-flamegraph if not already installed
cargo install flamegraph

# Generate flamegraph for a specific benchmark
./scripts/generate_flamegraph.sh single_option

# Generate flamegraphs for all benchmarks
./scripts/generate_flamegraph.sh --all
```

Flamegraphs will be saved to `target/flamegraphs/`.

## Detecting Performance Regressions

### Local Regression Detection

The `scripts/compare_benchmarks.sh` script helps detect performance regressions:

```bash
# First, save the current state as baseline
./scripts/compare_benchmarks.sh --save

# Make your changes, then compare against the baseline
./scripts/compare_benchmarks.sh

# Run specific benchmarks
./scripts/compare_benchmarks.sh --bench single_option --bench single_greeks

# Include batch benchmarks with required features
./scripts/compare_benchmarks.sh --bench batch_pricing --features bench_throughput
```

The script will show a table with performance changes, highlighting improvements and regressions.

### CI Integration

Performance regression detection is integrated into the CI pipeline. It will:

1. Run on pull requests to the main branch
2. Compare benchmark results against the base branch
3. Flag warnings for performance degradations of 5% or more
4. Fail the build for degradations of 10% or more
5. Generate a summary report in the GitHub PR

## Benchmark Categories

### 1. Single Option Operations

Benchmarks for individual option calculations:
- Black-Scholes pricing (`calc_price`)
- Rational pricing (`calc_rational_price`)
- Greeks calculation (Delta, Gamma, Vega, Theta, Rho)
- Implied volatility calculation

### 2. Batch Operations

Benchmarks for processing batches of options:
- Sequential batch processing
- Struct-of-arrays (SoA) format for potential SIMD optimization
- Varying batch sizes (10, 100, 1,000, etc.)

### 3. High-Volume Throughput

Benchmarks focused on operations per second:
- Large batch processing (10,000+ options)
- Throughput measurements with different processing strategies
- Scaling behavior as batch size increases

## Adding New Benchmarks

To add a new benchmark:

1. Create a new file in the appropriate directory
2. Add the benchmark to `Cargo.toml`:
   ```toml
   [[bench]]
   name = "my_benchmark"
   path = "benches/category/my_benchmark.rs"
   harness = false
   ```
3. If it requires features, add:
   ```toml
   required-features = ["bench_throughput"]
   ```

## Best Practices

1. **Consistent Environments**: Run comparative benchmarks on the same hardware
2. **Multiple Samples**: Don't rely on a single run; Criterion automatically runs multiple samples
3. **Focus on Critical Paths**: Prioritize benchmarking performance-critical operations
4. **Baseline Comparison**: Always compare against a baseline when evaluating optimizations
5. **Profile First**: Use flamegraphs to identify hotspots before trying to optimize
6. **Document Improvements**: Add notes about performance improvements to your PRs 