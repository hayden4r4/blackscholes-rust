# BlackScholes Performance Benchmarking

This document outlines the modern benchmarking infrastructure for the BlackScholes library, which leverages GitHub Actions and GitHub Pages for automated performance tracking and visualization.

## Overview

The benchmarking system is built around **Criterion.rs** and **github-action-benchmark**, providing:

- 📊 **Interactive Charts**: Professional visualizations hosted on GitHub Pages
- 🔍 **Automated Regression Detection**: CI integration that catches performance regressions
- 📈 **Historical Tracking**: Continuous monitoring of performance trends over time  
- 💬 **PR Integration**: Automatic benchmark comparison comments on pull requests

**Live Benchmark Results**: https://przemyslawolszewski.github.io/bs-rs/

## Automated Benchmarking Workflow

### Trigger Conditions

Benchmarks run automatically when:

- **Pull Requests**: Against `main`/`master` with changes to:
  - `src/**` (source code changes)
  - `benches/**` (benchmark changes) 
  - `Cargo.*` (dependency changes)
  - `.github/workflows/benchmark.yml` (workflow changes)

- **Pushes**: To `main`/`master` branch (updates the baseline and GitHub Pages)

- **Nightly**: Daily at 02:00 UTC to maintain continuous baseline data

### Workflow Features

- **Path Filtering**: Only runs on relevant file changes (performance optimization)
- **Concurrency Control**: Prevents overlapping benchmark runs
- **Pinned Environment**: Uses `ubuntu-22.04` for reproducible results
- **Artifact Preservation**: Criterion HTML reports available for 7 days
- **Security**: Minimal permissions with `contents:write` and `pull-requests:write`

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

## Running Benchmarks Locally

### Basic Usage

Run a specific benchmark:
```bash
cargo bench --bench single_option
```

Run all benchmarks:
```bash
cargo bench
```

### Generate CSV Output (GitHub Action Compatible)

To generate CSV output compatible with the GitHub Action:
```bash
cargo bench -- --output-format=csv
```

### View Local Results

Criterion generates HTML reports in `target/criterion/`. Open `target/criterion/report/index.html` in your browser to view detailed results.

## Performance Regression Detection

### Automated Detection

The GitHub Action automatically:
- ✅ **Comments on PRs**: Shows performance comparison between PR and base branch
- ⚠️ **Alert Threshold**: Warns when performance degrades by >10%
- ❌ **Fail on Alert**: Fails the CI check if regression threshold is exceeded
- 📊 **Summary Reports**: Provides structured benchmark comparison tables

### Configuration

Current thresholds:
- **Alert Threshold**: 110% (10% slower triggers warning)
- **Fail on Alert**: `true` (CI fails on regressions)
- **Comment on Alert**: `true` (PR comments enabled)

### Threshold Tuning

To adjust sensitivity, modify `.github/workflows/benchmark.yml`:
```yaml
alert-threshold: '105%'  # 5% threshold (more sensitive)
# or
alert-threshold: '120%'  # 20% threshold (less sensitive)
```

## Adding New Benchmarks

To add a new benchmark:

1. **Create the benchmark file** in the appropriate directory:
   ```rust
   // benches/new_feature/my_benchmark.rs
   use criterion::{criterion_group, criterion_main, Criterion};

   fn benchmark_my_feature(c: &mut Criterion) {
       c.bench_function("my_feature", |b| {
           b.iter(|| {
               // Your benchmark code here
           })
       });
   }

   criterion_group!(benches, benchmark_my_feature);
   criterion_main!(benches);
   ```

2. **Add to Cargo.toml**:
   ```toml
   [[bench]]
   name = "my_benchmark"
   path = "benches/new_feature/my_benchmark.rs"
   harness = false
   ```

3. **Commit and push** - the workflow will automatically include your new benchmark!

## System Dependencies

For local development, ensure you have:

- **Gnuplot**: Required for Criterion's HTML reports
  - Ubuntu/Debian: `sudo apt-get install gnuplot`
  - Fedora/RHEL: `sudo dnf install gnuplot`  
  - Arch Linux: `sudo pacman -S gnuplot`
  - macOS: `brew install gnuplot`

## Repository Configuration

### GitHub Pages Setup

**Repository Settings → Pages:**
- Source: "Deploy from branch"
- Branch: `gh-pages`
- Folder: `/` (root)

### Required Permissions

**Repository Settings → Actions → General:**
- Workflow permissions: "Read and write permissions"
- "Allow GitHub Actions to create and approve pull requests": ✅

## Monitoring and Maintenance

### Historical Data

- **Baseline Collection**: Nightly runs ensure continuous data even during quiet periods
- **Trend Analysis**: GitHub Pages shows performance trends over time
- **Version Tracking**: Each commit's performance is tracked and visualized

### Best Practices

1. **Monitor False Positives**: Adjust `alert-threshold` if too many false alerts occur
2. **Review PR Comments**: Use benchmark comparison comments to validate performance changes
3. **Check GitHub Pages**: Regularly review the trend charts for gradual performance changes
4. **Local Reproduction**: Use the CSV format locally to reproduce CI results

### Future Enhancements

Consider these additions as the project evolves:

- **Multi-Platform Testing**: Matrix builds for different Rust versions/platforms
- **Benchmark Categorization**: Separate different benchmark suites
- **External Integration**: Connect with performance monitoring tools
- **Custom Metrics**: Add domain-specific performance measurements

## Troubleshooting

### Common Issues

1. **No benchmark results**: Ensure Criterion outputs CSV with `--output-format=csv`
2. **Permission errors**: Verify repository settings have write permissions enabled
3. **Missing GitHub Pages**: Check that `gh-pages` branch exists and Pages is configured
4. **Threshold too sensitive**: Adjust `alert-threshold` in the workflow file

### Getting Help

- Check the [GitHub Action logs](../../actions) for detailed error messages
- Review the [github-action-benchmark documentation](https://github.com/benchmark-action/github-action-benchmark)
- Examine the `target/criterion/` directory structure locally 