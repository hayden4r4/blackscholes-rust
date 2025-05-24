use std::time::Duration;

use blackscholes::{Greeks, Inputs, OptionType, Pricing};
use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use rand::thread_rng;

#[path = "../common.rs"]
mod common;
use common::{generate_random_inputs, get_sample_config, BatchSize, InputsSoA};

// High-throughput benchmarking for option pricing
fn bench_throughput(c: &mut Criterion) {
    let mut group = c.benchmark_group("Option Pricing Throughput");

    // Configure the benchmark group
    group.warm_up_time(Duration::from_millis(500));

    // For throughput testing, we'll use larger batch sizes
    let batch_sizes = [
        (BatchSize::Medium as usize, "medium"),
        (BatchSize::Large as usize, "large"),
    ];

    let mut rng = thread_rng();

    for &(size, size_name) in batch_sizes.iter() {
        // Adjust sample count and measurement time based on batch size
        let (sample_count, measurement_time) = get_sample_config(size);
        group.sample_size(sample_count);
        group.measurement_time(measurement_time);

        // Generate random inputs for batch processing
        let inputs = generate_random_inputs(size, &mut rng);

        // Set throughput measurement - this will show operations per second
        group.throughput(Throughput::Elements(size as u64));

        // Benchmark price calculation throughput
        group.bench_function(BenchmarkId::new("price", size_name), |b| {
            b.iter(|| {
                let mut results = Vec::with_capacity(inputs.len());
                for input in black_box(&inputs) {
                    results.push(input.calc_price().unwrap());
                }
                black_box(results)
            })
        });

        // Benchmark rational price calculation throughput
        group.bench_function(BenchmarkId::new("rational_price", size_name), |b| {
            b.iter(|| {
                let mut results = Vec::with_capacity(inputs.len());
                for input in black_box(&inputs) {
                    results.push(input.calc_rational_price().unwrap());
                }
                black_box(results)
            })
        });

        // Benchmark delta calculation throughput
        group.bench_function(BenchmarkId::new("delta", size_name), |b| {
            b.iter(|| {
                let mut results = Vec::with_capacity(inputs.len());
                for input in black_box(&inputs) {
                    results.push(input.calc_delta().unwrap());
                }
                black_box(results)
            })
        });

        // This is a placeholder for future parallel implementation with Rayon
        // We'll implement a parallel version later with Rayon in Task 9
        // group.bench_function(BenchmarkId::new("parallel_price", size_name), |b| {
        //     b.iter(|| {
        //         // Future parallel implementation
        //     })
        // });
    }

    group.finish();
}

criterion_group!(
    name = benches;
    config = Criterion::default()
        .with_plots()
        .measurement_time(Duration::from_secs(10));
    targets = bench_throughput
);
criterion_main!(benches);
