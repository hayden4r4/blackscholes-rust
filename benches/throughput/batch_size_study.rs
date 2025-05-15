use blackscholes::{Greeks, Inputs, Pricing};
use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use rand::thread_rng;
use std::time::Duration;

#[path = "../common.rs"]
mod common;
use common::{generate_random_inputs, BatchSize};

// Simple batch size study to analyze scaling behavior
fn bench_batch_size_scaling(c: &mut Criterion) {
    let mut group = c.benchmark_group("Batch Size Scaling");

    // Configure timeout for large batch sizes
    group.measurement_time(Duration::from_secs(5));
    group.warm_up_time(Duration::from_millis(500));

    // Test a range of batch sizes to detect scaling behaviors
    let batch_sizes = [
        10,     // Tiny
        100,    // Small
        1_000,  // Medium
        5_000,  // Medium-large
        10_000, // Large (if time permits)
    ];

    let mut rng = thread_rng();

    for size in batch_sizes {
        // Generate random inputs
        let inputs = generate_random_inputs(size, &mut rng);

        // Set throughput measurement for proper ops/sec calculation
        group.throughput(Throughput::Elements(size as u64));

        // Standard price calculation
        group.bench_with_input(
            BenchmarkId::new("standard_price", size),
            &inputs,
            |b, inputs| {
                b.iter(|| {
                    let mut results = Vec::with_capacity(inputs.len());
                    for input in black_box(inputs) {
                        results.push(input.calc_price().unwrap());
                    }
                    black_box(results)
                })
            },
        );

        // Rational price calculation
        group.bench_with_input(
            BenchmarkId::new("rational_price", size),
            &inputs,
            |b, inputs| {
                b.iter(|| {
                    let mut results = Vec::with_capacity(inputs.len());
                    for input in black_box(inputs) {
                        results.push(input.calc_rational_price().unwrap());
                    }
                    black_box(results)
                })
            },
        );

        // Delta calculation
        group.bench_with_input(BenchmarkId::new("delta", size), &inputs, |b, inputs| {
            b.iter(|| {
                let mut results = Vec::with_capacity(inputs.len());
                for input in black_box(inputs) {
                    results.push(input.calc_delta().unwrap());
                }
                black_box(results)
            })
        });

        // Gamma calculation
        group.bench_with_input(BenchmarkId::new("gamma", size), &inputs, |b, inputs| {
            b.iter(|| {
                let mut results = Vec::with_capacity(inputs.len());
                for input in black_box(inputs) {
                    results.push(input.calc_gamma().unwrap());
                }
                black_box(results)
            })
        });

        // All greeks calculation (to test if calculating all at once is more efficient)
        group.bench_with_input(
            BenchmarkId::new("all_greeks", size),
            &inputs,
            |b, inputs| {
                b.iter(|| {
                    let mut results = Vec::with_capacity(inputs.len());
                    for input in black_box(inputs) {
                        results.push(input.calc_all_greeks().unwrap());
                    }
                    black_box(results)
                })
            },
        );
    }

    group.finish();
}

criterion_group!(
    name = benches;
    config = Criterion::default().with_plots();
    targets = bench_batch_size_scaling
);
criterion_main!(benches);
