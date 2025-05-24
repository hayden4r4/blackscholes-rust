use std::{hint::black_box, time::Duration};

use blackscholes::Pricing;
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use indicatif::{ProgressBar, ProgressStyle};
use rand::thread_rng;

#[path = "../common.rs"]
mod common;
use common::{generate_random_inputs, get_sample_config};

const BATCH_SIZES: [usize; 9] = [
    10,      // Tiny
    100,     // Small
    500,     // Medium-small
    1_000,   // Medium
    5_000,   // Medium-large
    10_000,  // Large
    25_000,  // Larger
    50_000,  // Very large
    100_000, // Huge
];

fn bench_scaling(c: &mut Criterion) {
    let mut group = c.benchmark_group("Batch Size Scaling");

    group.warm_up_time(Duration::from_millis(500));

    let mut rng = thread_rng();

    let progress_bar = ProgressBar::new(BATCH_SIZES.len() as u64);
    progress_bar.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} {msg}")
            .unwrap()
            .progress_chars("#>-"),
    );
    progress_bar.set_message("Running batch scaling benchmarks...");

    // Run benchmarks for each batch size
    for (idx, &size) in BATCH_SIZES.iter().enumerate() {
        progress_bar.set_position(idx as u64);
        progress_bar.set_message(format!("Benchmarking batch size: {}", size));

        let (sample_count, measurement_time) = get_sample_config(size);
        group.sample_size(sample_count);
        group.measurement_time(measurement_time);

        let inputs = generate_random_inputs(size, &mut rng);

        group.throughput(Throughput::Elements(size as u64));

        let id = BenchmarkId::new("price", size);
        group.bench_function(id, |b| {
            b.iter(|| {
                let mut results = Vec::with_capacity(inputs.len());
                for input in black_box(&inputs) {
                    results.push(input.calc_price().unwrap());
                }
                black_box(results)
            })
        });

        println!(
            "Completed benchmarking batch size: {} (approximate throughput varies)",
            size
        );
    }

    group.finish();
    progress_bar.finish_with_message("Benchmarks complete");

    println!("Benchmark completed. Use criterion's HTML reports for detailed analysis.");
}

criterion_group!(
    name = benches;
    config = Criterion::default().with_plots();
    targets = bench_scaling
);
criterion_main!(benches);
