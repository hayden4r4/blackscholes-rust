use std::{hint::black_box, time::Duration};

use blackscholes::{Greeks, Inputs, OptionType, Pricing};
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use rand::thread_rng;

#[path = "../common/mod.rs"]
mod common;
use common::{generate_random_inputs, get_sample_config, BatchSize, InputsSoA};

fn bench_throughput(c: &mut Criterion) {
    let mut group = c.benchmark_group("Option Pricing Throughput");

    group.warm_up_time(Duration::from_millis(500));

    let batch_sizes = [
        (BatchSize::Medium as usize, "medium"),
        (BatchSize::Large as usize, "large"),
    ];

    let mut rng = thread_rng();

    for &(size, size_name) in batch_sizes.iter() {
        let (sample_count, measurement_time) = get_sample_config(size);
        group.sample_size(sample_count);
        group.measurement_time(measurement_time);

        let inputs = generate_random_inputs(size, &mut rng);

        group.throughput(Throughput::Elements(size as u64));

        group.bench_function(BenchmarkId::new("price", size_name), |b| {
            b.iter(|| {
                let mut results = Vec::with_capacity(inputs.len());
                for input in black_box(&inputs) {
                    results.push(input.calc_price().unwrap());
                }
                black_box(results)
            })
        });

        group.bench_function(BenchmarkId::new("rational_price", size_name), |b| {
            b.iter(|| {
                let mut results = Vec::with_capacity(inputs.len());
                for input in black_box(&inputs) {
                    results.push(input.calc_rational_price().unwrap());
                }
                black_box(results)
            })
        });

        group.bench_function(BenchmarkId::new("delta", size_name), |b| {
            b.iter(|| {
                let mut results = Vec::with_capacity(inputs.len());
                for input in black_box(&inputs) {
                    results.push(input.calc_delta().unwrap());
                }
                black_box(results)
            })
        });
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
