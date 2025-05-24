use std::time::Duration;

use blackscholes::{Greeks, Inputs, OptionType, Pricing};
use criterion::{
    black_box, criterion_group, criterion_main, measurement::WallTime, BenchmarkId, Criterion,
};
use rand::thread_rng;

#[path = "../common/mod.rs"]
mod common;
use common::{generate_random_inputs, get_sample_config, BatchSize, InputsSoA};

// Batch Greeks calculation benchmark
fn bench_batch_greeks(c: &mut Criterion) {
    let mut group = c.benchmark_group("Batch Greeks Calculation");

    // Configure the benchmark group
    group.warm_up_time(Duration::from_millis(500));

    let batch_sizes = [
        (BatchSize::Tiny as usize, "tiny"),
        (BatchSize::Small as usize, "small"),
    ];

    let mut rng = thread_rng();

    for &(size, size_name) in batch_sizes.iter() {
        // Adjust sample count and measurement time based on batch size
        let (sample_count, measurement_time) = get_sample_config(size);
        group.sample_size(sample_count);
        group.measurement_time(measurement_time);

        // Generate random inputs for batch processing
        let inputs = generate_random_inputs(size, &mut rng);

        // Benchmark naive batch processing of individual Greeks
        group.bench_function(BenchmarkId::new("delta", size_name), |b| {
            b.iter(|| {
                let mut results = Vec::with_capacity(inputs.len());
                for input in black_box(&inputs) {
                    results.push(input.calc_delta().unwrap());
                }
                black_box(results)
            })
        });

        group.bench_function(BenchmarkId::new("gamma", size_name), |b| {
            b.iter(|| {
                let mut results = Vec::with_capacity(inputs.len());
                for input in black_box(&inputs) {
                    results.push(input.calc_gamma().unwrap());
                }
                black_box(results)
            })
        });

        group.bench_function(BenchmarkId::new("vega", size_name), |b| {
            b.iter(|| {
                let mut results = Vec::with_capacity(inputs.len());
                for input in black_box(&inputs) {
                    results.push(input.calc_vega().unwrap());
                }
                black_box(results)
            })
        });

        // Benchmark batch all_greeks calculation
        group.bench_function(BenchmarkId::new("all_greeks", size_name), |b| {
            b.iter(|| {
                let mut results = Vec::with_capacity(inputs.len());
                for input in black_box(&inputs) {
                    results.push(input.calc_all_greeks().unwrap());
                }
                black_box(results)
            })
        });

        // Convert to SoA format for future optimization
        let inputs_soa = InputsSoA::from_inputs(&inputs);

        // Placeholder for future SIMD-optimized batch Greeks calculation
        group.bench_function(BenchmarkId::new("delta_soa", size_name), |b| {
            b.iter(|| {
                let mut results = Vec::with_capacity(inputs_soa.len());
                for i in 0..inputs_soa.len() {
                    let input = Inputs::new(
                        inputs_soa.option_types[i],
                        inputs_soa.spots[i],
                        inputs_soa.strikes[i],
                        None,
                        inputs_soa.rates[i],
                        inputs_soa.dividends[i],
                        inputs_soa.times[i],
                        Some(inputs_soa.volatilities[i]),
                    );
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
    targets = bench_batch_greeks
);
criterion_main!(benches);
