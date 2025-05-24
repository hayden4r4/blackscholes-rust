use std::{hint::black_box, time::Duration};

use blackscholes::{Greeks, Inputs, OptionType};
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};

#[path = "../common.rs"]
mod common;
use common::{generate_standard_inputs, Moneyness, TimeToMaturity, VolatilityLevel};

fn bench_greeks(c: &mut Criterion) {
    let mut group = c.benchmark_group("Single Option Greeks");

    // Configure the benchmark group
    group.warm_up_time(Duration::from_millis(500));
    group.measurement_time(Duration::from_secs(2));

    // We'll test ATM call options with medium volatility at different maturities
    let configurations = [
        ("short_term", TimeToMaturity::ShortTerm),
        ("medium_term", TimeToMaturity::MediumTerm),
        ("long_term", TimeToMaturity::LongTerm),
    ];

    // First-order Greeks (most frequently used)
    for (name, maturity) in configurations.iter() {
        let inputs = generate_standard_inputs(
            OptionType::Call,
            Moneyness::AtTheMoney,
            *maturity,
            VolatilityLevel::Medium,
        );

        // Delta
        group.bench_with_input(
            BenchmarkId::new("delta", name),
            &inputs,
            |b, inputs: &Inputs| b.iter(|| black_box(inputs.calc_delta().unwrap())),
        );

        // Gamma
        group.bench_with_input(
            BenchmarkId::new("gamma", name),
            &inputs,
            |b, inputs: &Inputs| b.iter(|| black_box(inputs.calc_gamma().unwrap())),
        );

        // Theta
        group.bench_with_input(
            BenchmarkId::new("theta", name),
            &inputs,
            |b, inputs: &Inputs| b.iter(|| black_box(inputs.calc_theta().unwrap())),
        );

        // Vega
        group.bench_with_input(
            BenchmarkId::new("vega", name),
            &inputs,
            |b, inputs: &Inputs| b.iter(|| black_box(inputs.calc_vega().unwrap())),
        );

        // Rho
        group.bench_with_input(
            BenchmarkId::new("rho", name),
            &inputs,
            |b, inputs: &Inputs| b.iter(|| black_box(inputs.calc_rho().unwrap())),
        );
    }

    group.finish();
}

fn bench_all_greeks(c: &mut Criterion) {
    let mut group = c.benchmark_group("All Greeks Calculation");

    // Configure the benchmark group
    group.warm_up_time(Duration::from_millis(500));
    group.measurement_time(Duration::from_secs(2));

    // Compare individual vs. all_greeks calculation for different option types
    let call_atm = generate_standard_inputs(
        OptionType::Call,
        Moneyness::AtTheMoney,
        TimeToMaturity::MediumTerm,
        VolatilityLevel::Medium,
    );

    let put_atm = generate_standard_inputs(
        OptionType::Put,
        Moneyness::AtTheMoney,
        TimeToMaturity::MediumTerm,
        VolatilityLevel::Medium,
    );

    // Benchmark each individual Greek calculation separately (sum of times)
    for (name, inputs) in [("call", call_atm), ("put", put_atm)] {
        group.bench_with_input(
            BenchmarkId::new("individual", name),
            &inputs,
            |b, inputs: &Inputs| {
                b.iter(|| {
                    let delta = black_box(inputs.calc_delta().unwrap());
                    let gamma = black_box(inputs.calc_gamma().unwrap());
                    let theta = black_box(inputs.calc_theta().unwrap());
                    let vega = black_box(inputs.calc_vega().unwrap());
                    let rho = black_box(inputs.calc_rho().unwrap());
                    black_box((delta, gamma, theta, vega, rho))
                })
            },
        );

        // Benchmark all_greeks calculation (calculates all at once)
        group.bench_with_input(
            BenchmarkId::new("all_greeks", name),
            &inputs,
            |b, inputs: &Inputs| {
                b.iter(|| {
                    let all = black_box(inputs.calc_all_greeks().unwrap());
                    black_box(all)
                })
            },
        );
    }

    group.finish();
}

// Second-order Greeks benchmark
fn bench_second_order_greeks(c: &mut Criterion) {
    let mut group = c.benchmark_group("Second Order Greeks");

    // Configure the benchmark group
    group.warm_up_time(Duration::from_millis(500));
    group.measurement_time(Duration::from_secs(2));

    let inputs = generate_standard_inputs(
        OptionType::Call,
        Moneyness::AtTheMoney,
        TimeToMaturity::MediumTerm,
        VolatilityLevel::Medium,
    );

    group.bench_function("vanna", |b| {
        b.iter(|| black_box(inputs.calc_vanna().unwrap()))
    });

    group.bench_function("charm", |b| {
        b.iter(|| black_box(inputs.calc_charm().unwrap()))
    });

    group.bench_function("vomma", |b| {
        b.iter(|| black_box(inputs.calc_vomma().unwrap()))
    });

    group.bench_function("speed", |b| {
        b.iter(|| black_box(inputs.calc_speed().unwrap()))
    });

    group.bench_function("zomma", |b| {
        b.iter(|| black_box(inputs.calc_zomma().unwrap()))
    });

    group.finish();
}

criterion_group!(
    name = benches;
    config = Criterion::default().with_plots();
    targets = bench_greeks, bench_all_greeks, bench_second_order_greeks
);
criterion_main!(benches);
