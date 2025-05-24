use std::{hint::black_box, time::Duration};

use blackscholes::{Inputs, OptionType, Pricing};
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};

// Fix the module import using a path relative to the crate root
#[path = "../common/mod.rs"]
mod common;
use common::{generate_standard_inputs, Moneyness, TimeToMaturity, VolatilityLevel};

fn bench_option_pricing(c: &mut Criterion) {
    let mut group = c.benchmark_group("Single Option Pricing");

    // Configure the benchmark group
    group.warm_up_time(Duration::from_millis(500));
    group.measurement_time(Duration::from_secs(2));

    // Standard option configurations to test - comprehensive coverage
    let configurations = [
        (
            "call_atm_short",
            OptionType::Call,
            Moneyness::AtTheMoney,
            TimeToMaturity::ShortTerm,
            VolatilityLevel::Medium,
        ),
        (
            "put_atm_short",
            OptionType::Put,
            Moneyness::AtTheMoney,
            TimeToMaturity::ShortTerm,
            VolatilityLevel::Medium,
        ),
        (
            "call_itm_medium",
            OptionType::Call,
            Moneyness::InTheMoney,
            TimeToMaturity::MediumTerm,
            VolatilityLevel::Medium,
        ),
        (
            "put_itm_medium",
            OptionType::Put,
            Moneyness::InTheMoney,
            TimeToMaturity::MediumTerm,
            VolatilityLevel::Medium,
        ),
        (
            "call_otm_long",
            OptionType::Call,
            Moneyness::OutOfTheMoney,
            TimeToMaturity::LongTerm,
            VolatilityLevel::Medium,
        ),
        (
            "call_atm_short_highvol",
            OptionType::Call,
            Moneyness::AtTheMoney,
            TimeToMaturity::ShortTerm,
            VolatilityLevel::High,
        ),
    ];

    for (name, option_type, moneyness, maturity, vol_level) in configurations {
        let inputs = generate_standard_inputs(option_type, moneyness, maturity, vol_level);

        // Benchmark standard pricing
        group.bench_with_input(
            BenchmarkId::new("calc_price", name),
            &inputs,
            |b, inputs: &Inputs| {
                b.iter(|| {
                    let result = black_box(inputs).calc_price().unwrap();
                    black_box(result)
                })
            },
        );

        // Benchmark rational pricing
        group.bench_with_input(
            BenchmarkId::new("calc_rational_price", name),
            &inputs,
            |b, inputs: &Inputs| {
                b.iter(|| {
                    let result = black_box(inputs).calc_rational_price().unwrap();
                    black_box(result)
                })
            },
        );
    }

    group.finish();
}

// Skip the d1d2 calculation benchmark since it's using a private function
// We'll add a replacement that uses public functions later

criterion_group!(
    name = benches;
    config = Criterion::default().with_plots();
    targets = bench_option_pricing
);
criterion_main!(benches);
