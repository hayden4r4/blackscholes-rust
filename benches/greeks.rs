use std::hint::black_box;

use blackscholes::{Greeks, Inputs, OptionType};
use criterion::{criterion_group, criterion_main, Criterion};

const INPUTS: Inputs = Inputs {
    option_type: OptionType::Call,
    s: 51.03,
    k: 55.0,
    p: None,
    r: 0.05,
    q: 0.02,
    t: 30.0 / 365.25,
    sigma: Some(0.3),
};

fn criterion_benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("Greeks");

    group.bench_function("delta", |b| b.iter(|| black_box(INPUTS.calc_delta())));
    group.bench_function("gamma", |b| b.iter(|| black_box(INPUTS.calc_gamma())));
    group.bench_function("theta", |b| b.iter(|| black_box(INPUTS.calc_theta())));
    group.bench_function("vega", |b| b.iter(|| black_box(INPUTS.calc_vega())));
    group.bench_function("rho", |b| b.iter(|| black_box(INPUTS.calc_rho())));
    group.bench_function("all_greeks", |b| {
        b.iter(|| black_box(INPUTS.calc_all_greeks()))
    });

    group.finish();
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
