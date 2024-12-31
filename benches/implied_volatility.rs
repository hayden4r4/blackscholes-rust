use blackscholes::{OptionType, ImpliedVolatility, valuators::black_scholes::Inputs};
use criterion::{black_box, criterion_group, criterion_main, Criterion};

const INPUTS: Inputs = Inputs {
    option_type: OptionType::Call,
    s: 51.03,
    k: 55.0,
    p: Some(1.24),
    r: 0.0,
    q: 0.0,
    t: 45.0 / 365.25,
    sigma: None,
};

fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("calc_rational_iv", |b| {
        b.iter(|| black_box(INPUTS.calc_iv().unwrap()))
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
