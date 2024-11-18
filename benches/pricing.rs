use blackscholes::{Inputs, OptionType, Pricing};
use criterion::{black_box, criterion_group, criterion_main, Criterion};

const INPUTS: Inputs = Inputs {
    option_type: OptionType::Call,
    s: 51.03,
    k: 55.0,
    p: None,
    r: 0.0,
    q: 0.0,
    t: 25.0 / 360.0,
    sigma: Some(0.5),
};

fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("calc_price", |b| {
        b.iter(|| {
            let input = black_box(INPUTS);
            let result = input.calc_price().map_err(|e| black_box(e)).unwrap();
            black_box(result + 1.0);
        })
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
