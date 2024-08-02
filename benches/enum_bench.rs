use criterion::{black_box, Criterion, criterion_group, criterion_main};
use rand::Rng;

use blackscholes::OptionType;

pub enum Testable {
    Enum,
    String,
}

fn generate_random_inputs(n: usize) -> Vec<(OptionType, f64, f64, f64, f64, f64)> {
    let mut rng = rand::thread_rng();
    (0..n)
        .map(|_| {
            let option_type = if rng.gen_bool(0.5) {
                OptionType::Call
            } else {
                OptionType::Put
            };
            let s = rng.gen_range(50.0..150.0);
            let k = rng.gen_range(50.0..150.0);
            let r = rng.gen_range(0.0..0.1);
            let q = rng.gen_range(0.0..0.1);
            let t = rng.gen_range(0.0..1.0);
            (option_type, s, k, r, q, t)
        })
        .collect()
}

fn calc_price_enum(
    s: f64,
    k: f64,
    r: f64,
    t: f64,
) -> Result<f64, Testable> {
    if t == 0.0 {
        return Err(Testable::Enum);
    }
    Ok(s / r / k)
}

fn calc_price_string(
    s: f64,
    k: f64,
    r: f64,
    t: f64,
) -> Result<f64, String> {
    if t == 0.0 {
        return Err("Time to maturity is 0".to_string());
    }
    Ok(s / r / k)
}

fn enum_vs_string_benchmarks(c: &mut Criterion) {
    let test_cases = generate_random_inputs(1000);

    c.bench_function("calc_price_enum", |b| {
        b.iter(|| {
            for &(option_type, s, k, r, q, t) in test_cases.iter() {
                let _ = black_box(calc_price_enum(s, k, r, q * 0.0));
            }
        })
    });

    c.bench_function("calc_price_string", |b| {
        b.iter(|| {
            for &(option_type, s, k, r, q, t) in test_cases.iter() {
                let _ = black_box(calc_price_string(s, k, r, q * 0.0));
            }
        })
    });
}

criterion_group!(benches, enum_vs_string_benchmarks);
criterion_main!(benches);