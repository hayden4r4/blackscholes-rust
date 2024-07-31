use blackscholes::OptionType;
use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use rand::Rng;

#[allow(dead_code)]
fn abs_vs_equality_benchmarks(c: &mut Criterion) {
    let test_cases = vec![
        (1.0f64, 1.0f64),                // Identical values
        (1.0f64, 1.0f64 + f64::EPSILON), // Nearly identical values
        (1.0f64, 2.0f64),                // Significantly different values
    ];

    for (i, (x, y)) in test_cases.iter().enumerate() {
        let bench_name = format!("abs <= 0: case {}", i);
        c.bench_function(&bench_name, |b| {
            b.iter(|| {
                let result = black_box(*x).abs() <= 0.0;
                black_box(result);
            });
        });

        let bench_name = format!("equality: case {}", i);
        c.bench_function(&bench_name, |b| {
            b.iter(|| {
                let result = black_box(*x) == black_box(*y);
                black_box(result);
            });
        });
    }
}

fn generate_random_test_cases(n: usize) -> Vec<(OptionType, f64)> {
    let mut rng = rand::thread_rng();
    (0..n)
        .map(|_| {
            let option_type = if rng.gen_bool(0.5) {
                OptionType::Call
            } else {
                OptionType::Put
            };
            let value = rng.gen_range(-100.0..100.0);
            (option_type, value)
        })
        .collect()
}

fn mul_benchmarks(c: &mut Criterion) {
    let mut group = c.benchmark_group("OptionType Multiplication");
    let test_cases = generate_random_test_cases(100);

    group.bench_function("mul_impl", |b| {
        b.iter(|| {
            for &(x, y) in test_cases.iter() {
                black_box(x * y);
            }
        })
    });

    group.bench_function("cast_mul", |b| {
        b.iter(|| {
            for &(x, y) in test_cases.iter() {
                black_box(x as i8 as f64 * y);
            }
        })
    });

    for size in [10, 100, 1000].iter() {
        let test_cases = generate_random_test_cases(*size);

        group.bench_with_input(
            BenchmarkId::new("mul_impl_size", size),
            &test_cases,
            |b, tc| {
                b.iter(|| {
                    for &(x, y) in tc.iter() {
                        black_box(x * y);
                    }
                })
            },
        );

        group.bench_with_input(
            BenchmarkId::new("cast_mul_size", size),
            &test_cases,
            |b, tc| {
                b.iter(|| {
                    for &(x, y) in tc.iter() {
                        black_box(x as i8 as f64 * y);
                    }
                })
            },
        );
    }

    group.finish();
}

criterion_group!(benches, mul_benchmarks);
// criterion_group!(benches, abs_vs_equality_benchmarks);
criterion_main!(benches);
