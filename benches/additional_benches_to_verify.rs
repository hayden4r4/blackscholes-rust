use criterion::{black_box, criterion_group, criterion_main, Criterion};

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

criterion_group!(benches, abs_vs_equality_benchmarks);
criterion_main!(benches);
