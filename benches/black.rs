extern crate criterion;
use blackscholes::lets_be_rational::black::{
    asymptotic_expansion_of_normalised_black_call,
    asymptotic_expansion_of_normalised_black_call_old,
};
use criterion::{black_box, criterion_group, criterion_main, Criterion};

fn benchmark_asymptotic_expansion_of_normalised_black_call_changed(c: &mut Criterion) {
    let test_values = vec![
        (-12.0, 0.1),
        (-20.0, 0.05),
        (-15.0, 0.2),
        (-10.0, 0.05),
        (-30.0, 0.01),
    ];

    for (h, t) in test_values {
        c.bench_function(
            &format!(
                "asymptotic_expansion_of_normalised_black_call_old h={}, t={}",
                h, t
            ),
            |b| {
                b.iter(|| {
                    asymptotic_expansion_of_normalised_black_call_old(black_box(h), black_box(t))
                })
            },
        );
    }
}

fn benchmark_asymptotic_expansion_of_normalised_black_call(c: &mut Criterion) {
    let test_values = vec![
        (-12.0, 0.1),
        (-20.0, 0.05),
        (-15.0, 0.2),
        (-10.0, 0.05),
        (-30.0, 0.01),
    ];

    for (h, t) in test_values {
        c.bench_function(
            &format!(
                "asymptotic_expansion_of_normalised_black_call_new_one h={}, t={}",
                h, t
            ),
            |b| {
                b.iter(|| asymptotic_expansion_of_normalised_black_call(black_box(h), black_box(t)))
            },
        );
    }
}

criterion_group!(
    benches,
    benchmark_asymptotic_expansion_of_normalised_black_call_changed,
    benchmark_asymptotic_expansion_of_normalised_black_call
);
criterion_main!(benches);
