use blackscholes::{ImpliedVolatility, Inputs, OptionType, Pricing};
use criterion::{black_box, criterion_group, criterion_main, Criterion, BenchmarkId};
use std::time::Duration;

#[path = "../common.rs"]
mod common;
use common::{generate_standard_inputs, Moneyness, TimeToMaturity, VolatilityLevel};

fn bench_implied_volatility(c: &mut Criterion) {
    let mut group = c.benchmark_group("Implied Volatility");
    
    // Configure the benchmark group
    group.warm_up_time(Duration::from_millis(500));
    group.measurement_time(Duration::from_secs(3)); // IV calculation takes longer
    
    // Different price scenarios to test IV calculation
    // We'll only use ATM options to avoid convergence issues
    let configurations = [
        ("call_atm", OptionType::Call, Moneyness::AtTheMoney),
        ("put_atm", OptionType::Put, Moneyness::AtTheMoney),
    ];
    
    // Generate inputs for benchmarking IV calculation
    for (name, option_type, moneyness) in configurations {
        // Create standard inputs
        let mut inputs = generate_standard_inputs(
            option_type, 
            moneyness, 
            TimeToMaturity::MediumTerm, 
            VolatilityLevel::Medium
        );
        
        // Calculate option price using inputs.sigma
        let sigma = inputs.sigma.unwrap();
        let price = inputs.calc_price().unwrap();
        
        // Set price and clear sigma for IV calculation
        inputs.p = Some(price);
        inputs.sigma = None;
        
        // Bench standard IV calculation with different tolerances
        for (tolerance_name, tolerance) in [
            ("high_precision", 0.00001),
            ("medium_precision", 0.0001),
            ("low_precision", 0.001),
        ] {
            group.bench_with_input(
                BenchmarkId::new(format!("calc_iv_{}", tolerance_name), name), 
                &(inputs.clone(), tolerance), 
                |b, (inputs, tolerance): &(Inputs, f64)| {
                    b.iter(|| {
                        let iv = black_box(inputs).calc_iv(*tolerance).unwrap();
                        black_box(iv)
                    })
                }
            );
        }
        
        // Bench rational IV calculation
        group.bench_with_input(
            BenchmarkId::new("calc_rational_iv", name), 
            &inputs, 
            |b, inputs: &Inputs| {
                b.iter(|| {
                    let iv = black_box(inputs).calc_rational_iv().unwrap();
                    black_box(iv)
                })
            }
        );
        
        // Compare to known result for verification
        let iv_standard = inputs.calc_iv(0.00001).unwrap();
        let iv_rational = inputs.calc_rational_iv().unwrap();
        
        println!(
            "Verification for {}: Original sigma = {:.6}, IV standard = {:.6}, IV rational = {:.6}, Diff = {:.6}",
            name, sigma, iv_standard, iv_rational, (sigma - iv_standard).abs()
        );
    }
    
    group.finish();
}

criterion_group!(
    name = benches;
    config = Criterion::default().with_plots();
    targets = bench_implied_volatility
);
criterion_main!(benches); 