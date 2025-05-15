use blackscholes::{Pricing, Greeks, Inputs};
use criterion::{black_box, criterion_group, criterion_main, Criterion, BenchmarkId, Throughput};
use std::time::Duration;
use rand::thread_rng;
use std::fs::File;
use std::io::Write;
use indicatif::{ProgressBar, ProgressStyle};

#[path = "../common.rs"]
mod common;
use common::{generate_random_inputs, BatchSize, get_sample_config};

// Define batch sizes to test from very small to very large
const BATCH_SIZES: [usize; 9] = [
    10,         // Tiny
    100,        // Small
    500,        // Medium-small
    1_000,      // Medium
    5_000,      // Medium-large
    10_000,     // Large
    25_000,     // Larger
    50_000,     // Very large
    100_000,    // Huge
];

// Define a struct to hold benchmark results for analysis
struct ScalingResult {
    batch_size: usize,
    total_time_ns: f64,
    time_per_option_ns: f64,
    throughput_ops_per_sec: f64,
}

// Benchmark scaling behavior with different batch sizes
fn bench_scaling(c: &mut Criterion) {
    let mut group = c.benchmark_group("Batch Size Scaling");
    
    // Configure the benchmark group
    group.warm_up_time(Duration::from_millis(500));
    
    let mut rng = thread_rng();
    let mut results = Vec::new();
    
    // Setup progress bar for user feedback
    let progress_bar = ProgressBar::new(BATCH_SIZES.len() as u64);
    progress_bar.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} {msg}")
            .unwrap()
            .progress_chars("#>-")
    );
    progress_bar.set_message("Running batch scaling benchmarks...");
    
    // Run benchmarks for each batch size
    for (idx, &size) in BATCH_SIZES.iter().enumerate() {
        progress_bar.set_position(idx as u64);
        progress_bar.set_message(format!("Benchmarking batch size: {}", size));
        
        // Adjust sample count and measurement time based on batch size
        let (sample_count, measurement_time) = get_sample_config(size);
        group.sample_size(sample_count);
        group.measurement_time(measurement_time);
        
        // Generate random inputs
        let inputs = generate_random_inputs(size, &mut rng);
        
        // Set throughput measurement for proper ops/sec calculation
        group.throughput(Throughput::Elements(size as u64));
        
        // Benchmark price calculation
        let id = BenchmarkId::new("price", size);
        group.bench_function(id.clone(), |b| {
            b.iter(|| {
                let mut results = Vec::with_capacity(inputs.len());
                for input in black_box(&inputs) {
                    results.push(input.calc_price().unwrap());
                }
                black_box(results)
            })
        });
        
        // Extract results for analysis
        let id_str = id.as_str().to_string();
        if let Some(estimate) = group.benchmark_id(&id_str) {
            let ns_per_iter = estimate.median.as_nanos() as f64;
            let ns_per_option = ns_per_iter / size as f64;
            let ops_per_sec = 1_000_000_000.0 / ns_per_option;
            
            results.push(ScalingResult {
                batch_size: size,
                total_time_ns: ns_per_iter,
                time_per_option_ns: ns_per_option,
                throughput_ops_per_sec: ops_per_sec,
            });
            
            println!("Batch size: {}, Time per option: {:.2} ns, Throughput: {:.2} million ops/sec", 
                size, ns_per_option, ops_per_sec / 1_000_000.0);
        }
    }
    
    group.finish();
    progress_bar.finish_with_message("Benchmarks complete");
    
    // Save results to CSV for external analysis
    save_results_to_csv(&results);
}

// Save benchmark results to CSV file for further analysis
fn save_results_to_csv(results: &[ScalingResult]) {
    if results.is_empty() {
        return;
    }
    
    let file_path = "target/scaling_results.csv";
    let mut file = match File::create(file_path) {
        Ok(file) => file,
        Err(e) => {
            eprintln!("Failed to create results file: {}", e);
            return;
        }
    };
    
    // Write CSV header
    if let Err(e) = writeln!(file, "batch_size,total_time_ns,time_per_option_ns,throughput_ops_per_sec") {
        eprintln!("Failed to write CSV header: {}", e);
        return;
    }
    
    // Write data rows
    for result in results {
        if let Err(e) = writeln!(
            file, 
            "{},{:.2},{:.2},{:.2}", 
            result.batch_size, 
            result.total_time_ns,
            result.time_per_option_ns,
            result.throughput_ops_per_sec
        ) {
            eprintln!("Failed to write result row: {}", e);
            return;
        }
    }
    
    println!("Results saved to {}", file_path);
}

criterion_group!(
    name = benches;
    config = Criterion::default().with_plots();
    targets = bench_scaling
);
criterion_main!(benches); 