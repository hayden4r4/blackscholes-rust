use blackscholes::{Pricing, Greeks, Inputs};
use criterion::{black_box, criterion_group, criterion_main, Criterion};
use std::time::{Duration, Instant};
use rand::thread_rng;
use plotters::prelude::*;

#[path = "../common.rs"]
mod common;
use common::{generate_random_inputs, BatchSize};

// Define batch sizes for visualization
const BATCH_SIZES: [usize; 10] = [
    10, 50, 100, 500, 1_000, 5_000, 10_000, 25_000, 50_000, 100_000
];

// Colors for different operations
const PRICE_COLOR: RGBColor = RGBColor(0, 119, 182);
const RATIONAL_COLOR: RGBColor = RGBColor(0, 180, 216);
const DELTA_COLOR: RGBColor = RGBColor(144, 224, 239);
const GAMMA_COLOR: RGBColor = RGBColor(202, 240, 248);

// Simple measurement for a large number of batch sizes
// This doesn't use criterion as we want to collect a lot of data points quickly
fn measure_scaling() {
    println!("Running scaling visualization...");
    
    // Prepare data structures
    let mut rng = thread_rng();
    let mut time_per_option_results: Vec<(String, Vec<(f64, f64)>)> = Vec::new();
    let mut throughput_results: Vec<(String, Vec<(f64, f64)>)> = Vec::new();
    
    // First data series: Standard pricing
    let mut price_times = Vec::new();
    let mut price_throughput = Vec::new();
    
    // Second data series: Rational pricing
    let mut rational_times = Vec::new();
    let mut rational_throughput = Vec::new();
    
    // Third data series: Delta calculation  
    let mut delta_times = Vec::new();
    let mut delta_throughput = Vec::new();
    
    // Fourth data series: Gamma calculation
    let mut gamma_times = Vec::new();
    let mut gamma_throughput = Vec::new();
    
    // Run benchmark for each batch size
    for &size in BATCH_SIZES.iter() {
        let inputs = generate_random_inputs(size, &mut rng);
        
        // Measure pricing time
        let pricing_time = measure_operation_time(size, &inputs, |input| {
            input.calc_price().unwrap()
        });
        price_times.push((size as f64, pricing_time));
        price_throughput.push((size as f64, 1_000_000_000.0 / pricing_time));
        
        // Measure rational pricing time
        let rational_time = measure_operation_time(size, &inputs, |input| {
            input.calc_rational_price().unwrap()
        });
        rational_times.push((size as f64, rational_time));
        rational_throughput.push((size as f64, 1_000_000_000.0 / rational_time));
        
        // Measure delta calculation time
        let delta_time = measure_operation_time(size, &inputs, |input| {
            input.calc_delta().unwrap()
        });
        delta_times.push((size as f64, delta_time));
        delta_throughput.push((size as f64, 1_000_000_000.0 / delta_time));
        
        // Measure gamma calculation time
        let gamma_time = measure_operation_time(size, &inputs, |input| {
            input.calc_gamma().unwrap()
        });
        gamma_times.push((size as f64, gamma_time));
        gamma_throughput.push((size as f64, 1_000_000_000.0 / gamma_time));
        
        println!("Batch size: {}, price: {:.2} ns/op, rational: {:.2} ns/op, delta: {:.2} ns/op, gamma: {:.2} ns/op",
            size, pricing_time, rational_time, delta_time, gamma_time);
    }
    
    // Organize data for visualization
    time_per_option_results.push(("Standard Price".to_string(), price_times));
    time_per_option_results.push(("Rational Price".to_string(), rational_times));
    time_per_option_results.push(("Delta".to_string(), delta_times));
    time_per_option_results.push(("Gamma".to_string(), gamma_times));
    
    throughput_results.push(("Standard Price".to_string(), price_throughput));
    throughput_results.push(("Rational Price".to_string(), rational_throughput));
    throughput_results.push(("Delta".to_string(), delta_throughput));
    throughput_results.push(("Gamma".to_string(), gamma_throughput));
    
    // Generate visualizations
    match create_time_per_option_chart(&time_per_option_results) {
        Ok(_) => println!("Time per option chart created successfully"),
        Err(e) => eprintln!("Failed to create time per option chart: {}", e),
    }
    
    match create_throughput_chart(&throughput_results) {
        Ok(_) => println!("Throughput chart created successfully"),
        Err(e) => eprintln!("Failed to create throughput chart: {}", e),
    }
}

// Helper function to measure operation time
fn measure_operation_time<F>(size: usize, inputs: &[Inputs], operation: F) -> f64
where
    F: Fn(&Inputs) -> f64,
{
    const WARMUP_ITERS: usize = 5;
    const MEASUREMENT_ITERS: usize = 10;
    
    // Warmup phase
    for _ in 0..WARMUP_ITERS {
        let mut results = Vec::with_capacity(inputs.len());
        for input in inputs {
            results.push(operation(input));
        }
        black_box(results);
    }
    
    // Measurement phase
    let mut total_duration = Duration::ZERO;
    
    for _ in 0..MEASUREMENT_ITERS {
        let start = Instant::now();
        let mut results = Vec::with_capacity(inputs.len());
        for input in inputs {
            results.push(operation(input));
        }
        black_box(results);
        total_duration += start.elapsed();
    }
    
    // Calculate average time per option in nanoseconds
    let avg_time_ns = total_duration.as_nanos() as f64 / (MEASUREMENT_ITERS as f64 * size as f64);
    avg_time_ns
}

// Create a chart showing time per option for different batch sizes
fn create_time_per_option_chart(data: &[(String, Vec<(f64, f64)>)]) -> Result<(), Box<dyn std::error::Error>> {
    let root = BitMapBackend::new("target/time_per_option.png", (1024, 768)).into_drawing_area();
    root.fill(&WHITE)?;
    
    // Find y-axis range
    let min_y = data.iter()
        .flat_map(|(_, points)| points.iter().map(|(_, y)| *y))
        .fold(f64::INFINITY, |acc: f64, y| acc.min(y));
    
    let max_y = data.iter()
        .flat_map(|(_, points)| points.iter().map(|(_, y)| *y))
        .fold(0.0_f64, |acc: f64, y| acc.max(y));
    
    let y_range = (min_y * 0.8)..=(max_y * 1.2);
    
    // Create chart
    let mut chart = ChartBuilder::on(&root)
        .caption("Time per Option vs Batch Size", ("sans-serif", 30).into_font())
        .margin(10)
        .x_label_area_size(40)
        .y_label_area_size(60)
        .build_cartesian_2d(
            (10.0..100_000.0).log_scale(),
            y_range,
        )?;
    
    chart.configure_mesh()
        .x_desc("Batch Size (log scale)")
        .y_desc("Time per Option (ns)")
        .draw()?;
    
    // Define colors for each series
    let colors = [PRICE_COLOR, RATIONAL_COLOR, DELTA_COLOR, GAMMA_COLOR];
    
    // Plot each data series
    for (idx, (name, points)) in data.iter().enumerate() {
        let color = colors[idx % colors.len()];
        
        // Plot the line
        chart.draw_series(LineSeries::new(
            points.iter().copied(),
            color.stroke_width(3)
        ))?
        .label(name)
        .legend(move |(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], color.stroke_width(3)));
    }
    
    // Draw legend
    chart.configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .position(SeriesLabelPosition::UpperRight)
        .draw()?;
    
    println!("Time per option chart saved to target/time_per_option.png");
    Ok(())
}

// Create a chart showing throughput for different batch sizes
fn create_throughput_chart(data: &[(String, Vec<(f64, f64)>)]) -> Result<(), Box<dyn std::error::Error>> {
    let root = BitMapBackend::new("target/throughput.png", (1024, 768)).into_drawing_area();
    root.fill(&WHITE)?;
    
    // Find y-axis range
    let min_y = data.iter()
        .flat_map(|(_, points)| points.iter().map(|(_, y)| *y))
        .fold(f64::INFINITY, |acc: f64, y| acc.min(y));
    
    let max_y = data.iter()
        .flat_map(|(_, points)| points.iter().map(|(_, y)| *y))
        .fold(0.0_f64, |acc: f64, y| acc.max(y));
    
    let y_range = (min_y * 0.8)..=(max_y * 1.2);
    
    // Create chart
    let mut chart = ChartBuilder::on(&root)
        .caption("Throughput vs Batch Size", ("sans-serif", 30).into_font())
        .margin(10)
        .x_label_area_size(40)
        .y_label_area_size(60)
        .build_cartesian_2d(
            (10.0..100_000.0).log_scale(),
            y_range,
        )?;
    
    chart.configure_mesh()
        .x_desc("Batch Size (log scale)")
        .y_desc("Throughput (ops/sec)")
        .draw()?;
    
    // Define colors for each series
    let colors = [PRICE_COLOR, RATIONAL_COLOR, DELTA_COLOR, GAMMA_COLOR];
    
    // Plot each data series
    for (idx, (name, points)) in data.iter().enumerate() {
        let color = colors[idx % colors.len()];
        
        chart.draw_series(LineSeries::new(
            points.iter().map(|&(x, y)| (x, y / 1_000_000.0)), // Convert to M ops/sec
            color.stroke_width(3)
        ))?
        .label(format!("{} (M ops/sec)", name))
        .legend(move |(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], color.stroke_width(3)));
    }
    
    // Draw legend
    chart.configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .position(SeriesLabelPosition::LowerRight)
        .draw()?;
    
    println!("Throughput chart saved to target/throughput.png");
    Ok(())
}

// Simple benchmark function that just triggers our measurement
fn bench_visualization(_c: &mut Criterion) {
    measure_scaling();
}

criterion_group!(benches, bench_visualization);
criterion_main!(benches); 