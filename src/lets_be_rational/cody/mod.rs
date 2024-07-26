pub mod optimized;

#[cfg(test)]
mod tests {
    use std::time::Instant;

    use assert_approx_eq::assert_approx_eq;
    use statrs::function::erf::erfc;

    use crate::lets_be_rational::cody::optimized;

    fn measure_time<F>(f: F) -> f64
    where
        F: Fn(),
    {
        let start = Instant::now();
        f();
        start.elapsed().as_secs_f64()
    }

    #[test]
    fn test_accuracy_erfc() {
        let test_values = [-3.0, -1.5, -0.5, 0.0, 0.5, 1.5, 3.0];

        for &x in &test_values {
            let erfc_original = erfc(x);
            let erfc_optimized = optimized::erfc(x);
            println!(
                "x: {}, original: {}, optimized: {}, diff: {}",
                x,
                erfc_original,
                erfc_optimized,
                (erfc_original - erfc_optimized).abs()
            );
            assert_approx_eq!(erfc_original, erfc_optimized, 1e-10);
        }
    }

    #[test]
    fn test_accuracy_erfcx() {
        let test_values = [-3.0, -1.5, -0.5, 0.0, 0.5, 1.5, 3.0];

        for &x in &test_values {
            let erfcx_original = erfcx(x);
            let erfcx_optimized = optimized::erfcx(x);
            println!(
                "x: {}, original: {}, optimized: {}, diff: {}",
                x,
                erfcx_original,
                erfcx_optimized,
                (erfcx_original - erfcx_optimized).abs()
            );
            assert_approx_eq!(erfcx_original, erfcx_optimized, 1e-10);
        }
    }

    #[test]
    fn test_performance() {
        let iterations = 1_000_000;
        let x = 1.5;

        // erfc
        let time_original_erfc = measure_time(|| {
            for _ in 0..iterations {
                let _ = statrs::function::erf::erfc(x);
            }
        });
        let time_optimized_erfc = measure_time(|| {
            for _ in 0..iterations {
                let _ = optimized::erfc(x);
            }
        });
        println!(
            "erfc: original = {:.6} s, optimized = {:.6} s",
            time_original_erfc, time_optimized_erfc
        );

        // erfcx
        let time_original_erfcx = measure_time(|| {
            for _ in 0..iterations {
                let _ = erfcx(x);
            }
        });
        let time_optimized_erfcx = measure_time(|| {
            for _ in 0..iterations {
                let _ = optimized::erfcx(x);
            }
        });
        println!(
            "erfcx: original = {:.6} s, optimized = {:.6} s",
            time_original_erfcx, time_optimized_erfcx
        );

        assert!(true);
    }

    fn erfcx(x: f64) -> f64 {
        (x * x).exp() * erfc(x)
    }
}
