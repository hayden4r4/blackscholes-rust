// NOTE: To run current bench you need to make pub `mod black;` line in src/lets_be_rational/mod.rs
// and uncomment bench in toml file
// and uncomment the following code

// extern crate criterion;
//
//
// use criterion::{black_box, Criterion, criterion_group};
//
// use blackscholes::lets_be_rational::black::{
//     asymptotic_expansion_of_normalised_black_call, small_t_expansion_of_normalised_black_call,
// };
//
// fn benchmark_asymptotic_expansion_of_normalised_black_call(c: &mut Criterion) {
//     let test_values = vec![
//         (-12.0, 0.1),
//         (-20.0, 0.05),
//         (-15.0, 0.2),
//         (-10.0, 0.05),
//         (-30.0, 0.01),
//     ];
//
//     for (h, t) in test_values {
//         c.bench_function(
//             &format!(
//                 "asymptotic_expansion_of_normalised_black_call_new_one h={}, t={}",
//                 h, t
//             ),
//             |b| {
//                 b.iter(|| asymptotic_expansion_of_normalised_black_call(black_box(h), black_box(t)))
//             },
//         );
//     }
// }
//
// fn benchmark_small_t_expansion(c: &mut Criterion) {
//     let test_values = vec![
//         (0.1, 0.1),
//         (0.05, 0.05),
//         (0.15, 0.15),
//         (0.2, 0.2),
//         (0.0, 0.0),
//     ];
//
//     for (h, t) in test_values {
//         c.bench_function(
//             &format!(
//                 "small_t_expansion_of_normalised_black_call_new h={}, t={}",
//                 h, t
//             ),
//             |b| b.iter(|| small_t_expansion_of_normalised_black_call(black_box(h), black_box(t))),
//         );
//     }
// }
//
// criterion_group!(
//     benches,
//     benchmark_asymptotic_expansion_of_normalised_black_call,
//     benchmark_small_t_expansion
// );
