use blackscholes::{Greeks, Inputs, OptionType, Pricing, ImpliedVolatility};
// use std::time::{Duration, Instant};
// // use crossbeam_channel::bounded;
// use crossbeam_channel::bounded;
// use crossbeam_utils::thread;
// use rayon;
// use std::sync::Arc;

// fn test() {
//     let mut times: Vec<Duration> = Vec::new();

//     let n_workers: usize = 1;
//     let n_jobs: usize = 10_000_000;

//     let pool = rayon::ThreadPoolBuilder::new()
//         .num_threads(n_workers)
//         .build()
//         .unwrap();

//     // let pool = ThreadPool::new(n_workers);

//     let inputs: blackscholes::Inputs = blackscholes::Inputs {
//         option_type: blackscholes::OptionType::Call,
//         s: 3901.36,
//         k: 3900.0,
//         p: Some(74.10),
//         r: 0.0278,
//         q: 0.0137,
//         t: 12.0 / 360.0,
//         sigma: Some(0.254),
//     };
//     let inputs_arc = Arc::new(inputs);

//     // let (tx, rx) = bounded(n_jobs);

//     for _ in 0..n_jobs {
//         let now = Instant::now();

//         let inputs_clone = inputs_arc.clone();
//         // let tx_clone = tx.clone();

//         pool.spawn(move || {
//             let price: f64 = blackscholes::Greeks::calc_vega(&*inputs_clone);
//             // tx_clone.send(1).unwrap();
//         });
//         // println!("{:.1$}", price, 4);

//         let elapsed = now.elapsed();
//         times.push(elapsed);
//     }

//     // drop(tx);
//     // let count = rx.iter().count();
//     // println!("{}", count);

//     let mut sum_time = 0.0;
//     for time in &times {
//         sum_time += time.as_secs_f64();
//     }
//     let mut len: f64 = 0.0;
//     for _ in &times {
//         len += 1.0;
//     }
//     let mean_time: f64 = sum_time / len;
//     println!("{:.64}", sum_time);
//     println!("{} iters", n_jobs);
// }

fn main() {
    let inputs: blackscholes::Inputs = blackscholes::Inputs {
        option_type: blackscholes::OptionType::Call,
        s: 105.0,
        k: 100.0,
        p: Some(10.0),
        r: 0.05,
        q: 0.05,
        t: 30.0 / 365.25,
        sigma: None,
    };

    // let vanna: f64 = inputs.calc_price();
    // println!("{}", vanna);

    // let mut v: Vec<f64> = Vec::with_capacity(10_000_000);

    // for _ in 0..100_000_000 {
        let iv: f32 = inputs.calc_iv(0.0001).unwrap();
        println!("{}", iv);
    //     // v.push(price);
    // }


    // println!("{}", price);
    // //println!("{}", inputs)

    // let v: bool = inputs.option_type == blackscholes::OptionType::Call;
    // println!("{}", v)
}
