const FRAC_1_SQRT_PI: f64 = 0.564_189_583_547_756_3_f64;

const THRESHOLD: f64 = 0.46875;
const XNEG: f64 = -26.6287357137514;
const XBIG: f64 = 26.543;

const A: [f64; 5] = [
    3.1611237438705656,
    113.864_154_151_050_16,
    377.485_237_685_302,
    3_209.377_589_138_469_4,
    0.185_777_706_184_603_15,
];
const B: [f64; 4] = [
    23.601_290_952_344_122,
    244.024_637_934_444_17,
    1_282.616_526_077_372_3,
    2_844.236_833_439_171,
];
const C: [f64; 9] = [
    0.564_188_496_988_670_1,
    8.883_149_794_388_377,
    66.119_190_637_141_63,
    298.635_138_197_400_1,
    881.952_221_241_769,
    1_712.047_612_634_070_7,
    2_051.078_377_826_071_6,
    1_230.339_354_797_997_2,
    2.153_115_354_744_038_3E-8,
];
const D: [f64; 8] = [
    15.744_926_110_709_835,
    117.693_950_891_312_5,
    537.181_101_862_009_9,
    1_621.389_574_566_690_3,
    3_290.799_235_733_459_7,
    4_362.619_090_143_247,
    3_439.367_674_143_721_6,
    1_230.339_354_803_749_5,
];
const P: [f64; 6] = [
    0.305_326_634_961_232_36,
    0.360_344_899_949_804_45,
    0.125_781_726_111_229_26,
    0.016_083_785_148_742_275,
    6.587_491_615_298_378E-4,
    0.016_315_387_137_302_097,
];
const Q: [f64; 5] = [
    2.568_520_192_289_822,
    1.872_952_849_923_460_4,
    0.527_905_102_951_428_5,
    0.060_518_341_312_441_32,
    0.002_335_204_976_268_691_8,
];

fn d_int(x: f64) -> f64 {
    if x > 0.0 {
        x.floor()
    } else {
        -(-x).floor()
    }
}

fn smoothened_exponential_of_negative_square(y: f64) -> f64 {
    let y_tilde = d_int(y * 16.0) / 16.0;
    (-y_tilde * y_tilde).exp() * (-(y - y_tilde) * (y + y_tilde)).exp()
}

fn smoothened_exponential_of_positive_square(x: f64) -> f64 {
    let x_tilde = d_int(x * 16.0) / 16.0;
    (x_tilde * x_tilde).exp() * ((x - x_tilde) * (x + x_tilde)).exp()
}

fn ab(z: f64) -> f64 {
    ((((A[4] * z + A[0]) * z + A[1]) * z + A[2]) * z + A[3])
        / ((((z + B[0]) * z + B[1]) * z + B[2]) * z + B[3])
}

fn cd(y: f64) -> f64 {
    ((((((((C[8] * y + C[0]) * y + C[1]) * y + C[2]) * y + C[3]) * y + C[4]) * y + C[5]) * y
        + C[6])
        * y
        + C[7])
        / ((((((((y + D[0]) * y + D[1]) * y + D[2]) * y + D[3]) * y + D[4]) * y + D[5]) * y + D[6])
            * y
            + D[7])
}

fn pq(z: f64) -> f64 {
    z * (((((P[5] * z + P[0]) * z + P[1]) * z + P[2]) * z + P[3]) * z + P[4])
        / (((((z + Q[0]) * z + Q[1]) * z + Q[2]) * z + Q[3]) * z + Q[4])
}

pub fn erfc(x: f64) -> f64 {
    let y = x.abs();
    if y <= THRESHOLD {
        return 1.0 - x * ab(y * y);
    }
    let erfc_abs_x = if y >= XBIG {
        0.0
    } else if y <= 4.0 {
        cd(y)
    } else {
        (FRAC_1_SQRT_PI - pq(1.0 / (y * y))) / y
    } * smoothened_exponential_of_negative_square(y);
    if x < 0.0 {
        2.0 - erfc_abs_x
    } else {
        erfc_abs_x
    }
}

fn erfcx_above_threshold(y: f64) -> f64 {
    debug_assert!(!(y <= THRESHOLD));
    if y <= 4.0 {
        cd(y)
    } else {
        (FRAC_1_SQRT_PI - pq(1.0 / (y * y))) / y
    }
}

pub fn erfcx(x: f64) -> f64 {
    let y = x.abs();
    if y <= THRESHOLD {
        let z = y * y;
        return z.exp() * (1.0 - x * ab(z));
    }
    if x < XNEG {
        return f64::MAX;
    }
    let result = erfcx_above_threshold(y);
    if x < 0.0 {
        let expx2 = smoothened_exponential_of_positive_square(x);
        (expx2 + expx2) - result
    } else {
        result
    }
}
