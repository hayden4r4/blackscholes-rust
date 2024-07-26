const FRAC_1_SQRT_PI: f64 = 0.564189583547756286948079451560772586_f64;

const THRESHOLD: f64 = 0.46875;
const XNEG: f64 = -26.6287357137514;
const XBIG: f64 = 26.543;

const A: [f64; 5] = [
    3.1611237438705656,
    113.864154151050156,
    377.485237685302021,
    3209.37758913846947,
    0.185777706184603153,
];
const B: [f64; 4] = [
    23.6012909523441209,
    244.024637934444173,
    1282.61652607737228,
    2844.23683343917062,
];
const C: [f64; 9] = [
    0.564188496988670089,
    8.88314979438837594,
    66.1191906371416295,
    298.635138197400131,
    881.95222124176909,
    1712.04761263407058,
    2051.07837782607147,
    1230.33935479799725,
    2.15311535474403846E-8,
];
const D: [f64; 8] = [
    15.7449261107098347,
    117.693950891312499,
    537.181101862009858,
    1621.38957456669019,
    3290.79923573345963,
    4362.61909014324716,
    3439.36767414372164,
    1230.33935480374942,
];
const P: [f64; 6] = [
    0.305326634961232344,
    0.360344899949804439,
    0.125781726111229246,
    0.0160837851487422766,
    6.58749161529837803E-4,
    0.0163153871373020978,
];
const Q: [f64; 5] = [
    2.56852019228982242,
    1.87295284992346047,
    0.527905102951428412,
    0.0605183413124413191,
    0.00233520497626869185,
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
