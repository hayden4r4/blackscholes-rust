#include <cmath>

double erf_std(double x) {
    return std::erf(x);
}

double erfc_std(double x) {
    return std::erfc(x);
}

double erfcx_std(double x) {
    // For this function, the standard C++ library does not provide a direct equivalent.
    // However, we can simulate its operation using the definition erfcx(x) = exp(x^2) * erfc(x).
    return std::exp(x * x) * std::erfc(x);
}