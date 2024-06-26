#include "erf_cody.h"
#include <cmath>
#include <cassert>
#include <cstdio>
#include <chrono>
#include <vector>
#include <fstream>

using namespace std;
using namespace std::chrono;

bool test_erf_cody_for_value(double x) {
    if (abs(erf_cody(x) - erf(x)) < 1e-9) {
        return true;
    }
    return false;
}

vector<bool> run_unit_tests() {
    const double test_values[] = {0.0, 0.5, 1.0, -0.5, -1.0};
    const size_t num_test_values = sizeof(test_values) / sizeof(test_values[0]);

    vector<bool> unit_test_results;
    for (size_t i = 0; i < num_test_values; ++i) {
        unit_test_results.push_back(test_erf_cody_for_value(test_values[i]));
    }
    return unit_test_results;
}

struct PerformanceResult {
    float x;
    long long duration_cody;
    long long duration_std;
};

void print_summary(const vector<bool> &unit_test_results, const vector<PerformanceResult> &performance_results,
                   long long total_duration_cody, long long total_duration_std) {
    int passed_tests = count(unit_test_results.begin(), unit_test_results.end(), true);
    printf("Unit Test Summary: %d/%lu tests passed.\n", passed_tests, unit_test_results.size());

    int cody_wins = 0;
    int std_wins = 0;
    for (const auto &result: performance_results) {
        if (result.duration_cody < result.duration_std) {
            ++cody_wins;
        } else if (result.duration_cody > result.duration_std) {
            ++std_wins;
        }
    }

    printf("Performance Test Summary: Cody wins: %d, Std wins: %d\n", cody_wins, std_wins);
    printf("Total execution time (Cody): %lld nanoseconds, (Std): %lld nanoseconds\n", total_duration_cody,
           total_duration_std);
}

struct PerformanceTestResults {
    vector<PerformanceResult> results;
    long long total_duration_cody;
    long long total_duration_std;
};

PerformanceTestResults run_performance_tests() {
    const float start = -4.0f;
    const float end = 4.0f;
    const float step = 0.001f;

    PerformanceTestResults test_results;
    test_results.total_duration_cody = 0;
    test_results.total_duration_std = 0;

    try {
        for (float x = start; x <= end; x += step) {
            PerformanceResult result;
            result.x = x;

            auto start_time = high_resolution_clock::now();
            erf_cody(x);
            auto end_time = high_resolution_clock::now();
            result.duration_cody = duration_cast<nanoseconds>(end_time - start_time).count();
            test_results.total_duration_cody += result.duration_cody;

            start_time = high_resolution_clock::now();
            erf(x);
            end_time = high_resolution_clock::now();
            result.duration_std = duration_cast<nanoseconds>(end_time - start_time).count();
            test_results.total_duration_std += result.duration_std;

            test_results.results.push_back(result);
        }
    } catch (const exception &e) {
        printf("An error occurred: %s\n", e.what());
    }

    return test_results;
}

int main() {
    auto unit_test_results = run_unit_tests();
    auto performance_test_results = run_performance_tests();
    print_summary(unit_test_results, performance_test_results.results, performance_test_results.total_duration_cody,
                  performance_test_results.total_duration_std);
    return 0;
}
