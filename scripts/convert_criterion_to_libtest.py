#!/usr/bin/env python3
"""
Convert Criterion benchmark output to libtest format for github-action-benchmark
"""
import re
import sys

def convert_criterion_to_libtest(input_file, output_file):
    # Regex to match Criterion output format
    # Example: calc_price              time:   [26.608 ns 26.696 ns 26.858 ns]
    # Example: Greeks/delta            time:   [863.41 ps 873.28 ps 887.98 ps]
    # Example: Single Option Greeks/delta/short_term
    #                        time:   [26.970 ns 27.019 ns 27.083 ns]
    criterion_pattern = r'^(.+?)\s+time:\s+\[([0-9.]+)\s+(\w+)\s+([0-9.]+)\s+(\w+)\s+([0-9.]+)\s+(\w+)\]'
    
    with open(input_file, 'r') as f:
        content = f.read()
    
    libtest_output = []
    
    for line in content.split('\n'):
        match = re.match(criterion_pattern, line.strip())
        if match:
            bench_name = match.group(1).strip()
            min_time = float(match.group(2))
            min_unit = match.group(3)
            avg_time = float(match.group(4))
            avg_unit = match.group(5)
            max_time = float(match.group(6))
            max_unit = match.group(7)
            
            # Clean benchmark name - replace spaces and slashes with underscores for libtest compatibility
            clean_name = re.sub(r'[/\s]+', '_', bench_name).strip('_')
            
            # Convert to nanoseconds if needed and calculate range
            if avg_unit == 'ps':
                # For picoseconds, keep some precision to avoid 0 values
                value_ns = avg_time / 1000
                if value_ns < 1:
                    # Show as fractional nanoseconds for very fast benchmarks
                    value = f"{value_ns:.2f}".rstrip('0').rstrip('.')
                else:
                    value = int(value_ns)
                range_val = max(0.01, (max_time - min_time) / 1000 / 2)
                unit = 'ns/iter'
            elif avg_unit == 'ns':
                value = int(avg_time)
                range_val = max(1, int((max_time - min_time) / 2 * 100) / 100)  # Keep 2 decimal places, min 1
                unit = 'ns/iter'
            elif avg_unit == 'μs' or avg_unit == 'us':
                value = int(avg_time * 1000)
                range_val = max(1, int((max_time - min_time) * 1000 / 2))
                unit = 'ns/iter'
            elif avg_unit == 'ms':
                value = int(avg_time * 1_000_000)
                range_val = max(1, int((max_time - min_time) * 1_000_000 / 2))
                unit = 'ns/iter'
            elif avg_unit == 's':
                value = int(avg_time * 1_000_000_000)
                range_val = max(1, int((max_time - min_time) * 1_000_000_000 / 2))
                unit = 'ns/iter'
            else:
                # Default to ns
                value = int(avg_time)
                range_val = max(1, int((max_time - min_time) / 2 * 100) / 100)
                unit = 'ns/iter'
            
            # Format range value appropriately
            if range_val == int(range_val):
                range_str = str(int(range_val))
            else:
                range_str = f"{range_val:.2f}".rstrip('0').rstrip('.')
            
            # Format as libtest expects: test name ... bench: value unit (+/- range)
            libtest_line = f"test {clean_name} ... bench: {value} {unit} (+/- {range_str})"
            libtest_output.append(libtest_line)
    
    # Write the converted output
    with open(output_file, 'w') as f:
        f.write('\n'.join(libtest_output))
        if libtest_output:  # Add final newline if we have content
            f.write('\n')
    
    print(f"Converted {len(libtest_output)} benchmarks from Criterion to libtest format")
    return len(libtest_output)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 convert_criterion_to_libtest.py <input_file> <output_file>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    try:
        count = convert_criterion_to_libtest(input_file, output_file)
        if count == 0:
            print("Warning: No Criterion benchmarks found in input file")
            sys.exit(1)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1) 