use std::path::Path;

fn main() {
    let dir = Path::new("./lib");
    let build_dir = dir.join("build");

    println!(r"cargo:rerun-if-changed={}", dir.display());
    println!(r"cargo:rustc-env=NUM_JOBS=6");

    cc::Build::new()
        .files(&[
            dir.join("erf_cody.cpp"),
            dir.join("rationalcubic.cpp"),
            dir.join("normaldistribution.cpp"),
            dir.join("lets_be_rational.cpp"),
        ])
        .cpp(true)
        .flag("-fPIC")
        .flag("-DNDEBUG")
        .flag("-Ofast")
        .debug(false)
        .shared_flag(true)
        .opt_level(3)
        // .out_dir(&build_dir)
        .compile("liblets_be_rational");

    println!(
        r"cargo:rustc-link-search={}",
        build_dir.join("liblets_be_rational.a").display()
    );
}
