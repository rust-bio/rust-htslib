// Copyright 2014 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

#[cfg(feature = "serde")]
use bindgen;
use cc;
use fs_utils::copy::copy_directory;
use glob::glob;

use std::env;
use std::fs;
use std::path::PathBuf;
use std::process::Command;

fn sed_htslib_makefile(out: &PathBuf, patterns: &[&str], feature: &str) {
    for pattern in patterns {
        if !Command::new("sed")
            .current_dir(out.join("htslib"))
            .arg("-i")
            .arg("-e")
            .arg(pattern)
            .arg("Makefile")
            .status()
            .unwrap()
            .success()
        {
            panic!("failed to strip {} support", feature);
        }
    }
}

fn main() {
    let out = PathBuf::from(env::var("OUT_DIR").unwrap());
    let mut cfg = cc::Build::new();
    let want_static = cfg!(feature = "static") || env::var("HTS_STATIC").is_ok();

    if want_static {
        cfg.warnings(true).static_flag(true).pic(true);
    } else {
        cfg.warnings(true).static_flag(false).pic(true); 
    }

    if let Ok(z_inc) = env::var("DEP_Z_INCLUDE") {
        cfg.include(z_inc);
    }

    if !out.join("htslib").exists() {
        copy_directory("htslib", &out).unwrap();
    }

    let use_bzip2 = env::var("CARGO_FEATURE_BZIP2").is_ok();
    if !use_bzip2 {
        let bzip2_patterns = vec!["s/ -lbz2//", "/#define HAVE_LIBBZ2/d"];
        sed_htslib_makefile(&out, &bzip2_patterns, "bzip2");
    } else if let Ok(inc) = env::var("DEP_BZIP2_ROOT")
        .map(PathBuf::from)
        .map(|path| path.join("include"))
    {
        cfg.include(inc);
    }

    let use_lzma = env::var("CARGO_FEATURE_LZMA").is_ok();
    if !use_lzma && want_static {
        let lzma_patterns = vec!["s/ -llzma//", "/#define HAVE_LIBLZMA/d"];
        sed_htslib_makefile(&out, &lzma_patterns, "lzma");
    } else if let Ok(inc) = env::var("DEP_LZMA_INCLUDE").map(PathBuf::from) {
        cfg.include(inc);
    }

    let use_curl = env::var("CARGO_FEATURE_CURL").is_ok();
    if !use_curl {
        let curl_patterns = vec!["s/ -lcurl//", "/#define HAVE_LIBCURL/d"];
        sed_htslib_makefile(&out, &curl_patterns, "curl");
    } else if let Ok(inc) = env::var("DEP_CURL_INCLUDE").map(PathBuf::from) {
        cfg.include(inc);
    }

    let tool = cfg.get_compiler();
    let (cc_path, cflags_env) = (tool.path(), tool.cflags_env());
    let cc_cflags = cflags_env.to_string_lossy().replace("-O0", "");
    let host = env::var("HOST").unwrap_or_default();
    // autoreconf & ./configure (with no args) steps are necessary to include the htslib plugins (hfile_s3.o, hfile_s3_writer.o, etc...)
    if  // cleanup first
        // TODO: Have top level "cargo clean" do this instead of in here, see:
        // https://github.com/rust-lang/cargo/issues/572#issuecomment-632456478
        !Command::new("make")
        .current_dir(out.join("htslib"))
        .arg("clean")
        .status().unwrap().success()
    
        &&

        !Command::new("autoreconf")
        .current_dir(out.join("htslib"))
        .env("CFLAGS", &cc_cflags)
        .status().unwrap().success()
        
        &&
       
        !Command::new("./configure")
        .current_dir(out.join("htslib"))
        .env("CFLAGS", &cc_cflags)
        .arg(format!("--host={}", &host))
        .status().unwrap().success()
    {
        panic!("could not configure htslib nor any of its plugins")
    }

    if !Command::new("make")
        .current_dir(out.join("htslib"))
        .arg(format!("CC={}", cc_path.display()))
        .arg(format!("CFLAGS={}", &cc_cflags))
        .arg("lib-static")
        .status()
        .unwrap()
        .success()
    {
        panic!("failed to build htslib");
    }

    cfg.file("wrapper.c").compile("wrapper");

    // If bindgen is enabled, use it
    #[cfg(feature = "bindgen")]
    {
        bindgen::Builder::default()
            .header("wrapper.h")
            .generate_comments(false)
            .blacklist_function("strtold")
            .blacklist_type("max_align_t")
            .generate()
            .expect("Unable to generate bindings.")
            .write_to_file(out.join("bindings.rs"))
            .expect("Could not write bindings.");
    }

    // If no bindgen, use pre-built bindings
    #[cfg(all(not(feature = "bindgen"), target_os="macos"))]
    {
        fs::copy("osx_prebuilt_bindings.rs", out.join("bindings.rs"))
            .expect("couldn't copy prebuilt bindings");
        println!("cargo:rerun-if-changed=osx_prebuilt_bindings.rs");
    }

    #[cfg(all(not(feature = "bindgen"), target_os="linux"))]
    {
        fs::copy("linux_prebuilt_bindings.rs", out.join("bindings.rs"))
            .expect("couldn't copy prebuilt bindings");
        println!("cargo:rerun-if-changed=linux_prebuilt_bindings.rs");
    }

    let include = out.join("include");
    fs::create_dir_all(&include).unwrap();
    if include.join("htslib").exists() {
        fs::remove_dir_all(include.join("htslib")).expect("remove exist include dir");
    }
    copy_directory(out.join("htslib").join("htslib"), &include).unwrap();
    fs::copy(out.join("htslib").join("libhts.a"), out.join("libhts.a")).unwrap();

    println!("cargo:root={}", out.display());
    println!("cargo:include={}", include.display());
    println!("cargo:libdir={}", out.display());
    println!("cargo:rustc-link-lib=static=hts"); // XXX: only for static, adapt for dynamic?
    println!("cargo:rerun-if-changed=wrapper.c");
    println!("cargo:rerun-if-changed=wrapper.h");
    println!("cargo:rerun-if-changed=htslib/Makefile");
    let globs = std::iter::empty()
        .chain(glob("htslib/*.[ch]").unwrap())
        .chain(glob("htslib/cram/*.[ch]").unwrap())
        .chain(glob("htslib/htslib/*.h").unwrap())
        .chain(glob("htslib/os/*.[ch]").unwrap())
        .filter_map(Result::ok);
    for htsfile in globs {
        println!("cargo:rerun-if-changed={}", htsfile.display());
    }
}
