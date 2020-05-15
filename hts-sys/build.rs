// Copyright 2014 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use fs_utils::copy::copy_directory;
use glob::glob;

use std::env;
use std::fs;
use std::path::PathBuf;
use std::process::Command;

fn sed_htslib_makefile(out: &PathBuf, patterns: &[&str], feature: &str) {
    println!("SUBSTITUTING THINGS!");
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

    // if let Ok(z_inc) = env::var("DEP_Z_INCLUDE") {
    //     if !want_static {
    //         cfg.include(z_inc);
    //     }
    // }

    if !out.join("htslib").exists() {
        copy_directory("htslib", &out).unwrap();
    }
    
    // Nice flags to have
    let cflags_patterns = vec!["s/-g -Wall -O2 -fvisibility=hidden/-g -Wall -O2 -fvisibility=hidden -fPIC -pedantic -std=c90 -D_XOPEN_SOURCE=600/"];
    sed_htslib_makefile(&out, &cflags_patterns, "cflags");

    // let use_bzip2 = env::var("CARGO_FEATURE_BZIP2").is_ok();
    // if !use_bzip2 {
    //     let bzip2_patterns = vec!["s/ -lbz2//", "/#define HAVE_LIBBZ2/d"];
    //     sed_htslib_makefile(&out, &bzip2_patterns, "bzip2");
    // } else if let Ok(inc) = env::var("DEP_BZIP2_ROOT")
    //     .map(PathBuf::from)
    //     .map(|path| path.join("include"))
    // {
    //     cfg.include(inc);
    // }

    // let use_lzma = env::var("CARGO_FEATURE_LZMA").is_ok();
    // if !use_lzma && want_static {
    //     let lzma_patterns = vec!["s/ -llzma//", "/#define HAVE_LIBLZMA/d"];
    //     sed_htslib_makefile(&out, &lzma_patterns, "lzma");
    // } else if let Ok(inc) = env::var("DEP_LZMA_INCLUDE").map(PathBuf::from) {
    //     cfg.include(inc);
    // }

    // let use_curl = env::var("CARGO_FEATURE_CURL").is_ok();
    // if !use_curl {
    //     let curl_patterns = vec!["s/ -lcurl//", "/#define HAVE_LIBCURL/d"];
    //     sed_htslib_makefile(&out, &curl_patterns, "curl");
    // } else if let Ok(inc) = env::var("DEP_CURL_INCLUDE").map(PathBuf::from) {
    //     cfg.include(inc);
    // }

    let tool = cfg.get_compiler();
    let (cc_path, cflags_env) = (tool.path(), tool.cflags_env());
    let _cc_cflags = cflags_env.to_string_lossy().replace("-O0", "");
    let extra_cc_cflags = "-g -Wall -O2 -fvisibility=hidden -fPIC".to_string();
    // Some other flags which can be critical for cross-compiling to targets like MUSL
    //let cppflags = env::var("CPPFLAGS").unwrap_or_default();
    let _ldflags= env::var("LDFLAGS").unwrap_or_default();
    let host = env::var("HOST").unwrap_or_default();
    // Those two steps are necessary to include the htslib plugins in the resulting libhts.a (hfile_s3.o, hfile_s3_writer.o, etc...)
    if !Command::new("autoreconf")
        .current_dir(out.join("htslib"))
        .env("CFLAGS", &extra_cc_cflags)
        .status().unwrap().success()
        
        &&
       
        !Command::new("./configure")
        .current_dir(out.join("htslib"))
        .env("CFLAGS", &extra_cc_cflags)
        .arg(format!("--host={}", &host))
        .status().unwrap().success()
    {
        panic!("could not configure htslib nor any of its plugins")
    }

    if !Command::new("make")
        .current_dir(out.join("htslib"))
        .arg(format!("CC={}", cc_path.display()))
        .env("CFLAGS", &extra_cc_cflags)
        //.arg(format!("CPPFLAGS=\"{}\"", &cppflags))
        //.arg(format!("LDFLAGS=\"{}\"", &ldflags))
        .arg("lib-static")
        .arg("-j40")
        //.arg("-B")
        .status()
        .unwrap()
        .success()
    {
        panic!("failed to build htslib");
    }

    cfg.file("wrapper.c").compile("wrapper");

    bindgen::Builder::default()
        .header("wrapper.h")
        .generate_comments(false)
        .blacklist_function("strtold")
        .blacklist_type("max_align_t")
        .generate()
        .expect("Unable to generate bindings.")
        .write_to_file(out.join("bindings.rs"))
        .expect("Could not write bindings.");

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
    println!("cargo:rustc-link-lib=static=hts"); // XXX: only for static
    println!("cargo:rerun-if-changed=wrapper.c");
    println!("cargo:rerun-if-changed=wrapper.h");
    for htsfile in glob("htslib/**/*").unwrap() {
        let htsfile = htsfile.as_ref().unwrap().to_str().unwrap();
        println!("cargo:rerun-if-changed={}", htsfile);
    }
}
