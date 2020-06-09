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
use std::io::Write;
use std::path::PathBuf;

// these need to be kept in sync with the htslib Makefile
const FILES: &[&str] = &[
    "kfunc.c",
    "knetfile.c",
    "kstring.c",
    "bcf_sr_sort.c",
    "bgzf.c",
    "errmod.c",
    "faidx.c",
    "header.c",
    "hfile.c",
    "hfile_net.c",
    "hts.c",
    "hts_os.c",
    "md5.c",
    "multipart.c",
    "probaln.c",
    "realn.c",
    "regidx.c",
    "region.c",
    "sam.c",
    "synced_bcf_reader.c",
    "vcf_sweep.c",
    "tbx.c",
    "textutils.c",
    "thread_pool.c",
    "vcf.c",
    "vcfutils.c",
    "cram/cram_codecs.c",
    "cram/cram_decode.c",
    "cram/cram_encode.c",
    "cram/cram_external.c",
    "cram/cram_index.c",
    "cram/cram_io.c",
    "cram/cram_samtools.c",
    "cram/cram_stats.c",
    "cram/mFILE.c",
    "cram/open_trace_file.c",
    "cram/pooled_alloc.c",
    "cram/rANS_static.c",
    "cram/string_alloc.c",
];

fn main() {

    let out = PathBuf::from(env::var("OUT_DIR").unwrap());
    if !out.join("htslib").exists() {
        println!("copying...");
        copy_directory("htslib", &out).unwrap();
    }
    
    let mut cfg = cc::Build::new();
    
    // default files
    let out_htslib = out.join("htslib");
    let htslib = PathBuf::from("htslib");
    for f in FILES {
        let c_file = out_htslib.join(f);
        cfg.file(&c_file);
        println!("cargo:rerun-if-changed={}", htslib.join(f).display());
    }

    cfg.include(out.join("htslib"));

    let want_static = cfg!(feature = "static") || env::var("HTS_STATIC").is_ok();

    if want_static {
        cfg.warnings(true).static_flag(true).pic(true);
    } else {
        cfg.warnings(true).static_flag(false).pic(true); 
    }

    if let Ok(z_inc) = env::var("DEP_Z_INCLUDE") {
        cfg.include(z_inc);
    }


    // We build a config.h ourselves, rather than rely on Makefile or ./configure
    let mut config_lines = vec![
        "/* Default config.h generated by build.rs */",
        "#define HAVE_DRAND48 1",
    ];

    let use_bzip2 = env::var("CARGO_FEATURE_BZIP2").is_ok();
    if use_bzip2 {
         if let Ok(inc) = env::var("DEP_BZIP2_ROOT")
        .map(PathBuf::from)
        .map(|path| path.join("include"))
        {
            cfg.include(inc);
            config_lines.push("#define HAVE_LIBBZ2 1");
        }
    }

    let use_libdeflater = env::var("CARGO_FEATURE_LIBDEFLATER").is_ok();
    if use_libdeflater {
        if let Ok(inc) = env::var("DEP_LIBDEFLATE_INCLUDE").map(PathBuf::from) {
            cfg.include(inc);
            config_lines.push("#define HAVE_LIBDEFLATE 1");
        } else {
            panic!("no DEP_LIBDEFLATE_INCLUDE");
        }
    }

    let use_lzma = env::var("CARGO_FEATURE_LZMA").is_ok();
    if use_lzma {
        if let Ok(inc) = env::var("DEP_LZMA_INCLUDE").map(PathBuf::from) {
            cfg.include(inc);
            config_lines.push("#define HAVE_LIBBZ2 1");
            config_lines.push("#ifndef __APPLE__");
            config_lines.push("#define HAVE_LZMA_H 1");
            config_lines.push("#endif");
        }
    }

    let use_curl = env::var("CARGO_FEATURE_CURL").is_ok();
    if use_curl {
        if let Ok(inc) = env::var("DEP_CURL_INCLUDE").map(PathBuf::from) {
            cfg.include(inc);
            config_lines.push("#define HAVE_LIBCURL 1");
            cfg.file("htslib/hfile_libcurl.c");
            println!("cargo:rerun-if-changed=htslib/hfile_libcurl.c");

            #[cfg(target_os="macos")]
            {
                // Use builtin MacOS CommonCrypto HMAC
                config_lines.push("#define HAVE_COMMONCRYPTO 1");
            }

            #[cfg(not(target_os="macos"))]
            {
                if let Ok(inc) = env::var("DEP_OPENSSL_INCLUDE").map(PathBuf::from) {
                    // Must use hmac from libcrypto in openssl
                    cfg.include(inc);
                    config_lines.push("#define HAVE_HMAC 1");
                } else {
                    panic!("No OpenSSL dependency -- need OpenSSL includes");
                }
            }
        }
    }

    let use_gcs = env::var("CARGO_FEATURE_GCS").is_ok();
    if use_gcs {
        config_lines.push("#define ENABLE_GCS 1");
        cfg.file("htslib/hfile_gcs.c");
        println!("cargo:rerun-if-changed=htslib/hfile_gcs.c");
    }

    let use_s3 = env::var("CARGO_FEATURE_S3").is_ok();
    if use_s3 {
        config_lines.push("#define ENABLE_S3 1");
        cfg.file("htslib/hfile_s3.c");
        println!("cargo:rerun-if-changed=htslib/hfile_s3.c");
        cfg.file("htslib/hfile_s3_write.c");
        println!("cargo:rerun-if-changed=htslib/hfile_s3_write.c");
    }


    // write out config.h which controls the options htslib will use
    {
        let mut f = std::fs::File::create(out.join("htslib").join("config.h")).unwrap();
        for l in config_lines {
            writeln!(&mut f, "{}", l).unwrap();
        };
    }
    
    // write out version.h
    {
        let version = std::process::Command::new(out.join("htslib").join("version.sh"))
                             .output()
                             .expect("failed to execute process");
        let version_str = std::str::from_utf8(&version.stdout).unwrap().trim();

        let mut f = std::fs::File::create(out.join("htslib").join("version.h")).unwrap();
        writeln!(&mut f, "#define HTS_VERSION_TEXT \"{}\"", version_str).unwrap();
    }

    cfg.file("wrapper.c");
    cfg.compile("hts");

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

    println!("cargo:root={}", out.display());
    println!("cargo:include={}", include.display());
    println!("cargo:libdir={}", out.display());
    println!("cargo:rerun-if-changed=wrapper.c");
    println!("cargo:rerun-if-changed=wrapper.h");

    let globs = std::iter::empty()
        .chain(glob("htslib/*.[h]").unwrap())
        .chain(glob("htslib/cram/*.[h]").unwrap())
        .chain(glob("htslib/htslib/*.h").unwrap())
        .chain(glob("htslib/os/*.[h]").unwrap())
        .filter_map(Result::ok);
    for htsfile in globs {
        println!("cargo:rerun-if-changed={}", htsfile.display());
    }

    // Note: config.h is a function of the cargo features. Any feature change will 
    // cause build.rs to re-run, so don't re-run on that change.
    //println!("cargo:rerun-if-changed=htslib/config.h");
}