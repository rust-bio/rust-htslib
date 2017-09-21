// Copyright 2014 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

extern crate fs_utils;

use fs_utils::copy::copy_directory;

use std::process::Command;
use std::env;
use std::path::Path;

fn main() {
    let out = env::var("OUT_DIR").unwrap();
    let within_out = |p| format!("{}/{}", out, p);
    if !Path::new(&within_out("htslib")).exists() {
        copy_directory("htslib", &out).unwrap();
    }

    if Command::new("make").current_dir(within_out("htslib"))
                           .arg("CFLAGS=-g -Wall -O2 -fPIC")
                           .arg("lib-static")
                           .arg("-B")
                           .status().unwrap()
                           .success() != true {
        panic!("failed to build htslib");
    }
    if Command::new("gcc").arg("-I").arg(within_out("htslib"))
                          .arg("-c")
                          .arg("src/htslib/macros.c")
                          .arg("-o").arg(within_out("macros.o"))
                          .status().unwrap().success() != true {
        panic!("error building macro wrapper");
    }
    if Command::new("ar").arg("r")
                         .arg(within_out("libmacros.a"))
                         .arg(within_out("macros.o"))
                         .status().unwrap().success() != true {
        panic!("error buiding macro wrapper");
    }

    println!("cargo:rustc-flags=-L {out}/htslib -L {out} -l static=hts -l static=macros -l z",
             out=out);
}
