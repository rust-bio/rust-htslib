// Copyright 2014 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.


use std::process::Command;
use std::env;

fn main() {
    let root = env::var("CARGO_MANIFEST_DIR").unwrap();
    if Command::new("make").current_dir("htslib")
                           .arg("CFLAGS=-g -Wall -O2 -fPIC")
                           .arg("lib-static")
                           .status().unwrap()
                           .success() != true {
        panic!("failed to build htslib");
    }
    if Command::new("gcc").arg("-Ihtslib")
                          .arg("-c")
                          .arg("src/htslib/macros.c")
                          .status().unwrap().success() != true {
        panic!("error building macro wrapper");
    }
    if Command::new("ar").arg("r").arg("libmacros.a").arg("macros.o").status().unwrap().success() != true {
        panic!("error buiding macro wrapper");
    }

    println!("cargo:rustc-flags=-L {}/htslib -L . -l static=hts -l static=macros -l z", root);
}
