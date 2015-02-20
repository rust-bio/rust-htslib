#![feature(process)]
#![feature(env)]

use std::process::Command;
use std::env;

fn main() {
    let out_dir = env::var("OUT_DIR").unwrap();

    Command::new("make").current_dir("htslib")
                        .arg("CFLAGS=-g -Wall -O2 -fPIC")
                        .arg("lib-static")
                        .status().ok().expect("Failed to build htslib");

    println!("cargo:rustc-flags=-L htslib -l static=hts -l z");
}
