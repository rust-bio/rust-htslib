

use std::process::Command;

fn main() {

    Command::new("make").current_dir("htslib")
                        .arg("CFLAGS=-g -Wall -O2 -fPIC")
                        .arg("lib-static")
                        .status().ok().expect("Failed to build htslib");

    println!("cargo:rustc-flags=-L htslib -l static=hts -l z");
}
