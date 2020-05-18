#!/bin/sh -x

curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
. /root/.cargo/env
rustup target add x86_64-unknown-linux-musl

echo "Run source $HOME/.cargo/env to continue" 
