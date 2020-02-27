# cross rustembedded containers

Allows to compile (rust-)htslib in a variety of environments and architectures via [rustembedded cross](https://github.com/rust-embedded/cross).

## Quickstart

```shell
$ cd docker
$ docker build -t brainstorm/cross-x86_64-unknown-linux-musl:libcurl-openssl . -f dockerfile.musl
$ docker build -t brainstorm/cross-x86_64-unknown-linux-gnu:libcurl-openssl . -f dockerfile.gnu
```

Then to build and test rust-htslib with the above containers, proceed as you would with `cargo`, using `cross` instead, i.e:

```shell
$ cross build --target x86_64-unknown-linux-musl 
```
