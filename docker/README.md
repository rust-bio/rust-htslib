# cross rustembedded containers

Allows to compile (rust-)htslib in a variety of environments and architectures via [rustembedded cross](https://github.com/rust-embedded/cross).

## Quickstart

```shell
$ cd docker
$ docker build -t brainstorm/cross-x86_64-unknown-linux-musl:1.0.0 . -f Dockerfile.musl
$ docker build -t brainstorm/cross-x86_64-unknown-linux-gnu:1.0.0 . -f Dockerfile.gnu
```

Then, to build and test rust-htslib with the above containers, proceed as you would with `cargo`, using `cross` instead, i.e:

```shell
$ cross build --target x86_64-unknown-linux-musl 
```
