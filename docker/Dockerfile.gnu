FROM rustembedded/cross:x86_64-unknown-linux-gnu-0.1.16

#ENV LIBCLANG_PATH /usr/lib/llvm-9/lib/clang/9.0.0/lib/linux/
#ENV LLVM_CONFIG_PATH /usr/bin
RUN apt-get update && \
    apt-get install -y libcurl4-openssl-dev zlib1g-dev libbz2-dev libc6-dev-i386
#gcc-multilib libclang-common-9-dev llvm-9 clang-9 clang libclang-dev libc6-dev && \
#    ln -sf /usr/bin/llvm-config-9 /usr/bin/llvm-config
