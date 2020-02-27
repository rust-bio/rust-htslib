FROM rustembedded/cross:x86_64-unknown-linux-gnu

ENV LIBCLANG_PATH /usr/lib/llvm-9/lib/clang/9.0.0/lib/linux/
ENV LLVM_CONFIG_PATH /usr/bin
RUN apt-get update && \
    apt-get install -y libcurl4-openssl-dev zlib1g-dev libbz2-dev liblzma-dev gcc-multilib libc6-dev llvm-dev libclang-dev clang && \
    ln -sf /usr/bin/llvm-config-9 /usr/bin/llvm-config
