FROM rustembedded/cross:x86_64-unknown-linux-gnu

ENV LIBCLANG_PATH /usr/lib/x86_64-linux-gnu
ENV LLVM_CONFIG_PATH /usr/bin
RUN apt-get update && \
    apt-get install -y libcurl4-openssl-dev zlib1g-dev libbz2-dev liblzma-dev clang-8 && \
    ln -sf /usr/bin/llvm-config-8 /usr/bin/llvm-config && \
    ln -sf /usr/lib/x86_64-linux-gnu/libclang-8.so.1 /usr/lib/x86_64-linux-gnu/libclang.so.1
