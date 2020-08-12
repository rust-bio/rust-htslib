FROM rustembedded/cross:x86_64-unknown-linux-gnu

ENV LIBCLANG_PATH /usr/lib/llvm-10/lib
ENV LLVM_CONFIG_PATH /usr/bin
RUN apt-get update
RUN apt-get install -y build-essential wget gnupg lsb-release software-properties-common apt-transport-https ca-certificates # Otherwise LLVM bump below fails
RUN bash -c "$(wget -O - https://apt.llvm.org/llvm.sh)"
