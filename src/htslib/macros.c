#include <sys/types.h>
#include <htslib/bgzf.h>

int64_t bgzf_tell_func(BGZF *fp) {
    return (((fp)->block_address << 16) | ((fp)->block_offset & 0xFFFF));
}
