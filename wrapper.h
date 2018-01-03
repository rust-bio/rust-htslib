#include "htslib/htslib/hts.h"
#include "htslib/htslib/vcf.h"
#include "htslib/htslib/sam.h"
#include "htslib/htslib/bgzf.h"
#include "htslib/htslib/vcfutils.h"


/* Replacements for macros, as we cannot automatically generate bindings for
   them.
 */
int64_t bgzf_tell_func(BGZF *fp) {
    return (((fp)->block_address << 16) | ((fp)->block_offset & 0xFFFF));
}
