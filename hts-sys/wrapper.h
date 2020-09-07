#include "htslib/htslib/hts.h"
#include "htslib/htslib/vcf.h"
#include "htslib/htslib/sam.h"
#include "htslib/htslib/cram.h"
#include "htslib/htslib/bgzf.h"
#include "htslib/htslib/vcfutils.h"
#include "htslib/htslib/tbx.h"
#include "htslib/htslib/synced_bcf_reader.h"
#include "htslib/htslib/kbitset.h"
#include "htslib/htslib/faidx.h"
#include "htslib/htslib/thread_pool.h"

// The following functions have to be wrapped here because they are inline in htslib.

/**
 * <div rustbindgen replaces="kbs_init2"></div>
 */
kbitset_t *wrap_kbs_init2(size_t ni, int fill);

/**
 * <div rustbindgen replaces="kbs_init"></div>
 */
kbitset_t *wrap_kbs_init(size_t ni);


/**
 * <div rustbindgen replaces="kbs_insert"></div>
 */
void wrap_kbs_insert(kbitset_t *bs, int i);


/**
 * <div rustbindgen replaces="kbs_destroy"></div>
 */
void wrap_kbs_destroy(kbitset_t *bs);
