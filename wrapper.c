#include "wrapper.h"


kbitset_t *wrap_kbs_init2(size_t ni, int fill)
{
	return kbs_init2(ni, fill);
}

kbitset_t *wrap_kbs_init(size_t ni)
{
	return wrap_kbs_init2(ni, 0);
}

void wrap_kbs_insert(kbitset_t *bs, int i)
{
  kbs_insert(bs, i);
}

void wrap_kbs_destroy(kbitset_t *bs)
{
  kbs_destroy(bs);
}
