#ifndef SMC_H
#define SMC_H

typedef double (*real_rv) (void);
typedef void* (*sampler) (double*);
typedef double (*pdf) (double*);
typedef double (*metric) (const void *, const void *);

#endif
