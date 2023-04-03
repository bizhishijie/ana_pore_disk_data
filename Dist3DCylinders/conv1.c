// Compile: 
// >> mex -O -R2017b conv1.c
// >> mex -O -R2018a conv1.c

#include "mex.h"
#include "matrix.h"

#define FULL_SHAPE						0
#define SAME_SHAPE						1
#define VALID_SHAPE						2

#define CEIL(x)		ceil(x)
#define FLOOR(x)	floor(x)
#define ROUND(x)	round(x)

void conv1(const double *a, int na, const double *b, int nb, 
		 double *c, int *pint_nc, int shape)
/* Naive implementation, efficient when A or B is small */
{
	int nc, i, j, j1, j2, jmax, nanb;
	double s;
    double *c0;
	const double *pa, *pb;

    //c0 = c;
	switch (shape)
	{
	case VALID_SHAPE:
		j1 = nb - 1;
		j2 = na;
		break;
	case FULL_SHAPE:
		j1 = 0;
		j2 = na + nb - 1;
		break;
	case SAME_SHAPE:
		j1 = nb / 2;
		j2 = j1 + na;
		break;
	}

	if (pint_nc)
	{
		/* Compute length of the output of requested */
		nc = j2 - j1;
		if (nc < 0) nc = 0;
		*pint_nc = nc;
	}

	na--; nb--;
	jmax = min(na, nb);
	for (j = j1; j < jmax; j++)
	{
		i = j + 1;
		s = 0;
		pa = a;
		pb = b + j;
		while (i--) s += *(pa++) * *(pb--);
		*(c++) = s;
	}
	if (nb <= na)
		for (; j < nb; j++)
		{
			i = j + 1;
			s = 0;
			pa = a;
			pb = b + j;
			while (i--) s += *(pa++) * *(pb--);
			*(c++) = s;
		}
	else
	{
		jmax = min(nb, j2);
		na++;
		for (; j < jmax; j++)
		{
			i = na;
			s = 0;
			pa = a;
			pb = b + j;
			while (i--) s += *(pa++) * *(pb--);
			*(c++) = s;
		}
		na--;
	}
	b += nb;
	a -= nb;
	nb++;
	for (; j < na; j++)
	{
		i = nb;
		s = 0;
		pa = a + j;
		pb = b;
		while (i--) s += *(pa++) * *(pb--);
		*(c++) = s;
	}
	nanb = na + nb;
	jmax = min(nanb, j2);
	for (; j < jmax; j++)
	{
		i = nanb - j;
		s = 0;
		pa = a + j;
		pb = b;
		while (i--) s += *(pa++) * *(pb--);
		*(c++) = s;
	}
	for (; j < j2; j++) *(c++) = 0;

} /* conv1 */

#define A       prhs[0]
#define B       prhs[1]
#define SHAPE   prhs[2]

#define C       plhs[0]

void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[]) {
    
    double *a, *b, *c;
    int na, nb, nc;
    int shape;
    
#if MX_HAS_INTERLEAVED_COMPLEX
    a = mxGetDoubles(A);
    b = mxGetDoubles(B);
#else
    a = mxGetPr(A);
    b = mxGetPr(B);
#endif
    na = (int)mxGetN(A);
    nb = (int)mxGetN(B);
    c = mxMalloc(sizeof(double)*(na+nb-1));
    shape = FULL_SHAPE;
    if (nrhs >= 3)
        shape = (int)mxGetScalar(SHAPE);
    shape = max(shape,0);
    shape = min(shape,2);
    conv1(a, na, b, nb, c, &nc, shape);
    C = mxCreateDoubleMatrix(1,nc,mxREAL);
#if MX_HAS_INTERLEAVED_COMPLEX 
    memcpy(mxGetDoubles(C),c,sizeof(double)*nc);
#else
    memcpy(mxGetPr(C),c,sizeof(double)*nc);
#endif
    mxFree(c);
}
