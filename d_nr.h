/*
 * d_nr.h
 */

#ifndef _TOMO_D_NR_H_
#define _TOMO_D_NR_H_

void d_jacobi(double **a, int n, double d[], double **v, int *nrot);
void d_choldc(double **a, int n, double p[]);
void d_realft(double data[], unsigned long n, int isign);
void d_four1(double data[], unsigned long nn, int isign);

#endif /* _TOMO_D_NR_H_ */
