/*
 * File:    bch3.c
 * Title:   Encoder/decoder for binary BCH codes in C (Version 3.1)
 * Author:  Robert Morelos-Zaragoza
 * Date:    August 1994
 * Revised: June 13, 1997
 *
 * ===============  Encoder/Decoder for binary BCH codes in C =================
 *
 * Version 1:   Original program. The user provides the generator polynomial
 *              of the code (cumbersome!).
 * Version 2:   Computes the generator polynomial of the code.
 * Version 3:   No need to input the coefficients of a primitive polynomial of
 *              degree m, used to construct the Galois Field GF(2**m). The
 *              program now works for any binary BCH code of length such that:
 *              2**(m-1) - 1 < length <= 2**m - 1
 *
 * Note:        You may have to change the size of the arrays to make it work.
 *
 * The encoding and decoding methods used in this program are based on the
 * book "Error Control Coding: Fundamentals and Applications", by Lin and
 * Costello, Prentice Hall, 1983.
 *
 * Thanks to Patrick Boyle (pboyle@era.com) for his observation that 'bch2.c'
 * did not work for lengths other than 2**m-1 which led to this new version.
 * Portions of this program are from 'rs.c', a Reed-Solomon encoder/decoder
 * in C, written by Simon Rockliff (simon@augean.ua.oz.au) on 21/9/89. The
 * previous version of the BCH encoder/decoder in C, 'bch2.c', was written by
 * Robert Morelos-Zaragoza (robert@spectra.eng.hawaii.edu) on 5/19/92.
 *
 * NOTE:    
 *          The author is not responsible for any malfunctioning of
 *          this program, nor for any damage caused by it. Please include the
 *          original program along with these comments in any redistribution.
 *
 *  For more information, suggestions, or other ideas on implementing error
 *  correcting codes, please contact me at:
 *
 *                           Robert Morelos-Zaragoza
 *                           5120 Woodway, Suite 7036
 *                           Houston, Texas 77056
 *
 *                    email: r.morelos-zaragoza@ieee.org
 *
 * COPYRIGHT NOTICE: This computer program is free for non-commercial purposes.
 * You may implement this program for any non-commercial application. You may 
 * also implement this program for commercial purposes, provided that you
 * obtain my written permission. Any modification of this program is covered
 * by this copyright.
 *
 * == Copyright (c) 1994-7,  Robert Morelos-Zaragoza. All rights reserved.  ==
 *
 * m = order of the Galois field GF(2**m) 
 * n = 2**m - 1 = size of the multiplicative group of GF(2**m)
 * length = length of the BCH code
 * t = error correcting capability (max. no. of errors the code corrects)
 * d = 2*t + 1 = designed min. distance = no. of consecutive roots of g(x) + 1
 * k = n - deg(g(x)) = dimension (no. of information bits/codeword) of the code
 * p[] = coefficients of a primitive polynomial used to generate GF(2**m)
 * g[] = coefficients of the generator polynomial, g(x)
 * alpha_to [] = log table of GF(2**m) 
 * index_of[] = antilog table of GF(2**m)
 * data[] = information bits = coefficients of data polynomial, i(x)
 * bb[] = coefficients of redundancy polynomial x^(length-k) i(x) modulo g(x)
 * numerr = number of errors 
 * errpos[] = error positions 
 * recd[] = coefficients of the received polynomial 
 * decerror = number of decoding errors (in _message_ positions) 
 *
 */

#include <math.h>
#include <stdio.h>
#define noop (void)0
#define SIZE 131
#define SIZE2 22

#include "hal.h"
#include "simpleserial.h"
#include <stdint.h>
#include <stdlib.h>
#define M 7
// false or true
#define THROWS 0 

#define LEN 127
#define EFFORT 10
uint8_t             m, n, length, k, t, d;
uint8_t             p[14];
uint8_t             alpha_to[SIZE];
int8_t index_of[SIZE];
uint8_t g[SIZE];
uint8_t             recd[SIZE], data[SIZE], bb[SIZE];
int             seed;
uint8_t             numerr, errpos[SIZE], decerror = 0;
uint8_t mem[20];

uint32_t var_test = 0;
uint32_t temp_test = 0;




void 
read_p()
/*
 *	Read m, the degree of a primitive polynomial p(x) used to compute the
 *	Galois field GF(2**m). Get precomputed coefficients p[] of p(x). Read
 *	the code length.
 */
{
	uint8_t			i, ninf;

    m=M;
     if ( !(m>1) || !(m<21) ) exit(0);
	for (i=1; i<m; i++)
		p[i] = 0;
	p[0] = p[m] = 1;
	if (m == 2)			p[1] = 1;
	else if (m == 3)	p[1] = 1;
	else if (m == 4)	p[1] = 1;
	else if (m == 5)	p[2] = 1;
	else if (m == 6)	p[1] = 1;
	else if (m == 7)	p[1] = 1;
	else if (m == 8)	p[4] = p[5] = p[6] = 1;
	else if (m == 9)	p[4] = 1;
	else if (m == 10)	p[3] = 1;
    n = 1;
	for (i = 0; i < m; i++) {
        n *= 2;
        }
	n = n - 1;
	ninf = (n + 1) / 2 - 1;
	length = LEN;
    // pt[13]=n;
    // pt[14]=length;
    // pt[15]=ninf;
	if( !((length <= n)&&(length>ninf)) ) exit(0);
}


void 
generate_gf()
/*
 * Generate field GF(2**m) from the irreducible polynomial p(X) with
 * coefficients in p[0]..p[m].
 *
 * Lookup tables:
 *   index->polynomial form: alpha_to[] contains j=alpha^i;
 *   polynomial form -> index form:	index_of[j=alpha^i] = i
 *
 * alpha=2 is the primitive element of GF(2**m) 
 */
{
	register uint8_t    i, mask;

	mask = 1;
	alpha_to[m] = 0;
	for (i = 0; i < m; i++) {
		alpha_to[i] = mask;
		index_of[alpha_to[i]] = i;
		if (p[i] != 0)
			alpha_to[m] ^= mask;
		mask <<= 1;
	}
	index_of[alpha_to[m]] = m;
	mask >>= 1;
	for (i = m + 1; i < n; i++) {
		if (alpha_to[i - 1] >= mask)
		  alpha_to[i] = alpha_to[m] ^ ((alpha_to[i - 1] ^ mask) << 1);
		else
		  alpha_to[i] = alpha_to[i - 1] << 1;
		index_of[alpha_to[i]] = i;
	}
	index_of[0] = -1;
}


void 
gen_poly()
/*
 * Compute the generator polynomial of a binary BCH code. Fist generate the
 * cycle sets modulo 2**m - 1, cycle[][] =  (i, 2*i, 4*i, ..., 2^l*i). Then
 * determine those cycle sets that contain integers in the set of (d-1)
 * consecutive integers {1..(d-1)}. The generator polynomial is calculated
 * as the product of linear factors of the form (x+alpha^i), for every i in
 * the above cycle sets.
 */
{
	register uint8_t	ii, jj, ll, kaux;
	register uint8_t	test, aux, nocycles, root, noterms, rdncy;
	uint8_t             cycle[SIZE][11], size[SIZE], min[SIZE], zeros[SIZE];

	/* Generate cycle sets modulo n, n = 2**m - 1 */
	cycle[0][0] = 0;
	size[0] = 1;
	cycle[1][0] = 1;
	size[1] = 1;
	jj = 1;			/* cycle set index */
	do {
		/* Generate the jj-th cycle set */
		ii = 0;
		do {
			ii++;
			cycle[jj][ii] = (cycle[jj][ii - 1] * 2) % n;
			size[jj]++;
			aux = (cycle[jj][ii] * 2) % n;
		} while (aux != cycle[jj][0]);
		/* Next cycle set representative */
		ll = 0;
		do {
			ll++;
			test = 0;
			for (ii = 1; ((ii <= jj) && (!test)); ii++)	
			/* Examine previous cycle sets */
			  for (kaux = 0; ((kaux < size[ii]) && (!test)); kaux++)
			     if (ll == cycle[ii][kaux])
			        test = 1;
		} while ((test) && (ll < (n - 1)));
		if (!(test)) {
			jj++;	/* next cycle set index */
			cycle[jj][0] = ll;
			size[jj] = 1;
		}
	} while (ll < (n - 1));
	nocycles = jj;		/* number of cycle sets modulo n */

    t = EFFORT;
	d = 2 * t + 1;

	/* Search for roots 1, 2, ..., d-1 in cycle sets */
	kaux = 0;
	rdncy = 0;
	for (ii = 1; ii <= nocycles; ii++) {
		min[kaux] = 0;
		test = 0;
		for (jj = 0; ((jj < size[ii]) && (!test)); jj++)
			for (root = 1; ((root < d) && (!test)); root++)
				if (root == cycle[ii][jj])  {
					test = 1;
					min[kaux] = ii;
				}
		if (min[kaux]) {
			rdncy += size[min[kaux]];
			kaux++;
		}
	}
	noterms = kaux;
	kaux = 1;
	for (ii = 0; ii < noterms; ii++)
		for (jj = 0; jj < size[min[ii]]; jj++) {
			zeros[kaux] = cycle[min[ii]][jj];
			kaux++;
		}

	k = length - rdncy;

    if (k<0)
      {
         exit(0);
      }

	//("This is a (%d, %d, %d) binary BCH code\n", length, k, d);

	/* Compute the generator polynomial */
	g[0] = alpha_to[zeros[1]];
	g[1] = 1;		/* g(x) = (X + zeros[1]) initially */
	for (ii = 2; ii <= rdncy; ii++) {
	  g[ii] = 1;
	  for (jj = ii - 1; jj > 0; jj--)
	    if (g[jj] != 0)
	      g[jj] = g[jj - 1] ^ alpha_to[(index_of[g[jj]] + zeros[ii]) % n];
	    else
	      g[jj] = g[jj - 1];
	  g[0] = alpha_to[(index_of[g[0]] + zeros[ii]) % n];
	}
}


void 
encode_bch()
/*
 * Compute redundacy bb[], the coefficients of b(x). The redundancy
 * polynomial b(x) is the remainder after dividing x^(length-k)*data(x)
 * by the generator polynomial g(x).
 */
{
	register int    i, j;
	register int    feedback;

	for (i = 0; i < length - k; i++)
		bb[i] = 0;
	for (i = k - 1; i >= 0; i--) {
		feedback = data[i] ^ bb[length - k - 1];
		if (feedback != 0) {
			for (j = length - k - 1; j > 0; j--)
				if (g[j] != 0)
					bb[j] = bb[j - 1] ^ feedback;
				else
					bb[j] = bb[j - 1];
			bb[0] = g[0] && feedback;
		} else {
			for (j = length - k - 1; j > 0; j--)
				bb[j] = bb[j - 1];
			bb[0] = 0;
		}
	}
}


int8_t 
pompom()
/*
 * Simon Rockliff's implementation of Berlekamp's algorithm.
 *
 * Assume we have received bits in recd[i], i=0..(n-1).
 *
 * Compute the 2*t syndromes by substituting alpha^i into rec(X) and
 * evaluating, storing the syndromes in s[i], i=1..2t (leave s[0] zero) .
 * Then we use the Berlekamp algorithm to find the error location polynomial
 * elp[i].
 *
 * If the degree of the elp is >t, then we cannot correct all the errors, and
 * we have detected an uncorrectable error pattern. We output the information
 * bits uncorrected.
 *
 * If the degree of elp is <=t, we substitute alpha^i , i=1..n into the elp
 * to get the roots, hence the inverse roots, the error location numbers.
 * This step is usually called "Chien's search".
 *
 * If the number of errors located is not equal the degree of the elp, then
 * the decoder assumes that there are more than t errors and cannot correct
 * them, only detect them. We output the information bits uncorrected.
 */
{
	register uint8_t    i, j, u, q, t2, count = 0, syn_error = 0;
	int8_t             elp[SIZE2][SIZE2], d[SIZE], l[SIZE], u_lu[SIZE], s[SIZE];
	uint8_t             loc[SIZE], err[SIZE], reg[SIZE];
	t2 = 2 * t;
    /*
    var_test = var_test + temp_test;
    var_test = var_test + temp_test;
    var_test = var_test + temp_test;
    var_test = var_test + temp_test;
    var_test = var_test + temp_test;
    var_test = var_test + temp_test;
    var_test = var_test + temp_test;*/
	/* first form the syndromes */
	for (i = 1; i <= t2; i++) {
        //var_test = var_test + temp_test;
        //var_test = var_test + temp_test;
		s[i] = 0;
		for (j = 0; j < length; j++)
			s[i] ^= recd[j] * alpha_to[(i * j) % n];
		syn_error |= (s[i] != 0); /* set error flag if non-zero syndrome */
/*
 * Note:    If the code is used only for ERROR DETECTION, then
 *          exit program here indicating the presence of errors.
 */
		/* convert syndrome from polynomial form to index form  */
		s[i] = index_of[s[i]];
	}
    
    var_test = var_test + temp_test;
    var_test = var_test + temp_test;
    var_test = var_test + temp_test;
    var_test = var_test + temp_test;
    var_test = var_test + temp_test;
    var_test = var_test + temp_test;
    var_test = var_test + temp_test;
    
	if (syn_error) {	/* if there are errors, try to correct them */
		/*
		 * Compute the error location polynomial via the Berlekamp
		 * iterative algorithm. Following the terminology of Lin and
		 * Costello's book :   d[u] is the 'mu'th discrepancy, where
		 * u='mu'+1 and 'mu' (the Greek letter!) is the step number
		 * ranging from -1 to 2*t (see L&C),  l[u] is the degree of
		 * the elp at that step, and u_l[u] is the difference between
		 * the step number and the degree of the elp. 
		 */
		/* initialise table entries */
		d[0] = 0;			/* index form */
		d[1] = s[1];		/* index form */
		elp[0][0] = 0;		/* index form */
		elp[1][0] = 1;		/* polynomial form */
		for (i = 1; i < t2; i++) {
			elp[0][i] = -1;	/* index form */
			elp[1][i] = 0;	/* polynomial form */
		}
		l[0] = 0;
		l[1] = 0;
		u_lu[0] = -1;
		u_lu[1] = 0;
		u = 0;
 
		do {
			u++;
			if (d[u] == -1) {
				l[u + 1] = l[u];
				for (i = 0; i <= l[u]; i++) {
					elp[u + 1][i] = elp[u][i];
					elp[u][i] = index_of[elp[u][i]];
				}
			} else
				/*
				 * search for words with greatest u_lu[q] for
				 * which d[q]!=0 
				 */
			{
				q = u - 1;
				while ((d[q] == -1) && (q > 0))
					q--;
				/* have found first non-zero d[q]  */
				if (q > 0) {
				  j = q;
				  do {
				    j--;
				    if ((d[j] != -1) && (u_lu[q] < u_lu[j]))
				      q = j;
				  } while (j > 0);
				}
 
				/*
				 * have now found q such that d[u]!=0 and
				 * u_lu[q] is maximum 
				 */
				/* store degree of new elp polynomial */
				if (l[u] > l[q] + u - q)
					l[u + 1] = l[u];
				else
					l[u + 1] = l[q] + u - q;
 
				/* form new elp(x) */
				for (i = 0; i < t2; i++)
					elp[u + 1][i] = 0;
				for (i = 0; i <= l[q]; i++)
					if (elp[q][i] != -1)
						elp[u + 1][i + u - q] = 
                                   alpha_to[(d[u] + n - d[q] + elp[q][i]) % n];
				for (i = 0; i <= l[u]; i++) {
					elp[u + 1][i] ^= elp[u][i];
					elp[u][i] = index_of[elp[u][i]];
				}
			}
			u_lu[u + 1] = u - l[u + 1];
 
			/* form (u+1)th discrepancy */
			if (u < t2) {	
			/* no discrepancy computed on last iteration */
			  if (s[u + 1] != -1)
			    d[u + 1] = alpha_to[s[u + 1]];
			  else{
			    d[u + 1] = 0;
              }
              for (i = 1; i <= l[u + 1]; i++)
                if ((s[u + 1 - i] != -1) && (elp[u + 1][i] != 0))
                  d[u + 1] ^= alpha_to[(s[u + 1 - i] 
                              + index_of[elp[u + 1][i]]) % n];
			  /* put d[u+1] into index form */
			  d[u + 1] = index_of[d[u + 1]];	
			}
		} while ((u < t2) && (l[u + 1] <= t));
 
		u++;
		if (l[u] <= t) {/* Can correct errors */
			/* put elp into index form */
			for (i = 0; i <= l[u]; i++)
				elp[u][i] = index_of[elp[u][i]];

			/* Chien search: find roots of the error location polynomial */
			for (i = 1; i <= l[u]; i++)
				reg[i] = elp[u][i];
			count = 0;
			for (i = 1; i <= n; i++) {
				q = 1;
				for (j = 1; j <= l[u]; j++)
					if (reg[j] != -1) {
						reg[j] = (reg[j] + j) % n;
						q ^= alpha_to[reg[j]];
					}
				if (!q) {	/* store root and error
						 * location number indices */
					loc[count] = n - i;
					count++;
				}
			}
			if (count == l[u] || !THROWS)	{ // this if is harmful if you need to get a result every time
			/* no. roots = degree of elp hence <= t errors */
				for (i = 0; i < l[u]; i++)
					recd[loc[i]] ^= 1;
                //buffer[0] = count;
                return l[u];
            }
			else{
            }/* elp has degree >t hence cannot solve */
				
		}
	}
    else return 0;
    return -1;
}


//uint8_t info[128];

//uint8_t keep[64];
uint8_t diff = 0;
uint8_t helper = 0;

uint8_t
set_diff(uint8_t* m, uint8_t len) {
    diff = m[0];
    simpleserial_put('r', 1, m);
    return 0x00;
}
uint32_t base32[4];

void
preset(uint8_t* pt, uint8_t* buffer) {
    
    
  //   for (i = 0; i < k; i++) {
		// keep[i] = data[i] = random() & 1; //( random() & 65536 ) >> 16;
  //   }
  

	// encode_bch();
 //    for (i = 0; i < length - k; i++)
	// 	recd[i] = bb[i];
	// for (i = 0; i < k; i++)
	// 	recd[i + length - k] = data[i];

    // should be in handle
    uint8_t i;
	for(i = 0; i < 4; i++) {
		base32[i] = 0;
        // recd[i] ^= info[i];
	}

    
    
}



uint8_t set_var(uint8_t* pt, uint8_t len)
{
    var_test = 0;
	simpleserial_put('r', 4, pt);
    temp_test = ((uint32_t*)pt)[0];
	return 0x00;
}




void
handle(uint8_t* pt, uint8_t* buffer)
{
	uint16_t i;
    
    //for(i = 0; i < 4; i++) {
	//	base32[i] ^= pt32[i];
        // recd[i] ^= info[i];
    //}
	/* Randomly generate DATA */
    
	// for(i = 0; i < length; i++) {
	// 	recd[i] ^= mem[i];
 //        recd[i] ^= info[i];
	// }
	// to combat any evil compiler optimization
    // helper = 0;
    // helper = diff;
    //pt[15]=diff&1;
    //uint32_t* mem32 = mem;
    //uint32_t* pt32 = (uint32_t*)pt;
    /*
    var_test = var_test + temp_test;
    var_test = var_test + temp_test;
    var_test = var_test + temp_test;
    var_test = var_test + temp_test;
    var_test = var_test + temp_test;
    var_test = var_test + temp_test;
    var_test = var_test + temp_test;*/
    for(i = 0; i < 4; i++) {
		base32[i] = (((uint32_t*)mem)[i] ^ ((uint32_t*)pt)[i+4]);
        // recd[i] ^= info[i];
	}
    /*var_test = var_test + temp_test;
    var_test = var_test + temp_test;
    var_test = var_test + temp_test;
    var_test = var_test + temp_test;
    var_test = var_test + temp_test;
    var_test = var_test + temp_test;
    var_test = var_test + temp_test;*/
    for(i = 0; i < 4; i++) {
		base32[i] ^= ((uint32_t*)pt)[i];
        // recd[i] ^= info[i];
	}
    
    /*var_test = var_test + temp_test;
    var_test = var_test + temp_test;
    var_test = var_test + temp_test;
    var_test = var_test + temp_test;
    var_test = var_test + temp_test;
    var_test = var_test + temp_test;
    var_test = var_test + temp_test;*/
	/*
	 * recd[] are the coefficients of c(x) = x**(length-k)*data(x) + b(x)
	 */
    
    // numerr = 20;
	/*
	 * recd[] are the coefficients of r(x) = c(x) + e(x)
	 */
	// for (i = 0; i < numerr; i++)
	// 	errpos[i]=127-2*i;
	// if (numerr)
	// 	for (i = 0; i < numerr; i++)
	// 		recd[errpos[i]] ^= 1;
    
    //trigger_low();
    for (i = 0; i < length; i++){ // puts input in recd (to decode later)
		recd[i] = (((uint8_t*)base32)[i/8] >> (i%8)) & 1;
    }
    /*var_test = var_test + temp_test;
    var_test = var_test + temp_test;
    var_test = var_test + temp_test;
    var_test = var_test + temp_test;
    var_test = var_test + temp_test;
    var_test = var_test + temp_test;
    var_test = var_test + temp_test;*/
    //trigger_high();
	int8_t errs = pompom();             /* DECODE received codeword recv[] */

    
    /*var_test = var_test + temp_test;
    var_test = var_test + temp_test;
    var_test = var_test + temp_test;
    var_test = var_test + temp_test;
    var_test = var_test + temp_test;
    var_test = var_test + temp_test;
    var_test = var_test + temp_test;*/
    //pt[15]=errs;
    //pt[14]=buffer[0];
    //pt[13]=k;
}

void
finish(uint8_t* pt, uint8_t start, uint8_t* from) { // 
 //    for (int i = length - k; i < length; i++) { 
	// 	recd[i] ^= data[i-length+k];//info[i-length+k];
	// }
    int             i;
    uint8_t x = 0;
    int x_ = 0;
    for (i = start; i < length; i++) { 
		// if (data[i - length + k] != recd[i]) {
		// 	decerror++;
  //       }

        x >>= 1;
        x += from[i] << 7;

        // if(recd[i] != 0 && recd[i] != 1) {
        //     for(int i_=0; i_ < 16; i_++) pt[i_] = recd[i]+i_; // error
        // }
		
        if((++x_) % 8 == 0) {
            pt[x_/8-1] = x;
            x = 0;
        }
    }
    if(x_ % 8 != 0) {
        pt[x_/8] = x >> (8-(x_% 8));
    }
	/*
	 * print out original and decoded data
	 */

	/*
	 * DECODING ERRORS? we compare only the data portion
	 */
	// for (i = length - k; i < length; i++)
	// 	if (data[i - length + k] != recd[i])
	// 		decerror++;
    // if(pt[15] == -1) {
    //     pt[15] = 'L';
    //     return;
    // }
    
}



uint8_t get_mask(uint8_t* m, uint8_t len)
{
  aes_indep_mask(m, len);
  return 0x00;
}

uint8_t get_key(uint8_t* k, uint8_t len)
{

	int i=0;
    for (i = 0; i < 16; i++){
		mem[i] = k[i];//(k[i/8] >> (i%8)) & 1;
    }
	simpleserial_put('r', 16, k);
	return 0x00;
}

uint8_t get_pt(uint8_t* pt, uint8_t len)
{
    uint8_t buffer[16];
    preset(pt, buffer);
    noop;
    noop;
    noop;
    noop;
    noop;
    noop;
    noop;
    noop;
    noop;
    noop;
    noop;
    noop;
    noop;
    noop;
	trigger_high();
    
  #ifdef ADD_JITTER
  //for (volatile uint8_t k = 0; k < (*pt & 0x0F); k++);
  #endif

    handle(pt, buffer); /* encrypting the data block */

    trigger_low();
    uint8_t i;
 //    for(i = length-k; i < length; i++) {
 //        recd[i] ^= keep[i-length+k];
	// }
    finish(pt,length-k, &(recd[0]));


	simpleserial_put('r', 16, pt);
	return 0x00;
}

uint8_t set_seed(uint8_t* pt, uint8_t len) {
    seed =  * (int*) pt;
    srandom(seed);
    simpleserial_put('r', 4, pt);
    return 0x00;
}

uint8_t get_enc(uint8_t* pt, uint8_t len) {
    uint16_t i=0;
    for (i = 0; i < k; i++){ // puts input in recd (to decode later)
		data[i] = (pt[i/8] >> (i%8)) & 1;
    }
    for (i =0; i<length; i++) {
        pt[i]=0;
    }
    encode_bch();
    for (i = 0; i < length - k; i++)
		recd[i] = bb[i];
	for (i = 0; i < k; i++)
		recd[i + length - k] = data[i];
    finish(pt, 0, &(recd[0]));
    simpleserial_put('r', 16, pt);
    return 0x00;
}

uint8_t reset(uint8_t* x, uint8_t len)
{
    // Reset key here if needed
	return 0x00;
}
int main(void)
{


    platform_init();
    init_uart();
    trigger_setup();


    /* Uncomment this to get a HELLO message for debug */
	read_p();               /* Read m */
	generate_gf();          /* Construct the Galois Field GF(2**m) */
	gen_poly();             /* Compute the generator polynomial of BCH code */
	seed = 56214563;
	srandom(seed);
    // putch('h');
    // putch('e');
    // putch('l');
    // putch('l');
    // putch('o');
    // putch('\n');

	simpleserial_init();
    // #if SS_VER == SS_VER_2_1
    // simpleserial_addcmd(0x01, 16, aes);
    // #else
    simpleserial_addcmd('k', 16, get_key);
    simpleserial_addcmd('p', 32,  get_pt);
    simpleserial_addcmd('a', 4,  set_var);
    simpleserial_addcmd('e', 16,  get_enc);
    simpleserial_addcmd('t', 1,  set_diff);
    simpleserial_addcmd('s', 4,  set_seed);
    //simpleserial_addcmd('x',  0,   reset);
    simpleserial_addcmd_flags('m', 18, get_mask, CMD_FLAG_LEN);
    //simpleserial_addcmd('s', 2, enc_multi_setnum);
    //simpleserial_addcmd('f', 16, enc_multi_getpt);
    // #endif
    while(1)
        simpleserial_get();
}
