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




#include "hal.h"
#include "simpleserial.h"
#include <stdint.h>
//#include <stdlib.h>


#define LEN 8
#define WIDTH uint32_t

uint8_t get_pt(uint8_t* pt, uint8_t len)
{
    // masking of 3rd order
    WIDTH* arr = (WIDTH*)pt;
    WIDTH x0 = arr[0]; //*();
    WIDTH x1 = arr[1]; //*();
    WIDTH x2 = arr[2]; //*();
    WIDTH y0 = arr[3]; //*();
    WIDTH y1 = arr[4]; //*();
    WIDTH y2 = arr[5]; //*();
    for(int i=0; i<len;i++) pt[i]=0;
    WIDTH z0=0, z1=0, z2=0;
    trigger_high();
    z0 = x0 ^ y0; // no leakage
    z1 = x1 ^ y1;
    z2 = x2 ^ y2;
    // uint8_t z2 = z0 ^ z1;
    arr[0] = z0;
    arr[1] = z1;
    arr[2] = z2;
    trigger_low();
	simpleserial_put('r', 12, pt);
	return 0x00;
}



int main(void)
{

    platform_init();
    init_uart();
    trigger_setup();

	//aes_indep_init();
	//aes_indep_key(tmp);

    /* Uncomment this to get a HELLO message for debug */

    // putch('h');
    // putch('e');
    // putch('l');
    // putch('l');
    // putch('o');
    // putch('\n');
	simpleserial_init();
    simpleserial_addcmd('p', 24, get_pt);
    while(1)
        simpleserial_get();
}
