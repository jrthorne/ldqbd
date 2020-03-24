
/***********************
* File name = ldqbd.h  *
************************/
#ifndef _LDQBD_
#define _LDQBD_
/* Constants */
#define MAX_L 15 /* the maximum number of iterations allowable */
#define TOL 1.0e-8
/* types */
enum UD_type {R, G};
typedef enum UD_type UDtype;
/* Function templates */
/************************************************************************/
/* q_gen generates the elements of the matrix Q and returns a pointer
   to the first of the three elements. */
/* removed in matrix_ptr q_gen(int dim, int k);
matrix_ptr q_generate(int dim, int k); */
/************************************************************************/
/* UD_R calculates the UD pair for the calculation af R and returns
   a pointer to U(l,k) and pointer + 1 references D(l,k+2^l) */
matrix_ptr UD_R(matrix_ptr U, matrix_ptr D, int l, int k, int dim);
/***********************************************************************/
/* compR calculates the R matrix. Adapted from alg 3.4.4 */
matrix_ptr compR(float *diff, int l, matrix_ptr U, matrix_ptr D);
/************************************************************************/
/* calculate U for G seperately at l = 0 */
/* J1.10 matrix_ptr calcU_G(int l, int k, int dim);  */
/***********************************************************************/
/* Use equation 2.5.4 on page 15 (rearranged) to calculate R(k) from
   R(k+1). Iterate backeards */
matrix_ptr backR(matrix_ptr r, int k, int dim);
/************************************************************************/
/* UD_G calculates the UD pair for the calculation af G and returns
   a pointer to U(l,k-1+2^l) and pointer + 1 references D(l,k+2^(l+1)).
   For the method used in this function refer to  Bright P40
   */
matrix_ptr UD_G(matrix_ptr U, matrix_ptr D, int l, int k, int dim);
/***********************************************************************/
/* Use equation 2.7.7 & 2.7.8 on page 21 (rearranged) to calculate G(k) from
   G(k+1). Iterate backeards */
matrix_ptr backG(matrix_ptr g, int k, int dim);
/***********************************************************************/
/* compG calculates the G matrix. Adapted from alg 3.4.1 */
matrix_ptr compG(float *diff, int l, matrix_ptr U, matrix_ptr D);
/************************************************************************/
/* GfromR calculates the G matrix at k from Rk (see p21 eq 2.7.12 */
matrix_ptr GfromR(matrix_ptr Rk, int k, int dim);
/***********************************************************************/
/* RG process and compG/R are collectively the implementation of
algorithms 3.4.3 and 3.4.2 respectively */
matrix_ptr RGprocess(UDtype type, int k, int dim);
/*************************************************************************/
/* get_outer runs in conjunction with G/Rprocess to find the outer UD matricies
and store them in memory pointed to by the matrix _outer, where the index of the
matrix appertains to the level of l that the UD belong to. J1.11 */
void get_outer(int dim, UDtype type);
/*************************************************************************/
/* for use with get_outer, shuffleUD shuffles the UD pointers within the
   three sets of matricies */
void shuffleUD(int i);
/**************************************************************************/
/* chuck the correct matriceis in Uval and Dval for the calculation of the
   next UD */
void set_UDval(int i);

#endif
