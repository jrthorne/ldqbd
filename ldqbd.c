/*********************************************************************
* File name: ldqbd.c
* Date start/finish: 17-01-1997/dd-mm-yyyy
* By: Jason Thorne
* Name: Level Dependant Quasi-Birth Death routines
* Project: Matrix analytic methods in Quasi-stationary birth-death
*      processes Honours project 1996-1997
* Description: The routines for the calculations of R
*   and G in Les Bright's thesis are etc.
* Updates:
* JX.X:DD-MM-YYYY|XXX|Description
* J1.1:12-01-1997|JRT|Wrote compR, adapted from alg3.3.4 of L.Bright's
*                     thesis. Added functions calcU and calcD to calculate
*                     these values seperately.
* J1.2:13-01-1997|JRT|Changes at home.
* J1.3:14-01-1997|JRT|Changes at uni
* J1.4:15-01-1997|JRT|Changes at Uni
* J1.5:16-01-1997|JRT|Changes at uni. Code to chaeck my R matrix. Code G
* routines.
* J1.6:17-01-1997|JRT|Changes at uni, QA G code. Moved QBD routines to
* ldqbd.c & h
* J1.7:20-01-1997|JRT|At Uni. Changed UD_G to algorithm 3.4.3.
* J1.8:29-01-1997|JRT|Uni. Major change to incorporate new method, viz;
* algorithm 3.4.3.
* J1.9:30-01-1997|JRT|uni. Changed q_gen so that system could be tested with
* another model.
* J1.10:02-02-1997|JRT|Uni. tidied up the code a little. Checked that everything
* was being freed correctly and removed redundant routines.
* J1.11:04-02-1997|JRT|Uni. major change to Gprocess. implemented my new memory
* efficient algorithm.  It works for G with the example given.
* J1.12:05-02-1997|JRT|Uni. Changed J1.11 algorithm to be generic for R and G.
* J1.13:16-05-1997|JRT|Uni. Moved the q_gen routines to a seperate file.
* J1.15:19-05-1997|JRT|
* J1.16:20-05-1997|JRT|Debugging code.
* J2.0:02-11-2002|JRT| Found bug in freeing same matrix twice, and not freeing
*           one of them
***********************************************************************
* SerialNo:	20021102
* Version:	J2.0
***********************************************************************/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<float.h>
#include"matutil.h"
#include"ldqbd.h"
#include"qgen.h" /* J1.13 */
/* extermal variables whose scope is this source code and ldqbd.c */
extern FILE *report;
/************************************************************************/
/* calculate U seperately for R from page 35. J1.1 REMOVED IN J1.12 */
/************************************************************************/
/* calculate U for G seperately from page 31. J1.5
REDUNDANT CODE REMOVED IN J1.10 */
/************************************************************************/
/* calculate D seperately from page 35. J1.1 REMOVED IN J1.12  */
/************************************************************************/
/* calculate D for G seperately from page 31. J1.5
redundant code removed in J1.1.0 */
/************************************************************************/
/* J1.12 modified to be passed UD(l-1, k) to calculate the new UD.
UD_R calculates the UD pair for the calculation af R and returns
a pointer to U(l,k) and pointer + 1 references D(l,k+2^(l+1)).
For the method used in this function refer to  Bright P44 with particular
reference to equation 3.4.6 & 3.4.7. The passed UD need be;
U0 = Uk, U1 = U(k+2^(l-1)), U2 = U(k+2^l)
D0 = D(k+2^(l+1)), D1 = D(k+2^l), D2 = D(k+3*2^(l-1)) */
matrix_ptr UD_R(matrix_ptr U, matrix_ptr D, int l,  int k, int dim)
{
matrix_ptr Qk, Qk1, Qk2, negQ, negQinv;
matrix_ptr UD = salloc(2);
matrix_ptr UD_inv = salloc(0);
matrix_ptr I = salloc(0);
matrix_ptr A = salloc(0);   /* J1.8 */
matrix_ptr B = salloc(0);
matrix_ptr C = salloc(0);
int two_l = 1;
int p;
int adim;
/* J1.1 */
/* begin */
/* J1.12 initialise two_l. initialised this way because l is small and the
function pow() uses double precision arithmetic (overkill). */
for (p=0; p<l; p++) two_l *= 2;
if (l==0) {
    Qk = q_gen(dim, k);
    Qk1 = q_gen(dim, k+1);
/* calc U */
    negQ = matscalmult(-1.0, *(Qk1+1));
    negQinv = mat_invert(negQ);
    A = matmult(*Qk, *negQinv);
    mat_copy(UD, A);
    free_matrix(A);
    free_matrix(negQ);
    free_matrix(negQinv);
/* calc D */
/* J1.16 */
    free_matrix_array(Qk, 3);
    Qk = q_gen(dim, k+2);
    Qk2 = q_gen(dim, k+1);
    negQ = matscalmult(-1.0, *(Qk2+1));
    negQinv = mat_invert(negQ);
    A = matmult(*(Qk+2), *negQinv);
    mat_copy((UD+1), A);
    free_matrix(A);
    free_matrix_array(Qk, 3);
    free_matrix_array(Qk1, 3);
    free_matrix_array(Qk2, 3);
    free_matrix(negQ);
    free_matrix(negQinv);
/*testies
fprintf(report, "l=%i: k=%i\n", l, k);
write_matrix(report, UD+1);
fprintf(report, "\n");
*/
} else {
/*testies
fprintf(report, "begin...l=%i: k=%i\n", l, k);
write_matrix(report, U);
fprintf(report, "\n");
write_matrix(report, U+1);
fprintf(report, "\n");
write_matrix(report, U+2);
fprintf(report, "\n");
write_matrix(report, (D+2));
fprintf(report, "\n");
write_matrix(report, D);
fprintf(report, "\n");
write_matrix(report, D+1);
fprintf(report, "\n");
fprintf(report, "end the three\n");
*/
    A = matmult(*(U+2), *(D+1));
    B = matmult(*D, *(U+1));
    adim = A->row;    /* NEED ERROR HANDLE */
    I = identity_gen(adim);
    C = mattake(*I,*A);
    free_matrix(A);
    A = mattake(*C,*B);
    free_matrix(B);
    free_matrix(C);
    UD_inv = mat_invert(A);
    free_matrix(A);
   /* calculate U(l,k) */
    A = matmult(*U, *(U+1));
    B = matmult(*A, *UD_inv);
    mat_copy(UD, B);
    free_matrix(A);
    free_matrix(B);
       /* calculate D(l,k+2^(l+1)) */
    A = matmult(*(D+2), *(D+1));
    B = matmult(*A, *UD_inv);
    mat_copy((UD+1), B);
    free_matrix(A);
    free_matrix(B);
    free_matrix(UD_inv); /* J1.5 */
} /* end else */
return UD;
} /* end UD_R */
/************************************************************************/
/* UD_G calculates the UD pair for the calculation af G and returns
a pointer to U(l,k-1+2^l) and pointer + 1 references D(l,k-1+2^l).
For the method used in this function refer to  Bright algorithm 3.4.3 P42.
J1.8 made this a lot simpler
The values of U and D need be (J1.12 changed this order);
U0 = U(k-1+2^(l-1)), U1 = U(k-1+2^l), U2 = U(k-1+3*2^(l-1))
D0 = D(k-1+2^(l-1)), D1 = D(k-1+2^l), D2 = D(k-1+3*2^(l-1))  */
matrix_ptr UD_G(matrix_ptr U, matrix_ptr D, int l, int k, int dim)
{
matrix_ptr Qk, negQ, negQinv;
matrix_ptr UD = salloc(2);
matrix_ptr UD_inv = salloc(0);
matrix_ptr I = salloc(0);
matrix_ptr A = salloc(0);   /* J1.8 */
matrix_ptr B = salloc(0);
matrix_ptr C = salloc(0);
int adim;
/* J1.1 */
/* begin */
if (l==0) {
    Qk = q_gen(dim, k);
/* calc U */
    negQ = matscalmult(-1.0, *(Qk+1));
    negQinv = mat_invert(negQ);
    A = matmult(*negQinv, *Qk);
    mat_copy(UD, A);
    free_matrix(A);
/* calc D */
    A = matmult(*negQinv, *(Qk+2));
    mat_copy((UD+1), A);
    free_matrix(A);
    free_matrix_array(Qk, 3);
    free_matrix(negQ);
    free_matrix(negQinv);
} else {
    A = matmult(*(U+1), *(D+2));
    B = matmult(*(D+1), *U);
    adim = A->row;    /* NEED ERROR HANDLE */
    I = identity_gen(adim);
    C = mattake(*I,*A);
    free_matrix(A);
    A = mattake(*C,*B);
    free_matrix(B);
    free_matrix(C);
    UD_inv = mat_invert(A);
    free_matrix(A);
   /* calculate U(l,k-1+2^l) */
    A = matmult(*(U+1), *(U+2));
    B = matmult(*UD_inv, *A);
    mat_copy(UD, B);
    free_matrix(A);
    free_matrix(B);
       /* calculate D(l,k+2^(l+1)) */
    A = matmult(*(D+1), *D);
    B = matmult(*UD_inv, *A);
    mat_copy((UD+1), B);
    free_matrix(A);
    free_matrix(B);
    free_matrix(UD_inv); /* J1.5 */
} /* end else */
return UD;
} /* end UD_G */
/***********************************************************************/
/* compR calculates the R matrix. Adapted from alg 3.4.4 */
matrix_ptr compR(float *diff, int l, matrix_ptr U, matrix_ptr D)
{
static matrix_ptr PI;
static matrix_ptr Rnew; /* value is retained between calls */
matrix_ptr Rold = salloc(1);
matrix_ptr A, B; /* comes in handy */
int udim;
/* begin */
if (l==0) {
   Rnew = salloc(1);
   mat_copy(Rnew, U);
   PI = salloc(1);
   mat_copy(PI, D);
   *diff = FLT_MAX; /* the first iteration has no
		  previous value to compare,
		   so the difference is big  */
   free_matrix(Rold);
} else {
   mat_copy(Rold, Rnew);
   A =  matmult(*U, *PI);
   B = matadd(*Rold, *A);
   *diff = matabs(A);
   free_matrix(A);
   A = mattake(*Rold, *Rnew);
   mat_copy(Rnew, B);
   free_matrix(B);
   A = matmult(*D, *PI);
   mat_copy(PI, A);
   free_matrix(A);
} /* end if */
return Rnew;
}/* end compR */
/***********************************************************************/
/* Use equation 2.5.4 on page 15 (rearranged) to calculate R(k) from
R(k+1). Iterate backeards */
matrix_ptr backR(matrix_ptr r, int k, int dim)
{
matrix_ptr Qk1 = salloc(0);
matrix_ptr Qk2 = salloc(0);
matrix_ptr Qk  = salloc(0);
matrix_ptr result = salloc(0);
matrix_ptr A = salloc(0);
matrix_ptr B = salloc(0);
int rdim = r->col; /* R should be square */
/* begin */
    Qk = q_gen(dim, k);
    Qk1 = q_gen(dim, k+1);
    Qk2 = q_gen(dim, k+2);
    A = matmult(*r, *(Qk2+2));
    B = matadd(*(Qk1+1), *A);
    free_matrix(A);
    A = mat_invert(B);
    free_matrix(B);
    B = matscalmult(-1.0, *Qk);
    result = matmult(*B, *A);

    /* now do final house keeping */
    free_matrix_array(Qk, 3);
    free_matrix_array(Qk1, 3);
    free_matrix_array(Qk2, 3); /* J2.0 use the routine I created */
    free_matrix(A);
    free_matrix(B);
return result;
} /* end backR */
/***********************************************************************/
/* compG calculates the G matrix. Adapted from alg 3.4.2 P40 J1.8 lots
of changes. */
matrix_ptr compG(float *diff, int l, matrix_ptr U, matrix_ptr D)
{
float tol = 1.0e-3; /* the comparison value used in truncating G */
static matrix_ptr PI;
static matrix_ptr Gnew; /* value is retained between calls */
matrix_ptr Gold = salloc(1);
matrix_ptr A, B; /* comes in handy */
int udim;
/* begin */
if (l==0) 
{
    Gnew = salloc(1);
    mat_copy(Gnew, D);
    PI = salloc(1);
    mat_copy(PI, U);
    *diff = FLT_MAX; /* the first iteration has no
		  previous value to compare,
		   so the difference is big  */
    free_matrix(Gold); /* J1.12 */
} else {
    mat_copy(Gold, Gnew);
    A =  matmult(*PI, *D);
    *diff = matabs(A); /* J1.16 This is the diffeence in G(l) */
    /* testies  fprintf(report, "iteration %i: difference %5.2e.\n", l, *diff); */
    B = matadd(*Gold, *A);
    free_matrix(A);
    mat_copy(Gnew, B);
    free_matrix(B);
    A = matmult(*PI, *U);
    mat_copy(PI, A);
    free_matrix(A);
    free_matrix(Gold); /* J1.12 */
} /* end if */
return Gnew;
}/* end compG */
/***********************************************************************/
/* Use equation 2.7.14 on page 21 (rearranged) to calculate G(k) from
G(k+1). Iterate backwards */
matrix_ptr backG(matrix_ptr g, int k, int dim)
{
matrix_ptr Qk = salloc(0);
matrix_ptr result = salloc(0);
matrix_ptr A = salloc(0);
matrix_ptr B = salloc(0);
int rdim = g->col; /* G should be sqwuare */
/* begin */
    Qk = q_gen(dim, k);
    A = matmult(*Qk, *g);
    B = matadd(*(Qk+1), *A);
    free_matrix(A);
    A = matscalmult(-1.0, *B);
    free_matrix(B);
    B = mat_invert(A);
    result = matmult(*B, *(Qk+2));
    /* now do final house keeping J2.1 use free_marix_array */
    free_matrix_array(Qk, 3);
    free_matrix(A);
    free_matrix(B);
return result;
} /* end backG */
/************************************************************************/
/* GfromR calculates the G matrix at k from Rk (see p21 eq 2.7.12 */
matrix_ptr GfromR(matrix_ptr Rk, int k, int dim)
{
matrix_ptr A = salloc(0);
matrix_ptr B = salloc(0);
matrix_ptr C = salloc(0);
matrix_ptr result = salloc(0);
matrix_ptr Qk = salloc(0);
matrix_ptr Qk1 = salloc(0);
/* begin */
    Qk = q_gen(dim, k);
    Qk1 = q_gen(dim, k+1);
    A = matscalmult(-1.0, *(Qk+1));
    B = matmult(*Rk, *(Qk1+2));
    C = mattake(*A, *B);
    free_matrix(A);
    free_matrix(B);	
    A = mat_invert(C);
    free_matrix(C);
    result = matmult(*A, *(Qk+2));
    free_matrix(A);
    /* J2.0 use free matrix array function */
    free_matrix_array(Qk, 3);
    free_matrix_array(Qk1, 3);
return result;
} /* end GfromR */
/***************************************************************************/
/* J1.11 These variables have scope over the following functions */
matrix_ptr Uouter, Douter;  /* The LHS UD on l */
matrix_ptr Uinner, Dinner;  /* The centre UD on l*/
matrix_ptr Unew, Dnew;  /* The RHS UD becomes
				    the new LHS UD on l */
int noUD[MAX_L]; /* The number of UD stored for particular l */
int phasej, leveli; /* i & j were given names a little more elaborate as they
	      now have more scope and thus more likely to get confused
	      */
int two_l, two_i; /* two to the power of l and i respectively */
matrix_ptr A; /* helpful storage */
int l; /* Aw come on! Everyone knows what this is */

/***********************************************************************/
/* J1.11 RGprocess, compR/G and get_outer run hand in hand */
matrix_ptr RGprocess(UDtype type, int k, int dim)
{
float diff; /* The scalar difference between RG matricies in successive
       iterations */
matrix_ptr RGl = salloc(0); /* stores the iterative value of R or G */
int p;
/* begin */
l = 0;
two_l = 1;
Uouter = salloc(MAX_L);
Douter = salloc(MAX_L);
Uinner = salloc(MAX_L);
Dinner = salloc(MAX_L);
Unew = salloc(MAX_L);
Dnew = salloc(MAX_L);
if (type == G)  /* J1.12 */
    A = UD_G((matrix_ptr) NULL, (matrix_ptr) NULL, l, k, dim);
else
    A = UD_R((matrix_ptr) NULL, (matrix_ptr) NULL, l, k, dim);
mat_copy(Uouter, A);
mat_copy(Douter, (A+1));
free_matrix(A);

if (type == G) /*J1.12 */
    RGl = compG(&diff, l, Uouter, Douter);
else
    RGl = compR(&diff, l, Uouter, Douter);
while (diff > TOL && l < MAX_L) { /* J1.11 */
        leveli = 0;
        two_i = 1;
        for (p=0; p<MAX_L; p++)
            noUD[p] = 1; /* There will always be one outer UD on a level
                            before the next level UD need be calculated */
        phasej = k - 1 + 2*two_l;
        while (leveli <= l)
            get_outer(dim, type);
        l++;
        two_l *= 2;
        if (type == G)   /* J1.12 */
            RGl = compG(&diff, l, (Uouter+l), (Douter+l));
        else
            RGl = compR(&diff, l, (Uouter+l), (Douter+l));
    } /* end while */
if (type==G)
   printf("The number of iterations for G was %i at k = %i\n", l, k);
else
   printf("The number of iterations for R was %i at k = %i\n", l, k);
   
    free_matrix_array(Uouter, MAX_L);
    free_matrix_array(Douter, MAX_L);
    free_matrix_array(Uinner, MAX_L);
    free_matrix_array(Dinner, MAX_L);
    free_matrix_array(Unew, MAX_L);
    free_matrix_array(Dnew, MAX_L);
    return RGl;
} /* end RGprocess */

/*************************************************************************/
/* these variables are defined outside the process as it is recursive, and
   any variables local to the process are stored on the stack for each recursive
   call, so it was designed to utilise external variables. */
int m; /* counter for the number of matricies stored due to for j */
int stval, eval; /* the start and end values respectively of the j loop */
matrix_ptr Dval;
matrix_ptr Uval; /* the values of UD for the calculation of UD at
                                sthe next value l */
int freed; /* J2.0 this is non zero if Uval and Dval need to be freed */
/*************************************************************************/
/* get_outer runs in conjunction with G/Rprocess to find the outer UD matricies
and store them in memory pointed to by the matrix _outer, where the index of the
matrix appertains to the level of l that the UD belong to. J1.11 */
void get_outer(int dim, UDtype type)
{
/* begin */
Dval = salloc(3);
Uval = salloc(3);
freed = 0; /* J2.0 freed is false, as they have just been allocated memory */
    switch (noUD[leveli]) {
    case 1:
        if (leveli > 0) 
            set_UDval(leveli-1);
           if (type==G)
              A = UD_G(Uval, Dval, leveli, phasej, dim);
           else
              A = UD_R(Uval, Dval, leveli, phasej, dim); /* J1.16 */
        mat_copy(Uinner+leveli, A);
        mat_copy(Dinner+leveli, (A+1));
        free_matrix(A);
        noUD[leveli] = 2;
        phasej += two_i;
        if (leveli>0) {
            leveli -= 1;
            two_i /= 2;
            if (type == R) phasej += two_i; /* J1.12 */
            shuffleUD(leveli);
        }
        break;
    case 2:
        if (leveli > 0) 
            set_UDval(leveli-1);
           if (type==G)
              A = UD_G(Uval, Dval, leveli, phasej, dim);
           else
              A = UD_R(Uval, Dval, leveli, phasej, dim); /* J1.16 */
        mat_copy(Unew+leveli, A);
        mat_copy(Dnew+leveli, (A+1));
        free_matrix(A);
        if (leveli > 0)
            shuffleUD(leveli-1);
        if (type == G) /* J1.12 */
            phasej -= two_i;
        else
            phasej -= 2*two_i;
        leveli += 1;
        two_i *= 2;
        if (leveli <= l) {
            get_outer(dim, type); /* here is the one and only recursive call */
            if (leveli <= l && leveli > 0) {
                leveli -= 1;
                two_i /= 2;
                if (type == R) phasej += two_i; /* J1.12 */
            /* not sure that this is the correct amount to increment */
            } /* end if */
        } else {
           if (leveli > 0) 
               set_UDval(leveli-1);
           if (type==G)
              A = UD_G(Uval, Dval, leveli, phasej, dim);
           else
              A = UD_R(Uval, Dval, leveli, phasej, dim); /* J1.16 */
       	   shuffleUD(leveli-1); /* J1.16 */
            mat_copy(Uouter+leveli, A);
            mat_copy(Douter+leveli, (A+1));
            free_matrix(A);
        } /* end if */
        break;
    default:
        printf("ERROR!!!! in switch get_outer\n");
        exit(1);
        break;
    } /* end switch */
    /* J2.0 because of the recursive call, and the fact that Uval and Dval are 
       global, they tend to get freed twice. Brought in a boolean variable to detect
       if they have been freed or not. this is also global */
    if (!freed)
    {
        free_matrix_array(Uval, 3);
        free_matrix_array(Dval, 3); /* J2.0 bug here. Uval was freed twice */
        freed = 1;
    }
    return;
} /* end get_outer */
/*************************************************************************/
/* for use with get_outer, shuffleUD shuffles the UD pointers within the
   three sets of matricies */
void shuffleUD(int i)
{
    extern struct matrix_type null_mat; /* defined in matutil.c */
    mat_copy(Uouter+i, Unew+i);
    mat_copy(Douter+i, Dnew+i);
    noUD[i]=1;  /* J1.15 */
    return;
} /* end shuffleUD */
/**************************************************************************/
/* chuck the correct matriceis in Uval and Dval for the calculation of the
   next UD. J1.12 changed order. */
void set_UDval(int i)
{
    mat_copy(Uval, Uouter+i);
    mat_copy(Uval+1, Uinner+i);
    mat_copy(Uval+2, Unew+i);
    mat_copy(Dval, Douter+i);
    mat_copy(Dval+1, Dinner+i);
    mat_copy(Dval+2, Dnew+i);
    return;
} /* end set_UDval */
