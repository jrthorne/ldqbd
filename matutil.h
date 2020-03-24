
/**********************************************************************
* File name: matutil.h
* Date start/finish: 20-12-1996/dd-mm-yyyy
* By: Jason Thorne
* Name: Math Utilities
* Project: Matrix analytic methods in Quasi-stationary birth-death
*          processes Honours project 1996-1997
* Description: This is the header file forthe mathematic utilities
* needed to implement the methods as described in Les Brights thesis
* Updates: see matutil.c
*****************************************************************
* SerialNo:		20021101
***********************************************************************/
#ifndef _MATUTIL_
#define _MATUTIL_
#include<stdio.h>
#include<stdlib.h>
#include<stdarg.h>
#define MATTYPE float /* J2.1 the matrix type is either MATTYPE or float */
/***********************************************************************/
/* constants */
/* wrote some code to find that 2**53 = 1.11e-16 is smallest no that
   1.0 + CLOSEZERO == 1.0 for MATTYPE precision */
#define CLOSEZERO 5.0e-16 /* a very small number */
/***********************************************************************/
/* Type definitions */
typedef struct matrix_type *matrix_ptr;
struct matrix_type {
        int row, col; /* The dimention of the matrix */
        MATTYPE *element; /* The elements are stored in
                            typical C style as a
                            continuous linked list */
        };
/***********************************************************************/
/* routine templates */
/* This routine multiplies two matricies together and the resultant matrix
   is returned as the structure matrix */
matrix_ptr matmult(struct matrix_type, struct matrix_type);
/* This routine writes a matrix to a file. The file pointer is global */
void write_matrix(FILE *fp, matrix_ptr a);
/* This routine is the assignment operator for the matrix structure
   the number of rows followed by the number of columns is sent followed
   by a variable length argument list of the elements */
matrix_ptr mat_assign(int row, int col, ...);
/* free the element list belonging to the matrix */
void free_matrix(matrix_ptr a);
/* add matricies together */
matrix_ptr matadd(struct matrix_type, struct matrix_type);
/* take matricies from each other */
matrix_ptr mattake(struct matrix_type, struct matrix_type);
/* scalar multiply */
matrix_ptr matscalmult(MATTYPE, struct matrix_type);
/* This routine reads a matrix from a file. The structure used in the file
   is first the number of rows. followed by the number of columns, then the
   matrix.
*/
matrix_ptr mat_read(FILE *);
/* This function generates the identity matrix (allways square) */
matrix_ptr identity_gen(int dim);
/* Call the two nag routines f07adf and f07ajf to reduce the given
   matrix to LU form and then calculate the inverse respectively.
   The result overwrites the passed parameter *a */
matrix_ptr mat_invert(matrix_ptr a);
/***********************************************************************/
/* Copy the matrix structure */
void mat_copy(matrix_ptr to, matrix_ptr from);
/* return a scalar value representing the matrix. At this stage just
   sum all the absolute values of the elements of the matrix divided by the
   number of elements */
MATTYPE matabs(matrix_ptr);
/* J1.2 allocate memory and initialise the element pointer */
matrix_ptr salloc(int n);
/* J1.4 free array of matricies./ no_matricies is the same no as in salloc */
void free_matrix_array(matrix_ptr A, int no_matricies);
/**************************************************************************/
/* J1.6 return absolute value of the MATTYPEing point number */
MATTYPE myabs(MATTYPE x);
/***************************************************************/
/* J1.12 This functions performs an elementry arithmetic row operation 
   on a matrix. The sum is: A(lhs) = A(lhs) + coeff*A(rhs)
   where "A(i)" = row i of matrix A */
/***************************************************************/
void arith_elem_row_op(matrix_ptr A, int lhs, int rhs, MATTYPE coeff);
/***************************************************************/
/* J1.12 This function returns the given matrix augmented with 
   the identity matrix for finding the inverse matrix using
   row operations */
/***************************************************************/
matrix_ptr augment_with_I(matrix_ptr A);
/***************************************************************/
/* J1.13 This function returns the given matrix deaugmented 
   Used for finding the inverse matrix using
   row operations */
/***************************************************************/
matrix_ptr deaugment(matrix_ptr A);
/****************************************************************
* J1.13 this function is the elementary row operation swap.
* it changes the matrix passed to it and returns nothing
****************************************************************/
void swap_elem_row_op(matrix_ptr A, int xrow, int yrow);
/****************************************************************
* J2.1 this function is the elementary column operation swap
* col x for col y
* it changes the matrix passed to it and returns nothing
****************************************************************/
void swap_elem_col_op(matrix_ptr A, int x, int y);

#endif
