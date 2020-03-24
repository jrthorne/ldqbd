/**********************************************************************
* File name: matmain.c
* Date start/finish: 21-12-1996/dd-mm-yyyy
* By: Jason Thorne
* Name: Main Code
* Project: Matrix analytic methods in Quasi-stationary birth-death
*      processes Honours project 1996-1997
* Description: The main code to test the routines derived from Bright's thesis
* Updates:
* JX.X:DD-MM-YYYY|XXX|Description
* J1.1:12-01-1997|JRT|Wrote compR, adapted from alg3.3.4 of L.Bright's
*                     thesis. Added functions calcU and calcD to calculate
*                     these values seperately.
* J1.2:13-01-1997|JRT|Changes at home.
* J1.3:14-01-1997|JRT|Changes at uni
* J1.4:15-01-1997|JRT|Changes at Uni
* J1.5:16-01-1997|JRT|Changes at uni. Code to chaeck my R matrix. Code G
routines.
* J1.6:17-01-1997|JRT|Changes at uni, QA G code. Moved QBD routines to
* ldqbd.c & h
* J2.0:01-11-2002|JRT| Makes use of the new matutil which will calculate
*                      an inverse.
***********************************************************************
* SerialNo:	20021101
* Version:	J2.0
***********************************************************************/
#include<stdio.h>
#include<time.h>
#include<stdlib.h>
#include<math.h>
#include"matutil.h"
#include"ldqbd.h"  /* J1.6 */
#include"qgen.h"
    FILE *report;
    FILE *data;
/* external variables whose scope is this source code and ldqbd.c */
    matrix_ptr a, s, s0;
  /* insert from J1.8 */
    int dime = 8;
    int levk = 100;
    MATTYPE lambda = 1.0;
    MATTYPE theta = 0.05;
    MATTYPE rho = 0.5;
    MATTYPE mu;
/* J1.8 this satisfies lambda/(c*mu) < 1 (P97) */
/* the main code */
void main(void)
{
    time_t stm, ftm;
    MATTYPE takentm; /* for timing the routine */
    matrix_ptr r = salloc(0);
    matrix_ptr g = salloc(0);
    int l = 0;
    int b, i, j;
    matrix_ptr check = salloc(0);
    matrix_ptr checkR = salloc(0);
    matrix_ptr difference;
    MATTYPE diff;
    mu = lambda/((dime-1)*rho) * 1.1; 
/* the data and report files */
    report = fopen("tst.rep", "w");
    data   = fopen("qbd.dat", "r");
/* begin */
    stm = time(&stm);
    g = RGprocess(G, levk, dime);      
    ftm = time(&ftm);
    takentm = difftime(ftm, stm);
    fprintf(report, "The matrix G at k = %i by LR algorithm\n",levk);
    fprintf(report, "took %10.3f seconds to calculate as:\n", takentm);
    write_matrix(report, g); 
    fprintf(report, "\n");
    stm = time(&stm);
    r = RGprocess(R, levk, dime);      /* calculate R(100) */
    ftm = time(&ftm);
    takentm = difftime(ftm, stm);
    
    fprintf(report, "The matrix R at k = %i by LR algorithm\n",levk);
    fprintf(report, "took %10.3f seconds to calculate as:\n", takentm);
    write_matrix(report, r); 
    checkR = GfromR(r, levk, dime);
    fprintf(report, "The matrix G at k = %i calculated from R(%i) \n",
levk,levk);
    write_matrix(report, checkR);
    fprintf(report, "\n");
    difference = mattake(*g, *checkR);
    diff = matabs(difference);
    fprintf(report, " The error in G is %5.2e \n", diff);
    
    check = RGprocess(R, levk+1, dime);
    checkR = backR(check,levk,dime);
    difference = mattake(*r, *checkR);
    diff = matabs(difference);
    fprintf(report, " The error in R is %5.2e \n", diff);
    free_matrix(checkR);
    free_matrix(check);
    free_matrix(g);
/* finalise */
    fprintf(report, "****END OF REPORT****");
    fclose(report);
    fclose(data);
} /* end main */
