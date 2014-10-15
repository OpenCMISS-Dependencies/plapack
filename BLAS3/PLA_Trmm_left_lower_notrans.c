/*
   PLAPACK Release 3.1
   
   Copyright (C) 2001
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

void PLA_Trmm_left_lower_notrans( int diag, PLA_Obj alpha, PLA_Obj A, PLA_Obj B, int nb_alg )
{
  PLA_Obj     ATL=NULL, ATR=NULL,     A00=NULL, A01=NULL, A02=NULL,    BT=NULL,           B0=NULL,
              ABL=NULL, ABR=NULL,     A10=NULL, A11=NULL, A12=NULL,    BB=NULL,           B1=NULL,
                                      A20=NULL, A21=NULL, A22=NULL,                  B2=NULL,
              A11_dup=NULL, B1_mv=NULL, A21_dpmv=NULL, B1_dpmv=NULL,
              ONE=NULL;

  int         b, size_row;

  PLA_Create_constants_conf_to( A, NULL, NULL, &ONE );

  PLA_Part_2x2( A,  &ATL, /**/ &ATR,
                  /* ************** */
                    &ABL, /**/ &ABR,   0, 0,     /* submatrix */ PLA_BR );

  PLA_Part_2x1( B,  &BT, 
                   /***/
                    &BB,               0, /* length submatrix */ PLA_BOTTOM );

  while ( TRUE ){
    PLA_Obj_global_length( ATL, &size_row );
    b = min( size_row, nb_alg );
    if ( 0 == b ) break;

    PLA_Repart_2x2_to_3x3( ATL, /**/ ATR,         &A00, &A01, /**/ &A02,
                                /**/              &A10, &A11, /**/ &A12,
                        /* ************* */    /* ********************* */
                           ABL, /**/ ABR,         &A20, &A21, /**/ &A22,   
                           b, b, /* A11 from */ PLA_TL );

    PLA_Repart_2x1_to_3x1( BT,                   &B0,
                                                 &B1,
                          /**/                   /**/
                           BB,                   &B2,    
                           b, /* length B1 from */ PLA_TOP );

    /* ********************************************************************* */

    PLA_Obj_set_orientation( B1, PLA_PROJ_ONTO_ROW );

    PLA_Pmvector_create_conf_to( B2, PLA_PROJ_ONTO_ROW, PLA_ALL_ROWS, b, &B1_dpmv );

    PLA_Pmvector_create_conf_to( B2, PLA_PROJ_ONTO_COL, PLA_ALL_COLS, b, &A21_dpmv );

    PLA_Copy( B1, B1_dpmv );

    PLA_Copy( A21, A21_dpmv );

    PLA_Local_gemm( PLA_NO_TRANSPOSE, PLA_NO_TRANSPOSE, alpha, A21_dpmv, B1_dpmv, ONE, B2 );

    PLA_Mscalar_create_conf_to( A11, PLA_ALL_ROWS, PLA_ALL_COLS, &A11_dup );

    PLA_Copy( A11, A11_dup );

    PLA_Local_trmm( PLA_LEFT, PLA_LOWER_TRIANGULAR, PLA_NO_TRANSPOSE, diag, alpha, A11_dup, B1_dpmv );
    
    PLA_Copy( B1_dpmv, B1 );

    /* ********************************************************************* */

    PLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,         A00, /**/ A01, A02,
                            /* ************** */   /* ****************** */
                                    /**/               A10, /**/ A11, A12,
                              &ABL, /**/ &ABR,         A20, /**/ A21, A22, 
			      /* A11 added to */ PLA_BR );

    PLA_Cont_with_3x1_to_2x1( &BT,                     B0,
                             /***/                    /**/
                                                       B1,
                              &BB,                     B2,                 
			      /* B1  added to */ PLA_BOTTOM );
  }

  PLA_Obj_free( &ATL );
  PLA_Obj_free( &ATR );
  PLA_Obj_free( &A00 );
  PLA_Obj_free( &A01 );
  PLA_Obj_free( &A02 );
  PLA_Obj_free( &BT );
  PLA_Obj_free( &B0 );
  PLA_Obj_free( &ABL );
  PLA_Obj_free( &ABR );
  PLA_Obj_free( &A10 );
  PLA_Obj_free( &A11 );
  PLA_Obj_free( &A11_dup );
  PLA_Obj_free( &A12 );
  PLA_Obj_free( &BB );
  PLA_Obj_free( &B1 );
  PLA_Obj_free( &B1_mv );
  PLA_Obj_free( &B1_dpmv );
  PLA_Obj_free( &A20 );
  PLA_Obj_free( &A21 );
  PLA_Obj_free( &A21_dpmv );
  PLA_Obj_free( &A22 );
  PLA_Obj_free( &B2 );
  PLA_Obj_free( &ONE );
}
