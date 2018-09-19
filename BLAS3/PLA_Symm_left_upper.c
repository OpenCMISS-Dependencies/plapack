#include "PLA.h"

int PLA_Symm_left_upper ( PLA_Obj alpha, PLA_Obj A, PLA_Obj B, PLA_Obj C,
                          int nb_alg )
{
  PLA_Obj     ATL=NULL, ATR=NULL,     A00=NULL, A01=NULL, A02=NULL,  BT=NULL,  B0=NULL,  CT=NULL,  C0=NULL,
              ABL=NULL, ABR=NULL,     A10=NULL, A11=NULL, A12=NULL,  BB=NULL,  B1=NULL,  CB=NULL,  C1=NULL,
                                      A20=NULL, A21=NULL, A22=NULL,       B2=NULL,       C2=NULL,
             ONE=NULL, 
             B1_dpmv=NULL, C1_dpmv=NULL, A01_dpmv=NULL, A12_dpmv=NULL, A11_dmsc=NULL;
  int         b;

  PLA_Create_constants_conf_to( A, NULL, NULL, &ONE );

  PLA_Part_2x2( A,  &ATL, /**/ &ATR,
                  /* ************** */
                    &ABL, /**/ &ABR,   
	        0, 0,     /* submatrix */ PLA_TL );

  PLA_Part_2x1( B,  &BT, 
                   /***/
                    &BB,               
                0, /* length submatrix */ PLA_TOP );

  PLA_Part_2x1( C,  &CT, 
                   /***/
                    &CB,  
	        0, /* length submatrix */ PLA_TOP );

  while ( b = min( PLA_Obj_length( ABR ), nb_alg ) ){
    PLA_Repart_2x2_to_3x3( ATL, /**/ ATR,         &A00, /**/ &A01, &A02,
                        /* ************* */    /* ********************* */
                                /**/              &A10, /**/ &A11, &A12,
                           ABL, /**/ ABR,         &A20, /**/ &A21, &A22,   
                           b, b, /* A11 from */ PLA_BR );

    PLA_Repart_2x1_to_3x1( BT,                    &B0,
                          /**/                    /**/
                                                  &B1,
                           BB,                    &B2,                    
                           b, /* length B1 from */ PLA_BOTTOM );

    PLA_Repart_2x1_to_3x1( CT,                    &C0,
                          /**/                    /**/
                                                  &C1,
                           CB,                    &C2,                     
                           b, /* length C1 from */ PLA_BOTTOM );

    /* ********************************************************************* */

    PLA_Obj_set_orientation( A12, PLA_PROJ_ONTO_ROW );
    PLA_Obj_set_orientation( B1,  PLA_PROJ_ONTO_ROW );
    PLA_Obj_set_orientation( C1,  PLA_PROJ_ONTO_ROW );

    PLA_Pmvector_create_conf_to( C0, PLA_PROJ_ONTO_COL, PLA_ALL_COLS, 
				 b, &A01_dpmv );
    PLA_Pmvector_create_conf_to( C2, PLA_PROJ_ONTO_COL, PLA_ALL_COLS, 
				 b, &A12_dpmv );
    PLA_Pmvector_create_conf_to( C0, PLA_PROJ_ONTO_ROW, PLA_ALL_ROWS, 
				 b, &B1_dpmv );
    PLA_Pmvector_create_conf_to( C0, PLA_PROJ_ONTO_ROW, PLA_ALL_ROWS, 
				 b, &C1_dpmv );
    PLA_Mscalar_create_conf_to ( A11, PLA_ALL_ROWS, PLA_ALL_COLS, &A11_dmsc );

    PLA_Copy( A01, A01_dpmv );
    PLA_Copy( A12, A12_dpmv );
    PLA_Copy( B1,  B1_dpmv );
    PLA_Copy( C1,  C1_dpmv );
    PLA_Copy( A11, A11_dmsc );

    PLA_Local_gemm( PLA_NO_TRANSPOSE, PLA_NO_TRANSPOSE, 
		    alpha, A01_dpmv, B1_dpmv, ONE, C0 );

    PLA_Local_gemm( PLA_NO_TRANSPOSE, PLA_NO_TRANSPOSE, 
		    alpha, A12_dpmv, B1_dpmv, ONE, C2 );

    PLA_Local_symm( PLA_SIDE_LEFT, PLA_UPPER_TRIANGULAR,
		    alpha, A11_dmsc, B1_dpmv, ONE, C1_dpmv );

    PLA_Copy( C1_dpmv, C1 );

    /* ********************************************************************* */

    PLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,         A00, A01, /**/ A02,
                                    /**/               A10, A11, /**/ A12,
                            /* ************** */   /* ****************** */
                              &ABL, /**/ &ABR,         A20, A21, /**/ A22, 
                              /* A11 added to */ PLA_TL );

    PLA_Cont_with_3x1_to_2x1( &BT,                     B0,
                                                       B1,
                             /***/                    /**/
                              &BB,                     B2,                 
                              /* B1  added to */ PLA_TOP );

    PLA_Cont_with_3x1_to_2x1( &CT,                     C0,
                                                       C1,
                             /***/                    /**/
                              &CB,                     C2,                 
                              /* A1  added to */ PLA_TOP );
  }

  PLA_Obj_free( &ATL );
  PLA_Obj_free( &ATR );
  PLA_Obj_free( &A00 );
  PLA_Obj_free( &A01 );
  PLA_Obj_free( &A02 );
  PLA_Obj_free( &BT );
  PLA_Obj_free( &B0 );
  PLA_Obj_free( &CT );
  PLA_Obj_free( &C0 );
  PLA_Obj_free( &ABL );
  PLA_Obj_free( &ABR );
  PLA_Obj_free( &A10 );
  PLA_Obj_free( &A11 );
  PLA_Obj_free( &A12 );
  PLA_Obj_free( &BB );
  PLA_Obj_free( &B1 );
  PLA_Obj_free( &CB );
  PLA_Obj_free( &C1 );
  PLA_Obj_free( &A20 );
  PLA_Obj_free( &A21 );
  PLA_Obj_free( &A22 );
  PLA_Obj_free( &B2 );
  PLA_Obj_free( &C2 );
  PLA_Obj_free( &ONE );
  PLA_Obj_free( &B1_dpmv );
  PLA_Obj_free( &C1_dpmv );
  PLA_Obj_free( &A01_dpmv );
  PLA_Obj_free( &A12_dpmv );
  PLA_Obj_free( &A11_dmsc );

  return PLA_SUCCESS;
}
