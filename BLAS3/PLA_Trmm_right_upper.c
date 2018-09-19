/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Trmm_right_upper( int transa, int diag,
			  PLA_Obj alpha, PLA_Obj A, PLA_Obj B )
{
  PLA_Obj 
    A_TL = NULL, A_BR = NULL, A_1 = NULL, A_1_dup = NULL, A_11_dup = NULL, 
    B_L = NULL,  B_R = NULL,  B_1 = NULL, B_1_dup = NULL,
    one = NULL,  zero = NULL;

  int 
    nb_alg, size;

  PLA_Template
    templ;

  PLA_Local_scal( alpha, B );   
  PLA_Create_constants_conf_to( A, NULL, &zero, &one );
  
  PLA_Obj_global_length( B, &size );
  nb_alg = min( size, nb_alg );

  if ( transa == PLA_NO_TRANS ){
    PLA_Obj_template( A, &templ );
    PLA_Environ_nb_alg( PLA_OP_MAT_PAN, templ, &nb_alg );

    PLA_Obj_view_all( A, &A_TL );        PLA_Obj_view_all( B, &B_L );

    PLA_Pmvector_create_conf_to( 
  	  B, PLA_PROJ_ONTO_COL, PLA_ALL_COLS, nb_alg, &B_1_dup );
    PLA_Pmvector_create_conf_to( 
          B, PLA_PROJ_ONTO_ROW, PLA_ALL_ROWS, nb_alg, &A_1_dup );

    while ( TRUE ){
      PLA_Obj_global_length( A_TL, &size );
      if ( ( size = min( size, nb_alg ) ) == 0 ) break;

      PLA_Obj_vert_split_2( A_TL, -size, PLA_DUMMY, &A_1 );

      if ( size != nb_alg ){
	PLA_Obj_horz_split_2( A_1_dup, size, &A_1_dup,
                                             PLA_DUMMY );
	PLA_Obj_vert_split_2( B_1_dup, size, &B_1_dup, PLA_DUMMY );
      }

      PLA_Copy( A_1, A_1_dup );       

      PLA_Obj_vert_split_2( A_1_dup, -size, PLA_DUMMY, &A_11_dup );

      PLA_Set_triang_to_zero( PLA_LOWER_TRIANGULAR, diag, A_11_dup );

      PLA_Local_gemm( PLA_NO_TRANS, PLA_TRANS, one,  B_L,  A_1_dup,
                                               zero, B_1_dup );

      PLA_Obj_vert_split_2( B_L, -size, &B_L, &B_1 );

      PLA_Reduce( B_1_dup, MPI_SUM, B_1 );

      PLA_Obj_split_4( A_TL, -size, -size, &A_TL,      PLA_DUMMY,
                                            PLA_DUMMY, PLA_DUMMY );

      PLA_Obj_vert_split_2( A_1_dup, -size, &A_1_dup, PLA_DUMMY );

    }
  }
  else if ( transa == PLA_TRANS ){
    PLA_Obj_template( A, &templ );
    PLA_Environ_nb_alg( PLA_OP_MAT_PAN, templ, &nb_alg );

    PLA_Obj_view_all( A, &A_BR );        PLA_Obj_view_all( B, &B_R );

    PLA_Pmvector_create_conf_to( 
  	  B, PLA_PROJ_ONTO_COL, PLA_ALL_COLS, nb_alg, &B_1_dup );
    PLA_Pmvector_create_conf_to( 
          B, PLA_PROJ_ONTO_ROW, PLA_ALL_ROWS, nb_alg, &A_1_dup );

    while ( TRUE ){
      PLA_Obj_global_length( A_BR, &size );
      if ( ( size = min( size, nb_alg ) ) == 0 ) break;

      PLA_Obj_horz_split_2( A_BR, size, &A_1,
                                        PLA_DUMMY );

      PLA_Obj_set_orientation( A_1, PLA_PROJ_ONTO_ROW );

      if ( size != nb_alg ){
	PLA_Obj_horz_split_2( A_1_dup, size, &A_1_dup,
                                             PLA_DUMMY );
	PLA_Obj_vert_split_2( B_1_dup, size, &B_1_dup, PLA_DUMMY );
      }

      PLA_Copy( A_1, A_1_dup );       

      PLA_Obj_vert_split_2( A_1_dup, size, &A_11_dup, PLA_DUMMY );

      PLA_Set_triang_to_zero( PLA_UPPER_TRIANGULAR, diag, A_11_dup );

      PLA_Local_gemm( PLA_NO_TRANS, PLA_TRANS, one,  B_R,  A_1_dup,
                                               zero, B_1_dup );

      PLA_Obj_vert_split_2( B_R, size, &B_1, &B_R );

      PLA_Reduce( B_1_dup, MPI_SUM, B_1 );

      PLA_Obj_split_4( A_BR, size, size,    PLA_DUMMY, PLA_DUMMY,
                                            PLA_DUMMY, &A_BR );

      PLA_Obj_vert_split_2( A_1_dup, size, PLA_DUMMY, &A_1_dup );

    }
  }
  else 
    PLA_Abort("Trmm: transpose case not yet implemented", __LINE__, __FILE__ );

  PLA_Obj_free( &A_TL );
  PLA_Obj_free( &A_BR );
  PLA_Obj_free( &A_1 );
  PLA_Obj_free( &A_1_dup );
  PLA_Obj_free( &A_11_dup );
  PLA_Obj_free( &B_L );
  PLA_Obj_free( &B_R );
  PLA_Obj_free( &B_1 );
  PLA_Obj_free( &B_1_dup );
  PLA_Obj_free( &one );
  PLA_Obj_free( &zero );

  return PLA_SUCCESS;
}


