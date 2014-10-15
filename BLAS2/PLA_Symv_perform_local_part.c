/*
   PLAPACK Release 3.0
   
   Copyright (C) 2003
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Symv_perform_local_part( int uplo, 
				 PLA_Obj A, 
				 PLA_Obj x_dup_onto_rows, 
				 PLA_Obj x_dup_onto_cols, 
				 PLA_Obj y_dup_onto_rows, 
				 PLA_Obj y_dup_onto_cols ) 
{
  int 
    myrow, mycol, align_row, align_col,
    owner_top, owner_bottom, owner_left, owner_right,
    size, size_top,  size_bottom,  size_left,  size_right;

  PLA_Obj
    A_cur = NULL, A_11 = NULL, A_21 = NULL, A_12 = NULL,
    x_dup_onto_rows_cur = NULL, x_dup_onto_rows_1 = NULL,
    x_dup_onto_cols_cur = NULL, x_dup_onto_cols_1 = NULL,
    y_dup_onto_rows_cur = NULL, y_dup_onto_rows_1 = NULL,
    y_dup_onto_cols_cur = NULL, y_dup_onto_cols_1 = NULL,
    one = NULL;

  PLA_Template
    templ;

  PLA_Obj_set_to_zero( y_dup_onto_rows );
  PLA_Obj_set_to_zero( y_dup_onto_cols );

  PLA_Obj_template ( A, &templ );
  PLA_Temp_comm_col_rank( templ, &myrow );
  PLA_Temp_comm_row_rank( templ, &mycol );

  PLA_Obj_global_align_row( A, &align_row );
  PLA_Obj_global_align_col( A, &align_col );

  if ( align_row != align_col )
    PLA_Abort( "PLA_Symv_perform_local_part: only works when matrix is aligned on diagonal", __LINE__, __FILE__ );


  PLA_Create_constants_conf_to( A, NULL, NULL, &one );

  PLA_Obj_view_all( A, &A_cur );
  PLA_Obj_view_all( x_dup_onto_rows, &x_dup_onto_rows_cur );
  PLA_Obj_view_all( x_dup_onto_cols, &x_dup_onto_cols_cur );
  PLA_Obj_view_all( y_dup_onto_rows, &y_dup_onto_rows_cur );
  PLA_Obj_view_all( y_dup_onto_cols, &y_dup_onto_cols_cur );
    
  while ( TRUE ){
    PLA_Obj_split_size( A_cur, PLA_SIDE_TOP,  &size_top,  &owner_top );
    PLA_Obj_split_size( A_cur, PLA_SIDE_LEFT, &size_left, &owner_left );
    if ( 0 == ( size = min( size_top, size_left ) ) ) break;

    PLA_Obj_split_4( A_cur, size, size, &A_11, &A_12,
		                        &A_21, &A_cur );

    PLA_Obj_vert_split_2( x_dup_onto_rows_cur, size, 
			  &x_dup_onto_rows_1, &x_dup_onto_rows_cur );
    PLA_Obj_horz_split_2( x_dup_onto_cols_cur, size, &x_dup_onto_cols_1, 
			                             &x_dup_onto_cols_cur );
    PLA_Obj_vert_split_2( y_dup_onto_rows_cur, size, 
			  &y_dup_onto_rows_1, &y_dup_onto_rows_cur );
    PLA_Obj_horz_split_2( y_dup_onto_cols_cur, size, &y_dup_onto_cols_1, 
			                             &y_dup_onto_cols_cur );

    if ( mycol == owner_left ){
      if ( myrow == owner_top )
	PLA_Local_symv( uplo, one, A_11, x_dup_onto_rows_1, 
			one, y_dup_onto_cols_1 ); 

      if ( uplo == PLA_LOWER_TRIANGULAR ){
	PLA_Local_gemv( PLA_NO_TRANS, one, A_21, x_dup_onto_rows_1, 
			one, y_dup_onto_cols_cur );

	PLA_Local_gemv( PLA_TRANS, one, A_21, x_dup_onto_cols_cur, 
			one, y_dup_onto_rows_1 );
      }
      else{
	PLA_Local_gemv( PLA_NO_TRANS, one, A_12, x_dup_onto_rows_cur, 
			one, y_dup_onto_cols_1 );

	PLA_Local_gemv( PLA_TRANS, one, A_12, x_dup_onto_cols_1, 
			one, y_dup_onto_cols_cur );
      }
    }
  }

  PLA_Obj_free( &A_cur );
  PLA_Obj_free( &A_11 );
  PLA_Obj_free( &A_21 );
  PLA_Obj_free( &A_12 );
  PLA_Obj_free( &x_dup_onto_rows_cur );
  PLA_Obj_free( &x_dup_onto_rows_1 );
  PLA_Obj_free( &x_dup_onto_cols_cur );
  PLA_Obj_free( &x_dup_onto_cols_1 );
  PLA_Obj_free( &y_dup_onto_rows_cur );
  PLA_Obj_free( &y_dup_onto_rows_1 );
  PLA_Obj_free( &y_dup_onto_cols_cur );
  PLA_Obj_free( &y_dup_onto_cols_1 );
  PLA_Obj_free( &one );

  return PLA_SUCCESS;
}
	    




