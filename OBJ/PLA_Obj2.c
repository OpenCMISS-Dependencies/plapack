/*
   PLAPACK Release 3.1
   
   Copyright (C) 2001
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

/* ***************************************************************************

   PLA_Part_2x2( )

 *************************************************************************** */

int PLA_Part_2x2( PLA_Obj obj, 
                  PLA_Obj *A11, PLA_Obj *A12,
		  PLA_Obj *A21, PLA_Obj *A22, 
		  int    mb,   int nb,      int quadrant )
{
  int
    m, n;

  PLA_Obj_global_length( obj, &m );
  PLA_Obj_global_width ( obj, &n );

  switch ( quadrant ){
  case PLA_TL:
    PLA_Obj_split_4( obj,   mb,   nb,  A11, A12,
		                       A21, A22 );
    break;
  case PLA_TR:
    PLA_Obj_split_4( obj,   mb, n-nb,  A11, A12,
		                       A21, A22 );
    break;
  case PLA_BL:
    PLA_Obj_split_4( obj,  m-mb,  nb,  A11, A12,
		                       A21, A22 );
    break;
  case PLA_BR:
    PLA_Obj_split_4( obj, m-mb, n-nb,  A11, A12,
		                       A21, A22 );
    break;
  }

  return PLA_SUCCESS;
}

/* ***************************************************************************

   PLA_Repart_2x2_to_3x3( )

 *************************************************************************** */

int PLA_Repart_2x2_to_3x3( 
   PLA_Obj ATL, PLA_Obj ATR,  PLA_Obj *A00, PLA_Obj *A01, PLA_Obj *A02,
	 		      PLA_Obj *A10, PLA_Obj *A11, PLA_Obj *A12,
   PLA_Obj ABL, PLA_Obj ABR,  PLA_Obj *A20, PLA_Obj *A21, PLA_Obj *A22,
   int    mb,   int nb,      int quadrant )
{
  int
    m, n;

  switch ( quadrant ){
  case PLA_TL:
    PLA_Obj_global_length( ATL, &m );
    PLA_Obj_global_width ( ATL, &n );

    PLA_Obj_split_4( ATL, m-mb, n-nb,    A00, A01,
                       		         A10, A11 );

    PLA_Obj_horz_split_2( ATR, m-mb,                 A02, 
                                                     A12 );    

    PLA_Obj_vert_split_2( ABL, n-nb,    A20, A21 );

    PLA_Obj_view_all ( ABR,                          A22 );    

    break;

  case PLA_TR:
    PLA_Obj_global_length( ATR, &m );

    PLA_Obj_split_4( ATR, m-mb, nb,                A01, A02,
                       		                   A11, A12 );

    PLA_Obj_horz_split_2( ATL, m-mb,        A02, 
                                            A12 );    

    PLA_Obj_vert_split_2( ABR, nb,                 A21, A22 );

    PLA_Obj_view_all ( ABL,                 A20 );    

    break;

  case PLA_BL:
    PLA_Obj_global_width ( ABL, &n );

    PLA_Obj_vert_split_2( ATL, n-nb,    A00, A01 );

    PLA_Obj_view_all ( ATR,                          A02 );    

    PLA_Obj_split_4( ABL,   mb, n-nb,    A10, A11,
                       		         A20, A21 );

    PLA_Obj_horz_split_2( ABR,   mb,                 A12, 
                                                     A22 );    

    break;

  case PLA_BR:
    PLA_Obj_split_4( ABR, mb, nb,             A11, A12,
                       		              A21, A22 );

    PLA_Obj_horz_split_2( ABL, mb,      A10, 
                                        A20 );    

    PLA_Obj_vert_split_2( ATR, nb,            A01, A02 );

    PLA_Obj_view_all ( ATL,             A00 );    

    break;
  }
  
  return PLA_SUCCESS;
}

/* ***************************************************************************

   PLA_Cont_with_3x3_to_2x2( )

 *************************************************************************** */

int PLA_Cont_with_3x3_to_2x2( PLA_Obj *ATL, PLA_Obj *ATR,  PLA_Obj A00, PLA_Obj A01, PLA_Obj A02,
 		                                           PLA_Obj A10, PLA_Obj A11, PLA_Obj A12,
                               PLA_Obj *ABL, PLA_Obj *ABR, PLA_Obj A20, PLA_Obj A21, PLA_Obj A22,
			       int quadrant )
{
  int mb, nb;

  PLA_Obj_global_length( A11, &mb );
  PLA_Obj_global_width ( A11, &nb );

  PLA_Obj_view_all( A00, ATL );     PLA_Obj_view_all( A02, ATR );
  PLA_Obj_view_all( A20, ABL );     PLA_Obj_view_all( A22, ABR );

  switch ( quadrant ){
  case PLA_TL:
    PLA_Obj_view_shift( *ATL,       0,
		 	      0,        nb,
			          mb );

    PLA_Obj_view_shift( *ATR,       0,
		 	      0,        0,
			          mb );

    PLA_Obj_view_shift( *ABL,       0,
		 	      0,        nb,
			           0 );

    break;

  case PLA_BR:
    PLA_Obj_view_shift( *ATR,       0,
		 	    -nb,         0,
			           0 );

    PLA_Obj_view_shift( *ABL,      -mb,
		 	      0,        0,
			          0 );

    PLA_Obj_view_shift( *ABR,      -mb,
                             -nb,        0,
			           0 );

    break;

  default:
    printf("not yet implemented 1\n");
  }
}

/* ***************************************************************************

   PLA_Part_2x1( )

 *************************************************************************** */

int PLA_Part_2x1 ( PLA_Obj obj,   PLA_Obj *A1, 
                                  PLA_Obj *A2,
		   int    mb,  int side )
{
  int 
    m;

  PLA_Obj_global_length( obj, &m );

  if ( side == PLA_TOP )
    PLA_Obj_horz_split_2( obj, mb,   A1,
			             A2 );
  else
    PLA_Obj_horz_split_2( obj, m-mb, A1,
			             A2 );
  
  return PLA_SUCCESS;
}

/* ***************************************************************************

   PLAR_epart_2x1_to_3x1( )

 *************************************************************************** */

int  PLA_Repart_2x1_to_3x1( PLA_Obj AT,   PLA_Obj *A0,
				          PLA_Obj *A1,
			    PLA_Obj AB,   PLA_Obj *A2,
			    int    mb,    int side )
{
  int
    m;

  if ( side == PLA_TOP ){
    PLA_Obj_global_length( AT, &m );

    PLA_Obj_horz_split_2( AT, m-mb, A0,
                                    A1 );

    PLA_Obj_view_all( AB, A2 );
  }
  else {
    PLA_Obj_view_all( AT, A0 );

    PLA_Obj_horz_split_2( AB, mb,   A1,
                                    A2 );
  }

  return PLA_SUCCESS;
}

/* ***************************************************************************

   PLA_Cont_with_3x1_to_2x1( )

 *************************************************************************** */

int PLA_Cont_with_3x1_to_2x1( PLA_Obj *AT,   PLA_Obj A0,
                                             PLA_Obj A1,
			      PLA_Obj *AB,   PLA_Obj A2,
			      int side )
{
  int 
    mb;

  PLA_Obj_global_length( A1, &mb );

  PLA_Obj_view_all( A0, AT );
  PLA_Obj_view_all( A2, AB );

  if ( side == PLA_TOP ){
    PLA_Obj_view_shift( *AT,        0,
		 	      0,        0,
			          mb );
  }
  else {
    PLA_Obj_view_shift( *AB,      -mb,
		 	      0,        0,
			          0 );
  }

  return PLA_SUCCESS;
}

/* ***************************************************************************

   PLA_Part_1x2( )

 *************************************************************************** */

int PLA_Part_1x2( PLA_Obj obj, 
		  PLA_Obj *A1, PLA_Obj *A2,
		  int    nb,  int side )
{
  int 
    n;

  PLA_Obj_global_width( obj, &n );

  if ( side == PLA_LEFT )
    PLA_Obj_vert_split_2( obj, nb,   A1,  A2 );
  else
    PLA_Obj_vert_split_2( obj, n-nb, A1,  A2 );
  
  return PLA_SUCCESS;
}

/* ***************************************************************************

   PLA_Repart_1x2_to_1x3( )

 *************************************************************************** */

int  PLA_Repart_1x2_to_1x3( PLA_Obj  AL,              PLA_Obj  AR,
			    PLA_Obj *A0, PLA_Obj *A1, PLA_Obj *A2,
			    int    nb,    int side )
{
  int
    n;

  if ( side == PLA_LEFT ){
    PLA_Obj_global_width( AL, &n );

    PLA_Obj_vert_split_2( AL, n-nb, A0, A1 );

    PLA_Obj_view_all( AR, A2 );
  }
  else {
    PLA_Obj_view_all( AL, A0 );

    PLA_Obj_vert_split_2( AR, nb,   A1, A2 );
  }

  return PLA_SUCCESS;
}

/* ***************************************************************************

   PLA_Cont_with_1x3_to_1x2( )

 *************************************************************************** */

int  PLA_Cont_with_1x3_to_1x2( PLA_Obj *AL,           PLA_Obj *AR,
			       PLA_Obj A0, PLA_Obj A1, PLA_Obj A2,
			       int side )
{
  int 
    nb;

  PLA_Obj_global_width( A1, &nb );

  PLA_Obj_view_all( A0, AL );
  PLA_Obj_view_all( A2, AR );

  if ( side == PLA_LEFT ){
    PLA_Obj_view_shift( *AL,        0,
		 	      0,       nb,
			           0 );
  }
  else {
    PLA_Obj_view_shift( *AR,       0,
		 	   -nb,        0,
			          0 );
  }

  return PLA_SUCCESS;
}

