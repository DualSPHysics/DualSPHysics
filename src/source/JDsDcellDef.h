//HEAD_DSPH
/*
 <DUALSPHYSICS>  Copyright (c) 2021 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

 EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
 School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

 This file is part of DualSPHysics. 

 DualSPHysics is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License 
 as published by the Free Software Foundation; either version 2.1 of the License, or (at your option) any later version.
 
 DualSPHysics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details. 

 You should have received a copy of the GNU Lesser General Public License along with DualSPHysics. If not, see <http://www.gnu.org/licenses/>. 
*/

/// \file JDsDcellDef.h \brief Includes definitions for cell codificantion as unsigned (32 bits) value.

#ifndef _JDsDcellDef_
#define _JDsDcellDef_

//##############################################################################
// Codification of local cells according to the position (using 1+31 bits).
//##############################################################################
#define DCEL_CodeMapOut   0xffffffff  //-Special value for cells outside the limits.
#define DCEL_CodeSpecial  0x80000000  //-Mask value for special codes.
#define DCEL_CodeDomLeft  0x80000000  //-Special value for cells to the left of the limits.
#define DCEL_CodeDomRight 0x80000001  //-Special value for cells to the right of the limits.
//#define DCEL_GetCode(sx,sy,sz) (((sx)<<25)|(sy<<20)|(sz<<15)|((sy+sz)<<10)|((sx+sz)<<5)|(sx+sy))  //-Clave de codificacion (orden de valores: sx,sy,sz,sy+sz,sx+sz,sx+sy). | Encryption key (order of values: sx,sy,sz,sy+sz,sz+sx,sx+sy) (using 32bits).
#define DCEL_GetCode(sx,sy,sz) (((sx+1)<<25)|(sy<<20)|(sz<<15)|((sy+sz)<<10)|((sx+1+sz)<<5)|(sx+1+sy))  //-Clave de codificacion (orden de valores: sx,sy,sz,sy+sz,sx+sz,sx+sy). | Encryption key (order of values: sx,sy,sz,sy+sz,sz+sx,sx+sy) (extra 1 bit in sx).
#define DCEL_GetSx(dcc) (dcc>>25)       //-Numero de bits para coordenada X de celda. | Number of bits for X coordinate cell.
#define DCEL_GetSy(dcc) ((dcc>>20)&31)  //-Numero de bits para coordenada Y de celda. | Number of bits for Y coordinate cell.
#define DCEL_GetSz(dcc) ((dcc>>15)&31)  //-Numero de bits para coordenada Z de celda. | Number of bits for Z coordinate cell.
#define DCEL_Cellx(dcc,cel) ((*((unsigned*)&cel))>>((dcc>>10)&31))              //-Coordenada X de celda. | X coordinate of the cell.
#define DCEL_Celly(dcc,cel) (((*((unsigned*)&cel))<<(dcc>>25))>>((dcc>>5)&31))  //-Coordenada Y de celda. | Y coordinate of the cell.
#define DCEL_Cellz(dcc,cel) (((*((unsigned*)&cel))<<(dcc&31))>>(dcc&31))        //-Coordenada Z de celda. | Z coordinate of the cell.
#define DCEL_Cell(dcc,cx,cy,cz) ((cx<<((dcc>>10)&31))|(cy<<((dcc>>15)&31))|cz)  //-Valor de celda para cx, cy y cz. | Cell value for cx,cy and cz.
//#define DCEL_MaxCellx(dcc) ((0xffffffff>>((dcc>>10)&31)))            //-Coordenada X de celda maxima. | Maximum X coordinate of the cell (using 32bits).
#define DCEL_MaxCellx(dcc) ((0xffffffff>>((dcc>>10)&31))>>1)           //-Coordenada X de celda maxima. | Maximum X coordinate of the cell (ignores 1 bit).
#define DCEL_MaxCelly(dcc) ((0xffffffff<<(dcc>>25))>>((dcc>>5)&31))    //-Coordenada Y de celda maxima. | Maximum Y coordinate of the cell.
#define DCEL_MaxCellz(dcc) ((0xffffffff<<(dcc&31))>>(dcc&31))          //-Coordenada Z de celda maxima. | Maximum Z coordinate of the cell.


#endif


