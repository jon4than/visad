
//
// ShadowBarbTupleTypeJ2D.java
//

/*
VisAD system for interactive analysis and visualization of numerical
data.  Copyright (C) 1996 - 1998 Bill Hibbard, Curtis Rueden, Tom
Rink and Dave Glowacki.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 1, or (at your option)
any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License in file NOTICE for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/

package visad.bom;

import visad.*;
import visad.java2d.*;

import java.util.*;
import java.rmi.*;

/**
   The ShadowBarbTupleTypeJ2D class shadows the TupleType class for
   BarbRendererJ2D, within a DataDisplayLink, under Java2D.<P>
*/
public class ShadowBarbTupleTypeJ2D extends ShadowTupleTypeJ2D {

  public ShadowBarbTupleTypeJ2D(MathType t, DataDisplayLink link,
                                ShadowType parent)
         throws VisADException, RemoteException {
    super(t, link, parent);
  }

  public VisADGeometryArray[] makeFlow(int which, float[][] flow_values,
                float flowScale, float[][] spatial_values,
                byte[][] color_values, boolean[][] range_select)
         throws VisADException {

    DataRenderer renderer = getLink().getRenderer();
    boolean direct = renderer.getIsDirectManipulation();
    if (direct && renderer instanceof BarbManipulationRendererJ2D) {
      return ShadowBarbRealTupleTypeJ2D.staticMakeFlow(getDisplay(), which,
                 flow_values, flowScale, spatial_values, color_values,
                 range_select, renderer, true);
    }
    else {
      return ShadowBarbRealTupleTypeJ2D.staticMakeFlow(getDisplay(), which,
                 flow_values, flowScale, spatial_values, color_values,
                 range_select, renderer, false);
    }
  }

}

