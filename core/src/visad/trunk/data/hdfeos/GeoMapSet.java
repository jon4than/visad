//
// GeoMapSet.java
//

/*
VisAD system for interactive analysis and visualization of numerical
data.  Copyright (C) 1996 - 1999 Bill Hibbard, Curtis Rueden, Tom
Rink, Dave Glowacki, Steve Emmerson, Tom Whittaker, Don Murray, and
Tommy Jasmin.

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Library General Public
License as published by the Free Software Foundation; either
version 2 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Library General Public License for more details.

You should have received a copy of the GNU Library General Public
License along with this library; if not, write to the Free
Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
MA 02111-1307, USA
*/

package visad.data.hdfeos;

import java.lang.*;
import java.util.*;


public class GeoMapSet 
{
  Vector mapSet;
  
  public GeoMapSet() 
  {
    mapSet = new Vector(); 
  }

  public void add( GeoMap obj ) 
  {
    mapSet.addElement( obj );
  }

  public int getSize() 
  {
    int size = mapSet.size();
    return size;
  }

  public GeoMap getElement( int ii )  
  {
    if ( mapSet.size() == 0 ) 
    {
      return null;
    }
    else 
    {
      GeoMap obj = (GeoMap) mapSet.elementAt(ii);
      return obj;
    }
  }

  public GeoMap getGeoMap( NamedDimension obj ) 
  {
    String name = obj.getName();
    return getGeoMap( name );
  }

  public GeoMap getGeoMap( String name ) 
  {
    int size = this.getSize();

    if ( size == 0 ) 
    {
      return null;
    }
    else 
    {
      for ( int ii = 0; ii < size; ii++ ) 
      {
         GeoMap obj = (GeoMap) mapSet.elementAt(ii);

         if(( obj.toDim.equals( name ) ) || ( obj.fromDim.equals( name ) )) 
         {
           return obj;
         }
      }
      return null;
    }
  }
}
