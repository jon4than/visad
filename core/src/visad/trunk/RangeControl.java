//
// RangeControl.java
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

package visad;

import java.rmi.*;

/**
   RangeControl is the VisAD class for controlling SelectRange display scalars.<P>
*/
public class RangeControl extends Control {

  private boolean initialized = false;
  private float RangeLow = Float.NaN;
  private float RangeHi = Float.NaN;

  public RangeControl(DisplayImpl d) {
    super(d);
  }

  /** initialize the range of selected values as (range[0], range[1]) */
  public void initializeRange(float[] range)
    throws VisADException, RemoteException
  {
    changeRange(range, initialized);
  }

  /** set the range of selected values as (range[0], range[1]) */
  public void setRange(float[] range)
         throws VisADException, RemoteException {
    if (RangeLow != RangeLow || Math.abs(RangeLow - range[0]) > 0.0001 ||
        RangeHi != RangeHi || Math.abs(RangeHi - range[1]) > 0.0001)
    {
      changeRange(range, true);
    }
  }

  private void changeRange(float[] range, boolean notify)
    throws VisADException, RemoteException
  {
    RangeLow = range[0];
    RangeHi = range[1];
    initialized = (RangeLow == RangeLow && RangeHi == RangeHi);
    if (notify) {
      changeControl(true);
    }
  }

  /** return the range of selected values */
  public float[] getRange() {
    float[] range = new float[2];
    range[0] = RangeLow;
    range[1] = RangeHi;
    return range;
  }

  /** copy the state of a remote control to this control */
  public void syncControl(Control rmt)
    throws VisADException
  {
    if (rmt == null) {
      throw new VisADException("Cannot synchronize " + getClass().getName() +
                               " with null Control object");
    }

    if (!(rmt instanceof RangeControl)) {
      throw new VisADException("Cannot synchronize " + getClass().getName() +
                               " with " + rmt.getClass().getName());
    }

    RangeControl rc = (RangeControl )rmt;

    boolean changed = false;

    if (Math.abs(RangeLow - rc.RangeLow) > 0.0001) {
      changed = true;
      RangeLow = rc.RangeLow;
    }
    if (Math.abs(RangeHi - rc.RangeHi) > 0.0001) {
      changed = true;
      RangeHi = rc.RangeHi;
    }
    initialized = (RangeLow == RangeLow && RangeHi == RangeHi);

    if (changed) {
      try {
        changeControl(true);
      } catch (RemoteException re) {
        throw new VisADException("Could not indicate that control" +
                                 " changed: " + re.getMessage());
      }
    }
  }

  public boolean equals(Object o)
  {
    if (o == null || !(o instanceof RangeControl)) {
      return false;
    }

    RangeControl rc = (RangeControl )o;

    if (Math.abs(RangeLow - rc.RangeLow) > 0.0001) {
      return false;
    }
    if (Math.abs(RangeHi - rc.RangeHi) > 0.0001) {
      return false;
    }

    return true;
  }

  public String toString()
  {
    return "RangeControl[" + RangeLow + "," + RangeHi + "]";
  }
}

