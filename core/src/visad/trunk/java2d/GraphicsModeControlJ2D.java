//
// GraphicsModeControlJ2D.java
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

package visad.java2d;
 
import visad.*;

import java.rmi.*;

/**
   GraphicsModeControlJ2D is the VisAD class for controlling various
   mode settings for rendering.<P>

   A GraphicsModeControlJ2D is not linked to any DisplayRealType or
   ScalarMap.  It is linked to a DisplayImpl.<P>
*/
public class GraphicsModeControlJ2D extends GraphicsModeControl {

  private float lineWidth; // for LineAttributes; >= 1.0
  private float pointSize; // for PointAttributes; >= 1.0
  private boolean pointMode; // true => points in place of lines and surfaces
  private boolean textureEnable; // true => allow use of texture mapping
  private boolean scaleEnable; // true => display X, Y and Z scales

  private int transparencyMode = 0;
  private int projectionPolicy = 0;
  private int polygonMode = 0;

  private boolean missingTransparent = false;
  private int curvedSize = 10;

  public GraphicsModeControlJ2D(DisplayImpl d) {
    super(d);
    lineWidth = 1.0f;
    pointSize = 1.0f;
    pointMode = false;
    textureEnable = true;
    scaleEnable = false;
  }
 
  public boolean getMode2D() {
    return getDisplayRenderer().getMode2D();
  }

  public float getLineWidth() {
    return lineWidth;
  }

  public void setLineWidth(float width)
         throws VisADException, RemoteException {
    if (width < 1.0f) {
      throw new DisplayException("GraphicsModeControlJ2D." +
                                 "setLineWidth: width < 1.0");
    }
    lineWidth = width;
    changeControl(true);
    getDisplay().reDisplayAll();
  }

  public void setLineWidth(float width, boolean dummy) {
    if (width >= 1.0f) {
      lineWidth = width;
    }
  }

  public float getPointSize() {
    return pointSize;
  }

  public void setPointSize(float size)
         throws VisADException, RemoteException {
    if (size < 1.0f) {
      throw new DisplayException("GraphicsModeControlJ2D." +
                                 "setPointSize: size < 1.0");
    }
    pointSize = size;
    changeControl(true);
    getDisplay().reDisplayAll();
  }

  public void setPointSize(float size, boolean dummy) {
    if (size >= 1.0f) {
      pointSize = size;
    }
  }

  public boolean getPointMode() {
    return pointMode;
  }

  public void setPointMode(boolean mode)
         throws VisADException, RemoteException {
    pointMode = mode;
    changeControl(true);
    getDisplay().reDisplayAll();
  }

  public void setTextureEnable(boolean enable)
         throws VisADException, RemoteException {
    textureEnable = enable;
    changeControl(true);
    getDisplay().reDisplayAll();
  }

  public boolean getTextureEnable() {
    return textureEnable;
  }

  public void setScaleEnable(boolean enable)
         throws VisADException, RemoteException {
    scaleEnable = enable;
    getDisplayRenderer().setScaleOn(enable);
    changeControl(true);
  }
 
  public boolean getScaleEnable() {
    return scaleEnable;
  }

  public int getTransparencyMode() {
    return transparencyMode;
  }

  public void setTransparencyMode(int mode)
         throws VisADException, RemoteException {
    if (mode == 0) {
      transparencyMode = mode;
    }
    else {
      throw new DisplayException("GraphicsModeControlJ2D." +
                                 "setTransparencyMode: bad mode");
    }
  }

  public void setProjectionPolicy(int policy)
         throws VisADException, RemoteException {
    if (policy == 0) {
      projectionPolicy = policy;
    }
    else {
      throw new DisplayException("GraphicsModeControlJ2D." +
                                 "setProjectionPolicy: bad policy");
    }
  }

  public int getProjectionPolicy() {
    return projectionPolicy;
  }

  public void setPolygonMode(int mode)
         throws VisADException, RemoteException {
    if (mode == 0) {
      polygonMode = mode;
    }
    else {
      throw new DisplayException("GraphicsModeControlJ2D." +
                                 "setPolygonMode: bad mode");
    }
  }

  public int getPolygonMode() {
    return polygonMode;
  }

  public boolean getMissingTransparent() {
    return missingTransparent;
  }

  public void setMissingTransparent(boolean missing)
         throws VisADException, RemoteException {
    if (!missing) {
      missingTransparent = missing;
    }
    else {
      throw new DisplayException("GraphicsModeControlJ2D." +
                                 "setMissingTransparent: must be false");
    }
  }

  public int getCurvedSize() {
    return curvedSize;
  }

  public void setCurvedSize(int curved_size) {
    curvedSize = curved_size;
  }

  public Object clone() {
    GraphicsModeControlJ2D mode =
      new GraphicsModeControlJ2D(getDisplay());
    mode.lineWidth = lineWidth;
    mode.pointSize = pointSize;
    mode.pointMode = pointMode;
    mode.textureEnable = textureEnable;
    mode.scaleEnable = scaleEnable;
    mode.transparencyMode = transparencyMode;
    mode.projectionPolicy = projectionPolicy;
    mode.missingTransparent = missingTransparent;
    mode.polygonMode = polygonMode;
    mode.curvedSize = curvedSize;
    return mode;
  }

}

