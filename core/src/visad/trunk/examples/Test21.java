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

import java.awt.Component;

import java.rmi.RemoteException;

import visad.*;

import visad.java3d.DisplayImplJ3D;

import visad.util.SelectRangeWidget;

public class Test21
	extends UISkeleton
{
  public Test21() { }

  public Test21(String args[])
	throws VisADException, RemoteException
  {
    super(args);
  }

  DisplayImpl[] setupData()
	throws VisADException, RemoteException
  {
    RealType[] types = {RealType.Latitude, RealType.Longitude};
    RealTupleType earth_location = new RealTupleType(types);
    RealType vis_radiance = new RealType("vis_radiance", null, null);
    RealType ir_radiance = new RealType("ir_radiance", null, null);
    RealType[] types2 = {vis_radiance, ir_radiance};
    RealTupleType radiance = new RealTupleType(types2);
    FunctionType image_tuple = new FunctionType(earth_location, radiance);

    int size = 64;
    FlatField imaget1 = FlatField.makeField(image_tuple, size, false);

    DisplayImpl display1;
    display1 = new DisplayImplJ3D("display1", DisplayImplJ3D.APPLETFRAME);
    display1.addMap(new ScalarMap(RealType.Latitude, Display.YAxis));
    display1.addMap(new ScalarMap(RealType.Longitude, Display.XAxis));
    display1.addMap(new ScalarMap(vis_radiance, Display.ZAxis));
    display1.addMap(new ScalarMap(vis_radiance, Display.Green));
    display1.addMap(new ConstantMap(0.5, Display.Blue));
    display1.addMap(new ConstantMap(0.5, Display.Red));
    // display1.addMap(new ConstantMap(0.25, Display.Alpha));

    ScalarMap range1map = new ScalarMap(ir_radiance, Display.SelectRange);
    display1.addMap(range1map);

    GraphicsModeControl mode = display1.getGraphicsModeControl();
    mode.setPointSize(2.0f);
    mode.setPointMode(false);
    mode.setMissingTransparent(true);

    DataReferenceImpl ref_imaget1 = new DataReferenceImpl("ref_imaget1");
    ref_imaget1.setData(imaget1);
    display1.addReference(ref_imaget1, null);

    DisplayImpl[] dpys = new DisplayImpl[1];
    dpys[0] = display1;

    return dpys;
  }

  String getFrameTitle() { return "VisAD select range slider"; }

  Component getSpecialComponent(DisplayImpl[] dpys)
	throws VisADException, RemoteException
  {
    ScalarMap range1map = (ScalarMap )dpys[0].getMapVector().lastElement();
    return new SelectRangeWidget(range1map);
  }

  public String toString() { return ": SelectRange and SelectRangeWidget"; }

  public static void main(String args[])
	throws VisADException, RemoteException
  {
    Test21 t = new Test21(args);
  }
}
