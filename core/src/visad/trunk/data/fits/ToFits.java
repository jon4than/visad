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

package visad.data.fits;

import java.io.IOException;

import java.rmi.RemoteException;

import java.net.URL;

import visad.Data;
import visad.VisADException;

import visad.data.DefaultFamily;
import visad.data.DataNode;

public class ToFits
{
  public static void main(String args[])
	throws VisADException, RemoteException, IOException
  {
    DefaultFamily dflt = new DefaultFamily("default");

    if (args.length == 0) {
      args = new String[1];
      args[0] = "testdata/sseclogo.fits";
    }

    for (int i = 0; i < args.length; i++) {
      Data data = dflt.open(args[i]);

      try {
	System.out.println("ToFits " + args[i] + ": " + data.getType());
      } catch (Exception e) {
	System.err.println(args[i] + " print threw " + e.getMessage());
	e.printStackTrace(System.err);
	data = null;
	continue;
      }

      String name = "foo" + i;
      FitsForm form = new FitsForm();
      form.save(name, data, true);
      System.out.println("Wrote " + name);
    }

    System.exit(0);
  }
}
