/*
VisAD system for interactive analysis and visualization of numerical
data.  Copyright (C) 1996 - 2002 Bill Hibbard, Curtis Rueden, Tom
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

import java.awt.Dimension;
import java.awt.Toolkit;

import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;

import java.io.IOException;

import java.net.URL;

import java.rmi.RemoteException;

import javax.swing.BoxLayout;
import javax.swing.JFrame;
import javax.swing.JPanel;

import visad.CellImpl;
import visad.Data;
import visad.DataReferenceImpl;
import visad.Display;
import visad.DisplayImpl;
import visad.DisplayRenderer;
import visad.FieldImpl;
import visad.GraphicsModeControl;
import visad.Integer1DSet;
import visad.Real;
import visad.ScalarMap;
import visad.ShapeControl;
import visad.Tuple;
import visad.VisADException;

import visad.data.amanda.F2000Form;

import visad.java3d.DisplayImplJ3D;

import visad.util.LabeledColorWidget;
import visad.util.VisADSlider;

class DisplayFrame
  extends WindowAdapter
{
  private Display display;

  DisplayFrame(String title, Display display, JPanel panel,
               int width, int height)
    throws VisADException, RemoteException
  {
    this.display = display;

    JFrame frame = new JFrame(title);

    frame.addWindowListener(this);

    frame.getContentPane().add(panel);

    frame.setSize(width, height);

    Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
    frame.setLocation((screenSize.width - width)/2,
                      (screenSize.height - height)/2);

    frame.setVisible(true);
  }

  public void windowClosing(WindowEvent evt)
  {
    try { display.destroy(); } catch (Exception e) { }
    System.exit(0);
  }
}

class DisplayMaps
{
  ScalarMap trackmap, shapemap, letmap;

  DisplayMaps(F2000Form form, Display display)
    throws RemoteException, VisADException
  {
    // compute x, y and z ranges with unity aspect ratios
    final double xrange = form.getXMax() - form.getXMin();
    final double yrange = form.getYMax() - form.getYMin();
    final double zrange = form.getZMax() - form.getZMin();
    final double half_range = -0.5 * Math.max(xrange, Math.max(yrange, zrange));
    final double xmid = 0.5 * (form.getXMax() + form.getXMin());
    final double ymid = 0.5 * (form.getYMax() + form.getYMin());
    final double zmid = 0.5 * (form.getZMax() + form.getZMin());
    final double xmin = xmid - half_range;
    final double xmax = xmid + half_range;
    final double ymin = ymid - half_range;
    final double ymax = ymid + half_range;
    final double zmin = zmid - half_range;
    final double zmax = zmid + half_range;

    ScalarMap xmap = new ScalarMap(form.getX(), Display.XAxis);
    display.addMap(xmap);
    xmap.setRange(xmin, xmax);

    ScalarMap ymap = new ScalarMap(form.getY(), Display.YAxis);
    display.addMap(ymap);
    ymap.setRange(ymin, ymax);

    ScalarMap zmap = new ScalarMap(form.getZ(), Display.ZAxis);
    display.addMap(zmap);
    zmap.setRange(zmin, zmax);

    // ScalarMap eventmap = new ScalarMap(form.getEventIndex(),
    //                                    Display.SelectValue);
    // display.addMap(eventmap);

    this.trackmap = new ScalarMap(form.getTrackIndex(), Display.SelectValue);
    display.addMap(this.trackmap);

    // ScalarMap energymap = new ScalarMap(energy, Display.RGB);
    // display.addMap(energymap);

    this.shapemap = new ScalarMap(form.getAmplitude(), Display.Shape);
    display.addMap(this.shapemap);

    ScalarMap shape_scalemap = new ScalarMap(form.getAmplitude(),
                                             Display.ShapeScale);
    display.addMap(shape_scalemap);
    shape_scalemap.setRange(-20.0, 50.0);

    this.letmap = new ScalarMap(form.getLet(), Display.RGB);
    display.addMap(this.letmap);
  }
}

/** run 'java F2000App in_file' to display data.<br>
 *  try 'java F2000App 100events.r'
 */
public class F2000App
{
  public static void main(String args[])
         throws VisADException, RemoteException, IOException
  {
    if (args == null || args.length != 1) {
      System.out.println("to test read an F2000 file, run:");
      System.out.println("  'java visad.amanda.F2000Form in_file'");
      System.exit(1);
      return;
    }

    F2000Form form = new F2000Form();
    Data temp = null;
    if (args[0].startsWith("http://")) {
      // with "ftp://" this throws "sun.net.ftp.FtpProtocolException: RETR ..."
      URL url = new URL(args[0]);
      temp = form.open(url);
    }
    else {
      temp = form.open(args[0]);
    }

    final FieldImpl amanda = (FieldImpl) ((Tuple) temp).getComponent(0);
    final FieldImpl modules = (FieldImpl) ((Tuple) temp).getComponent(1);

    DisplayImpl display = new DisplayImplJ3D("amanda");

    DisplayMaps maps = new DisplayMaps(form, display);

    // GraphicsModeControl mode = display.getGraphicsModeControl();
    // mode.setScaleEnable(true);
    DisplayRenderer displayRenderer = display.getDisplayRenderer();
    displayRenderer.setBoxOn(false);

    ShapeControl scontrol = (ShapeControl) maps.shapemap.getControl();
    scontrol.setShapeSet(new Integer1DSet(form.getAmplitude(), 1));
    scontrol.setShapes(form.getCubeArray());

    final int nevents = amanda.getLength();

    // fixes track display?
    // SelectValue bug?
    // amanda = ((FieldImpl) amanda).getSample(99);

    final DataReferenceImpl amanda_ref = new DataReferenceImpl("amanda");
    // amanda_ref.setData(amanda);
    display.addReference(amanda_ref);

    final DataReferenceImpl modules_ref = new DataReferenceImpl("modules");
    modules_ref.setData(modules);
    display.addReference(modules_ref);

System.out.println("amanda MathType\n" + amanda.getType());
// visad.jmet.DumpType.dumpDataType(amanda, System.out);

    final DataReferenceImpl event_ref = new DataReferenceImpl("event");

    CellImpl cell = new CellImpl() {
      public void doAction() throws VisADException, RemoteException {
        Real r = (Real )event_ref.getData();
        if (r != null) {
          int index = (int )r.getValue();
          if (index < 0) {
            index = 0;
          } else if (index > nevents) {
            index = nevents;
          }
          amanda_ref.setData(amanda.getSample(index));
        }
      }
    };
    // link cell to hour_ref to trigger doAction whenever
    // 'hour' value changes
    cell.addReference(event_ref);

/*
    LabeledColorWidget energy_widget =
      new LabeledColorWidget(energymap);
    widget_panel.add(energy_widget);
    energy_widget.setMaximumSize(new Dimension(400, 250));
*/
    LabeledColorWidget let_widget =
      new LabeledColorWidget(maps.letmap);
    let_widget.setMaximumSize(new Dimension(400, 250));

    VisADSlider event_slider = new VisADSlider("event", 0, nevents, 0, 1.0,
                                               event_ref,
                                               form.getEventIndex());

    JPanel widget_panel = new JPanel();
    widget_panel.setLayout(new BoxLayout(widget_panel, BoxLayout.Y_AXIS));
    widget_panel.setMaximumSize(new Dimension(400, 600));

    widget_panel.add(let_widget);
    // widget_panel.add(new VisADSlider(eventmap));
    widget_panel.add(new VisADSlider(maps.trackmap));
    widget_panel.add(event_slider);

    JPanel display_panel = (JPanel) display.getComponent();
    Dimension dim = new Dimension(800, 800);
    display_panel.setPreferredSize(dim);
    display_panel.setMinimumSize(dim);

    // create JPanel in frame
    JPanel panel = new JPanel();
    panel.setLayout(new BoxLayout(panel, BoxLayout.X_AXIS));

    panel.add(widget_panel);
    panel.add(display_panel);

    new DisplayFrame("VisAD AMANDA Viewer", display, panel, 1200, 800);
  }
}
