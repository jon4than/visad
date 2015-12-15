//
// LAMBnav.java
//

/*
This source file is part of the edu.wisc.ssec.mcidas package and is
Copyright (C) 1998 - 2015 by Tom Whittaker, Tommy Jasmin, Tom Rink,
Don Murray, James Kelly, Bill Hibbard, Dave Glowacki, Curtis Rueden
and others.
 
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

package edu.wisc.ssec.mcidas;

/**
 * Navigation class for Lambert Conformal (LAMB) type nav. This code was 
 * modified from the original FORTRAN code (nvxlamb.dlm) on the McIDAS system. 
 * It only supports latitude/longitude to line/element transformations (LL) 
 * and vice/versa. Transform to 'XYZ' not implemented.
 * @see <A HREF="http://www.ssec.wisc.edu/mcidas/doc/prog_man.html">
 *      McIDAS Programmer's Manual</A>
 *
 * @author  Don Murray
 */
public final class LAMBnav extends AREAnav 
{

    int iwest;
    int ihem;
    double xrow;
    double xcol;
    double xlat1;
    double xlat2;
    double xspace;
    double xqlon;
    double xblat;
    double xfac;
    double xpole;

    /**
     * Set up for the real math work.  Must pass in the int array
     * of the LAMB nav 'codicil'.
     *
     * @param iparms  the nav block from the image file
     * @throws IllegalArgumentException
     *           if the nav block is not a LAMB type.
     */
    public LAMBnav (int[] iparms) 
        throws IllegalArgumentException
    {

        if (iparms[0] != LAMB ) 
            throw new IllegalArgumentException("Invalid navigation type" + 
                                                iparms[0]);
        xrow = iparms[1];
        xcol = iparms[2];
        int ipole = iparms[11];
        if (ipole == 0) ipole = 900000;
        ihem = 1;
        if (ipole < 0) ihem = -1;
        xpole = McIDASUtil.integerLatLonToDouble(ipole);

        xlat1 = McIDASUtil.integerLatLonToDouble(iparms[3]);
        xlat2 = McIDASUtil.integerLatLonToDouble(iparms[4]);
        xlat1 = (90. - ihem*xlat1)*DEGREES_TO_RADIANS;
        xlat2 = (90. - ihem*xlat2)*DEGREES_TO_RADIANS;

        xspace = iparms[5]/1000.;
        xqlon = McIDASUtil.integerLatLonToDouble(iparms[6]);
        double r = iparms[7]/1000.;

        iwest = iparms[10];
        if (iwest >= 0) iwest = 1;

        xfac = (Math.log(Math.sin(xlat1)) - Math.log(Math.sin(xlat2)))/
               (Math.log(Math.tan(.5*xlat1)) - Math.log(Math.tan(.5*xlat2)));
        xblat = r * Math.sin(xlat1)/
                (xspace*xfac*Math.pow(Math.tan(xlat1*.5), xfac));

    }

    /** converts from satellite coordinates to latitude/longitude
     *
     * @param  linele	  array of line/element pairs.  Where 
     *                     linele[indexLine][] is a 'line' and 
     *                     linele[indexEle][] is an element. These are in 
     *                     'file' coordinates (not "image" coordinates.)
     *
     * @return latlon[][]  array of lat/long pairs. Output array is 
     *                     latlon[indexLat][] of latitudes and 
     *                     latlon[indexLon][] of longitudes.
     *
     */
    public double[][] toLatLon(double[][] linele) 
    {

        double xldif;
        double xedif;
        double xlon;
        double xlat;
        double xrlon, radius;

        int number = linele[0].length;
        double[][] latlon = new double[2][number];

        // Convert array to Image coordinates for computations
        double[][] imglinele = areaCoordToImageCoord(linele);

        for (int point=0; point < number; point++) 
        {
            xldif = ihem*(imglinele[indexLine][point] - xrow)/xblat;
            xedif = -ihem*(imglinele[indexEle][point] - xcol)/xblat;
            xrlon = 0;
            if (!(xldif == 0 && xedif == 0)) xrlon = Math.atan2(xedif, xldif);
            xlon = ihem*xrlon/xfac/DEGREES_TO_RADIANS + xqlon;
            xlon = (xlon+900.)%360. - 180.0;
            radius = Math.sqrt(xldif*xldif + xedif*xedif);
            if (Math.abs(radius) < 1.e-10)
                xlat = ihem*90;
            else
                xlat = ihem*(90. - 2*Math.atan(
                          Math.exp(Math.log(radius)/xfac))/DEGREES_TO_RADIANS);

            latlon[indexLat][point] = xlat;
            latlon[indexLon][point] = (iwest == 1) ? -xlon  : xlon;

        } // end point for loop

        return latlon;

    }

    /**
     * toLinEle converts lat/long to satellite line/element
     *
     * @param  latlon	 array of lat/long pairs. Where latlon[indexLat][]
     *                    are latitudes and latlon[indexLon][] are longitudes.
     *
     * @return linele[][] array of line/element pairs.  Where
     *                    is an element.  These are in 'file' coordinates
     *                    (not "image" coordinates);
     */
    public double[][] toLinEle(double[][] latlon) 
    {
        double xlon;
        double xlat;
        double xrlon, xrlat, xclat;

        int number = latlon[0].length;
        double[][] linele = new double[2][number];

        for (int point=0; point < number; point++) 
        {
            xlat = latlon[indexLat][point];
            // transform to McIDAS (west positive longitude) coordinates
            xlon = (iwest == 1) 
                   ? -latlon[indexLon][point]
                   : latlon[indexLon][point];

            xrlon = ihem*(xlon-xqlon);
            xrlon = (xrlon+900.)%360. - 180.;
            xrlon = xrlon*xfac*DEGREES_TO_RADIANS;
            xclat = (90. - ihem*xlat)*DEGREES_TO_RADIANS*.5;
            if (xclat == 0.0)
                xrlat = 0.0;
            else
                xrlat =  xblat*Math.pow(Math.tan(Math.abs(xclat)), xfac);

            linele[indexLine][point] = xrow + ihem*(xrlat*Math.cos(xrlon));
            linele[indexEle][point]  = xcol - ihem*(xrlat*Math.sin(xrlon));

        } // end point loop

        // Return in 'File' coordinates
        return imageCoordToAreaCoord(linele, linele);
    }
}
