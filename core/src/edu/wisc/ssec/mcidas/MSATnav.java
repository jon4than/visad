//
// MSATnav.java
//

/*
This source file is part of the edu.wisc.ssec.mcidas package and is
Copyright (C) 1998 - 2017 by Tom Whittaker, Tommy Jasmin, Tom Rink,
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

import static java.lang.Math.PI;
import static java.lang.Math.abs;
import static java.lang.Math.acos;
import static java.lang.Math.asin;
import static java.lang.Math.atan;
import static java.lang.Math.atan2;
import static java.lang.Math.cos;
import static java.lang.Math.pow;
import static java.lang.Math.sin;
import static java.lang.Math.sqrt;
import static java.lang.Math.tan;

/**
 * Navigation class for Meteosat (MSAT) type nav. This code was modified
 * from the original FORTRAN code (nvxmsat.dlm) on the McIDAS system. It
 * only supports latitude/longitude to line/element transformations (LL) 
 * and vice/versa. Transform to 'XYZ' not implemented.
 * @see <A HREF="http://www.ssec.wisc.edu/mcidas/doc/prog_man.html">
 *      McIDAS Programmer's Manual</A>
 *
 * @author  Don Murray
 */
public final class MSATnav extends AREAnav {

    private boolean isEastPositive = true;

    final double NOMORB = 42164.;   // nominal radial distance of satellite (km)
    final double EARTH_RADIUS = 6378.155; // earth equatorial radius (km)

    int itype;
    double h;
    double a;
    double rp;
    double lpsi2;
    double deltax;
    double deltay;
    double rflon;                  // reference longitude;
    double sublon;
    int[] ioff = new int[3];

    /**
     * Set up for the real math work.  Must pass in the int array
     * of the MSAT nav 'codicil'.
     *
     * @param iparms the nav block from the image file
     * @throws IllegalArgumentException if the nav block is not a MSAT type.
     */
    public MSATnav(int[] iparms)
            throws IllegalArgumentException {

/* No longer needed.  Kept for consistency with nvxmsat.dlm
        if (ifunc != 1) 
        {
            if (iparms[0] == XY ) itype = 1;
            if (iparms[0] == LL ) itype = 2;
            return;
        }
*/

        if (iparms[0] != MSAT)
            throw new IllegalArgumentException("Invalid navigation type" +
                    iparms[0]);
        itype = 2;

        System.arraycopy(iparms, 3, ioff, 0, 3);
        h = (double) NOMORB - EARTH_RADIUS;
        a = 1. / 297.;
        rp = EARTH_RADIUS / (1. + a);
        lpsi2 = 1;
        deltax = 18. / 2500.;
        deltay = 18. / 2500.;
        rflon = 0.0;
        sublon = McIDASUtil.integerLatLonToDouble(iparms[6]);
    }

    /**
     * converts from satellite coordinates to latitude/longitude
     *
     * @param linele array of line/element pairs.  Where
     *               linele[indexLine][] is a 'line' and
     *               linele[indexEle][] is an element. These are in
     *               'file' coordinates (not "image" coordinates.)
     * @return latlon[][]  array of lat/long pairs. Output array is
     * latlon[indexLat][] of latitudes and
     * latlon[indexLon][] of longitudes.
     */
    public double[][] toLatLon(double[][] linele) {

        double xele, xlin;
        double xele2, xlin2;
        double xfi, xla;
        double x, y;
        double xr, yr;
        double tanx, tany;
        double val1, val2;
        double yk;
        double vmu;
        double cosrf, sinrf;
        double teta;
        double xt, yt, zt;
        double rs;

        int number = linele[0].length;
        double[][] latlon = new double[2][number];

        // Convert array to Image coordinates for computations
        double[][] imglinele = areaCoordToImageCoord(linele);

        for (int point = 0; point < number; point++) {
            xlin = imglinele[indexLine][point];
            xele = imglinele[indexEle][point];

            xele2 = xele / 2.;
            xlin2 = xlin / 2.;
            x = 1250.5 - xele2;
            y = ioff[2] - (xlin2 + ioff[1] - ioff[0]);
            xr = x;
            yr = y;
            x = xr * lpsi2 * deltax * DEGREES_TO_RADIANS;
            y = yr * lpsi2 * deltay * DEGREES_TO_RADIANS;
            rs = EARTH_RADIUS + h;
            tanx = tan(x);
            tany = tan(y);
            val1 = 1. + tanx * tanx;
            val2 = 1. + (tany * tany) * ((1. + a) * (1. + a));
            yk = rs / EARTH_RADIUS;
            if ((val1 * val2) > ((yk * yk) / (yk * yk - 1))) {
                latlon[indexLat][point] = Double.NaN;
                latlon[indexLon][point] = Double.NaN;
            } else {
                vmu = (rs - (EARTH_RADIUS * (sqrt((yk * yk) -
                        (yk * yk - 1) * val1 * val2)))) / (val1 * val2);
                cosrf = cos(rflon * DEGREES_TO_RADIANS);
                sinrf = sin(rflon * DEGREES_TO_RADIANS);
                xt = (rs * cosrf) + (vmu * (tanx * sinrf - cosrf));
                yt = (rs * sinrf) - (vmu * (tanx * cosrf + sinrf));
                zt = vmu * tany / cos(x);
                teta = asin(zt / rp);
                xfi = (atan(((tan(teta)) * EARTH_RADIUS) / rp)) *
                        RADIANS_TO_DEGREES;
                xla = -atan(yt / xt) * RADIANS_TO_DEGREES;

                // change longitude for correct subpoint
                xla = xla + sublon;

                //  put longitude into East Positive (form)
                if (isEastPositive) xla = -xla;

                latlon[indexLat][point] = xfi;
                latlon[indexLon][point] = xla;
            }  // end lat/lon point calculation 
        } // end point for loop

        return latlon;

    }

    /**
     * toLinEle converts lat/long to satellite line/element
     *
     * @param latlon array of lat/long pairs. Where latlon[indexLat][]
     *               are latitudes and latlon[indexLon][] are longitudes.
     * @return linele[][] array of line/element pairs.  Where
     * linele[indexLine][] is a line and linele[indexEle][]
     * is an element.  These are in 'file' coordinates
     * (not "image" coordinates);
     */
    public double[][] toLinEle(double[][] latlon) {
        double y;
        double x1, y1;
        double xfi, xla;
        double rom;
        double r1, r2;
        double coslo, sinlo;
        double teta;
        double xt, yt, zt;
        double px, py;
        double rs;
        double reph, rpph;
        double xr, yr;

        int number = latlon[0].length;
        double[][] linele = new double[2][number];

        for (int point = 0; point < number; point++) {

            x1 = latlon[indexLat][point];

            // expects positive East Longitude.
            y1 = isEastPositive
                    ? latlon[indexLon][point]
                    : -latlon[indexLon][point];

            // if in cartesian coords, transform to lat/lon
            if (itype == 1) {
                y = latlon[indexLon][point];
                // NXYZLL(x,y,z,zlat,zlon);
                y1 = -y1;
            }

            // correct for sublon
            y1 = y1 + sublon;
            xfi = x1 * DEGREES_TO_RADIANS;
            xla = y1 * DEGREES_TO_RADIANS;
            rom =
                    (EARTH_RADIUS * rp) /
                            sqrt(
                                    rp * rp * cos(xfi) * cos(xfi) +
                                            EARTH_RADIUS * EARTH_RADIUS * sin(xfi) * sin(xfi));
            y = sqrt(h * h + rom * rom - 2 * h * rom * cos(xfi) * cos(xla));
            r1 = y * y + rom * rom;
            r2 = h * h;
            if (r1 > r2)  // invalid point
            {
                linele[indexLine][point] = Double.NaN;
                linele[indexEle][point] = Double.NaN;
            } else          // calculate line an element
            {
                rs = EARTH_RADIUS + h;
                reph = EARTH_RADIUS;
                rpph = rp;
                coslo = cos(rflon * DEGREES_TO_RADIANS);
                sinlo = sin(rflon * DEGREES_TO_RADIANS);
                teta = atan((rpph / reph) * tan(xfi));
                xt = reph * cos(teta) * cos(xla);
                yt = reph * cos(teta) * sin(xla);
                zt = rpph * sin(teta);

                px = atan((coslo * (yt - rs * sinlo) - (xt - rs * coslo) * sinlo) /
                        (sinlo * (yt - rs * sinlo) + (xt - rs * coslo) * coslo));
                py = atan(zt * ((tan(px) * sinlo -
                        coslo) / (xt - rs * coslo)) * cos(px));
                px = px * RADIANS_TO_DEGREES;
                py = py * RADIANS_TO_DEGREES;
                xr = px / (deltax * lpsi2);
                yr = py / (deltay * lpsi2);
                xr = 1250.5 - xr;
                yr = yr + ioff[2] + ioff[1] - ioff[0];
                xr = xr * 2;
                yr = 5000 - yr * 2;
                linele[indexLine][point] = yr;
                linele[indexEle][point] = xr;

            }  // end calculations
        } // end point loop

        // Return in 'File' coordinates
        return imageCoordToAreaCoord(linele, linele);
    }

    /**
     * converts from satellite coordinates to latitude/longitude
     *
     * @param linele array of line/element pairs.  Where
     *               linele[indexLine][] is a 'line' and
     *               linele[indexEle][] is an element. These are in
     *               'file' coordinates (not "image" coordinates.)
     * @return latlon[][]  array of lat/long pairs. Output array is
     * latlon[indexLat][] of latitudes and
     * latlon[indexLon][] of longitudes.
     */
    public float[][] toLatLon(float[][] linele) {

        double xele, xlin;
        double xele2, xlin2;
        double xfi, xla;
        double x, y;
        double xr, yr;
        double tanx, tany;
        double val1, val2;
        double yk;
        double vmu;
        double cosrf, sinrf;
        double teta;
        double xt, yt, zt;
        double rs;

        int number = linele[0].length;
        float[][] latlon = new float[2][number];

        // Convert array to Image coordinates for computations
        float[][] imglinele = areaCoordToImageCoord(linele);

        for (int point = 0; point < number; point++) {
            xlin = imglinele[indexLine][point];
            xele = imglinele[indexEle][point];

            xele2 = xele / 2.;
            xlin2 = xlin / 2.;
            x = 1250.5 - xele2;
            y = ioff[2] - (xlin2 + ioff[1] - ioff[0]);
            xr = x;
            yr = y;
            x = xr * lpsi2 * deltax * DEGREES_TO_RADIANS;
            y = yr * lpsi2 * deltay * DEGREES_TO_RADIANS;
            rs = EARTH_RADIUS + h;
            tanx = tan(x);
            tany = tan(y);
            val1 = 1. + tanx * tanx;
            val2 = 1. + (tany * tany) * ((1. + a) * (1. + a));
            yk = rs / EARTH_RADIUS;
            if ((val1 * val2) > ((yk * yk) / (yk * yk - 1))) {
                latlon[indexLat][point] = Float.NaN;
                latlon[indexLon][point] = Float.NaN;
            } else {
                vmu = (rs - (EARTH_RADIUS * (sqrt((yk * yk) -
                        (yk * yk - 1) * val1 * val2)))) / (val1 * val2);
                cosrf = cos(rflon * DEGREES_TO_RADIANS);
                sinrf = sin(rflon * DEGREES_TO_RADIANS);
                xt = (rs * cosrf) + (vmu * (tanx * sinrf - cosrf));
                yt = (rs * sinrf) - (vmu * (tanx * cosrf + sinrf));
                zt = vmu * tany / cos(x);
                teta = asin(zt / rp);
                xfi = (atan(((tan(teta)) * EARTH_RADIUS) / rp)) *
                        RADIANS_TO_DEGREES;
                xla = -atan(yt / xt) * RADIANS_TO_DEGREES;

                // change longitude for correct subpoint
                xla = xla + sublon;

                //  put longitude into East Positive (form)
                if (isEastPositive) xla = -xla;

                latlon[indexLat][point] = (float) xfi;
                latlon[indexLon][point] = (float) xla;
            }  // end lat/lon point calculation 
        } // end point for loop

        return latlon;

    }

    /**
     * toLinEle converts lat/long to satellite line/element
     *
     * @param latlon array of lat/long pairs. Where latlon[indexLat][]
     *               are latitudes and latlon[indexLon][] are longitudes.
     * @return linele[][] array of line/element pairs.  Where
     * linele[indexLine][] is a line and linele[indexEle][]
     * is an element.  These are in 'file' coordinates
     * (not "image" coordinates);
     */
    public float[][] toLinEle(float[][] latlon) {
        double y;
        double x1, y1;
        double xfi, xla;
        double rom;
        double r1, r2;
        double coslo, sinlo;
        double teta;
        double xt, yt, zt;
        double px, py;
        double rs;
        double reph, rpph;
        double xr, yr;

        int number = latlon[0].length;
        float[][] linele = new float[2][number];

        for (int point = 0; point < number; point++) {

            x1 = latlon[indexLat][point];

            // expects positive East Longitude.
            y1 = isEastPositive
                    ? latlon[indexLon][point]
                    : -latlon[indexLon][point];

            // if in cartesian coords, transform to lat/lon
            if (itype == 1) {
                y = latlon[indexLon][point];
                // NXYZLL(x,y,z,zlat,zlon);
                y1 = -y1;
            }

            // correct for sublon
            y1 = y1 + sublon;
            xfi = x1 * DEGREES_TO_RADIANS;
            xla = y1 * DEGREES_TO_RADIANS;
            rom =
                    (EARTH_RADIUS * rp) /
                            sqrt(
                                    rp * rp * cos(xfi) * cos(xfi) +
                                            EARTH_RADIUS * EARTH_RADIUS * sin(xfi) * sin(xfi));
            y = sqrt(h * h + rom * rom - 2 * h * rom * cos(xfi) * cos(xla));
            r1 = y * y + rom * rom;
            r2 = h * h;
            if (r1 > r2)  // invalid point
            {
                linele[indexLine][point] = Float.NaN;
                linele[indexEle][point] = Float.NaN;
            } else          // calculate line an element
            {
                rs = EARTH_RADIUS + h;
                reph = EARTH_RADIUS;
                rpph = rp;
                coslo = cos(rflon * DEGREES_TO_RADIANS);
                sinlo = sin(rflon * DEGREES_TO_RADIANS);
                teta = atan((rpph / reph) * tan(xfi));
                xt = reph * cos(teta) * cos(xla);
                yt = reph * cos(teta) * sin(xla);
                zt = rpph * sin(teta);

                px = atan((coslo * (yt - rs * sinlo) - (xt - rs * coslo) * sinlo) /
                        (sinlo * (yt - rs * sinlo) + (xt - rs * coslo) * coslo));
                py = atan(zt * ((tan(px) * sinlo -
                        coslo) / (xt - rs * coslo)) * cos(px));
                px = px * RADIANS_TO_DEGREES;
                py = py * RADIANS_TO_DEGREES;
                xr = px / (deltax * lpsi2);
                yr = py / (deltay * lpsi2);
                xr = 1250.5 - xr;
                yr = yr + ioff[2] + ioff[1] - ioff[0];
                xr = xr * 2;
                yr = 5000 - yr * 2;
                linele[indexLine][point] = (float) yr;
                linele[indexEle][point] = (float) xr;

            }  // end calculations
        } // end point loop

        // Return in 'File' coordinates
        return imageCoordToAreaCoord(linele, linele);
    }

    /**
     * {@inheritDoc}
     */
    @Override public boolean canCalculateAngles() {
        return true;
    }

//    public static int iday = 0;
    public static final double r = 6371.221;

    /**
     * {@inheritDoc}
     */
    @Override public double[] angles(final int jday,
                                     final int jtime,
                                     final double xlat,
                                     final double xlon,
                                     final double gha,
                                     final double dec)
    {
        final double rdpdg = PI / 180.0;

        double[] xs = new double[3];

        int inorb = 0;

        // C
        // C  XS(1) is the length along the x-axis, that is a line from the
        // C	center of Earth through (0,0).  The y-axis (XS(2)) is from
        // C	the center of the Earth through (0,-90).  The z coordinate
        // C	is 0 because we're in the equatorial plane
        // C  XS(1) is 42164.365 (satellite height above Earth center) * cos 50.25

        // C          it is < 0 because the satellite is west of 90 E
        // C  XS(2) is 42164.365*sin 50.25...It's >0 because the Sat is over the E
        xs[1] = 42164.36499999999796273186802864074707031 * sin(39.75 * rdpdg) / 6378.136999999999716237653046846389770508;
        xs[0] = -(42164.36499999999796273186802864074707031 * cos(39.75 * rdpdg) / 6378.136999999999716237653046846389770508);
        // wtf!?
        xs[2] = (6378.38816 + 42164.0) / 6378.136999999999716237653046846389770508;
        xs[2] = 42164.0 / 6378.136999999999716237653046846389770508;
        xs[2] = 0.0;

//        if ((iday == jday)) {
//            Dummy.go_to("nvxmtst/Mtstang",1);
//        }
//        iday = jday;
        inorb = 0;
//        label1:
//        Dummy.label("nvxmtst/Mtstang",1);

        final double pictim = McIDASUtil.mcPackedIntegerToDouble(jtime);

        // determine satellite position
        final double xsat = xs[0] * 6378.136999999999716237653046846389770508;
        final double ysat = xs[1] * 6378.136999999999716237653046846389770508;
        final double zsat = xs[2] * 6378.136999999999716237653046846389770508;

        final double height = sqrt(pow(xsat, 2) + pow(ysat, 2) + pow(zsat, 2));
        final double ylat = AREAnav.geolat(rdpdg * xlat, 1);
        final double ylon = rdpdg * xlon;
        final double slat = sin(ylat);
        final double clat = cos(ylat);
        final double slon = sin(ylon);
        final double clon = cos(ylon);
        final double xsam = (r * clat) * clon;
        final double ysam = (r * clat) * slon;
        final double zsam = r * slat;

        // determine zenith angle of sun
        final double snlg = -(pictim * PI / 12.0) - rdpdg * gha;
        final double sndc = rdpdg * dec;
        final double cosdec = cos(sndc);
        final double us = cos(snlg) * cosdec;
        final double vs = sin(snlg) * cosdec;
        final double ws = sin(sndc);
        final double sunang = acos((us*xsam + vs * ysam + ws * zsam) / r) / rdpdg;

        // determine zenith angle of satellite
        final double xvec = xsat - xsam;
        final double yvec = ysat - ysam;
        final double zvec = zsat - zsam;
        final double xfact = sqrt(pow(xvec, 2) + pow(yvec, 2) + pow(zvec, 2));
        final double satang = acos((xvec * xsam + yvec * ysam + zvec * zsam) / (r * xfact)) / rdpdg;

        // determine relative angle
        final double x1 = clat * clon;
        final double y1 = clat * slon;
        final double z1 = slat;
        final double x2 = slon;
        final double y2 = -clon;
        final double x3 = -(slat * clon);
        final double y3 = -(slat * slon);
        final double z3 = clat;
        final double xc1 = us - x1;
        final double yc1 = vs - y1;
        final double zc1 = ws - z1;
        final double xc2 = xsat / height - x1;
        final double yc2 = ysat / height - y1;
        final double zc2 = zsat / height - z1;
        final double xan1 = xc1 * x3 + yc1 * y3 + zc1 * z3;
        final double xan2 = xc2 * x3 + yc2 * y3 + zc2 * z3;
        final double yan1 = xc1 * x2 + yc1 * y2;
        final double yan2 = xc2 * x2 + yc2 * y2;
        final double xan3 = xan1 * xan2 + yan1 * yan2;
        final double yan3 = -(yan1 * xan2) + xan1 * yan2;
        final double relang = abs(atan2(yan3, xan3) / rdpdg);

        return new double[] { satang, sunang, relang };
    }
}
