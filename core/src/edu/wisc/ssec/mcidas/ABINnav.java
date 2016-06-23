//
// ABINnav.java
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
 * ABINnav is used to provide {@literal "navigation"} for ABIN image data.
 *
 * This code is essentially a direct port of {@code nvxabin.dlm} from McIDAS-X.
 * Note: the variable naming convention has been retained from the original
 * McIDAS-X source code.
 */
public class ABINnav extends AREAnav {

    /** Radius of satellite orbit (kilometers). */
    private static final double dh = 42164.16000000000349245965480804443359375;

    /** Convert degrees to radians. */
    private static final double deg_to_rad = PI / 180.0;

    /** Radius of the Earth at the equator (kilometers). */
    private static final double r_eq =
        6378.136999999999716237653046846389770508;

    /** Radius of the Earth at the equator squared (kilometers). */
    // note: dreq2 and r_eq2 have been kept around to conform to McIDAS-X code
    private static final double dreq2 = r_eq * r_eq;

    /** Radius of the Earth at the poles (kilometers). */
    private static final double drpo =
        6356.752300000000104773789644241333007812;

    /** Radius of the Earth at the poles squared (kilometers). */
    private static final double drpo2 = drpo * drpo;

    /** Radius of the Earth at the equator squared (kilometers). */
    // note: dreq2 and r_eq2 have been kept around to conform to McIDAS-X code
    private static final double r_eq2 = dreq2;

    /**
     * Pre-calculated constant (value is dh^2 - r_eq^2).
     *
     * @see #dh
     * @see #r_eq
     */
    private static final double d = dh * dh - r_eq2;

    private static final double FP =
        1.006802999999999892466462370066437870264;

    private boolean isEastPositive = true;

    private int itype = 1;

    /** Line offset. */
    private double loff;

    /** Column offset. */
    private double coff;

    /** Line factor. */
    private double lfac;

    /** Column factor. */
    private double cfac;

    /** Subpoint. */
    private double plon;

    /** Base resolution. */
    private double bres;

    /**
     * Initializes the ABIN navigation code with the given set of navigation
     * parameters.
     *
     * @param navblock Navigation parameters from image file.
     *
     * @throws IllegalArgumentException if {@code navblock} is not ABIN.
     */
    public ABINnav(int[] navblock) {
        if (navblock[0] != ABIN) {
            throw new IllegalArgumentException("Invalid navigation type: " +
                                               navblock[0]);
        }
        loff = navblock[1] / 100000000.0;
        coff = navblock[2] / 100000000.0;
        lfac = navblock[3] / 100000000.0;
        cfac = navblock[4] / 100000000.0;
        plon = navblock[5] / 10.0;
        bres = navblock[6];
    }

    /**
     * Convert satellite lines/elements to latitude/longitude coordinates.
     *
     * @param linele Array of line/element pairs.
     *               Where {@code linele[indexLine]} are {@literal "lines"}
     *               and {@code linele[indexEle]} are {@literal "elements"}.
     *               These coordinates must be {@literal "file"} rather than
     *               {@literal "image"} coordinates.
     *
     * @return Array of latitude/longitude pairs. {@code latlon[indexLat]} are
     *         latitudes and {@code latlon[indexLon]} are longitudes.
     */
    public double[][] toLatLon(double[][] linele) {
        final double sub_lon_radians = plon * (PI / 180.0);

        double xlin;
        double xele;
        double lamda_goes;
        double theta_goes;
        double lamda_geos;
        double theta_geos;
        double cosx;
        double cosy;
        double sinx;
        double siny;
        double c1;
        double c2;
        double sd;
        double sdd;
        double sn;
        double s1;
        double s2;
        double s3;
        double sxy;

        int length = linele[indexLine].length;
        double[][] latLons = new double[2][length];
        double[][] imageLineElems = areaCoordToImageCoord(linele);

        for (int point = 0; point < length; point++) {
            double rlin = imageLineElems[indexLine][point];
            double rele = imageLineElems[indexEle][point];

            // start img_to_ll
            xlin = 0.0;
            xele = 0.0;
            lamda_goes = 0.0;
            theta_goes = 0.0;
            lamda_geos = 0.0;
            theta_geos = 0.0;
            cosx = 0.0;
            cosy = 0.0;
            sinx = 0.0;
            siny = 0.0;

            c1 = 0.0;
            c2 = 0.0;
            sd = 0.0;
            sdd = 0.0;
            sn = 0.0;
            s1 = 0.0;
            s2 = 0.0;
            s3 = 0.0;
            sxy = 0.0;

            double xlat;
            double xlon;

            // adjust using Base RESolution
            xlin = (rlin + bres - 1.0) / bres;
            xele = (rele + bres - 1.0) / bres;

            // Intermediate coordinates (coordinates will be radians)
            theta_goes = xlin * lfac + loff;
            lamda_goes = xele * cfac + coff;

            // convert GOES to GEOS
            theta_geos = asin(sin(theta_goes) * cos(lamda_goes));
            lamda_geos = atan(tan(lamda_goes) / cos(theta_goes));

            // SIN and COS for computations below
            cosx = cos(lamda_geos);
            cosy = cos(theta_geos);
            sinx = sin(lamda_geos);
            siny = sin(theta_geos);

            c1 = dh * cosx * cosy * dh * cosx * cosy;
            c2 = (cosy * cosy + FP * siny * siny) * d;

            sdd = c1 - c2;
            if ((sdd < 0.0))  {
                xlat = Double.NaN;
                xlon = Double.NaN;
            } else {
                sd = sqrt(sdd);

                sn = (dh * cosx * cosy - sd) / (cosy * cosy + FP * siny * siny);

                s1 = dh - sn * cosx * cosy;
                s2 = sn * sinx * cosy;
                s3 = -(sn * siny);

                sxy = sqrt(s1 * s1 + (s2 * s2));
                xlon = atan(s2 / s1) + sub_lon_radians;

                xlat = atan(-(FP * s3 / sxy));

                // convert radians to degrees
                xlon = xlon * (180.0 / PI);
                xlat = xlat * (180.0 / PI);

                // Longitudes in [-180,180]
                if ((xlon > 180)) {
                    xlon = xlon - 360.0;
                }
                if ((xlon < -180)) {
                    xlon = xlon + 360.0;
                }
            }
            // end img_to_ll

            latLons[indexLat][point] = xlat;
            latLons[indexLon][point] = xlon;
        }
        return latLons;
    }

    /**
     * Convert latitudes/longitudes to satellite lines/elements.
     *
     * @param latlon Array of latitude/longitude pairs.
     *               Where {@code latlon[indexLat]} are latitudes and
     *               {@code latlon[indexLon]} are longitudes.
     *
     * @return Array of line/element pairs. {@code linele[indexLine]} are lines
     *         and {@code linele[indexEle]} are elements. These coordinates are
     *         {@literal "file"} rather than {@literal "image"} coordinates.
     */
    public double[][] toLinEle(double[][] latlon) {
        final double d_geographic_ssl = plon * deg_to_rad;

        double rlin;
        double rele;
        double d_geographic_lat;
        double d_geocentric_lat;
        double d_geographic_lon;
        double r_earth;
        double r_1;
        double r_2;
        double r_3;
        double lamda;
        double theta;

        int length = latlon[indexLat].length;
        double[][] lineEles = new double[2][length];

        for (int point = 0; point < length; point++) {
            double rlat = latlon[indexLat][point];
            double rlon = latlon[indexLon][point];
            if (!isEastPositive) {
                rlon = -rlon;
            }

            // start ll_to_img
            rlin = 0.0;
            rele = 0.0;
            d_geographic_lat = 0.0;
            d_geocentric_lat = 0.0;
            d_geographic_lon = 0.0;
            r_earth = 0.0;
            r_1 = 0.0;
            r_2 = 0.0;
            r_3 = 0.0;
            lamda = 0.0;
            theta = 0.0;

            double xlin;
            double xele;

            // Earth (Geographic) Coordinates are converted to Radians
            d_geographic_lat = rlat * deg_to_rad;
            d_geographic_lon = rlon * deg_to_rad;

            d_geocentric_lat = atan(drpo2 / dreq2 * tan(d_geographic_lat));

            r_earth = drpo / sqrt(1.0 - (dreq2 - drpo2) / dreq2 * cos(d_geocentric_lat) * cos(d_geocentric_lat));

            r_1 = dh - r_earth * cos(d_geocentric_lat) * cos(d_geographic_lon - d_geographic_ssl);

            r_2 = -(r_earth * cos(d_geocentric_lat) * sin(d_geographic_lon - d_geographic_ssl));

            r_3 = r_earth * sin(d_geocentric_lat);

            if ((r_1 > dh))  {
                xlin = Double.NaN;
                xele = Double.NaN;
            } else {
                lamda = asin(-(r_2 / sqrt(r_1 * r_1 + r_2 * r_2 + r_3 * r_3)));
                theta = atan(r_3 / r_1);

                // image line and element
                rlin = (theta - loff) / lfac;
                rele = (lamda - coff) / cfac;

                // Adjust using Base RESolution
                xlin = rlin * bres - (bres - 1.0);
                xele = rele * bres - (bres - 1.0);
            }
            // end of ll_to_img

            lineEles[indexLine][point] = xlin;
            lineEles[indexEle][point] = xele;
        }
        return imageCoordToAreaCoord(lineEles, lineEles);
    }

    /**
     * {@inheritDoc}
     */
    @Override public boolean canCalculateAngles() {
        return true;
    }

    public static int iday = 0;
    public static float r = 6371.221f;

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
        final double xplat = 0.0;
        final double xplon = plon / 10.0;
        final double snlt = sin(xplat * rdpdg);
        final double cslt = cos(xplat * rdpdg);
        final double csln = cos(xplon * rdpdg);
        final double snln = sin(xplon * rdpdg);

//        xs[0] = 42164.36499999999796273186802864074707031 * cslt * csln / 6378.136999999999716237653046846389770508;
//        xs[1] = 42164.36499999999796273186802864074707031 * cslt * snln / 6378.136999999999716237653046846389770508;
//        xs[2] = 42164.36499999999796273186802864074707031 * snlt / 6378.136999999999716237653046846389770508;
        final double[] xs = {
            42164.36499999999796273186802864074707031 * cslt * csln / 6378.136999999999716237653046846389770508,
            42164.36499999999796273186802864074707031 * cslt * snln / 6378.136999999999716237653046846389770508,
            42164.36499999999796273186802864074707031 * snlt / 6378.136999999999716237653046846389770508,
        };

        double xsat = xs[0] * 6378.136999999999716237653046846389770508;

//        double ylat = 0.0;

        int inorb = 0;

//
//            if ((iday == jday)) {
//                Dummy.go_to("nvxabin/Abinang",1);
//            }
        iday = jday;
        inorb = 0;
//        label1:
//            Dummy.label("nvxabin/Abinang",1);
        final double pictim = McIDASUtil.mcPackedIntegerToDouble(jtime);

        // determine satellite position
        xsat = xs[0] * 6378.136999999999716237653046846389770508;
        final double ysat = xs[1] * 6378.136999999999716237653046846389770508;
        final double zsat = xs[2] * 6378.136999999999716237653046846389770508;

        final double height = sqrt(pow(xsat, 2) + pow(ysat, 2) + pow(zsat, 2));
        final double ylat = AREAnav.geolat(rdpdg * xlat, 1);
        final double ylon = rdpdg * xlon;
        final double slat = sin(ylat);
        final double clat = cos(ylat);
        final double slon = sin(ylon);
        final double clon = cos(ylon);
        final double xsam = r * clat * clon;
        final double ysam = r * clat * slon;
        final double zsam = r * slat;

        // determine zenith angle of sun
        final double snlg = -(pictim * PI / 12.0) - rdpdg * gha;
        final double sndc = rdpdg * dec;
        final double cosdec = cos(sndc);
        final double us = cos(snlg) * cosdec;
        final double vs = sin(snlg) * cosdec;
        final double ws = sin(sndc);
        final double sunang = acos((us * xsam + vs * ysam + ws * zsam) / r) / rdpdg;

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