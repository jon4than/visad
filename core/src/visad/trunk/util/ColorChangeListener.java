/*

@(#) $Id: ColorChangeListener.java,v 1.5 1999-07-14 21:57:37 dglo Exp $

VisAD Utility Library: Widgets for use in building applications with
the VisAD interactive analysis and visualization library
Copyright (C) 1998 Nick Rasmussen
VisAD is Copyright (C) 1996 - 1998 Bill Hibbard, Curtis Rueden, Tom
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

package visad.util;

import java.util.EventListener;

/** 
 * The interface that all objects must implement to recieve color change
 * events from the color widget 
 *
 * @author Nick Rasmussen nick@cae.wisc.edu
 * @version $Revision: 1.5 $, $Date: 1999-07-14 21:57:37 $
 * @since Visad Utility Library, 0.5
 */

public interface ColorChangeListener extends EventListener {

	/** The function called when the color widget has changed */
	void colorChanged(ColorChangeEvent e);

}
