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

$Id: DataVisitor.java,v 1.4 1999-07-14 22:01:39 dglo Exp $
*/

package visad.data;


import java.rmi.RemoteException;
import visad.FlatField;
import visad.Tuple;
import visad.VisADException;


/**
 * Abstract class for visiting a VisAD data object.  The derived,
 * concrete subclasses are data-form dependent.  The default action
 * upon visiting a VisAD data object is to do nothing and tell the caller
 * to continue.
 */
public abstract class
DataVisitor
{
    /**
     * Visit a VisAD Tuple.
     *
     * @param tuple	The VisAD Tuple being visited.
     * @precondition	<code>tuple</code> is non-null.
     * @postcondition	<code>tuple</code> has been visited.
     * @exception BadFormException	The Tuple doesn't fit the data model
     *					used by the visitor.
     * @exception VisADException	Core VisAD problem (probably couldn't
     *					create a VisAD object).
     * @see visad.data.DataNode
     */
    public boolean
    visit(Tuple tuple)
	throws BadFormException, VisADException, RemoteException
    {
	return true;
    }


    /**
     * Visit a VisAD FlatField.
     *
     * @param field	The VisAD FlatField being visited.
     * @precondition	<code>field</code> is non-null.
     * @postcondition	<code>field</code> has been visited.
     * @exception BadFormException	The Tuple doesn't fit the data model
     *					used by the visitor.
     * @exception VisADException	Core VisAD problem (probably couldn't
     *					create a VisAD object).
     * @see visad.data.DataNode
     */
    public boolean
    visit(FlatField field)
	throws BadFormException, VisADException, RemoteException
    {
	return true;
    }
}
