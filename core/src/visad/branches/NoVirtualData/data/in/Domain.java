/* Generated by Together */

package visad.data.in;

import visad.*;

/**
 * Instances are mutable.
 */
public class Domain
{
    /**
     * @supplierCardinality 1
     * @directed 
     */
    private VirtualSet		set;

    public Domain(VirtualSet set)
    {
	this.set = set;
    }

    public MathType getMathType()
    {
	return set.getMathType();
    }

    public SetType getSetType()
    {
	return set.getSetType();
    }

    public RealTupleType getRealTupleType()
    {
	return set.getRealTupleType();
    }

    /**
     * Returns the numer of points in the domain.
     */
    public int getLength()
	throws VisADException
    {
	return set.getLength();
    }

    public boolean isFactorable()
    {
	return set.isFactorable();
    }

    /**
     * {@link #isFactorable()} must be <code>true</code>.
     */
    public FactoredDomain factor()
	throws VisADException
    {
	return new FactoredDomain(this);
    }

    public VirtualSet getVirtualSet()
    {
	return set;
    }

    public SampledSet getSet()
    {
	return set.getSet();
    }

    public boolean equals(Object obj)
    {
	boolean	equals;
	if (!(obj instanceof Domain))
	{
	    equals = false;
	}
	else
	{
	    Domain	that = (Domain)obj;
	    equals = this == that || set.equals(that.set);
	}
	return equals;
    }

    public int hashCode()
    {
	return set.hashCode();
    }

    public SampledSet getData()
	throws VisADException
    {
	return set.getSet();
    }

    public static class FactoredDomain
    {
	Domain	innerDomain;
	Domain	outerDomain;

	protected FactoredDomain(Domain domain)
	    throws VisADException
	{
	    VirtualSet.FactoredSet	factoredSet = domain.set.factor();
	    innerDomain = new Domain(factoredSet.getInnerSet());
	    outerDomain = new Domain(factoredSet.getOuterSet());
	}

	public Domain getOuterDomain()
	{
	    return outerDomain;
	}

	public Domain getInnerDomain()
	{
	    return innerDomain;
	}
    }
}