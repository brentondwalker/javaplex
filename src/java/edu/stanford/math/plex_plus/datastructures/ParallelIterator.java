package edu.stanford.math.plex_plus.datastructures;

import java.util.Collection;
import java.util.Iterator;

import edu.stanford.math.plex_plus.datastructures.pairs.GenericPair;
import edu.stanford.math.plex_plus.utility.ExceptionUtility;

public class ParallelIterator<T, U> implements Iterator<GenericPair<T, U>> {
	private final Iterator<T> TIterator;
	private final Iterator<U> UIterator;
	
	public ParallelIterator(Collection<T> TCollection, Collection<U> UCollection) {
		ExceptionUtility.verifyNonNull(TCollection);
		ExceptionUtility.verifyNonNull(UCollection);
		ExceptionUtility.verifyEqual(TCollection.size(), UCollection.size());
		
		this.TIterator = TCollection.iterator();
		this.UIterator = UCollection.iterator();
		
	}
	
	public boolean hasNext() {
		return this.TIterator.hasNext();
	}

	public GenericPair<T, U> next() {
		T t = this.TIterator.next();
		U u = this.UIterator.next();
		return new GenericPair<T, U>(t, u);
	}

	public void remove() {
		this.TIterator.remove();
		this.UIterator.remove();
	}

}