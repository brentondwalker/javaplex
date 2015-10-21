/*
 * Vine.java
 *
 * Copyright (C) 2015 Brenton Walker
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package edu.stanford.math.plex4.homology.vineyards;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Vector;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Map.Entry;

import edu.stanford.math.plex4.utility.ComparisonUtility;
import edu.stanford.math.primitivelib.autogen.pair.ObjectObjectPair;
import edu.stanford.math.primitivelib.utility.Infinity;
import edu.stanford.math.plex4.homology.barcodes.*;

/**
 * This class implements the functionality of a single vine in a persistence vineyard.
 * A vine is essentially just a sequence of points (or Intervals in the barcode language use here)
 * in a sequence of persistence diagrams, along with their time values.
 * I'm trying to define the class in terms of the existing classes.
 * 
 * I don't think we can implement comparable for this, unless you can decide on a scalar comparison
 * function for vines.
 * * the comparison fn for intervals looks at infiniteness of endpoints first, then falls
 *   back to comparing the start points
 * * use Interval comparison for single-time values of the vines
 * * could compare the max values of the vines.
 * * could compare the mean values of the vines, or median
 * * could compare the start values of the vines (bad comparison in general)
 *  
 * @author Brenton Walker
 *
 * @param <T> the underlying type for the index set of the filtration- e.g. most likely Integer, Double, or Float
 * @param <K> the "time" parameter for the vineyard - also most likely Integer, Double, or Float
 * @param <S> the type for the basis. Generally Simplex or Cell
 * @param <G> the type for the generator. Generally SparseFormalSum of S or somthing similar
 */
//public class Interval<T extends Comparable<T>> implements Comparable<Interval<T>>, Serializable {

public class Vine<T extends Comparable<T>, K extends Comparable<K>, S, G> {

	/*
	 * The vine is a sequence of points in a sequence of persistence diagrams, along
	 * with the corresponding time values.
	 */
	protected Vector<ObjectObjectPair<K,AnnotatedInterval<T,G>>> vine = new Vector<ObjectObjectPair<K,AnnotatedInterval<T,G>>>();

	protected int dimension;
	
	/*
	 * Some accessor functions
	 */
	public int length() {
		return vine.size();
	}
		
	public K start() {
		return vine.firstElement().getFirst();
	}

	public K end() {
		return vine.lastElement().getFirst();
	}
	
	public int dimension() {
		return this.dimension;
	}

	/**
	 * Constructor
	 * 
	 * @param dim - the homologocal dimension of the Vine
	 */
	public Vine(int dim) {
		this.dimension = dim;
	}
	
	/*
	 * Append a time val to the vine
	 */
	public void append(K time, AnnotatedInterval<T,G> point) {
		if (vine.size()>0) {
			if (vine.lastElement().getFirst().compareTo(time) <= 0) {
				System.out.println("ERROR: Vine.append called on out-of-order time value!");
				return;
			}
		}
		vine.add(new ObjectObjectPair<K,AnnotatedInterval<T,G>>(time, point));
	}
	
	/*
	 * Insert a time val to the vine in the appropriate place
	 * Why would you ever want to do this??  Not sure.  Just adding the method.
	 */
	public void insert(K time, AnnotatedInterval<T,G> point) {
		if (vine.size()==0) {
			vine.add(new ObjectObjectPair<K,AnnotatedInterval<T,G>>(time, point));
		} else {
			for (int i=0; i<vine.size(); i++) {
				if (vine.get(i).getFirst().compareTo(time) < 0) {
					vine.insertElementAt(new ObjectObjectPair<K,AnnotatedInterval<T,G>>(time, point), i+1);
					return;
				}
			}
		}
	}
	
	/*
	 * This returns the interval/point in the vine at the specified time
	 * If the time is before start or after end it returns null.
	 * If the time is between values it returns the last defined interval /less or equal to/ the requested time.
	 */
	public AnnotatedInterval<T,G> slice(K time) {
		if (vine.size()==0) return null;
		if ((start().compareTo(time) > 0) || (end().compareTo(time) < 0)) return null;

		// should do a binary search...
		// for now, whatever...
		Iterator<ObjectObjectPair<K,AnnotatedInterval<T,G>>> vine_iterator = vine.iterator();
		
		while (vine_iterator.hasNext()) {
			ObjectObjectPair<K,AnnotatedInterval<T,G>> s = vine_iterator.next();
			K iter_time = s.getFirst();
			if (time.compareTo(iter_time) <= 0) {
				return s.getSecond();
			}
		}
		
		// something's gone horribly wrong!
		System.out.println("ERROR: Vine.slice() failed!");
		return null;
	}

	/*
	 * This returns the interval/point in the vine at the specified time
	 * Unlike slice() this will return a new Interval interpolated between the neighboring ones.
	 * If the time is before start or after end it returns null.
	 * If the time is between values it returns the defined interval /after/ the requested time.
	 */
	public AnnotatedInterval<T,G> interpolatedSlice(K time) {
		System.out.println("WARNING: Vine.interpolatedSlice() is not implemented yet!");
		return null;
	}


	
}
