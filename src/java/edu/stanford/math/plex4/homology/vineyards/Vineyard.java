/*
 * Vineyard.java
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
 * This class implements the functionality of a persistence vineyard.
 * Each Vine is essentially just a sequence of points (or Intervals in the barcode language use here)
 * in a sequence of persistence diagrams, along with their time values.  The Vineyard is a collection
 * of Vines.  It makes sense to subclass PersistenceInvariantDescriptor here, since that essentially just
 * pairs topological invariants with their generators.
 * I'm trying to define the class in terms of the existing classes.
 * However I don't think this can be an extension of PersistenceIntervalDescriptor.  I would like to make
 * each Vine be an entry in the that class, but it associates a single generator with each entry.  In a vineyard
 * the generator can change.  Which is too bad because that superclass provides a bunch of bookkeeping
 * and accessor functionality
 * 
 * I don't think we can properly implement comparable for this, unless you can decide on a scalar comparison
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
public class Vineyard<T extends Comparable<T>, K extends Comparable<K>, S, G> {
	/*
	 * The main data structure of a vineyard.
	 * It's just a dimension-indexed set of Lists of Vines
	 */
	protected final Map<Integer, List<Vine<T,K,S, G>>> vines = new HashMap<Integer, List<Vine<T,K,S, G>>>();

	/*
	 * Just trying to be like AnnotatedBarcodeCollection here...
	 */
	protected boolean useLeftClosedDefault = true;
	protected boolean useRightClosedDefault = false;

	/**
	 * Constructor
	 */
	public Vineyard() {}

	/** 
	 * method to add a vine to the vineyard
	 * 
	 * @param vine - the Vine to be added to the vineyard
	 */
	public void addVine(Vine<T,K,S,G> vine) {
		int dim = vine.dimension();
		if (!this.vines.containsKey(dim)) {
			this.vines.put(dim, new ArrayList<Vine<T,K,S,G>>());
		}
		vines.get(dim).add(vine);
	}

}


