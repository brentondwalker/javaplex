/*
 * BooleanPLHomologyVineyard.java
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

import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.Set;
import java.util.Vector;

import edu.stanford.math.plex4.streams.interfaces.AbstractFilteredStream;
import edu.stanford.math.plex4.streams.utility.FilteredComparator;
import edu.stanford.math.primitivelib.autogen.formal_sum.BooleanPrimitiveFreeModule;
import edu.stanford.math.primitivelib.autogen.formal_sum.BooleanSparseFormalSum;
import edu.stanford.math.primitivelib.autogen.pair.ObjectObjectPair;
import edu.stanford.math.plex4.homology.barcodes.AnnotatedInterval;
import edu.stanford.math.plex4.homology.vineyards.Vine;
import edu.stanford.math.plex4.homology.vineyards.Vineyard;
import gnu.trove.THashMap;
import gnu.trove.THashSet;


public class BooleanPLHomologyVineyard<U> {

	/**
	 * Whether or not to check our incremental updates against a fresh run of the
	 * persistence algorithm.
	 */
	private static final boolean sanity_check = false;
	
	/**
	 * The matrices for the persistence alg stored in sparse form
	 * These cannot be local to the pHcol() method as in the standard homology algs,
	 * because we need to update them as the filtration changes.
	 * 
	 * We are going to use D=R*W decomposition to match the vineyard alg description, instead
	 * of the standard R=D*V decomp of the PHcol alg.
	 * 
	 * Here R is a sparse representation of the columns, as usual, and W is a sparse
	 * representation of the /rows/.  Kind-of the opposite of the usual
	 */
	THashMap<U, BooleanSparseFormalSum<U>> R = new THashMap<U, BooleanSparseFormalSum<U>>();
	THashMap<U, BooleanSparseFormalSum<U>> W = new THashMap<U, BooleanSparseFormalSum<U>>();

	/**
	 * This maps a simplex to the set of columns containing the key as its low value.
	 */
	THashMap<U, THashSet<U>> lowMap = new THashMap<U, THashSet<U>>();
	
	/**
	 * This is an index of the vines by their birth simplices.
	 * Every interval has a birth simplex.  Not every interval has a death simplex.
	 * This will need to be updated whenever the pairing changes.
	 */
	THashMap<U, Vine<Double,Double,U,BooleanSparseFormalSum<U>>> vineIndex = new THashMap<U, Vine<Double,Double,U,BooleanSparseFormalSum<U>>>();
	
	/**
	 * It's possible that a single update with a new homotopy frame will cause
	 * several transpositions.  Since finalized streams (and their comparators)
	 * don't allow updates, we need to keep track of which simplices have swapped order
	 * relative to the "current" stream, so we can reverse the result of comparisons.
	 * 
	 * This is a /little/ bit hackish.  The ideal solution would be to create mutable
	 * streams that allow for swaps, and that update their comparator functions.  But
	 * I'm not sure even then that the comparators would accomodate arbitrary simplex
	 * swaps, since if two simplices have the same filtration value their order is
	 * lexicographic.
	 * 
	 * The key in this hash is intended to be "x:y" where x and y are the string
	 * representations of the two simplices in question.  It's a cheap way
	 * to do a 2D hash, since Java doesn't have a default Pair class.
	 */
	THashMap<String,Boolean> orderReversals = new THashMap<String,Boolean>();
	
	/**
	 * Keep track of the most recent frame of the homotopy.
	 * When the user passes in a new stream will need to compare to this one.
	 */
	AbstractFilteredStream<U> current_stream = null;
	
	/**
	 * This objects performs the chain computations.
	 */
	protected final BooleanPrimitiveFreeModule<U> chainModule;
	
	/**
	 * This comparator defines the ordering on the basis elements.
	 */
	protected final Comparator<U> basisComparator;
	
	/**
	 * This comparator provides the dictionary ordering on filtration value - basis element
	 * pairs.
	 */
	protected Comparator<U> filteredComparator = null;

	/**
	 * This stores the minimum dimension for which to compute (co)homology.
	 */
	protected int minDimension = 0;
	
	/**
	 * This stores the maximum dimension for which to compute (co)homology.
	 */
	protected int maxDimension = 2;	

	/**
	 * Keep track of whether or not the first frame of the homotopy has been passed in yet
	 */
	protected boolean initialized = false;
	
	/**
	 * A simple class for representing transpositions
	 * Really a transposition is just a pair of "U" type things bein swapped, but
	 * java has no Pair class, and arrays are not covariant, so we can't allocate U[2],
	 * so making this class is as sensible as any other solution.
	 */
	protected class Transposition<UV> implements Comparable<Transposition<UV>> {
		public UV a, b;
		public double t_index, t;
		
		public Transposition(double t, double t_index, UV a, UV b) {
			this.t = t;
			this.t_index = t_index;
			this.a = a;
			this.b = b;
		}
		
		public int compareTo(Transposition<UV> otr) throws ClassCastException {
			if (this.t_index < otr.t_index) return -1;
			if (this.t_index > otr.t_index) return 1;
			return 0;
		}
	}
	
	/**
	 * Constructor
	 */
	public BooleanPLHomologyVineyard(Comparator<U> basisComparator, int minDimension, int maxDimension) {
		this.chainModule = new BooleanPrimitiveFreeModule<U>();
		this.basisComparator = basisComparator;
		this.minDimension = minDimension;
		this.maxDimension = maxDimension;
	}
	
	/**
	 * A debugging tool. Print out the current pairing info.
	 */
	public void printPairingInfo() {
		System.out.println("\n\t"+lowMap+"\n");
	}
	
	/**
	 * This is the main entry point for people using this class.  Allows the user to repeatedly send in
	 * filtered streams.  On the first call this just does the normal decomposition.  On later calls 
	 * it figures out the transpositions needed and updates the decomposition and pairing.
	 * 
	 * @param stream
	 */
	public void addHomotopyFrame(AbstractFilteredStream<U> stream, double t) {
		System.out.println("addHomotopyFrame( "+stream+", "+t+")");
		// if this is the first frame, compute the initial D=R*U decomposition and return
		if (!this.initialized) {
			System.out.println("initializing...");
			this.filteredComparator = new FilteredComparator<U>(stream, this.basisComparator);
			this.current_stream = stream;
			DRUpHcol(stream);
			printPairingInfo();
			this.initialized = true;
			
			// After computing the initial pairing, create the first frame of the vineyard
			Vineyard<Double, Double, U, BooleanSparseFormalSum<U>> vineyard = new Vineyard<Double, Double, U, BooleanSparseFormalSum<U>>();
			for (U k : lowMap.keySet()) {
				Vine<Double, Double, U, BooleanSparseFormalSum<U>> v = new Vine<Double, Double, U, BooleanSparseFormalSum<U>>(stream.getDimension(k));
				vineyard.addVine(v);
				vineIndex.put(k, v);
				// now add the first homotopy frame info to the vine...
				AnnotatedInterval<Double,BooleanSparseFormalSum<U>> point = null;
				if (lowMap.get(k) == null) {
					// this is an unpaired simplex
					// how to get the generator?  Need to look up an example....

					point = AnnotatedInterval.makeRightInfiniteClosedInterval(stream.getFiltrationValue(k), null);
					// new AnnotatedInterval<Double,BooleanSparseFormalSum<U>>(stream.getFiltrationValue(k), 0.0, true, false, false, true, null);	
				} else {
					// it's a paired simplex
					//point = AnnotatedInterval.makeFiniteClosedInterval(stream.getFiltrationValue(k), stream.getFiltrationValue(lowMap.get(k)), null);
				}
				
				v.append(t, point);
			}

		}
		
		// keep track of what simplices have been swapped
		orderReversals.clear();
		
		// otherwise figure out the list of transpositions and perform updates on the decomp
		Vector<Transposition<U>> transpositions = generateTranspositionList(this.current_stream, stream);
		for (Transposition<U> tr : transpositions) {
			// each transposition is a pair of adjacent (they better be)
			// simplices that are being swapped
			updatePairingSwap(tr.a, tr.b);
			printPairingInfo();
		}
		this.filteredComparator = new FilteredComparator<U>(stream, this.basisComparator);
		this.current_stream = stream;		

		if (BooleanPLHomologyVineyard.sanity_check) {
			BooleanPLHomologyVineyard<U> vinetest = new BooleanPLHomologyVineyard<U>(this.basisComparator, 0, maxDimension);
	    	vinetest.addHomotopyFrame(stream, 0.0);
	    	if (vinetest.pairingsEqual(this.lowMap)) {
	    		System.out.println("*** SANITY CHECK SUCCESS! ***");
	    	} else {
	    		System.out.println("*** SANITY CHECK FAIL! ***");
	    	}
		}
	}


	/**
	 * A method to test if a pairing map is equal to our own.
	 * Used for sanity checking
	 * 
	 * @param omap
	 * @return
	 */
	boolean pairingsEqual(THashMap<U, THashSet<U>> omap) {
		
		if (! omap.keySet().equals(lowMap.keySet())) {
			return false;
		}
		for (U k : omap.keySet()) {
			if (! omap.get(k).equals(lowMap.get(k))) {
				return false;
			}
		}
		
		return true;
	}
	
	
	/**
	 * Generate the list of transpositions necessary to turn one stream into another.
	 * This is done in a naiive way, that isn't efficient for RAM or speed, but who cares.
	 * Dimitriy has an alg that uses CGAL's kinetic sort, which is more efficient.
	 * I'd like to do the kinetic sort later.
	 * 
	 * Two streams, S1 and S2.
	 * The idea is to compare every simplex in S1 to every simplex in S2, figure out
	 * - if they cross
	 * - when will they cross
	 * Then sort the list of crossings.  This is one way to ensure *consistency* of the
	 * intermediate streams.
	 */
	Vector<Transposition<U>> generateTranspositionList(AbstractFilteredStream<U> s1, AbstractFilteredStream<U> s2) {
		
		/*
		 * Unfortunately the filtration index recorded in javaplex streams isn't enough
		 * for this computation.  If two simplices have the same filt value, it gives them
		 * the same filtration index.  This makes it impossible to properly sort the
		 * transpositions so that we don't occasionally end up with an inconsistent
		 * filtration (i.e. simplices entering before their faces).
		 * Maybe it's possible to do this naturally by changing the way we iterate through, or
		 * by including teh simplex dimension in the list sort, but I don't see how without
		 * all this extra bookkeeping (and it really is a waste).
		 */
		THashMap<U,Integer> eindex1 = new THashMap<U,Integer>();
		THashMap<U,Integer> eindex2 = new THashMap<U,Integer>();
		int i=0;
		for (U x : s1) { eindex1.put(x, i++); }
		i=0;
		for (U x : s2) { eindex2.put(x, i++); }

		// as said above, these filtered comparators won't really help
		//Comparator<U> fc1 = new FilteredComparator<U>(s1, this.basisComparator);
		//Comparator<U> fc2 = new FilteredComparator<U>(s2, this.basisComparator);

		Vector<Transposition<U>> transpositions = new Vector<Transposition<U>>();
		
		//System.out.pln("generateTranspositionList():");
		//System.out.println("\nstream1:\n===========");
		//for (U x : s1) { System.out.print(" ("+x+","+s1.getFiltrationIndex(x)+","+eindex1.get(x)+")"); }
		//System.out.println("\nstream2:\n===========");
		//for (U x : s2) { System.out.print(" ("+x+","+s2.getFiltrationIndex(x)+","+eindex2.get(x)+")"); }
		
		for (U x : s1) {
			for (U y : s1) {
				int xi1 = eindex1.get(x);
				int yi1 = eindex1.get(y);
				int xi2 = eindex2.get(x);
				int yi2 = eindex2.get(y);
				//if (fc1.compare(x, y) > 0) {
				if (xi1 > yi1) {
					// check if the simplices cross between frames
					// if the product of compare() is negative then they must change order
					//if (fc1.compare(x, y)*fc2.compare(x, y) < 0) {
					if ((xi1-xi2)*(yi1-yi2) < 0) {
						// figure out *when* they cross
						// some time between 0 and 1
						// we know in stream s1 that x>y
						//System.out.println("\nSimplices must be transposed: ( "+x+" , "+y+" )");
						//System.out.println("Filtration indices are: ( "+s1.getFiltrationIndex(x)+" , "+s1.getFiltrationIndex(y)+
						//		" ) --> ( "+s2.getFiltrationIndex(x)+" , "+s2.getFiltrationIndex(y)+" )");
						double startspan = 1.0*(s1.getFiltrationValue(y) - s1.getFiltrationValue(x));
						double ctime = startspan/(1.0*(s2.getFiltrationValue(x)-s2.getFiltrationValue(y)+startspan));
						if ((ctime < 0.0) || (ctime > 1.0)) {
							System.out.println("ERROR: intersection time is "+ctime);
						}
						//double start_index_diff = 1.0*(s1.getFiltrationIndex(y) - s1.getFiltrationIndex(x));
						//double cindex_time = startspan/(1.0*(s2.getFiltrationIndex(x)-s2.getFiltrationIndex(y)+start_index_diff));
						double cindex_time = (1.0*(yi1-xi1)) / (1.0*(xi2-yi2+yi1-xi1));
						if ((cindex_time < 0.0) || (cindex_time > 1.0)) {
							System.out.println("ERROR: intersection time is "+cindex_time);
						}
						//System.out.println("Transposing:\t t="+ctime+"\t( "+x+" , "+y+" )");
						//System.out.println("Index: \t t="+cindex_time+"\t( "+x+" , "+y+" )");
						transpositions.add(new Transposition<U>(ctime, cindex_time, y, x));
					}
				}
			}
		}
		Collections.sort(transpositions);
		
		System.out.println("Final Transposition List:");
		for (Transposition<U> tr : transpositions) {
			System.out.println(tr.t_index+"\t"+tr.t+"\t"+tr.a+"\t"+tr.b);
		}
		System.out.println("");
		
		return transpositions;
	}

	
	/**
	 * Compute the sequence of transpositions necessary to turn current_stream stream into another
	 * For starters this will just be a simple bubble sort.
	 * In order to do it, we pretty much need the current_stream in an ordered, accessible, mutable structure
	 * Since we can't make an array, let's use a Vector
	 * This is annoying
	 */
	/*
	public Vector<Transposition<U>> streamTranspositionList(AbstractFilteredStream<U> new_stream) {
		// should do some sanity checks here, so make sure the streams have the same elements
		if (new_stream.getSize() != current_stream.getSize()) {
			System.out.println("streamTranspositionList(): streams aren't the same size!  Bailing out");
			return null;
		}
		if (new_stream.getSize() == 0) {
			System.out.println("streamTranspositionList(): new_stream is empty!  Bailing out");
			return null;
		}

		// turn current_stream into a Vector
		Iterator<U> csi = current_stream.iterator();
		Vector<U> csv = new Vector<U>(current_stream.getSize());
		while (csi.hasNext()) {
			csv.add(csi.next());
		}
		
		Vector<Transposition<U>> transpositions = new Vector<Transposition<U>>();
		Iterator<U> nsi = new_stream.iterator();

		// iterate through the entries in new_stream
		// for each entry, find the matching cell in the current stream, and "bubble it up" if necessary
		int match_len = 0;
		while (nsi.hasNext()) {
			U nsu = nsi.next();
			// find the position of nsu in the current stream
			int next_cell_pos = csv.indexOf(nsu, match_len);
			if (next_cell_pos < 0) {
				System.out.println("ERROR: streamTranspositionList(): mismatched streams!");
				System.exit(0);
			}
			if (next_cell_pos > match_len) {
				// "bubble up" the cell in the current stream
				for (int i=next_cell_pos; match_len<i; i--) {
					U a = csv.get(i-1);
					U b = csv.get(i);
					transpositions.add(new Transposition<U>(a, b));
					csv.add(i, a);
					csv.add(i-1, b);
					System.out.println("\ncsv:\t"+csv);
				}
			}
			match_len++;
		}
		return transpositions;
	}
	*/
	
	/**
	 * And this is the heart of vineyard computation.  Update the pairing by transposing the
	 * two (simplices)
	 * 
	 * One pitfall I'm seeing: if you have have several swaps in a single call to addHomotopyFrame(),
	 * then after the first swap the filteredComparator of the stream is no longer valid.  We have
	 * to either find a way to update that thing, or just keep track of and rely on (real)
	 * filtration indices.
	 * 
	 * @param i
	 * @param j
	 */
	public void updatePairingSwap(U a, U b) {
		// verify that the simplices are, in fact, adjacent
		//if ((current_stream.getFiltrationIndex(a)+1) != current_stream.getFiltrationIndex(b)) {
		//	System.out.println("ERROR: simplices "+a+" (index "+current_stream.getFiltrationIndex(a)+") and "+
		//						b+" (index "+current_stream.getFiltrationIndex(b)+") are NOT adjacent in the stream!");
		//	return;
		//}

		System.out.println("updatePairingSwap( "+a+" , "+b+" )");

		/*
		 * Not sure if this matters, but just to make everything match the diagrams 
		 * in my notes, ensure that a < b in the current ordering.
		 * Even if the filtration has had some permutations already, if two simplices are
		 * swapping positions, their relative positions before we do this operation
		 * should be consistent with the original stream.
		 */
		if (this.filteredComparator.compare(a, b) > 0) {
			U tmp = b;
			b = a;
			a = tmp;
			System.out.println("\tupdatePairingSwap( "+a+" , "+b+" )");
		}

		// if the simplices are of different dimensions nothing to do!
		if (current_stream.getDimension(a)  != current_stream.getDimension(b)) {
			// go ahead and do the swap (we're using hash tables and sparse representations...
			// nothing really to do...?)
			// need a way to update the stream storage
			System.out.println("Different dimension exception - returning.");
			return;
		}

		BooleanSparseFormalSum<U> R_a = R.get(a);
		BooleanSparseFormalSum<U> R_b = R.get(b);
		U low_a = low(R_a);
		U low_b = low(R_b);

		THashSet<U> low_a_hash = lowMap.get(a);
		THashSet<U> low_b_hash = lowMap.get(b);
		
		// sanity check
		// if R is in reduced form a and b can only be the low element for one column
		if ((low_a_hash != null) && (low_a_hash.size() > 1)) {
			System.out.println("ERROR: R matrix is not in reduced form!  Bailing out...");
			return;
		}
		if ((low_b_hash != null) && (low_b_hash.size() > 1)) {
			System.out.println("ERROR: R matrix is not in reduced form!  Bailing out...");
			return;
		}
		
		// Case 1
		if ((low_a==null) && (low_b==null)) {
			System.out.println("** Case 1 **");
			// since the column for /a/ in /R/ is zero, we simply set
			// W[a,b]=0 (if it's not already), and PWP is upper-triangular
			W.get(a).remove(b);
			
			// Case 1.1
			// there are cols x, y with lowR[x]=a, lowR[y]=b, and R[y,a]=1
			if ((low_a_hash != null) && (low_b_hash != null)
					&& (low_a_hash.size() > 0) && (low_b_hash.size() > 0)) {
				
				U x = low_a_hash.iterator().next();
				U y = low_b_hash.iterator().next();
				int orderReversalFactor = 1;
				if (orderReversals.containsKey(x+":"+y)) {
					orderReversalFactor = -1;
					System.out.println("About to compare "+x+" to "+y+" but their order has changed!");
				}
				
				// test if R[y,a]=1
				if (R.get(y).getCoefficient(a)) {

					System.out.println("  ** Case 1.1 **");

					// Case 1.1.1 (x < y)
					if (this.filteredComparator.compare(x, y)*orderReversalFactor < 0) {
						System.out.println("    ** Case 1.1.1 **");
						// add col x to col y in R, and update W
						boolean negative_c = R.get(y).getCoefficient(a); // this is unnecessary in GF2, right?
						this.chainModule.accumulate(R.get(y), R.get(x), negative_c);
						this.chainModule.accumulate(W.get(x), W.get(y), negative_c);
						// this preserves low(x)=a, low(y)=b
					}
					// Case 1.1.2 (x > y)
					else {
						System.out.println("    ** Case 1.1.2 **");
						// add col y to col x in R, and update W
						boolean negative_c = R.get(x).getCoefficient(a);
						this.chainModule.accumulate(R.get(x), R.get(y), negative_c);
						this.chainModule.accumulate(W.get(y), W.get(x), negative_c);
						// this swaps the pairing, low(x)=b, low(y)=a
						// update the reverse map
						lowMap.get(a).remove(x); lowMap.get(a).add(y);
						lowMap.get(b).remove(y); lowMap.get(b).add(x);
					}
				}
				
				// Case 1.2a (?????)
				else {
					// the matrices are still reduced, but the pairing /might/ change...

					// not sure if there's a simple shortcut
					// just recompute the lowMap for a and b!
					if (low(R.get(y)).equals(a)) {
						// the pairing changed!
						System.out.println("  ** Case 1.2a **");

						// this swaps the pairing, low(y)=a, low(x)=b
						// update the reverse map
						lowMap.get(a).remove(x); lowMap.get(a).add(y);
						lowMap.get(b).remove(y); lowMap.get(b).add(x);
					} else {
						// the pairing didn't change!
						System.out.println("  ** Case 1.2 (no pairing change) **");
					}
				}
			}

			// Possible Case 1.2b
			else {
				// Case 1.2b (Case 1b in Pairing Change Thm)
				// one of the simplices is positive, but unpaired
				// by necessity, since a<b, the only case we're interested in is when
				// a is the unpaired one, and R[y,a]=1
				// then the pairing is going to change

				if ((low_b_hash != null) && (low_b_hash.size() > 0)) {
					U y = low_b_hash.iterator().next();
					if (R.get(y).getCoefficient(a)) {
						System.out.println("  ** Case 1.2b **");

						// this swaps the pairing, low(y)=a
						// update the reverse map
						//lowMap.get(y).remove(b);  lowMap.get(y).add(a);
						lowMap.remove(b);
						lowMap.put(a, new THashSet<U>());
						lowMap.get(a).add(y);
					}
				}
			}
		}

		// Case 2
		else if ((low_a != null) && (low_b!= null)) {
			System.out.println("** Case 2 **");
			// sanity check
			if (((low_a_hash != null) && (low_a_hash.size() > 0))
					|| ((low_b_hash != null) && (low_b_hash.size() > 0))) {
				System.out.println("ERROR: R matrix is not in reduced form!  Bailing out...");
				return;
			}

			// Case 2.1 (W[a,b]=1)
			if (W.get(a).containsObject(b)) {
				System.out.println("  ** Case 2.1 **");
				// just in case the pairing changes and we need to update lowMap
				U x = low(R.get(a));
				U y = low(R.get(b));
				
				// add row b to row a in W, and update R
				boolean negative_c = W.get(a).getCoefficient(b);
				this.chainModule.accumulate(W.get(a), W.get(b), negative_c);
				this.chainModule.accumulate(R.get(b), R.get(a), negative_c);
				
				int orderReversalFactor = 1;
				if (orderReversals.containsKey(low_a+":"+low_b)) {
					orderReversalFactor = -1;
					System.out.println("About to compare "+low_a+" to "+low_b+" but their order has changed!");
				}
				
				// Case 2.1.1 (lowR[a] < lowR[b])
				if (this.filteredComparator.compare(low_a, low_b)*orderReversalFactor < 0) {
					// do nothing!
					// this preserves the pairing, low(x)=a, low(y)=b
					System.out.println("    ** Case 2.1.1 **");
				}
				
				// Case 2.1.2 (lowR[a] > lowR[b])
				else {
					System.out.println("    ** Case 2.1.2 **");
					// add col b to col a in R
					// add row a to row b in W
					negative_c = R.get(a).getCoefficient(low(R.get(a)));
					this.chainModule.accumulate(R.get(a), R.get(b), negative_c);
					this.chainModule.accumulate(W.get(b), W.get(a), negative_c);
					
					// this swaps the pairing, low(a)=y, low(b)=x
					// update the reverse map
					lowMap.get(x).remove(a); lowMap.get(x).add(b);
					lowMap.get(y).remove(b); lowMap.get(y).add(a);
				}
			}
			
			// Case 2.2 (W[a,b]=0)
			else {
				// do nothing!
				System.out.println("  ** Case 2.2 **");
			}
		}
		
		// Case 3
		else if ((low_a != null) && (low_b == null)) {
			System.out.println("** Case 3 **");
			// Case 3.1 (W[a,b]=1)
			if (W.get(a).containsObject(b)) {
				System.out.println("  ** Case 3.1 **");
				// will definitely change pairing, will need to update lowMap(x)
				U x = low(R.get(a));
				
				// add row b to row a in W, and update R
				boolean negative_c = W.get(a).getCoefficient(b);
				this.chainModule.accumulate(W.get(a), W.get(b), negative_c);
				this.chainModule.accumulate(R.get(b), R.get(a), negative_c);
				
				// now R is not reduced (R[b] was zero, now R[b]=R[a])
				// fix it
				negative_c = R.get(a).getCoefficient(low(R.get(a)));
				this.chainModule.accumulate(R.get(a), R.get(b), negative_c);  // should have R[a]=0 now
				this.chainModule.accumulate(W.get(b), W.get(a), negative_c);

				// this swaps the pairing, low(x)=b
				// update the reverse map
				lowMap.get(x).remove(a); lowMap.get(x).add(b);
			}
			
			// Case 3.2
			else {
				// do nothing!
				System.out.println("  ** Case 3.2 **");
			}
		}
		
		// Case 4
		else if ((low_a == null) && (low_b != null)) {
			System.out.println("** Case 4 **");
			// since the column for /a/ in /R/ is zero, we simply set
			// W[a,b]=0 (if it's not already), and PWP is upper-triangular
			// the pairing is preserved
			W.get(a).remove(b);
		}
		
		// record the fact that we've swapped the order of a and b
		this.orderReversals.put(a+":"+b, true);
		this.orderReversals.put(b+":"+a, true);
	}
	
	/**
	 * This function computes the operation low_A(j) as described in the paper. Note that if
	 * the chain is empty (for example the column contains only zeros), then this function
	 * returns null.
	 * 
	 * @param chain the chain to search
	 * @return  the lowest element of the chain
	 */
	protected U low(BooleanSparseFormalSum<U> chain) {

		//System.out.println("low("+chain+")");
		U maxObject = null;

		U current = null;
		for (Iterator<U> iterator = chain.iterator(); iterator.hasNext(); ) {
			if (this.filteredComparator == null) {
				System.out.println("ERROR: filteredComparator is null!!");
				System.exit(0);
			}
			current = iterator.next();
			if (maxObject == null || this.filteredComparator.compare(current, maxObject) > 0) {
				maxObject = current;
			}
		}

		//System.out.println("\t\t = "+maxObject);
		return maxObject;
	}

	
	/**
	 * Taken initially from BooleanpersistentHomology.java
	 * 
	 * This function implements the alternative pHcol algorithm described in Dmitriy's thesis. It computes the decomposition
	 * D = R * W, where D is the boundary matrix, R is reduced, and W=V^-1 is invertible and upper triangular.
	 * This function returns the pair (R, U). Note that in our implementation, we represent a matrix by
	 * a hash map which maps a generating object to a formal sum which corresponds to a column in the matrix.
	 * Note that this is simply a sparse representation of a linear transformation on a vector space with
	 * free basis consisting of elements of type U.
	 * 
	 * @param stream the filtered chain complex which provides elements in increasing filtration order
	 * @return a ObjectObjectPair containing the matrices R and V
	 */
	private ObjectObjectPair<THashMap<U, BooleanSparseFormalSum<U>>, THashMap<U, BooleanSparseFormalSum<U>>> DRUpHcol(AbstractFilteredStream<U> stream) {

		System.out.println("DRUpHcol( "+stream+" )");
		for (U i : stream) {
			/*
			 * Do not process simplices of higher dimension than maxDimension.
			 */
			if (stream.getDimension(i) < this.minDimension) {
				continue;
			}

			if (stream.getDimension(i) > this.maxDimension + 1) {
				continue;
			}

			// initialize W to be the identity matrix
			W.put(i, this.chainModule.createNewSum(i));

			// form the column R[i] which equals the boundary of the current simplex.
			// store the column as a column in R
			R.put(i, chainModule.createNewSum(stream.getBoundaryCoefficients(i), stream.getBoundary(i)));

			// compute low_R(i)
			U low_R_i = this.low(R.get(i));

			// if the boundary of i is empty, then continue to next iteration since there
			// is nothing to process
			if (low_R_i == null) {
				continue;
			}

			// note: first time through this loop will not execute.
			// only when it gets to a duplicate low_R value, this loop will resolve it.
			// It seems like matchingLowSimplices can only contain one thing on each run through this loop,
			// just by construction...
			// The whole THashSet and iterator thing is unnecessary
			THashSet<U> matchingLowSimplices = lowMap.get(low_R_i);
			while (matchingLowSimplices != null && !matchingLowSimplices.isEmpty()) {
				Iterator<U> iterator = matchingLowSimplices.iterator();
				/**
				 * TODO: Is this the right thing to do???
				 * Ie. should the iterator.next go at the end of the loop?
				 * A: I don't think you need an iterator at all.  You don't iterate it.  You
				 * clobber matchingLowSimplices on each iteration.
				 */
				U j = iterator.next();

				assert (R.get(j).getCoefficient(low_R_i) == true);
				boolean negative_c = R.get(i).getCoefficient(low_R_i);
				// this is different than pHcol because we're representing W=inv(V) as sparse rows
				this.chainModule.accumulate(R.get(i), R.get(j), negative_c);
				this.chainModule.accumulate(W.get(j), W.get(i), negative_c);

				// remove old low_R(i) entry
				//lowMap.get(low_R_i).remove(i);

				// recompute low_R(i)
				low_R_i = this.low(R.get(i));

				matchingLowSimplices = lowMap.get(low_R_i);
			}

			// store the low value in the map
			if (low_R_i != null) {
				if (!lowMap.containsKey(low_R_i)) {
					lowMap.put(low_R_i, new THashSet<U>());
				}
				lowMap.get(low_R_i).add(i);
			}
		}

		// at this point we have computed the decomposition D = R * W
		// we return the pair (R, W)

		return new ObjectObjectPair<THashMap<U, BooleanSparseFormalSum<U>>, THashMap<U, BooleanSparseFormalSum<U>>>(R, W);
	}


}

