package adnascreen;

import java.io.Closeable;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

/**
 * 
 *
 * @param <Key>
 * @param <Value>
 */
public class Cache<Key, Value>{
	/*
	 * This class is a size-limited cache. 
	 * It is implemented as a doubly-linked list plus a hash table. 
	 * Each list node contains a (key, value) pair. 
	 * The hash table maps keys to linked list nodes for fast access.  
	 */
	// doubly linked list, with:
	// head: least recently used
	// tail: most recently used
	private Node leastRecentlyUsed;
	private Node mostRecentlyUsed;
	private int maxValues;
	private Map<Key, Node> nodeKeys;
	private long hits;
	private long misses;
	private long forcedCloses;
	
	public Cache(int maxValues){
		this.maxValues = maxValues;
		nodeKeys = new HashMap<Key, Node>(maxValues);
	}
	
	/**
	 * Check whether value for key is present in cache. 
	 * This will not affect ordering of cache for evicting the least recently used element. 
	 * @param key
	 * @return
	 */
	public boolean inCache(Key key){
		return nodeKeys.containsKey(key);
	}
	
	public int size(){
		return nodeKeys.size();
	}
	
	public long getHits(){
		return hits;
	}
	
	public long getMisses(){
		return misses;
	}
	
	public long getForcedCloses(){
		return forcedCloses;
	}
	
	public Value get(Key key){
		Node found = nodeKeys.get(key);

		if(found == null){
			misses++;
			return null;
		} else {
			hits++;
			updateLastUsed(found);
			return found.value;
		}
	}
	
	public void put(Key key, Value value){
		// check for key in nodeKeys
		Node found = nodeKeys.get(key);
		
		if(found == null){
			found = new Node(key, value);
			while(nodeKeys.size() >= maxValues){
				remove(leastRecentlyUsed);
			}
			nodeKeys.put(key, found);
		} else {
			found.value = value;
		}
		
		updateLastUsed(found);
	}
	
	public Value remove(Key key){
		Node node = nodeKeys.get(key);
		if(node != null){
			remove(node);
			return node.value;
		} else {
			return null;
		}
	}
	
	public void clear(){
		try{
			close();
		} catch (IOException e){
			System.err.println(e);
		}
		leastRecentlyUsed = null;
		mostRecentlyUsed = null;
		nodeKeys.clear();
	}
	
	private Node remove(Node toRemove){
		// remove key, and close if 
		if(toRemove.value instanceof Closeable){
			forcedCloses++;
			try {
				((Closeable) toRemove.value).close();
			} catch (IOException e) {
				e.printStackTrace(System.err);
			}
		}
		
		nodeKeys.remove(toRemove.key);
		// remove node from list
		return removeFromListOnly(toRemove);
	}
	
	private Node removeFromListOnly(Node toRemove){
		if(toRemove != null){
			// end cases
			if(leastRecentlyUsed == toRemove){
				leastRecentlyUsed = leastRecentlyUsed.next;
			}
			if(mostRecentlyUsed == toRemove){
				mostRecentlyUsed = mostRecentlyUsed.previous;
			}
			
			// remove this node from the list, if it was in it
			if(toRemove.previous != null){
				toRemove.previous.next = toRemove.next;
				toRemove.previous = null;
			}
			if(toRemove.next != null){
				toRemove.next.previous = toRemove.previous;
				toRemove.next = null;
			}
			return toRemove;
		} else {
			return null;
		}
	}
	
	private void updateLastUsed(Node used){
		removeFromListOnly(used);
		// startup edge case
		if(leastRecentlyUsed == null){
			leastRecentlyUsed = used;
		}
		
		if(mostRecentlyUsed == null){ // startup edge case
			mostRecentlyUsed = used;
		} else { // non-empty, add to end of list
			mostRecentlyUsed.next = used;
			mostRecentlyUsed = used;
		}
	}
	
	public void close() throws IOException{
		while(leastRecentlyUsed != null){
			Node current = leastRecentlyUsed;
			Value value = current.value;
			if(value instanceof Closeable){
				((Closeable) value).close();
			}
			leastRecentlyUsed = current.next;
		}
		mostRecentlyUsed = null;
		nodeKeys.clear();
	}
	
	private class Node{
		Key key; // required to remove map entries
		Value value;
		Node previous;
		Node next;
		
		public Node(Key key, Value value){
			this.key = key;
			this.value = value;
			this.previous = null;
			this.next = null;
		}
	}
}
