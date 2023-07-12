########################################################################
# UniquePriorityQueue.py
#
# A wrapper around the python PriorityQueue class that only adds items
# to the queue if they haven't been added to the queue before.
########################################################################
from queue import PriorityQueue 
from heapq import heappop, heappush
from numpy import inf

class UniquePriorityQueue(PriorityQueue):
    """Subclass of queue.PriorityQueue to only keep track of unique items (items must be hashable)
    
    Attributes
    ----------
    _maxsize : int
        the maximize size that the queue can become, this is required by the python queue class
    _heapsize : int
        the maximum size that the heap can become
    _expansion_factor : int 
        the maximum number of elements to pop from the queue when get_m() is called 
    """
    def __init__(self, maxsize=0, heapsize=inf, expansion_factor=1):
        """Constructor which sets heapsize and calls PriorityQueue constructor
        
        Parameters
        ----------
        maxsize : int, optional
            the maximize size that the queue can become, this is required by the python queue class
        heapsize : int, optional
            the maximum size that the heap can become
        expansion_factor : int, optional 
            the maximum number of elements to pop from the queue when get_m() is called 
        """
        super().__init__(maxsize)
        self._heapsize = heapsize
        self._expansion_factor = expansion_factor

    def decrease_heapsize(self, size=1):
        """Decrements heap size by a particular value
        
        Parameters
        ----------
        size : int, optional
            how much to decrease the queue size by
        """
        self._heapsize = max(self._heapsize-size, 0)

    def get_m(self):
        """Gets a certain number of elements from the queue determined by the size of the expansion_factor,
        then the queue is rebalanced to maintain the heap property.
        
        Returns
        -------
        list
            a list of items from the queue, the size of the list return is determined by the expansion_factor
        """
        items = [heappop(self.queue) for _ in range(min(self._expansion_factor, len(self.queue)))]
        self._rebalance()
        return items

    def expansion_factor(self):
        """Getter for the expansion factor
        
        Returns
        -------
        int
            the expansion_factor value
        """
        return self._expansion_factor

    def _init(self, maxsize):
        """Override of _init function to create the queue and unique set"""
        self.uniques = set()
        self.queue = []

    def _get(self):
        """Gets a single element from the queue, then the queue is rebalanced to maintain the heap property.
        
        Returns
        -------
        object
            the smallest element from the queue
        """
        item = heappop(self.queue)
        self._rebalance()
        return item

    def _put(self, item):
        """Adds an item to the queue assuming we've never seen the hash of the item before
        
        Parameters
        ----------
        item : object
            an element to add to the queue
        """
        item_hash = hash(item)

        if item_hash not in self.uniques:
            self.uniques.add(item_hash)
            heappush(self.queue, item)

    def _rebalance(self):
        "Rebalances queue to it's heapsize while maintaining its heap property"
        if len(self.queue) >= (self._heapsize - self._expansion_factor):
            self.queue = [heappop(self.queue) for _ in range(min(self._heapsize-self._expansion_factor, len(self.queue)))]
