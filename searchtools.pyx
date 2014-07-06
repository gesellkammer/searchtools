#cython: boundscheck=False
#cython: embedsignature=True
#cython: wraparound=False
#cython: infer_types=True

cdef extern from "math.h":
    cdef double floor(double n)
    
def binsearch(seq, t):
    """
    Search the sorted seq for item t by using a binary search.
    It only works with 1-dimensional sequences.
    
    If found, return the index in seq so that seq[index] == t
    If not found, return a negative number where:
    
        i = binary_search(seq, t)
        if i < 0:
            assert seq[abs(i)-1] < t < seq[abs(i)]
        seq.insert(abs(i), t)   # the seq remains sorted
        assert is_sorted(seq)
    
    If t is in seq, then binsearch(seq, t) is the same as
        * numpy.searchsorted(seq, t)
        * bisect.bisect_left(seq, t)

    If t is NOT part of seq, then binsearch(seq, t) is:
        * -numpy.searchsorted(seq, t)
        * -bisect.bisect_left(seq, t)
        
    Thus, binsearch can be used to check if t is a part of seq 
    (it will be a positive number) and to check where t would
    be inserted in seq.

    """
    if isinstance(t, (float, int, long)):
        if isinstance(seq, list):
            return binsearch_floatlist(seq, t)
        else:
            return binsearch_float(seq, t)
    else:
        return binsearch(seq, float(t))

cpdef int binsearch_float(seq, double t):
    cdef int mininum, maximum, m
    cdef double seq_min, seq_max, seq_at_m, tmp
    cdef int seq_len
    minimum = 0
    seq_len = len(seq)
    maximum = seq_len - 1
    if maximum < 0:
        return -1
    seq_min = seq[0]
    seq_max = seq[seq_len - 1]
    if t > seq_max:
        return -(maximum + 1)
    if t < seq_min:
        return 0
    while 1:
        #m = min + int( (max - min)  * ((t - seq_min) / (seq_max - seq_min)) + 0.999999999999)
        tmp = (t - seq_min) / (seq_max - seq_min)
        tmp = tmp * (maximum - minimum)
        tmp = floor(tmp + 0.999999999999)
        m = minimum + <int>tmp
        seq_at_m = seq[m]
        if seq_at_m > t:
            maximum = m - 1
            seq_max = seq[maximum]
            if seq_max < t:
                return -m
        elif seq_at_m < t:
            minimum = m + 1
            seq_min = seq[minimum]
            if seq_min > t:
                return -minimum
        else:
            return m

cpdef int binsearch_floatlist(list seq, double t):
    cdef int mininum, maximum, m
    cdef double seq_min, seq_max, seq_at_m, tmp
    cdef int seq_len
    minimum = 0
    seq_len = len(seq)
    maximum = seq_len - 1
    if maximum < 0:
        return -1
    seq_min = seq[0]
    seq_max = seq[seq_len - 1]
    if t > seq_max:
        return -(maximum + 1)
    if t < seq_min:
        return 0
    while 1:
        #m = min + int( (max - min)  * ((t - seq_min) / (seq_max - seq_min)) + 0.999999999999)
        tmp = (t - seq_min) / (seq_max - seq_min)
        tmp = tmp * (maximum - minimum)
        tmp = floor(tmp + 0.999999999999)
        m = minimum + <int>tmp
        seq_at_m = seq[m]
        if seq_at_m > t:
            maximum = m - 1
            seq_max = seq[maximum]
            if seq_max < t:
                return -m
        elif seq_at_m < t:
            minimum = m + 1
            seq_min = seq[minimum]
            if seq_min > t:
                return -minimum
        else:
            return m

def binsearch2 (seq, double t, int column=0):
    """
    binsearch2 (seq, double t, int column)
    
    search the seq for item t at the column specified
    
    use if seq is a multidimensional seq
    
    if found, return the index in seq so that seq[index] = t
    if not found, return a negative number where
    i = abs(binary_search(seq, t))
    seq[i-1] < t < seq[i]
    
    example:
    
    seq = ((0,0), (41,5), (7,23), (9, 43))
    
    binsearch2 (seq, 6.2, 1) --> returns -2, because seq[2][1] > 6.2
    binsearch2 (seq, 23, 1) --> returns 2, because seq[2][1] = 23
    """
    cdef int minimum 
    cdef int maximum 
    cdef double seq_min 
    cdef double seq_max 
    cdef int m
    cdef double seq_at_m
    if not seq:
        return -1
    minimum = 0
    maximum = len(seq) - 1
    if maximum < 0:
        return -1
    
    seq_min = seq[0][column]
    seq_max = seq[maximum][column]
    
    if t > seq_max:
        return -(maximum + 1)
    if t < minimum:
        return 0.0
   
    while 1:
        m = minimum + int( (maximum - minimum)  * ((t - seq_min) / (seq_max - seq_min)) + 0.999999999999)
        seq_at_m = seq[m][column]
        if seq_at_m > t:
            maximum = m - 1
            seq_max = seq[maximum][column]
            if seq_max < t:
                return -m
        elif seq_at_m < t:
            minimum = m + 1
            seq_min = seq[minimum][column]
            if seq_min > t:
                return -minimum
        else:
            return m
            
def nearest(seq, x):
    """
    return the element of seq which is nearest to x
    """
    cdef int i = binsearch(seq, x)
    cdef int n = len(seq)
    if i >= 0:
        return seq[i]
    i = abs(i)
    if i >= n:
        return seq[n - 1]
    elif i == 0:
        return seq[0]
    else:
        obj0 = seq[i - 1]
        obj1 = seq[i]
        if abs(x - obj0) < abs(x - obj1):
            return obj0
        return obj1
    
def nearest_index(seq, x):
    """
    the same as nearest but returns the index of the element
    """
    cdef int i = binsearch(seq, x)
    cdef int n = len(seq)
    if i >= 0:
        return i
    i = abs(i)
    if i >= n:
        return n - 1
    elif i == 0:
        return 0
    else:
        obj0 = seq[i - 1]
        obj1 = seq[i]
        if abs(x - obj0) < abs(x - obj1):
            return i - 1
        return i
                    