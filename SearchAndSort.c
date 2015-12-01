// Searching and sorting.  Michael J. Fairchild.  http://mikef.org/

#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <math.h>
#include <time.h>

int linear_search(int *, int, int, int);
int binary_search(int *, int, int, int);

void insertion_sort(int *, int, int);
void bubble_sort(int *v, int, int);
void selection_sort(int *v, int, int);
void merge_sort(int *, int, int);
int* _merge_sort(int *, int *, int, int);
void _merge(int *, int *, int, int, int);

void heap_sort(int *, int, int);
typedef struct {
	int *v;
	int length;
	int heap_size;
} heap;
void build_max_heap(heap *);
void max_heapify(heap *, int);
void quicksort(int *, int, int);
int quicksort_partition(int *, int, int);


/* Input:
   v: an array
   l: the first (left) index within v to search
   r: the last (right) index within v to search; r>=l
   x: the search key
  
   Output:
   This function returns -1 if the x does not belong to the
   array v[l]..v[r].  Otherwise, it returns the index i such that
   l <= i <= r and v[i]==x.
  
   Performance: Worst=O(n), Best=O(1)
*/
int linear_search(int *v, int l, int r, int x)
{
	int i, j=-1;
	for (i = l; i <= r; ++i) {
		if (v[i] == x) {
			j=i;
			break;
		}
	}
	return j;
}

/* Input:
   v: an array
   l: the first (left) index within v to search
   r: the last (right) index within v to search; r>=l
   x: the search key
   The elements v[l]..v[r] MUST be sorted in increasing order.
  
   Output:
   This function returns -1 if the x does not belong to the
   array v[l]..v[r].  Otherwise, it returns the index i such that
   l <= i <= s and v[i]==x.
  
   Performance: Worst=O(lg(n)), Best=O(1)
*/
int binary_search(int *v, int l, int r, int x)
{
	int i, j=-1, c; // i=current index, j=return value, c=comparison
	while (l <= r) {
		i=(l+r)/2;
		c = v[i]-x;
		if (c > 0)
			r = i-1;
		else if (c < 0)
			l = i+1;
		else {
			j=i;
			break;
		}
	}
	return j;
}

/* Input:
   v: an array
   l: the first (left) index within the array to sort
   r: the last (right) index within the array to sort
  
   Output:
   The vector v is updated so that v[l]..v[r] is sorted in increasing order.
  
   Peformance: Worst=O(n^2), Best=O(n)
*/
void insertion_sort(int *v, int l, int r)
{
	int i, j, vi;
	for (i = l+1; i <= r; ++i) {
		vi = v[i];
		for (j = i-1; j >= l && v[j] > vi; --j)
			v[j+1]=v[j];
		v[j+1]=vi;
	}
}

/* Input:
   v: an array
   l: the first (left) index within the array to sort
   r: the last (right) index within the array to sort
  
   Output:
   The vector v is updated so that v[l]..v[r] is sorted in increasing order.
  
   Performance: Worst=O(n^2), Best=O(n)
*/
void bubble_sort(int *v, int l, int r)
{
	int i,j,k,t,swaps=1;
	for (k=r; swaps > 0; --k) {
		swaps=0;
		for (i=l; i<k; ++i) {
			j=i+1;
			if (v[i]>v[j]) {
				t = v[i];
				v[i] = v[j];
				v[j] = t;
				++swaps;
			}
		}
	}
}

/* Input:
   v: an array
   l: the first (left) index within the array to sort
   r: the last (right) index within the array to sort
  
   Output:
   The vector v is updated so that v[l]..v[r] is sorted in increasing order.
  
   Performance: Worst=Best=O(n^2)
*/
void selection_sort(int *v, int l, int r)
{
	int i,j,k,t,min;
	for (i=l; i<r; ++i) {
		k=i;
		min=v[i];
		for (j=i+1; j<=r; ++j) {
			if (v[j] < min) {
				k=j;
				min=v[j];
			}
		}
		t=v[i];
		v[i]=min;
		v[k]=t;
	}
}

/* Input:
   v: an array
   l: the first (left) index within the array to sort
   r: the last (right) index within the array to sort
  
   Output:
   The vector v is updated so that v[l]..v[r] is sorted in increasing order.
   This sort allocates O(n) additional memory storage
  
   Performance: Worst=Best=O(n*lg(n))
*/
void merge_sort(int *v, int l, int r)
{
	int *w = _merge_sort(v, NULL, l, r);
	if (w != NULL)
		free(w);
}

/* Input:
   v: an array
   w: a temporary array (either NULL or the size of v) for use within _merge
   l: the first (left) index within the array to sort
   r: the last (right) index within the array to sort
  
   Output:
   The vector v is updated so that v[l]..v[r] is sorted in increasing order.
   This sort allocates O(n) additional memory storage
  
   Performance: Worst=Best=O(n*lg(n))
*/
int* _merge_sort(int *v, int *w, int l, int r)
{
	int n = r-l;
	if (n == 0)
		return NULL;
	else if (n <= 64) { // Insertion sort is fast enough for small arrays
		insertion_sort(v, l, r);
		return NULL;
	}
	else {
		if (w == NULL) {
			if ((w = (int*)malloc(n*sizeof(int))) == NULL) {
				fprintf(stderr, "Cannot allocate memory.\n");
				exit(-1);
			}
		}
		int m = (r+l)/2;
		_merge_sort(v,w,l,m);
		_merge_sort(v,w,m+1,r);
		_merge(v,w,l,m,r);
		return w;
	}
}

/* Input:
   v: an array
   w: a temporary array for merging, at least as big as n:=r-l+1
   l: first (left) index
   m: middle (middle) index
   r: last (right) leftindex
   Must have 0 <= l <= m < r
   The subarrays v[l]..v[m] and v[m+1]..v[r] must already by sorted in
   increasing order.  This does *not* sort in place; it requires O(n)
   memory storage.
   
   Output:
   The subarray v[l]..v[r] is sorted in increasing order
   
   Performance: Worst=Best=O(n)
*/
void _merge(int *v, int *w, int l, int m, int r)
{
	// Check arguments and create temporary vector w
	int n = r-l+1;
	if (l > m || m >= r || l < 0 || m < 0 || r < 0) {
		fprintf(stderr, "Must have 0 <= l <= m < r\n");
		exit(-1);
	}
	
	// Put v[l]..v[m] into w from the left; put v[m+1]..v[r]
	// into w in reverse order from the right.
	int *vl = &v[l], *vm = &v[m], *vr = &v[r];
	int k = 0;
	while (vl <= vm)
		w[k++] = *vl++;
	while (vr > vm)
		w[k++] = *vr--;
	
	// Merge from w back into v so that v[l]..v[r] is sorted.
	int *wl = &w[0], *wm = &w[m-l], *wr = &w[n-1];
	vl = &v[l];
	k=0;
	for (; wl <= wm && wr > wm; ++k)
		v[l+k] = (*wl <= *wr ? *wl++ : *wr--);
	for (; wl <= wm; ++k)
		v[l+k] = *wl++;
	for (; wr > wm; ++k)
		v[l+k] = *wr--;

}


/* Input:
   v: an array
   l: the first (left) index within the array to sort
   r: the last (right) index within the array to sort
  
   Output:
   The vector v is updated so that v[l]..v[r] is sorted in increasing order.
   This is an in-place sort.
  
   Performance: Worst=Best=O(n*lg(n))
*/
void heap_sort(int *v, int l, int r)
{
	int n = r-l+1;
	if (n <= 64) {
		insertion_sort(v, l, r); // insertion sort is fast enough for small arrays
		return;
	}
	
	int i, t;
	heap h;
	h.v = (l == 0 ? v : &v[l]);
	if ((h.length = n) <= 1)
		return;
	build_max_heap(&h);
	for (i = h.length-1; i > 0; --i) {
		t = v[0];
		v[0] = v[i];
		v[i] = t;
		h.heap_size--;
		max_heapify(&h, 0);
	}
}

/* Input:
   h: a heap
  
   Output:
   The heap h is a satisfies the max-heap property.  That is,
   the value of each node is at most the value of its parent node.
  
   Performance: Worst=O(n)
*/
void build_max_heap(heap *h)
{
	h->heap_size = h->length;
	int i;
	for (i = h->length/2 - 1; i >= 0; --i)
		max_heapify(h, i);
}

/* Input:
   h: a heap
   i: an index into the heap
  
   The subtree rooted at node i need not be a max heap, but it is
   assumed that the subtrees rooted at the left child of node i and
   the right child of node i satisfy the max-heap property.
  
   Output: The subtree rooted at node i is a max heap.
  
   Performance: Worst=O(lg(n))
*/
void max_heapify(heap *h, int i)
{
	int fl = (h->heap_size)>>1; // index of first leaf is heap_size/2
	while (i < fl) {
		int l = (i << 1) | 0x01; // left child index l = 2*i+1
		int r = (i + 1) << 1; // right child index r = 2*i+2 = 2*(i+1)
		int t, si=i; // t = temporary, si = swap index
		if (l < h->heap_size && h->v[l] > h->v[i])
			si = l;
		if (r < h->heap_size && h->v[r] > h->v[si])
			si = r;
		if (si != i) {
			t = h->v[i];
			h->v[i] = h->v[si];
			h->v[si] = t;
			i = si;
		} else
			break;
	}
}

/* Input:
   v: an array
   l: the first (left) index within the array to sort
   r: the last (right) index within the array to sort
  
   Output:
   The vector v is updated so that v[l]..v[r] is sorted in increasing order.
   This algorithm sorts in place.
  
   Performance: Worst=O(n^2), Best=O(n*lg(n))
*/
void quicksort(int *v, int l, int r)
{
	if (l < r) {
		int q = quicksort_partition(v,l,r);
		quicksort(v,l,q-1);
		quicksort(v,q+1,r);
	}
}

/* Input:
   v: an array
   l: the first (left) index within the array to sort
   r: the last (right) index within the array to sort
  
   Output:
   The array is partitioned into two (possibly empty)
   subarrays v[l]..v[i-1] and v[i+1]..v[r] such that each
   element of v[l]..v[i-1] is less than or equal to v[i], and
   v[i] is less than or equal to each v[i+1]..v[r].  The
   partitioning index i is returned.
*/
int quicksort_partition(int *v, int l, int r)
{
	int x = v[r], i=l-1, j, t;
	for (j = l; j <= r-1; j++) {
		if (v[j] <= x) {
			i++;
			t = v[i];
			v[i] = v[j];
			v[j] = t;
		}
	}
	t = v[i+1];
	v[i+1]=v[r];
	v[r]=t;
	return i+1;
}

void print_vector(int *v, int n)
{
	int i;
	for (i = 0; i < n; ++i)
		printf("%d,",v[i]);
	printf("\n");
}

int main(void)
{
	// Array of function pointers to sort functions
	typedef void (*sort_fptr)(int*,int,int);
	sort_fptr fptrs[] = {&insertion_sort, &bubble_sort, &selection_sort, &heap_sort, &merge_sort, &quicksort, NULL};

	// Output file to store timings
	FILE *fp;
	if ((fp = fopen("SortTest.txt","w")) == NULL) {
		fprintf(stderr, "Cannot open output file.\n");
		exit(-1);
	}
	
	// Loop over different array sizes (powers of 2)
	int i, j, n, *v, *w;
	for (n = 2; n <= 1024*1024*8; n*=2) {
		// Initialize arrays v,w of size n with random elements
		v = (int*)malloc(n*sizeof(int));
		w = (int*)malloc(n*sizeof(int));
		if (v == NULL || w == NULL) {
			fprintf(stderr, "Cannot allocate memory.\n");
			exit(-1);
		}
		for (i = 0; i < n; ++i)
			v[i] = rand() % n;
		// w is a copy of v in order to restore v between successive sorts
		memcpy(w, v, n*sizeof(int));
		
		// Now loop over each function pointer (sort function) and invoke it
		printf("%10d",n);
		fprintf(fp, "%10d", n);
		clock_t start, stop, diff;
		sort_fptr fptr = NULL;
		j = 0;
		for (fptr = fptrs[j]; fptr != NULL; fptr = fptrs[++j]) {
			if (n > 32*1024 && j <= 2) {
				// Don't even bother with insertion, bubble, or selection sort if n>32k
				diff = 0;
			}
			else {
				start = clock();
				fptr(v,0,n-1); // Invoke the sort function
				stop = clock();
				diff = stop-start;
				memcpy(v, w, n*sizeof(int)); // Restore v from w for next sort
			}
			printf("%10f",(float)diff/CLOCKS_PER_SEC);
			fprintf(fp, "%10f", (float)diff/CLOCKS_PER_SEC);
		}
		printf("\n");
		fprintf(fp, "\n");
		
		free(v);
		free(w);
	}

	// Cleanup
	fclose(fp);
}
