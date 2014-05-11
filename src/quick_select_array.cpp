/*
- The QuickSelect algorithm quickly finds the k-th smallest element of an unsorted array of n elements.
- It is an O(n), worst-case linear time, selection algorithm. A typical selection by sorting method would need atleast O(n log n) time.
- This algorithm is identical to quick sort but it does only a partial sort, since we already know which partition our desired element lies as the pivot is in final sorted position.
*/

#include <iostream>
using namespace std;

// A simple print function
void print(int *input)
{
    for ( int i = 0; i < 5; i++ )
        cout << input[i] << " ";
    cout << endl;
}


/* Function to swap the elements */
void swap(float *x, float *y)
{
    float temp = *x;
    *x = *y ;
    *y = temp; 
}

int partition(int* input, int p, int r)
{
    int pivot = input[r];
    
    while ( p < r )
    {
        while ( input[p] < pivot )
            p++;
        
        while ( input[r] > pivot )
            r--;
        
        if ( input[p] == input[r] )
            p++;
        else if ( p < r ) {
            int tmp = input[p];
            input[p] = input[r];
            input[r] = tmp;
        }
    }
    
    return r;
}

int quick_select(int* input, int p, int r, int k)
{
    if ( p == r ) return input[p];   // if the list just contains one element.

    int j = partition(input, p, r);   // partition the elements with all elements

    int length = j - p + 1;

    if ( length == k ) return input[j];          // we've reached the element.

    else if ( k < length ) return quick_select(input, p, j - 1, k);

    else  return quick_select(input, j + 1, r, k - length);
}

int main()
{
    
    int a[] = { 1, 4, 3, 5, 2, 9,7,6};

    cout << "5th smallest element :  " << quick_select(a, 0, 7, 5) << endl;

    int b[] = {1};
    cout << "1st smallest element :  " << quick_select(b, 0, 0, 1) << endl;

    int c[] = {1,2,4};
    cout << "2nd smallest element :  " << quick_select(c, 0, 2 , 2) << endl;

    int d[] = {0,1,2,4,5,21,-1,-11,-22};                                      // handling negative queries as well. In vp tree, co-ordinates may be negative, in case geo-spatial queries.
    cout << " 4th smallest element :  " << quick_select(d, 0, 9 , 4) << endl;
}