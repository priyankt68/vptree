#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <stdint.h>
#include <string>
#include <string.h>
#include <math.h>
#include <vector>

#define left(i)   ((i << 1) +1)
#define right(i)  ((i << 1) +2)
#define parent(i) ((i - 1) >> 1)
#define floor_to_power_of_2(x) (1 << ((int) floor(log2(x))))
#define power_of_2(x) pow(2,x)

/*Basic structure for data-points*/
struct point     // CHANGE TO POINT
{
    float x;
    float y;
};

/* Type to represent a list of circles */
typedef struct {
    int n;
    float *x, *y, *r;
} circle_list_t;


/*                              UTILITY FUNCTIONS                          */
/* ________________________________________________________________________*/

int compare (const void * a, const void * b)
{
  return ( *(float*)a - *(float*)b );
}


/*
    Function to print the data-points

*/

void print(point pts[], int n)
{

    for (int i = 0; i < n; ++i)
    {
        std :: cout << pts[i].x << "," << pts[i].y << std::endl ;
    }

}

/* Function to swap the elements */

void swap(float *x, float *y)
{
    float temp = *x;
    *x = *y ;
    *y = temp; 
}


/*                          HELPER FUNCTIONS                                */
/*__________________________________________________________________________*/

void compute_distance( point pts[],int j, int s, float dist[])
{

     //distance array
    int start = j+1;
    int end   = j+s-1;


    for ( int index = start ; index <= end; ++index)
    {   
        dist[index] =  ((pts[index].x - pts[j].x)*(pts[index].x - pts[j].x)) + ((pts[index].y - pts[j].y)*(pts[index].y - pts[j].y)) ;    // dist array contains the  
    }   
}

int partition(float* input, int p, int r)
{
    float pivot = input[r];
    
    std :: cout << " Pivot element is : " << pivot << std :: endl;

    while ( p < r )
    {
        while ( input[p] < pivot )
            p++;
        
        while ( input[r] > pivot )
            r--;
        
        if ( input[p] == input[r] )
            p++;

        else if ( p < r ) 
        {                                             
            float tmp = input[p];                     
            input[p] = input[r];
            input[r] = tmp;
        }
    }
    

    return r;
}

float quick_select(float* input, int p, int r, int k)
{
    if ( p == r ) return input[p];   // if the list just contains one element.
    
    int j = partition(input, p, r);   // partition the elements with all elements

    int length = j - p + 1;     // changing the length 

    if ( length == k )  return input[j];          // we've reached the element.
    
    else if ( k < length ) return quick_select(input, p, j - 1, k);

    else  return quick_select(input, j + 1, r, k - length);
}


/*Partitions the points at indices [j+1, j+s-1] in pts so that the points at indices [j+1,j+pp-1] 
 have distance at most pr from (pts.x[j], pts.y[j]), while the points at indices [j+pp, j+s-1] 
 have distance at least pr. */
void partition_by_distance(point pts[],int j,int s, float pr)
{

    int p = j + 1;
    int r = j + s - 1;

    std :: cout << "Computing distance in partition by distance function" << std :: endl;

    float distance_vector[s-1];

    compute_distance( pts, j , s , distance_vector);   // computing the distance array which contains the distance from the vantage point in the current subtree.

    std :: cout << "COMPUTED DISTANCE ARRAY in partition_by_distance" << std :: endl;

    std :: cout << "dist\t" << "X\t" << "Y\t" << std :: endl;

     for (int index = p; index <= r ; ++index)
    {
        std::cout << distance_vector[index] << "\t" << pts[index].x << "\t" << pts[index].y << "\t" << std::endl ;
    }

    int index_of_radii;

    for (int index = p; index <= r; ++index)
    {
        if( pr == distance_vector[index])
        {
            index_of_radii = index;
        }
    }

    std :: cout << "index of radii : " << index_of_radii <<  std :: endl;
    
    // moving the pivot to the last position, before applying the partitioning algorithm

    swap(  &distance_vector[index_of_radii] , &distance_vector[r]);
    swap(&pts[index_of_radii].x , &pts[r].x);
    swap(&pts[index_of_radii].y , &pts[r].y);

    float pivot = distance_vector[r];

    std :: cout << "Pivot before partitioning" << pivot << std :: endl;
    
    std :: cout << "dist\t" << "X\t" << "Y\t" << std :: endl;

     for (int index = p; index <= r ; ++index)
     {
        std::cout << distance_vector[index] << "\t" << pts[index].x << "\t" << pts[index].y << "\t" << std::endl ;
     }
  
    std :: cout << " Pivot element is : " << pivot << std :: endl;

    while ( p < r )
    {
        while ( distance_vector[p] < pivot )
            p++;
        
        while ( distance_vector[r] > pivot )
            r--;
        
        if ( distance_vector[p] == distance_vector[r] )
            p++;

        else if ( p < r ) 
        {                                              
            swap(  &distance_vector[p] , &distance_vector[r]);
            swap(            &pts[p].x , &pts[r].x          );
            swap(            &pts[p].y , &pts[r].y          );
        }
    }
    
    std::cout << distance_vector[r] << "\t" << pts[r].x << "\t" << pts[r].y << "\t" << std::endl ;
    
    std :: cout << "dist\t" << "X\t" << "Y\t" << std :: endl;

     for (int index = j+1; index <= (j+s-1) ; ++index)
    {
        std::cout << distance_vector[index] << "\t" << pts[index].x << "\t" << pts[index].y << "\t" << std::endl ;
    }

}

/* Chooses a random index r in [j, j+s-1] and swaps pts.x[r] with pts.x[j] and pts.y[r] with 
pts.y[j] */
void random_point(point pts[],int j,int s)
{
    std :: cout << "j = " << j ;                                       // in the first pass, choose the random number between 0 and 15;
    std :: cout << " j+s-1 = " << j+s-1 << std::endl;
   
    int r = (random() % s ) + j;
       
    std :: cout << "random index : " << r << std :: endl;

    /* Swapping elements with the random element selected*/
    std :: cout << "X to be swapped" << pts[r].x << "," << pts[j].x << std::endl;
    std :: cout << "Y to be swapped" << pts[r].y << "," << pts[j].y << std::endl;

    swap(&pts[r].x,&pts[j].x);
    swap(&pts[r].y,&pts[j].y);

    std:: cout << "After the swapped element with the index = " << r << std::endl;
    std:: cout << "Random index selected " << r << std::endl;
  
}

/* 
Finds the (pp-1)st smallest distance among the distances from (pts.x[j], pts.y[j]) of the points 
in pts with indices in [j+1, j+s-1] and stores this distance in pr. Note that this is the same 
as finding the ppth smallest distance among the distances from (pts.x[j], pts.y[j]) among the 
points in pts with indices [j, j+s-1], so it chooses the desired pivot distance. 
*/

void select_by_distance(point pts[], int j, int s, int pivot, float *pr)  
{

  //  std :: cout << "Pivot elements here  : " << pivot << std::endl;
    std :: cout << " j+1  : " << j+1 << std::endl;
    std :: cout << " j+s-1  : " << j+s-1 << std::endl;
    std :: cout << " s-1:  " << s-1 << std::endl;
    //distance array
    int start = j+1;
    int end   = j+s-1;

    int dist_length =  (end - start) + 1;  // (s-1 )

    float dist[s-1];   // distance array which will hold the distances.
   
    compute_distance( pts, j , s , dist );   // computing the distance array which contains the distance from the vantage point in the current subtree.

    float temp_dist[s-1];
    
     for (int index = start ; index <= end; ++index)
    {
        temp_dist[index] = dist[index] ;
    }    

    std :: cout << "COMPUTED DISTANCE ARRAY" << std :: endl;

     for (int index = start ; index <= end; ++index)
    {
        std::cout << dist[index] << std::endl ;
    }

    int index_of_radii;   // index of the radii

   *pr = quick_select(temp_dist, j+1, j+s-1 , pivot);

    
    std :: cout << pivot << "th smallest radii " << *pr << std :: endl;

     for (int index = start ; index <= end; ++index)
    {
        if(*pr == dist[index])
        {
            index_of_radii = index;
        }
    }

    std :: cout << "Partitioned distance array : " << std :: endl;

     for (int index = start; index <= end; ++index)
    {
        std::cout << temp_dist[index] << std::endl ;
    }

    std :: cout << index_of_radii-1 << " <- index of the smallest radii ->" << *pr << std :: endl;


   

}

 /*Partitions the points at indices [j+1, j+s-1] in pts so that the points at indices [j+1,j+pp-1] 
 have distance at most pr from (pts.x[j], pts.y[j]), while the points at indices [j+pp, j+s-1] 
 have distance at least pr. */
 


void build_vp_tree (point pts[], int n)
{
    int sum=0;
    int fringe=0;  
    
    /* The VP tree structure */
    circle_list_t vp;

    /* n rounded down to the closest power of 2, also the number of leaves
       in the left subtree of the root if that subtree is in fact full */
    int t = floor_to_power_of_2(n);

    /* Number of leaves in the "fringe" (last level) of the tree. */
    int f = n - t + 1;

    /* Displaying the fring-size */
    std :: cout << "fringe size for " << n << " data-points is "<< f << std :: endl ;

    /* Index of current node in vp tree data structure with the structure type as circle_list_t*/
    int i, ii;

    /* Difference between i and index of leftmost node on current level */
    int di;
    
    /* Starting position of subarray of pts holding the points of current subtree */
    int j;

    /* Index of first fringe leaf of current subtree */
    int l;
    
    /* Size of current subtree */
    int s;
    
    /* Pivot position in sorted subarray */
    int pp;
    
    /* Pivot point + radius */
    float px, py, pr;

    /* Allocate space for the arrays in the VP tree */

    vp.n = n;
    vp.x = (float*) malloc(((vp.n << 1) + (vp.n >> 1)) * sizeof(float)); // |____________(space fot x and y i.e.(vp.n * 2))_________|________space to store floor(n/2)__________________|
    vp.y = vp.x + vp.n;                                                 // pointing to the point after space for x has been given                   
    vp.r = vp.y + vp.n;                                                 // pointing to the point after the space allocated for the data points x,y have been covered followed by floor(n/2)
    

    for (i = 0; i < n ;) // check for n. (TO DO)                      // i and ii tracks the current point in the vantage point tree.
    {

        std :: cout << "Starting for the next node" << std::endl;

        /* Size of full subtree at current level */
        s = (t << 1) - 1;                                             //   n = 7 , t = 4 , s = 7 

        //std :: cout << "Size of the full subtree: " << s << std :: endl;

        /* Pivot position within a full subtree */
        pp = t - 1;                                                   // pp = 3

        std::cout << "Pivot inside build process :"  << pp << std::endl ;

        /* We're starting at the leftmost leaf */
        l = 0;                                                        // l = 0 

        /* We're starting at the leftmost node */
        di = 0;                                                       // di = 0
        
        /* Calculate starting position of points array of leftmost subtree */
        for (j = 0, ii = i; ii; ++j, ii = ii >> 1);                          // ii gets multiplied by two every time the iteration proceeds.       
        
        std :: cout << "Value of J before while loop even begins" << j << std :: endl; 


        /* Loop over subtrees at current level */
        while (j < n) 
        {

            /* The transition subtree that may not be full needs special treatment */
            if (l < f && (l + t) >= f)
             {

                /* Calculate subtree size */
                s = t - 1 + f - l;

                /* Reduce fringe size for the remainder of this level (and for the next) */
                t = t >> 1;
                
                /* Calculate position of pivot */
                if ((f - l )< t) 
                    pp = t + f - l-1 ; //todo  (t+f-l-1 )

            }

            /* Partitioning is needed only if we have stuff to partition */
            if (s > 1) 
            {
            
                
                random_point(pts, j, s);                /* Chooses a random index r in [j, j+s-1] and swaps pts.x[r] with pts.x[j] and pts.y[r] with pts.y[j] */

                                                                           
                
                /*
                   Finds the (pp-1)st smallest distance among the distances from (pts.x[j], pts.y[j])
                   of the points in pts with indices in [j+1, j+s-1] and stores this distance in pr.
                   Note that this is the same as finding the ppth smallest distance among the distances
                   from (pts.x[j], pts.y[j]) among the points in pts with indices [j+1, j+s-1], so it
                   chooses the desired pivot distance. */

                select_by_distance(pts, j, s, pp-1, &pr);                    // TODO (priyank) : check if its pp-1 or pp to be the pivotal element.           
                                                        
               
                                                                        

                /* 
                   Partitions the points at indices [j+1, j+s-1] in pts so that the points at indices
                   [j+1,j+pp-1] have distance at most pr from (pts.x[j], pts.y[j]), while the points
                   at indices [j+pp, j+s-1] have distance at least pr. */

                partition_by_distance(pts, j, s, pr);  

                /* Store the pivot radius in radius array of vp */
                vp.r[i] = pr;
                
            }

            //Once the partition has been done, then save the elements in the vantage-point tree data structure.

            /* Store current pivot in vp */
           // vp.x[i] = pts.x[j];                  // Check for the assignments which are made here and way they are made here.
           // vp.y[i] = pts.y[j];                  // Same here
            vp.x[i] = pts[j].x;
            vp.y[i] = pts[j].y;
            ++i;
            ++di;

            /* Calculate index in pts array for next subtree */
            for (j += s, ii = di; !(ii & 1); ++j, ii = ii >> 1);
            
            /* If we're in the transition subtree, reset s and pp to default values for the
               remainder of this level and for the next */
            if (l < f && (l + t) >= f) 
            {

                s = (t << 1) - 1;
                pp = t - 1;

            }

            /* Update fringe leaf index.  This index is meaningless for subtrees after the
               transition subtree, but the tests still succeed */
            l += t;
            
        }
    }
}



/*       MAIN()                   */
/*________________________________*/

int main()
{


    /* Extracting data from the .csv files to enable the build process */

    std::vector<point> points; // vector to temporarily store points
    std :: ifstream infile;
    infile.open("../data/data.txt"); // open file
    if(infile)
    {
        std :: string s="";
        point temp_point;
        while(infile)
        {
        getline(infile,s);
        char* pEnd;
        temp_point.x = (strtod (s.c_str(), &pEnd)) ;
        temp_point.y =(strtod (pEnd, NULL));
        points.push_back(temp_point);
        }
                
    }

    /*Total number of points in the data*/
    int n = points.size() - 1;  

    std :: cout << "Number of data points to process" << n << std::endl ;  

    /*Declaring the structure of arrays for the data points to be stored.*/
    point pts[n];
    int i=0;
    for(std::vector<point>::iterator it=points.begin(); it!=points.end(); ++it)
        {   
            pts[i].x =  it->x;
            pts[i].y =  it->y;
            i++;
        }

    /* Printing data */
    print(pts,n);

    /* Calling the build process of the tree */
    build_vp_tree(pts,n);



    return 0;
}