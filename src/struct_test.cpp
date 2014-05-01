#include <iostream>
#include <fstream>
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
struct point_list
{
    std::string coordinates;
    int n;
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

void print(point_list pts[], int n)
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


/*Function to calculate the squared distances between two points */


float distance( const point_list& p1, const point_list& p2 )
{
    float a = (p1.x-p2.x);
    float b = (p1.y-p2.y);
    return sqrt(a*a+b*b);
}


/*                          HELPER FUNCTIONS                                */
/*__________________________________________________________________________*/

/* Chooses a random index r in [j, j+s-1] and swaps pts.x[r] with pts.x[j] and pts.y[r] with 
pts.y[j] */

void random_point(point_list pts[],int j,int s)
{
    std :: cout << "j = " << j ;                                       // in the first pass, choose the random number between 0 and 15;
    std :: cout << " j+s-1 = " << j+s-1 << std::endl;
   // int r = rand() % ( 2*j + s  );    // selecting a random element between
   //int r=0;    //TODO check using random()

    int r = (random() % s ) + j;
    // (rand % s) + j  TO DO && CORRECT  (DONE)
    
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
points in pts with indices [j, j+s-1], so it chooses the desired pivot distance. */


void select_by_distance(point_list pts[], int j, int s, int pivot, float *pr)  // check the last parameter here and      UPDATE: (Seems to be correct as of now, but give a thorough check)
                                                                             //  the way it's going to be accessed.  UPDATE : It's correct and is ready to go.
{

    std :: cout << "Pivot elements here  : " << pivot << std::endl;
    std :: cout << " j+1  : " << j+1 << std::endl;
    std :: cout << " j+s-1  : " << j+s-1 << std::endl;
    //distance array
    float dist[s];
    int i=0;
    for (int index = j+1; index <= j+s-1; ++index)
    {   
       
        dist[i++] = sqrt( ( pow(pts[index].x - pts[j].x,2) ) + ( pow(pts[index].y - pts[j].y,2) ) ) ;    // dist array contains the 
                 
    }

    float temp[s];
    /* Checking for the distance */

    for (int k = 0; k < i; ++k)
    {
        temp[k] =  dist[k];
        //std::cout << dist[k] << std::endl ;
    }

    std::sort (temp, temp+i);

    std::cout << "After sorting the distances" << std::endl ;

    for (int k = 0; k < i; ++k)
    {
        std::cout << temp[k] << std::endl ;
    }
    
    std::cout << "Original distances" << std::endl ;

      for (int k = 0; k < i; ++k)
    {
        std::cout << dist[k] << std::endl ;
    }
    std::cout << "\n" << pivot << "th smallest element is : " << temp[pivot-1] << std::endl ;
   
    *pr = temp[pivot-1];

    //free(temp);
    
    int start = j+1;
    int end = j+s-1;

    std::cout << "start for partitioning : " <<  start << "end for partitioning" << end << std::endl;

    while (dist[start] <= *pr) 
        start += 1;
    while (*pr < dist[end]) 
        end -= 1 ;

    while (start < end)
    {
        swap ( &pts[start].x , &pts[end].x );
        swap ( &pts[start].y , &pts[end].y );
        swap ( &dist[start]  , &dist[end]  );

        start += 1;
        end   -= 1;

        while (dist[start] <= *pr && start < end) 
            start += 1;
        while ( *pr < dist[end] && start < end)
            end   -= 1;


    }

    std :: cout << "partitioning index : (end) "<< end << std :: endl;
    std :: cout << "pivot elemtn" << pivot << std::endl;

    std :: cout << "Distance\t" << "X\t" << "Y\t" << std :: endl;
    for (int k = 0; k < i; ++k)
    {
        std::cout << dist[k] <<"\t" << pts[k].x <<"\t" << pts[k].y << std::endl ;
    }
    

}

 /*Partitions the points at indices [j+1, j+s-1] in pts so that the points at indices [j+1,j+pp-1] 
 have distance at most pr from (pts.x[j], pts.y[j]), while the points at indices [j+pp, j+s-1] 
 have distance at least pr. */
void partition_by_distance(point_list pts[],int j,int s, float pr)
{

    //std::cout << "inside partition_by_distance : "  << pr << std::endl ;

}

void build_vp_tree (point_list pts[], int n)
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
    std :: cout << "fringe size for " << n << "data-points is "<< f << std :: endl ;

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
    

    for (i = 0; i < n - 1;) // check for n. (TO DO)                      // i and ii tracks the current point in the vantage point tree.
    {

        std :: cout << "Starting for the next node" << std::endl;

        /* Size of full subtree at current level */
        s = (t << 1) - 1;                                             //   n = 15 , t = 8 , s = 15 

        /* Pivot position within a full subtree */
        pp = t - 1;                                                   // pp = 7
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
            if (l < f && l + t >= f)
             {

                /* Calculate subtree size */
                s = t - 1 + f - l;

                /* Reduce fringe size for the remainder of this level (and for the next) */
                t = t >> 1;
                
                /* Calculate position of pivot */
                if (f - l < t) pp = t + f - l ; //todo  (t+f-l-1 )

            }

            /* Partitioning is needed only if we have stuff to partition */
            if (s > 1) 
            {
            
                
                random_point(pts, j, s);                /* Chooses a random index r in [j, j+s-1] and swaps pts.x[r] with pts.x[j] and pts.y[r] with pts.y[j] */

//                                                                            In the first pass, pts array (untampered, j=0, s= 15)
                
                /*
                   Finds the (pp-1)st smallest distance among the distances from (pts.x[j], pts.y[j])
                   of the points in pts with indices in [j+1, j+s-1] and stores this distance in pr.
                   Note that this is the same as finding the ppth smallest distance among the distances
                   from (pts.x[j], pts.y[j]) among the points in pts with indices [j+1, j+s-1], so it
                   chooses the desired pivot distance. */
                select_by_distance(pts, j, s, pp-1, &pr);                    // TODO (priyank) : check if its pp-1 or pp to be the pivotal element.           
                                                        /*
                                                                            In the first pass, pts array, j=0, s=15, pp = 7
                                                                        */


                /* This function needs to be implemented
                   Partitions the points at indices [j+1, j+s-1] in pts so that the points at indices
                   [j+1,j+pp-1] have distance at most pr from (pts.x[j], pts.y[j]), while the points
                   at indices [j+pp, j+s-1] have distance at least pr. */

//                partition_by_distance(pts, j, s, pr);   // this is implemented in select_by_distance().

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
            if (l < f && l + t >= f) 
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

    std::vector<point_list> points; // vector to temporarily store points

    FILE* file = fopen("../data/data.csv", "rt");
    for (;;) 
    {
       
        char buffer[100];
        point_list point;
        if ( !fgets(buffer, 100, file ) ) 
            {
            fclose( file );
            break;
            }
        
        point.coordinates =  buffer;
        size_t comma = point.coordinates.rfind(","); // performing regex to find comma
        sscanf(buffer + comma + 1, "%f", &point.y);
        comma = point.coordinates.rfind(",", comma-1);
        sscanf(buffer + comma + 1, "%f", &point.x);
       // printf("%lg, %lg\n", point.x, point.y);
        points.push_back(point);
                
    }

    /*Total number of points in the data*/
    //    point_list *pointer_to_data = points.data;

    int n= points.size();  

    std :: cout << "Number of data points to process" << n << std::endl ;  

    /*Declaring the structure of arrays for the data points to be stored.*/
    point_list pts[n];
    int i=0;
    for(std::vector<point_list>::iterator it=points.begin(); it!=points.end(); ++it)
        {   
            pts[i].x =  it->x;
            pts[i].y =  it->y;
            pts[i].n = n;
            i++;
        }

    /* Printing data */
    print(pts,n);

    /* Calling the build process of the tree */
    build_vp_tree(pts,n);

    return 0;
}