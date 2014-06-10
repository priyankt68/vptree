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
#include <omp.h>

#include <new>
#include <vector>

#define left(i)   ((i << 1) +1)
#define right(i)  ((i << 1) +2)
#define parent(i) ((i - 1) >> 1)
#define floor_to_power_of_2(x) (1 << ((int) floor(log2(x))))    // t = 2 raised to power H, h=height of sub-tree at that level
#define power_of_2(x) pow(2,x)

/*Basic structure for data-points*/
struct point     
{
    float x;
    float y;
};

/* Type to represent a list of circles */
struct circle_list_t
{
    //int n;
    float *x, *y;
    float *r;
};


/*                              UTILITY FUNCTIONS                          */
/* ________________________________________________________________________*/

float compare (const void * a, const void * b)
{
  return ( *(float*)a - *(float*)b );
}


/* Function to print the data-points */

void print_input_data(point pts[], int n)
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
    int start = j+1;
    int end   = j+s-1;

    for ( int index = start ; index <= end; ++index)
    {   
        dist[index] =  sqrt(((pts[index].x - pts[j].x)*(pts[index].x - pts[j].x)) + ((pts[index].y - pts[j].y)*(pts[index].y - pts[j].y))) ;    // dist array contains the  
    }   
 }


int partition(float* input, point pts[], int p, int r)
{
    float pivot = input[r];
   
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
            swap(            &pts[p].x , &pts[r].x          );
            swap(            &pts[p].y , &pts[r].y          );

        }
    }
   
    return r;
}


float quick_select(float* input, point pts[], int p, int r, int k)
{
    if ( p == r ) return input[p];   // if the list just contains one element.
    
    int j = partition(input,pts, p, r);   // partition the elements with all elements

    int length = j - p + 1;     // changing the length 

    if ( length == k )  return input[j];          // we've reached the element.
    
    else if ( k < length ) return quick_select(input,pts, p, j - 1, k);

    else  return quick_select(input, pts,j + 1, r, k - length);
}


/*Partitions the points at indices [j+1, j+s-1] in pts so that the points at indices [j+1,j+pp-1] 
 have distance at most pr from (pts.x[j], pts.y[j]), while the points at indices [j+pp, j+s-1] 
 have distance at least pr. */
void partition_by_distance(point pts[], int j, int s, float pr, float distance[])
{

    int p = j + 1;
    int r = j + s - 1;
  /*  std::cout << " p = " << p << std :: endl;
    std::cout << " r = " << r << std :: endl;

    std::cout << " pr = " << pr << std :: endl;
    std::cout << " Distance" << std :: endl;

for (int index = p; index <= r; ++index)
    {
      std::cout << index << ": "<<distance[index] << std :: endl;
    }
*/
    int index_of_radii;
    for (int index = p; index <= r; ++index)
    {
        if( pr == distance[index])
        {
            index_of_radii = index;
        }
    }

  //  std::cout << "INDEX OF RADII "<<index_of_radii << std :: endl;
    

    // moving the pivot to the last position, before applying the partitioning algorithm

    //std::cout << "distance[index_of_radii]" << distance[index_of_radii] <<std :: endl;

    //std::cout << "distance[r]" << distance[r] <<std :: endl; 

    float temp = distance[index_of_radii];
    distance[index_of_radii] =  distance[r];
    distance[r] = temp;
    
    swap(&pts[index_of_radii].x , &pts[r].x);
    swap(&pts[index_of_radii].y , &pts[r].y);

    float pivot = distance[r];
      
    while ( p < r )
    {
        while ( distance[p] < pivot )
            p++;
        
        while ( distance[r] > pivot )
            r--;
        
        if ( distance[p] == distance[r] )
            p++;

        else if ( p < r ) 
        {                                              
            swap(  &distance[p] , &distance[r]);
            swap(            &pts[p].x , &pts[r].x          );
            swap(            &pts[p].y , &pts[r].y          );
        }
    }
    
}

/* Chooses a random index r in [j, j+s-1] and swaps pts.x[r] with pts`.x[j] 
and pts.y[r] with  pts.y[j] */
void random_point(point pts[],int j,int s)
{
   int r = (rand() % s ) + j;
       
   swap(&pts[r].x,&pts[j].x);
   swap(&pts[r].y,&pts[j].y);
}


/* 
Finds the (pp-1)st smallest distance among the distances from (pts.x[j], pts.y[j]) of the points 
in pts with indices in [j+1, j+s-1] and stores this distance in pr. Note that this is the same 
as finding the ppth smallest distance among the distances from (pts.x[j], pts.y[j]) among the 
points in pts with indices [j, j+s-1], so it chooses the desired pivot distance. 
*/
 void select_by_distance(point pts[], int j, int s, int pivot, float *pr, float distance[])  
{

    int start = j+1;
    int end   = j+s-1;

    int index_of_radii;   // index of the radii

    *pr = quick_select(distance,pts, j+1, j+s-1 , pivot);

  /*  std :: cout << "Printing points order" << std::endl;
            for (int index = start ; index <= end; ++index)
    {
       std :: cout << distance[index]  << "::" <<pts[index].x << "," << pts[index].y << std:: endl;
    }      
    */

}

void print_vp_tree(circle_list_t *vp, int n)
{
    std :: cout << "Printing vp-tree" << std :: endl;

    for (int i = 0; i < n; ++i)
    {
        std :: cout  << vp->x[i] <<" " << vp->y[i] << " "<< vp->r[i] << " " << std :: endl; 
    }

}

int read_input_data(point pts[])
{   
    std :: ifstream infile;
    infile.open("../data/data.txt"); // open file
    int i=0;

    if(infile)
        {  
        
        std :: string s="";
        
            while(infile)
            {
                getline(infile,s);
                char* pEnd;
                pts[i].x = (strtod (s.c_str(), &pEnd)) ;
                pts[i].y = (strtod (pEnd, NULL));
                i++;    
            }
        
        }

    infile.close();
    
      
    return i-1;
}

void save_vp_tree(circle_list_t *vp, int n)
{
    std::ofstream MyFile;
    MyFile.open ("../data/vp_tree_build.txt",std::ofstream::out | std::ofstream::trunc) ;
    
    for (int i=0;i<n;i++)
    {
        
        MyFile << vp->x[i] << " " << vp->y[i] << " " << vp->r[i] << "\n";
    }
    
MyFile.close();

}

/*Partitions the points at indices [j+1, j+s-1] in pts so that the points at indices [j+1,j+pp-1] 
 have distance at most pr from (pts.x[j], pts.y[j]), while the points at indices [j+pp, j+s-1] 
 have distance at least pr. */
 
circle_list_t* build_vp_tree (point pts[], int n, circle_list_t *vp)
{
    int sum=0;
    int fringe=0;  
    

    /* Global distance array*/

    float distance_vector[n];  // this distance array is accessed every time distance is calculated. No other distance vector is calculated again.


    /* n rounded down to the closest power of 2, also the number of leaves
       in the left subtree of the root if that subtree is in fact full */
    int t = floor_to_power_of_2(n);

    /* Number of leaves in the "fringe" (last level) of the tree. */
    int f = n - t + 1;
  
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

    for (i = 0; i < n ;) 
    {

        /* Size of full subtree at current level */
        s = (t << 1) - 1;                                             
     
        /* Pivot position within a full subtree */
        pp = t - 1;                                                   

        /* We're starting at the leftmost leaf */
        l = 0;                                                        
        di = 0;  
        
       /* Calculate starting position of points array of leftmost subtree */
        for (j = 0, ii = i; ii; ii = ii >> 1,j++);                          // ii gets divided by two every time the iteration proceeds.       
        
        /* Loop over subtrees at current level */
        while (j < n) 
        {
            /* The transition subtree that may not be full needs special treatment */
            if (l < f && (l + t) >= f)
             {

                /* Calculate subtree size */
                s = t - 1 + f - l;  

                /* Reduce fringe size for the remainder of this level (and for the next) */
                
                /* Calculate position of pivot */
                if ((f - l ) < t>>1)
                { 
                    pp = (t/2 - 1) + (f - l) ; //todo  (t+f-l-1 )

                }

            }
            /*
std :: cout << "---------------------------"  << std :: endl;
        std :: cout << "j = " << j << std :: endl;

                  std :: cout << "s = " << s << std :: endl;
                    std :: cout << "t = " << t << std :: endl;
                    std :: cout << "pp = " << pp << std :: endl;
                    std :: cout << "l = " << l << std :: endl;
                    std :: cout << "i = " << i << std :: endl;
                     std :: cout << "di = " << di << std :: endl;
*/
                    
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

               compute_distance(pts, j, s, distance_vector);

               select_by_distance(pts, j, s, pp, &pr, distance_vector);                    // TODO (priyank) : check if its pp-1 or pp to be the pivotal element.           
             
             /*       std :: cout << "-----" << std:: endl;
                    for (int index = j+1 ; index <= j+s-1; ++index)
    {
       std :: cout << distance_vector[index]  << "::" <<pts[index].x << "," << pts[index].y << std:: endl;
    }                                          
               */                                                              

                /* 
                   Partitions the points at indices [j+1, j+s-1] in pts so that the points at indices
                   [j+1,j+pp-1] have distance at most pr from (pts.x[j], pts.y[j]), while the points
                   at indices [j+pp, j+s-1] have distance at least pr. */

                partition_by_distance(pts, j, s, pr, distance_vector);  
             
                 
   //             std :: cout << "Printing points order" << std::endl;
          //      print(pts,18);
                /* Store the pivot radius in radius array of vp */
                vp->r[i] = pr;
             
            }

          //s  std :: cout << "Values to be put in the vantage point tree" << std::endl;
       
                
          
          
            vp->x[i] = pts[j].x;
            vp->y[i] = pts[j].y;
           
            //std :: cout << vp.x[i] << "," << vp.y[i] << std :: endl;
            ++i;
            ++di;


            //std :: cout << "Value of j" << j<<std :: endl;
           // std :: cout << "Value of s" << s<<std :: endl;
           /* Calculate index in pts array for next subtree */
            for (j += (s), ii = di; !(ii & 1) ; ii = ii >> 1,j++);


            // std :: cout << "Value of j later" << j<<std :: endl;
        ///std :: cout << "Value of i,  later" << i<<std :: endl;
          //  std :: cout << "Value of di later" << di <<std :: endl;
                
             if ((l < f) && (l + t >= f) )
             {    //( l < f && (l + 2t >= f ) )
           //         std :: cout << "Inside the transition tree" << std:: endl;

                  // std :: cout << "s = " << s << std :: endl;
                   // std :: cout << "t = " << t << std :: endl;
                    //std :: cout << "pp = " << pp << std :: endl;
                    //std :: cout << "l = " << l << std :: endl;
                s = t - 1 ;
                t >>= 1;
                pp = t - 1;
                l = f;    

/*
                   std :: cout << "s = " << s << std :: endl;
                    std :: cout << "t = " << t << std :: endl;
                    std :: cout << "pp = " << pp << std :: endl;
                    std :: cout << "l = " << l << std :: endl;
                    std :: cout << "di = " << di << std :: endl;
*/
            } else l += t;

           
        }
    }
  
}


//===== Main program ========================================================
int main()
{



    point *pts = new point[100];

    if (pts == NULL)
    {
        std :: cout << "points not allocated" << std :: endl;
    }

    int n = read_input_data(pts);

    print_input_data(pts,n);                    /* Printing input data */
//    circle_list_t *vp = new circle_list_t();  // this initialises the struct 
      //values to zero.
  
/*
    if (vp == NULL)
    {
        std :: cout << "points for vptree  not allocated" << std :: endl;   
    }
    vp->x = (float *)malloc(((n*2) + (n)) * sizeof(float));
    vp->y = vp->x + n;
    vp->r = vp->y + n;

    
    double startp = omp_get_wtime();

    build_vp_tree(pts,n,vp);       /* Calling the build process of the tree */

//    double endp = omp_get_wtime();

  //  delete pts;
   // std :: cout << "Build Process for " << n << " :: " << (endp - startp) << std :: endl;

    //print_vp_tree(vp,n);           /* Printing VP tree */

    //save_vp_tree(vp,n);            /* Saving VP tree to file */

    //free(vp);

   return 0;
}
