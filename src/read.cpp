#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <stdint.h>
#include <string>
#include <string.h>
#include <math.h>
#include <vector>
#include <omp.h>
#include "vptree.h"

#define NoP 100
//ing namespace std;

/*
void QueryPerformanceCounter( uint64_t* val )
{
    timeval tv;
    struct timezone tz = {0, 0};
    gettimeofday( &tv, &tz );
    *val = tv.tv_sec * 1000000 + tv.tv_usec;
}
*/

struct HeapItem {
        HeapItem( int index, double dist) :
        index(index), dist(dist) {}
        int index;
        double dist;
        bool operator<(const HeapItem& o) const {
            return dist < o.dist;
        }
    };

struct Point
{
    /* data */
    std::string coordinates;
    double x;
    double y;
};

double distance( const Point& p1, const Point& p2 )
{
    double a = (p1.x-p2.x);
    double b = (p1.y-p2.y);
    return sqrt(a*a+b*b);
}

void linear_search( const std::vector<Point>& items, const Point& target, int k, std::vector<Point>* results, std::vector<double>* distances) 
{
    std::priority_queue<HeapItem> heap;
    for ( int i = 0; i < items.size(); i++ ) 
    {
        double dist = distance( target, items[i] );
        if ( heap.size() < k || dist < heap.top().dist ) 
        {
            if (heap.size() == k) heap.pop();
            heap.push( HeapItem( i, dist ) );
        }
    }

    results->clear();
    distances->clear();

    while( !heap.empty() ) 
    {
        results->push_back( items[heap.top().index] );
        distances->push_back( heap.top().dist );
        heap.pop();
    }

    std::reverse( results->begin(), results->end() );
    std::reverse( distances->begin(), distances->end() );
}

int main()
{
std::vector<Point> points;
printf("Reading coordinates database...\n");
FILE* file = fopen("../data/data.csv", "rt");
    for(;;) 
    {
        char buffer[100];
        Point point;
        if ( !fgets(buffer, 100, file ) ) 
            {
            fclose( file );
            break;
            }
        point.coordinates = buffer;
        size_t comma = point.coordinates.rfind(",");
        sscanf(buffer + comma + 1, "%lg", &point.y);
        comma = point.coordinates.rfind(",", comma-1);
        sscanf(buffer + comma + 1, "%lg", &point.x);
    //xcz    printf("%lg, %lg\n", point.x, point.y);
        points.push_back(point);
        //if(points.size()>50000)break;

    }


/*
printf("My points are \n");
    for(std::vector<Point>::iterator it=points.begin(); it!=points.end(); ++it)  
        printf("%2.1f, %2.1f \n", (it->x),(it->y));
  */

/*Starting of process of construction

                            in vantage point trees */


    VpTree<Point, distance> tree;    //declaring a tree
   //uint64_t start, end;
    double start, end;
 // QueryPerformanceCounter( &start );   //starting to measure the time in which the tree is created

    std::ofstream file1;
    
    //file1.open ("../data/timings.csv") ; //ios::out | ios::ate | ios::app
    //file1<< "NoP" << "\t" <<"create" << "\t" <<"serial_nn_search" << "\t" <<"par_nn_search" << "\t" << "linear_search" << std::endl ;

    start = omp_get_wtime();

    tree.create( points )  ; // so in vector -> Point , every point contains: double x , double y

    // printf("%d\n", tree._root.index);

    end = omp_get_wtime();
  // QueryPerformanceCounter( &end );    // end performance to create a tree

 /*Creation of the tree is over at this point*/

    printf("create took %f \n",(end- start));
    file1<< NoP << "\t" << (end-start);


/* Visualising the tree */

   // tree.visualise();    



/*Starting of search process 

                            in vantage point trees      */

    //Randomly selecting 10% of the dataset to query


    Point point;    //query point
    point.x = 3.1 ;  // query point being: ( 3.1 , 1.1 )
    point.y = 1.1;

    std::vector<Point> results;
    std::vector<double> distances;

    start = omp_get_wtime();

    tree.search( point, 8, &results, &distances );

    end = omp_get_wtime();

    printf("Search took %f\n", (end-start));

    for( int i = 0; i < results.size(); i++ ) 
        printf("%s %lg\n", results[i].coordinates.c_str(), distances[i]);

    printf("---\n");
    //QueryPerformanceCounter( &start );

    start = omp_get_wtime();

      
    linear_search( points, point, 8, &results, &distances );

    end = omp_get_wtime();
   // QueryPerformanceCounter( &end );
    printf("Linear search took %f\n", (end-start));

    for( int i = 0; i < results.size(); i++ )
        printf("%s %lg\n", results[i].coordinates.c_str(), distances[i]);
    


    return 0;


}