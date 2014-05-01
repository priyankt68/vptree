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
//void printdata(std::vector<Point>& temp);

struct HeapItem 
    {
        HeapItem( int index, double dist) :
        index(index), dist(dist) {}
        int index;
        double dist;
        bool operator<(const HeapItem& o) const   // operator overloading

        {
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

void printdata(std::vector<Point>& temp)
{
    for(std::vector<Point>::iterator it=temp.begin(); it!=temp.end(); ++it)  
        printf("%2.1f, %2.1f \n", (it->x),(it->y));
}

void build_vp_tree(std::vector<Point>& tem, int s)
{


    /*Segregating data for the sake of memory colaescing in gpus and openmp thread memory accesses.*/
    std::vector<float> pointx;
    std::vector<float> pointy;
    Point point;
        for(std::vector<Point>::iterator it=tem.begin(); it!=tem.end(); ++it)
        {
            pointx.push_back(it->x);
            pointy.push_back(it->y);
        }

        /*Segregating data for the sake of colaescing memeory accesses in GPUs*/
        for(std::vector<float>::iterator itx=pointx.begin(),ity=pointy.begin(); itx!=pointx.end(),ity!=pointy.end(); ++itx,++ity)
            std::cout << *itx <<","<<*ity<< std::endl;
        


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
       // printf("%lg, %lg\n", point.x, point.y);
        points.push_back(point);
        //if(points.size()>50000)break;
        
    }

        //printdata(points);

    std::cout << "Array Size: " << points.size() << std::endl;
    int pointssize = points.size();

      // segregatedata(points);

        /*Function to build vantage point tree.   Returns a permuted structure representing vantage point tree built in-place procedure.*/
    build_vp_tree(points, pointssize);




/*
printf("My points are \n");
    for(std::vector<Point>::iterator it=points.begin(); it!=points.end(); ++it)  
        printf("%2.1f, %2.1f \n", (it->x),(it->y));

*/

/* Starting of process of construction

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


    end = omp_get_wtime();
  // QueryPerformanceCounter( &end );    // end performance to create a tree

 /*Creation of the tree is over at this point*/

    printf("create took %f \n",(end- start));
//    file1<< NoP << "\t" << (end-start);
int c=0;
for (std::vector<Point>::iterator it=points.begin(); it!=points.end(); ++it)
{
   // std::cout << ' ' << it->x << ',' << it->y << std::endl;
  //  c++;     // counting the number of points in total.
}

    //std::cout << "Number of points:" << c << std::endl;

/*Starting of search process 

                            in vantage point trees      */

    //Randomly selecting 10% of the dataset to query

    file = fopen("../data/data.csv", "rt");

    //Randomly shuffling to
    std::random_shuffle ( points.begin(), points.end() );    // #include<algorithm>
std::vector<Point>::iterator iter = points.begin();
    std::vector<Point> queryset;
    for (int i=0;i< ((points.size())/10);i++)
    {
        
        Point point;    //query point
        point.x = iter->x ;  // query point being: ( 3.1 , 1.1 )
        point.y = iter->y;
        queryset.push_back(point);
        ++iter;      
    }
    /*for (std::vector<Point>::iterator iter = points.begin(); iter!= points.end(); ++iter)
    {
        Point point;    //query point
        point.x = iter->x ;  // query point being: ( 3.1 , 1.1 )
        point.y = iter->y;
        queryset.push_back(point);
    }*/

   // for (std::vector<Point>::iterator iter = queryset.begin(); iter!= queryset.end(); ++iter)
    //{
       // std::cout << ' ' << iter->x << ',' << iter->y << std::endl;
      //  c++;
    //}

   // printf("%d\n Queryset %ld",c, queryset.size() );
   /* Point point;    //query point
    point.x = 3.1 ;  // query point being: ( 3.1 , 1.1 )
    point.y = 1.1;
*/


    std::vector<Point> results;
    std::vector<double> distances;

    //std::ofstream file1;

    file1.open ("../data/serial_vs_parallel.csv") ; //ios::out | ios::ate | ios::app

    double startp = omp_get_wtime();
#pragma omp parallel for
    for (int i = 0; i < queryset.size(); ++i)
    {
        tree.search(queryset[i],8,&results,&distances);
    }

    double endp = omp_get_wtime();

    printf("Parallel Search took %f\n", (endp-startp));


double starts = omp_get_wtime();
//#pragma omp parallel for 
    for (int i = 0; i < queryset.size(); ++i)
    {
        tree.search(queryset[i],8,&results,&distances);
    }

double ends = omp_get_wtime();

    printf("serial Search took %f\n", (ends - starts));

file1<< queryset.size() << "," << (endp-startp) << "," << (ends-starts) << std::endl ;


    for( int i = 0; i < results.size(); i++ ) 
        printf("%s %lg\n", results[i].coordinates.c_str(), distances[i]);

    printf("---\n");
    //QueryPerformanceCounter( &start );

  //  start = omp_get_wtime();

      
    //linear_search( points, point, 8, &results, &distances );

    //end = omp_get_wtime();
   // QueryPerformanceCounter( &end );
   // printf("Linear search took %f\n", (end-start));

    //for( int i = 0; i < results.size(); i++ )
      //  printf("%s %lg\n", results[i].coordinates.c_str(), distances[i]);
    


    return 0;


}