#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <stdint.h>
#include <string>
#include <string.h>
#include <math.h>
#include <vector>
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
int main()
{
std::vector<Point> points;
printf("Reading coordinates database...\n");
FILE* file = fopen("data.txt", "rt");
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
        printf("%lg, %lg\n", point.x, point.y);
        points.push_back(point);
        //if(points.size()>50000)break;
    }

printf("My points are \n");
    for(std::vector<Point>::iterator it=points.begin(); it!=points.end(); ++it)  
        printf("%2.1f, %2.1f \n", (it->x),(it->y));
  
}