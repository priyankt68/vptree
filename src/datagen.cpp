#include <iostream>
#include <fstream>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
using namespace std;


int main()
{

 float a = 5.0;
 float b = 2.0;
 std::ofstream MyFile;
 MyFile.open ("../data/data.csv", ios::out | ios::ate | ios::app) ;
    for (int i=0;i<100000000;i++)
    {
    	float x =((float)rand()/(float)(RAND_MAX)) * a;
		float y =((float)rand()/(float)(RAND_MAX)) * b;
		MyFile<<x<<","<<y<<"\n";
    }
    
MyFile.close();

    return 0;
}

