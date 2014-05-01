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
 MyFile.open ("../data/data.csv") ;
 	//MyFile << "X" << "," << "Y" << endl;
    for (int i=0,j=0;i<16,j<15;i++,j++)
    {
    	float x =((float)rand()/(float)(RAND_MAX)) * a;
		float y =((float)rand()/(float)(RAND_MAX)) * b;
		MyFile<<x<<","<<y<<"\n";
    }
    
MyFile.close();

    return 0;
}

