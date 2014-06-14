#include <stdlib.h>
#include <algorithm>
#include <iostream>
#include <stdlib.h>
#include <sstream>
#include <fstream>
#include <vector>
#include <omp.h>
#include <queue>
#include <math.h>

#define left(i)   ((2*i) + 1)
#define right(i)  ((2*i) + 2)
#define parent(i) ((i - 1) >> 1)

/* Type to represent a list of circles */
typedef struct {
    float x, y, r;
} circle_list_t;

/* Type to represent a list of query data-points */
typedef struct {
    float x, y;
} query_list;

    // An item on the intermediate result queue
    struct HeapItem {
        HeapItem( float index, float dist, float x, float y) :
        index(index), dist(dist) ,x(x), y(y){}
        float index;
        float dist;
        float x;
        float y;
        bool operator<(const HeapItem& o) const {
            return dist > o.dist;
        }
    };
    
float compute_distance(float *vp_x,float *vp_y,float *qx,float *qy)
{

	return  sqrt( (*vp_x- *qx)*(*vp_x - *qx) + (*vp_y- *qy)*(*vp_y - *qy) ) ;
	
}

 void print_vp_tree_data(std::vector<circle_list_t> vp_tree,int n)
 {
	
	for (int i=0;i<n; ++i)
	{
		std :: cout << vp_tree[i].x << " " <<vp_tree[i].y  << " " << " " << vp_tree[i].r << std::endl;
	}

 }

void print_query_data(std::vector<query_list>& query,int n)
 {
	
	for (int i=0;i<n; ++i)
	{
		std :: cout << query[i].x << " " <<query[i].y << std::endl;
	}

 }


int read_query_data(std::vector<query_list>& query)
{	

    std :: ifstream query_infile;
    query_infile.open("../data/query_data.txt"); // open file
    int i=0;
    if(query_infile)
    {  
        std :: string s="";

        while(query_infile)
        {
        getline(query_infile,s);
        char* pEnd;
        query[i].x = (strtod (s.c_str(), &pEnd)) ;
        query[i].y =(strtod (pEnd, NULL));

        i++;
        }
                
    }

    return i-1;
}
int read_query_size()
{
    std :: ifstream infile;
    infile.open("../data/query.txt"); // open file
    int n;
    if(infile)
    {   
         std :: string s="";
        getline(infile,s);
        
        n = atoi(s.c_str());
        
    }

    return n;
}

int read_vp_tree_size()
{
    std :: ifstream infile;
    infile.open("../data/data.txt"); // open file
    int n;
    if(infile)
    {   
        std :: string s="";
        getline(infile,s);
        
        n = atoi(s.c_str());
        
    }

    return n;
}
int read_vp_tree_data(std::vector<circle_list_t> vp_tree)
{
	std :: ifstream infile;
    infile.open("../data/vp_tree_build.txt"); // open file
   
	int i=0;
 	int c = 0;
    int l = 0;
	while(!infile.eof() )
	{
		if (c%3 == 0) 	infile >> vp_tree[l].x;
	
		else if (c%3 == 1) infile >> vp_tree[l].y;
		
		else
			{
			infile >> vp_tree[l].r;
			l++;
			}
	c++;
    }

  	infile.close();
	return l;
}

void nn_search(std::vector<circle_list_t>& vp, float qx, float qy,  int n, int k,std::priority_queue<HeapItem>& heap)
{
	
    /* Index of current node */
    int i = 0;

    /* Index of previous node */
    int j = -1;

    /* Distance of currently checked point */
    float d;
 	    	    	    
    /* Best distance found so far, corresponding point stored in (nnx, nny) */
	// /float t = compute_distance(&(vp[0].x), &(vp[0].y), &qx, &qy);
	float t = 1000.99;
	//std :: cout << t << std::endl;
    // heap.push(HeapItem(t, 0, vp[0].x, vp[0].y)); 
	//std :: cout <<"NN_x = " << *nnx << " NN_y = " << *nny << std:: endl; 
    while (i >= 0 )  
    {
    	if ( i >= n) j=i , i = parent(i);
    	else
    	{
        /* Calculate distance of query from node's point */
		d = compute_distance(&vp[i].x, &vp[i].y, &qx, &qy);    // Compute distance between target and current node
        //std :: cout << " i =  "<< i << " d = " << d << " " << " t = " << t << 
        //" r = "<< vp[i].r <<  std :: endl;

        /* Do this if we came from the parent */
        if (j < i) 
        {

            /* First check the node's point's distance from (qx, qy) */
            if (d < t) {
                
                if(heap.size() == k) heap.pop();                 // remove furthest node from result list (if we already have k results)
           		heap.push(HeapItem(vp[i].r, d, vp[i].x, vp[i].y));           // add current node to result list
            	if(heap.size() == k) t = heap.top().dist;  // t tracks the farthest distance
                                
                          //     std :: cout << "Change in the value of nearest neighbour ";
            }

            /* Proceed to left or right child depending on whether q is inside
               or outside the circle */
            if (d <= vp[i].r) j = i, i = left(i); 
            else              j = i, i = right(i);

        }
      
        /* Do this if we came from a child */
        else {

		//	std :: cout << " \nComing from a child" << std :: endl;

            if (j == left(i)  && (d <= (vp[i].r)) && (t >= (vp[i].r - d)) )
            	{
            		j = i;
            	    i = right(i);
            	}

            	
            else if (j == right(i) && d >  vp[i].r && t >= d - vp[i].r)
	            { 
	            	j = i;
	            	i = left(i);
	            }
            else 
	            {
	             j = i;
	             i = parent(i);
	            }
       
	    	}
	}	}
}



int main(int argc, char const *argv[])
{

    int  k;
    if (argc != 2)
    {
        std :: cout << "Usage " << argv[0] << " " << "k" << std :: endl;
        exit(0);
    }
    else
         k = atoi(argv[1]);

    /*Reading Query data*/
	std :: ifstream infile;
    infile.open("../data/query.txt"); // open file
    int query_size;
    if(infile)
    {   
        std :: string s="";
        getline(infile,s);
        
        
        query_size = atoi(s.c_str());
     //   std :: cout << query_size << std :: endl;

    }

    std::vector<query_list> query(query_size); 
        int i=0;
        while(i<query_size)
        { 
        std :: string s="";
         
        getline(infile,s);

        char* pEnd;

        query[i].x = (strtod (s.c_str(), &pEnd)) ;
        query[i].y =(strtod (pEnd, NULL));
        i++;
        }
    

    infile.close();


  //  print_query_data(query,query_size);



    /*Reading VP Tree data*/

	std :: ifstream vp_infile;
    vp_infile.open("../data/vp_tree_build.txt"); // open file
    int vp_data_size;
    if(vp_infile)
    {   
        std :: string s="";
        getline(vp_infile,s);
        
        vp_data_size = atoi(s.c_str());
        //std :: cout << vp_data_size << std :: endl;

    }
  //  std :: cout << vp_data_size << std :: endl;
    std::vector<circle_list_t> vp(vp_data_size); // vector to temporarily store points   

    
    
   
    i=0;
    int c = 0;
    int l = 0;
    while(!vp_infile.eof() )
    {
        if (c%3 == 0)   vp_infile >> vp[l].x;
    
        else if (c%3 == 1) vp_infile >> vp[l].y;
        
        else
            {
            vp_infile >> vp[l].r;
            l++;
            }
    c++;
    }

    
    

    vp_infile.close();
    

    //print_vp_tree_data(vp,vp_data_size);


    int iteration = query_size;
    	
	double startp = omp_get_wtime();
	while (query_size--)
	{
		 // Use a priority queue to store intermediate results on
    std::priority_queue<HeapItem> heap;

    

    nn_search(vp,query[i].x,query[i].y,vp_data_size,k,heap);

    




//		std :: cout << query[i].x << "," << query[i].y << std :: endl;
		/*
		while (!heap.empty())
  		{
  			std :: cout << heap.top().dist  << " :: " << heap.top().x << " , " << heap.top().y << std :: endl;
    	 	heap.pop();
    	 	std :: cout << "\n";
  		}
*/
		i++;
//		std::cout << '\n';
	}
	double endp = omp_get_wtime();

    std :: cout << iteration << " " << endp - startp << std :: endl;

	return 0;
}





