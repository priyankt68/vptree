#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>

#define left(i)   ((2*i) +1)
#define right(i)  ((2*i) +2)
#define parent(i) ((i - 1) >> 1)

/* Type to represent a list of circles */
typedef struct {
    float x, y, r;
} circle_list_t;


float compute_distance(float *vp_x,float *vp_y,float *qx,float *qy)
{

 float dista =  ( (*vp_x- *qx)*(*vp_x - *qx) + (*vp_y- *qy)*(*vp_y - *qy) ) ;
  return dista;
}

 void print_data(circle_list_t vp_tree[],int n)
 {
	
	for (int i=0;i<n; ++i)
	{
		std :: cout << vp_tree[i].x << " " <<vp_tree[i].y  << " " << " " << vp_tree[i].r << std::endl;
	}

 }

int read_data(circle_list_t vp_tree[])
{
	std :: ifstream infile;
    infile.open("../data/vp_tree_build.txt"); // open file
    std::string line;
    
    int i=0;

    
	int c = 0;
    int l = 0;

	while(!infile.eof() )
	{
	
		if (c%3 == 0)
			{
				infile >> vp_tree[l].x;
			}
		else if (c%3 == 1)
			{
				infile >> vp_tree[l].y;
			}
		else
			{
				infile >> vp_tree[l].r;
				l++;
			}
	
	c++;
    }
  
//  std:: cout << c << std::endl;      

	infile.close();
	
	return l;
}


 /* Updated NN search procedure for above VP tree */

void nn_search(circle_list_t vp[], float qx, float qy, float *nnx, float *nny, int n)
{
    /* Index of current node */
    int i = 0;

    /* Index of previous node */
    int j = -1;

    /* Distance of currently checked point */
    float d;
 	    	    	    
    /* Best distance found so far, corresponding point stored in (nnx, nny) */
		float t = compute_distance(&(vp[0].x), &(vp[0].y), &qx, &qy);
		std :: cout << t << std::endl;
  //  float t;
//		return;
		std :: cout << vp[0].x << std::endl;
		std :: cout << vp[0].y << std::endl;
    *nnx = vp[0].x;
    
   *nny = vp[0].y;

//std :: cout << "nn_x = "<< nnx << "nn_y" << nny << std::endl;
 
    while (i >= 0 && i <= n) 
    {


        /* Calculate distance of query from node's point */
        d = compute_distance(&vp[i].x, &vp[i].y, &qx, &qy);    // Compute distance between target and current node

       std :: cout << d << " " << i << std :: endl;

        /* Do this if we came from the parent */
        if (j < i) 
        {

            /* First check the node's point's distance from (qx, qy) */
            if (d < t) {
                t = d;
                *nnx = vp[i].x;
                *nny = vp[i].y;
            }

 
            /* Proceed to left or right child depending on whether q is inside
               or outside the circle */
            if (d <= vp[i].r) j = i, i = left(i); 
            else              j = i, i = right(i);


        }
      
        /* Do this if we came from a child */
        else {

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

	}	
}
int main(int argc, char const *argv[])
{

	if (argc != 3)
	{
		std :: cout << "Usage" << argv[0] << "<space>" << " query_x_cordinate <space>" << "query_y_cordinate" << std :: endl;
	}
	else
	{
		float qx = atof(argv[1]);
		float qy =  atof(argv[2]);
	}
	circle_list_t *vp_tree = new circle_list_t[50] ;
	
	int n = read_data(vp_tree);  // returns

	//print_data(vp_tree,n);

	float nnx,nny;
	nn_search(vp_tree,2.5,1.6,&nnx,&nny,n);

	std :: cout << "nn_x = "<< nnx << " nn_y" << nny << std::endl;
	
	
	return 0;
}





