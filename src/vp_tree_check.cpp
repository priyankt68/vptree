#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <math.h>

#define left(i)   ((i << 1) +1)
#define right(i)  ((i << 1) +2)
#define parent(i) ((i - 1) >> 1)

typedef struct {
    float x, y, r;
} circle_list_t;

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

void print_vp_tree(circle_list_t vp[], int n)
{
    std :: cout << "Printing vp-tree" << std :: endl;

    for (int i = 0; i < n; ++i)
    {   
    	std :: cout  << vp[i].x <<" " << vp[i].y << " "<< vp[i].r << " " << std :: endl; 
    }
}

void position_check(float x,float y, float p, float q, float r, int flag)
{	

	/*	If flag is equal to one, we check for the less than or euqal to condition, 
		else we're in the right sub-tree and we need to check for the greater than equal to condition.
	*/
    float d = sqrt((x-p)*(x-p) + (y-q)*(y-q));
	
	if (flag) 
	{

		if (d==0 || r==0)
		{
			
		}
		
		// we're in the right sub-tree, check for greater than equal to condition
		if (d<=r && (x!=p && y!=q) )
		{
			std :: cout << "WRONG POSITION" << std:: endl;
			std:: cout << "Distance : " << d << "|| " << x <<" " << y << std :: endl;			
		}
		
	}
	else
	{
		if (d==0 || r==0)
		{
			
		}
		if(d==r)
		{

		}
	
		if (d > r )
		{
			std :: cout << "WRONG POSITION neechr" << std:: endl;
			std:: cout << "Distance : " << d << "|| " << x <<" " << y <<  std :: endl;			
		}	
		

	}
/*
	std:: cout << "Distance : " << d << std:: endl;
	 std:: cout << x <<" " << y << " "<< "with root being" << p << " " << q << std :: endl;
*/	
}

void post_order_vp_tree(std::vector<circle_list_t>& vp, int nodeindex, int &i, int n, int root, int flag)
{
	if (nodeindex > n)
      	return ;	
    

    post_order_vp_tree(vp, left(nodeindex), i, n, root,flag) ;
    
    if (vp[nodeindex].x == vp[root].x && vp[nodeindex].y == vp[root].y)
    {
    	//std:: cout << "\n\t|| RIGHT SUBTREE || "  << std:: endl;
    	flag = 1;  // when the flag is one, we say that, yes, we've entered the right side of my sub-tree.
    }
    
    position_check(vp[nodeindex].x , vp[nodeindex].y, vp[root].x, vp[root].y, vp[root].r,flag);
  
    post_order_vp_tree(vp,  right(nodeindex), ++i, n, root,flag);
    
    
}


int main(int argc, char const *argv[])
{
	/*Reading VP Tree data*/

	std :: ifstream vp_infile;
    vp_infile.open("../data/vp_tree_build_500.txt"); // open file

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
   
   
   int i=0;
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
  
	int n = vp_data_size;

	std :: cout <<  n << std :: endl;
	
	int j=0;
	for (int i = 0; i < n-1; ++i)
	{	
		std :: cout << "-----------------------\n" << std:: endl;
		std:: cout << "\t\tRoot: "  << vp[i].x <<" " << vp[i].y << " "<< vp[i].r << " " << std :: endl;
		post_order_vp_tree(vp,i,j,n-1,i,0);
	}
	
	return 0;
}