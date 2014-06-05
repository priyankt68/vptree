#include <iostream>
#include <fstream>
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
		if (d<r && (d!=0 || r!=0))
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
		if (d > r && (d!=0 || r!=0))
		{
			std :: cout << "WRONG POSITION" << std:: endl;
			std:: cout << "Distance : " << d << "|| " << x <<" " << y <<  std :: endl;			
		}	
		

	}
/*
	std:: cout << "Distance : " << d << std:: endl;
	 std:: cout << x <<" " << y << " "<< "with root being" << p << " " << q << std :: endl;
*/	
}

void post_order_vp_tree(circle_list_t vp[], int nodeindex, int &i, int n, int root, int flag)
{
	if (nodeindex > n)
      	return ;	
    

    post_order_vp_tree(vp, left(nodeindex), i, n, root,flag) ;
    
    if (vp[nodeindex].x == vp[root].x && vp[nodeindex].y == vp[root].y)
    {
    	std:: cout << "\n\t|| RIGHT SUBTREE || "  << std:: endl;
    	flag = 1;  // when the flag is one, we say that, yes we've entered the right side of my sub-tree.
    }
    
    position_check(vp[nodeindex].x , vp[nodeindex].y, vp[root].x, vp[root].y, vp[root].r,flag);
  
    post_order_vp_tree(vp,  right(nodeindex), ++i, n, root,flag);
    
    
}


void vp_tree_check(circle_list_t vp[], int n)
{
		
}


int main(int argc, char const *argv[])
{
	circle_list_t *vp_tree = new circle_list_t[50] ;
	int n = read_data(vp_tree);  
	//vp_tree_check(vp_tree,n);
	int j=0;
	for (int i = 0; i < n-1; ++i)
	{	
		std :: cout << "-----------------------\n" << std:: endl;
		std:: cout << "\t\tRoot: "  << vp_tree[i].x <<" " << vp_tree[i].y << " "<< vp_tree[i].r << " " << std :: endl;
		post_order_vp_tree(vp_tree,i,j,n-1,i,0);
	}
	
	return 0;
}