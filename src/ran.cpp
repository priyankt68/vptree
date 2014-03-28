#include<iostream>
#include<vector>
#include<algorithm>


int main()
{
 std::vector<int> v(10,0);

std::cout<<"size of v"<<v.size();

  for (int i = 0; i < 10; ++i)
  	 {
 	 	v.push_back(i);
 			
 	}


  for (int i = 0; i < 10; ++i)
  	 {
 	 
 	 try
 	 {
 	 	std::cout<<v.at(i)<<" "<<std::endl;
 	 }
 	 catch (std::exception& e)
 	 {
		std::cout<<"array index out of bounds"<<i<<std::endl;
 	 }
 			
 	}
return 0;
}