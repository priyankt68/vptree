#include<iostream>
#include<vector>
using namespace std;
struct test{

	int f;
	char ff;
};

int main()
{
	std::vector<int> vector1(10,0);
	std::vector<test> vector2(10);
	std::vector<char> v(10,'a');

	//display the original size
	cout<<"size="<<v.size()<<endl;


	for (int i = 0; i < v.size(); ++i)
	{
		//v[i]
		cout<<v[i]+2<<endl;
		//cout<<vector1[i]+1<<endl;
	}

	vector1[1]=1;
	
    

	//std::cout<<vector1[1];
	//std::cout<<vector2[0];
	
	

return 0;

}