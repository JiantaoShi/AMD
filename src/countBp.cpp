#include "stdafx.h"
#include "countBp.h"
long countBp(const char *filename)
{
	ifstream input;input.open(filename,ios::in);
	if(input.fail()) {cout<<"Error! Can not open "<<filename<<"! Please check it!"<<endl;exit(0);}
	string head,buff;
	long count1=0;
	int signal=0;
	while(getline(input,head))
	{
		if(head[0]=='>')
		{
			if(signal==1)  
			{ 
				size_t len=buff.size();
				for(size_t j=0;j<len;j++)
				{
					if(buff[j]=='N')continue;
					  count1++;
				}
				buff.clear();
			}
			else { buff.clear();signal=1;}
		}
		else buff+=head;
	}
	size_t len=buff.size();
	for(size_t j=0;j<len;j++)
	{
		if(buff[j]=='N')continue;
		  count1++;
	}
	buff.clear();
    input.close();
   return count1;
}
