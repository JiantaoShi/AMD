#include "stdafx.h"
#include "Correlation.h"
void Translate(char *motif,char *Trans)
{ 
   size_t lens=strlen(motif);
   for(size_t i=0;i<lens;i++)
   {
	   switch(motif[lens-1-i])
	   {
	   case 'A':
		   Trans[i]='T';break;
	   case 'G':
		   Trans[i]='C';break;
	   case 'C':
		   Trans[i]='G';break;
	   case 'T':
		   Trans[i]='A';break;
	   case 'M':
		   Trans[i]='K';break;
	   case 'K':
	       Trans[i]='M';break;
	   case 'R':
		   Trans[i]='Y';break;
	   case 'Y':
		   Trans[i]='R';break;
	   case 'W':
		   Trans[i]='W';break;
	   case 'S':
		   Trans[i]='S';break;
	   case 'B':
		   Trans[i]='V';break;
	   case 'V':
		   Trans[i]='B';break;
	   case 'D':
		   Trans[i]='H';break;
	   case 'H':
		   Trans[i]='D';break;
	   case 'N':
		   Trans[i]='N';break;
	   case 'X':
		   Trans[i]='X';break;
	   default:
		   break;
	   }
   }
  Trans[lens]=NULL;
}
void countMatrix(double matrix[4][50],char *motif,double beta,int sign)
{
   int i=0;
   size_t lens=strlen(motif);
   while(motif[i]!=NULL)
   {
	   switch(motif[i])
	   {
	   case 'A':
		   matrix[0][i]++;break;
	   case 'G':
		   matrix[1][i]++;break;
	   case 'C':
		   matrix[2][i]++;break;
	   case 'T':
		   matrix[3][i]++;break;
	   case 'M':
		   matrix[0][i]++;matrix[2][i]++;break;
	   case 'K':
	       matrix[1][i]++;matrix[3][i]++;break;
	   case 'R':
		   matrix[0][i]++;matrix[1][i]++;break;
	   case 'Y':
		    matrix[2][i]++;matrix[3][i]++;break;
	   case 'W':
		   matrix[0][i]++;matrix[3][i]++;break;
	   case 'S':
		   matrix[1][i]++;matrix[2][i]++;break;
	   case 'B':
		   matrix[1][i]++;matrix[2][i]++;matrix[3][i]++;break;
	   case 'V':
		   matrix[0][i]++;matrix[1][i]++;matrix[2][i]++;break;
	   case 'D':
		   matrix[0][i]++;matrix[1][i]++;matrix[3][i]++;break;
	   case 'H':
		   matrix[0][i]++;matrix[2][i]++;matrix[3][i]++;break;
	   case 'N':
		   matrix[0][i]++;matrix[1][i]++;matrix[2][i]++;matrix[3][i]++;break;

	   default:
		   break;
	   }
	   i++;
   }
 if(sign==1)
 {
   for(size_t i=0;i<lens;i++)
   {
	   double buff=matrix[0][i]+matrix[1][i]+matrix[2][i]+matrix[3][i];
	   for(int j=0;j<4;j++) matrix[j][i]=(matrix[j][i]/buff+beta*bg[j])/(1.0+beta);	   
   }
 }
}
void countMatrixE(double matrix[4][50],char *motif,double &maxscore)
{
   int i=0;
   size_t lens=strlen(motif);
   while(motif[i]!=NULL)
   {
	   switch(motif[i])
	   {
	   case 'A':
		   matrix[0][i]+=1.1;maxscore+=1.1;break;
	   case 'G':
		   matrix[1][i]+=1.1;maxscore+=1.1;break;
	   case 'C':
		   matrix[2][i]+=1.1;maxscore+=1.1;break;
	   case 'T':
		   matrix[3][i]+=1.1;maxscore+=1.1;break;
	   case 'M':
		   matrix[0][i]++;matrix[2][i]++;maxscore+=1;break;
	   case 'K':
	       matrix[1][i]++;matrix[3][i]++;maxscore+=1;break;
	   case 'R':
		   matrix[0][i]++;matrix[1][i]++;maxscore+=1;break;
	   case 'Y':
		    matrix[2][i]++;matrix[3][i]++;maxscore+=1;break;
	   case 'W':
		   matrix[0][i]++;matrix[3][i]++;maxscore+=1;break;
	   case 'S':
		   matrix[1][i]++;matrix[2][i]++;maxscore+=1;break;
	   case 'B':
		   matrix[1][i]++;matrix[2][i]++;matrix[3][i]++;maxscore+=1;break;
	   case 'V':
		   matrix[0][i]++;matrix[1][i]++;matrix[2][i]++;maxscore+=1;break;
	   case 'D':
		   matrix[0][i]++;matrix[1][i]++;matrix[3][i]++;maxscore+=1;break;
	   case 'H':
		   matrix[0][i]++;matrix[2][i]++;matrix[3][i]++;maxscore+=1;break;
	   case 'N':
		   matrix[0][i]++;matrix[1][i]++;matrix[2][i]++;matrix[3][i]++;maxscore+=1;break;

	   default:
		   break;
	   }
	   i++;
   }
}
void Trans(double matrix1[4][50],char *motif,int start,int end,double beta)
{
   double matrix[4][50];
   memset(matrix,0,sizeof(matrix));
   size_t lens=strlen(motif);
   for(int j=0;j<start;j++)
   {
	   for(int i=0;i<4;i++)
		   matrix[i][j]=bg[i];
   }
   for(int j=lens+start;j<lens+start+end;j++)
   {
	   for(int i=0;i<4;i++)
		   matrix[i][j]=bg[i];
   }
   int i=0;
   while(motif[i]!=NULL)
   {
	   switch(motif[i])
	   {
	   case 'A':
		   matrix[0][i+start]++;break;
	   case 'G':
		   matrix[1][i+start]++;break;
	   case 'C':
		   matrix[2][i+start]++;break;
	   case 'T':
		   matrix[3][i+start]++;break;
	   case 'M':
		   matrix[0][i+start]++;matrix[2][i+start]++;break;
	   case 'K':
	       matrix[1][i+start]++;matrix[3][i+start]++;break;
	   case 'R':
		   matrix[0][i+start]++;matrix[1][i+start]++;break;
	   case 'Y':
		    matrix[2][i+start]++;matrix[3][i+start]++;break;
	   case 'W':
		   matrix[0][i+start]++;matrix[3][i+start]++;break;
	   case 'S':
		   matrix[1][i+start]++;matrix[2][i+start]++;break;
	   case 'B':
		   matrix[1][i+start]++;matrix[2][i+start]++;matrix[3][i+start]++;break;
	   case 'V':
		   matrix[0][i+start]++;matrix[1][i+start]++;matrix[2][i+start]++;break;
	   case 'D':
		   matrix[0][i+start]++;matrix[1][i+start]++;matrix[3][i+start]++;break;
	   case 'H':
		   matrix[0][i+start]++;matrix[2][i+start]++;matrix[3][i+start]++;break;
	   case 'N':
		   matrix[0][i+start]++;matrix[1][i+start]++;matrix[2][i+start]++;matrix[3][i+start]++;break;
	   default:
		   break;
	   }
	   i++;
   }
   for(size_t i=0;i<lens+start+end;i++)
   {
	   double buff=matrix[0][i]+matrix[1][i]+matrix[2][i]+matrix[3][i];
	   for(int j=0;j<4;j++) matrix[j][i]=(matrix[j][i]/buff+beta*bg[j])/(1.0+beta);	   
   }
   for(size_t i=0;i<lens+start+end;i++)
   {
	   for(int j=0;j<4;j++) matrix1[j][i]+=matrix[j][i];
   }

}
double Cor(int sign,double matrix[3][4][50],size_t lens,int& signal)
{
	double r=0,temp=0;double buff[2][5];
	int pos=sign;
	if(pos<0) pos=0;
	for(size_t i=pos;i<(lens+pos);i++)
	{
		memset(buff,0,sizeof(buff));
		for(int j=0;j<4;j++)
		{
			buff[0][0]+=matrix[0][j][i]*matrix[1][j][i-sign];
			buff[1][0]+=matrix[0][j][i]*matrix[2][j][i-sign];
			buff[0][1]+=matrix[0][j][i];buff[1][1]+=matrix[0][j][i];
			buff[0][2]+=matrix[1][j][i-sign];buff[1][2]+=matrix[2][j][i-sign];
			buff[0][3]+=pow(matrix[0][j][i],2);buff[1][3]+=pow(matrix[0][j][i],2);
			buff[0][4]+=pow(matrix[1][j][i-sign],2);buff[1][4]+=pow(matrix[2][j][i-sign],2);
		}
		r+=(4*buff[0][0]-buff[0][1]*buff[0][2])/(sqrt(4*buff[0][3]-pow(buff[0][1],2))*sqrt(4*buff[0][4]-pow(buff[0][2],2)));
		temp+=(4*buff[1][0]-buff[1][1]*buff[1][2])/(sqrt(4*buff[1][3]-pow(buff[1][1],2))*sqrt(4*buff[1][4]-pow(buff[1][2],2)));
	}
	if(temp>r) {r=temp;signal=1;}
	r=r/(double)lens;
      return r;
}
inline double inmin(double a,double b,double c,double d)
{
	double temp=a;
	if(b<temp) temp=b;
	if(c<temp) temp=c;
	if(d<temp) temp=d;
	return temp;
}
inline double getscore(double matrix[4][50],char *motif)
{
	double score=0;
    size_t lens=strlen(motif);
	for(size_t i=0;i<lens;i++)
	{
	 switch(motif[i])
	   {
	   case 'A':
		   score+=matrix[0][i];break;
	   case 'G':
		   score+=matrix[1][i];break;
	   case 'C':
		   score+=matrix[2][i];break;
	   case 'T':
		   score+=matrix[3][i];break;
	   case 'N':
		   score=score+inmin(matrix[0][i],matrix[1][i],matrix[2][i],matrix[3][i]);break;
	   default:
		   break;
	   }
   }
	return score;
}
void scanSe(string seq,double matrix[4][50],int k,int sign,int &C,int &N,double maxscore)
{
  char motif[50],Trans[50];
  size_t lens=seq.size();
  double temp,score;
  for(size_t i=0;i<lens-k-1;)
  {
	  int j=0;
	  for(j=0;j<k;j++)
	  {
		  if(seq[i+j]=='N') { i=i+j+1;break;}
		  motif[j]=seq[i+j];
	  }
      motif[j]=NULL;
	  i++;
	  if(strlen(motif)!=k) continue;
	  Translate(motif,Trans);
      score=getscore(matrix,motif);
	  temp=getscore(matrix,Trans);
	  if(temp>score) score=temp;
	  if((score/maxscore)>=0.8)
	  {
		  if(sign==0) C++;
		  else N++;
	  }
  }
}
double Correlation(char *motif1,char *motif2,double beta,int& site)
{
   size_t lens=strlen(motif1);
   char* motif3=(char*)new char[lens+1];
   int signal=0;
   Translate(motif2,motif3);
   double matrix[3][4][50];double buff=0,temp[2];
   memset(matrix,0,sizeof(matrix));
   countMatrix(matrix[0],motif1,beta,1);
   countMatrix(matrix[1],motif2,beta,1);
   countMatrix(matrix[2],motif3,beta,1);  
   for(size_t i=0;i<lens-5;i++)
   {
	   int sign[2]={0,0},sitebuff=0-i;
       temp[0]=Cor(i,matrix,lens-i,sign[0]);
	   temp[1]=Cor(0-i,matrix,lens-i,sign[1]);
	   if(temp[0]>temp[1]){temp[1]=temp[0];sign[1]=sign[0];sitebuff=i;}
	   if(buff<temp[1]) {buff=temp[1];signal=sign[1];site=sitebuff;}
   }//end for
  if(signal==1) strcpy(motif2,motif3);
  delete [] motif3;
   return buff;
}
void merge(struct Result *head,int start,int end,double beta,double matrix[4][50])
{
    struct Result *cur=head;
	size_t lens=strlen(cur->motif);
	while(cur!=NULL)
	{
		Trans(matrix,cur->motif,cur->site-start,end-cur->site,beta);
		cur=cur->next;
	}
	for(size_t i=0;i<lens-start+end;i++)
   {
	   double buff=matrix[0][i]+matrix[1][i]+matrix[2][i]+matrix[3][i];
	   for(int j=0;j<4;j++) matrix[j][i]=matrix[j][i]/buff;	   
   }
}
char Getbp(int *sign)
{
	char seq[5];int count=0;
	char *bp[15]={"A","C","G","T","AC","AG","AT","GC","CT","GT","GCT","AGT","ACT","AGC","AGCT"};
	char buff[15]={'A','C','G','T','M','R','W','S','Y','K','B','D','H','V','N'};
	for(int i=0;i<4;i++)
	{
		if(sign[i]==1)
		{
			switch(i)
			{
			case 0:
				seq[count]='A';break;
			case 1:
				seq[count]='G';break;
			case 2:
				seq[count]='C';break;
			case 3:
				seq[count]='T';break;
			default:
				break;

			}
			count++;
		}
	}
	seq[count]=NULL;
	char re='X';
	for(int i=0;i<15;i++)
	{
		if(strcmp(seq,bp[i])==0)
		{
			re=buff[i];
			break;
		}
	}
   return re;
}
void Getconsensus(double matrix[4][50],char *consensus,int lens)
{
   double buff1=0;
   for(int i=0;i<4;i++)
	   buff1+=log(bg[i])/log(double(2))*bg[i];
   for(int i=0;i<lens;i++)
   {
	   int sign[4]={0,0,0,0};
	   double buff=0,temp;
	   for(int j=0;j<4;j++)
	   {
	     temp=log(matrix[j][i])/log(double(2))*matrix[j][i];
		 if((log(matrix[j][i]/bg[j])/log(double(2)))>0) sign[j]=1; 
		 buff+=temp;
	   }
	   if((buff-buff1)<0.25) consensus[i]='N';
	   else consensus[i]=Getbp(sign);	   
   }
   consensus[lens]=NULL;
}
void Tranmatrix(double matrix[4][50],int length,double transmatrix[3][4][50])
{
     for(int i=0;i<length;i++)
	 {
		 for(int j=0;j<4;j++)
			 transmatrix[1][j][i]=matrix[j][i];
		 transmatrix[2][3][i]=matrix[0][length-1-i];
         transmatrix[2][2][i]=matrix[1][length-1-i];
		 transmatrix[2][1][i]=matrix[2][length-1-i];
		 transmatrix[2][0][i]=matrix[3][length-1-i];
	 }
}
void normalMatrix(double matrix[4][50],double beta,int length)
{
	for( int i=0;i<length;i++)
   {
	   double buff=matrix[0][i]+matrix[1][i]+matrix[2][i]+matrix[3][i];
	   if(buff==0) 
		   for(int j=0;j<4;j++)
			   matrix[j][i]=bg[j];
	   else 
	      for(int j=0;j<4;j++) matrix[j][i]=(matrix[j][i]/buff+beta*bg[j])/(1.0+beta);	   
   }
}

