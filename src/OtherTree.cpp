#include "stdafx.h"
#include "OtherTree.h"
void constructTree(tree* root,int len,int count,int sign)
{
  
     switch(sign)
     {
      case 0:
          root->Bp='N';break;
      case 1:
          root->Bp='A';break;
      case 2:
  	      root->Bp='G';break;
      case 3:
    	  root->Bp='C';break;
      case 4:
    	  root->Bp='T';break;
      default:
    	 break;
     }
     root->countB=0;
     root->countF=0;
	 root->sign=0;
     if(count<2)
     {
		      root->child1[0]=new tree;
              root->child2[0]=new tree;
              root->child3[0]=new tree;
              root->child4[0]=new tree;
              constructTree(root->child1[0],len,count+1,1);
              constructTree(root->child2[0],len,count+1,2);
              constructTree(root->child3[0],len,count+1,3);
              constructTree(root->child4[0],len,count+1,4);
         for(int i=1;i<15;i++)
        {
              root->child1[i]=NULL;
              root->child2[i]=NULL;
              root->child3[i]=NULL;
              root->child4[i]=NULL;
                 
        }
    }
	else if(count<4)
	{
		for(int i=0;i<(15-len);i++)
        {
              root->child1[i]=new tree;
              root->child2[i]=new tree;
              root->child3[i]=new tree;
              root->child4[i]=new tree;
              constructTree(root->child1[i],i,count+1,1);
              constructTree(root->child2[i],i,count+1,2);
              constructTree(root->child3[i],i,count+1,3);
              constructTree(root->child4[i],i,count+1,4);       
        }
		   for(int i=15-len;i<15;i++)
        {
              root->child1[i]=NULL;
              root->child2[i]=NULL;
              root->child3[i]=NULL;
              root->child4[i]=NULL;
                 
        }

	 }
	else if(count<6)
	{
		      root->child1[0]=new tree;
              root->child2[0]=new tree;
              root->child3[0]=new tree;
              root->child4[0]=new tree;
              constructTree(root->child1[0],len,count+1,1);
              constructTree(root->child2[0],len,count+1,2);
              constructTree(root->child3[0],len,count+1,3);
              constructTree(root->child4[0],len,count+1,4);
       
         for(int i=1;i<15;i++)
        {
              root->child1[i]=NULL;
              root->child2[i]=NULL;
              root->child3[i]=NULL;
              root->child4[i]=NULL;
                 
        }

	}
	else
	{
          for(int i=0;i<15;i++)
        {
              root->child1[i]=NULL;
              root->child2[i]=NULL;
              root->child3[i]=NULL;
              root->child4[i]=NULL;

        }
	}
}

inline tree * getadd(char Bp,int pos,tree *root)
{
	 struct tree *curbuff=NULL;
    switch(Bp)
	{
	case 'A':
		curbuff=root->child1[pos];break;
	case 'G':
		curbuff=root->child2[pos];break;
	case 'C':
		curbuff=root->child3[pos];break;
	case 'T':
		curbuff=root->child4[pos];break;
	default:
		break;
	}//end switch
	return curbuff;

}
tree * getroot(char *motif,tree *root)
{
     struct tree *curbuff;
	 curbuff=getadd(motif[0],0,root);
     curbuff=getadd(motif[1],0,curbuff);
	 curbuff=getadd(motif[2],0,curbuff);
	 size_t lens=strlen(motif);
	 curbuff=getadd(motif[lens-3],lens-6,curbuff);
	 curbuff=getadd(motif[lens-2],0,curbuff);
	 curbuff=getadd(motif[lens-1],0,curbuff);
	 return curbuff;
	
}
void TranNumber(char bp,int& i)
{
	 switch(bp)
	 {
	 case 'A':
		 i=0;break;
	 case 'G':
		 i=1;break;
	 case 'C':
		 i=2;break;
	 case 'T':
         i=3;break;
	 case 'X':
		 i=4;break;
	 case 'N':
		 i=5;break;
	 default:
		 break;
	 }
}
 void GetvalueAll(string sequence,struct Bgpro &bglog)
{
	size_t lens=sequence.size();
	Bseqscur->seq=new char[lens+38];
    for(size_t j=0;j<7;j++) Bseqscur->seq[j]='N';
    for(size_t j=7;j<(lens+7);j++)
	{
		Bseqscur->seq[j]=sequence[j-7];
	}
    for(size_t j=7+lens;j<(lens+37);j++) Bseqscur->seq[j]='N';
    Bseqscur->seq[lens+37]=NULL;
	Bseqscur->len=lens+37;
	for(size_t i=0;i<lens;i++)
	{
		int j,k,m,n;
		j=k=m=n=6;
		if(i==0) 
		{
			TranNumber(sequence[i],j);
			if(j<4) 
			{
				bglog.Obp[j]++;
				if(!MICRO)bglog.Obp[3-j]++;
			}
		}
		else if(i==1)
		{
			TranNumber(sequence[i],j);
		    TranNumber(sequence[i-1],k);
            if(j<4)
			{
				bglog.Obp[j]++;
				if(!MICRO)bglog.Obp[3-j]++;
			}
		    if(j<4&&k<4) 
			{
				bglog.Tbp[k][j]++;
				if(!MICRO)bglog.Tbp[3-j][3-k]++;
			}
		}
		else if(i==2)
		{
            TranNumber(sequence[i],j);
			TranNumber(sequence[i-1],k);
			TranNumber(sequence[i-2],m);
            if(j<4) 
			{
				bglog.Obp[j]++;
				if(!MICRO)bglog.Obp[3-j]++;
			}
			if(j<4&&k<4)
			{
				bglog.Tbp[k][j]++;
                if(!MICRO)bglog.Tbp[3-j][3-k]++;
			}
			if(j<4&&k<4&&m<4)
			{
				bglog.Thbp[m][k][j]++;
				if(!MICRO)bglog.Thbp[3-j][3-k][3-m]++;
			}
		}
		else 
		{
			TranNumber(sequence[i],j);
			TranNumber(sequence[i-1],k);
			TranNumber(sequence[i-2],m);
			TranNumber(sequence[i-3],n);
			if(j<4) 
			{
				bglog.Obp[j]++;
                if(!MICRO)bglog.Obp[3-j]++;
			}
			if(j<4&&k<4)
			{
				bglog.Tbp[k][j]++;
                if(!MICRO)bglog.Tbp[3-j][3-k]++;
			}
			if(j<4&&k<4&&m<4)
			{
				bglog.Thbp[m][k][j]++;
				if(!MICRO)bglog.Thbp[3-j][3-k][3-m]++;
			}
			if(j<4&&k<4&&m<4&&n<4) 
			{
				bglog.Fbp[n][m][k][j]++;
				if(!MICRO)bglog.Fbp[3-j][3-k][3-m][3-n]++;
			}
		}
	}

}

 void Getvalue(char *sequence,struct tree *root,int sign)
{
	int posi=0;
	if(sequence[posi]!='N'&&sequence[posi+1]!='N'&&sequence[posi+2]!='N')
	{
		struct tree *curbuff;
		curbuff=getadd(sequence[posi],0,root);
		curbuff=getadd(sequence[posi+1],0,curbuff);
		curbuff=getadd(sequence[posi+2],0,curbuff);
		for(int j=posi+3;j<(posi+18);j++)
		 {
			if(sequence[j]!='N'&&sequence[j+1]!='N'&&sequence[j+2]!='N')
			{
				struct tree *curbuff1=curbuff;
				curbuff1=getadd(sequence[j],j-posi-3,curbuff);
				curbuff1=getadd(sequence[j+1],0,curbuff1);
				curbuff1=getadd(sequence[j+2],0,curbuff1);
				if(sign==0) curbuff1->countF+=1;
				else curbuff1->countB+=1;
			}
		  }//end for r
	}
}
void DistroyTree(tree *root,int len)
{
	 int signal=0;
	 for(int i=0;i<(15-len);i++)
     { 
             if(root->child1[i]!=NULL&&root->child2[i]!=NULL&&root->child3[i]!=NULL&&root->child4[i]!=NULL)
             {
                DistroyTree(root->child1[i],i);
                DistroyTree(root->child2[i],i);
                DistroyTree(root->child3[i],i);
                DistroyTree(root->child4[i],i);
             }
             else
             {
                 signal+=1;
                 }
     }
	 if(signal!=(15-len)){ delete root;root=NULL;}
	 if(signal==15)
	 {
		 delete root;
		 root=NULL;
	 }
}
void TranslateE(char *motif,char *Trans)
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
void getZscore(char* seq,tree *root,int number[2],int count,char *seqbuff,int sign)
{
	if(count<6)
	{
		int posi=0,posi1=count,len=strlen(seq);
		if(count==3) posi=len-6;
		if(count>2) posi1=len-3+count-3;
		switch(seq[posi1])
		{
		   case 'A':
			   seqbuff[count]='A';getZscore(seq,root->child1[posi],number,count+1,seqbuff,sign);break;
		   case 'G':
			   seqbuff[count]='G';getZscore(seq,root->child2[posi],number,count+1,seqbuff,sign);break;
		   case 'C':
			   seqbuff[count]='C';getZscore(seq,root->child3[posi],number,count+1,seqbuff,sign);break;
		   case 'T':
			   seqbuff[count]='T';getZscore(seq,root->child4[posi],number,count+1,seqbuff,sign);break;
		   case 'M':
			   seqbuff[count]='A';getZscore(seq,root->child1[posi],number,count+1,seqbuff,sign);seqbuff[count]='C';getZscore(seq,root->child3[posi],number,count+1,seqbuff,sign);break;
		   case 'K':
			   seqbuff[count]='G';getZscore(seq,root->child2[posi],number,count+1,seqbuff,sign);seqbuff[count]='T';getZscore(seq,root->child4[posi],number,count+1,seqbuff,sign);break;
		   case 'R':
			   seqbuff[count]='A';getZscore(seq,root->child1[posi],number,count+1,seqbuff,sign);seqbuff[count]='G';getZscore(seq,root->child2[posi],number,count+1,seqbuff,sign);break;
		   case 'Y':
			   seqbuff[count]='C';getZscore(seq,root->child3[posi],number,count+1,seqbuff,sign);seqbuff[count]='T';getZscore(seq,root->child4[posi],number,count+1,seqbuff,sign);break;
		   case 'W':
			   seqbuff[count]='A';getZscore(seq,root->child1[posi],number,count+1,seqbuff,sign);seqbuff[count]='T';getZscore(seq,root->child4[posi],number,count+1,seqbuff,sign);break;
		   case 'S':
			   seqbuff[count]='G';getZscore(seq,root->child2[posi],number,count+1,seqbuff,sign);seqbuff[count]='C';getZscore(seq,root->child3[posi],number,count+1,seqbuff,sign);break;
		   case 'B':
			   seqbuff[count]='C';getZscore(seq,root->child3[posi],number,count+1,seqbuff,sign);seqbuff[count]='G';getZscore(seq,root->child2[posi],number,count+1,seqbuff,sign);seqbuff[count]='T';getZscore(seq,root->child4[posi],number,count+1,seqbuff,sign);break;
		   case 'V':
			   seqbuff[count]='A';getZscore(seq,root->child1[posi],number,count+1,seqbuff,sign);seqbuff[count]='C';getZscore(seq,root->child3[posi],number,count+1,seqbuff,sign);seqbuff[count]='G';getZscore(seq,root->child2[posi],number,count+1,seqbuff,sign);break;
		   case 'D':
			   seqbuff[count]='A';getZscore(seq,root->child1[posi],number,count+1,seqbuff,sign);seqbuff[count]='G';getZscore(seq,root->child2[posi],number,count+1,seqbuff,sign);seqbuff[count]='T';getZscore(seq,root->child4[posi],number,count+1,seqbuff,sign);break;
		   case 'H':
			   seqbuff[count]='A';getZscore(seq,root->child1[posi],number,count+1,seqbuff,sign);seqbuff[count]='C';getZscore(seq,root->child3[posi],number,count+1,seqbuff,sign);seqbuff[count]='T';getZscore(seq,root->child4[posi],number,count+1,seqbuff,sign);break;
		   case 'N':
			   seqbuff[count]='A';getZscore(seq,root->child1[posi],number,count+1,seqbuff,sign);seqbuff[count]='G';getZscore(seq,root->child2[posi],number,count+1,seqbuff,sign);seqbuff[count]='C';getZscore(seq,root->child3[posi],number,count+1,seqbuff,sign);seqbuff[count]='T';getZscore(seq,root->child4[posi],number,count+1,seqbuff,sign);break;
		   default:
			   break;
		 }
	}//end if
	else 
	{
		root->sign=2;
		seqbuff[count]=NULL;
		size_t lens=strlen(seq);
		char motif[50],trans[50];
		for(int i=0;i<3;i++) motif[i]=seqbuff[i];
		for(size_t i=3;i<(3+lens-6);i++) motif[i]='N';
		for(int i=lens-3;i<lens;i++) motif[i]=seqbuff[i-lens+3+3];
		motif[lens]=NULL;
		TranslateE(motif,trans);
		tree *cur=getroot(trans,curroot);
		if(sign==0)
		{
			
			if(cur->sign!=2&&!MICRO)
			{
				number[0]=number[0]+root->countF+cur->countF;
				number[1]=number[1]+root->countB+cur->countB;
			}
			else
			{
				number[0]+=root->countF;
				number[1]+=root->countB;
			}
			root->sign=1;
		}
		else {
			root->sign=3;
			if(!MICRO) cur->sign=3;
		}
	}
}