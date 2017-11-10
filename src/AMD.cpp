// AMD.cpp : The main program
//

#include "stdafx.h"
#include "OtherTree.h"
#include "math.h"
#include <iostream>
#include <fstream>
#include <string>
#include "Correlation.h"
#include "countBp.h"
#include "time.h"

using namespace std;
#define MARK1 3*1024
#define MARK2 3*1024*4
#define MARK3 3*1024*16
#define SMALL 0.00000025


struct Gsize NC;
struct Seqstru
{
	char *seq;
	double *seqlog;
	double *transeqlog;
	size_t len;
	Seqstru *next;
};
int pos=0,Topcut=50,heap=1,Nnumber=6,markov=0,sortsign=0,MICRO=0;
double Foldcut=1.1,CO=0.6;
struct Result *choose,*Rhead;
struct Bgpro bglog;
const char *filename2,*filename1;
//double bg[4]={0.31,0.19,0.19,0.31};
tree *curroot;
struct Seqstru *seqHead,*seqscur,*seqsdel;
struct Bseqstru *BseqHead,*Bseqscur,*Bseqsdel;
double bg[4]={0.295,0.205,0.205,0.295};

void TrantoPoint(Result *choose,int size);
void Redundancy(Result *head);
void GetLong(Result *choose);
void getMscore(Result *cur);
void GetLongCore(Result *head,Gsize NC);
void getBglog(char *seq,double *score);
void CalBglog(struct Bgpro &bglog,int &markov);
void travalTree(tree *root,int len,char* seq,int sign,Gsize NC);
inline void sortbuff(double buffm[6],int sign,int max);
inline void getMbglo(double *slog,int i,int k,double *bglo);
void sort(struct Result *Top,int count,int maxpos);
void scanSe(string seq,struct tree *root,int sign);
void seqTranto(struct Seqstru *Head,string seq);
inline void Getmatrix(Result *cur,double score);
void updataCore(Result *choose,tree *root,Gsize NC);
Result* reRedun(Result *head);
void resort(Result *head);
void getOut(Result *head);
void filledMatrix(double matrixbuff[4][50],double matrix[4][50],int left,int right,int len);
void GetFontseq(Result *head);
int main(int argc, char* argv[])
{
   time_t procstart,procend;
   procstart=time(NULL);
   register int argF,argT,argFC,argB,argCO,argMI;
   const char *filename;
   struct tree *root=new tree;
   memset(&bglog,0,sizeof(bglog));
   seqscur=new Seqstru;
   seqHead=seqscur;
   Bseqscur=new Bseqstru;
   BseqHead=Bseqscur;

   ifstream fgfp,bgfp;
   char seq[20];
   string Gsequence,buffseq;
   int Mlens=0,signal=0;

   //read in arguments from command line
   argF=argT=argFC=argB=argCO=argMI=50;
   for(int i=1;i<argc;i++)
   {
	   if(strcmp(argv[i],"-F")==0)//Fg file
		   argF=i;
	   if(strcmp(argv[i],"-B")==0)//Background file
		   argB=i;
	   if(strcmp(argv[i],"-T")==0)//Top number
		   argT=i;
	   if(strcmp(argv[i],"-CO")==0)//Correlation
		   argCO=i;
	   if(strcmp(argv[i],"-FC")==0)//Cutoff of Fold Change
		   argFC=i;
	   if(strcmp(argv[i],"-MI")==0)//Micro RNA
		   argMI=i;
   }//end for
   if(argF==50||argB==50)
   {
	   cout<<"Usage:\n"<<"\tRequired arguments:\n"<<"\t[-F] [Foreground file]\n"<<"\t[-B] [Background file]\n"<<endl;
	   cout<<"\tOptional arguments:\n"<<"\t[-MI] [Mirco RNA switch]\n"<<"\t[-T] [Number of top k-mers considered,50 is default]\n"<<"\t[-CO] [Similarity cutoff to eliminate redundancy,0.6 is default]\n"<<"\t[-FC] [Fold change cutoff,1.2 is default]"<<endl;
	   exit(0);
   }



   //Foreground File
   if(argF!=50)
        filename2=argv[argF+1];
   //correlation
    if(argCO!=50)
		CO=atof(argv[argCO+1]);

   //Background File
   if(argB!=50)
	   filename1=argv[argB+1];

   //cutoffs
   if(argFC!=50)
	   Foldcut=atof(argv[argFC+1]);
   if(argT!=50)
	   Topcut=atoi(argv[argT+1]);
   if(argMI!=50)
	   MICRO=1;

   // Construct Foreground File
   NC.Fsize=countBp(filename2);
   NC.Bsize=countBp(filename1);


   //Construct tree

   constructTree(root,14,0,0);
   curroot=root;
   fgfp.open(filename2,ios::in);
   if(fgfp.fail()) {cout<<"Error!Can not open the Foreground result file!"<<endl;exit(0);}
   bgfp.open(filename1,ios::in);
   if(bgfp.fail()) {cout<<"Error!Can not open the Background result file!"<<endl;exit(0);}
   cout<<"Scanning background sequences!"<<endl;
   while(getline(bgfp,buffseq))
   {
		if(buffseq[0]=='>')
		{
			if(signal==1)  {
                           if(Gsequence.size()>40)
                             scanSe(Gsequence,root,1);
                            Gsequence.clear();
                            }
			else { Gsequence.clear();signal=1;}
		}
		else{
			Gsequence+=buffseq;}
   }
   if(Gsequence.size()>40)
    scanSe(Gsequence,root,1);
   Gsequence.clear();
    Bseqsdel->next=NULL;
   delete Bseqscur;Bseqscur=NULL;
   //calcuate bglog
    CalBglog(bglog,markov);

    for(int i=0;i<4;i++)
	{
		bg[i]=exp(bglog.Obp[i]);
	}


  cout<<"Scanning foreground sequences!"<<endl;
   signal=0;
   while(getline(fgfp,buffseq))
   {
		if(buffseq[0]=='>')
		{
			if(signal==1)  {
                           if(Gsequence.size()>40)
                             scanSe(Gsequence,root,0);
                           Gsequence.clear();
                           }
			else { Gsequence.clear();signal=1;}
		}
		else{
			Gsequence+=buffseq;}
   }
   if(Gsequence.size()>40)
    scanSe(Gsequence,root,0);
    Gsequence.clear();
   fgfp.close();bgfp.close();
   seqsdel->next=NULL;
   delete seqscur;seqscur=NULL;

   //Judge mode
   char Rfile1[500],Rfile2[500];
   strcpy(Rfile1,filename2);
   strcpy(Rfile2,filename2);
   strcat(Rfile1,".Matrix");
   strcat(Rfile2,".Details");
   FILE *Moti=fopen(Rfile1,"wt"),*Tmoti=fopen(Rfile2,"wt");
   choose=new Result[Topcut+1];
   Rhead=choose;
   cout<<"Get top kmers!"<<endl;
   cout<<"Number\tMotif Consensus\tScore"<<endl;

   travalTree(root,14,seq,0,NC);//get top
   if(pos<Topcut)
        sort(choose,0,pos);

   for(int i=0;i<pos;i++) cout<<"Top--"<<i<<"\t"<<choose[i].motif<<"\t"<<choose[i].Zscore<<endl;
   if(pos==0) { cout<<"Cannot find !"<<endl;exit(0);}

   TrantoPoint(choose,Topcut);
   Redundancy(choose);
   Result *curbuff=Rhead;
   cout<<"Updata cores!"<<endl;
   cout<<"Initial cores\t"<<"Degenerate cores"<<endl;
   while(curbuff!=NULL)
   {
      cout<<curbuff->motif<<"\t";
	  updataCore(curbuff,root,NC);
	  cout<<curbuff->motif<<endl;
	  curbuff=curbuff->next;
   }

   Redundancy(Rhead);
    cout<<"Extending cores!"<<endl;
   cout<<"Degenerate cores\t"<<"Extended cores"<<endl;
   curbuff=Rhead;
   while(curbuff!=NULL)
   {
	  int number[2];
	  char seqbuff[50];
	  getZscore(curbuff->motif,root,number,0,seqbuff,1);
	  curbuff=curbuff->next;
   }

   GetFontseq(Rhead);

   GetLongCore(Rhead,NC);
   //////////////////////////////

   ////////////////////////////
    cout<<"Refining seed motif!"<<endl;
   GetLong(Rhead);
   DistroyTree(root,14);
   int count=1;
   fprintf(Tmoti,"Motif\tCore\tConsensus\tSegmentNum\tMAPscore\n");
   double buff1=0;
	 for(int i=0;i<4;i++)
	   buff1+=log(bg[i])/log(double(2))*bg[i];
 //////////////////////////////////////////////////
  Rhead=reRedun(Rhead);

  resort(Rhead);

  cout<<"Get output!"<<endl;
  getOut(Rhead);

//////////////////////////////////
  fprintf(Moti,"Motif\tA\tC\tG\tT\n");
   while(Rhead!=NULL)
   {
    if(count>=10) break;
    char cons[50];
	if(Rhead->segment==0){Rhead=Rhead->next; continue;}
	Getconsensus(Rhead->matrix[0],cons,Rhead->len);
	if(strlen(cons)>30||strlen(cons)<6) {Rhead=Rhead->next;continue;}
	fprintf(Moti,"Motif%d:\n",count);
	for(int i=0;i<Rhead->len;i++)
	{
		fprintf(Moti,"%d\t",i+1);
		fprintf(Moti,"%f\t",Rhead->matrix[0][0][i]);fprintf(Moti,"%f\t",Rhead->matrix[0][2][i]);
		fprintf(Moti,"%f\t",Rhead->matrix[0][1][i]);fprintf(Moti,"%f\n",Rhead->matrix[0][3][i]);
		double temp=Rhead->matrix[0][1][i];
		Rhead->matrix[0][1][i]=Rhead->matrix[0][2][i];
		Rhead->matrix[0][2][i]=temp;
	}
	fprintf(Tmoti,"Motif%d\t%s\t%s\t%d\t%f\n",count,Rhead->motif,cons,int(Rhead->segment),Rhead->Mscore);

	count++;
	cout<<cons<<endl;

	Rhead=Rhead->next;
   }

   delete [] choose;

   //free memory


   //Consenses

  seqsdel=seqHead;
  Bseqsdel=BseqHead;
  while(seqsdel!=NULL)
  {
	  Seqstru *delcur=seqsdel;
	  seqsdel=seqsdel->next;
	  delete [] delcur->seq;
	  delete [] delcur->seqlog;
	  delete [] delcur->transeqlog;
      delete delcur;
  }
  while(Bseqsdel!=NULL)
  {
	  Bseqstru *delcur=Bseqsdel;
	  Bseqsdel=Bseqsdel->next;
	  delete [] delcur->seq;
      delete delcur;
  }
  procend=time(NULL);
  fprintf(Tmoti,"Runtime:%f",difftime(procend,procstart));
  fclose(Moti);fclose(Tmoti);
  return 0;
}
char tranTobp(int bpos)
{
	char bp='N';
	switch(bpos)
	{
	case 0:
		bp='A';break;
	case 1:
		bp='G';break;
	case 2:
		bp='C';break;
	case 3:
		bp='T';break;
	default:break;
	}
	return bp;
}
void tranToseq(char motif[50],int pos1,int pos2)
{
	for(int i=0;i<4;i++)
	{
		for(int j=0;j<4;j++)
		{
			for(int k=0;k<4;k++)
			{
				int sum=i*16+j*4+k;
				if(sum==pos1)
				{
                   motif[0]=tranTobp(i);motif[1]=tranTobp(j);motif[2]=tranTobp(k);
				}
				if(sum==pos2)
				{
					motif[3]=tranTobp(i);motif[4]=tranTobp(j);motif[5]=tranTobp(k);
				}
			}
		}
	}
	motif[6]=NULL;
}

inline void computeSum(double sum[2][4][50],Gsize NC,char *seq,Result *cur,int i)
{
	char motif[50];
	strcpy(motif,seq);
	char buff[15]={'A','G','C','T','R','M','W','S','K','Y','V','D','H','B','N'};
	double sumary1[14]={sum[0][0][i],sum[0][1][i],sum[0][2][i],sum[0][3][i],sum[0][0][i]+sum[0][1][i],sum[0][0][i]+sum[0][2][i],sum[0][0][i]+sum[0][3][i],
		sum[0][1][i]+sum[0][2][i],sum[0][1][i]+sum[0][3][i],sum[0][2][i]+sum[0][3][i],sum[0][0][i]+sum[0][1][i]+sum[0][2][i],sum[0][0][i]+sum[0][1][i]+sum[0][3][i],
		sum[0][0][i]+sum[0][2][i]+sum[0][3][i],sum[0][1][i]+sum[0][2][i]+sum[0][3][i]};
    double sumary2[14]={sum[1][0][i],sum[1][1][i],sum[1][2][i],sum[1][3][i],sum[1][0][i]+sum[1][1][i],sum[1][0][i]+sum[1][2][i],sum[1][0][i]+sum[1][3][i],
		sum[1][1][i]+sum[1][2][i],sum[1][1][i]+sum[1][3][i],sum[1][2][i]+sum[1][3][i],sum[1][0][i]+sum[1][1][i]+sum[1][2][i],sum[1][0][i]+sum[1][1][i]+sum[1][3][i],
		sum[1][0][i]+sum[1][2][i]+sum[1][3][i],sum[1][1][i]+sum[1][2][i]+sum[1][3][i]};

   for(int j=0;j<14;j++)
   {
	   char motifbuff[50];
	   strcpy(motifbuff,seq);

	   double E=NC.Fsize/(double)NC.Bsize*sumary2[j],F;
       double p=sumary2[j]/(double)NC.Bsize;
	   if(E!=0)
	   {
		   F=sumary1[j]/E;
		   double score=(sumary1[j]-E)/sqrt(E);
		   if(score>cur->Zscore&&F>cur->Fold)
		   {
			   cur->Zscore=score;
			   cur->F=F;
			   motifbuff[i]=buff[j];
			   strcpy(motif,motifbuff);
		   }
	   }
   }
   strcpy(seq,motif);
}

inline void computeSum1(int sum1[2][14][14][4][4],Gsize NC,char *seq,Result *cur,int i,int k)
{
	 char motif[50];
	 strcpy(motif,seq);
	  char buff[15]={'A','G','C','T','R','M','W','S','K','Y','V','D','H','B','N'};
	   int sumary1[4][14],sumary2[4][14];
	   for(int j=0;j<4;j++)
	   {
		   int sumbuff1[14]={sum1[0][i][k-i-1][j][0],sum1[0][i][k-i-1][j][1],sum1[0][i][k-i-1][j][2],sum1[0][i][k-i-1][j][3],sum1[0][i][k-i-1][j][0]+sum1[0][i][k-i-1][j][1],sum1[0][i][k-i-1][j][0]+sum1[0][i][k-i-1][j][2],sum1[0][i][k-i-1][j][0]+sum1[0][i][k-i-1][j][3],
			   sum1[0][i][k-i-1][j][1]+sum1[0][i][k-i-1][j][2],sum1[0][i][k-i-1][j][1]+sum1[0][i][k-i-1][j][3],sum1[0][i][k-i-1][j][2]+sum1[0][i][k-i-1][j][3],sum1[0][i][k-i-1][j][0]+sum1[0][i][k-i-1][j][1]+sum1[0][i][k-i-1][j][2],sum1[0][i][k-i-1][j][0]+sum1[0][i][k-i-1][j][1]+sum1[0][i][k-i-1][j][3],
			   sum1[0][i][k-i-1][j][0]+sum1[0][i][k-i-1][j][2]+sum1[0][i][k-i-1][j][3],sum1[0][i][k-i-1][j][1]+sum1[0][i][k-i-1][j][2]+sum1[0][i][k-i-1][j][3]};
		   memcpy(sumary1[j],sumbuff1,sizeof(sumbuff1));
           int sumbuff2[14]={sum1[1][i][k-i-1][j][0],sum1[1][i][k-i-1][j][1],sum1[1][i][k-i-1][j][2],sum1[1][i][k-i-1][j][3],sum1[1][i][k-i-1][j][0]+sum1[1][i][k-i-1][j][1],sum1[1][i][k-i-1][j][0]+sum1[1][i][k-i-1][j][2],sum1[1][i][k-i-1][j][0]+sum1[1][i][k-i-1][j][3],
			   sum1[1][i][k-i-1][j][1]+sum1[1][i][k-i-1][j][2],sum1[1][i][k-i-1][j][1]+sum1[1][i][k-i-1][j][3],sum1[1][i][k-i-1][j][2]+sum1[1][i][k-i-1][j][3],sum1[1][i][k-i-1][j][0]+sum1[1][i][k-i-1][j][1]+sum1[1][i][k-i-1][j][2],sum1[1][i][k-i-1][j][0]+sum1[1][i][k-i-1][j][1]+sum1[1][i][k-i-1][j][3],
			   sum1[1][i][k-i-1][j][0]+sum1[1][i][k-i-1][j][2]+sum1[1][i][k-i-1][j][3],sum1[1][i][k-i-1][j][1]+sum1[1][i][k-i-1][j][2]+sum1[1][i][k-i-1][j][3]};
		    memcpy(sumary2[j],sumbuff2,sizeof(sumbuff2));
	   }
       int tsum1[14][14],tsum2[14][14];
	   for(int j=0;j<14;j++)
	   {
		   int sumbuff1[14]={sumary1[0][j],sumary1[1][j],sumary1[2][j],sumary1[3][j],sumary1[0][j]+sumary1[1][j],sumary1[0][j]+sumary1[2][j],sumary1[0][j]+sumary1[3][j],
			   sumary1[1][j]+sumary1[2][j],sumary1[1][j]+sumary1[3][j],sumary1[2][j]+sumary1[3][j],sumary1[0][j]+sumary1[1][j]+sumary1[2][j],sumary1[0][j]+sumary1[1][j]+sumary1[3][j],
			   sumary1[0][j]+sumary1[2][j]+sumary1[3][j],sumary1[1][j]+sumary1[2][j]+sumary1[3][j]};
		   memcpy(tsum1[j],sumbuff1,sizeof(sumbuff1));
		   int sumbuff2[14]={sumary2[0][j],sumary2[1][j],sumary2[2][j],sumary2[3][j],sumary2[0][j]+sumary2[1][j],sumary2[0][j]+sumary2[2][j],sumary2[0][j]+sumary2[3][j],
			   sumary2[1][j]+sumary2[2][j],sumary2[1][j]+sumary2[3][j],sumary2[2][j]+sumary2[3][j],sumary2[0][j]+sumary2[1][j]+sumary2[2][j],sumary2[0][j]+sumary2[1][j]+sumary2[3][j],
			   sumary2[0][j]+sumary2[2][j]+sumary2[3][j],sumary2[1][j]+sumary2[2][j]+sumary2[3][j]};
		   memcpy(tsum2[j],sumbuff2,sizeof(sumbuff2));
	   }
	   for(int j=0;j<14;j++)
	   {
		   for(int m=0;m<14;m++)
		   {
			   char motifbuff[50];
			   strcpy(motifbuff,seq);
              double E=NC.Fsize/(double)NC.Bsize*(double)tsum2[j][m],F;
			  double p=(double)tsum2[j][m]/(double)NC.Bsize;
			   if(E!=0)
			   {
				   F=double(tsum1[j][m])/E;
				   double score=(tsum1[j][m]-E)/sqrt(E);
				   if(score>cur->Zscore&&F>cur->Fold)
				   {
					   cur->Zscore=score;
					   cur->F=F;
					   motifbuff[i]=buff[m];
					   motifbuff[k]=buff[j];
					   strcpy(motif,motifbuff);
				   }//end if
			   }
		   }//end for k
	   }//end for j
	   strcpy(seq,motif);
}
void Getcore(int matrix2[2][14][14][4][4],double matrix1[2][4][50],Result *head,Gsize NC,char *seq,int sign)
{
	size_t lens=strlen(head->motif);

	char seqbuff1[50],seqbuff2[50],motif[50];
	int gap=lens-6,mlens=13,start,end;
	if(gap%2==0) mlens=14;
	start=0;end=mlens;
	//int loop=(mlens-gap)/2;
	//if(loop>3)
	//{start=(mlens-gap)/2-3;end=(mlens-gap)/2+gap+3;}
    strcpy(seqbuff2,seq);
	for(int i=start;i<end;i++)
	{
		strcpy(motif,seqbuff1);strcpy(seqbuff1,seq);
		computeSum(matrix1,NC,seqbuff1,head,i);
		if(strcmp(seqbuff1,motif)!=0) strcpy(seqbuff2,seqbuff1);
	}
	if(sign)
	{
		for(int i=start;i<(end-1);i++)
		{
		   for(int j=i+1;j<end;j++)
		   {
			   strcpy(seqbuff1,seq);strcpy(motif,seqbuff1);
			   computeSum1(matrix2,NC,seqbuff1,head,i,j);
			   if(strcmp(seqbuff1,motif)!=0) strcpy(seqbuff2,seqbuff1);
		   }
		}
	}
	strcpy(seq,seqbuff2);
}
void ScanSeq(int matrix2[2][14][14][4][4],double matrix1[2][4][50],Result *head,double matrix[4][50])
{
	size_t lens=strlen(head->motif),slens=head->Fseq.size();
	int gap=lens-6,mlens=13;
	if(gap%2==0) mlens=14;
	for(size_t p=0;p<slens;)
	{
		int lenbuff=0;
		char motif[20];
		for(lenbuff=0;lenbuff<mlens;lenbuff++)
			motif[lenbuff]=head->Fseq[p+lenbuff];
		motif[mlens]=NULL;
		double score=getscore(matrix,motif);
		if(score!=double(mlens)){p+=mlens; continue;}
		countMatrix(matrix1[0],motif,0.001,0);
		for(int i=0;i<mlens-1;i++)
		{
			for(int j=i+1;j<mlens;j++)
			{
			   int x=6,y=6;
			   TranNumber(motif[i],x);
               TranNumber(motif[j],y);
			   if(x<4&&y<4) matrix2[0][i][j-i-1][x][y]++;
			}
		}
		p=p+mlens;
	}
	slens=head->Bseq.size();
   for(size_t p=0;p<slens;)
	{
		int lenbuff=0;
		char motif[20];
		for(lenbuff=0;lenbuff<mlens;lenbuff++)
			motif[lenbuff]=head->Bseq[p+lenbuff];
		motif[mlens]=NULL;
		double score=getscore(matrix,motif);
		if(score!=double(mlens)){p+=mlens; continue;}
		countMatrix(matrix1[1],motif,0.001,0);
		for(int i=0;i<mlens-1;i++)
		{
			for(int j=i+1;j<mlens;j++)
			{
			   int x=6,y=6;
			   TranNumber(motif[i],x);
               TranNumber(motif[j],y);
			   if(x<4&&y<4) matrix2[1][i][j-i-1][x][y]++;
			}
		}
		p=p+mlens;
   }
}
void GetLongCore(Result *head,Gsize NC)
{
	Result *cur=head;
	while(cur!=NULL)
	{
      cout<<cur->motif<<"\t";
	 char seqbuff[50],seq[50]="NNNNNNNNNNNNN",seqfont[50];
	 strcpy(seqfont,seq);
	 int gap=strlen(cur->motif)-6;
	 int mlens=13,start,end;
	 if(gap%2==0){ mlens=14;seq[13]='N';seq[14]=NULL;}
	 int left=int((mlens-gap)/2),right=mlens-left-gap;
	 int matrix2[2][14][14][4][4];double matrix1[2][4][50];
	 memset(cur->buffmatrix,0,sizeof(cur->buffmatrix));
	 countMatrix(cur->buffmatrix,seq,0.001,0);
	 memset(matrix1,0,sizeof(matrix1));memset(matrix2,0,sizeof(matrix2));
	 ScanSeq(matrix2,matrix1,cur,cur->buffmatrix);
	 cur->Fold=cur->F;
	 Getcore(matrix2,matrix1,cur,NC,seq,1);
	 cur->Fold=cur->F;
	 int signal=0;
	 if(NC.Fsize>MARK3) signal=1;
	 if(strcmp(seqfont,seq)!=0)
	 {
		 memset(cur->buffmatrix,0,sizeof(cur->buffmatrix));
		 countMatrix(cur->buffmatrix,seq,0.001,0);
		 memset(matrix1,0,sizeof(matrix1));memset(matrix2,0,sizeof(matrix2));
		 ScanSeq(matrix2,matrix1,cur,cur->buffmatrix);
		 Getcore(matrix2,matrix1,cur,NC,seq,signal);
	 }
	 cur->Fseq.clear();
	 cur->Bseq.clear();
	 for(int i=0;i<left;i++) seqbuff[i]=seq[i];
	 for(int i=left;i<(left+3);i++) seqbuff[i]=cur->motif[i-left];
	 for(int i=left+3;i<(left+3+gap);i++) seqbuff[i]=seq[i-left-3+left];
	 for(int i=left+3+gap;i<(left+3+gap+3);i++) seqbuff[i]=cur->motif[i-left-3-gap+3+gap];
	 for(int i=left+3+gap+3;i<(left+3+gap+3+right);i++) seqbuff[i]=seq[i-left-3-gap-3+left+gap];
	 seqbuff[left+3+gap+3+right]=NULL;
	 int newlens=left+3+gap+3+right;
	 start=-1;end=newlens+1;
	 char motifbuff[50];
	 for(int i=0;i<newlens;i++)
	 {
		 if(seqbuff[i]!='N') {start=i;break;}
	 }
	 for(int i=newlens-1;i>=0;i--)
	 {
		 if(seqbuff[i]!='N') {end=i;break;}
	 }
	 if(start!=-1&&end!=(1+newlens)&&start!=end)
	 {

		 for(int i=start;i<=end;i++)
		 {
			 motifbuff[i-start]=seqbuff[i];
		 }
		 motifbuff[end-start+1]=NULL;
		 strcpy(cur->motif,motifbuff);
	  }
      cout<<cur->motif<<endl;
	 cur=cur->next;
	}
}
void WriteSeq(Result *head,char *motif,char *seq,int sign)
{
	Result *cur=head;
	char trans[50],seqtrans[15];
	Translate(motif,trans);
    Translate(seq,seqtrans);
	while(cur!=NULL)
	{
		if(strlen(motif)!=strlen(cur->motif)) {cur=cur->next;continue;}
		double score1=getscore(cur->buffmatrix,motif),score2=getscore(cur->buffmatrix,trans),len=double(strlen(cur->motif));
		if(score1==score2) score2=-1;
		if(score1==len&&score1!=score2)
		{
			if(sign==0)
			{
				cur->Bseq+=seq;
			}
			if(sign==1)
			{
				cur->Fseq+=seq;
			}
		}
		if(score2==len&&score2!=score1&&!MICRO)
		{
			if(sign==0)
			{
				cur->Bseq+=seqtrans;
			}
			if(sign==1)
			{
				cur->Fseq+=seqtrans;
			}
		}
		cur=cur->next;

	}
}
inline void GetvalueFont(char *sequence,Result *head,int sign)
{
	int posi=7;
	if(sequence[posi]!='N'&&sequence[posi+1]!='N'&&sequence[posi+2]!='N')
	{
		char buff[15],motifbuff[30];
		struct tree *curbuff;
		curbuff=getadd(sequence[posi],0,curroot);
		curbuff=getadd(sequence[posi+1],0,curbuff);
		curbuff=getadd(sequence[posi+2],0,curbuff);
		for(int i=0;i<3;i++) motifbuff[i]=sequence[posi+i];
		for(int j=posi+3;j<(posi+18);j++)
		 {

			if(sequence[j]!='N'&&sequence[j+1]!='N'&&sequence[j+2]!='N')
			{
				struct tree *curbuff1=curbuff;
				curbuff1=getadd(sequence[j],j-posi-3,curbuff);
				curbuff1=getadd(sequence[j+1],0,curbuff1);
				curbuff1=getadd(sequence[j+2],0,curbuff1);
				if(curbuff1->sign==3)
				{
					int gap=j-posi-3,left=int((14-gap)/2),right=14-gap-left;
					if(gap%2!=0) right=13-gap-left;
					for(int i=3;i<(3+gap);i++) motifbuff[i]='N';
					for(int i=3+gap;i<(6+gap);i++) motifbuff[i]=sequence[j+i-3-gap];
					motifbuff[6+gap]=NULL;
					for(int i=0;i<left;i++) buff[i]=sequence[posi-left+i];
					for(int i=0;i<gap;i++) buff[i+left]=sequence[posi+3+i];
					for(int i=0;i<right;i++) buff[i+left+gap]=sequence[j+3+i];
					buff[left+right+gap]=NULL;
					WriteSeq(head,motifbuff,buff,sign);
				}
			}
		  }//end for r
	}

}
void GetFontseq(Result *head)
{
	Bseqstru *Bcur=BseqHead;
	Seqstru *Fcur=seqHead;
	Result *cur=head;
	while(cur!=NULL)
	{
		memset(cur->buffmatrix,0,sizeof(cur->buffmatrix));
		countMatrix(cur->buffmatrix,cur->motif,0.001,0);
		cur=cur->next;
	}
	while(Bcur!=NULL)
	{
		  int lens=int(Bcur->len);
		  char motif[50];
		  for(int i=0;i<lens-34;i++)
		  {
			  int j=0;
			  for(j=0;j<34;j++)
			  {
				  motif[j]=Bcur->seq[i+j];
			  }
			  motif[j]=NULL;
			  GetvalueFont(motif,head,0);
          }
		  Bcur=Bcur->next;
	}
	while(Fcur!=NULL)
	{
		  int lens=int(Fcur->len);
		  char motif[50];
		  for(int i=3;i<lens-37;i++)
		  {
			  int j=0;
			  for(j=0;j<34;j++)
			  {
				  motif[j]=Fcur->seq[i+j];
			  }
			  motif[j]=NULL;
			  GetvalueFont(motif,head,1);
          }
		  Fcur=Fcur->next;
	}

}
void filledMatrix(double matrixbuff[4][50],double matrix[4][50],int left,int right,int len)
{
	for(int i=0;i<left;i++)
	{
		for(int j=0;j<4;j++)
			matrixbuff[j][i]=bg[j];
	}
	for(int i=0;i<len;i++)
	{
		for(int j=0;j<4;j++)
			matrixbuff[j][i+left]=matrix[j][i];
	}
	for(int i=0;i<right;i++)
		for(int j=0;j<4;j++)
			matrixbuff[j][i+len+left]=bg[j];
}
void getOut(Result *head)
{
	if(head!=NULL)
	{
		Result *cur=head->next,*curbuff=cur;
		int start,end,lensbuff;
		while(cur!=NULL)
		{
		    double matrixbuff[3][4][50];
			start=6-cur->len;
			end=head->len-6;
			Tranmatrix(cur->matrix[0],cur->len,matrixbuff);
	        double matrix[4][50];
			memcpy(matrix,matrixbuff[2],sizeof(matrixbuff[2]));
			memset(matrixbuff,0,sizeof(matrixbuff));
			double rbuff=0,temp;
			for(int i=start;i<=end;i++)
			 {
				 int left[3],right[3];
				 if(i<0)
				 {
					 left[0]=0-i;left[1]=left[2]=0;
				 }
				 else
				 {
					 left[0]=0;left[1]=left[2]=i;
				 }
				  int gap=head->len-i-cur->len;
				 if(gap<0)
				 {
					 right[0]=0-gap;right[1]=right[2]=0;
				 }
				 else
				 {
					 right[0]=0;right[1]=right[2]=gap;
				 }
				 filledMatrix(matrixbuff[0],head->matrix[0],left[0],right[0],head->len);
                 filledMatrix(matrixbuff[1],cur->matrix[0],left[1],right[1],cur->len);
				 filledMatrix(matrixbuff[2],matrix,left[2],right[2],cur->len);
				 lensbuff=head->len+left[0]+right[0];
				 int sign=0;
			     temp=Cor(0,matrixbuff,lensbuff,sign);
			     if(rbuff<temp) rbuff=temp;
			 }//end for
			if(rbuff>CO)
			{
				if(cur==head->next) {head->next=cur->next;cur=head->next;}
				else
				{
					curbuff->next=cur->next;
					cur=curbuff->next;
				}//end if
			}
			else
			{
				curbuff=cur;
				cur=cur->next;
			}//end if
		}//end while
		getOut(head->next);
	}
}
void resort(Result *head)
{
	    if(head!=NULL)
		{
			Result *cur=head->next;
			while(cur!=NULL)
			{
				if(cur->Mscore>head->Mscore)
				{
					double Mscore=head->Mscore,segment=head->segment,matrix[2][4][50];
					memcpy(matrix,head->matrix,sizeof(head->matrix));
					int len=head->len;
					char motif[50];
					strcpy(motif,head->motif);
					head->Mscore=cur->Mscore;cur->Mscore=Mscore;head->segment=cur->segment;cur->segment=segment;
					head->len=cur->len;cur->len=len;memcpy(head->matrix,cur->matrix,sizeof(cur->matrix));
					memcpy(cur->matrix,matrix,sizeof(matrix));
					strcpy(head->motif,cur->motif);strcpy(cur->motif,motif);
				}
				cur=cur->next;
			}
			resort(head->next);
		}
}
Result *reRedun(Result *head)
{

	  Result *cur=head,*curbuff=cur;
	  while(cur!=NULL)
	  {
		  if(cur->len<5||cur->len>25)
		  {
			  if(cur==head){ head=head->next;cur=head;}
			  else
			  {
				  curbuff->next=cur->next;
				  cur=curbuff->next;
			  }
		  }
		  else
		  {
			  curbuff=cur;
			  cur=cur->next;
		  }
	   }
      return head;
}
inline void getBpbuff(char bp,char *buff)
{
	   switch(bp)
		{
		case 'A':
			strcpy(buff,"AMRWDHVN");break;//DHV
		case 'G':
			strcpy(buff,"GRSKBDVN");break;//BDV
		case 'C':
			strcpy(buff,"CMSYBHVN");break;//BHV
		case 'T':
			strcpy(buff,"TWYKBDHN");break;//BDH
		default:
			break;
		}
}
void Getseq(char *seq,int count,int k,tree *root,Gsize NC,Result *cur,int sign)
{
	if(count<k)
	{
		int len=8,signal=sign;
		if(sign==4) len=1;
		int posi=count;
		char buff[10];
		size_t lens=strlen(cur->motif);
		if(count>2) posi=lens-3+count-3;
		getBpbuff(cur->motif[posi],buff);
		for(int i=0;i<len;i++)
		{
			if(i>0) signal=sign+1;
			seq[count]=buff[i];
			Getseq(seq,count+1,k,root,NC,cur,signal);
		}
	}
	else
	{
	   seq[count]=NULL;
	   char seqbuff[50];
	   char motif[50];
	   size_t lens=strlen(cur->motif);
	   for(int i=0;i<3;i++) motif[i]=seq[i];
	   for(size_t i=3;i<(3+lens-6);i++) motif[i]='N';
	   for(int i=lens-3;i<lens;i++) motif[i]=seq[i-lens+3+3];
	   motif[lens]=NULL;
	   int number[2];
	   memset(number,0,sizeof(number));
	   getZscore(motif,root,number,0,seqbuff,0);
	   double E=NC.Fsize/(double)NC.Bsize*(double)number[1],F;
	   double p=(double)number[1]/(double)NC.Bsize;
	   if(E!=0)
		{
			 double score=(double(number[0])-E)/sqrt(E);
			 F=double(number[0])/E;
			 if(score>cur->Zscore&&F>cur->Fold)
			 {
				cur->Zscore=score;
				cur->F=F;
				strcpy(cur->motif,motif);
			}
		}//end if
	}
}
void updataCore(Result *choose,tree *root,Gsize NC)
{
	Result *cur=choose;
	char seq[50];
	cur->Fold=cur->F;
    Getseq(seq,0,6,root,NC,cur,0);

}
void seqTranto(struct Seqstru *Head,string seq)
{
	size_t lens=seq.size();
	Head->len=lens+40;
	Head->seq=new char[lens+42];
	Head->seqlog=new double[lens+42];
	Head->transeqlog=new double[lens+42];
	char *seqbuff=new char[lens+42],*Transeq=new char [lens+42];
	for(size_t j=0;j<10;j++) seqbuff[j]='N';
	for(size_t j=10;j<(lens+10);j++)
	{
		seqbuff[j]=seq[j-10];
	}
	for(size_t j=10+lens;j<(lens+40);j++) seqbuff[j]='N';
	seqbuff[lens+40]=NULL;
	lens=lens+40;
	Translate(seqbuff,Transeq);
	getBglog(seqbuff,Head->seqlog);
	getBglog(Transeq,Head->transeqlog);
	strcpy(Head->seq,seqbuff);
	delete [] seqbuff;
	delete [] Transeq;
}

void scanSe(string seq,tree *root,int sign)
{
  char motif[40];
  int lens=seq.size();
  for(int i=0;i<lens;i++)
  {
      if(!isupper(seq[i])) seq[i]=toupper(seq[i]);
      if(seq[i]!='A'&&seq[i]!='G'&&seq[i]!='C'&&seq[i]!='T'&&seq[i]!='N') seq[i]='N';
  }
  if(sign==1)
  {
	  GetvalueAll(seq,bglog);
	  Bseqsdel=Bseqscur;
	  Bseqscur->next=new Bseqstru;
	  Bseqscur=Bseqscur->next;
  }
  if(sign==0)
  {
	  seqTranto(seqscur,seq);
	  seqsdel=seqscur;
	  seqscur->next=new Seqstru;
	  seqscur=seqscur->next;
  }
  char buff[30];
  for(int i=0;i<20;i++)  buff[i]='N';
  buff[20]=NULL;
  seq+=buff;
  for(int i=0;i<lens;i++)
  {
	  int j=0;
	  for(j=0;j<20;j++)
	  {
		  motif[j]=seq[i+j];
	  }
      motif[j]=NULL;
	  Getvalue(motif,root,sign);
  }
 }


void sort(struct Result *Top,int count,int maxpos)
{
		if(count<maxpos)
		{
			for(int i=count;i<maxpos;i++)
			{
				struct Result temp=Top[count];
				if(Top[count].Zscore<Top[i].Zscore) {Top[count]=Top[i];Top[i]=temp;}
			}
			sort(Top,count+1,maxpos);
		}
}
void travalTree(tree *root,int len,char* seq,int sign,Gsize NC)
{
     int signal=0;
     if(sign!=0)
     {
         if(sign!=1)
          {
             for(int i=0;i<len;i++)
               seq[sign-2-i]='N';
           }
          seq[sign-1]=root->Bp;


     }
     for(int i=0;i<(15-len);i++)
     {

             if(root->child1[i]!=NULL&&root->child2[i]!=NULL&&root->child3[i]!=NULL&&root->child4[i]!=NULL)
             {
                travalTree(root->child1[i],i,seq,sign+1+i,NC);
                travalTree(root->child2[i],i,seq,sign+1+i,NC);
                travalTree(root->child3[i],i,seq,sign+1+i,NC);
                travalTree(root->child4[i],i,seq,sign+1+i,NC);
             }
             else
             {
                 signal+=1;
                 }
     }
	 if(signal==15&&root->sign==0)
     {
		        root->sign=1;
                seq[sign]=NULL;
                strcpy(choose[pos].motif,seq);
				char Trans[50];
				Translate(seq,Trans);
				size_t lens=strlen(seq);
				tree *cur=getroot(Trans,curroot);
				if(cur->sign==0&&!MICRO)
				{
				choose[pos].countF=root->countF+cur->countF;
				choose[pos].countN=root->countB+cur->countB;
				cur->sign=1;
				}
				else
				{
					choose[pos].countF=root->countF;
				    choose[pos].countN=root->countB;
				}
                double p=choose[pos].countN/(double)NC.Bsize;
        		double E=NC.Fsize/(double)NC.Bsize*(double)choose[pos].countN,F;
        		if(E!=0)
        		{
					//choose[pos].Tscore=choose[pos].countF/(double)NC.Fsize*log(1/p);
        			  F=choose[pos].countF/E;
        			  choose[pos].F=F;
					  choose[pos].Zscore=(double(choose[pos].countF)-E)/sqrt(E);
        			  if(choose[pos].Zscore>0&&F>Foldcut)
        			  {
        				 pos++;
        			  }//end if
        		}
        		if(E!=0&&pos==Topcut+1)
        		{
					if(sortsign==0){sort(choose,0,pos);sortsign=1;}
					else if(choose[pos-1].Zscore>choose[pos-2].Zscore)
        		      sort(choose,0,pos);
        	  	    pos--;
        		}//end if
     }
}
void TrantoPoint(Result *choose,int size)
{
	if(pos<size)
	{
		for(int i=0;i<pos-1;i++)
		{
			choose[i].next=&choose[i+1];
		}
		choose[pos-1].next=NULL;
	}
	else{
		for(int i=0;i<size-1;i++)
		{
			choose[i].next=&choose[i+1];
		}
		choose[size-1].next=NULL;
	}
}
void getMscore(Result *cur)
{
	double bgl=0,fgl=0;
	Aseq *curbuff=cur->head;
	int lens=cur->end-cur->start+1,start=cur->start,end=cur->end;
	while(curbuff!=NULL)
	{
		for(int i=start;i<=end;i++)
          	bgl+=curbuff->bglo[i];
		curbuff=curbuff->next;
	}
	cur->bg=bgl;
	bgl=bgl/cur->segment;
	for(int i=0;i<4;i++)
	{
		for(int j=0;j<lens;j++)
			fgl+=cur->matrix[0][i][j]*log(cur->matrix[0][i][j]);
	}
	cur->Mscore=(cur->segment)/lens*(fgl-bgl);
}
void updataMatrix(Result *cur)
{
	Aseq *curbuff=cur->head,*curbuff1;
	double matrix[2][4][50];
	memset(matrix,0,sizeof(matrix));
	memcpy(matrix[1],cur->matrix[1],sizeof(cur->matrix[1]));
	int  lens=cur->end-cur->start+1,start=cur->start,end=cur->end;
	while(curbuff!=NULL)
	{
		if(cur->segment<=1) break;
		memset(matrix[0],0,sizeof(matrix[0]));
		countMatrix(matrix[0],curbuff->seq,0.001,0);
		double fg=0,bgbuff=0,bgl;
		for(int i=start;i<=end;i++)
		{
			for(int j=0;j<4;j++)
				matrix[1][j][i]=matrix[1][j][i]-matrix[0][j][i];
		}
		memset(matrix[0],0,sizeof(matrix[0]));
        for(int i=start;i<=end;i++)
		{
			for(int j=0;j<4;j++)
				matrix[0][j][i-start]=matrix[1][j][i];
		}
		normalMatrix(matrix[0],0.001,lens);
        for(int i=start;i<=end;i++)
			   bgbuff+=curbuff->bglo[i];
        for(int i=0;i<4;i++)
	    {
			for(int j=0;j<lens;j++)
				fg+=matrix[0][i][j]*log(matrix[0][i][j]);
	    }
		bgl=(cur->bg-bgbuff)/(cur->segment-1);
		double scorebuff=log(cur->segment-1)/lens*(fg-bgl);
		if(scorebuff>cur->Mscore)
		{
			cur->segment+=-1;
			cur->Mscore=scorebuff;
			cur->bg=cur->bg-bgbuff;
			memcpy(cur->matrix[1],matrix[1],sizeof(matrix[1]));
			if(curbuff==cur->head){ cur->head=curbuff->next;curbuff1=cur->head;}
			else curbuff1->next=curbuff->next;
			Aseq *del=curbuff;
			curbuff=curbuff->next;
			delete del;del=NULL;
		}
		else
		{
            memcpy(matrix[1],cur->matrix[1],sizeof(cur->matrix[1]));
			curbuff1=curbuff;
		    curbuff=curbuff->next;
		}
	}
        memset(cur->matrix[0],0,sizeof(cur->matrix[0]));
        for(int i=start;i<=end;i++)
		{
			for(int j=0;j<4;j++)
				cur->matrix[0][j][i-start]=cur->matrix[1][j][i];
		}
	   normalMatrix(cur->matrix[0],0.001,lens);
}
void GetLong(Result *choose)
{
   Result *cur=choose;
   double buff1=0,matrixbuff[4][50];
   for(int i=0;i<4;i++)
	   buff1+=log(bg[i])/log(double(2))*bg[i];
	while(cur!=NULL)
	{
		size_t lens=strlen(cur->motif);
		int start=lens+10,end=-1;
		memset(cur->matrix,0,sizeof(cur->matrix));
		memset(matrixbuff,0,sizeof(matrixbuff));
		double maxscore=0;
		countMatrixE(matrixbuff,cur->motif,maxscore);
		//countMatrix(matrixbuff,cur->motif,0.001,0);
		for(size_t i=5;i<(5+lens);i++)
		{
			for(int j=0;j<4;j++)
			   cur->matrix[1][j][i]=matrixbuff[j][i-5];
		}
		cur->current=new Aseq;
		cur->head=cur->current;
		cur->segment=0;

		Getmatrix(cur,maxscore-1.1);
		if(cur->segment<2){cur->len=0;cur=cur->next;continue;}
		struct Aseq *curbuff=cur->current;
		cur->font->next=NULL;
		cur->current=cur->font;
		delete curbuff;
		curbuff=NULL;
		memcpy(cur->matrix[1],cur->matrix[0],sizeof(cur->matrix[0]));
	  ////////////////////////////////////////////////////////
		size_t lensbuff=lens+10;
		cur->start=0;
		cur->end=lens+9;
		normalMatrix(cur->matrix[0],0.001,lens+10);


		getMscore(cur);
		cur->Mscore=cur->Mscore/cur->segment*log(cur->segment);
        updataMatrix(cur);

		memset(matrixbuff,0,sizeof(matrixbuff));
		for(size_t i=0;i<lensbuff;i++)
		{
		   double buff=0,temp;
		   int sign=0;
		   for(int j=0;j<4;j++)
		   {
			 temp=log(cur->matrix[0][j][i])/log(double(2))*cur->matrix[0][j][i];
			 if((log(cur->matrix[0][j][i]/bg[j])/log(double(2)))>0) sign++;
			 buff+=temp;
		   }
		   if((buff-buff1)>0.25&&sign>0&&sign<3) {start=i;break;}
		}//end for

		for(int i=lensbuff-1;i>=0;i--)
		{
           double buff=0,temp;
		   int sign=0;
		   for(int j=0;j<4;j++)
		   {
			 temp=log(cur->matrix[0][j][i])/log(double(2))*cur->matrix[0][j][i];
			 if((log(cur->matrix[0][j][i]/bg[j])/log(double(2)))>0) sign++;
			 buff+=temp;
		   }
		   if((buff-buff1)>0.25&&sign>0&&sign<3) {end=i;break;}
		}
		if(start==lensbuff||end==-1) {cur->len=0;cur=cur->next;continue;}
		for(size_t i=start;i<=end;i++)
		   for(int j=0;j<4;j++)
			   matrixbuff[j][i-start]=cur->matrix[0][j][i];

		memset(cur->matrix[0],0,sizeof(cur->matrix[0]));
		memcpy(cur->matrix[0],matrixbuff,sizeof(matrixbuff));
		cur->start=start;
		cur->end=end;
		getMscore(cur);
		cur->Mscore=cur->Mscore/cur->segment*log(cur->segment)*log(double(end-start+1));
		while(cur->head!=NULL)
		{
			curbuff=cur->head;
			cur->head=cur->head->next;
			delete curbuff;
			curbuff=NULL;
		}
		cur->current=NULL;
		cur->len=end+1-start;
		cur=cur->next;
	}

}
inline void sortbuff(double buffm[5],int sign,int max)
{
	if(sign<max)
	{
		double temp=buffm[sign];
		for(int i=sign;i<max;i++)
		{
			if(buffm[sign]<buffm[i]){ buffm[sign]=buffm[i];buffm[i]=temp;}
		}
		sortbuff(buffm,sign+1,max);
	}

}
inline void getMbglo(double *slog,int i,int k,double *bglo)
{
	for(int j=i;j<(k+i);j++)
		bglo[j-i]=slog[j];
}
inline void Getmatrix(Result *cur,double score)
{
	Seqstru *curbuff=seqHead;
	size_t k=strlen(cur->motif);
    double buffmatrix[4][50];
	while(curbuff!=NULL)
	{
         char motifbuff[50],Transbuff[100];
		 for(int i=0;i<int(curbuff->len-k-9);i++)
		  {
			  int j=0;
			  for(j=0;j<(k+10);j++)
			  {
				 motifbuff[j]=curbuff->seq[i+j];
			  }
			  motifbuff[j]=NULL;
			  if(strlen(motifbuff)!=k+10) continue;
			  Translate(motifbuff,Transbuff);
			  double scorebuff1=getscore(cur->matrix[1],motifbuff),scorebuff2=getscore(cur->matrix[1],Transbuff);
			  memset(buffmatrix,0,sizeof(buffmatrix));
			  if(scorebuff2==scorebuff1) scorebuff2=-1;
			   if(scorebuff2>=score&&scorebuff2>scorebuff1&&!MICRO)
			   {
				   strcpy(cur->current->seq,Transbuff);
				   getMbglo(curbuff->transeqlog,curbuff->len-i-k-10,k+10,cur->current->bglo);
				   cur->font=cur->current;
				   cur->current->next=new Aseq;
				   cur->current=cur->current->next;
					  countMatrix(buffmatrix,Transbuff,0.001,0);
					   cur->segment++;
					   for(int m=0;m<(10+k);m++)
						  for(int j=0;j<4;j++)
							 cur->matrix[0][j][m]+=buffmatrix[j][m];

			    }//end if
				if(scorebuff1>=score&&scorebuff1>scorebuff2)
				{
					strcpy(cur->current->seq,motifbuff);
					getMbglo(curbuff->seqlog,i,k+10,cur->current->bglo);
					cur->font=cur->current;
					cur->current->next=new Aseq;
				    cur->current=cur->current->next;
						countMatrix(buffmatrix,motifbuff,0.001,0);
						cur->segment++;
						for(int m=0;m<(10+k);m++)
							for(int j=0;j<4;j++)
								cur->matrix[0][j][m]+=buffmatrix[j][m];
				 }//end if

		 }//end for
		 curbuff=curbuff->next;
	}//end while

}
void getBglog(char* seq,double *score)
{
	size_t len=strlen(seq);
	for(size_t j=0;j<len;j++)
	{
		  int n=6,m=6,k=6,l=6;
		switch(markov)
		{
		case 0:
			TranNumber(seq[j],n);
			if(n>3) score[j]=0;
			else score[j]=bglog.Obp[n];
			break;
		case 1:
			if(j<10) score[j]=0;
			else if(j==10)
			{
			   TranNumber(seq[j],n);
			   if(n>3) score[j]=0;
			   else score[j]=bglog.Obp[n];
			}
			else
			{
				TranNumber(seq[j-1],n);
				TranNumber(seq[j],m);
				if(n>3||m>3) score[j]=0;
				else score[j]=bglog.Tbp[n][m];
			}
			break;
		case 2:
			if(j<10) score[j]=0;
			else if(j==10)
			{
			   TranNumber(seq[j],n);
			   if(n>3) score[j]=0;
			   else score[j]=bglog.Obp[n];
			}
			else if(j==11)
			{
				TranNumber(seq[j-1],n);
				TranNumber(seq[j],m);
				if(n>3||m>3) score[j]=0;
				else score[j]=bglog.Tbp[n][m];
			}
			else
			{
				TranNumber(seq[j-2],n);
                TranNumber(seq[j-1],m);
				TranNumber(seq[j],k);
				if(n>3||m>3||k>3) score[j]=0;
				else score[j]=bglog.Thbp[n][m][k];
			}
			break;
		case 3:
			if(j<10) score[j]=0;
			else if(j==10)
			{
			   TranNumber(seq[j],n);
			   if(n>3) score[j]=0;
			   else score[j]=bglog.Obp[n];
			}
			else if(j==11)
			{
				TranNumber(seq[j-1],n);
				TranNumber(seq[j],m);
				if(n>3||m>3) score[j]=0;
				else score[j]=bglog.Tbp[n][m];
			}
			else if(j==12)
			{
				TranNumber(seq[j-2],n);
                TranNumber(seq[j-1],m);
				TranNumber(seq[j],k);
				if(n>3||m>3||k>3) score[j]=0;
				else score[j]=bglog.Thbp[n][m][k];
			}
			else
			{
				TranNumber(seq[j-3],n);
                TranNumber(seq[j-2],m);
                TranNumber(seq[j-1],k);
				TranNumber(seq[j],l);
				if(n>3||m>3||k>3||l>3) score[j]=0;
				else score[j]=bglog.Fbp[n][m][k][l];
			}
			break;
		default:
			break;
		}//end switch
	}//end for
}

void Redundancy(Result *head)
{
	if(head!=NULL)
	{
		Result *cur=head,*curnext=cur->next;
	  while(curnext!=NULL)
	  {
		  Result *curbuff=curnext->next;
		  char motif[50];
		  Translate(curnext->motif,motif);
		  if(strcmp(head->motif,motif)==0||strcmp(head->motif,curnext->motif)==0)
		  {
			  cur->next=curbuff;
			  curnext=curbuff;
		  }
		  else{
		  cur=curnext;
		  curnext=curbuff;
		  }
	   }
	  Redundancy(head->next);
	}
}
void CalBglog(struct Bgpro &bglog,int &markov)
{
   int sum1=0;
   for(int i=0;i<4;i++)
	   sum1+=bglog.Obp[i];
   if(sum1>MARK3) markov=3;
   else if(sum1>MARK2) markov=2;
   else if(sum1>MARK1) markov=1;
   else markov=0;
   for(int i=0;i<4;i++)
   {
	   bglog.Obp[i]=log((bglog.Obp[i]+SMALL)/sum1);
	   if(markov==0) continue;
	   int sum2=0;
	   for(int j=0;j<4;j++) sum2+=bglog.Tbp[i][j];
	   for(int j=0;j<4;j++)
	   {
		   bglog.Tbp[i][j]=log((bglog.Tbp[i][j]+SMALL)/sum2);
		   if(markov==1) continue;
		   int sum3=0;
		   for( int k=0;k<4;k++) sum3+=bglog.Thbp[i][j][k];
		   for(int k=0;k<4;k++)
		   {
			   bglog.Thbp[i][j][k]=log((bglog.Thbp[i][j][k]+SMALL)/sum3);
			   if(markov==2) continue;
			   int sum4=0;
			   for(int m=0;m<4;m++) sum4+=bglog.Fbp[i][j][k][m];
			   for(int m=0;m<4;m++)
			   {
				   bglog.Fbp[i][j][k][m]=log((bglog.Fbp[i][j][k][m]+SMALL)/sum4);
			   }//end for m
		   }//end for k
	   }//end for j
   }//end for i

}
