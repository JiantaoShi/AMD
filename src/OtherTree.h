#include <string>
#include <iostream>
using namespace std;
struct tree
{
	char  Bp;
	int countF;
	int countB;
	int sign;
	tree* child1[15];
	tree* child2[15];
	tree* child3[15];
	tree* child4[15];  
};
struct Bgpro
{
	double Obp[4];
	double Tbp[4][4];
	double Thbp[4][4][4];
	double Fbp[4][4][4][4];
};
struct Bseqstru
{
	char *seq;
	size_t len;
	Bseqstru *next;
};
struct Gsize
{
	long Fsize;
	long Bsize;
};//bp size
extern tree *curroot;
extern int MICRO;

extern struct Gsize NC;
extern struct Bseqstru *BseqHead,*Bseqscur,*Bseqsdel;
void constructTree(tree* root,int len,int count,int sign);
tree * getroot(char *motif,tree *root);
void TranNumber(char bp,int& i);
void GetvalueAll(string sequence,struct Bgpro &bglog);
void Getvalue(char *sequence,struct tree *root,int sign);
void DistroyTree(tree *root,int len);
void getZscore(char* seq,tree *root,int number[2],int count,char *seqbuff,int sign);
inline tree * getadd(char Bp,int pos,tree *root);


