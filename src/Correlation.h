#include "math.h"
#include <string>
#include <iostream>
#include <fstream>

using namespace std;
struct Aseq
{
	char seq[35];
	double bglo[35];
	Aseq *next;
};

struct Result
{ 
	int countF;
	int countN;
	char motif[50];
    double Zscore;
	double Mscore;
	double Tscore;
	string Fseq;
	string Bseq;
	struct Aseq *head;
	struct Aseq *current;
	struct Aseq *font;
	int start;
	int end;
	double bg;
	double segment;
	double matrix[2][4][50];
	double buffmatrix[4][50];
	int len;
	double F;
	double Fold;
	int site;
	struct Result *next;
};
extern const char *filename2,*filename1;
extern double bg[4];
double Correlation(char *motif1,char *motif2,double beta,int &site);
void merge(struct Result* head,int start,int end,double beta,double matrix[4][50]);
void Translate(char *motif,char *Trans);
void Getconsensus(double matrix[4][50],char *consensus,int lens);
double Cor(int sign,double matrix[3][4][50],size_t lens,int& signal);
void normalMatrix(double matrix[4][50],double beta,int length);
void countMatrix(double matrix[4][50],char *motif,double beta,int sign);
void Tranmatrix(double matrix[4][50],int length,double transmatrix[3][4][50]);
inline double getscore(double matrix[4][50],char *motif);
void countMatrixE(double matrix[4][50],char *motif,double &maxscore);
