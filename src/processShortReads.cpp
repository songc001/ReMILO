#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<unistd.h>
#include<iostream>
#include<fstream>
#include<vector>
using namespace std;

#define  percent  0.20
#define   maxBuflen  1000000
#define   MAPRATE  0.4
#define   PERRATE  0.6 
#define   MaxNum   1000000
typedef  struct  Read{
         int  pos;
         char name[100];
         int  direction;
         float maprate;
         float perfectrate;
         int length;
}Read;
typedef struct  SubContig{
    string name;
    vector<Read>   readlist1;
    vector<Read>   readlist2;
    int    length;
    int    flag;
    string seq;
    int    pos;
}SubContig;
typedef  struct  MyContig{
        //   char  cigar[100000];
           vector<Read>   read1;
           vector<Read>   read2;
           int   length;
           char  name[100];
           vector<SubContig>  sublist;
}Contig;
typedef   struct  Miscontig{
          Contig  contig1;
          Contig  contig2;
          Contig  thirdcontig1;
          Contig  thirdcontig2;
          int     support[3];
}Miscontig;
typedef   struct  Index{
          char    str[1000];
          int     num;
          }Index;

//char  qname[1000000],q 

//FILE  *fp,*fq,*fq1,*fp1;
/*
Contig  g_contig1,g_contig2;
//Mispos   mispos[10];
int convertInt(char str[])
{
   int length;
   int sum=0;

   length=strlen(str);
   for(int i=0; i<length; i++)
   {
      sum+=(str[i]-'0')*pow(10,(length-i-1));
   }
   return sum;
} 
*/
/*
   1   represent  reverse  
   0   represent  positive 
*/

typedef struct  haContig{
    string   name;
    vector<Read>    readlist;
}haContig;
typedef struct  ArcNode{
   int  adjvex;
   struct ArcNode *nextarc;
}ArcNode;
typedef struct  VNode{
   string str;
   int  conpos;
   int  chrpos;
   int  coverage;
   int color;
   ArcNode *firstarc;
}VNode;
typedef struct
{
   VNode  vertex[MaxNum];
   int   vexnum,arcnum;
}DeGraph;

vector<haContig>  haContiglist;
vector<Contig>    contiglist;
double    oveConCov,oveChrCov;

vector<string>  split(string str,const char *delim)
{
    char *p;
    char buffer[maxBuflen];
    vector<string>  data;
     
    strcpy(buffer,str.c_str());
    p=strtok(buffer,delim);
    while(p)
    {
        data.push_back(p);
        p=strtok(NULL,delim);
    }
    return data;
}
void  calOveConCov()
{
    int len=0;
    int coverage=0;
    for(int i=0;i<contiglist.size();i++)
    {
        for(int j=0;j<contiglist[i].sublist.size();j++)
        {
           for(int k=0;k<contiglist[i].sublist[j].readlist1.size();k++)
           {
               len+=contiglist[i].sublist[j].readlist1[k].length;
           }
           for(int k=0;k<contiglist[i].sublist[j].readlist2.size();k++)
           {
               len+=contiglist[i].sublist[j].readlist2[k].length;
           }
        }
    }

   
}
void  calChrCov(int connum,int subnum)
{   
    int  tolen=0; 
    for(int i=0;i<contiglist[connum].sublist[subnum].readlist1.size();i++)
    {
        tolen+=contiglist[connum].sublist[subnum].readlist1[i].length;
    }
}
int   calCov(int connum,int subnum)
{
    int tolen=0,coverage=0;
    for(int i=0;i<contiglist[connum].sublist[subnum].readlist1.size();i++)
    {
        tolen+=contiglist[connum].sublist[subnum].readlist1[i].length;
    
    }
    for(int i=0;i<contiglist[connum].sublist[subnum].readlist2.size();i++)
    {
        tolen+=contiglist[connum].sublist[subnum].readlist2[i].length;
    }
    coverage=tolen/(contiglist[connum].sublist[subnum].length);

    return  coverage; 
}

void  getCfile(ifstream  &c)
{
    string  buffer;
    Contig  contig;
    SubContig  subcontig;
    haContig   hacontig;
    vector<string>  data;
    int  order;
    
    subcontig.length=0;
    while(c.good())
    {
       getline(c,buffer);
       if(buffer[0]=='>')
       {
          if(subcontig.length!=0)
          {
             contig.sublist.push_back(subcontig);
             subcontig.length=0;
             subcontig.seq="";
          }
          buffer.erase(0,1);
          data=split(buffer,"_");
          subcontig.name=buffer;
          hacontig.name=buffer;
          haContiglist.push_back(hacontig);
          if(data.size()==3)
          {
             subcontig.flag=1;
          }
          else
          {
             subcontig.flag=0;
          }
          if(contig.name!=data[0])
          {
             if(contig.name[0]=='c')  
             {   
                contiglist.push_back(contig);
                contig.sublist.clear();
             }
             strcpy(contig.name,data[0].c_str());
          }
       }
       else
       {
          subcontig.length+=buffer.size();
          subcontig.seq+=buffer;
       }
          
    }
    contig.sublist.push_back(subcontig);
    contiglist.push_back(contig);    
}
       
void initGraph(DeGraph &g,int kmer,int connum,int subnum)
{
    int tmp,tmp1;
    string str;
    int count=0;
    g.vexnum=0;
    g.arcnum=0;
    for(int i=0;i<4;i++)
    {
       tmp=g.vexnum;
       g.vexnum+=contiglist[connum].sublist[subnum+i].length-kmer+1;
       for(int j=tmp;j<g.vexnum;j++)
       {
           g.vertex[j].str="";
           for(int k=count;k<count+kmer;k++)
           {
              g.vertex[j].str+=contiglist[connum].sublist[subnum].seq[k];
           }
           g.vertex[j].conpos=count+1;
           g.vertex[j].chrpos=count+1+contiglist[connum].sublist[subnum+i].pos;
           g.vertex[j].coverage=calCov(connum,subnum);
           count++;
       }
       tmp1=g.arcnum;
       g.arcnum+=contiglist[connum].sublist[subnum].length-kmer;
       for(int j=tmp1;j<g.arcnum;j++)
       {
           g.vertex[j].firstarc=NULL;
       }
    }
}    
void createGraph(DeGraph  &g,int kmer,int connum,int subnum)
{
    int i,s,e,w;
    ArcNode *p;
    initGraph(g,kmer,connum,subnum);
    for (i=1;i<g.arcnum;i++)
    {
      p=(ArcNode*)malloc(sizeof(ArcNode));
      p->adjvex=i;
      p->nextarc=g.vertex[i+1].firstarc;
      g.vertex[i+1].firstarc=p;
    }
}
void color(DeGraph  &g,int kind)
{
    int count=0,count1=0,color=0;
    for(int i=0;i<g.vexnum;i++)
    {
       if(kind==1&&(g.vertex[i].coverage<oveConCov-10||g.vertex[i].coverage>oveConCov+10))
       {
          g.vertex[i].color=1;
       }
       else
       {
          if(kind==0&&(g.vertex[i].coverage<oveConCov-4||g.vertex[i].coverage>oveConCov+4))
          {
            g.vertex[i].color=1;
          }
          else
          {
            g.vertex[i].color=0;
          }
        }
     }
     for(int i=0;i<g.vexnum;i++)
     {
        if(g.vertex[i].color==1)
        {
           count++;
        }
        else
        {
           count1++;
        }
      }
      if(count1>count)
      {
         color=1;
      }
      else
      {
        color=0;
      }
      for(int i=0;i<g.vexnum;i++)
      {
        g.vertex[i].color=color;
      }
}
void detectError(DeGraph &g,int  vexnum[4])
{
   int flag[4];
   int iferror=0;
   flag[0]=0;
   for(int i=0;i<vexnum[0];i++)
   {
      if(g.vertex[i].color==1)
      {
         flag[0]=1;
      }
   }
   for(int i=vexnum[0];i<vexnum[1];i++)
   {
      if(g.vertex[i].color==1)
      {
         flag[1]=1;
      }
   }
   
   for(int i=vexnum[1];i<vexnum[2];i++)
   {
      if(g.vertex[i].color==1)
      {
         flag[2]=1;
      }
   }
   
   for(int i=vexnum[2];i<vexnum[3];i++)
   {
      if(g.vertex[i].color==1)
      {
         flag[3]=1;
      }
   }
    if((flag[0]==1&&flag[3]==1)||(flag[1]==1&&flag[2]==1))
    {
        iferror=1;
    }
   
}
void calMapRate(string cigar,Read  &read)
{
    int  i;
    int mc=0,sc=0,hc=0;
    int length=0;
    //vector<int>   num;
    string  num;
    int   figure=0;
    
   // strcpy(cigar,strcigar.c_str());
    while(cigar[i]!='\0')
    {
       if(cigar[i]>=48&&cigar[i]<=57)
       {
           num+=cigar[i];
       } 
       else
       {
           length++;
           figure=atoi(num.c_str());
           if(cigar[i]=='M')
           {
              mc+=figure;
           }
           if(cigar[i]=='S')
           {
              sc+=figure;
           }
           if(cigar[i]=='H')
           {
              hc+=figure;
           }
           num="";
           figure=0;
       }
       i++;
    }
 //   cout<<mc<<hc<<sc<<endl;
    read.length=length;
    read.maprate=(length-sc-hc)/(length*1.0);
    read.perfectrate=mc/( (length-sc-hc)*1.0 );    
}
        
void  assignRead(vector<string> data,Read &read)
{ 
    strcpy(read.name,data[0].c_str());
    read.pos=atoi(data[3].c_str());
    if( (atoi(data[1].c_str())/16)%2==0 )
    {
        read.direction=0;
    }
    else
    {
        read.direction=1;
    }
    calMapRate(data[5],read);
}

void dealRead(vector<string>  data)
{
    Read  read;
    cout<<data[0]<<endl;
    assignRead(data,read);
    if(read.maprate>MAPRATE  &&  read.perfectrate>PERRATE)
    {
       for(int i=0;i<haContiglist.size();i++)
       {
          if(data[2]==haContiglist[i].name)
          {
             haContiglist[i].readlist.push_back(read);
          }
       }
    }
  //  getchar();

}
    
void  getRfile(ifstream &r1)
{
    string buffer;
    vector<string>  data;
    Read   read;
    getline(r1,buffer);
    while(r1.good()&&buffer[0]=='@')
    {
        getline(r1,buffer);
    }
    while(r1.good())
    {
        data=split(buffer,"	");
        getline(r1,buffer);
        if(data[2]!="*")
        {
           dealRead(data);
        }
        data.clear();
    } 
}
void  linkR1toC()
{    
    for(int i=0;i<haContiglist.size();i++)
    {
       for(int j=0;j<contiglist.size();j++)
       {
          for(int k=0;k<contiglist[j].sublist.size();k++)
          {
             if(contiglist[j].sublist[k].name==haContiglist[i].name)
             {
                contiglist[j].sublist[k].readlist1=haContiglist[i].readlist;
                haContiglist[i].readlist.clear();
                goto  part1;
             }
          }
       }
       part1:   true;//cout<<"sucesss"<<endl;
           
    }          
}
void  linkR2toC()
{    
    for(int i=0;i<haContiglist.size();i++)
    { 
       for(int j=0;j<contiglist.size();j++)
       {
          for(int k=0;k<contiglist[j].sublist.size();k++)
          {
             if(contiglist[j].sublist[k].name==haContiglist[i].name)
             {
                contiglist[j].sublist[k].readlist2=haContiglist[i].readlist;
                haContiglist[i].readlist.clear();
                goto  part1;
             }
          }
       }
       part1:  true;
    }          
}
void   formalOutput(ofstream &outmis)
{
    for(int i=0;i<contiglist.size();i++)
    {
       outmis<<'>'<<contiglist[i].name<<'	'<<endl;
       for(int j=0;j<contiglist[i].sublist.size();j++)
       {
          outmis<<'	'<<contiglist[i].sublist[j].name<<": "<<contiglist[i].sublist[j].readlist1.size()<<'	'<<contiglist[i].sublist[j].readlist2.size()<<endl;
       }
    }
}
int judgeDirect(int direct)
{
    direct=direct/16;
    direct=direct%2;
    if(direct==1)
    {
         return  1;
    }
    else
    {
        return   0;
    }
}/*
int    dealSeq(char seq[],char (&moseq)[10000000])
{
 {
    int i;
    if(contig.direction==0)
    {
     contig.refpos2=contig.refpos1+contig.tolen;
    }
    else
    {
      contig.refpos2=contig.refpos1+contig.tolen;
    }
  //  cout<<contig.refpos1<<" "<<contig.refpos2<<" "<<length<<endl;
    return 0;
 }
  */          
//output  subcontig
 
int main(int argc,char *argv[])
{ 
  FILE *fp,*fp1,*fp2,*fq,*fq1,*fq2;
  char ch;
  char chflag[10];
  char chpos[20];
  char str[1000000];
  char tempnum[1000];
  char qname[100], qname1[100],refname[100];
  char cigar[100000],sequence[10000000];
  int  flag;
  int  pos;
  int  con_length;
  int  k,thousand=0;;
  int  count, mark, connum=0,misnum=0; 
  int tempcount=0; 
  int  precigar[2000],becigar[2000],lencigar[2000];
  vector<char>  chromosome; 
  int chrcount=0,chrnum=0;
  Index   index[3000];
  Contig contig[200000];
  Miscontig  miscontig[20000];
  int        miscount=0,refcount=0,misassemblynum=0;
  char    tempname[100];
  Read    read;
  int direction;
  char  tempstr[100],tempstr1[100];
  if((fp=fopen(argv[1],"r"))==NULL)
  {
     printf("read1.sam open error");
  }
  if((fp1=fopen(argv[2],"r"))==NULL)
  {
     printf("read2.sam  pen error");
  }
  if((fp2=fopen(argv[3],"r"))==NULL)
  {
     printf("contig file open error");
  }
  if((fq=fopen(argv[4],"w"))==NULL)
  {
     printf("fq open error");
  }
  if((fq1=fopen(argv[5],"w"))==NULL)
  {
     printf("fq open error");
  }
  if((fq2=fopen(argv[6],"w"))==NULL)
  {
     printf("fq open error");
  }
/*  if((fq1=fopen("subcontig","w"))==NULL)
  {
     printf("fq1 open error");
  }
*/
//  cout<<"1"<<endl;
  if(argc==1)
  {
     printf("usage:");
     printf("./shortread fq1 fq2 ");
  }
  ch=fgetc(fp2);
    tempstr[0]='\0';
  while(ch!=EOF)
  {
    if(ch=='>')
    {
         ch=fgetc(fp2);
         count=0;
         while(ch!='_')
         {
            tempstr1[count]=ch;
            miscontig[miscount].contig1.name[count++]=ch;
            ch=fgetc(fp2);
         }
         while(ch!='\n')
         {
            miscontig[miscount].contig1.name[count++]=ch;
            ch=fgetc(fp);
         }
         miscontig[miscount].contig1.name[count++]='\0';
         if(strcmp(tempstr,tempstr1)==0)
         {
            strcpy( miscontig[miscount].thirdcontig1.name,miscontig[miscount].contig1.name);
            miscontig[miscount].contig1=miscontig[miscount-1].contig1;
            while(ch!='>')
            {
                if(ch!='\n')    chromosome.push_back(ch);
                ch=fgetc(fp2);
            }
            miscontig[miscount].thirdcontig1.length=chromosome.size();
         }
         else
         {
              while(ch!='>')
              {
                  if(ch!='\n')    chromosome.push_back(ch);
                  ch=fgetc(fp2);
              }
              miscontig[miscount].contig1.length=chromosome.size();
          
              ch=fgetc(fp2);
              count=0;
              while(ch!='\n')
              {
                  miscontig[miscount].thirdcontig1.name[count++]=ch;
                  ch=fgetc(fp2);
              }
              miscontig[miscount].thirdcontig1.name[count++]='\0';
              while(ch!='>')
              {
                if(ch!='\n')    chromosome.push_back(ch);
                ch=fgetc(fp2);
              }
              miscontig[miscount].thirdcontig1.length=chromosome.size();
         }

  //    cout<<miscontig[miscount].thirdcontig1.name<<endl; 
       
         ch=fgetc(fp2);
         count=0;
         while(ch!='\n')
         {
         	 miscontig[miscount].thirdcontig2.name[count++]=ch;
           	ch=fgetc(fp2);
         }
      		miscontig[miscount].thirdcontig2.name[count]='\0';
     	 while(ch!='>')
         {
        	 if(ch!='\n')    chromosome.push_back(ch);
         	ch=fgetc(fp2);
      	 }
         miscontig[miscount].thirdcontig2.length=chromosome.size();
         ch=fgetc(fp2);
   
         count=0;
         while(ch!='\n')
         {
                 misassemblynum=5;
        	 miscontig[miscount].contig2.name[count++]=ch;
     		 ch=fgetc(fp2);
         }
         miscontig[miscount].contig2.name[count]='\0';
      
         while(ch!='>'&&ch!=EOF)
         {
        	 if(ch!='\n')    chromosome.push_back(ch);
        	 ch=fgetc(fp2);
         }
         miscontig[miscount].contig2.length=chromosome.size();
         miscount++;
         strcpy(tempstr,tempstr1);
     }
    
   }
//  cout<<miscount<<endl;
//  cout<<"load "<<miscount<<" misassembled contig"<<endl;
 // getchar();
  ch=fgetc(fp);
  while(ch!=EOF)
  {
//read  read1.sam
      if(ch=='@')
      {
        fgets(str,1000,fp);
        ch=fgetc(fp);
      }
      else
      {
          k=0;    
          while(ch!='	')
          {
            qname[k++]=ch;
            ch=fgetc(fp);
          }
          qname[k]='\0';

          k=0;
          ch=fgetc(fp);
          while(ch!='	')
          {
            chflag[k++]=ch;
            ch=fgetc(fp);
          }
          chflag[k]='\0';
          flag=atoi(chflag);
         // read  refname
          k=0;
          ch=fgetc(fp);
          while(ch!='	')
          {
            refname[k++]=ch;
            ch=fgetc(fp);
          }
          refname[k]='\0';
         //contig's pos on  the reference  genome
          
          k=0;
          ch=fgetc(fp);
          while(ch!='	')
          {
            chpos[k++]=ch;
            ch=fgetc(fp);
          }
          chpos[k]='\0';
          pos=atoi(chpos);

          ch=fgetc(fp);
          while(ch!='	')
          {
            ch=fgetc(fp);
          }
        // cigar        
          k=0;
          ch=fgetc(fp);
          while(ch!='	')
          {
             cigar[k++]=ch;
             ch=fgetc(fp);
          }
          cigar[k]='\0';  
          /*
          save info
*/
          direction=judgeDirect(flag);
//deal  reads
          if(strcmp(refname,"*")!=0)
          {
                refcount++;
                strcpy(read.name,qname);
                read.pos=pos;
                read.direction=direction;
                for(int i=0;i<miscount;i++)
                {
                   if(( strcmp(refname,miscontig[i].contig1.name) )==0
                        &&(  ((miscontig[i].contig1.length-pos)<300) || (pos<300)  ) )
                   {
                      miscontig[i].contig1.read1.push_back(read);
                   } 
       //     cout<<qname<<endl;
                   if(( strcmp(refname,miscontig[i].contig2.name) )==0
                        &&(  ((miscontig[i].contig2.length-pos)<300)|| (pos<300)  ) )
                   {
                      miscontig[i].contig2.read1.push_back(read);
                   } 
                   if(( strcmp(refname,miscontig[i].thirdcontig1.name) )==0
                        &&(  ((miscontig[i].thirdcontig1.length-pos)<300)|| (pos<300)  )  )
                   {
                      miscontig[i].thirdcontig1.read1.push_back(read);
                   } 
                   if(( strcmp(refname,miscontig[i].thirdcontig2.name) )==0
                        &&(  ((miscontig[i].thirdcontig2.length-pos)<300)|| (pos<300)  )  )
                   {
                      miscontig[i].thirdcontig2.read1.push_back(read);
                   } 
                }
          }
        
       //  cout<<qname<<endl; 
          if(strcmp(qname,qname1)!=0)
          {
          	connum++;
          	if(connum==1000000)
          	{
          	   thousand++;
             	   cout<<"deal"<<connum*thousand<<"reads"<<endl;
                //   cout<<miscontig[0].contig1.read1[thousand].pos<<endl;
               //    cout<<miscontig[0].contig1.read1[thousand].name<<endl;
                   connum=0;
                }
           
          }          
          strcpy(qname1,qname);
          fgets(str,1000000,fp);   
          ch=fgetc(fp);
       
      }
  }  
  //cout<<refcount<<endl;
  ch=fgetc(fp1);
//read2.sam
  refcount=0;
  while(ch!=EOF)
  {
//read file
      if(ch=='@')
      {
        fgets(str,1000,fp1);
        ch=fgetc(fp1);
      }
/*
read qname 
*/
      else
      {
          k=0;    
          while(ch!='	')
          {
            qname[k++]=ch;
            ch=fgetc(fp1);
          }
          qname[k]='\0';
/*
read flag  and convert
*/
          k=0;
          ch=fgetc(fp1);
          while(ch!='	')
          {
            chflag[k++]=ch;
            ch=fgetc(fp1);
          }
          chflag[k]='\0';
          flag=atoi(chflag);
/* 
 read  refname
*/ 
          k=0;
          ch=fgetc(fp1);
          while(ch!='	')
          {
            refname[k++]=ch;
            ch=fgetc(fp1);
          }
          refname[k]='\0';
/*
  contig's pos on  the reference  genome
*/          
          k=0;
          ch=fgetc(fp1);
          while(ch!='	')
          {
            chpos[k++]=ch;
            ch=fgetc(fp1);
          }
          chpos[k]='\0';
          pos=atoi(chpos);

          ch=fgetc(fp1);
          while(ch!='	')
          {
            ch=fgetc(fp1);
          }
  //        cout<<"load cigar"<<endl;
/*
 cigar
*/          
          k=0;
          ch=fgetc(fp1);
          while(ch!='	')
          {
             cigar[k++]=ch;
             ch=fgetc(fp1);
          }
          cigar[k]='\0';  
    
//deal read2
          direction=judgeDirect(flag);
          if(strcmp(refname,"*")!=0)
          {
                refcount++;
                strcpy(read.name,qname);
                read.pos=pos;
                read.direction=direction;
                for(int i=0;i<miscount;i++)
                {
                   if(( strcmp(refname,miscontig[i].contig1.name) )==0
                        &&(  ((miscontig[i].contig1.length-pos)<300) || (pos<300)  ) )
                   {
                      miscontig[i].contig1.read2.push_back(read);
                   } 
       //     cout<<qname<<endl;
                   if(( strcmp(refname,miscontig[i].contig2.name) )==0
                        &&(  ((miscontig[i].contig2.length-pos)<300)|| (pos<300)  ) )
                   {
                      miscontig[i].contig2.read2.push_back(read);
                   } 
                   if(( strcmp(refname,miscontig[i].thirdcontig1.name) )==0
                        &&(  ((miscontig[i].thirdcontig1.length-pos)<300)|| (pos<300)  )  )
                   {
                      miscontig[i].thirdcontig1.read2.push_back(read);
                   } 
                   if(( strcmp(refname,miscontig[i].thirdcontig2.name) )==0
                        &&(  ((miscontig[i].thirdcontig2.length-pos)<300)|| (pos<300)  )  )
                   {
                      miscontig[i].thirdcontig2.read2.push_back(read);
                   } 
                }
          }
           //  getchar();
       //  cout<<qname<<endl; 
          if(strcmp(qname,qname1)!=0)
          {
             connum++;
             if(connum==1000000)
             {
               thousand++;
               cout<<"deal"<<connum*thousand<<"reads"<<endl;
               cout<<miscontig[0].contig1.read2[thousand].pos<<endl;
               cout<<miscontig[0].contig1.read2[thousand].name<<endl;
               connum=0;
             }
           
          }
          
        //  cout<<"remilo have dealed "<<thousand*1000000+connum<<" reads"<<endl; 
          strcpy(qname1,qname);
          fgets(str,1000000,fp);   
          ch=fgetc(fp);
       
      }
  }  
  for(int i=0;i<miscount;i++)
  {
      for(int j=0;j<miscontig[i].contig1.read1.size();j++)
      {
        fputs(miscontig[i].contig1.read1[j].name,fq);
        fputc(' ',fq);
        sprintf(chpos,"%d",miscontig[i].contig1.read1[j].pos);
        fputs(chpos,fq);
        fputc('	',fq);
      }
      fputc('\n',fq);
      for(int j=0;j<miscontig[i].contig2.read1.size();j++)
      {
        fputs(miscontig[i].contig2.read1[j].name,fq);
        fputc(' ',fq);
        sprintf(chpos,"%d",miscontig[i].contig2.read1[j].pos);
        fputs(chpos,fq);
        fputc('	',fq);
      }
      fputc('\n',fq);
      for(int j=0;j<miscontig[i].thirdcontig1.read1.size();j++)
      {
        fputs(miscontig[i].thirdcontig1.read1[j].name,fq);
        fputc(' ',fq);
        sprintf(chpos,"%d",miscontig[i].thirdcontig1.read1[j].pos);
        fputs(chpos,fq);
        fputc('	',fq);
      }
      fputc('\n',fq);
      for(int j=0;j<miscontig[i].thirdcontig2.read1.size();j++)
      {
        fputs(miscontig[i].thirdcontig2.read1[j].name,fq);
        fputc(' ',fq);
        sprintf(chpos,"%d",miscontig[i].thirdcontig2.read1[j].pos);
        fputs(chpos,fq);
        fputc('	',fq);
      }
      fputc('\n',fq);
  }
  for(int i=0;i<miscount;i++)
  {
      for(int j=0;j<miscontig[i].contig1.read2.size();j++)
      {
        fputs(miscontig[i].contig1.read2[j].name,fq);
        fputc(' ',fq);
        sprintf(chpos,"%d",miscontig[i].contig1.read2[j].pos);
        fputs(chpos,fq);
        fputc('	',fq);
      }
      fputc('\n',fq);
      for(int j=0;j<miscontig[i].contig2.read2.size();j++)
      {
        fputs(miscontig[i].contig2.read2[j].name,fq);
        fputc(' ',fq);
        sprintf(chpos,"%d",miscontig[i].contig2.read2[j].pos);
        fputs(chpos,fq);
        fputc('	',fq);
      }
      fputc('\n',fq);
      for(int j=0;j<miscontig[i].thirdcontig1.read2.size();j++)
      {
        fputs(miscontig[i].thirdcontig1.read2[j].name,fq);
        fputc(' ',fq);
        sprintf(chpos,"%d",miscontig[i].thirdcontig1.read2[j].pos);
        fputs(chpos,fq);
        fputc('	',fq);
      }
      fputc('\n',fq);
      for(int j=0;j<miscontig[i].thirdcontig2.read2.size();j++)
      {
        fputs(miscontig[i].thirdcontig2.read2[j].name,fq);
        fputc(' ',fq);
        sprintf(chpos,"%d",miscontig[i].thirdcontig2.read2[j].pos);
        fputs(chpos,fq);
        fputc('	',fq);
      }
      fputc('\n',fq);
  }
  fclose(fp);
  fclose(fp1);
  fclose(fq);  
//   cout<<"detect "<<misassemblynum<<" misassembly error"<<endl;;
//  fclose(fq1);
//  cout<<tempcount<<endl;
  
  return 0;
}

