#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<unistd.h>
#include<iostream>
#include<vector>
using namespace std;
typedef  struct  Read{
         int  pos;
         char name[100];
         int  direction;
}Read;
typedef  struct  MyContig{
        //   char  cigar[100000];
           vector<Read>   read1;
           vector<Read>   read2;
           int   length;
           char  name[100];
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
#define  percent  0.20
//char  qname[1000000],qname1[1000000]; 

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
   int count=0;
   int length;

   length=strlen(seq);
   for(int i=length-1;i>=0;i--)
   {
      switch(seq[i])
      {
        case 'A':moseq[length-i-1]='T';break;
        case 'T':moseq[length-i-1]='A';break;
        case 'G':moseq[length-i-1]='C';break;
        case 'C':moseq[length-i-1]='G';break;
      }
   }
   moseq[length]='\0';
     // cout<<seq[0]<<endl;
     // cout<<moseq[length-1]<<endl;
     // getchar();
   return 0;
}
 
 
int judgeBypos(Contig (&contig)[200000],int num,Miscontig  &(miscontig),int misnum)
{
   int lo_flag=0;
   int
 temp1,temp2;
   int temp[4];
   int distance=0;
   int lo_dis;
   int  j=0;
   if(classify(contig,num,miscontig)>1)
   {
     for(int i=0; i<miscontig.kind-1; i++)
     {      
        j=i+1;
    //    cout<<contig[i].refpos1<<contig[i].refpos2<<endl;
        if(strcmp(contig[i].refname,contig[j].refname)==0)
        {
           if(contig[i].flag!=contig[j].flag&&contig[j].flag!=0&&contig[j].flag!=0)
           {
              if((contig[i].refpos1<contig[j].refpos1)&&(contig[i].refpos2<contig[j].refpos2))
              { 
                  if((contig[i].direction==0)&&(contig[j].direction==0))
                  {
              // cout<<i<<" "<<contig[i].direction<<" "<<j<<" "<<contig[j].direction<<endl;
                    distance=abs((contig[i].refpos2-contig[j].refpos1));
               //cout<<"mark"<<endl;
                  }
                  if((contig[i].direction==0)&&(contig[j].direction==1))
                  {
                    distance=abs((contig[i].refpos2-contig[j].refpos2));
                  }
                  if((contig[i].direction==1)&&(contig[j].direction==0))
                  {
                    distance=abs((contig[i].refpos1-contig[j].refpos1));
                  }
                  if((contig[i].direction==1)&&(contig[j].direction==1))
                  {
                    distance=abs((contig[i].refpos1-contig[j].refpos2));
                  }
              }
              if((contig[i].refpos1>contig[j].refpos1)&&(contig[i].refpos2>contig[j].refpos2))
              {
                  if((contig[j].direction==0)&&(contig[i].direction==0))
                  {
             // cout<<i<<" "<<contig[i].direction<<" "<<j<<" "<<contig[j].direction<<endl;
                    distance=abs((contig[j].refpos2-contig[i].refpos1));
          //  cout<<"mark"<<endl;
                  }
                  if((contig[j].direction==0)&&(contig[i].direction==1))
                  {
                    distance=abs((contig[j].refpos2-contig[i].refpos2));
                  }
                  if((contig[j].direction==1)&&(contig[i].direction==0))
                  {
                    distance=abs((contig[j].refpos1-contig[i].refpos1));
                  }
                  if((contig[j].direction==1)&&(contig[i].direction==1))
                  {
                    distance=abs((contig[j].refpos1-contig[i].refpos2));
                  }
              }
              if(distance>85)
              {
                lo_flag=1;
                lo_dis=distance;
              }
           }
        }
        else
        {
           lo_flag=1;
           lo_dis=distance;
        }
  
  // cout<<contig[i].direction<<endl;  
   }

int   getPosOfgenome(int &length,Contig &contig)
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
//output  subcontigs  
int main(int argc,char *argv[])
{ 
  FILE *fp,*fp1,*fp2,*fq;
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
  int        miscount=0,refcount=0;
  char    tempname[100];
  Read    read;
  int direction;
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
  if((fq=fopen("misassembiles","w"))==NULL)
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
  while(ch!=EOF)
  {
    if(ch=='>')
    {
         ch=fgetc(fp2);
         count=0;
         while(ch!='\n')
         {
            miscontig[miscount].contig1.name[count++]=ch;
            ch=fgetc(fp2);
         }
         miscontig[miscount].contig1.name[count++]='\0';
         while(ch!='>')
         {
            if(ch!='\n')    chromosome.push_back(ch);
            ch=fgetc(fp2);
         }
       miscontig[miscount].contig1.length=chromosome.size();
   
//      cout<<miscontig[miscount].contig1.name<<endl; 
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
    }
    cout<<miscount<<endl;
  }
//  cout<<miscount<<endl;
  cout<<"load "<<miscount<<" contig"<<endl;
 // getchar();
//  strcpy(contig[connum++].qname,qname)




  ch=fgetc(fp);
  while(ch!=EOF)
  {
      if(ch=='@')
      {
        fgets(str,1000,fp);
        ch=fgetc(fp);
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
            ch=fgetc(fp);
          }
          qname[k]='\0';
/*
read flag  and convert
*/
          k=0;
          ch=fgetc(fp);
          while(ch!='	')
          {
            chflag[k++]=ch;
            ch=fgetc(fp);
          }
          chflag[k]='\0';
          flag=atoi(chflag);
/* 
 read  refname
*/ 
          k=0;
          ch=fgetc(fp);
          while(ch!='	')
          {
            refname[k++]=ch;
            ch=fgetc(fp);
          }
          refname[k]='\0';
/*
  contig's pos on  the reference  genome
*/          
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
  //        cout<<"load cigar"<<endl;
/*
 cigar
*/          
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
          ch=fgetc(fp);
          while(ch!='	')
          {
             ch=fgetc(fp);
          }
          
          ch=fgetc(fp);
          while(ch!='	')
          {
             ch=fgetc(fp);
          }
       
          ch=fgetc(fp);
          while(ch!='	')
          {
             ch=fgetc(fp);
          }

          k=0;
          ch=fgetc(fp);
          while(ch!='	')
          {
             ch=fgetc(fp);           
          }
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
           //  getchar();
       //  cout<<qname<<endl; 
          if(strcmp(qname,qname1)!=0)
          {
             connum++;
             if(connum==1000000)
             {
               thousand++;
               cout<<"deal"<<connum*thousand<<"reads"<<endl;
               cout<<miscontig[0].contig1.read1[thousand].pos<<endl;
               cout<<miscontig[0].contig1.read1[thousand].name<<endl;
               connum=0;
             }
           
          }
          
          strcpy(qname1,qname);
          fgets(str,1000000,fp);   
          ch=fgetc(fp);
       
      }
  }  
  cout<<refcount<<endl;
  getchar();

  ch=fgetc(fp1);
  while(ch!=EOF)
  {
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
          /*
          save info
*/
          ch=fgetc(fp1);
          while(ch!='	')
          {
             ch=fgetc(fp1);
          }
          
          ch=fgetc(fp1);
          while(ch!='	')
          {
             ch=fgetc(fp1);
          }
       
          ch=fgetc(fp1);
          while(ch!='	')
          {
             ch=fgetc(fp1);
          }

          k=0;
          ch=fgetc(fp1);
          while(ch!='	')
          {
             ch=fgetc(fp1);           
          }
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
//  fclose(fq1);
//  cout<<tempcount<<endl;
  
  return 0;
}
