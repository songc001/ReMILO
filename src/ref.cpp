#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<unistd.h>
#include<iostream>
#include<vector>
using namespace std;
typedef  struct  MyContig{
        //   char  cigar[100000];
           char  qname[100];
           char  refname[100];
           int   flag; 
           int   direction;//the direction of contig,1 represent reverse
           int   refpos1;
           int   refpos2;
        //   int   conpos;//0 represent pre,1 represent behind
           int   conBegin;// the begin of subcontig in contig
           int   conEnd;//the end of subcontig in contig
           int   tolen;
           int   cigarLen;
}Contig;
typedef   struct  Miscontig{
          Contig  contig[2000];
          int     flagnum[10];
          int     kind;
          int     num;
          char    seq[10000000];
}Miscontig;
typedef   struct  Index{
          char    str[1000];
          int     num;
          }Index;
#define  percent  0.20
//char  qname[1000000],qname1[1000000]; 

FILE  *fp,*fq,*fq1,*fp1;
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
}
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
     
     
   
//int sort(
int classify(Contig (&contig)[200000],int num,Miscontig  &(miscontig))
{
  int kind=0;
  int count[10000];
  int beginsum[10000],endsum[10000];
  int a[10000];
  int best;
  Contig tempcontig;
  for(int i=0;i<num;i++)
  {
    contig[i].flag=0;
  }
  for(int i=num-1;i>0;i--)
  {
    for(int j=0;j<i;j++)
    {
      if(abs(contig[j].conBegin-contig[j].conEnd)<abs(contig[j+1].conBegin-contig[j+1].conEnd))
      {
        tempcontig=contig[j];
        contig[j]=contig[j+1];
        contig[j+1]=tempcontig;
      }
    }
  }    
  for(int i=0;i<num;i++)
  {
    if(contig[i].flag==0)
    {
      contig[i].flag=kind+1;
      for(int j=i+1;j<num;j++)
      {
        if(contig[j].flag==0)
        {
          if(contig[i].conEnd<contig[j].conBegin)
          {
            contig[j].flag=0;
          }
          if(contig[i].conBegin>contig[j].conEnd)
          {
            contig[j].flag=0;

          }
          if((contig[i].conBegin<=contig[j].conBegin)&&(contig[i].conEnd<=contig[j].conEnd)&&(contig[i].conEnd>=contig[j].conBegin))
          {
             
             if( ( abs(contig[i].conEnd-contig[j].conBegin)<(contig[i].cigarLen*percent) )&&( abs(contig[i].conEnd-contig[j].conBegin)<(contig[i].cigarLen*percent) ) )
             {
               contig[j].flag=0;
             }
             else
             {
               contig[j].flag=kind+1;
             }
          }
          if((contig[i].conBegin>=contig[j].conBegin)&&(contig[i].conEnd>=contig[j].conEnd)&&(contig[i].conBegin<contig[j].conEnd)) 
          {
          
             if( (abs(contig[i].conBegin-contig[j].conEnd)<(contig[i].cigarLen*percent) )&&( abs(contig[i].conBegin-contig[j].conEnd)<(contig[i].cigarLen*percent) ) )
             {
               contig[j].flag=0;
             }
             else
             {
               contig[j].flag=kind+1;
             }
          }
          if((contig[i].conBegin>contig[j].conBegin)&&(contig[i].conEnd<contig[j].conEnd))
          {
             contig[j].flag=kind+1;
          }
          if((contig[i].conBegin<contig[j].conBegin)&&(contig[i].conEnd>contig[j].conEnd))
          {
             contig[j].flag=kind+1;
          }
        }
      }
    kind++;
    }
  }
 /* for(int j=0;j<num;j++)
  {
    cout<<contig[j].flag<<endl;
    cout<<contig[j].conBegin<<" "<<contig[j].conEnd<<endl;
  }
*/
  for(int i=0;i<=10000;i++)
  {
     count[i]=0;
     beginsum[i]=0;
     endsum[i]=0;
     a[i]=0;
  }
  miscontig.kind=0;
  for(int i=0;i<num;i++)
  {
     for(int j=1;j<=kind;j++)
     {
       if(contig[i].flag==j&&a[j]!=1)
       {
          a[j]=1;
          miscontig.contig[miscontig.kind++]=contig[i];
       }
     }
  }
if(kind>1)
{
  for(int i=miscontig.kind-1;i>0;i--)
  {
     for(int j=0;j<i;j++)
     {
        if(miscontig.contig[j].conBegin>miscontig.contig[j+1].conBegin)
        {
           tempcontig=miscontig.contig[j];
           miscontig.contig[j]=miscontig.contig[j+1];
           miscontig.contig[j+1]=tempcontig;
          //cout<<'1'<<endl;
        }
     }
   }
}
return kind;
}  
 
int judgeBypos(Contig (&contig)[200000],int num,Miscontig  &(miscontig),int misnum)
{
   int lo_flag=0;
   int temp1,temp2;
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
     }
  // cout<<contig[i].direction<<endl;  
   }
return lo_flag;     
}

int  dealCigar( int &pre, int  &be, int &length,char  cigar[],Contig &contig)  //length represent length of real seq
  {
      char  str[100];
      int   num;
      int   sum;
      int   i=0;//  cigar
      int   k=0;//  str
      int   count=0;;// number of H and S
      int   preflag=0;
      int   beflag=0;
      int   flag=0;
      int   temp;
      contig.cigarLen=0;
      be=pre=0;
      contig.conBegin=0;
      contig.conEnd=0;
      contig.tolen=0;
     // cout<<cigar<<endl;
      str[0]='\0';
      while(cigar[i]!='\0')
      {
	if(cigar[i]>=48&&cigar[i]<=57)
        {
          str[k++]=cigar[i];
        }
        else
        {
          str[k]='\0';
          k=0;
          num=convertInt(str);
       //   cout<<num;
        //  cout<<cigar[i];
          if(cigar[i]=='H')
          {  
             if(cigar[i+1]!='\0')
             {
               contig.conBegin+=num;
             }
             else
             {
               contig.conEnd+=num;
             }
          }
          if(cigar[i]=='S')
          {
             if(cigar[i+1]!='\0')
             {
               contig.conBegin+=num;
             }
             else
             {
               contig.conEnd+=num;
             }
          }
          if(cigar[i]=='M'||cigar[i]=='I'||cigar[i]=='='||cigar[i]=='X')
          {
             contig.cigarLen+=num;
          } 
          if(cigar[i]!='D')
          { 
            contig.tolen+=num;
          }
          str[0]='\0';
        }
        i++;
      }
      if(contig.direction==0)
      {
          contig.conEnd=contig.tolen-contig.conEnd-1;
      }
      if(contig.direction==1)
      {
          temp=contig.conBegin;
          contig.conBegin=contig.conEnd;
          contig.conEnd=contig.tolen-temp-1;
      }
     // cout<<contig.qname<<endl;
     // cout<<contig.conBegin<<" "<<contig.conEnd<<endl;
    //  getchar();
      return 0;
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
            
//output  subcontigs  
int main(int argc,char *argv[])
{ 
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
  int  k;
  int  count, mark, connum=0,misnum=0; 
  int tempcount=0; 
  int  precigar[2000],becigar[2000],lencigar[2000];
  vector<char>  chromosome[100]; 
  int chrcount=0,chrnum=0;
  Index   index[3000];
  Contig contig[200000];
  Miscontig  miscontig;

  if((fp=fopen(argv[1],"r"))==NULL)
  {
       cout<<"don't exit sam file:"<<argv[1]<<endl;
  }
  if((fp1=fopen(argv[2],"r"))==NULL)
  {
       cout<<"don't exit fasta file:"<<argv[2]<<endl;
  }
  if((fq=fopen(argv[3],"w"))==NULL)
  {
     printf("fq open error");
  }
  if((fq1=fopen(argv[4],"w"))==NULL)
  {
     printf("fq1 open error");
  }
  cout<<"load  sequence"<<endl;
  ch=fgetc(fp1);
  while(ch!=EOF)
  {
    if(ch=='>')
    {
       ch=fgetc(fp1);
       count=0;
       while(ch!='\n')
       {
          index[chrcount].str[count++]=ch;
          ch=fgetc(fp1);
       }
       index[chrcount].str[count]='\0';
  //      fgets(str,10000,fp1);
       index[chrcount].num=chrcount;
       cout<<chromosome[abs(chrcount-1)].size()<<endl;
       chrcount++;
    }
    if(ch!='\n')
    {
       chromosome[chrcount-1].push_back(ch);
    }
    ch=fgetc(fp1);
  }
  ch=fgetc(fp);
 cout<<"already load   sequence "<<endl;
//getchar();
//  strcpy(contig[connum++].qname,qname)
  while(ch!=EOF)
  {
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
          flag=convertInt(chflag);
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
          pos=convertInt(chpos);

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
             sequence[k++]=ch;
             ch=fgetc(fp);           
          }
          sequence[k]='\0';
            cout<<strlen(sequence)<<endl;
      //      getchar();
       //   cout<<"analyze"<<endl;
          if(strcmp(qname,qname1)==0)
          {
             if( (connum>=2)&&(judgeBypos(contig,connum,miscontig,misnum)==1) )
             {
           //    cout<<"input"<<endl; 
               tempcount++;
               fputc('>',fq);        
               fputs(contig[0].qname,fq);
               fputc('\n',fq);
               for(int i=0;i<(miscontig.kind-1)&&miscontig.contig[i].flag!=0;i++)
               {
               sprintf(tempnum,"%d",miscontig.contig[i].conBegin);
               fputs(tempnum,fq);
               fputc('	',fq);
               sprintf(tempnum,"%d",miscontig.contig[i].conEnd);
               fputs(tempnum,fq);
               fputc('	',fq);
               sprintf(tempnum,"%d",miscontig.contig[i+1].conBegin);
               fputs(tempnum,fq);
               fputc('	',fq);
               sprintf(tempnum,"%d",miscontig.contig[i+1].conEnd);
               fputs(tempnum,fq);
               fputc('\n',fq);
               }
               
               for(int i=0;i<miscontig.kind;i++)
               {
               //    subcontig1
                 fputc('>',fq1);
                 fputs(contig[0].qname,fq1);
                 fputc('_',fq1);
                 sprintf(tempnum,"%d",i);
                 fputs(tempnum,fq1);
            //     fputs("",fq1);
                 fputc('\n',fq1);
               //    cout<<strlen(miscontig.seq)<<endl;
               
                 for(int j=miscontig.contig[i].conBegin;j<=miscontig.contig[i].conEnd;j++)
                 {
                   //   cout<<miscontig.seq[j]<<endl;
                  //  getchar()
                    fputc(miscontig.seq[j],fq1);
                    if(((j-miscontig.contig[i].conBegin+1)%60)==0&&(j!=miscontig.contig[i].conBegin))
                    {
                       fputc('\n',fq1);
                    }
                 }
                 if(((miscontig.contig[i].conEnd-miscontig.contig[i].conBegin+1)%60)!=0)
                 {
                    fputc('\n',fq1);
                 }
            //       third contig
                 if(i<miscontig.kind-1)
                 {
             //       cout<<"tts"<<endl; 
                    fputc('>',fq1);
                    fputs(contig[0].qname,fq1);
          //          cout<<"tts"<<endl;
                    fputc('_',fq1); 
                    sprintf(tempnum,"%d",i);
                    fputs(tempnum,fq1);
                    fputs("_b",fq1);
                    fputc('\n',fq1);
             //       cout<<'F';
                    for(int k=0;k<chrcount;k++)
                    {
                       if(strcmp(index[k].str,miscontig.contig[i].refname)==0)
                       {
                          chrnum=index[k].num;
                          break;
                       }
                    }
                    if(miscontig.contig[i].direction==0)
                    {
                       if(miscontig.contig[i].refpos2+1000<chromosome[chrnum].size())
                       {
                          for(int j=miscontig.contig[i].refpos2;j<miscontig.contig[i].refpos2+1000;j++)
                          {
                   
                            fputc(chromosome[chrnum].at(j),fq1);
                            if(((j-miscontig.contig[i].refpos2+1)%60)==0&&(j!=miscontig.contig[i].refpos2+1000))
                            {
                              fputc('\n',fq1);
                            }
                          }
                       }
                       else
                       {
                          for(int j=miscontig.contig[i].refpos2;j<chromosome[chrnum].size();j++)
                          {
                            fputc(chromosome[chrnum].at(j),fq1);
                            if(((j-miscontig.contig[i].refpos2+1)%60)==0&&(j!=chromosome[chrnum].size()))
                            {
                              fputc('\n',fq1);
                            }
                          }
                       }
                          
                    }
                    else
                    {
                       if(miscontig.contig[i].refpos1>1000)
                       {  
                         cout<<miscontig.contig[i].refname<<endl; 
                         for(int j=miscontig.contig[i].refpos1-1000;j<miscontig.contig[i].refpos1;j++)
                         {
                            fputc(chromosome[chrnum].at(j),fq1);
                            if((j-miscontig.contig[i].refpos1+1001)%60==0&&(j!=(miscontig.contig[i].refpos1)))
                            {
                              fputc('\n',fq1);
                            }
                         }
                       }
                       else
                       {
                 ///        cout<<"H"<<endl;
                         for(int j=1;j<miscontig.contig[i].refpos1;j++)
                         {
                            fputc(chromosome[chrnum].at(j),fq1);
                            if((j+1)%60==0&&(j!=(miscontig.contig[i].refpos1)))
                            {
                              fputc('\n',fq1);
                            }
                         }
                       }
                         
                    }
                    fputc('\n',fq1); 
                    cout<<"second"<<endl;
                    fputc('>',fq1);
                    fputs(contig[0].qname,fq1);
                    fputc('_',fq1);
                    sprintf(tempnum,"%d",i+1);
                    fputs(tempnum,fq1);
                    fputs("_f",fq1);
                    fputc('\n',fq1);
                        
                    for(int k=0;k<chrcount;k++)
                    {
                       if(strcmp(index[k].str,miscontig.contig[i+1].refname)==0)
                       {
                          chrnum=index[k].num;
                          break;
                       }
                    }
            //        cout<<"cao"<<endl;
                    if(miscontig.contig[i+1].direction==0)
                    {
                       if(miscontig.contig[i+1].refpos1>1000)
                       {
                         cout<<"L"<<miscontig.contig[i+1].refpos1<<endl;
                         cout<<chromosome[chrnum].size()<<endl;
                          for(int j=miscontig.contig[i+1].refpos1-1000;j<=miscontig.contig[i+1].refpos1;j++)
                          {
                            fputc(chromosome[chrnum].at(j),fq1);
                            if(((j-miscontig.contig[i+1].refpos1+1001)%60)==0&&(j!=(miscontig.contig[i+1].refpos1)))
                            {
                              fputc('\n',fq1);
                            }
                          }
                          
                       }
                       else
                       {
                          cout<<"L2"<<endl;
                        
                          for(int j=0;j<miscontig.contig[i+1].refpos1;j++)
                          {
                            fputc(chromosome[chrnum].at(j),fq1);
                            if(((j+1)%60)==0&&(j!=miscontig.contig[i].refpos1))
                            {
                              fputc('\n',fq1);
                            }
                          }
                       }
                          
                    }
                    else
                    {                  
              
                       if(miscontig.contig[i+1].refpos2+1000<chromosome[chrnum].size())
                       {
                          for(int j=miscontig.contig[i+1].refpos2;j<=miscontig.contig[i+1].refpos2+1000;j++)
                          {
                            fputc(chromosome[chrnum].at(j),fq1);
                            if(((j-miscontig.contig[i+1].refpos2+1)%60)==0&&(j!=(miscontig.contig[i+1].refpos2+1000)))
                            {
                              fputc('\n',fq1);
                            }
                          }
                       }
                       else
                       {
                          for(int j=miscontig.contig[i+1].refpos2;j<chromosome[chrnum].size();j++)
                          {
                            fputc(chromosome[chrnum].at(j),fq1);
                            if(((j-miscontig.contig[i+1].refpos2+1)%60)==0&&(j!=chromosome[chrnum].size()))
                            {
                              fputc('\n',fq1);
                            }
                          }
                       }
                          
                    }
                    fputc('\n',fq1); 
 
                 }
                 
               }     
                
               
             }
             connum=0;
          }
          contig[connum].direction=judgeDirect(flag);
          strcpy(contig[connum].qname,qname);
          strcpy(contig[connum].refname,refname);
    //      cout<<"dealcigar"<<endl;
          dealCigar(precigar[connum],becigar[connum],lencigar[connum],cigar,contig[connum]);
      //    cout<<"getpos"<<endl;
          contig[connum].refpos1=pos;
          getPosOfgenome(lencigar[connum],contig[connum]);
   //       cout<<contig[connum].refpos1<<contig[connum].refpos2<<endl;
          if(connum==0)
          {
          //  cout<<contig[0].qname<<' '<<strlen(sequence)<<' '; 
            if(contig[connum].direction==1)
            {
               dealSeq(sequence,miscontig.seq);
            }
            else
            {
               strcpy(miscontig.seq,sequence);
            }
       //     cout<<"tolen "<<contig[0].tolen<<endl;
        //    cout<<"selen"<<strlen(sequence)<<endl;
       //     getchar();
          }
          if(contig[connum].cigarLen>100)
          {
            connum++;
          }
     //    cout<<qname<<endl; 
          strcpy(qname1,qname);
          fgets(str,1000000,fp);   
          ch=fgetc(fp);
          cout<<connum<<endl;
          cout<<contig[connum-1].qname<<endl;
      }
  }  
  fclose(fp);
  fclose(fq);  
  fclose(fq1);
//  cout<<tempcount<<endl;
  return 0;
}
