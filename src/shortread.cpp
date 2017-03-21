#include<iostream>
#include<fstream>
#include<cstring>
#include<vector>
#include<cstdlib>
#include<ctime>
#include<cmath>

#define   maxBuflen  1000000
#define   MAPRATE  0.4
#define   PERRATE  0.6 
using namespace std;

typedef struct  Read{
    string name;
    int pos;
    int direction;
    float maprate;
    float perfectrate;
}Read;

typedef struct  SubContig{
    string name;
    Read   read1;
    Read   read2;
    int    length;
}SubContig;
    
typedef struct  Contig{
    string   name;
    vector<SubContig>  sublist;
}Contig;

vector<Contig>    contiglist;

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
    read.maprate=(length-sc-hc)/(length*1.0);
    read.perfectrate=mc/( (length-sc-hc)*1.0 );    
}
        
void  assignRead(vector<string> data,Read &read)
{ 
    read.name=data[0];
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

void linkReadToContig(vector<string>  data)
{
    Read  read;
    cout<<data[0]<<endl;
    assignRead(data,read);
    if(read.maprate>MAPRATE  &&  read.perfectrate>PERRATE)
    {
      cout<<read.maprate<<endl;
    }
    getchar();

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
           linkReadToContig(data);
        }
        data.clear();
    }    
}
        
int main(int argc,char *argv[])
{
    ifstream  r1(argv[1]);
    ifstream  r2(argv[2]);
    ifstream  c(argv[3]);
    ofstream outmis(argv[4]);
    time_t   start,end;
    start=time(NULL); 
    if(r1.is_open())
    {
        cout<<"file1"<<endl;
        getRfile(r1);         
    }
    else
    {
        cout<<"can't find file "<<argv[1]<<endl;
        exit(0);
    }
    if(!r2.is_open())
    {
        getRfile(r2);
    }
    {
        cout<<"can't find file "<<argv[2]<<endl;
        exit(0);
    }
    if(!outmis.is_open())
    {
        cout<<"can't find file "<<argv[3]<<endl;
        exit(0);
    }
    end=time(NULL);
    cout<<"program  run about  "<<end-start<<"seconds"<<endl; 
  
}
