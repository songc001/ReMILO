#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<iostream>
#include<fstream>
#include<vector>
#include<string>
#include<algorithm>
using namespace std;
 
#define maxBuflen  10000000
#define distance   85
#define maxQuality 1
#define percent   0.3
typedef struct Chr {
     string name;
     string sequence;
}Chr;
typedef struct Op{
     char type;
     int length;
}Op;//storage the type of cigar;
typedef struct subContig{
     int refpos;
     int conpos;
     int direction;
     int length;
     string chrName;
     string sequence;
}subContig;
typedef struct Contig{
     string  name;
     int     length;
     vector<subContig>  subContigList;
}Contig;
typedef struct Seg{
     string refname;
     int begin;
     int end;
     int refbegin;
     int refend;
     int direction;
//     bool type;
     double quality;
     bool  flag;//indels(true) or struction variation(0);
     int kind;
}Seg;
typedef struct MisSeg{
     Seg   seg1;
     Seg   seg2;
}MisSeg;
/*
typedef struct Segment{
     int pos1;
     int pos2;
     int refpos1;
     int refpos2;
     int direction;
     double quality;
}Segment;*/
typedef struct TSegment{
     vector<Seg>  seglist;
     vector<Seg>  indelseglist;
     double quality;
}TSegment;
vector<Chr>  reference;  
vector<Contig>   contiglist;
vector<Op> oplist;
//vector<Segment> segmentlist;
vector<TSegment>  tsegmentlist; 
int COUNT=0;

string GetChrName(string str)
{
   string name;
   if(str[0]!='>')    
   {
      cout<<"The name of  reference is not correct"<<endl;
      exit(1);
   }
   for(int i=1;i<str.size();i++)
   {
      if(str[i]=='	'||str[i]==' ')
      {
         break;
      }
      name+=str[i];
   }
   return name;
}
void getFfile(ifstream  &f)
{
    string buffer;
    string sequence;
    Chr    chr;
    while(f.good())
    {
       getline(f,buffer);
      // cout<<'1'<<endl;
       
       if(buffer[0]=='>')
       {
          if(chr.sequence!="")
          {
             reference.push_back(chr);
             chr.name.clear();
             chr.sequence.clear();
          }
          chr.name=GetChrName(buffer);
       }
       else
       {
          chr.sequence+=buffer;
          buffer.clear();
       }
    }
    reference.push_back(chr);
}

int JudgeDirection(int direction)
{
    return (direction>>4)%2;
}
vector<string>  split(string str,const char * delim)
{
    vector<string> data;
    char *p;
    char buffer[maxBuflen];
    strcpy(buffer,str.c_str());
    p=strtok(buffer,delim);
    while(p)
    {
       data.push_back(p);
       p=strtok(NULL,delim);
    }
    return data;
}
void AssignContig(Contig  &contig,vector<string> data)
{
   						
}
void ClearContig(Contig   &contig)
{
    contig.name="";
    contig.length=0;
    contig.subContigList.clear();
}
void AssignSubContig(subContig  &subcontig,vector<string> data)
{
    subcontig.direction=JudgeDirection(atoi(data[1].c_str()));
    subcontig.chrName=data[2];
    subcontig.refpos=atoi(data[3].c_str());
    
}
void ClearSubContig(subContig  &subcontig)
{
    subcontig.refpos=0;
    subcontig.conpos=0;
    subcontig.direction=0;
    subcontig.length=0;
    subcontig.chrName="";
    subcontig.sequence="";
}

double  judge_1(double a)
{
    if(a==1.0)   return 1.0;
    else
       return 0;   
}
void   AssignDirection(string direc,vector<Seg>  &seglist)
{
    int direction;
    direction=JudgeDirection( atoi( direc.c_str() ) );
    for(int i=0;i<seglist.size();i++)
    {
       seglist[i].direction=direction;
    }
}
double CalcuQua(int M,int indel)
{
    return  M/( (M+indel)*1.0 );
}
vector<int>  CalcuCigarInfo(vector<string> data)
{
    vector<int> counts(10,0);//counts[0]: total length of  contig;counts[1]:the length of refgenome which is aligned by contig;counts[3]the number of 'I','D','=',
    //vector<char> clip;//counts[2]:the number of 'M';
    Seg  seg;
    char ch; 
    int length;
    int refpos;
    TSegment tsegment;
    Op  op;
    vector<Seg>  seglist;
    vector<Seg>  indelseglist;
    refpos=atoi(data[3].c_str());
    seg.begin=counts[0];
    seg.quality=0;
    counts[1]=refpos;
    seg.refbegin=counts[1];
    tsegment.quality=0;
    seg.flag=false;
    for(int i=0;i<oplist.size()-1;i++)//deal n-1 var;
    {
       ch=oplist[i].type;
       length=oplist[i].length;
       if(ch=='M')
       {
          counts[0]+=length;
          counts[1]+=length;
          counts[2]+=length;
       }
       if(ch=='I')
       {
          counts[0]+=length;
          counts[3]+=length;
          if(length>distance)
          {
             seg.end=counts[0]-length;
             seg.quality=CalcuQua(counts[2],counts[3]-length);
             seg.flag=true;
             seg.refend=counts[1];
             seg.refname=data[2];
             indelseglist.push_back(seg);
             counts[2]=counts[3]=0;
             seg.begin=counts[0];
             seg.flag=true;
             seg.refbegin=counts[1]+1;
            // seg.type=true;
          }
       }
       if(ch=='D')
       {
          // data[0]+=length;
          counts[1]+=length;
          counts[3]+=length;
          if(length>distance)
          {
           //  cout<<length<<endl;
             seg.end=counts[0];
             seg.quality=CalcuQua(counts[2],counts[3]-length);
             seg.flag=true;
             seg.refend=counts[1]-length;
             seg.refname=data[2];
             indelseglist.push_back(seg);
             counts[2]=counts[3]=0;
             seg.begin=counts[0]+1;
             seg.flag=true;
             seg.refbegin=counts[1];
          }
       }
       if(ch=='S'||ch=='H')
       {
          counts[0]+=length;
          if(counts[0]==length)
          {
             seg.begin=counts[0];
            // seg.quality=CalcuQua(counts[2],counts[3]-length);
             seg.flag=false;
             seg.refbegin=counts[1];
          }
          else
          {
             seg.end=counts[0]-length;
             seg.quality=CalcuQua(counts[2],counts[3]-length);
             if(seg.flag!=true)   seg.flag=false;
             seg.refend=counts[1];
             seg.refname=data[2];
             if(seg.flag==true)   
             {
                indelseglist.push_back(seg);
             }
             else
             {
                seglist.push_back(seg);
             }
             counts[2]=counts[3]=0;
             seg.begin=counts[0];
             seg.flag=false;
             seg.refbegin=counts[1]+1;
          }
       }
       if(ch=='=')
       {
          counts[0]+=length;
          counts[1]+=length;
          counts[2]+=length;
       }          
       if(ch=='X')
       {
          counts[0]+=length;
          counts[1]+=length;
          counts[3]+=length;
          cout<<data[0]<<endl;///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
          cout<<"X"<<endl;////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
          getchar();
       }   
    } 
    op=oplist[oplist.size()-1];
    if(op.type=='H'||op.type=='S')
    {
       seg.end=counts[0];
       seg.quality=CalcuQua(counts[2],counts[3]);
       if(seg.flag!=true)      seg.flag=false;
       seg.refend=counts[1];
       seg.refname=data[2];
       if(seg.flag==true)
       {
          indelseglist.push_back(seg);
       }
       else
       {
          seglist.push_back(seg);
       }
       counts[2]=counts[3]=0;
       counts[0]+=op.length;
    }
    else
    {
       counts[0]+=op.length;
       counts[1]+=op.length;
       counts[2]+=op.length;
       seg.end=counts[0];
       seg.quality=CalcuQua(counts[2],counts[3]);
       seg.refend=counts[1];
       seg.refname=data[2];
       if(seg.flag==true)
       {
          indelseglist.push_back(seg);
       }
       else
       {
          seglist.push_back(seg);
       }
       counts[2]=counts[3]=0;
       counts[0]+=op.length;
    }
 //   if(tsegment.quality!=1.0) tsegment.quality=CalcuQua(data[1],data[2]);
    if(!seglist.empty())    AssignDirection(data[1],seglist);
    if(!indelseglist.empty())  AssignDirection(data[1],indelseglist);
    if(!seglist.empty()||!indelseglist.empty())
    {
       tsegment.seglist=seglist;
       tsegment.indelseglist=indelseglist;
 //      tsegment.quality=0;//need to be modified;
       tsegmentlist.push_back(tsegment);
    }
    oplist.clear();
    seglist.clear();
    indelseglist.clear();
    return  counts;      
}
bool a_less_b(const Seg &a,const Seg &b)
{
    return a.begin<b.begin;
}
void ConvertCigar(string cigar)
{
    vector<int>  data;
    int i=0;
    string str;
    Op  op;
   //      cout<<"w"<<endl;
    while(i<cigar.size())
    {
       if(cigar[i]>=48&&cigar[i]<=57)
       {
          str+=cigar[i];
       }
       else
       {
          op.type=cigar[i];
          op.length=atoi(str.c_str());
          oplist.push_back(op);
          str="";
       }
       i++; 
    }
   // cout<<oplist.size()<<endl;
}
int ClassifySeg(Seg  seg1,Seg  seg2)
{
   int  samelen=0;
   int  seg1_len=0;
   int  seg2_len=0;
   samelen=abs(seg1.end-seg2.begin);
   seg1_len=abs(seg1.end-seg1.begin)*percent;
   seg2_len=abs(seg2.end-seg2.begin)*percent;
   if(seg2.kind!=0)    return seg2.kind;
   if(seg1.end<seg2.begin) return 0;
   if(seg1.begin>seg2.end)   return 0;
   if( (seg1.begin<=seg2.begin)&&(seg1.end<=seg2.end)&&(seg1.end>=seg2.begin) )
   {
      if(samelen<seg1_len&&samelen<seg2_len)   return 0;
      else
      {
         return seg1.kind;
      }
   }
   if( (seg1.begin>=seg2.begin)&&(seg1.end>=seg2.end)&&(seg1.end<seg2.begin) )
   {  
      if(samelen<seg1_len&&samelen<seg2_len)   return 0;
      else
      {
         return seg1.kind;
      }
   }
   if( (seg1.begin>seg2.begin)&&(seg1.end<seg2.end) )
   {
      return  seg1.kind;
   }
   if( (seg1.begin<seg2.begin)&&(seg1.end>seg2.end) )
   {
      return  seg1.kind;
   }
}
vector<Seg>  MergeSeg(vector<Seg>  sortseglist,int kind)
{
    vector<Seg>  fseglist;
    Seg  tempseg;
    for(int i=1;i<=kind;i++)
    {
       tempseg.begin=tempseg.end=0;
       for(int j=0;j<sortseglist.size();j++)
       {
          if(sortseglist[j].kind==i && abs(sortseglist[j].end-sortseglist[j].begin)>abs(tempseg.end-tempseg.begin) )
          {
             tempseg=sortseglist[j];
          }
  //        cout<<'d'<<endl;
       }
       fseglist.push_back(tempseg);
    }
//    cout<<'e'<<endl;
    sort(fseglist.begin(),fseglist.end(),a_less_b);
//    cout<<'s'<<endl;
    return fseglist;
}
    
vector<MisSeg>  JudgeMis(vector<Seg>  fseglist)
{
   MisSeg  misseg;
   vector<MisSeg>  misseglist;
   for(int i=0;i<fseglist.size()-1;i++)
   {
      if(fseglist[i].refname!=fseglist[i+1].refname)   
      {
         misseg.seg1=fseglist[i];
         misseg.seg2=fseglist[i+1];
         misseglist.push_back(misseg);
      }
      else
      {
         if(abs(fseglist[i].refend-fseglist[i+1].refbegin)>distance)  
         {
            misseg.seg1=fseglist[i];
            misseg.seg2=fseglist[i+1];
            misseglist.push_back(misseg);
         }
      }
   }
   return misseglist;
}
void FormatOutputBreakPoints(vector<MisSeg> misseglist,ofstream &out)
{
   for(int i=0;i<misseglist.size();i++)
   {
      out<<misseglist[i].seg1.begin<<'	'<<misseglist[i].seg1.end<<'	'<<misseglist[i].seg2.begin<<'	'<<misseglist[i].seg2.end<<endl;
   }
}
Chr FindChr(string name)
{
   for(int i=0;i<reference.size();i++)
   {
      if(reference[i].name==name)
      {
         return  reference[i];
      }
   }
   cout<<"the sam file doesn't match fasta file"<<endl;
   exit(1);
}
string  CutRefSeq(int begin,int end,string  chrname,int direction) 
{
    Chr chr;
    string  seq="";
    int i=0;
    chr=FindChr(chrname);
    if(direction==0)
    {
       begin=end;
       end=end+1000;
    }
    else
    {
       end=begin;
       if( (begin-1000)>0 )   
       {
          begin=(begin-1000);
       }
       else
       {
          begin=0;
       }
    }
    i=begin;
    while(i<end&&i<chr.sequence.size())
    {
       seq+=chr.sequence[i];
       i++;
    }
    return  seq;
}
string CutConSeq(int begin,int end,string sequence)
{ 
    string seq="";
    for(int i=abs(begin-1);i<end;i++)
    {
       seq+=sequence[i];
    }
    return seq;        
}
void FormatOutput(string str,ofstream  &output)
{
    string buffer="";
    for(int i=0;i<str.size();i++)
    {
       buffer+=str[i];
       if(buffer.size()==60)
       {
          output<<buffer<<endl;
          buffer.clear();
       }
    }
    if( (str.size()%60)!=0 )
    {
       output<<buffer<<endl;
    }  
}
void OutputSegment(vector<MisSeg>  misseglist,string sequence,string conname,ofstream  &output)
{
    string  seq;;
    string  con;
    Seg  seg1,seg2;
    int count=0;
  //  cout<<misseglist.size()<<endl;
//    cout<<conname<<endl;
    for(int i=0;i<misseglist.size();i++)
    {
       seg1=misseglist[i].seg1;
       seg2=misseglist[i].seg2;
       if(seg2.direction==0)    seg2.direction=1;
       else
       {
          seg2.direction=0;
       }
     //  FormatOutput(sequence);
       if(i>=1)
       {  
          seq=CutRefSeq(seg1.refbegin,seg1.refend,seg1.refname,seg1.direction);
          output<<'>'<<conname<<"_"<<count<<"_b"<<endl;
          FormatOutput(seq,output);
          count++;
       }
       else
       {
          seq=CutConSeq(seg1.begin,seg1.end,sequence);
          output<<'>'<<conname<<"_"<<count<<endl;
          FormatOutput(seq,output);
          seq=CutRefSeq(seg1.refbegin,seg1.refend,seg1.refname,seg1.direction);
          output<<'>'<<conname<<"_"<<count<<"_b"<<endl;
          FormatOutput(seq,output);
          count++;
       }       
    //   cout<<i<<endl;
       seq=CutRefSeq(seg2.refbegin,seg2.refend,seg2.refname,seg2.direction);
       output<<'>'<<conname<<"_"<<count<<"_f"<<endl;
       FormatOutput(seq,output);
       output<<'>'<<conname<<"_"<<count<<endl;
       seq=CutConSeq(seg2.begin,seg2.end,sequence);
       FormatOutput(seq,output);
    //   getchar();
    } 
}
void InferSegment(string sequence,string  conname,ofstream &out,ofstream &output)
{
    vector<Seg> sortsegList;
    vector<Seg> fseglist;
    vector<MisSeg>  misseglist;
    int kind=0;
    for(int i=0;i<tsegmentlist.size();i++)
    {
      // if(tsegmentlist[i].segmentlist.size()<2)  continue;
       for(int j=0;j<tsegmentlist[i].indelseglist.size();j++)
       {
          sortsegList.push_back(tsegmentlist[i].indelseglist[j]);
          
       }
       for(int j=0;j<tsegmentlist[i].seglist.size();j++)
       {
          sortsegList.push_back(tsegmentlist[i].seglist[j]);
       }
    }
    sort(sortsegList.begin(),sortsegList.end(),a_less_b);
/*    for(int i=0;i<sortsegList.size();i++)
    {
       out<<sortsegList[i].begin<<' '<<sortsegList[i].end<<endl;
       out<<sortsegList[i].refbegin<<' '<<sortsegList[i].refend<<endl;
       out<<sortsegList[i].quality<<endl;
       out<<endl;
       
       //getchar();
    }
*/
    for(int i=0;i<sortsegList.size();i++)
    {
       sortsegList[i].kind=0;
    }
    for(int i=0;i<sortsegList.size();i++)
    {
       if(sortsegList[i].kind==0)
       {
          kind++;
          sortsegList[i].kind=kind; 
       }
       else    continue;
       for(int j=i+1;j<sortsegList.size();j++)
       {
          sortsegList[j].kind=ClassifySeg(sortsegList[i],sortsegList[j]);
       }
    }
    if(kind>1) 
    {
       fseglist=MergeSeg(sortsegList,kind);
       misseglist=JudgeMis(fseglist);
       FormatOutputBreakPoints(misseglist,out);
       OutputSegment(misseglist,sequence,conname,output);
    }
    misseglist.clear();
    tsegmentlist.clear();
}
void getCfile(ifstream  &c,ofstream  &out,ofstream  &output)															
{
    Contig  contig;
    subContig  subcontig; 
    string  buffer;
    string  name;  
    vector<string>  data; 
    string  sequence;
    while(c.good())
    {
       getline(c,buffer);
       if(buffer[0]!='@'&&buffer!="")
       {
          data=split(buffer,"	");
        //  cout<<data[0]<<endl; 
          if(data[2]=="*")
          {
              continue;
          }
          if(data[0]!=contig.name)
          {
             if(contig.name!="")
             {
                //contiglist.push_back(contig);
              //  judgeMis(contig);   
                if(tsegmentlist.size()>0)   
                {
                   out<<'>'<<contig.name<<endl;
                   InferSegment(sequence,contig.name,out,output);
                }
                ClearContig(contig);
             }
             contig.name=data[0];
             sequence=data[9];
          }
       //   cout<<data[5]<<endl;
          ConvertCigar(data[5]);
          CalcuCigarInfo(data);
  //        AssignSubContig(subcontig,data);

   //       contig.subContigList.push_back(subcontig);
   //       ClearSubContig(subcontig);
          data.clear();         
       }
       
    }
}
int main(int argc,char *argv[])
{
    ifstream  f(argv[1]);//contig file
    ifstream  c(argv[2]);//sam file
    ofstream  out(argv[3]);
    ofstream  output(argv[4]);
   // ifstream  mis(argv[3]);
  //  ifstream  split(argv[4]);
    //test();
    if(f.is_open())
    {
       getFfile(f);
  //     cout<<reference[0].sequence.length()<<endl;
   //    getchar(); 
    }
    if(c.is_open())
    {
       getCfile(c,out,output);
       cout<<"sam file loaded"<<endl;
       cout<<COUNT<<endl;
    }
    
}
  
