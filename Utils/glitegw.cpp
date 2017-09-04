/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) <year>  <name of author>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

*/

#include "glitegw.h"
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include "visivoutils.h"

#ifdef WIN32
  #include <time.h>
#endif

//---------------------------------------------------------------------
std::string gLiteGw::readVF (std::map<std::string, std::string> arguments,int idOp)
//---------------------------------------------------------------------
{
  std::string localFilename;
  std::map<std::string, std::string>::iterator iter;
   std::multimap<std::string, std::string>::iterator iter1;
 

  std::stringstream sstmp1;
  std::string  randat;
  time_t rawtime;
  struct tm * timeinfo;
  char buffer [80];
  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  strftime (buffer,80,"%Y%m%d%H%M%S",timeinfo);
  sstmp1<<"_"<<rand()<<buffer<<"_";  
  randat=sstmp1.str();

  switch(idOp)
{
/*** Merge **/
  case 3:
  case 8:
  {
    std::multimap<std::string, std::string> newFileList;
    iter=arguments.find("filelist");
    if(iter==arguments.end()) return "";
    std::string fileList=iter->second;

    std::ifstream fileInput(fileList.c_str());
    if(!fileInput.is_open()) return "";

    while(!fileInput.eof())
	{
		std::string tab, col;
		fileInput >> tab;
		fileInput >> col;
		if(tab.substr(0,6) =="lfn://")
		{
		    if(tab.find(".bin") == std::string::npos)
   			tab.append(".bin");
		    iter=arguments.find("VO");
		    if( iter == arguments.end())
		    {
			std::cerr<<"VO parameter is not given. Reading from gLite catalogue aborted."<<std::endl;
			fileInput.close();
			return "NogLiteSuccess";
		    }  
		    std::string vo(iter->second);
		    localFilename=tab.substr(7,tab.size()-1);
		    localFilename=getName(localFilename);
		    localFilename="gL"+randat+localFilename;
		    iter=arguments.find("out");
		    if( iter == arguments.end())
		    {
		      std::string dir;
		      dir=getenv("PWD");
		      localFilename=dir+"/"+localFilename;
		    }else{
		      std::string dir=getDir(iter->second);
		      if(dir=="")
		      {
			dir=getenv("PWD");
			localFilename=dir+"/"+localFilename;
		      } else
			localFilename=dir+localFilename;
		    }  
		    newFileList.insert(make_pair(localFilename,col));
		    iter=m_lfnDownloadedFile.find(tab);
		    if(iter==m_lfnDownloadedFile.end())
		    {
		      if(!remoteDownloadFiles(tab,vo,localFilename)) 
			return "NogLiteSuccess";
		      m_lfnDownloadedFile.insert(make_pair(tab,localFilename));
		      tab=tab+".head";
		      localFilename=localFilename+".head";
		      if(!remoteDownloadFiles(tab,vo,localFilename)) 
			 return "NogLiteSuccess";
		      m_lfnDownloadedFile.insert(make_pair(tab,localFilename));
		    }
		} else 
		    newFileList.insert(make_pair(tab,col));
	}  
	fileInput.close();

	std::string filegLiteGw="case3"+randat+".txt";
        std::ofstream fileOut(filegLiteGw.c_str());
	if(!fileOut.is_open())
	{
	  std::cerr<<"No write permission for "<<filegLiteGw<<" Operation aborted."<<std::endl;
	  return "NogLiteSuccess";
	}  
	for(iter1=newFileList.begin();iter1!=newFileList.end();iter1++)
	{ 
	  if(iter1->first !="")
	  {
	   fileOut<<iter1->first;
	   fileOut<<"         ";
	   fileOut<<iter1->second;
	   fileOut<<std::endl;
	  }
	}
        localFilename=filegLiteGw; //used for return value 
    break;
  }

/*** Append **/
  case 4:
  {
    std::multimap<std::string, std::string> newFileList;
    iter=arguments.find("filelist");
    if(iter==arguments.end()) return "";
    std::string fileList=iter->second;

    std::ifstream fileInput(fileList.c_str());
    if(!fileInput.is_open()) return "";

    while(!fileInput.eof())
	{
		std::string tab;
		fileInput >> tab;
		if(tab.substr(0,6) =="lfn://")
		{
		    if(tab.find(".bin") == std::string::npos)
   			tab.append(".bin");
		    iter=arguments.find("VO");
		    if( iter == arguments.end())
		    {
			std::cerr<<"VO parameter is not given. Reading from gLite catalogue aborted."<<std::endl;
			fileInput.close();
			return "NogLiteSuccess";
		    }  
		    std::string vo(iter->second);
		    localFilename=tab.substr(7,tab.size()-1);
		    localFilename=getName(localFilename);
		    localFilename="gL"+randat+localFilename;
		    iter=arguments.find("out");
		    if( iter == arguments.end())
		    {
		      std::string dir;
		      dir=getenv("PWD");
		      localFilename=dir+"/"+localFilename;
		    }else{
		      std::string dir=getDir(iter->second);
		      if(dir=="")
		      {
			dir=getenv("PWD");
			localFilename=dir+"/"+localFilename;
		      } else
			localFilename=dir+localFilename;
		    }  
		    newFileList.insert(make_pair(localFilename,""));
		    iter=m_lfnDownloadedFile.find(tab);
		    if(iter==m_lfnDownloadedFile.end())
		    {
		      if(!remoteDownloadFiles(tab,vo,localFilename)) 
			return "NogLiteSuccess";
		      m_lfnDownloadedFile.insert(make_pair(tab,localFilename));
		      tab=tab+".head";
		      localFilename=localFilename+".head";
		      if(!remoteDownloadFiles(tab,vo,localFilename)) 
			 return "NogLiteSuccess";
		      m_lfnDownloadedFile.insert(make_pair(tab,localFilename));
		    }
		} else 
		    newFileList.insert(make_pair(tab,""));
	}  
	fileInput.close();

	std::string filegLiteGw="case4"+randat+".txt";
        std::ofstream fileOut(filegLiteGw.c_str());
	if(!fileOut.is_open())
	{
	  std::cerr<<"No write permission for "<<filegLiteGw<<" Operation aborted."<<std::endl;
	  return "NogLiteSuccess";
	}  
	for(iter1=newFileList.begin();iter1!=newFileList.end();iter1++)
	{ 
	  if(iter1->first !="")
	  {
	   fileOut<<iter1->first;
	   fileOut<<std::endl;
	  }
	}
        localFilename=filegLiteGw; //used for return value 
    break;
  }
  case 17:
  {
    std::stringstream ssInfileparameters;
    iter=arguments.find("infiles");
    ssInfileparameters.str(iter->second);
    std::string file1[2];
    ssInfileparameters>>file1[0];
    ssInfileparameters>>file1[1];
    if(file1[0]=="" || file1[1]=="")
    {
      localFilename=iter->second;
      break;
    }

    std::string returnString="";
   
    for(int i=0;i<2;i++)
	{
		if(i==1) returnString=returnString+" ";
		if(file1[i].substr(0,6) =="lfn://")
		{
		    if(file1[i].find(".bin") == std::string::npos)
   			file1[i].append(".bin");
		    iter=arguments.find("VO");
		    if( iter == arguments.end())
		    {
			std::cerr<<"VO parameter is not given. Reading from gLite catalogue aborted."<<std::endl;
			return "NogLiteSuccess";
		    }  
		    std::string vo(iter->second);
		    localFilename=file1[i].substr(7,file1[i].size()-1);
		    localFilename=getName(localFilename);
		    localFilename="gL"+randat+localFilename;
		    iter=arguments.find("out");
		    if( iter == arguments.end())
		    {
		      std::string dir;
		      dir=getenv("PWD");
		      localFilename=dir+"/"+localFilename;
		    }else{
		      std::string dir=getDir(iter->second);
		      if(dir=="")
		      {
			dir=getenv("PWD");
			localFilename=dir+"/"+localFilename;
		      } else
			localFilename=dir+localFilename;
		    }  
		    if(!remoteDownloadFiles(file1[i],vo,localFilename)) 
			return "NogLiteSuccess";
		    returnString=returnString+file1[i];
		    m_lfnDownloadedFile.insert(make_pair(file1[i],localFilename));
		    file1[i]=file1[i]+".head";
		    localFilename=localFilename+".head";
		    if(!remoteDownloadFiles(file1[i],vo,localFilename)) 
			 return "NogLiteSuccess";
		    m_lfnDownloadedFile.insert(make_pair(file1[i],localFilename));
		} else 
		    returnString=returnString+file1[i];
	}  


        localFilename=returnString; //used for return value 
    break;
  }
/*** Grid2Point ***/

case 10000:
  {
    break;
  }

/*** Default **/
  case -1:
  default:
  {
    iter=arguments.find("file");
    if( iter == arguments.end())
    {
      std::cerr<<"file parameter is not given. Reading from gLite aborted."<<std::endl;
      return "";
    }  
    std::string lfnFilename(iter->second);    
    if(lfnFilename.substr(0,6) !="lfn://")
    {                
      std::cerr<<"A localfile is given. VisIVO read from local filesystem"<<std::endl;
      return lfnFilename;
    }  
    
    iter=arguments.find("VO");
    if( iter == arguments.end())
    {
      std::cerr<<"VO parameter is not given. Reading from gLite catalogue aborted."<<std::endl;
      return "NogLiteSuccess";
    }  
    std::string vo(iter->second);

    localFilename=lfnFilename.substr(7,lfnFilename.size()-1); 
    localFilename=getName(localFilename);
//    std::clog<<"t1: localFilename="<<localFilename<<std::endl;        
    iter=arguments.find("out");
    if( iter == arguments.end())
    {
       std::string dir;
       dir=getenv("PWD");
       localFilename=dir+"/"+localFilename;
//       std::clog<<"t2: localFilename="<<localFilename<<" dir="<<dir<<std::endl;        
    }else{
      std::string dir=getDir(iter->second);
      if(dir=="")
      {
       dir=getenv("PWD");
       localFilename=dir+"/"+localFilename;	
      } else
        localFilename=dir+"/"+localFilename;
//       std::clog<<"t3: localFilename="<<localFilename<<" dir="<<dir<<std::endl;        
    }  
//     std::clog<<"t4: localFilename="<<localFilename<<std::endl;        
   
    if(!remoteDownloadFiles(lfnFilename,vo,localFilename)) 
      return "NogLiteSuccess";
    std::string header=localFilename+".head";
    lfnFilename=lfnFilename+".head";
//     std::clog<<"t5: lfnFilename="<<lfnFilename<<" header:"<<header<<std::endl;        
    if(!remoteDownloadFiles(lfnFilename,vo,header))
    {
      remove(localFilename.c_str());
      return "NogLiteSuccess";
    }
    break;
  }
}
  return localFilename;
}
//---------------------------------------------------------------------
bool gLiteGw::writeVF (std::map<std::string, std::string> arguments,int idOp)
//---------------------------------------------------------------------
{  
std::map<std::string, std::string>::iterator iter;
switch(idOp)
{
/*** Default **/
  case 10000:
  {
    break;
  }
  case -1:
  default:
  { 
    bool isavbt=true;
    iter=arguments.find("realOutFilename");
    if( iter == arguments.end())
    {
      std::cerr<<"file output name is not available. Writing into gLite catalogue aborted."<<std::endl;
      return false;
    }  
    std::string realOutFilename=iter->second;
    iter=arguments.find("lfnout");
    if( iter == arguments.end())
    {
      std::cerr<<"lfn not available. Writing into gLite catalogue aborted."<<std::endl;
      return false;
    }  
    std::string lfnFilename=iter->second;
    iter=arguments.find("VO");
    if( iter == arguments.end())
    {
      std::cerr<<"VO not available. Writing into gLite catalogue aborted."<<std::endl;
      return false;
    }  
    std::string vo=iter->second;
    std::string se;
    iter=arguments.find("se");
    if( iter != arguments.end()) se=iter->second;
    if(idOp!=10 && idOp!=11 )  // --out is a VBT
    {  
    if(lfnFilename.find(".bin") == std::string::npos)
   			lfnFilename.append(".bin");
    } else
      isavbt=false;
    if(!remoteLfn(realOutFilename,se, vo,lfnFilename,isavbt))
      return false;
    break;
  }
}
  return true;
}

//---------------------------------------------------------------------
void gLiteGw::rmv ()
//---------------------------------------------------------------------
// remove downloaded temporary files
{
std::map<std::string, std::string>::iterator iter;

for(iter=m_lfnDownloadedFile.begin();iter!=m_lfnDownloadedFile.end();iter++)
 remove(iter->second.c_str());

return;
}