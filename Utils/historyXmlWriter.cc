//
//  historyXmlWriter.cpp
//  VisIVOServer_New072012
//
//  Created by Fabio Vitello on 27/06/13.
//
//


#include "historyXmlWriter.h"

#include "tinyxml.h"
#include "tinystr.h"

#include <cstdlib>
#include <cstring>

#include <iostream>
#include <fstream>
#include <sstream>



const float m_missing=-1.0918273645e+23;
const float m_text=-1.4536271809e+15;




//---------------------------------------------------------------------
HistoryXmlWriter::HistoryXmlWriter(const char* histFile,const char* format,const char* out,const char* tableOrVolume,double comput[],double size[],
                                   const char* login, const char* binaryHeader, float missing,float text, const char* endian,
                                   const char* type, long unsigned int points,
                                   const char* vo, const char* se, const char* lfnout, const char* inputFile)
//---------------------------------------------------------------------
{
    
    
    
    historyData["command"]="VisIVOImporter";
    historyData["fformat"]=format;
    historyData["file"]=inputFile;
    if (strcmp(out,"nopath")==0)
        historyData["out"]="VisIVOServerBinary.bin";
    else
        historyData["out"]=out;
    if (strcmp(tableOrVolume,"volume")==0)
    {
        historyData["volume"]="";

        
/*       
    historyData["compx"]=boost::lexical_cast<std::string>(comput[0]);
    historyData["compy"]=boost::lexical_cast<std::string>(comput[1]);
    historyData["compz"]=boost::lexical_cast<std::string>(comput[2]);
    historyData["sizex"]=boost::lexical_cast<std::string>(size[0]);
    historyData["sizey"]=boost::lexical_cast<std::string>(size[1]);
    historyData["sizez"]=boost::lexical_cast<std::string>(size[2]);
*/
        
        std::stringstream strs;
        strs << comput[0];
        historyData["compx"]=strs.str();
        
        strs.clear();
        strs << comput[1];
        historyData["compy"]=strs.str();
        
        
        strs.clear();
        strs << comput[2];
        historyData["compz"]=strs.str();
        
        strs.clear();
        strs << size[0];
        historyData["sizex"]=strs.str();
       
        strs.clear();
        strs << size[1];
        historyData["sizey"]=strs.str();
       
        strs.clear();
        strs << size[2];
        historyData["sizez"]=strs.str();
    }
    if (strcmp(login,"nouser")!=0)
    {
        historyData["userpwd"]=login;
    }
    if (strcmp(binaryHeader,"noheader")!=0)
    {
        historyData["binaryheader"]=binaryHeader;
    }
    if (missing != m_missing)
    {
        std::stringstream strs;
        strs << missing;
        historyData["missingvalue"]= strs.str();
    }
    if (text != m_text)
    {
        std::stringstream strs;
        strs << text;

        historyData["textvalue"]= strs.str();
    }
    if (strcmp(endian,"little")!=0)
    {
        historyData["bigendian"]= "";
    }
    if (strcmp(type,"float")!=0)
    {
        historyData["double"]= "";
    }
    if (points!=0)
    {
        std::stringstream strs;
        strs << points;

        historyData["npoints"]= strs.str();
    }
    if (strcmp(vo,"")!=0)
    {
        historyData["vo"]= vo;
    }
    if (strcmp(se,"")!=0)
    {
        historyData["se"]= se;
    }
    if (strcmp(lfnout,"")!=0)
    {
        historyData["lfnout"]= lfnout;
    }
    
    std::vector <std::string> outFilename;
    outFilename.insert(outFilename.begin(), historyData["out"].c_str());

    
    WriteXml(histFile,historyData,outFilename);
    
}

HistoryXmlWriter::HistoryXmlWriter(const char* histFile,const char* opName, std::map<std::string,std::string> appParameter,  std::vector <std::string> outFilename)
{
    appParameter["command"]="VisIVOFilters";
    appParameter["op"]=opName;
   // dump_map(appParameter);

    WriteXml(histFile,appParameter,outFilename);

}

HistoryXmlWriter::HistoryXmlWriter(const char* histFile, std::map<std::string,std::string> viewParameter,std::vector <std::string> outFilename)
{
    viewParameter["command"]="VisIVOViewer";
    // dump_map(viewParameter);
    
    
    
    WriteXml(histFile,viewParameter,outFilename);
    
}




void HistoryXmlWriter::WriteXml(const char* histFile,std::map<std::string, std::string> h_data,std::vector <std::string> outFilename)
{
    //dump_map(h_data);
    
    TiXmlDocument doc;
    doc.LoadFile (histFile);
    TiXmlElement* root = doc.FirstChildElement( "history" );
    if (! root )
    {

       // TiXmlDocument doc;
        TiXmlDeclaration * decl = new TiXmlDeclaration( "1.0", "", "" );
        doc.LinkEndChild( decl );
        
        //  TiXmlElement * root = new TiXmlElement( "history" );
        root = new TiXmlElement( "history" );
        doc.LinkEndChild( root );
    }
    
    TiXmlElement * command = new TiXmlElement( "command" );
    root->LinkEndChild( command );
    command ->SetAttribute("value",(h_data.find("command") -> second).c_str() );
    
    TiXmlElement * input = new TiXmlElement( "input" );
    command->LinkEndChild( input );
    input ->SetAttribute("value",(h_data.find("file") -> second).c_str() );
    
    
    for (std::vector<std::string>::iterator it = outFilename.begin() ; it != outFilename.end(); ++it)
    {
        TiXmlElement * output = new TiXmlElement( "output" );
        command->LinkEndChild( output );
       // output ->SetAttribute("value",(h_data.find("out") -> second).c_str() );
        std::string f_name=*it;
        output ->SetAttribute("value", f_name.c_str());
    }
   
    
    
    for ( std::map<std::string,std::string>::const_iterator it = h_data.begin(); it != h_data.end(); it++) {
        if(strcmp(it->first.c_str(), "history")== 0)
        {
            
        }
        else if(strcmp(it->first.c_str(), "command")!= 0  && strcmp(it->first.c_str(), "file")!= 0 && strcmp(it->first.c_str(), "out")!= 0   )
        {
            TiXmlElement * param = new TiXmlElement( "param" );
            command->LinkEndChild( param );
            param ->SetAttribute("name", it->first.c_str());
            param ->SetAttribute("value", it->second.c_str());
            
          //  std::cout<<"******************* "<< it->first.c_str()<<" - "<<it->second.c_str()<<" **********************"<<std::endl;

        }
        
    }
    
    doc.SaveFile( histFile );

}

void HistoryXmlWriter::dump_map(const std::map<std::string, std::string>& map) {
    for ( std::map<std::string,std::string>::const_iterator it = map.begin(); it != map.end(); it++) {
        std::cout << "Key: " << it->first << std::endl;
        std::cout << "Value: " << it->second <<std::endl;
    }
}