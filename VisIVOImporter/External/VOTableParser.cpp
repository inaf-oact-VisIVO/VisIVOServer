/*
 * VisIVO, multidimensional visualization software
 * Copyright (C) 2007 Marco Comparato
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 */

/////////////////////////////////////////////////////////////////////////////
// Name:        VOTableParser.cpp
//
// Author:      Marco Comparato
// Date:        29/03/2007
/////////////////////////////////////////////////////////////////////////////
#include <cstdlib>
#include <cstring>

#include "VOTableParser.h"
#include "string.h"
#include <iostream>
#include <sstream>
#include <fstream>

#include "votable1.1.hxx"

#include <xercesc/dom/DOM.hpp>
#include <xercesc/framework/Wrapper4InputSource.hpp>

#include <xsd/cxx/xml/string.hxx>
#include <xsd/cxx/xml/dom/elements.hxx>

#include <xsd/cxx/xml/dom/serialization.hxx>
#include <xsd/cxx/xml/dom/bits/error-handler-proxy.hxx>

#include <xsd/cxx/tree/exceptions.hxx>
#include <xsd/cxx/tree/error-handler.hxx>

using namespace std;

//-------------------------------------------------------------------------
VOTableParser::VOTableParser()
//-------------------------------------------------------------------------
{
  m_source = "";
  m_outFile = "";

  m_resourceCount = 0;
  m_currentResource = 0;
  m_currentTable = 0;
  m_currentRow = 0;

  m_dataType = UNKNOWN;;
}

//-------------------------------------------------------------------------
VOTableParser::~VOTableParser()
//-------------------------------------------------------------------------
{
}

//------------------------------------------------------------------------
int VOTableParser::Parse(string source)
//-------------------------------------------------------------------------
{
  m_source = source;

  try
  {
    auto_ptr<votable::VOTABLE> vot(votable::VOTABLE_(m_source, xml_schema::flags::dont_validate));
    m_vot = vot;
  }
  catch(xsd::cxx::tree::parsing<char> e)
  {
    cerr << e.what() << endl;
    return 1;
  }
  catch(xsd::cxx::tree::unexpected_element<char> e)
  {
    cerr << e.what() << endl;
    return 1;
  }
  catch(xsd::cxx::tree::expected_attribute<char> e)
  {
    cerr << e.what() << endl;
    return 1;
  }
  catch(xsd::cxx::tree::unexpected_enumerator<char> e)
  {
    cerr << e.what() << endl;
    return 1;
  }

  m_resourceCount = m_vot->RESOURCE().size();

  if(m_resourceCount<1)
	  return 1;
  else
  {

  res      = m_vot->RESOURCE().begin();

  if(res->TABLE().size() > 0)
  {
    tab      = res->TABLE().begin();
    dat      = tab->DATA();
    tabData  = dat->TABLEDATA();
  }

  return 0;
}
}
//----------------------------------------------------------------------
int VOTableParser::Write(char *outFile, float *ddata, int nRows, int nCols)
//----------------------------------------------------------------------
{
  int i = 0;
  int j = 0;

  if(outFile)
  {
    cerr << "Writing..." << endl;

    SetOutFile(outFile);

    if(m_fieldNames.size() < nCols)
    {
      m_fieldNames.clear();

      for(i = 0; i < nCols; i++)
      {
        stringstream nameStream;
        nameStream << "field" << i;
        m_fieldNames.push_back(nameStream.str());
      }
    }

    using namespace xercesc;
    namespace xml = xsd::cxx::xml;

    XMLPlatformUtils::Initialize ();

    {
      xml::dom::auto_ptr<DOMDocument> doc(
        createEmptyDOMDocument("VOTABLE",
        "http://www.ivoa.net/xml/VOTable/v1.1",
        ""));

      DOMElement* root(doc->getDocumentElement ());

      root->setAttributeNS(
        xml::string ("http://www.w3.org/2001/XMLSchema-instance").c_str (),
        xml::string ("xsi:schemaLocation").c_str (),
        xml::string ("http://www.ivoa.net/xml/VOTable/v1.1 http://www.ivoa.net/xml/VOTable/v1.1").c_str ());

      auto_ptr<votable::VOTABLE> vot(votable::VOTABLE_(*doc, xml_schema::flags::dont_validate));

      votable::VOTABLE::RESOURCE::container resCont;
      votable::RESOURCE::TABLE::container tabCont;
      votable::TABLE::FIELD::container fieldCont;
      votable::TABLE::DATA::container dataCont;
      votable::DATA::TABLEDATA::container tableDataCont;
      votable::TABLEDATA::TR::container trCont;
      votable::RESOURCE res;
      votable::TABLE tab;
      votable::DATA data;
      votable::TABLEDATA tableData;

      for(i = 0; i < nCols; i++)
      {
        votable::FIELD::datatype::type fieldDataType("float");
        votable::FIELD field(fieldDataType, m_fieldNames[i].c_str());
        fieldCont.push_back(field);
      }

      for(i = 0; i < nRows; i++)
      {
        votable::TR tr;
        votable::TR::TD::container tdCont;

        for(j = 0; j < nCols; j++)
        {
          stringstream sValue;
          sValue << ddata[i * nCols + j]; //sValue << (((double) rand() / RAND_MAX) * 50) - 25;

          votable::TD td(sValue.str());
          tdCont.push_back(td);
        }

        tr.TD(tdCont);
        trCont.push_back(tr);
      }

      tableData.TR(trCont);

      tableDataCont.set(tableData);
      data.TABLEDATA(tableDataCont);

      dataCont.set(data);

      tab.FIELD(fieldCont);
      tab.DATA(dataCont);

      tabCont.push_back(tab);
      res.TABLE(tabCont);

      resCont.push_back(res);

      vot->RESOURCE(resCont);

      std::ofstream ofs(outFile);

      votable::VOTABLE_(*doc, *vot, xml_schema::flags::dont_validate);

      serializeDOMDocument(ofs, *doc);
    }

    XMLPlatformUtils::Terminate ();

    return 0;
  }

  return 1;
}

//-------------------------------------------------------------------------
string VOTableParser::GetTableDescription(int resource, int table)
//-------------------------------------------------------------------------
{
  SetCurrentResource(resource);

  return GetTableDescription(table);
}

//-------------------------------------------------------------------------
string VOTableParser::GetTableDescription(int table)
//-------------------------------------------------------------------------
{
//  tab = res->TABLE().begin() + table;
//  
//  return tab->DESCRIPTION().get();
  
  return "";
}

//-------------------------------------------------------------------------
int VOTableParser::GetResourceTRCount(int resource, int table)
//-------------------------------------------------------------------------
{
  SetCurrentResource(resource);

  return GetResourceTRCount(table);
}


//-------------------------------------------------------------------------
int VOTableParser::GetResourceTRCount(int table)
//-------------------------------------------------------------------------
{
  if(m_vot.get())
  {
    tab = res->TABLE().begin() + table;
    dat = tab->DATA();
    tabData = dat->TABLEDATA();

    return tabData->TR().size();
  }

  return 0;
}
  
//-------------------------------------------------------------------------
int VOTableParser::GetResourceTableFieldsCount(int resource, int table)
//-------------------------------------------------------------------------
{
  SetCurrentResource(resource);
  SetCurrentTable(table);

  return GetResourceTableFieldsCount();
}

//-------------------------------------------------------------------------
int VOTableParser::GetResourceTableFieldsCount(int table)
//-------------------------------------------------------------------------
{
  SetCurrentTable(table);

  return GetResourceTableFieldsCount();
}

//-------------------------------------------------------------------------
int VOTableParser::GetResourceTableFieldsCount()
//-------------------------------------------------------------------------
{
  if(m_vot.get())
    return tab->FIELD().size();

  return 0;
}


//-------------------------------------------------------------------------
string VOTableParser::GetResourceTableFieldName(int resource, int table, int field)
//-------------------------------------------------------------------------
{
  SetCurrentResource(resource);
  SetCurrentTable(table);

  return GetResourceTableFieldName(field);
}

//-------------------------------------------------------------------------
string VOTableParser::GetResourceTableFieldName(int table, int field)
//-------------------------------------------------------------------------
{
  SetCurrentTable(table);

  return GetResourceTableFieldName(field);
}

//-------------------------------------------------------------------------
string VOTableParser::GetResourceTableFieldName(int field)
//-------------------------------------------------------------------------
{
  votable::TABLE::FIELD::iterator tabField = tab->FIELD().begin() + field;

  return tabField->name();
}


//-------------------------------------------------------------------------
string VOTableParser::GetResourceTableFieldUcd(int resource, int table, int field)
//-------------------------------------------------------------------------
{
  SetCurrentResource(resource);
  SetCurrentTable(table);

  return GetResourceTableFieldUcd(field);
}

//-------------------------------------------------------------------------
string VOTableParser::GetResourceTableFieldUcd(int table, int field)
//-------------------------------------------------------------------------
{
  SetCurrentTable(table);

  return GetResourceTableFieldUcd(field);
}

//-------------------------------------------------------------------------
string VOTableParser::GetResourceTableFieldUcd(int field)
//-------------------------------------------------------------------------
{
  votable::TABLE::FIELD::iterator tabField = tab->FIELD().begin() + field;

  if(tabField->ucd())
    return tabField->ucd().get();

  return "";
}


//-------------------------------------------------------------------------
string VOTableParser::GetResourceTableFieldDataType(int resource, int table, int field)
//-------------------------------------------------------------------------
{
  SetCurrentResource(resource);
  SetCurrentTable(table);

  return GetResourceTableFieldDataType(field);
}

//-------------------------------------------------------------------------
string VOTableParser::GetResourceTableFieldDataType(int table, int field)
//-------------------------------------------------------------------------
{
  SetCurrentTable(table);

  return GetResourceTableFieldDataType(field);
}

//-------------------------------------------------------------------------
string VOTableParser::GetResourceTableFieldDataType(int field)
//-------------------------------------------------------------------------
{
  votable::TABLE::FIELD::iterator tabField = tab->FIELD().begin() + field;

  return tabField->datatype();
}


//-------------------------------------------------------------------------
string VOTableParser::GetResourceTableFieldId(int resource, int table, int field)
//-------------------------------------------------------------------------
{
  SetCurrentResource(resource);
  SetCurrentTable(table);

  return GetResourceTableFieldId(field);
}

//-------------------------------------------------------------------------
string VOTableParser::GetResourceTableFieldId(int table, int field)
//-------------------------------------------------------------------------
{
  SetCurrentTable(table);

  return GetResourceTableFieldId(field);
}

//-------------------------------------------------------------------------
string VOTableParser::GetResourceTableFieldId(int field)
//-------------------------------------------------------------------------
{
  votable::TABLE::FIELD::iterator tabField = tab->FIELD().begin() + field;

  if(tabField->ID())
    return tabField->ID().get();

  return "";
}


//-------------------------------------------------------------------------
string VOTableParser::GetResourceTableFieldArraySize(int resource, int table, int field)
//-------------------------------------------------------------------------
{
  SetCurrentResource(resource);
  SetCurrentTable(table);

  return GetResourceTableFieldArraySize(field);
}

//-------------------------------------------------------------------------
string VOTableParser::GetResourceTableFieldArraySize(int field)
//-------------------------------------------------------------------------
{
  votable::TABLE::FIELD::iterator tabField = tab->FIELD().begin() + field;

  if(tabField->arraysize())
    return tabField->arraysize().get();

  return "";
}

//-------------------------------------------------------------------------
string VOTableParser::GetResourceTableFieldRef(int resource, int table, int field)
//-------------------------------------------------------------------------
{
  SetCurrentResource(resource);
  SetCurrentTable(table);

  return GetResourceTableFieldRef(field);
}

//-------------------------------------------------------------------------
string VOTableParser::GetResourceTableFieldRef(int field)
//-------------------------------------------------------------------------
{
  votable::TABLE::FIELD::iterator tabField = tab->FIELD().begin() + field;

  if(tabField->ref())
    return tabField->ref().get();

  return "";
}


//-------------------------------------------------------------------------
string VOTableParser::GetResourceTableFieldUnit(int resource, int table, int field)
//-------------------------------------------------------------------------
{
  SetCurrentResource(resource);
  SetCurrentTable(table);

  return GetResourceTableFieldUnit(field);
}

//-------------------------------------------------------------------------
string VOTableParser::GetResourceTableFieldUnit(int field)
//-------------------------------------------------------------------------
{
  votable::TABLE::FIELD::iterator tabField = tab->FIELD().begin() + field;

  if(tabField->unit())
    return tabField->unit().get();

  return "";
}


//-------------------------------------------------------------------------
string VOTableParser::GetResourceTableFieldPrecision(int resource, int table, int field)
//-------------------------------------------------------------------------
{
  SetCurrentResource(resource);
  SetCurrentTable(table);

  return GetResourceTableFieldPrecision(field);
}

//-------------------------------------------------------------------------
string VOTableParser::GetResourceTableFieldPrecision(int field)
//-------------------------------------------------------------------------
{
  votable::TABLE::FIELD::iterator tabField = tab->FIELD().begin() + field;

  if(tabField->precision())
    return tabField->precision().get();

  return "";
}

//-------------------------------------------------------------------------
string VOTableParser::GetResourceTableFieldDescription(int resource, int table, int field)
//-------------------------------------------------------------------------
{
  SetCurrentResource(resource);
  SetCurrentTable(table);

  return GetResourceTableFieldDescription(field);
}

//-------------------------------------------------------------------------
string VOTableParser::GetResourceTableFieldDescription(int field)
//-------------------------------------------------------------------------
{
  //TODO
  return "";
}


//-------------------------------------------------------------------------
unsigned int VOTableParser::GetResourceTableFieldWidth(int resource, int table, int field)
//-------------------------------------------------------------------------
{
  SetCurrentResource(resource);
  SetCurrentTable(table);

  return GetResourceTableFieldWidth(field);
}

//-------------------------------------------------------------------------
unsigned int VOTableParser::GetResourceTableFieldWidth(int field)
//-------------------------------------------------------------------------
{
  votable::TABLE::FIELD::iterator tabField = tab->FIELD().begin() + field;

  if(tabField->width())
    return tabField->width().get();

  return 0;
}

//-------------------------------------------------------------------------
string VOTableParser::GetVOTableDescription()
//-------------------------------------------------------------------------
{
  //TODO
  return "";
}

//-------------------------------------------------------------------------
string VOTableParser::GetResourceName(int resource)
//-------------------------------------------------------------------------
{
  SetCurrentResource(resource);

  return GetResourceName();
}

//-------------------------------------------------------------------------
string VOTableParser::GetResourceName()
//-------------------------------------------------------------------------
{
  if(res->name())
    return res->name().get();

  return "";
}

//-------------------------------------------------------------------------
string VOTableParser::GetResourceDescription(int resource)
//-------------------------------------------------------------------------
{
  SetCurrentResource(resource);

  return GetResourceDescription();
}

//-------------------------------------------------------------------------
string VOTableParser::GetResourceDescription()
//-------------------------------------------------------------------------
{
  //TODO

  return "";
}



//-------------------------------------------------------------------------
int VOTableParser::GetResourceCount()
//-------------------------------------------------------------------------
{
  return m_resourceCount;
}

//-------------------------------------------------------------------------
int VOTableParser::GetResourceTableCount(int resource)
//-------------------------------------------------------------------------
{
  SetCurrentResource(resource);

  return GetResourceTableCount();
}

//-------------------------------------------------------------------------
int VOTableParser::GetResourceTableCount()
//-------------------------------------------------------------------------
{
    return res->TABLE().size();
}

//-------------------------------------------------------------------------
int VOTableParser::GetResourceParamsCount(int resource)
//-------------------------------------------------------------------------
{
  SetCurrentResource(resource);

  return GetResourceParamsCount();
}

//-------------------------------------------------------------------------
int VOTableParser::GetResourceParamsCount()
//-------------------------------------------------------------------------
{
  return res->PARAM().size();
}

//-------------------------------------------------------------------------
int VOTableParser::GetResourceParamVOptionCount(int resource, int param)
//-------------------------------------------------------------------------
{
  SetCurrentResource(resource);

  return GetResourceParamVOptionCount(param);
}

//-------------------------------------------------------------------------
int VOTableParser::GetResourceParamVOptionCount(int param)
//-------------------------------------------------------------------------
{
  if(res->PARAM().size() > 0)
  {
    votable::RESOURCE::PARAM::iterator resParam = res->PARAM().begin() + param;
    votable::PARAM::VALUES::container parValue = resParam->VALUES();
 
    return parValue->OPTION().size();
  }

  return 0;
}


//-------------------------------------------------------------------------
string VOTableParser::GetResourceParamVOption(int resource, int param, int option)
//-------------------------------------------------------------------------
{
  SetCurrentResource(resource);

  return GetResourceParamVOption(param, option);
}

//-------------------------------------------------------------------------
string VOTableParser::GetResourceParamVOption(int param, int option)
//-------------------------------------------------------------------------
{
  votable::RESOURCE::PARAM::iterator resParam = res->PARAM().begin() + param;
  votable::PARAM::VALUES::container parValue = resParam->VALUES();

  if(parValue->OPTION().size())
    return (parValue->OPTION().begin() + option)->value();

  return "";
}

//-------------------------------------------------------------------------
string VOTableParser::GetData(int resource, int table, int i, int j)
//-------------------------------------------------------------------------
{
  SetCurrentResource(resource);
  //SetCurrentTRSet(trSet);
  SetCurrentTable(table);
  SetCurrentTR(i);

  return GetData(j);
}

//-------------------------------------------------------------------------
string VOTableParser::GetData(int tdSet)
//-------------------------------------------------------------------------
{
  votable::TABLEDATA::TR::iterator row = tabData->TR().begin() + m_currentRow;

  if(row->TD().size())
    return (row->TD())[tdSet];

  return "";
}

//-------------------------------------------------------------------------
void VOTableParser::SetCurrentResource(int resource)
//-------------------------------------------------------------------------
{
  m_currentResource = resource;
  m_currentTable = 0;
  m_currentRow = 0;

  if(m_vot.get())
  {
    res = m_vot->RESOURCE().begin() + m_currentResource;

    if(res->TABLE().size() > 0)
    {
      tab = res->TABLE().begin();
      dat = tab->DATA();
      tabData = dat->TABLEDATA();
    }
  }

  return;
}

//-------------------------------------------------------------------------
void VOTableParser::SetCurrentTable(int table)
//-------------------------------------------------------------------------
{
  m_currentTable = table;
  m_currentRow = 0;

  if(m_vot.get())
  {
    res = m_vot->RESOURCE().begin() + m_currentResource;      
    tab = res->TABLE().begin() + table;
    dat = tab->DATA();

    if(tabData = dat->TABLEDATA())
    {
      m_dataType = VOTABLE_TABLEDATA;
    }
    else if(dat->BINARY())
    {
      m_dataType = VOTABLE_BINARY;
    }
    else
    {
      m_dataType = UNKNOWN;
    }
  }

  return;
}

//-------------------------------------------------------------------------
void VOTableParser::SetCurrentTR(int i)
//-------------------------------------------------------------------------
{
  m_currentRow = i;

  return;
}

//-------------------------------------------------------------------------
int VOTableParser::GetDataType()
//-------------------------------------------------------------------------
{
  return m_dataType;
}

//-------------------------------------------------------------------------
xsd::cxx::xml::dom::auto_ptr<xercesc::DOMDocument>
VOTableParser::createEmptyDOMDocument(const std::string& name,
                                      const std::string& ns,
                                      const std::string& prefix)
//-------------------------------------------------------------------------
{
  using namespace xercesc;
  namespace xml = xsd::cxx::xml;

  const XMLCh ls_id [] = {chLatin_L, chLatin_S, chNull};

  // Get an implementation of the Load-Store (LS) interface.
  //
  DOMImplementation* impl (
    DOMImplementationRegistry::getDOMImplementation (ls_id));

  xml::dom::auto_ptr<DOMDocument> doc (
    impl->createDocument (
    (ns.empty () ? 0 : xml::string (ns).c_str ()),
    xml::string ((prefix.empty () ? name : prefix + ':' + name)).c_str (),
    0));

  return doc;
}

//-------------------------------------------------------------------------
void
VOTableParser::serializeDOMDocument(std::ostream& os,
           const xercesc::DOMDocument& doc,
           const std::string& encoding /*= "UTF-8"*/)
//-------------------------------------------------------------------------
{
  using namespace xercesc;
  namespace xml = xsd::cxx::xml;
  namespace tree = xsd::cxx::tree;

  const XMLCh ls_id [] = {chLatin_L, chLatin_S, chNull};

  // Get an implementation of the Load-Store (LS) interface.
  //
  DOMImplementation* impl (
    DOMImplementationRegistry::getDOMImplementation (ls_id));

  // Create a DOMWriter.
  //
  xml::dom::auto_ptr<DOMWriter> writer (impl->createDOMWriter ());

  // Set error handler.
  //
  tree::error_handler<char> eh;
  xml::dom::bits::error_handler_proxy<char> ehp (eh);
  writer->setErrorHandler (&ehp);

  // Set encoding.
  //
  writer->setEncoding(xml::string (encoding).c_str ());

  // Set some generally nice features.
  //
  writer->setFeature (XMLUni::fgDOMWRTDiscardDefaultContent, true);
  writer->setFeature (XMLUni::fgDOMWRTFormatPrettyPrint, true);

  // Adapt ostream to format target and serialize.
  //
  xml::dom::ostream_format_target oft (os);

  writer->writeNode (&oft, doc);

  eh.throw_if_failed<tree::parsing<char> > ();
}
