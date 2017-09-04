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

/*=========================================================================
  Language:			C++
  Date:				  29/03/2007
  Authors:			Marco Comparato mcomp@ct.astro.it
  Description:  class to access the VOTables
=========================================================================*/
 
#ifndef __VOTableParser_h
#define __VOTableParser_h

#include <vector>
#include <string>
#include <memory>
#include "votable1.1.hxx"

enum
{
  VOTABLE_TABLEDATA,
  VOTABLE_BINARY,
  UNKNOWN
};



class VOTableParser
{
  std::auto_ptr<votable::VOTABLE> m_vot;

  std::string m_source;
  std::string m_outFile;

  std::vector<std::string> m_fieldNames;

  int m_resourceCount;
  int m_currentResource;
  int m_currentTable;
  int m_currentRow;

  int m_dataType;

  votable::VOTABLE::RESOURCE::iterator res;
  votable::RESOURCE::TABLE::iterator tab;
  votable::TABLE::DATA::container dat;
  votable::DATA::TABLEDATA::container tabData;

  void SetSource(std::string source) { m_source = source; };
  void SetOutFile(std::string outFile) { m_outFile = outFile; };

  xsd::cxx::xml::dom::auto_ptr<xercesc::DOMDocument>
    createEmptyDOMDocument(const std::string& root_element_name,
    const std::string& root_element_namespace = "",
    const std::string& root_element_namespace_prefix = "");

  void
    serializeDOMDocument(std::ostream& os,
    const xercesc::DOMDocument& doc,
    const std::string& encoding = "UTF-8");

public:

  VOTableParser();
  ~VOTableParser();

  int Init() { return 0; };
  int Parse(std::string source);
  int Write(char *outFile, float *ddata, int nRows, int nCols);

  std::string GetResourceName(int resource);
  std::string GetResourceName();
  std::string GetResourceDescription(int resource);
  std::string GetResourceDescription();

  int GetResourceCount();
  int GetResourceTableCount(int resource);
  int GetResourceTableCount();
  int GetResourceParamsCount(int resource);
  int GetResourceParamsCount();
  int GetResourceParamVOptionCount(int resource, int param);
  int GetResourceParamVOptionCount(int param);
  int GetResourceTRCount(int resource, int table);
  int GetResourceTRCount(int table);
  int GetResourceTableFieldsCount(int resource, int table);
  int GetResourceTableFieldsCount(int table);
  int GetResourceTableFieldsCount();

  std::string  GetResourceTableFieldName(int resource, int table, int field);
  std::string  GetResourceTableFieldName(int table, int field);
  std::string  GetResourceTableFieldName(int field);
  std::string  GetResourceTableFieldDataType(int resource, int table, int field);
  std::string  GetResourceTableFieldDataType(int table, int field);
  std::string  GetResourceTableFieldDataType(int field);
  std::string  GetResourceTableFieldUcd(int resource, int table, int field);
  std::string  GetResourceTableFieldUcd(int table, int field);
  std::string  GetResourceTableFieldUcd(int field);
  std::string  GetResourceTableFieldId(int resource, int table, int field);
  std::string  GetResourceTableFieldId(int table, int field);
  std::string  GetResourceTableFieldId(int field);
  std::string  GetResourceTableFieldArraySize(int resource, int table, int field);
  std::string  GetResourceTableFieldArraySize(int field);
  std::string  GetResourceTableFieldRef(int resource, int table, int field);
  std::string  GetResourceTableFieldRef(int field);
  std::string  GetResourceTableFieldUnit(int resource, int table, int field);
  std::string  GetResourceTableFieldUnit(int field);
  std::string  GetResourceTableFieldPrecision(int resource, int table, int field);
  std::string  GetResourceTableFieldPrecision(int field);
  std::string  GetResourceTableFieldDescription(int resource, int table, int field);
  std::string  GetResourceTableFieldDescription(int field);
  unsigned int GetResourceTableFieldWidth(int resource, int table, int field);
  unsigned int GetResourceTableFieldWidth(int field);

  std::string GetResourceParamVOption(int resource, int param, int option);
  std::string GetResourceParamVOption(int param, int option);
  
  std::string GetData(int resource, int table, int i, int j);
  std::string GetData(int tdSet);
  std::string GetVOTableDescription();
  std::string GetTableDescription(int resource, int table);
  std::string GetTableDescription(int table);
  
  int GetDataType();
  
  void SetCurrentResource(int resource);
  void SetCurrentTable(int table);
  void SetCurrentTR(int i);

  void SetFieldNames(std::vector<std::string> fieldNames) { m_fieldNames = fieldNames; };
};

#endif
