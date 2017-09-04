/*
 * Persistent Systems Private Limited, Pune, India (website - http://www.pspl.co.in)
 *
 * Permission is hereby granted, without written agreement and without
 * license or royalty fees, to use, copy, modify, and distribute this
 * software and its documentation for any purpose, provided that this
 * copyright notice and the following paragraph appears in all
 * copies of this software.
 *
 * DISCLAIMER OF WARRANTY.
 * -----------------------
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND.
 * IN NO EVENT SHALL THE RIGHTHOLDER BE LIABLE FOR ANY CLAIM, DAMAGES OR
 * OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 *
 */

/*
* This class represents the voFieldParam in the MetaData.  
* 
* A field contains description, a max of 2 values and a
* list of links.
*
* Date created - 03 May 2002
*
*/

//date:					20/09/2005
//modified by:	Marco Comparato (mcomp@oact.inaf.it)
//description:	modified to be used as class for fields
#include <cstdlib>
#include <cstring>
#include "voFieldParam.h"
//#include "VOUtils.h"

voFieldParam::voFieldParam()
{
	m_description = NULL;

	m_ID = NULL;
	m_unit = NULL;
	m_precision = NULL;
	m_width = 0;
	m_ref = NULL;
	m_utype = NULL;
	m_name = NULL;
	m_UCD = NULL;
	m_arraySize = NULL;
	m_datatype = 0;
}


void voFieldParam::init()
{
	m_description = NULL;

	m_ID = NULL;
	m_unit = NULL;
	m_precision = NULL;
	m_width = 0;
	m_ref = NULL;
	m_utype = NULL;
	m_name = NULL;
	m_UCD = NULL;
	m_arraySize = NULL;
	
}

void voFieldParam::makecopy(const voFieldParam &f)
{

	m_datatype = f.m_datatype;
	m_width = f.m_width;
//	m_values = f.m_values;
//	m_linkList = f.m_linkList;


	//try 
	{
		if (NULL != f.m_ID)
		{
			m_ID = new char[strlen(f.m_ID) + 1];
			strcpy(m_ID, f.m_ID);
		}
		if (NULL != f.m_utype)
		{
			m_utype = new char[strlen(f.m_utype) + 1];
			strcpy(m_utype, f.m_utype);
		}

		if (NULL != f.m_unit)
		{
			m_unit = new char[strlen(f.m_unit) + 1];
			strcpy(m_unit, f.m_unit);
		}

		if (NULL != f.m_name)
		{
			m_name = new char[strlen(f.m_name) + 1];
			strcpy(m_name, f.m_name);
		}

		if (NULL != f.m_precision)
		{
			m_precision = new char[strlen(f.m_precision) + 1];
			strcpy(m_precision, f.m_precision);
		}

		if (NULL != f.m_ref)
		{
			m_ref = new char[strlen(f.m_ref) + 1];
			strcpy(m_ref, f.m_ref);
		}

		if (NULL != f.m_arraySize)
		{
			m_arraySize = new char[strlen(f.m_arraySize) + 1];
			strcpy(m_arraySize, f.m_arraySize);
		}

		if (NULL != f.m_UCD)
		{
			m_UCD = new char[strlen(f.m_UCD) + 1];
			strcpy(m_UCD, f.m_UCD);
		}

		if (NULL != f.m_description)
		{
			m_description = new char[strlen(f.m_description) + 1];
			strcpy(m_description, f.m_description);
		}
	
	} 
	//catch (bad_alloc ex) 
	//{
		// Ignore ?
	//}


}

void voFieldParam::cleanup(void)
{
	delete[] m_ID;
	delete[] m_utype;
	delete[] m_unit;
	delete[] m_precision;
	delete[] m_ref;
	delete[] m_name;
	delete[] m_UCD;
	delete[] m_arraySize;
	delete[] m_description;

}

//*
//* To be implemented in the subclasses.
//*
//int voFieldParam::setValue(char *s, int *status)
//{
	// does nothing
	// do not expect this to be called.
//	return 0;
//}

//*
//* To be implemented in the subclasses.
//*
//int voFieldParam::setType(field_type type, int *status)
//{
	// does nothing
	// do not expect this to be called.
//	return 0;
//}

//modified by MComp. 03-11-2005
int voFieldParam::setDescription(char * desc, int *status)
{
  free(m_description);
  m_description = strdup(desc);

  //memory leaks???

//	delete[] m_description;
//	VOUtils::copyString(m_description, desc, status);

	return 0;
}

int voFieldParam::setID(char * ID, int *status)
{
  free(m_ID);
  m_ID = strdup(ID);

//	delete[] m_ID;
//	VOUtils::copyString(m_ID, ID, status);

	return 0;
}

int voFieldParam::setUtype(char * utype, int *status)
{
  free(m_utype);
  m_utype = strdup(utype);

//	delete[] m_utype;
//	VOUtils::copyString(m_utype, utype, status);

	return 0;
}

int voFieldParam::setUnit(char * unit, int *status)
{
  free(m_unit);
  m_unit = strdup(unit);

//	delete[] m_unit;
//	VOUtils::copyString(m_unit, unit, status);

	return 0;
}

//modified by MComp. 20-9-2005
int voFieldParam::setDatatype(char *datatype, int *status)
{
  free(m_datatype);
  m_datatype = strdup(datatype);
	
  return 0;
}

//modified by MComp. 03-11-2005
int voFieldParam::setPrecision(char * precision, int *status)
{
  free(m_precision);
  m_precision = strdup(precision);

//	delete[] m_precision;
//	VOUtils::copyString(m_precision, precision, status);

	return 0;
}

int voFieldParam::setWidth(int width, int *status)
{
	m_width = width;
	return 0;
}

int voFieldParam::setRef(char * ref, int *status)
{
  free(m_ref);
  m_ref = strdup(ref);

//	delete[] m_ref;
//	VOUtils::copyString(m_ref, ref, status);

	return 0;
}

int voFieldParam::setName(char * name, int *status)
{
  free(m_name);
  m_name = strdup(name);

//	delete[] m_name;
//	VOUtils::copyString(m_name, name, status);

  return 0;
}

int voFieldParam::setUCD(char * ucd, int *status)
{
  free(m_UCD);
  m_UCD = strdup(ucd);

//	delete[] m_UCD;
//	VOUtils::copyString(m_UCD, ucd, status);

	return 0;
}

int voFieldParam::setArraySize(char * arraySize, int *status)
{
	free(m_arraySize);
	m_arraySize = strdup(arraySize);

	return 0;
}



//int voFieldParam::setLinks(vector<Link> link, int *status)
//{
//	m_linkList = link;
//	return 0;
//}

//int voFieldParam::setValues(Values v[], int numOfValues, int *status)
//{
//	for (int i = 0; i < numOfValues && i < 2; i++)
//	{
//		m_values.push_back (v[i]);
//	}
//	return 0;
//}


//int voFieldParam::replaceValues(Values v, int *status)
//{
//	m_values.pop_back();
//	m_values.push_back(v);
//	return 0;
//}

/*
* Get description
*/ 
int voFieldParam::getDescription(char * &desc, int * status)
{
  desc = strdup(m_description);

//	VOUtils::copyString(desc, m_description, status);	

	return *status;
}

//*
//* Get 'Values', given an index.
//*
//int voFieldParam::getValues(Values &v, int index, int *status)
//{
//	*status = VOERROR;		
//	if (index >= 0 && index < m_values.size() )
//	{
		// assignment operator overloaded.
//		try
//		{
//			v = m_values[index];
//			*status = 0;
//		} 
//		catch (bad_alloc ex)
//		{
//			*status = INSUFFICIENT_MEMORY_ERROR;
//		}
//	} 

//	return *status;
//}

//*
//* Get number of values.
//*
//int voFieldParam::getNumOfValues (int &numOfValues, int *status)
//{
//	numOfValues = m_values.size();
//	*status = 0;	
//	return *status;
//}

//*
//* Get number of links.
//*
//int voFieldParam::getNumOfLinks(int  &nLinks, int *status)
//{
//	nLinks = m_linkList.size();
//	*status = 0;
//	return *status;
//}


//*
//* Get the link, given the index.
//*
//int voFieldParam::getLink(Link &link, int linkNum, int *status)
//{
//	*status = VOERROR;
//	if (linkNum >= 0 && linkNum < m_linkList.size()) 
//	{
//		try
//		{
//			link = m_linkList[linkNum];
//			*status = 0;
//		} 
//		catch (bad_alloc ex)
//		{
//			*status = INSUFFICIENT_MEMORY_ERROR;
//		}

//	}
//	return *status;
//}


int voFieldParam::getID(char * &ID, int *status)
{
  ID = strdup(m_ID);

//	VOUtils::copyString(ID, m_ID, status);

	return *status;
}

int voFieldParam::getUtype(char * &utype, int *status)
{
  utype = strdup(m_utype);

//  VOUtils::copyString(utype, m_utype, status);	

	return *status;
}

int voFieldParam::getUnit(char * &unit, int *status)
{
  unit = strdup(m_unit);
	
//  VOUtils::copyString(unit, m_unit, status);	

  return *status;
}

//modified by MComp. 20-9-2005
int voFieldParam::getDatatype(char * &datatype, int *status)
{
	datatype = strdup(m_datatype);
	*status = 0;
	return 0;
}

int voFieldParam::getPrecision(char * &precision, int *status)
{
  precision = strdup(m_precision);

//  VOUtils::copyString(precision, m_precision, status);	

  return *status;
}

int voFieldParam::getWidth(int &width, int *status)
{
	width = m_width;
	*status = 0;
	return 0;
}

int voFieldParam::getRef(char * &ref, int *status)
{
  ref = strdup(m_ref);

//  VOUtils::copyString(ref, m_ref, status);

  return *status;
}

int voFieldParam::getName(char * &name, int *status)
{
  name = strdup(m_name);
	
//  VOUtils::copyString(name, m_name, status);	

  return *status;
}

int voFieldParam::getUCD(char * &ucd, int *status)
{
  ucd = strdup(m_UCD);

//  VOUtils::copyString(ucd, m_UCD, status);

  return *status;
}

int voFieldParam::getArraySize(char * &arraySize, int *status)
{
  arraySize = strdup(m_arraySize);

//	VOUtils::copyString(arraySize, m_arraySize, status);

  return *status;
}

int voFieldParam::isVariableType(bool &b, int *status)
{
	b = false;
	char c = '*';

	if (m_arraySize != NULL && strchr(m_arraySize, c) != NULL) 
	{
		b = true;
	}

	*status = 0;
	return 0;
}
