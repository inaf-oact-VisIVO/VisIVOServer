/*******************************************************************************************************
 *  ###################### README #####################################################    						 *
 * VARIABLE
 *       m_fArray[][]: array where table data are stored
 *       m_nOfCol: number of selected columns given as arguments in --field;
 *       m_nOfAllocatedCol: number of columns allocated to store data in new table.
 *      		                If --allcolumns is given m_nOfAllocatedCol=allcolumns of
 *                          original table, otherwise m_nOfAllocatedCol=m_nOfCol;
 *       m_nOfRow: number of rows             
 *       m_goodEventList[]: array of row index of original table corresponding to events
 *                        passing the n_sigma cut;
 *       m_nOfSigma: number of sigma required in the event selection;
 *       m_nOfEle: maximum number of rows of m_fArray that can be allocated
 *       colId[]: array with IDs of selected columns (e.g. --field X Z Vz ==> colId[3]={0,2,5});
 *       outcolId[]: array with IDs of columns to be stored in the new table
 *                   (e.g if --allcolumns is given outcolId[6]={0,1,2,3,4,5}, otherwise 
 *                    outcolId[3]={0,1,2});
 *       allcolId[]: array with IDs of allocated columns to be read from original table
 *                   (e.g if --allcolumns is given allcolId[]=outcolId[], otherwise 
 *                    allcolId[]=colId[]); 
 *       
 *  FUNCTION USAGE
 *      getColumn(colId,colNumber,StartRow,EndRow,Array): colId is the id of the columns to be read from the
 *                     original table, starting from row StartRow to EndRow (with respect to the original table).
 *                     These columns are stored in the array sequencially (e.g. if colId[3]={0,2,4} array[][0,1,2]
 *                             
 ********************************************************************************************************/


#include "vssigmacontoursop.h"
#include "vstable.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <algorithm>


const unsigned int VSSigmaContoursOp::MIN_NUMBER_OF_ROW = 10000;
const unsigned int VSSigmaContoursOp::MAX_NUMBER_TO_REDUCE_ROW = 10000;


VSSigmaContoursOp::VSSigmaContoursOp(){

	m_fArray=NULL;
	m_goodEventList=NULL;
  m_nOfCol= 0;
  m_nOfAllocatedCol= 0;
  m_nOfSigma= 1;
  m_exclude=false;
  m_nOfRow= 0;
  m_nOfEle= 0;
  
}//close constructor

VSSigmaContoursOp::~VSSigmaContoursOp(){
	
	
	if(m_fArray!=NULL) {
		for(int i=0;i<m_nOfAllocatedCol;i++){
			if(m_fArray[i]!=NULL) delete [] m_fArray[i];	
		}
		delete [] m_fArray;
	}
	
	if(m_goodEventList!=NULL) delete [] m_goodEventList;
	

}//close destructor

//--------------------------------------------------------------------
bool VSSigmaContoursOp::execute(){
//--------------------------------------------------------------------
	//#############################
	//####    PARSE ARGUMENTS  ####
	//#############################
	//nsigma
	if(!isParameterPresent("nsigma")) std::cout<<"No nsigma value is given...using default value (1 sigma)"<<std::endl;
	else if(getParameterAsString("nsigma").empty()||getParameterAsString("nsigma")=="unknown" ){
		std::cerr<<"VSSigmaContoursOp:: No sigma value is given"<<std::endl;
		return false;
	}
	else if(getParameterAsFloat("nsigma")<=0){
		std::cerr<<"VSSigmaContoursOp:: Negative or zero nsigma given...insert positive values!"<<std::endl;
		return false;		
	}
	else m_nOfSigma= getParameterAsFloat("nsigma");

	if(isParameterPresent("exclude") ) m_exclude=true;
	
	//fields
	bool fieldIsPresent= true;
	if(!isParameterPresent("field") ) fieldIsPresent= false;
	
	if((getParameterAsString("field").empty()||getParameterAsString("field")=="unknown")&&fieldIsPresent ){
		std::cerr<<"VSSigmaContoursOp:: No field is given"<<std::endl;
		return false;
	}
	
	std::stringstream ssField;
	
	if(fieldIsPresent){
		ssField << getParameterAsString("field");
	
		//Count all fields given as arguments of --field
		while(!ssField.eof()){
			std::string tmpfield;
			ssField>> tmpfield;
		
			unsigned int tmpId= m_tables[0]->getColId(tmpfield);;
			if(tmpId==-1){
				std::cerr<<"Invalid field name"<<tmpfield<<std::endl;
			} 
			else m_nOfCol++;		
		}//close while
	
	}//close if
	else m_nOfCol= m_tables[0]->getNumberOfColumns();
	if(m_nOfCol==0)
	{
		 std::cerr<<"No valid field names are given"<<std::endl;
		 return false;
	}
	
	//reset string buffer to initial position
	ssField.clear();
	ssField.seekg(0,std::ios_base::beg);

	//Store id of selected columns
	std::string * col;

	col = new std::string[m_nOfCol];

	try
	{
		col = new std::string[m_nOfCol];
	}
	catch(std::bad_alloc &e)
	{
		col=NULL;
	}

	if (col==NULL)
		return false;

	//unsigned int colId[m_nOfCol];//id of selected columns
	unsigned int * colId;
	try{
	colId = new unsigned int [m_nOfCol];
	}
	catch (std::bad_alloc &e) 
	{
	colId = NULL;
	}
	if (colId==NULL)
	{
		delete [] col;
		return false;
	}

	if(fieldIsPresent)
		for(int i=0;i<m_nOfCol;i++){
			ssField>>col[i];
			colId[i]= m_tables[0]->getColId(col[i]);
			std::cout<<"COL NAME= "<<col[i]<<" ==> Id="<<colId[i]<<std::endl;	
		}//close for
	else 
	 for(int i=0;i<m_nOfCol;i++){
			col[i]= m_tables[0]->getColName(i);
			colId[i]= i;
			std::cout<<"COL NAME= "<<col[i]<<" ==> Id="<<colId[i]<<std::endl;	
		}//close for
	
		
	//output file name
	std::string outName;
	if(!isParameterPresent("out")||getParameterAsString("out")=="unknown") outName="VSSigmaContoursOp.bin";
	else {
		outName= getParameterAsString("out");
		if(outName.find(".bin") == std::string::npos) outName.append(".bin");
	}
	
	m_realOutFilename.push_back(outName);

	//#############################
	//####   DEFINE NEW TABLE  ####
	//#############################
	//get number of table rows
	m_nOfRow= m_tables[0]->getNumberOfRows();
	
	VSTable newTable;
	newTable.setType("float");
	newTable.setIsVolume(m_tables[0]->getIsVolume());
	newTable.setCellSize(m_tables[0]->getCellSize()[0],m_tables[0]->getCellSize()[1],m_tables[0]->getCellSize()[2]);
	newTable.setCellNumber(m_tables[0]->getCellNumber()[0],m_tables[0]->getCellNumber()[1],m_tables[0]->getCellNumber()[2]);
	newTable.setEndiannes(m_tables[0]->getEndiannes());
	
	//Check if only the cut columns have to be stored to a new table 
	//otherwise write all columns
	unsigned int* allcolId;//id of allocated columns
	unsigned int* outcolId;//id of allocated columns to be stored in the new table
	allcolId=NULL;
	outcolId=NULL;
	
	if(!isParameterPresent("allcolumns")){	
		m_nOfAllocatedCol= m_nOfCol;
		allcolId= new unsigned int[m_nOfAllocatedCol];
		outcolId= new unsigned int[m_nOfAllocatedCol];
		for(int i=0;i<m_nOfAllocatedCol;i++){
			newTable.addCol(col[i]);
			allcolId[i]= colId[i];
			outcolId[i]= i;
		}//close for	
					
	}//close if
	else{
		unsigned int nColumns= m_tables[0]->getNumberOfColumns();
		m_nOfAllocatedCol= nColumns;
		allcolId= new unsigned int[m_nOfAllocatedCol];
		outcolId= new unsigned int[m_nOfAllocatedCol];
		for(unsigned int i=0;i<m_nOfAllocatedCol;i++){
			std::string colname= m_tables[0]->getColName(i);
			newTable.addCol(colname);			
			allcolId[i]= i;
			outcolId[i]= i;			
		}//close for
	}//close else
		
	newTable.setLocator(outName);
  
	
	if(m_nOfRow>getMaxNumberInt()) m_nOfEle= getMaxNumberInt();
	else m_nOfEle= m_nOfRow;
	
	std::cout<<"***************************************************************"<<std::endl;
	std::cout<<"*******     PARAMETERS                  ***********************"<<std::endl;
	std::cout<<"***************************************************************"<<std::endl;
	std::cout<<"Number of sigma ==> "<<m_nOfSigma<<std::endl;
	std::cout<<"Number of selected columns ==> "<<m_nOfCol<<std::endl;
	for(int i=0;i<m_nOfCol;i++) std::cout<<"                      col Id "<<colId[i]<<" ==> "<<col[i]<<std::endl;	
	std::cout<<"Number of stored columns ==> "<<m_nOfAllocatedCol<<std::endl;
	std::cout<<"New table stored in file ==> "<<outName<<std::endl;
	std::cout<<"***************************************************************"<<std::endl;
	
	//#############################
	//####    ALLOCATE ARRAY   ####
	//#############################	
	bool goodAllocation= allocateArray();
	
	if(!goodAllocation){
	  std::cerr<<"Invalid allocation"<<std::endl;
	  if(col!=NULL) delete [] col;
	  if(colId!=NULL) delete [] colId;
	  if(allcolId!=NULL) delete [] allcolId;
	  if(outcolId!=NULL) delete [] outcolId;	
	  return false;	
	} 	
	
	unsigned long long int index=0;
	unsigned long long int index_sel=0;
	float *mean, *sigma;
	try {
		mean = new float[m_nOfCol];
	}
	catch(std::bad_alloc &e)
	{
	mean = NULL;
	}
	if (mean == NULL)
	{
	  if(col!=NULL) delete [] col;
	  if(colId!=NULL) delete [] colId;
	  if(allcolId!=NULL) delete [] allcolId;
	  if(outcolId!=NULL) delete [] outcolId;	
	return false;
	}
	try {
		sigma = new float[m_nOfCol];
	}
	catch(std::bad_alloc &e)
	{
	sigma = NULL;
	}
	if (sigma == NULL)
	{
	  if(mean!=NULL) delete [] mean;
	  if(col!=NULL) delete [] col;
	  if(colId!=NULL) delete [] colId;
	  if(allcolId!=NULL) delete [] allcolId;
	  if(outcolId!=NULL) delete [] outcolId;	
	return false;
	}
	float sum= 0.;
	int m_nOfRowsToBeRead= m_nOfEle;
	
	//#############################
	//####    CALCULATE MEAN   ####
	//#############################
	std::cout<<"*** Calculate mean ***"<<std::endl;
	while(index< m_nOfRow){	
		m_tables[0]->getColumn(colId,m_nOfCol,index,index+m_nOfEle-1,m_fArray);
					
		for(int j=0;j<m_nOfCol;j++) {
			sum= 0.;//initialize to zero
			for(int i=0;i<m_nOfEle;i++) {			
				sum+= m_fArray[j][i];				
			}//close for i
			mean[j]+= sum;
		}//close for j
						
		index+= m_nOfEle;
		if(index+m_nOfEle>= m_nOfRow) m_nOfEle= m_nOfRow-index;
			
	}//close while	
	
	//divide for number of elements to get the mean
	for(int i=0;i<m_nOfCol;i++) mean[i]/= m_nOfRow;
	
	
	//#############################
	//####    CALCULATE SIGMA   ####
	//#############################
	index=0;
	m_nOfEle= m_nOfRowsToBeRead;
	std::cout<<"*** Calculate sigma ***"<<std::endl;
	while(index< m_nOfRow){
		m_tables[0]->getColumn(colId,m_nOfCol,index,index+m_nOfEle-1,m_fArray);				
		for(int j=0;j<m_nOfCol;j++) {
			sum= 0.;//initialize to zero
			for(int i=0;i<m_nOfEle;i++) {			
				sum+= (m_fArray[j][i]-mean[j])*(m_fArray[j][i]-mean[j]);				
			}//close for i
			sigma[j]+= sum;
		}//close for j
						
		index+= m_nOfEle;
		if(index+m_nOfEle>= m_nOfRow) m_nOfEle= m_nOfRow-index;
			
	}//close while
		
	for(int i=0;i<m_nOfCol;i++) sigma[i]= sqrt(sigma[i]/(m_nOfRow-1));
	
	//Print mean and sigma
	std::cout<<std::endl;
	for(int i=0;i<m_nOfCol;i++) std::cout<<std::setprecision(9)<<"variable "<<col[i]<<" MEAN= "<<mean[i]<<"  SIGMA="<<sigma[i]<<std::endl;
	std::cout<<std::endl;
		
	//#########################################
	//####    CALCULATE GOOD EVENT NUMBER  ####
	//#########################################
	unsigned long long int goodev_counter=0;
	index=0;
	m_nOfEle= m_nOfRowsToBeRead;
	std::cout<<"*** Selecting events ***"<<std::endl;
	while(index< m_nOfRow){
		m_tables[0]->getColumn(colId,m_nOfCol,index,index+m_nOfEle-1,m_fArray);				
		for(int i=0;i<m_nOfEle;i++)
		{
			int j;			
			for(j=0;j<m_nOfCol;j++) 
			{						
				if( (m_fArray[j][i]-mean[j])> m_nOfSigma*sigma[j]|| (m_fArray[j][i]-mean[j])< -m_nOfSigma*sigma[j])
				{ 
					if(!m_exclude) break;
				} else
				{
					if(m_exclude) break;

				}		
			}//close for j	
			if(j==m_nOfCol) goodev_counter++;		
		}//close i						
		index+= m_nOfEle;
		if(index+m_nOfEle>= m_nOfRow) m_nOfEle= m_nOfRow-index;
	}//close while

	std::cout<<"Number of selected events ==> "<<goodev_counter<<std::endl;

  if(goodev_counter==0) return true;
	newTable.setNumberOfRows(goodev_counter);
	newTable.writeHeader();
	
	
	//###################################
	//####    STORE EVENTS IN TABLE  ####
	//###################################
	index=0;
	m_nOfEle= m_nOfRowsToBeRead;
	std::cout<<"*** Storing events in new table ***"<<std::endl;
	while(index< m_nOfRow){
		m_tables[0]->getColumn(colId,m_nOfCol,index,index+m_nOfEle-1,m_fArray);
		
		int Nev_sel=0;								
		for(int i=0;i<m_nOfEle;i++)
		{		
			int j;
			for(j=0;j<m_nOfCol;j++) 
			{
					if( (m_fArray[j][i]-mean[j])> m_nOfSigma*sigma[j]|| (m_fArray[j][i]-mean[j])< -m_nOfSigma*sigma[j])
					{ 
						if(!m_exclude) break;
					} else
					{
						if(m_exclude) break;
					}	
			}//close for j	
			if(j==m_nOfCol)
			{ 
				m_goodEventList[Nev_sel]= index+i;
 				Nev_sel++;
			}				
		}//close i
								
		//write only if Nev>0						
		if(Nev_sel>0){						
			m_tables[0]->getColumnList(allcolId,m_nOfAllocatedCol,m_goodEventList,Nev_sel,m_fArray);			
			//for(int i=0;i<Nev_sel;i++) m_goodEventList[i]= index_sel+i;		
			newTable.putColumn(outcolId,m_nOfAllocatedCol,index_sel,index_sel+Nev_sel-1,m_fArray);													
		}//close if	
							
		index+= m_nOfEle;
		index_sel+= Nev_sel;
		if(index+m_nOfEle>= m_nOfRow) m_nOfEle= m_nOfRow-index;
			
	}//close while
 
	  if(mean!=NULL) delete [] mean;
	  if(col!=NULL) delete [] col;
	  if(colId!=NULL) delete [] colId;
	  if(allcolId!=NULL) delete [] allcolId;
	  if(outcolId!=NULL) delete [] outcolId;	

	return true;
}//close function


//---------------------------------------------------------------------
bool VSSigmaContoursOp::allocateArray(){
//---------------------------------------------------------------------
unsigned long long int tempLL=getMaxNumberInt();
if(((unsigned long long int)m_nOfEle*m_nOfAllocatedCol)>tempLL) 
	m_nOfEle=(int)tempLL/m_nOfAllocatedCol;


	try {
		m_fArray= new float*[m_nOfAllocatedCol];
	}
	
	catch(std::bad_alloc &e) {
		m_fArray=NULL;
	}
	

	if(m_fArray==NULL) return false;


	bool goodAllocation=false;
	while(!goodAllocation){
		goodAllocation=true;
		
		try{
			m_goodEventList= new unsigned long long int[m_nOfEle];
		}
		catch(std::bad_alloc &e){
			m_goodEventList= NULL;
		}
		
		for(unsigned int i=0;i<m_nOfAllocatedCol;i++){
			try{
				m_fArray[i]= new float[m_nOfEle];
			}
			catch(std::bad_alloc &e){
				m_fArray[i]= NULL;
			}

			if(m_fArray[i]==NULL||m_goodEventList==NULL) {	
				goodAllocation=false;
				//DELETE MEMORY PREVIOUSLY ALLOCATED
				delete [] m_goodEventList;
				for(unsigned int j=0;j<i;j++) delete [] m_fArray[j];				
				if(m_nOfEle==MIN_NUMBER_OF_ROW){ 
					delete [] m_fArray;
					m_fArray=NULL;
					m_goodEventList=NULL;
					return false;
				}//close if 
				m_nOfEle= m_nOfEle-MAX_NUMBER_TO_REDUCE_ROW;
				if(m_nOfEle<=MIN_NUMBER_OF_ROW) m_nOfEle=MIN_NUMBER_OF_ROW;
				break;
			}//close if

		}//close for i
		
	}//close while
	
	

	return true;
		
}
//---------------------------------------------------------------------


void VSSigmaContoursOp::printHelp(){

	std::cout<<"Create a new data table, where selected columns have values within (or outside) N sigma contours. NOTE: this filter can be executed only for columns that have a Gaussian distribution"<<std::endl<<std::endl;
	std::cout<<"Usage: VisIVOFilters --op sigmacontours [--nsigma nOfSigma] [--field column_list] [--allcolumns] [--out filename_out.bin] [--help] [--file] filename.bin"<<std::endl;

	std::cout<<"Example: VisIVOFilters --op sigmacontours --nsigma 2 --field F1 F2 F5 --allcolumns --out ncontours.bin --file example.bin"<<std::endl;

	std::cout<<"Note: A New table is produced. The command produce a new table where the columns F1, F2 and F5 have (all of them) values included in 2 sigma contours."<<std::endl;
	std::cout<<"      --nsigma. Number of sigma used in the variable selection. Default: 1 sigma countour "<<std::endl;
	std::cout<<"      --field. Selected column variables that must have  values within N sigma contours, or outside  N sigma contours if --exclude option is given. Default value: all the columns in the table"<<std::endl;
	std::cout<<"      --exclude. The output table reports fields outside of N sigma contours."<<std::endl;
	std::cout<<"      --allcolumns. All columns of the original tables are stored in the new table. But only the corresponding rows for the --field selected columns are reported"<<std::endl;
	std::cout<<"      --out. Name of the new table.  Default name is given."<<std::endl;
	std::cout<<"      --file. Input filename"<<std::endl;
	std::cout<<"      --help produce this output "<<std::endl;


	return;

}//close function



