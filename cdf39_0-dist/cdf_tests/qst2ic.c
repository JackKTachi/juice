/******************************************************************************
* Copyright 1996-2014 United States Government as represented by the
* Administrator of the National Aeronautics and Space Administration.
* All Rights Reserved.
******************************************************************************/
/******************************************************************************
*
*  NSSDC/CDF                  Quick Start Test Program (INTERNAL interface/C).
*
*  Version 1.11, 13-Jun-96, Hughes STX.
*
*  Modification history:
*
*   V1.0  27-May-91, J Love     Original version (for CDF V2.1).
*   V1.1  25-Jun-91, J Love     Renamed CDF for portability.
*   V1.2   2-Aug-91, J Love     Use 'Exit'/'ExitBAD'.  Use 'CDFlib'.
*   V1.3  24-Oct-91, J Love     Modified for IBM-PC port.
*   V1.4  21-Apr-92, J Love     CDF V2.2.
*   V1.5   2-Sep-92, J Love     CDF V2.3 (shareable/NeXT/zVar).
*   V1.6  13-Nov-92, J Love     Borland C 3.0 (partial bracketing).
*   V1.7  14-Jan-94, J Love     CDF V2.4.
*   V1.8  20-Dec-94, J Love     CDF V2.5.
*   V1.8a  9-Jan-95, J Love     Indentation on UNIX machines.
*   V1.8b  8-Mar-95, J Love     Fixed wrong data type for rEntry.
*   V1.9   6-Apr-95, J Love     POSIX.
*   V1.9a 18-Apr-95, J Love     More POSIX.
*   V1.10  2-Jan-95, J Love     Shortened argument lists for Think C.
*   V1.11 13-Jun-96, J Love     CDF V2.6.
*   V1.12 13-Mar-23, M Liu      Only handle zVariables only.
*
******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "cdf.h"

#if defined(vms)
#include <ssdef>
#define EXIT_SUCCESS_   SS$_NORMAL
#define EXIT_FAILURE_   SS$_ABORT
#else
#define EXIT_SUCCESS_   0
#define EXIT_FAILURE_   1
#endif

/******************************************************************************
* Increased stack size for Borland C on IBM PC.
******************************************************************************/

#if defined(BORLANDC)
extern unsigned _stklen = 12000u;
#endif

/******************************************************************************
* Macros/prototypes.
******************************************************************************/

#define zN_DIMSa        1
#define zDIM_0_SIZEa    5
#define zNUM_ELEMSa     8

void QuitCDF PROTOARGs((char *where, CDFstatus status));

/******************************************************************************
* Main.
******************************************************************************/

int main () {
CDFid id;
CDFstatus status;
int dim_n;
static long dimSizes[1] = { 0L };
static long encoding = IBMPC_ENCODING;
static long actual_encoding = IBMPC_ENCODING;
static long majority = COL_MAJOR;
static long checkSum = MD5_CHECKSUM, checkSum_out;
static long zNumDimsA = zN_DIMSa;
static long zDimSizesA[zN_DIMSa] = { zDIM_0_SIZEa };
static long zVarAdataType = CDF_CHAR;
static long zVarAdataTypeNew = CDF_UCHAR;
long zVarAdataType_out;
static long zVarAnumElements = zNUM_ELEMSa;
static long zVarAnumElementsNew = zNUM_ELEMSa;
long zVarAnumElements_out;
long zVarAnum_out, varNum_out2;
static char zVarAvalues[zDIM_0_SIZEa][zNUM_ELEMSa] = {
  {'1','1','1','1','1','1','1','1'},
  {'2','2','2','2','2','2','2','2'},
  {'3','3','3','3','3','3','3','3'},
  {'4','4','4','4','4','4','4','4'},
  {'5','5','5','5','5','5','5','5'}
};
static char zVarAvalue_out[zNUM_ELEMSa];
static long recNum = 0;
static long recStart = 0;
static long recCount = 1;
static long recInterval = 1;
static long zRecNum = 0;
static long zRecStart = 0;
static long zRecCount = 1;
static long zRecInterval = 1;
long zIndicesA[zN_DIMSa];
static long zCounts[zN_DIMSa] = { zDIM_0_SIZEa };
static long zIntervals[zN_DIMSa] = { 1 };
char zVarAbuffer_out[zDIM_0_SIZEa][zNUM_ELEMSa];
long attrNum_out;
static long entryNum = 2;
long maxEntry_out;
static long attrScope = GLOBAL_SCOPE;
static long attrScope2 = VARIABLE_SCOPE;
static long attrScope3 = VARIABLE_SCOPE;
long attrScope_out;
static long entryDataType = CDF_INT2;
static long entryDataTypeNew = CDF_UINT2;
long entryDataType_out;
static long entryNumElems = 1;
long entryNumElems_out;
static short entryValue = 1;
short entryValue_out;
long encoding_out;
long majority_out;
long zNumDimsA_out;
long zDimSizesA_out[zN_DIMSa];
long maxRec_out;
long numAttrs_out;
long version_out;
long release_out;
long increment_out;
char subincrement_out;
int i, x0, x;
static long zVarArecVariance = VARY;
static long zVarArecVarianceNew = NOVARY;
long zVarArecVariance_out;
static long zVarAdimVariances[zN_DIMSa] = { VARY };
static long zVarAdimVariancesNew[zN_DIMSa] = { NOVARY };
long zVarAdimVariances_out[zN_DIMSa];
static char zVarAname[] = "zVARa1";
static char new_zVarAname[] = "zVARa2";
char zVarAname_out[CDF_VAR_NAME_LEN256+1];
static char attrName[] = "ATTR1";
static char attrName2[] = "ATTR2";
static char attrName3[] = "ATTR3";
static char new_attrName[] = "ATTR1a";
char attrName_out[CDF_ATTR_NAME_LEN256];
char CopyrightText[CDF_COPYRIGHT_LEN+1];
char errorText[CDF_STATUSTEXT_LEN+1];
static double zEntryValue = 4.0;
double zEntryValueOut;
long numZvars, maxGentry, numGentries,
     maxZentry, numZentries, numGattrs, numVattrs;
static char pad3[zNUM_ELEMSa+1] = { "********" };
static char pad3out[zNUM_ELEMSa+1] = { "        " };
static long blockingfactor3 = 5;
long blockingfactorOut3;
long recStartOut, recCountOut, recIntervalOut, recNumOut;
long indicesOut[CDF_MAX_DIMS],
     countsOut[CDF_MAX_DIMS],
     intervalsOut[CDF_MAX_DIMS];
int dimN;
long entryNumOut1, entryNumOut3;
long maxAllocOut3;
long maxRecOut3, maxRecOut;
long nIndexRecsOut3;
long nIndexEntriesOut3;
static long allocRecs3 = 8;
static long nZvars = 1;
static long zVarNs[1] = { 0 };
static char zVarsRecBuffer[zDIM_0_SIZEa][zNUM_ELEMSa] = {
  {'%','%','%','%','%','%','%','%'},
  {'%','%','%','%','%','%','%','%'},
  {'%','%','%','%','%','%','%','%'},
  {'%','%','%','%','%','%','%','%'},
  {'%','%','%','%','%','%','%','%'}
};
static char zVarsRecBufferOut[zDIM_0_SIZEa][zNUM_ELEMSa];
long lastUpdated;
static long lastUpdatedO = 20170101;

/******************************************************************************
* Display title.
******************************************************************************/

printf ("Testing Internal/C interface...\n");

/* CDFsetFileBackward(BACKWARDFILEon); */
/* CDFsetValidate(VALIDATEFILEon); */

/******************************************************************************
* Create CDF.
******************************************************************************/

status = CDFlib (CREATE_, CDF_, "TEST", 0L, dimSizes, &id,
		 PUT_, CDF_MAJORITY_, majority,
		 NULL_);

if (status < CDF_OK) {
  if (status == CDF_EXISTS) {
    status = CDFlib (OPEN_, CDF_, "TEST", &id,
		     NULL_);
    if (status < CDF_OK) QuitCDF ("1.0", status);

    status = CDFlib (DELETE_, CDF_,
		     NULL_);
    if (status < CDF_OK) QuitCDF ("1.1", status);

    status = CDFlib (CREATE_, CDF_, "TEST", 0L, dimSizes, &id,
		     PUT_, CDF_MAJORITY_, majority,
		     NULL_);
    if (status < CDF_OK) QuitCDF ("1.2", status);
  }
  else
    QuitCDF ("1.3", status);
}

/******************************************************************************
* Create variables.
******************************************************************************/

status = CDFlib (CREATE_, zVAR_, zVarAname, zVarAdataType, zVarAnumElements,
				 zNumDimsA, zDimSizesA, zVarArecVariance,
				 zVarAdimVariances, &zVarAnum_out,
		 PUT_, zVAR_PADVALUE_, pad3,
                 NULL_);
if (status < CDF_OK) QuitCDF ("2.0", status);
status = CDFlib (SELECT_, CDF_, id,
		          zVAR_, zVarAnum_out,
		 GET_, zVAR_PADVALUE_, pad3out,
                 NULL_);
if (status < CDF_OK) QuitCDF ("2.2", status);

if (strcmp(pad3out,pad3)) QuitCDF ("2.6", status);

/******************************************************************************
* Modify variables.
******************************************************************************/

status = CDFlib (SELECT_, zVAR_, 0L,
		 PUT_, zVAR_DATASPEC_, zVarAdataTypeNew, zVarAnumElementsNew,
		       zVAR_RECVARY_, zVarArecVarianceNew,
		       zVAR_DIMVARYS_, zVarAdimVariancesNew,
		       zVAR_INITIALRECS_, 1L,
		 NULL_);
if (status < CDF_OK) QuitCDF ("2b.04", status);

/******************************************************************************
* Check leap second last upadted.
******************************************************************************/

status = CDFlib (GET_, CDF_LEAPSECONDLASTUPDATED_, &lastUpdated,
                 NULL_);
if (status < CDF_OK) QuitCDF ("2.7", status);
if ((int)lastUpdated != lastUpdatedO) QuitCDF("2.8", status);

/******************************************************************************
* Close CDF.
******************************************************************************/

status = CDFlib (CLOSE_, CDF_,
		 NULL_);
if (status < CDF_OK) QuitCDF ("3.0", status);

/******************************************************************************
* Reopen CDF.
******************************************************************************/

status = CDFlib (OPEN_, CDF_, "TEST", &id,
		 SELECT_, CDF_DECODING_, HOST_DECODING,
		 NULL_);
if (status < CDF_OK) QuitCDF ("4.0", status);

/******************************************************************************
* Delete CDF.
******************************************************************************/

status = CDFlib (DELETE_, CDF_,
		 NULL_);
if (status < CDF_OK) QuitCDF ("5.0", status);

/******************************************************************************
* Create CDF again (previous delete will allow this).
******************************************************************************/

status = CDFlib (CREATE_, CDF_, "TEST", 0L, dimSizes, &id,
		 PUT_, CDF_MAJORITY_, majority,
		 SELECT_, CDF_DECODING_, HOST_DECODING,
		 NULL_);
if (status < CDF_OK) QuitCDF ("6.0", status);

/******************************************************************************
* Set the checksum (to MD5).
******************************************************************************/

status = CDFlib (PUT_, CDF_CHECKSUM_, checkSum,
                 NULL_);
if (status < CDF_OK) QuitCDF ("6.5", status);

/******************************************************************************
* Create variables.
******************************************************************************/

status = CDFlib (CREATE_, zVAR_, zVarAname, zVarAdataType, zVarAnumElements,
				 zNumDimsA, zDimSizesA, zVarArecVariance,
				 zVarAdimVariances, &zVarAnum_out,
		 PUT_, zVAR_ALLOCATERECS_, allocRecs3,
		       zVAR_BLOCKINGFACTOR_, blockingfactor3,
		 NULL_);
if (status < CDF_OK) QuitCDF ("7.0", status);

/******************************************************************************
* PUT to variables.
******************************************************************************/

status = CDFlib (SELECT_, zVAR_, zVarAnum_out,
			  zVAR_RECNUMBER_, zRecNum,
		 NULL_);
if (status < CDF_OK) QuitCDF ("8.2", status);

for (x0 = 0; x0 < zDIM_0_SIZEa; x0++) {
   zIndicesA[0] = x0;
   status = CDFlib (SELECT_, zVAR_DIMINDICES_, zIndicesA,
		    PUT_, zVAR_DATA_, zVarAvalues[x0],
		    NULL_);
   if (status < CDF_OK) QuitCDF ("8.3", status);
}

/******************************************************************************
* GET from the variables.
******************************************************************************/

for (x0 = 0; x0 < zDIM_0_SIZEa; x0++) {
   zIndicesA[0] = x0;
   status = CDFlib (SELECT_, zVAR_DIMINDICES_, zIndicesA,
		    GET_, zVAR_DATA_, zVarAvalue_out,
		    NULL_);
   if (status < CDF_OK) QuitCDF ("9.3", status);

   for (i = 0; i < zNUM_ELEMSa; i++) {
      if (zVarAvalue_out[i] != zVarAvalues[x0][i]) QuitCDF ("9.4", status);
   }
}

/******************************************************************************
* HyperPUT to the variables.
******************************************************************************/

for (x0 = 0; x0 < zDIM_0_SIZEa; x0++)
   for (i = 0; i < zNUM_ELEMSa; i++) {
      zVarAvalues[x0][i]++;
   }

zIndicesA[0] = 0;

status = CDFlib (SELECT_, zVAR_RECNUMBER_, zRecStart,
			  zVAR_RECCOUNT_, zRecCount,
			  zVAR_RECINTERVAL_, zRecInterval,
			  zVAR_DIMINDICES_, zIndicesA,
			  zVAR_DIMCOUNTS_, zCounts,
			  zVAR_DIMINTERVALS_, zIntervals,
		 PUT_, zVAR_HYPERDATA_, zVarAvalues,
		 NULL_);
if (status < CDF_OK) QuitCDF ("10.1", status);

/******************************************************************************
* HyperGET from variables.
******************************************************************************/

status = CDFlib (GET_, zVAR_HYPERDATA_, zVarAbuffer_out,
		 NULL_);
if (status < CDF_OK) QuitCDF ("11.3", status);

for (x0 = 0; x0 < zDIM_0_SIZEa; x0++)
   for (i = 0; i < zNUM_ELEMSa; i++) {
      if (zVarAbuffer_out[x0][i] != zVarAvalues[x0][i])
	QuitCDF ("11.4", status);
   }

/******************************************************************************
* Confirm hyper parameters for a zVariable.
******************************************************************************/

status = CDFlib (CONFIRM_, zVAR_RECNUMBER_, &recStartOut,
			   zVAR_RECCOUNT_, &recCountOut,
			   zVAR_RECINTERVAL_, &recIntervalOut,
			   zVAR_DIMINDICES_, indicesOut,
			   zVAR_DIMCOUNTS_, countsOut,
			   zVAR_DIMINTERVALS_, intervalsOut,
		 NULL_);
if (status < CDF_OK) QuitCDF ("11a.0", status);

if (recStartOut != zRecStart) QuitCDF ("11a.1", status);
if (recCountOut != zRecCount) QuitCDF ("11a.2", status);
if (recIntervalOut != zRecInterval) QuitCDF ("11a.3", status);
for (dimN = 0; dimN < zN_DIMSa; dimN++) {
   if (indicesOut[dimN] != zIndicesA[dimN]) QuitCDF ("11a.4", status);
   if (countsOut[dimN] != zCounts[dimN]) QuitCDF ("11a.5", status);
   if (intervalsOut[dimN] != zIntervals[dimN]) QuitCDF ("11a.6", status);
}

/******************************************************************************
* Set/confirm sequential access position for a zVariable (and read/write a
* value).
******************************************************************************/

status = CDFlib (SELECT_, zVAR_SEQPOS_, zRecStart, zIndicesA,
		 GET_, zVAR_SEQDATA_, zVarAvalue_out,
		 PUT_, zVAR_SEQDATA_, zVarAvalue_out,
		 CONFIRM_, zVAR_SEQPOS_, &recNumOut, indicesOut,
		 NULL_);
if (status < CDF_OK) QuitCDF ("11b.0", status);

if (recNumOut != zRecStart) QuitCDF ("11b.1", status);
if (indicesOut[0] != zIndicesA[0] + 2) QuitCDF ("11b.2", status);

/******************************************************************************
* Create attributes.
******************************************************************************/

status = CDFlib (CREATE_, ATTR_, attrName, attrScope, &attrNum_out,
			  ATTR_, attrName2, attrScope2, &attrNum_out,
			  ATTR_, attrName3, attrScope3, &attrNum_out,
		 NULL_);
if (status < CDF_OK) QuitCDF ("12.0", status);

/******************************************************************************
* PUT to attributes.
******************************************************************************/

status = CDFlib (SELECT_, ATTR_, 0L,
			  gENTRY_, entryNum,
		 PUT_, gENTRY_DATA_, entryDataType, entryNumElems,
				     &entryValue,
		 SELECT_, ATTR_, 1L,
			  zENTRY_NAME_, zVarAname,
		 PUT_, zENTRY_DATA_, CDF_REAL8, 1L, &zEntryValue,
		 NULL_);
if (status < CDF_OK) QuitCDF ("13.0", status);

/******************************************************************************
* Confirm entry numbers.
******************************************************************************/

status = CDFlib (CONFIRM_, gENTRY_, &entryNumOut1,
			   zENTRY_, &entryNumOut3,
		 NULL_);
if (status < CDF_OK) QuitCDF ("13a.0", status);
if (entryNumOut1 != 2) QuitCDF ("13a.1", status);
if (entryNumOut3 != 0) QuitCDF ("13a.3", status);

/******************************************************************************
* GET from attributes.
******************************************************************************/

status = CDFlib (SELECT_, ATTR_, 0L,
			  gENTRY_, entryNum,
		 CONFIRM_, CURgENTRY_EXISTENCE_,
		 GET_, gENTRY_DATA_, &entryValue_out,
		 SELECT_, ATTR_, 1L,
			  zENTRY_, 0L,
		 CONFIRM_, CURzENTRY_EXISTENCE_,
		 GET_, zENTRY_DATA_, &zEntryValueOut,
		 NULL_);
if (status < CDF_OK) QuitCDF ("14.0", status);

if (entryValue_out != entryValue) QuitCDF ("14.1", status);
if (zEntryValue != zEntryValueOut) QuitCDF ("14.3", status);

/******************************************************************************
* Confirm existence of variables/attributes/entries.
******************************************************************************/

status = CDFlib (CONFIRM_, zVAR_EXISTENCE_, zVarAname,
			   ATTR_EXISTENCE_, attrName3,
		 NULL_);
if (status < CDF_OK) QuitCDF ("14a.0", status);

status = CDFlib (SELECT_, ATTR_, 0L,
		 CONFIRM_, gENTRY_EXISTENCE_, entryNum,
		 SELECT_, ATTR_, 1L,
		 CONFIRM_, zENTRY_EXISTENCE_, 0L,
		 NULL_);
if (status < CDF_OK) QuitCDF ("14a.1", status);

/******************************************************************************
* Get CDF documentation.
******************************************************************************/

status = CDFlib (GET_, LIB_VERSION_, &version_out,
		       LIB_RELEASE_, &release_out,
		       LIB_INCREMENT_, &increment_out,
		       LIB_subINCREMENT_, &subincrement_out,
		       LIB_COPYRIGHT_, CopyrightText,
		 NULL_);
if (status < CDF_OK) QuitCDF ("15.0", status);

/******************************************************************************
* Inquire CDF.
******************************************************************************/

status = CDFlib (GET_, CDF_ENCODING_, &encoding_out,
		       CDF_MAJORITY_, &majority_out,
		       CDF_NUMzVARS_, &numZvars,
		       CDF_NUMATTRS_, &numAttrs_out,
                       CDF_CHECKSUM_, &checkSum_out,
		 NULL_);
if (status < CDF_OK) QuitCDF ("16.0", status);

/* if (encoding_out != actual_encoding) QuitCDF ("16.4", status); */
if (majority_out != majority) QuitCDF ("16.5", status);
/* if (maxRec_out != 0) QuitCDF ("16.6", status); */
if (numZvars != 1) QuitCDF ("16.8", status);
if (numAttrs_out != 3) QuitCDF ("16.9", status);
if (checkSum != checkSum_out) QuitCDF ("16.10", status);

/******************************************************************************
* Inquire numbers.
******************************************************************************/

status = CDFlib (GET_, ATTR_NUMBER_, attrName3, &attrNum_out,
		       zVAR_NUMBER_, zVarAname, &varNum_out2,
		 NULL_);
if (status < CDF_OK) QuitCDF ("16a.0", status);

if (attrNum_out != 2) QuitCDF ("16a.1", status);
if (varNum_out2 != 0) QuitCDF ("16a.3", status);

/******************************************************************************
* Rename variables.
******************************************************************************/

status = CDFlib (SELECT_, zVAR_NAME_, zVarAname,
		 PUT_, zVAR_NAME_, new_zVarAname,
		 NULL_);
if (status < CDF_OK) QuitCDF ("17.0", status);

/******************************************************************************
* Read/write multiple variable data.
******************************************************************************/

status = CDFlib (SELECT_, zVARs_RECNUMBER_, 2L,
		 PUT_, zVARs_RECDATA_, nZvars, zVarNs, zVarsRecBuffer,
		 GET_, zVARs_RECDATA_, nZvars, zVarNs, zVarsRecBufferOut,
		 NULL_);
if (status < CDF_OK) QuitCDF ("17a.0", status);

if (memcmp(zVarsRecBufferOut,zVarsRecBuffer,
	   sizeof(zVarsRecBuffer))) QuitCDF ("17a.2", status);

/******************************************************************************
* Inquire variables.
******************************************************************************/

status = CDFlib (SELECT_, zVAR_, zVarAnum_out,
		 GET_, zVAR_NAME_, zVarAname_out,
		       zVAR_DATATYPE_, &zVarAdataType_out,
		       zVAR_NUMELEMS_, &zVarAnumElements_out,
		       zVAR_RECVARY_, &zVarArecVariance_out,
		       zVAR_DIMVARYS_, zVarAdimVariances_out,
		       zVAR_NUMDIMS_, &zNumDimsA_out,
		       zVAR_DIMSIZES_, zDimSizesA_out,
		       zVAR_BLOCKINGFACTOR_, &blockingfactorOut3,
		       zVAR_MAXallocREC_, &maxAllocOut3,
		       zVAR_MAXREC_, &maxRecOut3,
		       zVAR_nINDEXRECORDS_, &nIndexRecsOut3,
		       zVAR_nINDEXENTRIES_, &nIndexEntriesOut3,
		 CONFIRM_, zVAR_, &zVarAnum_out,
		 NULL_);
if (status < CDF_OK) QuitCDF ("18b.0", status);

if (strcmp(zVarAname_out,new_zVarAname) != 0) QuitCDF ("18b.1", status);
if (zVarAdataType_out != zVarAdataType) QuitCDF ("18b.2", status);
if (zVarAnumElements_out != zVarAnumElements) QuitCDF ("18b.3", status);
if (zVarArecVariance_out != zVarArecVariance) QuitCDF ("18b.4", status);
if (zNumDimsA_out != zNumDimsA) QuitCDF ("18b.5", status);
if (zVarAnum_out != 0L) QuitCDF ("18b.6", status);
if (blockingfactorOut3 != blockingfactor3) QuitCDF ("18b.7", status);
if (maxAllocOut3 + 1 != allocRecs3) QuitCDF ("18b.8", status);
if (maxRecOut3 != 2L) QuitCDF ("18b.9", status);

for (dim_n = 0; dim_n < zNumDimsA; dim_n++) {
   if (zDimSizesA_out[dim_n] != zDimSizesA[dim_n]) {
     QuitCDF ("18b.10", status);
   }
   if (zVarAdimVariances_out[dim_n] != zVarAdimVariances[dim_n]) {
     QuitCDF ("18b.11", status);
   }
}

/******************************************************************************
* Rename attribute.
******************************************************************************/

status = CDFlib (SELECT_, ATTR_NAME_, attrName,
		 PUT_, ATTR_NAME_, new_attrName,
		 NULL_);
if (status < CDF_OK) QuitCDF ("20.0", status);

/******************************************************************************
* Inquire attribute.
******************************************************************************/

status = CDFlib (GET_, ATTR_NAME_, attrName_out,
		       ATTR_SCOPE_, &attrScope_out,
		       ATTR_MAXgENTRY_, &maxEntry_out,
		 CONFIRM_, ATTR_, &attrNum_out,
		 NULL_);
if (status < CDF_OK) QuitCDF ("22.0", status);

if (strcmp(attrName_out,new_attrName) != 0) QuitCDF ("22.1", status);
if (attrScope_out != attrScope) QuitCDF ("22.2", status);
if (maxEntry_out != entryNum) QuitCDF ("22.3", status);
if (attrNum_out != 0L) QuitCDF ("22.4", status);

/******************************************************************************
* Inquire attribute entries.
******************************************************************************/

status = CDFlib (SELECT_, ATTR_, 0L,
			  gENTRY_, entryNum,
		 GET_, gENTRY_DATATYPE_, &entryDataType_out,
		       gENTRY_NUMELEMS_, &entryNumElems_out,
		 NULL_);
if (status < CDF_OK) QuitCDF ("23.0", status);

if (entryDataType_out != entryDataType) QuitCDF ("23.1", status);
if (entryNumElems_out != entryNumElems) QuitCDF ("23.2", status);

status = CDFlib (SELECT_, ATTR_, 1L,
			  zENTRY_, 0L,
		 GET_, zENTRY_DATATYPE_, &entryDataType_out,
		       zENTRY_NUMELEMS_, &entryNumElems_out,
		 NULL_);
if (status < CDF_OK) QuitCDF ("23.7", status);

if (entryDataType_out != CDF_REAL8) QuitCDF ("23.8", status);
if (entryNumElems_out != 1L) QuitCDF ("23.9", status);

/******************************************************************************
* Get error text.
******************************************************************************/

status = CDFlib (SELECT_, CDF_STATUS_, CDF_OK,
		 GET_, STATUS_TEXT_, errorText,
		 NULL_);
if (status < CDF_OK) QuitCDF ("24.0", status);

/******************************************************************************
* Select zMode and inquire CDF.
******************************************************************************/

status = CDFlib (SELECT_, CDF_zMODE_, zMODEon2,
		 SELECT_, ATTR_, 0L,
		 GET_, CDF_NUMgATTRS_, &numGattrs,
		       CDF_NUMvATTRS_, &numVattrs,
		       CDF_NUMzVARS_, &numZvars,
		       ATTR_MAXgENTRY_, &maxGentry,
		       ATTR_NUMgENTRIES_, &numGentries,
		       zVARs_MAXREC_, &maxRecOut,
		 SELECT_, ATTR_, 1L,
		 GET_, ATTR_MAXzENTRY_, &maxZentry,
		       ATTR_NUMzENTRIES_, &numZentries,
		 SELECT_, CDF_zMODE_, zMODEoff,
		 NULL_);
if (status < CDF_OK) QuitCDF ("25.0", status);

if (numGattrs != 1) QuitCDF ("25.1", status);
if (numVattrs != 2) QuitCDF ("25.2", status);
if (numZvars != 1) QuitCDF ("25.4", status);
if (maxGentry != entryNum) QuitCDF ("25.5", status);
if (numGentries != 1) QuitCDF ("25.6", status);
if (maxZentry != 0) QuitCDF ("25.8", status);
if (numZentries != 1) QuitCDF ("25.9", status);
if (maxRecOut != 2L) QuitCDF ("25.10", status);

/******************************************************************************
* Modify entries/attribute.
******************************************************************************/

status = CDFlib (SELECT_, ATTR_, 0L,
			  gENTRY_, entryNum,
		 PUT_, gENTRY_DATASPEC_, entryDataTypeNew, entryNumElems,
		 SELECT_, ATTR_, 1L,
			  zENTRY_, 0L,
		 PUT_, zENTRY_DATASPEC_, CDF_EPOCH, 1L,
		 SELECT_, ATTR_, 0L,
		 PUT_, ATTR_SCOPE_, VARIABLE_SCOPE,
		       ATTR_SCOPE_, GLOBAL_SCOPE,
		 NULL_);
if (status < CDF_OK) QuitCDF ("27.0", status);

/******************************************************************************
* Delete entries/attribute/variables and variable record.
******************************************************************************/

status = CDFlib (SELECT_, ATTR_, 0L,
			  gENTRY_, entryNum,
		 DELETE_, gENTRY_, 
		 SELECT_, ATTR_, 1L,
			  zENTRY_, 0L,
		 DELETE_, zENTRY_, 
		 NULL_);
if (status < CDF_OK) QuitCDF ("28.0", status);

status = CDFlib (SELECT_, ATTR_, 0L,
		 DELETE_, ATTR_, 
		 SELECT_, zVAR_, 0L,
		 DELETE_, zVAR_RECORDS_, 1L, 1L,
		 NULL_);
if (status < CDF_OK) QuitCDF ("28.1", status);

/******************************************************************************
* Close CDF.
******************************************************************************/

status = CDFlib (CLOSE_, CDF_,
		 NULL_);
if (status < CDF_OK) QuitCDF ("29.0", status);

/******************************************************************************
* Reopen and close CDF to validate the file.
******************************************************************************/

status = CDFlib (OPEN_, CDF_, "TEST", &id,
                 NULL_);
if (status < CDF_OK) QuitCDF ("36.0", status);

status = CDFlib (CLOSE_, CDF_,
                 NULL_);
if (status < CDF_OK) QuitCDF ("37.0", status);

/******************************************************************************
* Successful completion.
******************************************************************************/

return EXIT_SUCCESS_;
}


/******************************************************************************
* QuitCDF.
******************************************************************************/

void QuitCDF (where, status)
char *where;
CDFstatus status;
{
  char text[CDF_STATUSTEXT_LEN+1];
  printf ("Program failed at %s...\n", where);
  if (status < CDF_OK) {
    CDFlib (SELECT_, CDF_STATUS_, status,
	    GET_, STATUS_TEXT_, text,
	    NULL_);
    printf ("%s\n", text);
  }
  CDFlib (CLOSE_, CDF_,
	  NULL_);
  printf ("...test aborted.\n");
  exit (EXIT_FAILURE_);
}
