/******************************************************************************
 *
 *  SPDF/CDF                  CDF read functions to acquire full information.
 *
 *  Version 1.0, 10-Mar-23, Hughes STX.
 *
 *  Modification history:
 *
 *   V1.0  10-Mar-23, M Liu      Original version 
 *                               Show how CDF read functions to retrieve 
 *                               attributes and varibable data are called.
 *
 ******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "cdf.h"

int Do_reads (CDFid id, int i);
int Read_GlobalAttributes (CDFid id, int i);
int Read_VariablesAttributes (CDFid id, int i);
int Read_VariablesData (CDFid id, int i);

/* Read all global attributes */

int Read_GlobalAttributes (CDFid id, int ix) {
    CDFstatus status;
    int i, j, k, count;
    long numAttrs, numgAttrs, scope, maxEntry, dataType, numElems;
    char attrName[CDF_ATTR_NAME_LEN256+1];
    CDFdata *entryValues;
    CDFstatus *statusAll;
    long *dataTypes, *numElemsAll;

    status = CDFgetNumgAttributes (id, &numgAttrs);
    if (status != 0) {
      CDFcloseCDF (id);
      return 1;
    }
    if (numgAttrs > 0) {
      status = CDFgetNumAttributes (id, &numAttrs);

      for (j = 0; j < (int) numAttrs; ++j) {
        status = CDFattrInquire (id, (long) j, attrName, &scope,
                                 &maxEntry);
        if (scope != GLOBAL_SCOPE) continue;
        count = maxEntry + 1;
        entryValues = (CDFdata *) malloc(sizeof(CDFdata) * count);
        statusAll = (CDFstatus *) malloc(sizeof(CDFstatus) * count);
        dataTypes = (long *) malloc(sizeof(long) * count);
        numElemsAll = (long *) malloc(sizeof(long) * count);
        for (k = 0; k < count; ++k) {
          statusAll[k] = CDFreadAttrgEntry (id, (long) j, (long) k,
                                            &dataTypes[k], &numElemsAll[k], 
                                            &entryValues[k]);
        }
        for (k = 0; k < (int) count; ++k) {
          if (statusAll[k] != CDF_OK) continue;
          if (dataTypes[k] == CDF_CHAR || dataTypes[k] == CDF_UCHAR) { 
/*
            printf("cdf:%d g_attr(%d):%s entry:(%d) => [%s] (datatype:%ld "
                   "numelms:%ld)\n", ix, j, attrName,
                   k, (char *)entryValues[k], dataTypes[k],
                   numElemsAll[k]);
            fflush(stdout);
*/
          } else {
            int n; char text[81];
/*
            n = EncodeValuesFormat (dataTypes[k], numElemsAll[k],
                                    (void *)entryValues[k], text, NULL, 0, 0, 3,
                                    80);
            printf("cdf:%d g_attr(%d):%s entry:(%ld) => [%s] (datatype:%ld "
                   "numelms:%ld)\n", ix, j, attrName, 
                   k, text, dataTypes[k], numElemsAll[k]);
            fflush(stdout);
*/
          }
          CDFdataFree (entryValues[k]);
        }
        free (entryValues);
        free (statusAll);
        free (dataTypes);
        free (numElemsAll);
      }
    }
    return 0;
}

/* Read all variables' attributes */

int Read_VariablesAttributes (CDFid id, int ix) {
    int i, j, k;
    long numAttrs, numvAttrs, numrVars, numzVars, maxEntry, scope, dataType,
         numElems;
    char attrName[CDF_ATTR_NAME_LEN256+1], varName[CDF_VAR_NAME_LEN256+1];
    CDFstatus status;
    long *dataTypes, *numElemsAll, *attrNums;
    CDFstatus *statusAll;
    CDFdata *entryValues;

    status = CDFgetNumAttributes (id, &numAttrs);

    status = CDFgetNumvAttributes (id, &numvAttrs);

    status = CDFgetNumzVars (id, &numzVars);

    if (numzVars > 0 && numvAttrs > 0) {
      entryValues = (CDFdata *) malloc (sizeof(CDFdata) * numvAttrs);
      attrNums = (long *) malloc (sizeof(long) * numvAttrs);
      dataTypes = (long *) malloc (sizeof(long) * numvAttrs);
      numElemsAll = (long *) malloc (sizeof(long) * numvAttrs);
      statusAll = (CDFstatus *) malloc (sizeof(CDFstatus) * numvAttrs);
      for (k = 0; k < (int) numzVars; ++k) {
        status = CDFgetzVarName(id, k, varName);
        for (i = j = 0; j < numAttrs; ++j) {
            status = CDFattrInquire (id, (long) j, attrName, &scope,
                                     &maxEntry);
            if (scope == GLOBAL_SCOPE) continue;
            statusAll[i] = CDFreadAttrzEntry (id, (long) j, (long) k,
                                              &dataTypes[i], &numElemsAll[i],
                                              &entryValues[i]);
            if (statusAll[i] == CDF_OK) {
              attrNums[i] = j;
              i++;
            }
        }
        for (j = 0; j < i; ++j) {
            status = CDFgetAttrName (id, attrNums[j], attrName);
            if (dataTypes[j] == CDF_CHAR || dataTypes[j] == CDF_UCHAR) {
/*
              printf("cdf:%d v_attr(%ld):%s var:%s (%d) => [%s] (datatype:%ld "
                     "numelms:%ld)\n", ix, attrNums[j], 
                     attrName, varName, k, (char *) entryValues[j],
                     dataTypes[j], numElemsAll[j]);
*/
            } else {
              int n; char text[81];
/*
              n = EncodeValuesFormat (dataTypes[j], numElemsAll[j],
                                      (void *)entryValues[j], text, NULL, 0, 0,
                                      3, 80);
              printf("cdf:%d v_attr(%ld):%s var:%s (%d) => [%s] (datatype:%ld "
                     "numelms:%ld)\n", ix, attrNums[j],
                     attrName, varName, k, text,
                     dataTypes[j], numElemsAll[j]);
*/
            }
            CDFdataFree (entryValues[j]);
        }
      }
      free (statusAll);
      free (dataTypes);
      free (numElemsAll);
      free (attrNums);
    }
    return 0;
}

/* Read all variables' data */

int Read_VariablesData (CDFid id, int ix) {
    CDFstatus status;
    long numrVars, numzVars, *numRecs;
    char varName[CDF_VAR_NAME_LEN256+1];
    int i, j, k;
    CDFdata *dataValues;

    status = CDFgetNumzVars (id, &numzVars);
    if (numzVars > 0) {
        dataValues = (CDFdata *) malloc (sizeof(CDFdata) * numzVars);
        numRecs = (long *) malloc (sizeof(long) * numzVars);
        for (k = 0; k < numzVars; ++k) {
          status = CDFreadzVarDataByVarID (id, (long) k, &numRecs[k],
                                           &dataValues[k]);
        }
/*
        for (k = 0; k < numzVars; ++k) {
          status = CDFgetzVarName(id, k, varName);
          printf("cdf:%d [var:%s (%d) => %ld rec%s]\n", ix, varName,
                 k, numRecs[k],(numRecs[k]>1?"s":""));
          printf("\n");
*/
          CDFdataFree (dataValues[k]);
        }
        free (dataValues);
        free (numRecs);
    }
    return 0;
}

/* Each CDF is broken into 3 steps, reading global attribute, variable
 * attributes and data. */

int Do_reads (CDFid id, int ix) {
    int res;
 
    res = Read_GlobalAttributes (id, ix);
    if (res != 0) return res;
    res = Read_VariablesAttributes (id, ix);
    if (res != 0) return res;
    res = Read_VariablesData (id, ix);
    return res;
}

int main(int argc, char **argv) {
    CDFstatus status;
    int res;
    CDFid id;
    int i, j, k;

    for (i = 1; i < argc; ++i) {    
      status = CDFopenCDF (argv[i], &id);
      if (status != CDF_OK) {
        return 1;
      }
      status = CDFsetzMode (id, zMODEon2);
      res = Do_reads (id, i); 
      status = CDFcloseCDF (id);
    }
    return 0;
}

