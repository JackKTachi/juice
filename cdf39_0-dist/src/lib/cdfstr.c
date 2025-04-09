/******************************************************************************
* Copyright 1996-2014 United States Government as represented by the
* Administrator of the National Aeronautics and Space Administration.
* All Rights Reserved.
******************************************************************************/
/******************************************************************************
*
*  NSSDC/CDF				      CDF character string handling.
*
*  Version 1.7a, 28-Feb-97, Hughes STX.
*
*  Modification history:
*
*   V1.0  16-May-91, J Love	Original version (for CDF V2.1).  This is a
*				combination of nulterminate.c, isdescr.c,
*				descrtoref.c, and descrtorefnul.c.
*   V1.1  31-Jul-91, J Love	Added 'array' class to checking for a
*				descriptor (DSC$K_CLASS_A).
*   V1.2  23-Sep-91, J Love	Modified for IBM-PC port.
*   V1.3   5-Apr-92, J Love	CDF V2.2.
*   V1.4   7-Jul-92, J Love	CDF V2.3 (shareable/NeXT/zVar).
*   V1.5   1-Mar-93, J Love	CDF V2.4.
*   V1.5a  5-Feb-94, J Love	DEC Alpha/OpenVMS port.
*   V1.6  13-Dec-94, J Love	CDF V2.5.
*   V1.6a 10-Jan-95, J Love	Added `Trailer'.
*   V1.6b  6-Mar-95, J Love	Added `IsGraphChr'.  Solaris 2.3 IDL i/f.
*				`strstr' moved here (for SunOS before 4.1).
*   V1.7  22-Aug-96, J Love	CDF V2.6.
*   V1.7a 28-Feb-97, J Love	Windows NT for MS Visual C/C++ on an IBM PC.
*   V1.8  12-Aug-10, M Liu 	Nothing is performed if negative max length
*                               is passed in in strcpyX, strcatX, catchrX,
*                               PstrcatX and PstrcpyX.
*   V1.9  12-Aug-21, M Liu 	Added a few UTF-8 encoding handling functions.
*
******************************************************************************/

/******************************************************************************
 *  NOTE:
 *   Never assume that a string may be NUL-terminated already.  If it is,
 *   it might be accidental and the NUL might go away.
 *****************************************************************************/

#if defined(unix) || defined(__CYGWIN__)
#  define _XOPEN_SOURCE
#endif
#define _XOPEN_SOURCE_EXTENDED

#include <stdlib.h>
#include <wchar.h>
#include "cdflib.h"

#if defined(Fif_DESCR)
#if defined(vms)
#include <descrip.h>
#endif
#if defined(MPW_C)
#define MPW_FORTRAN_STRING	18
struct DESCRstruct {
  void *ptr;
  short length;
  Byte1 filler;
  Byte1 type;
};
#endif
#endif

#define IsASCII(c) (0 <= c && c <= 127)

/******************************************************************************
 * NULterminate.  NUL-terminate a character string.
 *****************************************************************************/

STATICforIDL char *NULterminate (string, length, ssh)
char	*string;			/* String to NUL-terminate. */
size_t	length;				/* Length of string (not including
					   the NUL terminator. */
struct STRINGstruct **ssh;		/* Pointer to head of STRINGstruct
					   linked list. */
{
  struct STRINGstruct *ss;

  if (*ssh == NULL) {	/* First string. */
    *ssh = (struct STRINGstruct *)
		 cdf_AllocateMemory ((size_t)sizeof(struct STRINGstruct),
				     NULL);
    if (*ssh == NULL) return NULL;
    ss = *ssh;
  }
  else {
    ss = *ssh;
    while (ss->next != NULL) ss = ss->next;
    ss->next = (struct STRINGstruct *)
	       cdf_AllocateMemory ((size_t)sizeof(struct STRINGstruct), NULL);
    if (ss->next == NULL) return NULL;
    ss = ss->next;
  }

  ss->next = NULL;

  ss->string = (char *) cdf_AllocateMemory ((size_t)length + 1, NULL);
  if (ss->string == NULL) return NULL;

  memmove (ss->string, string, length);
  ss->string[length] = NUL;

  return ss->string;
}

/******************************************************************************
 * FreeStrings.  Free all strings (and control structures) used during the
 * call.
 *****************************************************************************/

STATICforIDL void FreeStrings (ssh)
struct STRINGstruct *ssh;		/* Head of STRINGstruct chain. */
{
  struct STRINGstruct *ss, *sst;
  ss = ssh;
  while (ss != NULL) {
    cdf_FreeMemory (ss->string, NULL);
    sst = ss;
    ss = ss->next;
    cdf_FreeMemory (sst, NULL);
  }
  return;
}

#if defined(Fif_DESCR)
/******************************************************************************
 * isDESCR.  Is a descriptor being pointed to?
 *
 * NOTE:
 *     There is no checking on the length of the string because all lengths
 * are valid to VMS and to CDF (for attributes with character data types).  It
 * is assumed that all lengths will also be valid on the Macintosh using MPW
 * Fortran.
 *
 *****************************************************************************/

STATICforIDL Logical isDESCR (ptr, aptr, len)
void	*ptr;		/* could point to a descriptor or to the
			   string (if already passed by reference) */
char	**aptr;		/* pointer to string */
size_t	*len;		/* length of string */
{
#if defined(vms)
  struct dsc$descriptor *tPtr = (struct dsc$descriptor *) ptr;
  if (tPtr->dsc$b_dtype == DSC$K_DTYPE_T &&
      (tPtr->dsc$b_class == DSC$K_CLASS_S ||
       tPtr->dsc$b_class == DSC$K_CLASS_A ||
       tPtr->dsc$b_class == DSC$K_CLASS_NCA)) {
    *len = tPtr->dsc$w_length;
    *aptr = tPtr->dsc$a_pointer;
    return TRUE;
  }
  else
    return FALSE;
#endif
#if defined(MPW_C)
  struct DESCRstruct *tPtr = (struct DESCRstruct *) ptr;
  if (tPtr->type == MPW_FORTRAN_STRING && 0 <= tPtr->length) {
    *len = (size_t) tPtr->length;
    *aptr = (char *) tPtr->ptr;
    return TRUE;
  }
  else
    return FALSE;
#endif
}
#endif

#if defined(Fif_DESCR)
/******************************************************************************
 * DESCRtoREF.  Convert from passing by descriptor to passing by reference.
 *****************************************************************************/

STATICforIDL char *DESCRtoREF (ptr)
void	*ptr;		/* could point to a descriptor or to the string (if
			   already passed by reference) */
{
  char *aptr;		/* pointer to string */
  size_t len;		/* length of string */

  if (isDESCR (ptr, &aptr, &len))
    return aptr;
  else
    return (char *) ptr;	 	   /* already passed by reference */
}
#endif

#if defined(Fif_DESCR)
/******************************************************************************
 * DESCRtoREFnul.  Convert from passing by descriptor to passing by reference
 *		  and NUL-terminate.
 *  NOTE:
 *   If a descriptor is being pointed to, always use a temporary string to pass
 *   by reference - never assume that the string may already be NUL-terminated
 *   (in case its NUL-terminated by chance and the NUL might go away).
 *
 *****************************************************************************/

STATICforIDL char *DESCRtoREFnul (ptr, maxREFlen, ssh)
void	*ptr;				/* Could point to a descriptor or to
					   the string (if already passed by
					   reference). */
size_t	maxREFlen;			/* Maximum length of the string if
					   passed by reference. */
struct STRINGstruct **ssh;		/* Pointer to head of STRINGstruct
					   linked list. */
{
  char *aptr;		/* pointer to string */
  size_t len;		/* length of string */
  if (isDESCR (ptr, &aptr, &len))
    return NULterminate (aptr, len, ssh);
  else
    return NULterminate (ptr, maxREFlen, ssh);
}
#endif

#if defined(Fif_NOLEN)
/******************************************************************************
 * FindEndNUL.  Find the end of the string and NUL-terminate.
 *
 *  NOTE:
 *   Never assume that the string may already be NUL-terminated (in case its
 *   NUL-terminated by chance and the NUL might go away).
 *
 *****************************************************************************/

STATICforIDL char *FindEndNUL (ptr, maxREFlen, ssh)
char	*ptr;				/* Character string. */
size_t	maxREFlen;			/* Maximum length of the string. */
struct STRINGstruct **ssh;		/* Pointer to head of STRINGstruct
					   linked list. */
{
  size_t i;
  for (i = 0; i < maxREFlen; i++)
     if (!Printable(ptr[i])) return NULterminate (ptr, i, ssh);
  return NULterminate (ptr, maxREFlen, ssh);
}
#endif

/******************************************************************************
 * Trailer.
 *   Note that if the trailer is a null-string (length of zero), TRUE will be
 * returned.
 *****************************************************************************/

VISIBLE_PREFIX Logical Trailer (string, trailer)
char *string;
char *trailer;
{
  size_t stringL = strlen (string);
  size_t trailerL = strlen (trailer);
  if (trailerL > stringL) return FALSE;
  if (strcmp(&string[stringL-trailerL],trailer)) return FALSE;
  return TRUE;
}

/******************************************************************************
 * strcpyX.
 *    Copies from the source to the destination but only up to the maximum
 * number of characters specified.  Then NUL-terminates the destination
 * either at the actual end of the copied source string or after the maximum
 * number of characters.  If the maximum number of characters is zero, then
 * a normal `strcpy' is done. If the maximum number of characters is less than
 * zero, nothing is copied. 'max', if not zero, could be the destination's
 * actual size, or the maximum allowed size for the destination item defined
 * in CDF (when the actual size is not available in the calling function), not
 * including the terminating NUL character.
 *****************************************************************************/

VISIBLE_PREFIX char *strcpyX (dst, src, max)
char *dst;
char *src;
size_t max;
{
  if (max > 0 && src[0] != NUL) {
    size_t j, len;
    len = strlen(src);
    j = MINIMUM(max, len);
    memcpy (dst, src, j);
    dst[(int)j] = NUL;
  } else if (max == 0) {
    strcpy (dst, src);
  } else if (src[0] == NUL) {
    dst[0] = NUL;
  }
  return dst;
}

/******************************************************************************
 * strcatX.
 *    Concatenates from the source to the destination but only up to the
 * maximum number of characters, 'max', specified. In other words, 'max' is
 * used as a mechanism for controlling how many characters are allowed/stored
 * in the destination. Then NUL-terminates the destination either at the actual
 * end of the concatenated source string or after the maximum number of
 * characters.  If the maximum number of characters is zero, then a normal
 * `strcat' is done. If the maximum number of characters is less than zero, 
 * nothing is concatenated. 'max', if not zero, could be the destination's 
 * actual size, or the maximum allowed size for the destination item defined 
 * in CDF (when the actual size is not available in the calling function), not
 *including the terminating NUL character.
 *****************************************************************************/

VISIBLE_PREFIX char *strcatX (dst, src, max)
char *dst;
char *src;
size_t max;
{
  if (max > 0 && src[0] != NUL) {
    size_t i, j;
    i = strlen(dst);
    j = strlen(src);
    if (max > i) {
      size_t k = MINIMUM(max - i, j);
      if (k > 0) 
        strncat (dst, src, k);
    }
  } else if (max == 0) {
    strcat (dst, src);
  }
  return dst;
}

/******************************************************************************
 * catchrX.
 *    Concatenates one character (the source) to the destination string but
 * only up to the maximum number of characters specified.  If the maximum
 * number of characters is zero, the character is always concatenated.
 * If the maximum number of characters is less than zero, nothing is 
 * concatenated. 
 *****************************************************************************/

VISIBLE_PREFIX char *catchrX (dst, src, max)
char *dst;
int src;
size_t max;
{
  size_t l = strlen(dst);
  if (max > 0) {
    if (l < max) {
      dst[l] = (char) src;
      dst[l+1] = NUL;
    }
  } else if (max == 0) {
    dst[l] = (char) src;
    dst[l+1] = NUL;
  }
  return dst;
}

/******************************************************************************
 * catchrY.
 *    Concatenates a number of bytes (the source) to the destination string 
 * but only up to the maximum number of characters specified.  If the maximum
 * number of characters is zero, the character is always concatenated.
 * If the maximum number of characters is less than zero, nothing is 
 * concatenated. 
 *****************************************************************************/

VISIBLE_PREFIX char *catchrY (dst, src, nbytes, max)
char *dst;
char *src;
size_t nbytes;
size_t max;
{
  int ix;
  size_t l = strlen(dst);
  if (max > 0) {
    if ((l+nbytes) < max) {
      for (ix = 0; ix < nbytes; ++ix) {
        if ((l+ix) < max) 
          dst[l+ix] = (char) *(src+ix);
        else {
          dst[l+ix] = NUL;
          return dst;
        }
      }
      dst[l+nbytes] = NUL;
    }
  } else if (max == 0) {
    for (ix = 0; ix < nbytes; ++ix)
      dst[l+ix] = (char) *(src+ix);
    dst[l+nbytes] = NUL;
  }
  return dst;
}

/******************************************************************************
 * prependX.
 *    Place a string at the beginning of another string but don't let the
 * destination string exceed a maximum length.
 *****************************************************************************/

VISIBLE_PREFIX char *prependX (dst, src, max)
char *dst;
char *src;
size_t max;
{
  size_t lenDst = strlen(dst);
  size_t maxLen = max - lenDst;
  size_t lenSrc = MINIMUM (maxLen, strlen(src));
  memmove (&dst[lenSrc], dst, lenDst + 1);
  memmove (dst, src, lenSrc);
  return dst;
}

/******************************************************************************
 * PstrcpyX.
 *    Copy a Pascal-style string to a C-style string.
 *****************************************************************************/

#if defined(mac)
STATICforIDL char *PstrcpyX (dst, src, max)
char *dst;
char *src;
size_t max;
{
  if (max < 0) return dst;
  int len = (int) BOO(max > 0,MINIMUM(max,(int)src[0]),src[0]);
  strcpyX (dst, &src[1], len);
  dst[len] = NUL;
  return dst;
}
#endif

/******************************************************************************
* PstrcatX.
*    Concatinate a Pascal-style string to a C-style string.
******************************************************************************/

#if defined(mac)
STATICforIDL char *PstrcatX (dst, src, max)
char *dst;
char *src;
size_t max;
{
  if (max > 0) {
    int i, j;
    for (i = strlen(dst), j = 1; j <= src[0] && i < max; i++, j++) {
       dst[i] = src[j];
    }
    dst[i] = NUL;
  } else if (max == 0) {
     strncat (dst, &src[1], (size_t) src[0]);
  }
  return dst;
}
#endif

/******************************************************************************
 * PprependX.
 *    Place a Pascal-style string at the beginning of a C-style string but
 * don't let the C-style string exceed a maximum length.
******************************************************************************/

#if defined(mac)
STATICforIDL char *PprependX (dst, src, max)
char *dst;
char *src;
size_t max;
{
  size_t lenDst = strlen(dst);
  size_t maxLen = max - lenDst;
  size_t lenSrc = MINIMUM (maxLen, src[0]);
  memmove (&dst[lenSrc], dst, lenDst + 1);
  memmove (dst, &src[1], lenSrc);
  return dst;
}
#endif

/******************************************************************************
 * CtoFORTstring.
 *   Copies a NUL-terminated string (C-style) to a FORTRAN-style string.
 *****************************************************************************/

STATICforIDL void CtoFORTstring (Cstring, FORTstring, length)
char *Cstring;
void *FORTstring;
int length;
{
#if defined(Fif_DESCR)
  char *ptr; size_t lenT; int i;
  if (!isDESCR(FORTstring,&ptr,&lenT)) {
    ptr = FORTstring;
    lenT = (size_t) length;
  }
  strcpyX (ptr, Cstring, lenT);
  for (i = strlen(Cstring); i < lenT; i++) ptr[i] = ' ';
  return;
#else
  int i;
  strcpyX (FORTstring, Cstring, (size_t) length);
  for (i = strlen(Cstring); i < length; i++) ((char *) FORTstring)[i] = ' ';
  return;
#endif
}

/******************************************************************************
 * Printable.
 *****************************************************************************/

VISIBLE_PREFIX int Printable (c)
int c;
{
  int ret;
  return BOO(IsASCII(c),isprint(c),FALSE); 
/*
  if (IsASCII(c)) {
    ret = isprint(c);
    return (ret!=0?ret:FALSE);
  } else {
    wint_t wc = (wint_t) c;
    ret = iswprint(c);
    return (ret!=0?ret:FALSE);
  }
*/
}

/******************************************************************************
 * UpperCase.
 *****************************************************************************/

VISIBLE_PREFIX int UpperCase (c)
int c;
{
  int ret;
  if (IsASCII(c)) {
    ret = isupper(c);
    return (ret!=0?ret:FALSE);
  } else {
    wint_t wc = (wint_t) c;
    ret = iswupper(c);
    return (ret!=0?ret:FALSE);
  }
}

/******************************************************************************
 * LowerCase.
 *****************************************************************************/

VISIBLE_PREFIX int LowerCase (c)
int c;
{
  int ret;
  if (IsASCII(c)) {
    ret = islower(c);
    return (ret!=0?ret:FALSE);
  } else {
    wint_t wc = (wint_t) c;
    ret = iswlower(c);
    return (ret!=0?ret:FALSE);
  }
}

/******************************************************************************
 * Alphabetic.
 *****************************************************************************/

VISIBLE_PREFIX int Alphabetic (c)
int c;
{
  int ret;
  if (IsASCII(c)) {
    ret = isalpha(c);
    return (ret!=0?ret:FALSE);
  } else {
    wint_t wc = (wint_t) c;
    ret = iswalpha(c);
    return (ret!=0?ret:FALSE);
  }
}

/******************************************************************************
 * Decimal.
 *****************************************************************************/

VISIBLE_PREFIX int Decimal (c)
int c;
{
  int ret;
  if (IsASCII(c)) {
    ret = isdigit(c);
    return (ret!=0?ret:FALSE);
  } else {
    wint_t wc = (wint_t) c;
    ret = iswdigit(c);
    return (ret!=0?ret:FALSE);
  }
}

/******************************************************************************
 * Spacing.
 *****************************************************************************/

VISIBLE_PREFIX int Spacing (c)
int c;
{
  int ret;
  if (IsASCII(c)) {
    ret = isspace(c);
    return (ret!=0?ret:FALSE);
  } else {
    wint_t wc = (wint_t) c;
    ret = iswspace(c);
    return (ret!=0?ret:FALSE);
  }
}

/******************************************************************************
 * Visible.
 *****************************************************************************/

VISIBLE_PREFIX int Visible (c)
int c;
{
  int ret;
  if (IsASCII(c)) {
    ret = isgraph(c);
    return (ret!=0?ret:FALSE);
  } else {
    wint_t wc = (wint_t) c; 
    ret = iswgraph(c);
    return (ret!=0?ret:FALSE);
  }
}

/******************************************************************************
 * MakeLower.
 *****************************************************************************/

VISIBLE_PREFIX int MakeLower (c)
int c;
{
  if (IsASCII(c)) 
    return tolower(c);
  else {
    wint_t wc = (wint_t) c;
    return towlower(wc);
  }
}

/******************************************************************************
 * MakeUpper.
 *****************************************************************************/

VISIBLE_PREFIX int MakeUpper (c)
int c;
{
  if (IsASCII(c)) 
    return toupper(c);
  else {
    wint_t wc = (wint_t) c;
    return towupper(wc);
  }
}

/******************************************************************************
 * CDFstrcmpIgCase.
 * It returns  1 if compared strings are identical with the same length
 *             0 if compared strings are not same or there is a NULL string
 *****************************************************************************/

int CDFstrcmpIgCase (string1, string2)
char *string1;
char *string2;
{
  size_t string1L, string2L;
  int i;

  if (string1 == NULL || string2 == NULL) return 0;
  string1L = strlen(string1);
  string2L = strlen(string2);
  if (string1L != string2L) return 0;
  for (i = 0; i < (int) string1L; i++) {
     if (MakeLower(string1[i]) != MakeLower(string2[i]))
       return 0;
  }
  return 1;
}

/******************************************************************************
 * CDFstrstrIgCase.
 *****************************************************************************/

char *CDFstrstrIgCase (string, substring)
char *string;
char *substring;
{
  size_t stringL, substringL;
  int i, j;
  
  if (string == NULL || substring == NULL) return NULL;
  stringL = strlen(string);
  substringL = strlen(substring);
  if (stringL == 0 || substringL == 0) return NULL;
  if (stringL < substringL) return NULL;
  for (i = 0; i < (int) stringL; i++) {
     for (j = 0; j < (int) substringL; j++) {
        if (MakeLower(substring[j]) != MakeLower(string[i+j])) break;
     }
     if (j == substringL) return &string[i];
  }
  return NULL;
}

/******************************************************************************
 * UTF8StrChars - Find the number of UTF-8 characters in a given byte length
 * string. Returns -1 if this is not a valid UTF-8 string.
 *****************************************************************************/

VISIBLE_PREFIX int UTF8StrChars (const char *string, int byteCnt)
{
  int cnt; unsigned char *p;
  p = UTF8Check ((unsigned char *)string, byteCnt, &cnt);
  if (p != NULL)
    return -1;
  else
    return cnt;
}

/******************************************************************************
 * UTF8CharBytes - Get the byte count of a character. The argument should be a
 * pointer pointing to the starting byte location of a character in a
 * Unicode of UTF8 encoding. Otherwise, it returns 0.
 *****************************************************************************/
VISIBLE_PREFIX int UTF8CharBytes (const char *c)
{
  int cnt = 0;
  size_t len;
  if (c == NULL) return 0;
  len = strlen(c);
  if (len > 0) {
    cnt++;
    if ((*c & 0xc0) == 0xc0) { /* binary is 11xxxxxx */
      while( (c[cnt] & 0xc0) == 0x80  ) /* binary code is 10xxxxxx */
        cnt++;
    }
  }
  return cnt;
}

/******************************************************************************
 * UTF8CharAtBytes - Get the byte count of a given character in the 
 * null-terminating string. The character number is 0-based.  -1 is returned 
 * if the character number is greater than the total number of characters in 
 * the string.
 *****************************************************************************/
VISIBLE_PREFIX int UTF8CharAtBytes (const char *string, int charNum)
{
  int ix, iy, iz, loc;
  iz = (int) strlen (string);
  ix = UTF8StrChars (string, iz);
  if (charNum > ix) return -1;
  ix = -1;
  iy = loc = 0;
  while (ix < charNum) {
    iy = UTF8CharBytes (string+loc);
    ++ix;
    loc += iy;
  }
  return iy;
}

/******************************************************************************
 * UTF8CharAtOffset - find the byte offset from the beginning of the 
 * null-terminating string for a given character number (0-based). For
 * an ASCII string, it is the same as the character number. For multi-byte
 * character string, it is not. if the total number of characters is less than
 * the given character number, -1 is returned.
 *****************************************************************************/
VISIBLE_PREFIX int UTF8CharAtOffset (const char *string, int charNum)
{
  int ix, iy, i, j, len;
  iy = (int) strlen (string);
  ix = UTF8StrChars (string, iy);
  j = 0;
  if (charNum > ix) return -1;
  for (i = 0; i < charNum; ++i) {
    len = UTF8CharBytes (string+j);
    j += len;
  }
  return j;
}

/******************************************************************************
 * UTF8CharAtString - find the pointer location in the null-terminating string 
 * for a given character number (0-based). For an ASCII string, it is the same
 * as the character number. For multi-byte character string, it is not. if the
 * total number of characters is less than the given character number, NULL is 
 * returned.
******************************************************************************/
VISIBLE_PREFIX char *UTF8CharAtString (const char *string, int charNum)
{
  int ix, iy, i, j, len;
  iy = (int) strlen (string);
  ix = UTF8StrChars (string, iy);
  j = 0;
  if (charNum > ix) return NULL;
  for (i = 0; i < charNum; ++i) {
    len = UTF8CharBytes (string+j);
    j += len;
  }
  return ((char *)string)+j;
}

/******************************************************************************
 * UTF8ByteLengthFromCharAt - Get the remaining byte length in the 
 * null-terminating string from the given character number (0-based).
 *****************************************************************************/
VISIBLE_PREFIX int UTF8ByteLengthFromCharAt (const char *string, int charNum)
{
  int ix, iy, iz, loc;
  iz = strlen(string);
  if (charNum <= 0) return iz;
  ix = UTF8StrChars (string, iz);
  if (charNum > ix) return 0;
  ix = 0;
  iy = loc = 0;
  while (ix < charNum) {
    iy = UTF8CharBytes (string+loc);
    loc += iy;
    ++ix;
  }
  return (iz - loc);
}

/******************************************************************************
 * UTF8StrCharsFromCharAt - Find the remaining number of charatcers in the
 * null-terminating string from the given character number (0-based),
 * inclusive.
 *****************************************************************************/
VISIBLE_PREFIX int UTF8StrCharsFromCharAt (const char *string, int charNum)
{
   int ix, iy, i;
   iy = (int) strlen(string);
   ix = UTF8StrChars (string, iy);
   i = ix - charNum;
   if (i < 0) i = 0;
   return i;
}

/******************************************************************************
 * Encode Unicode to UTF8 string
 *  @param utf8 [Output] UTF8 character string (pre-allocated) passed in from
 *                       the caller. Make sure the size is big enough. It is
 *                       NUL-terminated when returns.
 *  @param codepoints [Input] An array of Unicode codepoints
 *  @param numCP [Input] The number of codepoints to convert
 *  @return The number of converted bytes in the output string (not including
 *          the NUL byte). If there is an error, a negative number is returned,
 *          indicating which codepoint is the bad one. The output will not be
 *          properly set.
 *
 * Modified from the original code der_encode_utf8_string.c by Tom St Denis
 * @file der_encode_utf8_string.c
 * ASN.1 DER, encode a UTF8 STRING, Tom St Denis
 * https://github.com/libtom/libtomcrypt
 * SPDX-License-Identifier: Unlicense
 *****************************************************************************/

VISIBLE_PREFIX int UnicodetoUTF8 (char *utf8, int *codepoints, int numCP)
{
  int i, bytePtr, size, cp;

  bytePtr = 0;
  /* store UTF8 */
  for (i = 0; i < numCP; i++) {
     cp = codepoints[i];
     if (cp <= 0x7F) size = 1;
     else if (cp <= 0x7FF) size = 2;
     else if (cp <= 0xFFFF) size = 3;
     else if (cp <= 0x10FFFF) size = 4;
     else size = -1;
     switch (size) {
        case 1:
           utf8[bytePtr++] = (char) cp;
           break;
        case 2:
           utf8[bytePtr++] = (char) (0xC0 | ((cp >> 6) & 0x1F));
           utf8[bytePtr++] = (char) (0x80 | (cp & 0x3F));
          break;
        case 3:
          utf8[bytePtr++] = (char) (0xE0 | ((cp >> 12) & 0x0F));
          utf8[bytePtr++] = (char) (0x80 | ((cp >> 6) & 0x3F));
          utf8[bytePtr++] = (char) (0x80 | (cp & 0x3F));
          break;
        case 4:
          utf8[bytePtr++] = (char) (0xF0 | ((cp >> 18) & 0x07));
          utf8[bytePtr++] = (char) (0x80 | ((cp >> 12) & 0x3F));
          utf8[bytePtr++] = (char) (0x80 | ((cp >> 6) & 0x3F));
          utf8[bytePtr++] = (char) (0x80 | (cp & 0x3F));
          break;
        default:
          return (-1-i);
     }
  }
  utf8[bytePtr] = (char) '\0';
  return bytePtr;
}

/******************************************************************************
 *  Decode a UTF8 string to Unicode codepoints.
 *  @param codepoint [out] The pre-allocated array for unicode codepoints
 *  @param utf8      [in] The encoded UTF8 STRING
 *  @param byteCnt   [in] The byte size of the UTF8 STRING
 *  @return          The number of Unicode codepoints if successful
 *                   Otherwise a negative number is returned, indicating which
 *                   byte position in the utf8 array is the bad one. The output
 *                   will not be properly set.
 *  Modified from the original code from der_decode_utf8_string.c
 *  by Tom St Denis
 *  @file der_decode_utf8_string.c
 *  ASN.1 DER, encode a UTF8 STRING, Tom St Denis
 *  https://github.com/libtom/libtomcrypt
 *  SPDX-License-Identifier: Unlicense
 *****************************************************************************/

int UTF8toUnicode(int *codepoints, const char *utf8, int byteCnt)
{
  int x, y, z, kk, tmp;

  if (byteCnt == 0) return 0;
  x = 0;
  for (y = 0; x < byteCnt; ) {
     /* read first byte */
     kk = x;
     tmp = utf8[x++];

     /* A unicode character is recovered from a sequence of 1 to 4 utf8 bytes.
      * the form of those bytes must match a row in the following table:
      *    0xxxxxxx
      *    110xxxxx 10xxxxxx
      *    1110xxxx 10xxxxxx 10xxxxxx
      *    11110xxx 10xxxxxx 10xxxxxx 10xxxxxx
      * the number of leading 1s in the first byte (0,2,3,4) determines the
      * number of remaining bytes to read (0,1,2,3)
      */

     /* determine z, the number of leading 1s.
        this is done by left-shifting tmp, which clears the ms-bits */
     for (z = 0; (tmp & 0x80) && (z <= 4); z++, tmp = (tmp << 1) & 0xFF);

     /* z should be in {0,2,3,4} */
     if (z == 1 || z > 4) return (-1-kk);

     /* right-shift tmp to restore least-sig bits */
     tmp >>= z;

     /* now update z so it equals the number of additional bytes to read */
     if (z > 0) --z;

     if ((x + z) > byteCnt) return (-1-kk);

     /* read remaining bytes */
     while (z-- != 0) {
        if ((utf8[x] & 0xC0) != 0x80) return (-1-kk);
        tmp = (tmp << 6) | ((int)utf8[x++] & 0x3F);
     }

     codepoints[y++] = tmp;
  }

  return y;
}

/******************************************************************************
 *  Checks if a character string is a UTF8 encoding.
 *  @param s         [in] The terminating character string.
 *  @param numBytes  [in] The number of bytes in the string to check. Make
 *                        sure the last bytes cover a full UTF-8 character.
 *                        If 0, all bytes are checked (make sure that the
 *                        string is NUL-terminated).
 *  @param charCnt   [out] The number of UTF-8 characters in the string if the
 *                         count is requested. Otherwise, pass in NULL.
 *  @return          NULL if the string passes the check. Otherwise, a pointer 
 *                   to the postion in the string that fails the check. 
 *
 * Modified from the original code from Markus Kuhn.
 *   Markus Kuhn <http://www.cl.cam.ac.uk/~mgk25/> -- 2005-03-30
 *   License: http://www.cl.cam.ac.uk/~mgk25/short-license.html
 *
 * Added the number of bytes to check in the string. If it is 0, the full
 * string is checked. 
 *****************************************************************************/
VISIBLE_PREFIX unsigned char *UTF8Check(unsigned char *s, int numBytes,
	                         	int* charCnt)
{
  int i;
  i = 0;
  if (charCnt != NULL) *charCnt = 0;
  while (*s) {
    if ((int) (unsigned char)*s < 0x80) {
      /* 1 byte: 0xxxxxxx */
      s++;
      i++;
      if (charCnt != NULL) *charCnt = *charCnt + 1;
      if (numBytes > 0 && i == numBytes) break;
    } else if ((s[0] & 0xe0) == 0xc0) {
      /* 2 bytes: 110XXXXx 10xxxxxx */
      if ((numBytes > 0 && (numBytes-i)<= 1) ||
          ((int) (unsigned char)(s[1] & 0xc0) != 0x80 ||
           (int) (unsigned char)(s[0] & 0xfe) == 0xc0)) { /* overlong? */
        return s;
      } else {
        s += 2;
        i += 2;
	if (charCnt != NULL) *charCnt = *charCnt + 1;
      }
      if (numBytes > 0 && i == numBytes) break;
    } else if ((int) (unsigned char)(s[0] & 0xf0) == 0xe0) {
      /* 3 bytes: 1110XXXX 10Xxxxxx 10xxxxxx */
      if ((numBytes > 0 && (numBytes-i)<= 2) ||
          ((int) (unsigned char)(s[1] & 0xc0) != 0x80 ||
           (int) (unsigned char)(s[2] & 0xc0) != 0x80 ||
           ((int) (unsigned char) s[0] == 0xe0 &&
            (int) (unsigned char)(s[1] & 0xe0) == 0x80) || /* overlong? */
           ((int) (unsigned char) s[0] == 0xed &&
            (int) (unsigned char)(s[1] & 0xe0) == 0xa0) || /* surrogate? */
           ((int) (unsigned char) s[0] == 0xef &&
            (int) (unsigned char) s[1] == 0xbf &&
            (int) (unsigned char)(s[2] & 0xfe) == 0xbe))) {/*U+FFFE or U+FFFF? */
        return s;
      } else {
        s += 3;
        i += 3;
	if (charCnt != NULL) *charCnt = *charCnt + 1;
      }
      if (numBytes > 0 && i == numBytes) break;
    } else if ((s[0] & 0xf8) == 0xf0) {
      /* 4 bytes: 11110XXX 10XXxxxx 10xxxxxx 10xxxxxx */
      if ((numBytes > 0 && (numBytes-i)<= 3) ||
          (int) (unsigned char)(s[1] & 0xc0) != 0x80 ||
          (int) (unsigned char)(s[2] & 0xc0) != 0x80 ||
          (int) (unsigned char)(s[3] & 0xc0) != 0x80 ||
          ((int) (unsigned char) s[0] == 0xf0 &&
           (int) (unsigned char)(s[1] & 0xf0) == 0x80) ||  /* overlong? */
          ((int) (unsigned char) s[0] == 0xf4 &&
           (int) (unsigned char) s[1] > 0x8f) ||
          (int) (unsigned char) s[0] > 0xf4) {             /* > U+10FFFF? */
        return s;
      } else {
        s += 4;
        i += 4;
	if (charCnt != NULL) *charCnt = *charCnt + 1;
      }
      if (numBytes > 0 && i == numBytes) break;
    } else {
      return s;
    }
  }
  return NULL;
}

/******************************************************************************
 *  Checks if a char value is a UTF8 encoding.
 *  @param value     [in] The value to check
 *  @return          1 if passed the test
 *                   0 otherwise
 *
 * Modified from the original code from Markus Kuhn.
 *   Markus Kuhn <http://www.cl.cam.ac.uk/~mgk25/> -- 2005-03-30
 *   License: http://www.cl.cam.ac.uk/~mgk25/short-license.html
 *****************************************************************************/
VISIBLE_PREFIX int IsUTF8Char(const char *s)
{
  if ((int)(unsigned char)*s < 0x80) {
    /* 0xxxxxxx */
    if (Printable((int)*s))
      return 1;
    else 
      return 0;
  } else if ((int) (unsigned char) (s[0] & 0xe0) == 0xc0) {
    /* 110XXXXx 10xxxxxx */
    if ((int) (unsigned char) (s[1] & 0xc0) != 0x80 ||
        (int) (unsigned char) (s[0] & 0xfe) == 0xc0)   /* overlong? */
      return 0;
    else
      return 1;
  } else if ((int) (unsigned char) (s[0] & 0xf0) == 0xe0) {
    /* 1110XXXX 10Xxxxxx 10xxxxxx */
    if ((int) (unsigned char) (s[1] & 0xc0) != 0x80 ||
        (int) (unsigned char) (s[2] & 0xc0) != 0x80 ||
        ((int) (unsigned char) s[0] == 0xe0 &&
         (int) (unsigned char) (s[1] & 0xe0) == 0x80) ||/* overlong? */
        ((int) (unsigned char) s[0] == 0xed &&
         (int) (unsigned char) (s[1] & 0xe0) == 0xa0) ||/* surrogate? */
        ((int) (unsigned char) s[0] == 0xef &&
         (int) (unsigned char) s[1] == 0xbf &&
         (int) (unsigned char) (s[2] & 0xfe) == 0xbe))  /* U+FFFE or U+FFFF? */
      return 0;
    else
      return 1;
  } else if ((int) (unsigned char) (s[0] & 0xf8) == 0xf0) {
    /* 11110XXX 10XXxxxx 10xxxxxx 10xxxxxx */
    if ((int) (unsigned char) (s[1] & 0xc0) != 0x80 ||
        (int) (unsigned char) (s[2] & 0xc0) != 0x80 ||
        (int) (unsigned char) (s[3] & 0xc0) != 0x80 ||
        ((int) (unsigned char) s[0] == 0xf0 &&
         (int) (unsigned char) (s[1] & 0xf0) == 0x80) ||/* overlong? */
        ((int) (unsigned char) s[0] == 0xf4 &&
         (int) (unsigned char) s[1] > 0x8f) ||
        (int) (unsigned char) s[0] > 0xf4)              /* > U+10FFFF? */
      return 0;
    else
      return 1;;
  } else
    return 0;
}

/******************************************************************************
* `strstr' (for SunOS before 4.1).  This should no longer be needed.  If it
* is, STATICforIDL may or may not need to be present.
******************************************************************************/

#if defined(SunOSpre41)
/*STATICforIDL*/ char *strstr (string, substring)
char *string;
char *substring;
{
  int stringL = strlen(string);
  int substringL = strlen(substring);
  int i, j;
  if (stringL == 0 || substringL == 0) return NULL;
  if (stringL < substringL) return NULL;
  for (i = 0; i < stringL; i++) {
     for (j = 0; j < substringL; j++)
	if (substring[j] != string[i+j]) break;
     if (j == substringL) return &string[i];
  }
  return NULL;
}
#endif
