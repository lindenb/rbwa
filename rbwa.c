/*
The MIT License (MIT)

Copyright (c) 2014 Pierre Lindenbaum PhD.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

*/
#include <math.h>
#include <stdio.h>
#include <zlib.h>
#include <string.h>
#include <errno.h>
#include <assert.h>
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include "bwamem.h"
#include "kseq.h"

KSEQ_DECLARE(gzFile)

#define WHERENL fprintf(stderr,"[%s:%d] ",__FILE__,__LINE__)
#define DIE_FAILURE(FormatLiteral,...) do { WHERENL; fprintf (stderr,"Exiting: " FormatLiteral "\n", ##__VA_ARGS__); exit(EXIT_FAILURE);} while(0)
#define DEBUG(FormatLiteral, ...) do { fputs("[DEBUG]",stderr); WHERENL; fprintf (stderr,"" FormatLiteral "\n", ##__VA_ARGS__);} while(0)


/**
 * Structure holding the BWA index and the options
 */ 
typedef struct rbwa_t
	{
	bwaidx_t *idx;
	mem_opt_t *opt;
	}RBwaHandler,*RBwaHandlerPtr;

/**
 * Close resources associated to bwa
 */
SEXP RBwaClose(SEXP handle)
	{	
	void *p = R_ExternalPtrAddr(handle);
	if(p!=NULL)
		{
		RBwaHandlerPtr h=(RBwaHandlerPtr)p;
		bwa_idx_destroy(h->idx);
		free(h->opt);
		free(h);
		}
	R_ClearExternalPtr(handle);
	return ScalarLogical(1);
}

/**
 * destructor
 * 
 */
static void _RBwaFinalizer(SEXP handle)
	{
	RBwaClose(handle);
	}


/**
 * Opeb bwa resources
 */
SEXP RBwaOpen(SEXP Rfilename)
	{
	RBwaHandlerPtr handler;
	const char* filename= CHAR(STRING_ELT(Rfilename, 0));
	if(filename==NULL) DIE_FAILURE("NULL Filename");
	handler=(RBwaHandlerPtr)calloc(1,sizeof(RBwaHandler));
	if(handler==NULL)
		{
		DIE_FAILURE("out of memory");
		}
	/* load BWA index */
	handler->idx = bwa_idx_load(filename, BWA_IDX_ALL); 
	if(handler->idx==NULL)
		{
		DIE_FAILURE("Cannot load index \"%s\"",filename);
		}
	/* load BWA option */
	handler->opt = mem_opt_init();
	/* wrap pointer in SEXP */
	SEXP ext = PROTECT(R_MakeExternalPtr(handler, R_NilValue, R_NilValue));
	/* register destructor */
	R_RegisterCFinalizerEx(ext,_RBwaFinalizer, TRUE);
     	UNPROTECT(1);
	return ext;
	}
/**
 * map a DNA sequence with BWA
 * http://stackoverflow.com/questions/23547625
 */
SEXP RBwaMap(SEXP handle,SEXP readseqR)
	{
	SEXP res, vChrom,vPos,vStrand,vMapq,vSecondary,vNM,cls,rownam;
	int nprotect=0;
	RBwaHandlerPtr handler;
	mem_alnreg_v ar;
	int i;
	kseq_t ks;
	const char* seq= CHAR(STRING_ELT(readseqR, 0));
	/* check DNA sequence */
	if(seq==NULL) DIE_FAILURE("NULL seq");
	/* retrieve and cast handler */
	handler=(RBwaHandlerPtr)R_ExternalPtrAddr(handle);
	
	/* initialize short read */	
	memset((void*)&ks,0,sizeof(kseq_t));
	ks.seq.l=strlen(seq);
	ks.seq.s=(char*)seq;
	
	/* run allign */
	ar = mem_align1(handler->opt,
		handler->idx->bwt,
		handler->idx->bns,
		handler->idx->pac,
		ks.seq.l,
		ks.seq.s
		); // get all the hits
	

	/* prepare the table and its columns */
	PROTECT(res = Rf_allocVector(VECSXP,6));nprotect++;
	PROTECT(vChrom = Rf_allocVector(VECSXP,ar.n)); nprotect++;
   	PROTECT(vPos = Rf_allocVector(VECSXP, ar.n)); nprotect++;
  	PROTECT(vStrand = Rf_allocVector(VECSXP,ar.n));nprotect++;
	PROTECT(vMapq = Rf_allocVector(VECSXP, ar.n));nprotect++;
	PROTECT(vNM = Rf_allocVector(VECSXP, ar.n));nprotect++;
	PROTECT(vSecondary = Rf_allocVector(VECSXP, ar.n));nprotect++;
	
	/** set the table as a data.frame */
	PROTECT(cls = allocVector(STRSXP, 1)); nprotect++;
  	SET_STRING_ELT(cls, 0, mkChar("data.frame"));
   	classgets(res, cls);
	
	

	/* set the columns of the table */
	SET_VECTOR_ELT(res, 0, vChrom);
  	SET_VECTOR_ELT(res, 1, vPos);
	SET_VECTOR_ELT(res, 2, vStrand);
	SET_VECTOR_ELT(res, 3, vMapq);
	SET_VECTOR_ELT(res, 4, vNM);
	SET_VECTOR_ELT(res, 5, vSecondary);

	/* set the columns' name of the table */
	SEXP sNames = PROTECT(allocVector(STRSXP, 6));nprotect++;
	SET_STRING_ELT(sNames, 0, mkChar("chrom"));
	SET_STRING_ELT(sNames, 1, mkChar("pos"));
	SET_STRING_ELT(sNames, 2, mkChar("strand"));
	SET_STRING_ELT(sNames, 3, mkChar("mapq"));
	SET_STRING_ELT(sNames, 4, mkChar("NM"));
	SET_STRING_ELT(sNames, 5, mkChar("secondary"));
	setAttrib(res, R_NamesSymbol, sNames);

	/* loop over the hits */
	for (i = 0; i < ar.n; ++i)
		{
		mem_aln_t a = mem_reg2aln(handler->opt, handler->idx->bns, handler->idx->pac, ks.seq.l, ks.seq.s, &ar.a[i]);
		SET_VECTOR_ELT(vChrom, i , mkString(handler->idx->bns->anns[a.rid].name));
		SET_VECTOR_ELT(vPos, i ,ScalarInteger((long)a.pos));
		SET_VECTOR_ELT(vStrand, i ,ScalarInteger((int)a.is_rev));
		SET_VECTOR_ELT(vMapq, i ,ScalarInteger((int)a.mapq));
		SET_VECTOR_ELT(vNM, i ,ScalarInteger((int)a.NM));
		SET_VECTOR_ELT(vSecondary, i ,ScalarInteger((int)(ar.a[i].secondary >= 0)));
		free(a.cigar);
		}	
	free(ar.a);
		
	
	/** set the row names */
	PROTECT(rownam = allocVector(STRSXP, ar.n));nprotect++; // row.names attribute
	for (i = 0; i < ar.n; ++i)
		{
		char rname[20];
		sprintf(rname,"%d",i+1);
		SET_STRING_ELT(rownam, i ,mkChar(rname));
		}
	setAttrib(res, R_RowNamesSymbol, rownam);

	

	UNPROTECT(nprotect);
	return res;
	}

