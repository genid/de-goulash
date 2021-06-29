/*  plugins/dosage.c -- prints genotype dosage.

    Copyright (C) 2014-2018 Genome Research Ltd.

    Author: Petr Danecek <pd3@sanger.ac.uk>

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
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#include <stdio.h>
#include <stdlib.h>
#include <htslib/vcf.h>
#include <math.h>
#include <getopt.h>
#include <inttypes.h>
#include "bcftools.h"


/*
    This short description is used to generate the output of `bcftools plugin -l`.
*/
const char *about(void)
{
    return "Prints genotype dosage determined from tags requested by the user.\n";
}

const char *usage(void)
{
    return 
        "\n"
        "About: Print genotype dosage\n"
        "Usage: bcftools +dosage [General Options] -- [Plugin Options]\n"
        "Options:\n"
        "   run \"bcftools plugin\" for a list of common options\n"
        "\n"
        "Plugin options:\n"
        "   -t, --tags <list>   VCF tags to determine the dosage from [PL,GL,GT]\n"
        "\n"
        "Example:\n"
        "   bcftools +dosage in.vcf -- -t GT\n"
        "\n";
}

bcf_hdr_t *in_hdr = NULL;
int pl_type = 0, gl_type = 0;
uint8_t *buf = NULL;
int nbuf = 0;   // NB: number of elements, not bytes
char **tags = NULL;
int ntags = 0;
float *vals = NULL, *dsg = NULL;
int mvals, mdsg;

typedef int (*dosage_f) (bcf1_t *);
dosage_f *handlers = NULL;
int nhandlers = 0;


int calc_dosage_PL(bcf1_t *rec)
{
    int i, j, nret = bcf_get_format_values(in_hdr,rec,"PL",(void**)&buf,&nbuf,pl_type);
    if ( nret<0 ) return -1;

    nret /= rec->n_sample;
    if ( nret != rec->n_allele*(rec->n_allele+1)/2 ) return -1;     // not diploid
    hts_expand(float, nret, mvals, vals);
    hts_expand(float, rec->n_allele, mdsg, dsg);
    #define BRANCH(type_t,is_missing,is_vector_end) \
    { \
        type_t *ptr = (type_t*) buf; \
        for (i=0; i<rec->n_sample; i++) \
        { \
            float sum = 0; \
            for (j=0; j<nret; j++) \
            { \
                if ( is_missing || is_vector_end ) break; \
                vals[j] = pow(10,-0.1*ptr[j]); \
                sum += vals[j]; \
            } \
            if ( j<nret ) \
                for (j=0; j<rec->n_allele; j++) dsg[j] = -1; \
            else \
            { \
                if ( sum ) for (j=0; j<nret; j++) vals[j] /= sum; \
                vals[0] = 0; \
                memset(dsg, 0, sizeof(float)*rec->n_allele); \
                int k, l = 0; \
                for (j=0; j<rec->n_allele; j++) \
                { \
                    for (k=0; k<=j; k++) \
                    { \
                        dsg[j] += vals[l]; \
                        dsg[k] += vals[l]; \
                        l++; \
                    } \
                } \
            } \
            for (j=1; j<rec->n_allele; j++) \
                printf("%c%f",j==1?'\t':',',dsg[j]); \
            ptr += nret; \
        } \
    }
    switch (pl_type)
    {
        case BCF_HT_INT:  BRANCH(int32_t,ptr[j]==bcf_int32_missing,ptr[j]==bcf_int32_vector_end); break;
        case BCF_HT_REAL: BRANCH(float,bcf_float_is_missing(ptr[j]),bcf_float_is_vector_end(ptr[j])); break;
    }
    #undef BRANCH
    return 0;
}

int calc_dosage_GL(bcf1_t *rec)
{
    int i, j, nret = bcf_get_format_values(in_hdr,rec,"GL",(void**)&buf,&nbuf,gl_type);
    if ( nret<0 ) return -1;

    nret /= rec->n_sample;
    if ( nret != rec->n_allele*(rec->n_allele+1)/2 ) return -1;     // not diploid
    hts_expand(float, nret, mvals, vals);
    hts_expand(float, rec->n_allele, mdsg, dsg);
    #define BRANCH(type_t,is_missing,is_vector_end) \
    { \
        type_t *ptr = (type_t*) buf; \
        for (i=0; i<rec->n_sample; i++) \
        { \
            float sum = 0; \
            for (j=0; j<nret; j++) \
            { \
                if ( is_missing || is_vector_end ) break; \
                vals[j] = pow(10,ptr[j]); \
                sum += vals[j]; \
            } \
            if ( j<nret ) \
                for (j=0; j<rec->n_allele; j++) dsg[j] = -1; \
            else \
            { \
                if ( sum ) for (j=0; j<nret; j++) vals[j] /= sum; \
                vals[0] = 0; \
                memset(dsg, 0, sizeof(float)*rec->n_allele); \
                int k, l = 0; \
                for (j=0; j<rec->n_allele; j++) \
                { \
                    for (k=0; k<=j; k++) \
                    { \
                        dsg[j] += vals[l]; \
                        dsg[k] += vals[l]; \
                        l++; \
                    } \
                } \
            } \
            for (j=1; j<rec->n_allele; j++) \
                printf("%c%f",j==1?'\t':',',dsg[j]); \
            ptr += nret; \
        } \
    }
    switch (gl_type)
    {
        case BCF_HT_INT:  BRANCH(int32_t,ptr[j]==bcf_int32_missing,ptr[j]==bcf_int32_vector_end); break;
        case BCF_HT_REAL: BRANCH(float,bcf_float_is_missing(ptr[j]),bcf_float_is_vector_end(ptr[j])); break;
    }
    #undef BRANCH
    return 0;
}

int calc_dosage_GT(bcf1_t *rec)
{
    int i, j, nret = bcf_get_genotypes(in_hdr,rec,(void**)&buf,&nbuf);
    if ( nret<0 ) return -1;

    nret /= rec->n_sample;
    int32_t *ptr = (int32_t*) buf;
    hts_expand(float, rec->n_allele, mdsg, dsg);
    for (i=0; i<rec->n_sample; i++)
    {
        memset(dsg, 0, sizeof(float)*rec->n_allele);
        for (j=0; j<nret; j++)
        {
            if ( ptr[j]==bcf_int32_vector_end || bcf_gt_is_missing(ptr[j]) ) break;
            int idx = bcf_gt_allele(ptr[j]);
            if ( idx > rec->n_allele ) error("The allele index is out of range at %s:%"PRId64"\n", bcf_seqname(in_hdr,rec),(int64_t) rec->pos+1);
            dsg[idx] += 1;
        }
        if ( !j )
            for (j=0; j<rec->n_allele; j++) dsg[j] = -1;
        for (j=1; j<rec->n_allele; j++)
            printf("%c%.1f",j==1?'\t':',',dsg[j]);
        ptr += nret;
    }
    return 0;
}


char **split_list(char *str, int *nitems)
{
    int n = 0, done = 0;
    char *ss = strdup(str), **out = NULL;
    while ( !done && *ss )
    {
        char *se = ss;
        while ( *se && *se!=',' ) se++;
        if ( !*se ) done = 1;
        *se = 0;
        n++;
        out = (char**) realloc(out,sizeof(char*)*n);
        out[n-1] = ss;
        ss = se+1;
    }
    *nitems = n;
    return out;
}

int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out)
{
    int i, id, c;
    char *tags_str = "PL,GL,GT";

    static struct option loptions[] =
    {
        {"tags",1,0,'t'},
        {0,0,0,0}
    };
    while ((c = getopt_long(argc, argv, "t:?h",loptions,NULL)) >= 0)
    {
        switch (c) 
        {
            case 't': tags_str = optarg; break;
            case 'h':
            case '?':
            default: fprintf(stderr,"%s", usage()); exit(1); break;
        }
    }
    tags = split_list(tags_str, &ntags);

    in_hdr = in;
    for (i=0; i<ntags; i++)
    {
        if ( !strcmp("PL",tags[i]) )
        {
            id = bcf_hdr_id2int(in_hdr,BCF_DT_ID,"PL");
            if ( bcf_hdr_idinfo_exists(in_hdr,BCF_HL_FMT,id) )
            {
                pl_type = bcf_hdr_id2type(in_hdr,BCF_HL_FMT,id);
                if ( pl_type!=BCF_HT_INT && pl_type!=BCF_HT_REAL )
                {
                    fprintf(stderr,"Expected numeric type of FORMAT/PL\n");
                    return -1;
                }
                handlers = (dosage_f*) realloc(handlers,(nhandlers+1)*sizeof(*handlers));
                handlers[nhandlers++] = calc_dosage_PL;
            }
        }
        else if ( !strcmp("GL",tags[i]) )
        {
            id = bcf_hdr_id2int(in_hdr,BCF_DT_ID,"GL");
            if ( bcf_hdr_idinfo_exists(in_hdr,BCF_HL_FMT,id) )
            {
                gl_type = bcf_hdr_id2type(in_hdr,BCF_HL_FMT,id);
                if ( gl_type!=BCF_HT_INT && gl_type!=BCF_HT_REAL )
                {
                    fprintf(stderr,"Expected numeric type of FORMAT/GL\n");
                    return -1;
                }
                handlers = (dosage_f*) realloc(handlers,(nhandlers+1)*sizeof(*handlers));
                handlers[nhandlers++] = calc_dosage_GL;
            }
        }
        else if ( !strcmp("GT",tags[i]) )
        {
            handlers = (dosage_f*) realloc(handlers,(nhandlers+1)*sizeof(*handlers));
            handlers[nhandlers++] = calc_dosage_GT;
        }
        else
        {
            fprintf(stderr,"No handler for tag \"%s\"\n", tags[i]);
            return -1;
        }
    }
    free(tags[0]);
    free(tags);

    printf("#[1]CHROM\t[2]POS\t[3]REF\t[4]ALT");
    for (i=0; i<bcf_hdr_nsamples(in_hdr); i++) printf("\t[%d]%s", i+5,in_hdr->samples[i]);
    printf("\n");

    return 1;
}


bcf1_t *process(bcf1_t *rec)
{
    int i,j, ret;

    printf("%s\t%"PRId64"\t%s", bcf_seqname(in_hdr,rec),(int64_t) rec->pos+1,rec->d.allele[0]);
    if ( rec->n_allele == 1 ) printf("\t.");
    else for (i=1; i<rec->n_allele; i++) printf("%c%s", i==1?'\t':',', rec->d.allele[i]);
    if ( rec->n_allele==1 )
    {
        for (j=0; j<rec->n_sample; j++) printf("\t0.0");
    }
    else
    {
        for (i=0; i<nhandlers; i++)
        {
            ret = handlers[i](rec);
            if ( !ret ) break;  // successfully printed
        }
        if ( i==nhandlers )
        {
            // none of the annotations present
            for (i=0; i<rec->n_sample; i++) printf("\t-1.0");
        }
    }
    printf("\n");

    return NULL;
}


void destroy(void)
{
    free(vals);
    free(dsg);
    free(handlers);
    free(buf);
}


