// sequtils: useful sequence things

#include "Python.h"

#include <stdio.h>


//Translation Table
const char DNA2AMINO[] = {
    'F',   //TTT
    'F',   //TTC
    'L',   //TTA
    'L',   //TTG

    'S',   //TCT
    'S',   //TCC
    'S',   //TCA
    'S',   //TCG

    'Y',   //TAT
    'Y',   //TAC
    '*',   //TAA
    '*',   //TAG

    'C',   //TGT
    'C',   //TGC
    '*',   //TGA
    'W',   //TGG

    'L',   //CTT
    'L',   //CTC
    'L',   //CTA
    'L',   //CTG

    'P',   //CCT
    'P',   //CCC
    'P',   //CCA
    'P',   //CCG

    'H',   //CAT
    'H',   //CAC
    'Q',   //CAA
    'Q',   //CAG

    'R',   //CGT
    'R',   //CGC
    'R',   //CGA
    'R',   //CGG

    'I',   //ATT
    'I',   //ATC
    'I',   //ATA
    'M',   //ATG

    'T',   //ACT
    'T',   //ACC
    'T',   //ACA
    'T',   //ACG

    'N',   //AAT
    'N',   //AAC
    'K',   //AAA
    'K',   //AAG

    'S',   //AGT
    'S',   //AGC
    'R',   //AGA
    'R',   //AGG

    'V',   //GTT
    'V',   //GTC
    'V',   //GTA
    'V',   //GTG

    'A',   //GCT
    'A',   //GCC
    'A',   //GCA
    'A',   //GCG

    'D',   //GAT
    'D',   //GAC
    'E',   //GAA
    'E',   //GAG

    'G',   //GGT
    'G',   //GGC
    'G',   //GGA
    'G',   //GGG
};

const char GENERIC_DNA[] = "gatcrywsmkhbvdn";
const char GENERIC_RNA[] = "gaucrywsmkhbvdn";
const char GENERIC_AMINO[] = "acdefghiklmnpqrstvwy";

const unsigned int DNA_LEN = 15;
const unsigned int RNA_LEN = 15;
const unsigned int AMINO_LEN = 20;


inline int baseAsInt(char base, unsigned int* out)
{
    if(base < 'a')
        base = base - 'A' + 'a';

    if(base == 't')
    {
        *out = 0;
        return 1;
    }
    else if(base == 'c')
    {
        *out = 1;
        return 1;
    }
    else if(base == 'a')
    {
        *out = 2;
        return 1;
    }
    else if(base == 'g')
    {
        *out = 3;
        return 1;
    }
    else 
    {
        unsigned int i = 0;
        for(i = 0; i < DNA_LEN; i++)
        {
            if(base == GENERIC_DNA[i])
            {
                *out = 4;
                return 1;
            }
        }
    }
    
    return 0;
}

void clear(char** oseq)
{
    unsigned int i;
    for(i = 0; i < 6; i++)
        free(oseq[i]);

    free(oseq);
}

void reverse(char* c, unsigned int len)
{
    char tmp;
    unsigned int i;
    for(i = 0; i < len/2; i++)
    {
        tmp = c[i];
        c[i] = c[len-i-1];
        c[len-i-1] = tmp;
    }
}

inline void 
translate_codon(unsigned int a, unsigned int  b, unsigned int c, char* out)
{
    if( (a > 3) || (b > 3) || (c > 3))
    {
        *out = 'X';
        return;
    }
    *out = DNA2AMINO[16*a + 4 * b + c];
}

void unknown_base(char base, unsigned int position)
{
    char* err = malloc(256);
    sprintf(err, "unknown base \'%c\' at position %d", base, position);
    PyErr_SetString(PyExc_ValueError, err);
}

static PyObject *
translate(PyObject *self, PyObject *args)
{
    const char * iseq = NULL;
    unsigned int ilen = 0;

    if(!PyArg_ParseTuple(args, "s#:translate", &iseq, &ilen))
    {
        return NULL;
    }

    if(ilen < 3)
    {
        PyErr_SetString(PyExc_ValueError, "sequence too short");
        return NULL;
    }

    unsigned int olen = ilen / 3;
    unsigned int i;

    char * oseq = malloc(olen * sizeof(char)+1);
    oseq[olen] = '\0';

    unsigned int a = 0, b = 0, c = 0;
    

    for(i = 0; i < ilen-2; i+=3)
    {
        if(!baseAsInt(iseq[i], &a))
        {
            unknown_base(iseq[i], i);
            free(oseq);
            return NULL;
        }
        if(!baseAsInt(iseq[i + 1], &b))
        {
            unknown_base(iseq[i+1], i+1);
            free(oseq);
            return NULL;
        }
        if(!baseAsInt(iseq[i + 2], &c))
        {
            unknown_base(iseq[i+2], i+2);
            free(oseq);
            return NULL;
        }

        translate_codon(a, b, c, oseq + i/3);

    }

    return Py_BuildValue("s", oseq);
}

static PyObject *
sixFrameTranslation(PyObject *self, PyObject *args)
{
    const char * iseq = NULL;
    unsigned int ilen = 0;

    if(!PyArg_ParseTuple(args, "s#:sixFrameTranslation", &iseq, &ilen))
    {
        return NULL;
    }

    if(ilen < 3)
    {
        PyErr_SetString(PyExc_ValueError, "sequence too short");
        return NULL;
    }

    unsigned int olen = ilen / 3;
    unsigned int i;

    char ** oseq = malloc(6 * sizeof(char*));

    for(i = 0; i < 6; i++)
    {
        oseq[i] = malloc(olen * sizeof(char) + 1);
        oseq[i][olen] = '\0';
    }


    //perform the translations
    unsigned int f = 0, j = 0;
    unsigned int a = 0, b = 0, c = 0;
    if(!baseAsInt(iseq[0], &a))
    {
        unknown_base(iseq[0], 0);
        clear(oseq);
        return NULL;
    }
    if(!baseAsInt(iseq[1], &b))
    {
        unknown_base(iseq[1], 1);
        clear(oseq);
        return NULL;
    }

    for(i = 0; i < ilen-2; i++)
    {
        if(!baseAsInt(iseq[i + 2], &c))
        {
            unknown_base(iseq[i+2], i+2);
            clear(oseq);
            return NULL;
        }
        //forward frame
        translate_codon(a, b, c, oseq[f] + j);
        //reverse frame
        translate_codon((c+2)%4, (b+2)%4, (a+2)%4, oseq[f+3] + j);

        f++;
        if(f >= 3)
        {
            f = 0;
            j++;
        }

        a = b;
        b = c;
    }

    //make sure the terminating nulls are in the right place
    for(f = 0; f < 3; f++)
    {
        oseq[f  ][(ilen-f)/3] = '\0';
        oseq[f+3][(ilen-f)/3] = '\0';
    }

    //reverse the direction of the three reverse frames
    reverse(oseq[3], olen);
    reverse(oseq[4], (ilen-1)/3);
    reverse(oseq[5], (ilen-2)/3);

    switch(ilen % 3)
    {
        case 0:
            return Py_BuildValue("ssssss", oseq[0], oseq[1], oseq[2], 
                                oseq[3], oseq[5], oseq[4]);
        case 1:
            return Py_BuildValue("ssssss", oseq[0], oseq[1], oseq[2], 
                                oseq[4], oseq[3], oseq[5]);
        case 2:
            return Py_BuildValue("ssssss", oseq[0], oseq[1], oseq[2], 
                                oseq[5], oseq[4], oseq[3]);
    }
    return NULL;
}

static PyObject *
seq_type(PyObject *self, PyObject *args)
{
    const char * iseq = NULL;
    unsigned int ilen = 0;

    if(!PyArg_ParseTuple(args, "s#:seq_type", &iseq, &ilen))
    {
        return NULL;
    }

    int DNA = 1, RNA = 1, AMINO = 1, TMP;
    unsigned int i, j;
    char c;
    for(i = 0; i < ilen; i++)
    {
        c = iseq[i];
        
        if(c < 'a')
        {
            c = c - 'A' + 'a';
        }

        if(DNA)
        {
            TMP = 0;
            for(j = 0; j < DNA_LEN; j++)
            {
                if(GENERIC_DNA[j] == c)
                    TMP = 1;
            }
            if(TMP == 0)
            {
                DNA = 0;
            }
        }
        if(RNA)
        {
            TMP = 0;
            for(j = 0; j < RNA_LEN; j++)
            {
                if(GENERIC_RNA[j] == c)
                    TMP = 1;
            }
            if(TMP == 0)
            {
                RNA = 0;
            }
        }
        if(AMINO)
        {
            TMP = 0;
            for(j = 0; j < AMINO_LEN; j++)
            {
                if(GENERIC_AMINO[j] == c)
                    TMP = 1;
            }
            if(TMP == 0)
            {
                AMINO = 0;
            }
        }
    }

    if(DNA)
        return Py_BuildValue("s", "DNA");
    if(RNA)
        return Py_BuildValue("s", "RNA");
    if(AMINO)
        return Py_BuildValue("s", "AMINO");

    PyErr_SetString(PyExc_ValueError, "Unknown Alphabet");
    return NULL;
}


// Module functions table.

static PyMethodDef
module_functions[] = {
    { "seq_type", seq_type, METH_VARARGS,
        "Attempt to guess the type of the sequence. Return \"DNA\", \"RNA\" "
        "or \"AMINO\""},
    { "translate", translate, METH_VARARGS, 
        "Convert a DNA sequence to its (+1 frame) Amino Acid sequence" },
    { "sixFrameTranslation", sixFrameTranslation, METH_VARARGS, 
        "Convert a DNA sequence to its six possible Amino Acid sequences" },
    { NULL }
};

// This function is called to initialize the module.

void
initsequtils(void)
{
    Py_InitModule3("sequtils", module_functions, 
        "Various sequence tools");
}