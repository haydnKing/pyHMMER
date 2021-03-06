// sequtils: useful sequence things

#include "Python.h"

#include <stdio.h>
#include <math.h>

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
const char GENERIC_AMINO[] = "acdefghiklmnpqrstvwy*";

const unsigned int DNA_LEN = 15;
const unsigned int RNA_LEN = 15;
const unsigned int AMINO_LEN = 21;


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

char intAsBase(int in){
    switch(in){
        case 0:
            return 't';
        case 1:
            return 'c';
        case 2:
            return 'a';
        case 3:
            return 'g';
    }
    return 'x';
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

char* _translate(const char* iseq, int ilen, int frame){

    int fwd = (frame > 0);
    frame = abs(frame) - 1;

    unsigned int olen = (ilen - frame) / 3;
    unsigned int i, j;

    char * oseq = malloc(olen * sizeof(char)+1);
    oseq[olen] = '\0';

    unsigned int tmp[] = {0,0,0};

    int di = fwd ? 1 : -1;
    int offset = fwd ? 0 : ilen-1;

    for(i = frame; i < ilen-2; i+=3)
    {
        for(j = 0; j < 3; j++){
            if(!baseAsInt(iseq[offset + di*(i+j)], tmp+j))
            {
                unknown_base(iseq[offset + di*(i+j)], i);
                free(oseq);
                return NULL;
            }
            if(fwd==0)
                tmp[j] = (tmp[j]+2)%4;
        }
        translate_codon(tmp[0],tmp[1],tmp[2], oseq + i/3);
    }

    return oseq;
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

    char* oseq = _translate(iseq, ilen, 1);

    if(oseq == NULL){
        return NULL;
    }

    return Py_BuildValue("s", oseq);
}

/*
 * sixFrameTranslation: iterator to return the translations one by one
 *
 */

typedef struct {
    PyObject_HEAD
    int pos;
    int ilen;
    const char* iseq;
} sequtils_SFT;

int SFT_frames[] = {1,2,3,-1,-2,-3};

/*
 * __iter__: returns self
 */
PyObject* sequtils_SFT_iter(PyObject *self)
{
    Py_INCREF(self);
    return self;
}

PyObject* sequtils_SFT_iternext(PyObject *self)
{
    sequtils_SFT *p = (sequtils_SFT *)self;
    if (p->pos < 6) {
        char* oseq = _translate(p->iseq, p->ilen, SFT_frames[p->pos]);
        (p->pos)++;
        if(oseq == NULL)
            return NULL;
        return Py_BuildValue("s", oseq);
    } else {
        /* Raising of standard StopIteration
         * exception with empty value. */
        PyErr_SetNone(PyExc_StopIteration);
        return NULL;
    }
}

static PyTypeObject sequtils_SFTType = {
    PyObject_HEAD_INIT(NULL)
        0,                         /*ob_size*/
    "sequtils._SFT",            /*tp_name*/
    sizeof(sequtils_SFT),       /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    0,                         /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_ITER,
    /* tp_flags: Py_TPFLAGS_HAVE_ITER tells python 
     * to use tp_iter and tp_iternext fields. */
    "Internal SFT iterator object.",           /* tp_doc */
    0,  /* tp_traverse */
    0,  /* tp_clear */
    0,  /* tp_richcompare */
    0,  /* tp_weaklistoffset */
    sequtils_SFT_iter,  /* tp_iter: __iter__() method */
    sequtils_SFT_iternext  /* tp_iternext: next() method */
};

PyObject * sixFrameTranslation(PyObject *self, PyObject *args)
{
    const char * iseq = NULL;
    unsigned int ilen = 0;
    sequtils_SFT *p;

    if(!PyArg_ParseTuple(args, "s#:sixFrameTranslation", &iseq, &ilen))
    {
        return NULL;
    }

    if(ilen < 3)
    {
        PyErr_SetString(PyExc_ValueError, "sequence too short");
        return NULL;
    }

    /* I don't need python callable __init__() method for this iterator,
     *      so I'll simply allocate it as PyObject and initialize it by hand. */

    p = PyObject_New(sequtils_SFT, &sequtils_SFTType);
    if (!p) return NULL;

    /* I'm not sure if it's strictly necessary. */
    if (!PyObject_Init((PyObject *)p, &sequtils_SFTType)) {
        Py_DECREF(p);
        return NULL;
    }

    p->pos = 0;
    p->iseq = iseq;
    p->ilen = ilen;
    return (PyObject *)p;
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
    //record which char ruled out each alphabet
    char notd = ' ', notr = ' ', nota = ' ';
    unsigned int i, j;
    char c;
    for(i = 0; i < ilen; i++)
    {
        c = iseq[i];
        
        //if c is upper case
        if(c >= 'A' && c <= 'Z')
        {
            //convert to lower case
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
                notd = c;
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
                notr = c;
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
                nota = c;
            }
        }

        //If only one has survived
        if((DNA + RNA + AMINO) == 1)
            break;
    }

    if(DNA)
        return Py_BuildValue("s", "DNA");
    if(RNA)
        return Py_BuildValue("s", "RNA");
    if(AMINO)
        return Py_BuildValue("s", "AMINO");

    char err_string[256];
    sprintf(err_string, "Unknown Alphabet: not DNA (%c), RNA (%c) or Amino(%c)",
            notd, notr, nota);
    PyErr_SetString(PyExc_ValueError, err_string);
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
        "Return an iterator to the six possible Amino Acid sequences" },
    { NULL }
};

// This function is called to initialize the module.

void
initsequtils(void)
{
    PyObject* m;

    sequtils_SFTType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&sequtils_SFTType) < 0)  return;

    m = Py_InitModule3("sequtils", module_functions, 
        "Various sequence tools");

    Py_INCREF(&sequtils_SFTType);
    PyModule_AddObject(m, "_MyIter", (PyObject *)&sequtils_SFTType);
}
