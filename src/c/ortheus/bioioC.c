#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <float.h>

#include "fastCMaths.h"
#include "commonC.h"
#include "ctype.h"
#include "bioioC.h"

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
//integer reader / writer
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////

void readIntegers(FILE *file, int intNumber, INT_32 *iA) { 
    INT_32 i;
    
    for(i=0; i<intNumber; i++) {
        fscanf(file, INT_STRING, iA + i);
    }
}

void writeIntegers(FILE *file, int intNumber, INT_32 *iA) {
    INT_32 i;
    
    for(i=0; i<intNumber; i++) {
        fprintf(file, INT_STRING "\n", iA[i]);
    }
}

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
//fasta reader
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////

char *addSeqToList(char *seq, INT_32 *length, INT_32 *maxLength, struct List *seqs, struct List *seqLengths) {
    INT_32 i;
    seq = arrayPrepareAppend(seq, maxLength, *length+1, sizeof(char));
    seq[*length] = '\0';
    i = sizeof(char)*((*length) + 1);
    listAppend(seqs, memcpy(mallocLocal(i), seq, i));
    listAppend(seqLengths, constructInt(*length));
    (*length)  = 0;
    return seq;
}  

void fastaRead(FILE *fastaFile, struct List *seqs, struct List *seqLengths) {
    //reads in group of sequences INT_32o lists
    char j;
    INT_32 k;
    static char *seq;
    static INT_32 seqLength;
    
    k=0;
    while((j = getc(fastaFile)) != EOF) { //initial termininating characters
        if(j == '>') { //fasta start
            fastaStart:
            while(TRUE) {
                j = getc(fastaFile); 
                if(j == EOF) {
                    seq = addSeqToList(seq, &k, &seqLength, seqs, seqLengths); //lax qualification for a sequence
                    return;
                }
                if(j == '\n')  {
                    break;
                }
            }
            while(TRUE) { //start of sequence
                j = getc(fastaFile);
                if(j == EOF) {
                    seq = addSeqToList(seq, &k, &seqLength, seqs, seqLengths);
                    return;
                }
                if(j != '\n') {
                    if(j == '>') {
                        //end of seq
                        seq = addSeqToList(seq, &k, &seqLength, seqs, seqLengths);
                        goto fastaStart;
                    }
                    else { //valid char
                        seq = arrayPrepareAppend(seq, &seqLength, k+1, sizeof(char));
                        seq[k++] = j;
                    }
                }
            }    
        }
    }
}

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
//newick tree parser
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////


inline char *eatWhiteSpace(char *string) {
    while(*string != '\0' && isspace(*string)) {
        string++;
    }
    return string;
}

inline char* eatString(char *string, char **newString) {
    char *i;
  
    i = string;
    while(*string != '\0' && !isspace(*string)) {
           string++; 
    }
    (*newString) = memcpy(mallocLocal((string-i + 1)*sizeof(char)), i, (string-i)*sizeof(char));
    (*newString)[string - i] = '\0'; 
    return eatWhiteSpace(string);
}

char *replaceString(char *oldString, char old, char *new, INT_32 newLength) {
    char *i;
    INT_32 j;
    char *k;
    char *newString;
    j=0;
    k=oldString;
    for(i=oldString; *i != '\0'; i++) {
        if(*i == old) {
            j++;
        }
    }
    newString = mallocLocal(sizeof(char)*(j*(newLength-1) + (i-k) + 1));
    k=newString;
    for(i=oldString; *i != '\0'; i++) {
        if(*i == old) {
            for(j=0; j<newLength; j++) {
                *k++ = new[j];
            }
        }
        else {
            *k++ = *i;
        }
    } 
    *k = '\0';
    return newString;
}

char *replaceAndFreeString(char *oldString, char old, char *new, INT_32 newLength) {
    char *i;
    
    i = replaceString(oldString, old, new, newLength);
    free(oldString);
    return i;
}

inline char *parseFloat(char *string, FLOAT_32 *f) {
    float f2;
    
    f2 = FLT_MIN;
    sscanf(string, "%f", &f2);
    assert(f2 != FLT_MIN);
    *f = f2;
    while(*string != '\0' && !isspace(*string)) {
        //assert(isdigit(*string) || *string == '.'); --thats wrong
        string++;
    }
    return eatWhiteSpace(string);
}

char *newickTreeParser_fn(char *newickTreeString, FLOAT_32 *distance) {
    if(*newickTreeString != '\0') {
        if (*newickTreeString == ':') {
            newickTreeString++;
            return parseFloat(eatWhiteSpace(newickTreeString), distance);
        }
    }
    return newickTreeString;
}
    
char *newickTreeParser_fn2(char *newickTreeString, FLOAT_32 defaultDistance, struct BinaryTree **binaryTree, struct List *strings) {
    struct BinaryTree *temp1;
    struct BinaryTree *temp2;
    FLOAT_32 f;
    
    temp1 = NULL;
    temp2 = NULL;
    
    if(*newickTreeString == '(') {
        temp1 = NULL;
        newickTreeString = eatWhiteSpace(++newickTreeString);
        assert(*newickTreeString != ')');
        do {
            newickTreeString = newickTreeParser_fn2(newickTreeString, defaultDistance, &temp2, strings);
            if(temp1 != NULL) {
                //merge node
                temp1 = constructBinaryTree(0.0f, TRUE, temp1, temp2); //default to zero distance for nodes of 
            }
            else {
                temp1 = temp2;
            }
        } while (*newickTreeString != ')');
        newickTreeString = eatWhiteSpace(++newickTreeString);
    }
    else {
        listAppend(strings, NULL);
        newickTreeString = eatString(newickTreeString, (char **)&(strings->list[strings->length-1]));
        temp1 = constructBinaryTree(0.0f, FALSE, NULL, NULL); 
    }
    f = defaultDistance;
    newickTreeString = newickTreeParser_fn(newickTreeString, &f);
    temp1->distance += f;
    (*binaryTree) = temp1;
    return newickTreeString;
}
 
struct BinaryTree *newickTreeParser(char *newickTreeString, FLOAT_32 defaultDistance, struct List **strings) {
    struct BinaryTree *binaryTree;
    char *i;
    *strings = constructEmptyList(0, free);
    //lax newick tree parser
    newickTreeString = replaceString(newickTreeString, '(', " ( ", 3);
    newickTreeString = replaceAndFreeString(newickTreeString, ')', " ) ", 3);
    newickTreeString = replaceAndFreeString(newickTreeString, ':', " : ", 3);
    newickTreeString = replaceAndFreeString(newickTreeString, ',', " ", 1);
    newickTreeString = replaceAndFreeString(newickTreeString, ';', "", 0);
    i = newickTreeString;
    newickTreeString = eatWhiteSpace(newickTreeString);
    newickTreeParser_fn2(newickTreeString, defaultDistance, &binaryTree, *strings);
    free(i);

    return binaryTree;
}

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
//read in fasta file, and turn into a column alignment
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////

struct CharColumnAlignment *multiFastaRead(char *fastaFile) {
    struct List *seqs;
    struct List *seqLengths;
    FILE *fileHandle;
    INT_32 alignmentLength = 0;
    struct CharColumnAlignment *charColumnAlignment;
    INT_32 i;
    INT_32 j;
    INT_32 k;
    
    seqs = constructEmptyList(0, free);
    seqLengths = constructEmptyList(0, free);
    fileHandle = fopen(fastaFile, "r");
    fastaRead(fileHandle, seqs, seqLengths);
    fclose(fileHandle); 
   
    alignmentLength = 0;
    if(seqLengths->length != 0) {
        alignmentLength = listGetInt(seqLengths, 0);
    }
    for(i=0; i<seqLengths->length; i++) {
        assert(alignmentLength == listGetInt(seqLengths, 0)); 
    }
    charColumnAlignment = mallocLocal(sizeof(struct CharColumnAlignment));
    charColumnAlignment->columnNo = alignmentLength;
    charColumnAlignment->seqNo = seqLengths->length;
    charColumnAlignment->columnAlignment = mallocLocal(sizeof(char)*(charColumnAlignment->columnNo)*(charColumnAlignment->seqNo));
    k=0;
    for(i=0; i<alignmentLength; i++) { 
        for(j=0; j<seqLengths->length; j++) {
            charColumnAlignment->columnAlignment[k++] = ((char *)seqs->list[j])[i];
        }
    }
    destructList(seqs);
    destructList(seqLengths);
    return charColumnAlignment;
}

char *charColumnAlignment_getColumn(struct CharColumnAlignment *charColumnAlignment, INT_32 col) {
    return &charColumnAlignment->columnAlignment[col*charColumnAlignment->seqNo];
}

void destructCharColumnAlignment(struct CharColumnAlignment *charColumnAlignment) {
    free(charColumnAlignment->columnAlignment);
    free(charColumnAlignment);
}
