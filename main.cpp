/* 
 * File:   main.cpp
 * Author: Sergey
 *
 * Created on February 7, 2013, 12:37 PM
 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <string.h>
#include "utils_unsorted.h"
#include "IntervalTreeMy.h"

#define DEFAULT_SEP '\t'
#define ID_SEP "|"
#define ANNO_SEP ","

//using namespace std;

/*
 * Program arguments example:
 * 
 * -a file_with_annotation -ai 5 -ai 9 -al 6 -ah 7 -aa 11 \
 * -i file_to_annotate -ii 0 -ii 1 -il 2 -ih 2 \
 * -o output_file \
 * [-asep sep_char] [-isep sep_char]
 * 
 * All columns are zero-bases (first column is 0)!
 * Input file should be tab-separated!
 * 
 * -a    : file with annotation
 * -ai   : id columns, like chromosome and strand, is case of multiple "ai", id 
 *         is a concatenation of them separated by '|'
 * -al   : column for annotation interval's low
 * -ah   : column for annotation interval's high
 * -aa   : column with annotation
 * -i    : file that needs to be annotated
 * -ii   : id columns, to intersect with -ai..., concatenated same as -ai
 * -il   : column for input file interval's low 
 * -ih   : column for input file interval's high
 * -o    : copy of -i file with annotation appended to the end of lines for 
 *         intervals that overlap between -i and -ai files
 * -asep : separator between columns in -a file
 * -isep : separator between columns in -i file
 */
// FIXME: no checks for many things
// ADDON: add support using SQL bases for annotation

int main(int argc, char** argv)
{    
    //cout << "Hello World!" << endl;
    char* fname_anno;
    char* fname_input;
    char* fname_output;

    vector<unsigned int> col_ai;
    int col_al = -1, col_ah = -1, col_aa = -1;
    vector<unsigned int> col_ii;
    int col_il = -1, col_ih = -1;
    
    char aSEP = DEFAULT_SEP;
    char iSEP = DEFAULT_SEP;

    // TODO: use config library and/or classes to read config/arguments
    for ( int i = 0; i < argc; i++ )
    {
        if ( strcmp(argv[i], "-a") == 0 )
        {
            fname_anno = argv[i + 1];
        }
        else if ( strcmp(argv[i], "-ai") == 0 )
        {
            col_ai.push_back(atoi(argv[i + 1]));
        }
        else if ( strcmp(argv[i], "-al") == 0 )
        {
            col_al = atoi(argv[i + 1]);
        }
        else if ( strcmp(argv[i], "-ah") == 0 )
        {
            col_ah = atoi(argv[i + 1]);
        }
        else if ( strcmp(argv[i], "-aa") == 0 )
        {
            col_aa = atoi(argv[i + 1]);
        }
        else if ( strcmp(argv[i], "-i") == 0 )
        {
            fname_input = argv[i + 1];
        }
        else if ( strcmp(argv[i], "-ii") == 0 )
        {
            col_ii.push_back(atoi(argv[i + 1]));
        }
        else if ( strcmp(argv[i], "-il") == 0 )
        {
            col_il = atoi(argv[i + 1]);
        }
        else if ( strcmp(argv[i], "-ih") == 0 )
        {
            col_ih = atoi(argv[i + 1]);
        }
        else if ( strcmp(argv[i], "-o") == 0 )
        {
            fname_output = argv[i + 1];
        }
        else if ( strcmp(argv[i], "-asep") == 0 )
        {
            escape_hex(argv[i + 1], aSEP);
        }
        else if ( strcmp(argv[i], "-isep") == 0 )
        {
            escape_hex(argv[i + 1], iSEP);
        }
    }
    
    if ( col_ai.size() == 0 || col_ii.size() == 0 ||
            col_al < 0 || col_ah < 0 || col_aa < 0 ||
            col_il < 0 || col_ih < 0 )
    {
        cout << "Invalid arguments." << endl;
        exit(-1);
    }

    ifstream fsanno(fname_anno);
    ifstream fsin(fname_input);
    ofstream fsout(fname_output, ios::trunc);

    if ( !fsanno.is_open() )
    {
        cout << "Unable to open file \"" << fname_anno << "\"" << endl;
        exit(-2);
    }
    if ( !fsin.is_open() )
    {
        cout << "Unable to open file \"" << fname_input << "\"" << endl;
        exit(-2);
    }
    if ( !fsout.is_open() )
    {
        cout << "Unable to open file \"" << fname_output << "\"" << endl;
        exit(-2);
    }

    // collection of interval trees for each ID in annotation/input files
    // e.g. a tree for each chrom and strand (chrX|+ etc.) 
    map<string, IntervalTree> trees_map;
    map<string, IntervalTree>::iterator trees_iter;

    // stuff for reading files
    string line;
    vector<string> data;

    // TODO: make something like an annotation reader classes
    cout << "Reading annotation intervals...\n" << flush;
    while (fsanno.good())
    {
        getline(fsanno, line);
        if ( line.length() == 0 ) continue;
        split(data, line, aSEP);

        string low = data[col_al];
        string high = data[col_ah];
        //if ( !is_number(low) || !is_number(high) ) continue;

        // make ID
        std::vector<unsigned int>::iterator it = col_ai.begin();
        string id = data[*it];
        ++it;
        for (; it != col_ai.end(); ++it )
        {
            id += ID_SEP + string(data[*it]);
        }

        // try to find a tree for the ID
        trees_iter = trees_map.find(id);
        // new interval to be added to the tree with the ID
        Interval<string> *p_ival = new Interval<string > (atoi(low.c_str()), atoi(high.c_str()), data[col_aa]);
        // Add new interval to the existing tree or to the new tree for the ID
        if ( trees_iter != trees_map.end() )
        {
            trees_iter->second.Insert(p_ival);
        }
        else
        {
            IntervalTree tree;
            tree.Insert(p_ival);
            trees_map[id] = tree;
        }
    }
    fsanno.close();
    cout << "Annotation intervals read\n" << flush;

    // TODO: make something input reader helper classes
    cout << "Reading input file, writing output...\n" << flush;
    while (fsin.good())
    {
        string stres = "";
        getline(fsin, line);
        if ( line.length() == 0 ) continue;
        split(data, line, iSEP);

        string low = data[col_il];
        string high = data[col_ih];
        //if ( !is_number(low) || !is_number(high) ) continue;

        // make ID
        std::vector<unsigned int>::iterator it = col_ii.begin();
        string id = data[*it];
        ++it;
        for (; it != col_ii.end(); ++it )
        {
            id += ID_SEP + string(data[*it]);
        }
        
        trees_iter = trees_map.find(id);
        if ( trees_iter != trees_map.end() )
        {
            vector<void *> * results = trees_iter->second.findOverlapping(atoi(low.c_str()), atoi(high.c_str()));
            if ( results->size() > 0 )
            {
                unsigned int i = 0;

                for ( i = 0; i < results->size() - 1; i++ )
                {
                    Interval<string>* currentInterval = (Interval<string> *) (*results)[i];
                    stres = stres + currentInterval->value + ANNO_SEP;
                }
                Interval<string>* currentInterval = (Interval<string> *) (*results)[i];
                stres = stres + currentInterval->value;
            }
        }
        fsout << line << "\t" << stres << "\n";
    }
    fsin.close();
    fsout.close();
    cout << "Annotation done\n" << flush;

    return 0;
}
