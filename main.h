#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <ctime>
#include <memory>

#ifndef MAIN_H

using namespace std;

typedef struct fragment_Node{
    char* fragment;

    fragment_Node* next_Node;

}fragment_Node;


bool read_Qurey(const string& file_Name);
void initialize_Hash();
void hash_Constructor(string fragment_String);
fragment_Node* new_Node(int hash_Index, string fragment_String);
void hash_Deconstructor();
void free_LL( fragment_Node* current_Ptr);
unsigned int radix_Noation(string fragment_Str);
int character_Value(char letter);

#endif 
