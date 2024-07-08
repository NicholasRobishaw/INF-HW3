// includes
#include "main.h"
using namespace std;

class Queries_HT{
    // hash size will be a configurable parameter

    public:
    int hash_Size = 12;
    int hash_mod = hash_Size;
    int radix_Base = 5; // A, C, G, T, N
    int radix_Character_Num[5] = {0, 1, 2, 3, 4};
    int fragment_Size = 33;

    fragment_Node** hash_Table;

    const long MAX_FRAGMENTS = 1000; 
    const long MAX_SCAFFOLD_COUNT = 608;

    unsigned int collision_Count;


    // file reader for query file
    // read query dataset file function
    bool read_Qurey(const string& file_Name){
        ifstream file(file_Name);
        string current_Line;
        int query_Num = 0;
        time_t stop_Watch = 0;
        
        // check if the file can be opened
        if (!file.is_open()) {
            // display error and exit
            cerr << "Error: Unable to open file " << file_Name << "\n";
            return false;
        }
        
        // iterate through file contents line by line
        while ( query_Num < MAX_FRAGMENTS && getline(file, current_Line)) {
            // check if we are at a valid fragment
            if(current_Line[0] != '>' && !current_Line.empty()){
                // increment query count
                query_Num++;

                // attempt to add fragment to hash table
                // if false is returned then there was a problem and program will clean up memory and end
                hash_Constructor(current_Line);
                
            }
        }

        // close the file
        file.close();

        // return sucess
        return true;
    }


    // file reader for genome file
    

    // initialize hash table
    void initialize_Hash(){
        // create an array of hash table size
        // each index of the array will store a LL
        hash_Table[hash_Size];

        collision_Count = 0;

        // iterate through the hash table and set each index to NULL
        for(int index=0; index < hash_Size; index++){
            hash_Table[index] = nullptr;
        }

    }


    // chaining method for duplicate fragments/hash indexes
    // constructor ( LL? )
        // this will add a node to the LL at the hash index
    void hash_Constructor(string fragment_String){
        unsigned int fragment_Radix;
        
        // get the hash index of the fragment
        fragment_Radix = radix_Noation(fragment_String);

        int hash_Index = fragment_Radix % hash_mod; 

        if(hash_Table[hash_Index] != nullptr){
            collision_Count++;
        }

        // crete new node for fragment at hash index ( will account for chaining )
        hash_Table[hash_Index] = new_Node(hash_Index, fragment_String);
    }


    fragment_Node* new_Node(int hash_Index, string fragment_String){
        fragment_Node* new_Node_Ptr = new fragment_Node;

        // create fragment string
        new_Node_Ptr->fragment = new char[fragment_Size];

        // copy fragment string into new node
        for(int index=0; index < fragment_Size; index++){
            new_Node_Ptr->fragment[index] = fragment_String[index];
        }

        // set next pointers to previous head ptrs
        new_Node_Ptr->next_Node = hash_Table[hash_Index];

        return new_Node_Ptr;
    }


    // deconstructor
    void hash_Deconstructor(){
        // iterate through the hash table
        for(int index = 0; index < hash_Size; index++){
            // destroy the LL
            free_LL(hash_Table[index]);
        }

        delete[] hash_Table;
    }

    
    // LL destroyer
    void free_LL( fragment_Node* current_Ptr){
        // check if we are not at the end of the LL
        if( current_Ptr != nullptr ){
            // check if we are at the end of the LL
            if( current_Ptr->next_Node != nullptr){
                // recurse to next node
                free_LL(current_Ptr->next_Node);
            }

            // check if there is data in the fragment string
            if( current_Ptr->fragment != nullptr){
                // free the character fragment array
                delete[] current_Ptr->fragment;
            }

            // free the current node
            delete current_Ptr;
        }
    }


    // search hash table function ( return bool )

    // insert into hash function( at start of LL )


    // We will be converting the 32 character fragments to radix notation so we can store it in the has table
        // the values for the proteins will be
        // A -> 0, C -> 1, G -> 2, T -> 3, N -> 4 

    // convert sequence to radix notation function
    unsigned int radix_Noation(string fragment_Str){

        int string_Size = fragment_Size - 2;
        int index;
        unsigned int radix_Conversion;
        int character_Num;

        // loop through the string from right to left
        for(index = string_Size; index >= 0; index--){
            // get the character number 
            character_Value(fragment_Str[index]);
            
            // convert string to radix notation
            radix_Conversion += (character_Num * radix_Base^(string_Size - index));
        }

        // return radix converted string
        return radix_Conversion;
    }
    
    // returns the radix number for the letter
    int character_Value(char letter){
        int index;

        // loop through the alphabet of values (A, C, G, T, and N)
        for(index = 0; index < radix_Base; index++){
            // check which letter is a match and return the index which is also the number for that letter
            if(letter == radix_Character_Num[index]){
                return index;
            }
        }

        return index;
    }

};

int main(int argc, char* argv[]){
    Queries_HT my_Hash;
    time_t stop_Watch = 0;

    // check if the program has the appropriate number of parameter
    // create error handler for no input file argument
    if( argc < 4 ){
        cout << "Error please input the correct command in\n Program End";
        return 1;
    }

    // Subproblem A -  Assess the impact of the hash table size
    // set hash table to fixed size ( 1 million, 10 million, 30 million and 60 million )

    time(&stop_Watch);
    cout << "Program Start at: " << ctime(&stop_Watch) << endl;

    // initialize the hash table
    my_Hash.initialize_Hash();

    time(&stop_Watch);
    cout << "Starting hash table population at: " << ctime(&stop_Watch) << endl;

    // populate hash table with the sequence fragments from the query dataset
    if(my_Hash.read_Qurey(argv[2])){
        // for each of your 4 hash table sizes, how many collisions did you observe while populating the hash? ( KEEP IN MIND I WILL NEED TO RUN THIS 4 DIFFERENT TIMES WITH THE 4 SIZES)

        //For each of your 4 hash table sizes, how long did it take you to populate the hash table? Do the timing results make sense?
        time(&stop_Watch);
        cout << "End of hash table population at: " << ctime(&stop_Watch) << endl;
    }

    
    // Subproblem B -  Searching Speed
    // set hash table to size 60 million
    // populate the hash table with the sequence fragments from the query dataset.
    // Read in genome sequence

    // implement a search function which would search for 16-character fragments of the subject sequence within the Queries_HT object.
    // Iterate through all 16-character long fragments of the subject dataset, searching for each one in the query dataset.

        // How long did it take to search for every possible 16 character long fragment of the subject dataset within the query dataset?
        // How many such fragments did you find
        // Print the first 10 fragments of the subject dataset that you found within the Query_HT
        
    // clean up hash table
    my_Hash.hash_Deconstructor();

    cout << "Program End\n";
    return 0;
}