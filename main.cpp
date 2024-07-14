// includes
#include "main.h"
using namespace std;

class Queries_HT{

    public:
    int hash_Size;
    int hash_Mod;
    int radix_Base = 5; // A, C, G, T, N
    int radix_Character_Num[5] = {0, 1, 2, 3, 4};
    int collision_Count;
    int fragment_Size = 33;
    
    const long MAX_FRAGMENTS = 125000000; 
    const long MAX_SCAFFOLD_COUNT = 608;
    long genome_Size = 0;
    long allocated_Genome_Size = 3000000000;
    long fragments_Found;
    long allocated_Frag_Size = 3057186554; 
    long scaffold_Count = 0;

    char* genome_Data = nullptr;
    char** fragments;
    fragment_Node** hash_Table;
    

    // file reader for query file
    // read query dataset file function
    bool read_Qurey(const string& file_Name){
        ifstream file(file_Name);
        string current_Line;
        int query_Num = 0;
        
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
    bool file_reader(const string& file_Name) {
        // variables
        ifstream file(file_Name);
        string current_Scaffold_Name;
        string current_Line;
        string temp_Genome = "";
        long prev_genome_size;
        unsigned long index;

        // check if the file can be opened
        if (!file.is_open()) {
            // display error and exit
            cerr << "Error: Unable to open file " << file_Name << "\n";
            return false;
        }
    
        // initialize starting array, set to 3 billion to save on runtime complexity
        genome_Data = new char[allocated_Genome_Size];
        
        // iterate through file contents line by line
        while ( scaffold_Count < MAX_SCAFFOLD_COUNT && getline(file, current_Line)) {
            // check if we hit a header
            if(current_Line[0] == '>' && !temp_Genome.empty()){
                
                // increment scaffold count by 1
                scaffold_Count++;
                
                // calculate new size for genome array
                prev_genome_size = genome_Size;
            
                // resize the genome array to fit last scaffold and copy data into new array
                genome_Constructor(temp_Genome, prev_genome_size);

                // reset the temp genome string to empty
                temp_Genome = "";
            }

            else if (!current_Line.empty() && current_Line[0] != '>') {
                // Iterate through current line and append characters to temp_Genome
                for (index = 0; index < current_Line.length() && index < 80; index++) {
                    if(                  current_Line[index] == 'A' 
                                      || current_Line[index] == 'C' 
                                      || current_Line[index] == 'G' 
                                      || current_Line[index] == 'T' 
                                      || current_Line[index] == 'N'){
                    
                    
                        // Append character to temp_Genome
                        temp_Genome += current_Line[index];
                    }
                }
            }
        }

        // post append information of last scaffold
        if (!temp_Genome.empty()) {
            // calculate new genome array size
            prev_genome_size = genome_Size;
            
            // resize the genome array to fit last scaffold and copy data into new array
            genome_Constructor(temp_Genome, prev_genome_size);
        }

        // close the file
        file.close();

        // return sucess
        return true;
        
    }// end of file reader function


    // deconstructor for deallocating the genome character array
    void genome_Deconstructor(char* genome_arr){
        delete[] genome_arr;
        genome_arr=nullptr;
    }


    // function to resize the genome character array
    void genome_Constructor(const string to_Add, long cat_Index) {
        long index;
        long length = to_Add.length();
        
        // check if the array needs to be resized
        if(genome_Size + length >= allocated_Genome_Size){
            // calculate the new size of the genome
            long new_Size = genome_Size + length + 1;
            
            // allocate new array size
            char* new_Genome_Arr = new char[new_Size];
            
            // copy data from old array to new one
            for(index = 0; index < genome_Size; index++){
                new_Genome_Arr[index] = genome_Data[index];
            }
            
            // deallocate memory from old array
            delete[] genome_Data;
            
            // set genome pointer to new array
            genome_Data = new_Genome_Arr;
            
            // set the new allocated size for genome array
            allocated_Genome_Size = new_Size;
        }
        
            
        for(index = 0; index < length; index++){
            //add in the new data 
            genome_Data[cat_Index + index] = to_Add[index];
        }
        
        // set the new genome size
        genome_Size = cat_Index + to_Add.length();
    }


    // initialize hash table
    void initialize_Hash(int size){
        // create an array of hash table size
        // each index of the array will store a LL
        hash_Size = size;
        hash_Mod = size;
        
        hash_Table = new fragment_Node*[hash_Size];

        collision_Count = 0;

        fragments_Found = 0;

        // iterate through the hash table and set each index to NULL
        for(int index=0; index < hash_Size; index++){
            hash_Table[index] = nullptr;
        }

    }


    // function to populate the hash table ( Avoids adding duplicate fragments to help runtime)
    void hash_Constructor(string fragment_String) {
        unsigned int fragment_Radix = radix_Noation(fragment_String);
        int hash_Index = fragment_Radix % hash_Mod;
        
        // Traverse the linked list at hash_Index to check for existing fragment
        fragment_Node* current = hash_Table[hash_Index];
        
        if(current != nullptr){
            // Found a collision
            collision_Count++;
        }
        
        while (current != nullptr) {
            if (current->fragment == fragment_String) {
                
                // end of function since a duplicate was found
                return; 
            }
            current = current->next_Node;
        }
        
        // Create new node and insert at hash_Index
        hash_Table[hash_Index] = new_Node(hash_Index, fragment_String);
    }


    // function to add a new node that Link to the hash table
    fragment_Node* new_Node(int hash_Index, string fragment_String){
        // create a new fragment node
        fragment_Node* new_Node_Ptr = new fragment_Node;

        // create fragment string
        new_Node_Ptr->fragment = new char[fragment_Size];

        // copy fragment string into new node
        for(int index=0; index < fragment_Size; index++){
            new_Node_Ptr->fragment[index] = fragment_String[index];
        }

        // set next pointers to previous head ptrs
        new_Node_Ptr->next_Node = hash_Table[hash_Index];

        // return ptr to new node
        return new_Node_Ptr;
    }


    // deconstructor
    void hash_Deconstructor(){
        // iterate through the hash table
        for(int index = 0; index < hash_Size; index++){
            // destroy the LL
            free_LL(hash_Table[index]);
        }

        // destroy the main array
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


    // function to convert the fragment string to a unsigned radix number
    unsigned int radix_Noation(string fragment_Str) {
        unsigned int radix_Conversion = 0;
        int character_Num;
        
        // iterate thorugh the string from 
        for (int index = 0; index < fragment_Size; index++) {
            // get the character value
            character_Num = character_Value(fragment_Str[index]);
            
            //calculate the radix conversion for the character
            radix_Conversion += (character_Num * pow(radix_Base, fragment_Size - 1 - index));
        }
    
        return radix_Conversion;
    }
    
    // returns the radix number for the letter
    int character_Value(char letter){
        int index;
        
        switch (letter) {
            case 'A':
                index = 0;
                break;
            case 'C':
                index = 1;
                break;
            case 'G':
                index = 2;
                break;
            case 'T':
                index = 3;
                break;
            default:
                index = 4;
        }
        
        return index;
    }


    // iterate through the genome array and attempt to locate every 32 character fragment in the hash table
    void search_Function(){
        long index = 0, high_Index = 0, genome_Index;
        long inner_Index;
        char test_Frag[fragment_Size-1];

        // create a n array to store the found fragments from the genome
        fragments = new char*[allocated_Frag_Size];

        // loop through the genome array
        while( index < genome_Size && high_Index < genome_Size ){

            // create fragment ( index to index+31 )
            high_Index = index + (fragment_Size - 1);

            if(high_Index < genome_Size){
                
                
                // if( index % 500000 == 0 ){
                //     cerr << "On index: " << index << endl;
                // }

                // reset the genome index back to the start
                genome_Index = 0;
    
                // read into character test frag array
                for( inner_Index = index; inner_Index < high_Index; inner_Index++){
                    test_Frag[genome_Index] = genome_Data[inner_Index];
                    genome_Index++;
                    
                }
                
                //set end of fragment
                test_Frag[fragment_Size-1] = '\0';
    
                // pass to helper function to look through the hash table ( return type unknown right now )
                if( lookup_Hash(test_Frag) ){
                    // increment found frag counter
                    fragments_Found++;
    
                    // if fragment was found add to fragment found list ( Possible 2d array of character fragments)
                    found_Frags(test_Frag);
                }
            
            }
        
            index++;
        }   

    }


    // helper function for search_Function to look for fragment in hash table
    bool lookup_Hash( char test_Fragment[] ){
        fragment_Node* temp_Ptr = nullptr;
        unsigned int fragment_Radix;

        // convert string to radix notation
        fragment_Radix = radix_Noation(test_Fragment);

        // get the hash index
        int test_Hash = fragment_Radix % hash_Mod;

        // set the temp pointer to the hash index
        temp_Ptr = hash_Table[test_Hash];

        // loop thorugh the Linked list until hitting null or match
        while( temp_Ptr != nullptr ){
            
            // if the current nodes fragment matches the test fragment
            if( strcmp(test_Fragment, temp_Ptr->fragment) == 0){
                // return success
                return true;
            }

            // otherwise move to next node
            temp_Ptr = temp_Ptr->next_Node;
        }
        
        // otherwise return failure
        return false;
    }


    // function for adding new found fragment to storage
    void found_Frags( char fragment[]){
        long index;
        int frag_Index;
        char **new_Arr = nullptr;
        char *new_Frag = nullptr;
        
        if( fragments_Found >= allocated_Frag_Size ){
            
            allocated_Frag_Size += 10000;
            
            // make a newly resized array to house the extra found fragment
            new_Arr = new char*[allocated_Frag_Size];
        
            // copy existing pointers from old array to the new resized array
            for( index = 0; index < fragments_Found; index++){
                new_Arr[index] = fragments[index];
            }
    
            // delete old fragment array
            delete[] fragments;
            
            // set head fragments ptr to the newly created array
            fragments = new_Arr;
            
        }

        // create a character array for the new sport
        new_Frag = new char[fragment_Size];

        // copy new fragment to the new spot
        for( frag_Index = 0; frag_Index < fragment_Size-1; frag_Index++){
            new_Frag[frag_Index] = fragment[frag_Index];
        }

        // add in the new fragment at the end
        fragments[fragments_Found-1] = new_Frag;
    }


    // function to deallocate the fragment array
    void found_Frags_Deconstructor(){
        long index;

        // loop through the array of pointers
        for( index = 0; index < fragments_Found; index++){
            // free each index
            delete[] fragments[index];
        }

        // free main array
        delete[] fragments;
    }
};

int main(int argc, char* argv[]){
    Queries_HT my_Hash;
    time_t stop_Watch = 0;
    int index, inner_Index, h_Size=0;
    bool part_B = false;

    // ./my_program human.txt human_reads_125_32.fa <HASH SIZE HERE> <-Subpart>

    // check if the program has the appropriate number of parameter
    // create error handler for no input file argument
    if( argc < 4 ){
        cout << "Error please input the correct command in\n Program End";
        return 1;
    }

    // check and set hash size
    if     (strcmp("60000000", argv[3]) ==0 ){
        // set the hash size
        h_Size = 60000000;
    } 
    else if (strcmp("30000000", argv[3]) ==0) {
        // set the hash size
        h_Size = 30000000;
    } 
    else if (strcmp("10000000", argv[3]) ==0) {
        // set the hash size
        h_Size = 10000000;
    } 
    else if (strcmp("1000000" , argv[3]) ==0 ){
        // set the hash size
        h_Size = 1000000;
    }
    // otherwise assunme bad input and return error
    else{
        cout << "Error please input the correct hash size in\n Program End";
        return 1;   
    }
    
    // check if a subproblem was incorrectly inputted
    if(strcmp("-A", argv[4]) !=0 && strcmp("-B", argv[4]) !=0){
        cout << "Error please input the correct subpart in\n Program End";
        return 1; 
    }
    
    // check if sub problem B will be executed
    if( strcmp("-B", argv[4]) ==0 ){
        part_B = true;
    }
    
    // Subproblem A -  Assess the impact of the hash table size
    // set hash table to fixed size ( 1 million, 10 million, 30 million and 60 million )
    time(&stop_Watch);
    cout << "Program Start at: " << ctime(&stop_Watch) << endl;

    // initialize the hash table
    my_Hash.initialize_Hash(h_Size);

    time(&stop_Watch);
    cout << "Starting hash table population at: " << ctime(&stop_Watch) << endl;

    // populate hash table with the sequence fragments from the query dataset
    if(my_Hash.read_Qurey(argv[2])){
        // for each of your 4 hash table sizes, how many collisions did you observe while populating the hash? ( KEEP IN MIND I WILL NEED TO RUN THIS 4 DIFFERENT TIMES WITH THE 4 SIZES)

        //For each of your 4 hash table sizes, how long did it take you to populate the hash table? Do the timing results make sense?
        time(&stop_Watch);
        cout << "End of hash table population at: " << ctime(&stop_Watch) << endl;
    }

    // display info
    cout << "Hash Size           : " << my_Hash.hash_Size << endl;
    cout << "Number of collisions: " << my_Hash.collision_Count << endl;

    // Subproblem B -  Searching Speed
    // set hash table to size 60 million
    // populate the hash table with the sequence fragments from the query dataset.
    // Read in genome sequence
    if( part_B ){
        if(my_Hash.file_reader(argv[1])){
            
            cout << "\nGenome size: " << my_Hash.genome_Size << endl;
            
            time(&stop_Watch);
            cout << "Starting Search at: " << ctime(&stop_Watch) << endl;
    
            // call to main driver search function
            my_Hash.search_Function();
                
            // How long did it take to search for every possible 32 character long fragment of the subject dataset within the query dataset?
            time(&stop_Watch);
            cout << "Finished Search at: " << ctime(&stop_Watch) << endl;   
    
            // How many such fragments did you find
            cout << "There were " << my_Hash.fragments_Found << " fragments found\n";
            
            // Print the first 10 fragments of the subject dataset that you found within the Query_HT
            for( index = 0; index < 10; index++){
                cout << " ";
    
                for( inner_Index = 0; inner_Index < my_Hash.fragment_Size; inner_Index++){
                    cout << my_Hash.fragments[index][inner_Index];
                }
    
                cout << endl;
            }
    
            cout << endl;
        }
        
        // clean up part B memory
        my_Hash.genome_Deconstructor(my_Hash.genome_Data);
        my_Hash.found_Frags_Deconstructor();
    }
    
    // clean up hash table
    my_Hash.hash_Deconstructor();

    cout << "Program End\n";
    return 0;
}
