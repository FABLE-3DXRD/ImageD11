

// From http://www.codetoad.com/forum/24_24964.asp

#include <iostream>

// Needed for the ifstream class
#include <fstream>

// Needed for the setw stream manipulator
// #include <iomanip>

#include <vector>

using namespace std;


int main()
{
// STEP1: OPEN THE FILE TO BE READ

// ifstream is an object reprenting an input file stream
// it can be used just like you use cin but reads a file 
// instead of the standard input stream
// To open a file with it, specify the name of the file in the parenthesis
// alternatively, you could write:
// ifstream file_to_read;
// file_to_read.open("data.txt");
// The file does not have to be .txt, any file name will do, 
// like "hello", as long as its ASCII
// (Binary file can also be accepted, it's just that the program 
// I wrote does not process binary
  ifstream file_to_read("peaks.out_merge_t200");
 const int max_num_of_char_in_a_line = 4096*8,
// Maximum number of characters expected in a single line in the header
num_of_header_lines = 1; // Number of header files to skip

// STEP2: SKIP ALL THE HEADER LINES

for(int j = 0; j < num_of_header_lines; ++j)
// The ignore member function ignores the number of characters in its first parameter
// or until the character in the 2nd paramter is reached. In this case, the newline character.
file_to_read.ignore(max_num_of_char_in_a_line, '\n');

// STEP3: READ THE FILE AND STORE THE DATA READ INTO ARRAYS

// These are arrays to store the processed data.
// vector objects are being used instead of the traditional array
// because vector objects can grow dynamically in size at runtime
// while traditional arrays have a fixed size.
// (ex: Time time_array[50] can hold 50 Time structures at max
// Note: The type of elements stored in a vector is specified in the
// template argument (that is between < and >)
//vector<Time> time_array;
 vector<float> v_array;
// Keep reading while the End Of File character is NOT reached
  float v;
  int i;
while(!file_to_read.eof()) {
  //Time time_read;
  // Read 2 characters and store them in time_read.hours
  // P.S.: the setw stream manupulator tells the stream to read only 2 characters
  // file_to_read >> setw(2) >> time_read.hours;

  // Ignore one character. In our case, the : character
  //file_to_read.ignore(1);
  //file_to_read >> setw(2) >> time_read.minutes;

  //file_to_read.ignore(1);
  //file_to_read >> time_read.seconds;

  // Store the Time structure read, time_read, into the end of time_array
  //time_array.push_back(time_read);


  // Read the volt values and store them in volt1, volt2 and volt3
  // Remeber that the default behaviour of a stream is to read until the space character is reached
  file_to_read >> v ;
  i += 1;
// Store each volt value into the appropriate array
v_array.push_back(v);
}

 cout << i;

return 0;
}
