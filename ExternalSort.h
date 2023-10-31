/* ----------------------------------------------------------------------------
File - MyExternalSort.h
Date - 10/28/2023

External Sort Implementation

Developed by: Shaik, Rajmal Basha

Purpose:
This algorithm efficiently sorts a large dataset using external sort and the PolyPhase Merge technique.
----------------------------------------------------------------------------*/

#ifndef _MY_EXTERNAL_SORT_
#define _MY_EXTERNAL_SORT_
#include <vector>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <queue>
#include <assert.h>
#include <string>
#include <fstream>
#include <limits.h>
#include <math.h>
#include <chrono>
#include <cstdio>
#include <sys/stat.h>
#include <unistd.h>

#pragma region Global Variables

// total number of tapes you have,
// note that you should leave one tape for output during merging
const int TOTAL_TAPES = 16;

// the number of elements that can be loaded in to the main memory each time
const int BATCH_SIZE = 32;

// a list of kth order Fibonacci numbers
static std::vector<size_t> fibo_series;

#pragma endregion Global Variables

#pragma region File Handling Functions

/*------------------------------------------------------------------------------
  Function: CreateDirectory
    A Generic function for creating new directory

  Variables:
    dir_name - Name of the directory to be created
------------------------------------------------------------------------------*/
void createDirectory(std::string &dir_name)
{
  if (mkdir(dir_name.c_str(), 0777) == 0)
  {
    std::cout << "Successfully ctreated the directory: " << dir_name << std::endl;
  }
  else
  {
    std::cerr << "Failed to create the directory: " << dir_name << std::endl;
  }
}

/*------------------------------------------------------------------------------
  Function: deleteDirectory
    A Generic function for deleting directory

  Variables:
    dir_name - Name of the directory to be deleted

  Note: This function works only when the directory is empty
------------------------------------------------------------------------------*/
void deleteDirectory(std::string &dir_name)
{
  if (rmdir(dir_name.c_str()) == 0)
  {
    std::cout << "Successfully deleted the directory: " << dir_name << std::endl;
  }
  else
  {
    std::cerr << "Error removing directory: " << dir_name << std::endl;
  }
}

/*------------------------------------------------------------------------------
  Function: CreateFile
    A Generic function for creating new files.

  Variables:
    file_name - File name along with the path and file extension
    example - dir_name/filename.extension

------------------------------------------------------------------------------*/
void createFile(std::string &file_name)
{
  std::ofstream file(file_name);

  if (file.is_open())
  {
    file.close();
    std::cout << "Successfully created the file: " << file_name << std::endl;
  }
  else
  {
    std::cerr << "Error creating file: " << file_name << std::endl;
  }
}

/*------------------------------------------------------------------------------
  Function: deleteFile
    A Generic function for deleting a file

  Variables:
    dir_name - File name along with the path and file extension
------------------------------------------------------------------------------*/
void deleteFile(std::string &file_name)
{
  if (std::remove(file_name.c_str()) != 0)
  {
    std::cerr << "Error deleting file: " << file_name << std::endl;
  }
  else
  {
    std::cout << "File deleted: " << file_name << std::endl;
  }
}

#pragma endregion File Handling Functions

#pragma region Tape Handling Functions
/*------------------------------------------------------------------------------
  Function: createTapes
    Creates tapes to store sorted blocks of size as defined by BATCH_SIZE

  Returns:
    A vector of strings pointing the tapes newly created
------------------------------------------------------------------------------*/
std::vector<std::string> createTapes()
{
  // Holds the newly created tape names.
  std::vector<std::string> tapes;

  std::string dir_name = "Tapes";
  createDirectory(dir_name);

  for (int i = 1; i <= TOTAL_TAPES; ++i)
  {
    std::string file_name = dir_name + "/tape" + std::to_string(i) + ".txt";
    createFile(file_name);
    tapes.push_back(file_name);
  }
  return tapes;
}

/*------------------------------------------------------------------------------
  Function: deleteTapes
    Delete all the tapes which are created for storing the sorted input blocks
------------------------------------------------------------------------------*/
void deleteTapes()
{
  std::string dir_name = "Tapes";

  for (int i = 1; i <= TOTAL_TAPES; i++)
  {
    std::string file_name = dir_name + "/tape" + std::to_string(i) + ".txt";
    deleteFile(file_name);
  }

  deleteDirectory(dir_name);
}

#pragma endregion Tapes Handling Functions

#pragma region Data Handling Functions
/*------------------------------------------------------------------------------
   Function: writeSortedBlockToTape
    Writes the sorted block of data to a tape

   Variables:
    sorted_block: Sorted block of size BATCH_SIZE
    tape: The file name to which the sorted block must be written
------------------------------------------------------------------------------*/
void writeSortedBlockToTape(std::vector<int> sorted_block, std::string tape)
{
  std::ofstream of;

  // Open file in append mode
  of.open(tape, std::ios::app);

  // Check if the file opened successfully
  if (of)
  {
    // Write each integer to a new line in the file
    for (int i = 0; i < sorted_block.size(); ++i)
    {
      of << sorted_block[i] << std::endl;
    }

    // Close the file
    of.close();
  }
  else
  {
    std::cerr << "Error opening " << tape << " for writing." << std::endl;
  }
}

/*------------------------------------------------------------------------------
  Function: writeToOutputFile
    Writes the sorted elements in tape to an output file

  Variables:
    input_file_name - Tape containing the sorted elements
    output_file_name - Passed as an argument during command execution and where the final sorted elements are to be copied
    number_of_elements - Computed during initial run. Required to ignore the INT_MAX elements from being copied to final output file
    number_of_blocks - Number of blocks intially created. Required to ignore the INT_MAX elements from being copied to final output file
    reverse - if set true then the elements are in descending order else in ascending order
------------------------------------------------------------------------------*/
void writeToOutputFile(std::string input_file, std::string output_file, int number_of_elements, int number_of_blocks, bool reverse = false)
{
  std::ifstream input(input_file);
  std::ofstream output(output_file, std::ios::trunc);

  if (!input.is_open() || !output.is_open())
  {
    std::cerr << "Error opening files." << input_file << std::endl;
    return;
  }

  int lines_to_copy = number_of_elements;
  int line_count = 0;
  std::string line;

  if (reverse)
  {
    // Skips the INT_MAX entries
    int numberOfLinesToSkip = (number_of_blocks * BATCH_SIZE) - number_of_elements;
    while (numberOfLinesToSkip--)
    {
      getline(input, line);
    }
  }

  while (line_count < lines_to_copy && std::getline(input, line))
  {
    output << line << std::endl;
    line_count++;
  }

  input.close();
  output.close();
}

/*------------------------------------------------------------------------------
  Function: RemoveDataFromFile
    Delete the elements in the tape

  Variables:
    file_name - File name with complete path address
------------------------------------------------------------------------------*/
void removeDataFromFile(std::string &file_name)
{
  // Open the file in truncation mode to clear its contents.
  std::ofstream file(file_name, std::ios::trunc);
  std::ifstream input(file_name);
  if (file.is_open())
  {
    file.close();
  }
  else
  {
    std::cout << "Error opening the file." << file_name << std::endl;
  }
}

#pragma endregion Data Handling Functions

#pragma region Computing Fibonnaci Series
/*------------------------------------------------------------------------------
  Function: ComputeKthFibonacci
    Compute K-th order Fibonacci number for computing number of blocks should store in each tape.

  Note: you may call this function as many times as needed.
------------------------------------------------------------------------------*/
void computeKthFibonacci(void)
{
  /* formula:
    F(k)(0, 1, ..., k-1) = 0, 0, ..., 1
    F(k)(n) = F(k)(n − 1) + F(k)(n − 2) + ··· + F(k)(n − k), n >= k

    when n == k:
      F(k)(k) = F(k)(k − 1) + F(k)(k − 2) + ··· + F(k)(k − k)
              = 1           + 0           + ... + 0 = 1
    when n > k:
      F(k)(n)     = F(k)(n − 1) + F(k)(n − 2) + ··· + F(k)(n − k)
      F(k)(n)     = F(k)(n − 1) + [F(k)(n − 2) + ··· + F(k)(n − k) + F(k)(n − k-1)] - F(k)(n − k-1)
      F(k)(n − 1) = F(k)(n − 2) + F(k)(n − 3) + ··· + F(k)(n − k) + F(k)(n − k-1)
      F(k)(n)     = 2*F(k)(n − 1) - F(k)(n − k-1)
  */
  int k = TOTAL_TAPES - 1;
  if (fibo_series.empty()) // initial fibo_series
  {
    for (int i = 0; i < k - 1; i++)
      fibo_series.push_back(0);
    fibo_series.push_back(1);
    fibo_series.push_back(1); // this is fibo_series[k]
  }

  for (int i = 0; i < 100; i++) // compute next 100 items for usage
  {
    int n = fibo_series.size();
    size_t new_item = 2 * fibo_series[n - 1] - fibo_series[n - k - 1];
    fibo_series.push_back(new_item);
  }
}

/*------------------------------------------------------------------------------
  Function: GetBlocksAt
    get the number of blocks you should write to at current pass.
  Vairables:
    pass      - current pass of PolyPhaseMergePass.
    blocks_at - where you want to store the block counts of eath tape to.

  Note:
  For the t-th input tape at pass n, it should hold blocks_at[t] blocks:
    blocks_at[t] = F(k)(n+k-2) + F(k)(n+k-3) + ... + F(k)(n+t-2), t = 1, 2, ..., k
------------------------------------------------------------------------------*/
void getBlocksAt(int pass, std::vector<int> &blocks_at)
{
  if (pass < 1)
    return;
  int k = TOTAL_TAPES - 1;
  blocks_at = std::vector<int>(k, 0);
  for (int t = 0; t < k; t++)
  {
    for (int tmp_i = pass + t - 1; tmp_i <= pass + k - 2; tmp_i++)
    {
      if (tmp_i >= fibo_series.size())
        computeKthFibonacci(); // compute more fibo_series
      blocks_at[t] += fibo_series[tmp_i];
    }
  }
}
#pragma endregion Computing Fibonacci Series

#pragma region Insertion Sort
/*------------------------------------------------------------------------------
  The insertion sort algorithm.

  Vairables:
   a              - the input array
   left and right - the left and end indexes of the range of the elements to be sorted, inclusive
   reverse        - if set true, sort in descending order. Default: false
------------------------------------------------------------------------------*/
template <typename Comparable>
void insertionSort(std::vector<Comparable> &a, int left, int right, bool reverse = false)
{
  // CODE BEGINS
  for (int i = left + 1; i <= right; ++i)
  {
    Comparable temp = a[i];
    int j = i;
    if (!reverse) // Ascending order
    {
      while (j > left && temp < a[j - 1])
      {
        a[j] = a[j - 1];
        j--;
      }
    }

    else // Descending order
    {
      while (j > left && temp > a[j - 1])
      {
        a[j] = a[j - 1];
        j--;
      }
    }
    a[j] = temp;
  }
}
#pragma endregion Insertion Sort

#pragma region Verify Sorted Result
/*------------------------------------------------------------------------------
  Function: IsSorted
    Checks if the content in input_file is sorted

  Vairables:
    input_file - the file name that you want to check if it's sorted
    reverse    -  if set true, sort in descending order; otherwise in ascending order

  Note: you can use this to check if your result is correctly sorted.
------------------------------------------------------------------------------*/
bool IsSorted(const std::string &input_file, bool reverse = false)
{

  bool sorted = true;
  std::ifstream in(input_file.c_str());
  if (!in.is_open())
  {
    std::cout << input_file << " doesn't exist!\n";
    return false;
  }
  else
  {
    std::string buffer;
    int prev = INT_MIN, curr;
    if (reverse)
      prev = INT_MAX;
    while (!in.eof())
    {
      in >> curr;
      if ((curr < prev && !reverse) || (curr > prev && reverse))
      {
        sorted = false;
        std::cout << "Out of order: " << prev << ", " << curr << std::endl;
        break;
      }
      prev = curr;
    }
  }
  in.close();
  return sorted;
}

#pragma endregion Verify Sorted Result

#pragma region Load Sorted Blocks To Tapes
/*------------------------------------------------------------------------------
   Function: loadSortInputBlocks
     (1) Loads input data as blocks;
     (2) sort each block internally;
     (3) distribute blocks to the tapes for PolyPhaseMerge;
     (4) return how many passes you need to perform the PolyPhaseMerge.

   Vairables:
     input_file  - the file name that contains the inputs
     ext_arrays  - names of the files that serve as the external tapes
     reverse     - if set true, sort in descending order; otherwise in ascending order
------------------------------------------------------------------------------*/
std::vector<int> loadSortInputBlocks(
    const std::string &input_file,
    const std::vector<std::string> &ext_arrays,
    bool reverse = false)
{
  int curr_pass = 1;
  int total_number_of_elements = 0;
  int total_number_of_blocks = 0;
  int tape_number = 0; // 0 based indexing is used for tapes

  std::ifstream input(input_file);
  std::string line;

  // int maxElement = 0;

  // Holds the number of blocks currently placed in each tape before reading further inputs.
  std::vector<int> current_block_count_in_each_tape(TOTAL_TAPES - 1, 0);

  // Compute the number of blocks that should be present in each tape in that particular pass
  std::vector<int> blocks_count_in_current_pass(TOTAL_TAPES - 1, 0);

  // Holds the difference between newly computed blocks in each tape and current number of blocks in each tape.
  std::vector<int> blocks_to_be_inserted_each_tape(TOTAL_TAPES - 1, 0);

  // Read all inputs as blocks of BATCH_SIZE, sort and place the sorted block in the corresponding tape which is computed by kth order fibonacci series
  do
  {
    // Get blocks for the new pass
    getBlocksAt(curr_pass, blocks_count_in_current_pass);

    int total_blocks_to_be_inserted = 0;

    // Compute total blocks that needs to be inserted in the current pass as well as number of blocks that needs to be inserted in a particular tape.
    for (int i = 0; i < TOTAL_TAPES - 1; i++) // O(15)  Constant time
    {
      blocks_to_be_inserted_each_tape[i] = blocks_count_in_current_pass[i] - current_block_count_in_each_tape[i];
      total_blocks_to_be_inserted += blocks_count_in_current_pass[i] - current_block_count_in_each_tape[i];
      current_block_count_in_each_tape[i] = blocks_count_in_current_pass[i];
    }

    // Total number of blocks is equivalent to summation of blocks to be inserted at each pass
    total_number_of_blocks += total_blocks_to_be_inserted;

    // Insert the new blocks that are computed above
    do
    {
      // Create block of size as defined by BATCH_SIZE
      std::vector<int> block;

      for (int i = 0; i < BATCH_SIZE; i++)
      {
        // Read each line and place it in block for sorting
        if (std::getline(input, line))
        {
          total_number_of_elements++;
          block.push_back(std::stoi(line));
          // maxElement = std::max(maxElement, stoi(line));
        }

        else
        {
          // Insert INT_MAX into block as special character when the input file is exhausted
          block.push_back(INT_MAX);
        }
      }

      // Sort the new block
      insertionSort(block, 0, BATCH_SIZE - 1, reverse);

      int tape_counter = 0;

      // Write the block to the corresponding tape
      while (true && tape_counter < 15)
      {
        if (blocks_to_be_inserted_each_tape[(tape_number % (TOTAL_TAPES - 1))])
        {
          writeSortedBlockToTape(block, ext_arrays[(tape_number % (TOTAL_TAPES - 1))]);

          // Reduce the number of blocks that needs to be inserted by 1
          total_blocks_to_be_inserted--;

          // Reduce the number of block to be inserted in the particular tape by 1.
          blocks_to_be_inserted_each_tape[(tape_number % (TOTAL_TAPES - 1))]--;
          break;
        }
        else
        {
          tape_counter++;
          tape_number++;
        }
      }

      tape_number++;

    } while (total_blocks_to_be_inserted > 0);

    curr_pass++;

  } while (!input.eof());

  return {curr_pass - 1, total_number_of_elements, total_number_of_blocks};
}

#pragma endregion Load Sorted Blocks To Tapes

#pragma region Sort & Merge Blocks
/*------------------------------------------------------------------------------
   Function: sortIthPass
     (1) Merge one block from each tape into one block, and write the merged block into an external tape.
     (2) Repeat step 1 until all the blocks with second least block count is completely merged
     (2) Clean up intermediate file.
     (3) Reopen the output file in input mode again for next merge phase

   Vairables:
     in[]   - An array of input file stream objects.
     output_index  - Output file index (0 based indexing is used for tapes)
     blocks_at_ith_merge - A reference to a vector of integers representing blocks at the i-th merge stage
     ext_arrays  - names of the files that serve as the external tapes
     number_of_elements_in_block - A reference to a vector of integers storing the number of elements in each block
     reverse     - if set true, sort in descending order; otherwise in ascending order
------------------------------------------------------------------------------*/
void sortIthPass(std::ifstream in[], int output_index, std::vector<int> &blocks_at_ith_merge, std::vector<std::string> &ext_arrays, std::vector<int> &number_of_elements_in_block, bool reverse)
{
  // Remaining elements in that particular block that are yet to be merged to final output
  std::vector<int> remaining_elements_in_block = number_of_elements_in_block;

  int firstmin_block_size = INT_MAX;
  int secondmin_block_size = INT_MAX;

  // Computing the second minimum number of blocks
  // The first minimum will be always 0 and it is the output tape
  // We need to merge all the blocks in the tape which has the second least number of blocks
  for (int block_size : blocks_at_ith_merge)
  {
    if (block_size < firstmin_block_size)
    {
      secondmin_block_size = firstmin_block_size;
      firstmin_block_size = block_size;
    }
    else if (block_size < secondmin_block_size && block_size > firstmin_block_size)
    {
      secondmin_block_size = block_size;
    }
  }

  int temp_secondmin_block_size = secondmin_block_size;

  // Open the output file in append mode
  std::ofstream of;
  of.open(ext_arrays[output_index % TOTAL_TAPES], std::ios::app);

  // Ascending Order
  if (!reverse)
  {
    while (temp_secondmin_block_size--)
    {
      remaining_elements_in_block = number_of_elements_in_block;

      std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<std::pair<int, int>>> pq;

      // Initially we load the first element in all the tapes to the priority queue
      for (int k = 0; k < TOTAL_TAPES; k++)
      {
        if ((output_index % TOTAL_TAPES) != k)
        {
          std::string line;
          if (getline(in[k], line))
          {
            pq.push({std::stoi(line), k});
            remaining_elements_in_block[k]--;
          }
        }
      }

      // We pop the minimum element from priority queue and insert new element from that particular tape into priority queue
      while (!pq.empty())
      {

        auto it = pq.top();
        pq.pop();
        int value = it.first;
        int tape_number = it.second;
        of << value << std::endl;

        if (remaining_elements_in_block[tape_number])
        {
          std::string line;
          if (getline(in[tape_number], line))
          {
            pq.push({std::stoi(line), tape_number});
            remaining_elements_in_block[tape_number]--;
          }
        }
      }
    }

    // Cleaning up all the tapes which contains second minimum number of blocks
    for (int i = 0; i < TOTAL_TAPES; i++)
    {
      if (blocks_at_ith_merge[i] == secondmin_block_size)
      {
        removeDataFromFile(ext_arrays[i]);
      }
    }

    // Updating the number of blocks in each tape after the merging is completed
    for (int i = 0; i < TOTAL_TAPES; i++)
    {
      blocks_at_ith_merge[i] = abs(blocks_at_ith_merge[i] - secondmin_block_size);
    }

    // Close output file
    of.close();

    // Reopen the output file in input mode for merging in next pass
    in[output_index % TOTAL_TAPES].close();
    in[output_index % TOTAL_TAPES].open(ext_arrays[output_index % TOTAL_TAPES].c_str());
  }

  // Descending Order
  else
  {

    while (temp_secondmin_block_size--)
    {

      remaining_elements_in_block = number_of_elements_in_block;

      std::priority_queue<std::pair<int, int>> pq;

      // Initially we load the first element in all the tapes to the priority queue
      for (int k = 0; k < TOTAL_TAPES; k++)
      {
        if ((output_index % TOTAL_TAPES) != k)
        {
          std::string line;
          if (getline(in[k], line))
          {
            pq.push({std::stoi(line), k});
            remaining_elements_in_block[k]--;
          }
        }
      }

      // We pop the highest value element from priority queue and insert new element from that particular tape into priority queue
      while (!pq.empty())
      {

        auto it = pq.top();
        pq.pop();
        int value = it.first;
        int tape_number = it.second;
        of << value << std::endl;

        if (remaining_elements_in_block[tape_number])
        {
          std::string line;
          if (getline(in[tape_number], line))
          {
            pq.push({std::stoi(line), tape_number});
            remaining_elements_in_block[tape_number]--;
          }
        }
      }
    }

    // Cleaning up all the tapes which contains second minimum number of blocks
    for (int i = 0; i < TOTAL_TAPES; i++)
    {
      if (blocks_at_ith_merge[i] == secondmin_block_size)
      {
        removeDataFromFile(ext_arrays[i]);
      }
    }

    // Updating the number of blocks in each tape after the merging is completed
    for (int i = 0; i < TOTAL_TAPES; i++)
    {
      blocks_at_ith_merge[i] = abs(blocks_at_ith_merge[i] - secondmin_block_size);
    }

    // Close output file
    of.close();

    // Reopen the output file in input mode for merging in next pass
    in[output_index % TOTAL_TAPES].close();
    in[output_index % TOTAL_TAPES].open(ext_arrays[output_index % TOTAL_TAPES].c_str());
  }
}
#pragma endregion Sort &Merge Blocks

#pragma region PolyPhaseMerge
/*------------------------------------------------------------------------------
   Function: polyPhaseMerge
     (1) load data from k external tapes
     (2) repeat (3) until you have merged all blocks in any one of the tapes:
     (3) merge one block from each tape into one block, and write the merged block into an external tape.
     (4) clean up intermediate files.

   Vairables:
     cnt_pass    - total number of passes you need to perform PolyPhaseMerges,
                   this is computed from LoadSortInputBlocks
     ext_arrays  - the names of the list of files that serve as the external arrays
     reverse     - if set true, sort in descending order; otherwise in ascending order
------------------------------------------------------------------------------*/
int polyPhaseMerge(const int cnt_pass, std::vector<std::string> &ext_arrays, bool reverse)
{
  // Initial output index is 15 which is Tape16 (0 based indexing used)
  int output_index = TOTAL_TAPES - 1;

  // In initial case the size of each block is BATCH_SIZE except for output tape which has no elements and block size is 0.
  // After the merging, number of elements in each block of output tape in that pass is (TOTAL_TAPES-1)*(BATCH SIZE) (This is only for initial case)
  // In first merge the number of elements in each block of output tape is 15 * 32 = 480 elements.
  std::vector<int> number_of_elements_in_block_after_merging(TOTAL_TAPES - 1, BATCH_SIZE);
  number_of_elements_in_block_after_merging.push_back(0);

  int curr_pass = cnt_pass;

  std::ifstream in[TOTAL_TAPES];

  // Open all the tapes in input mode for reading
  for (int k = 0; k < TOTAL_TAPES; k++)
  {
    in[k].open(ext_arrays[k].c_str());
  }

  // To store blocks information
  std::vector<int> blocks_in_each_tape_at_ith_merge(TOTAL_TAPES - 1, 0);

  // Get blocks size at that particular pass
  getBlocksAt(cnt_pass, blocks_in_each_tape_at_ith_merge);
  blocks_in_each_tape_at_ith_merge.push_back(0);

  // Iteratively Merge blocks by computing the number of blocks that needs to merged at that particular pass
  while (curr_pass)
  {
    sortIthPass(in, output_index, blocks_in_each_tape_at_ith_merge, ext_arrays, number_of_elements_in_block_after_merging, reverse);
    int elements_after_merging = 0;

    // Compute the number of elements in each block of the output tape of that particular pass.
    for (int i = 0; i < TOTAL_TAPES; i++)
    {
      elements_after_merging += number_of_elements_in_block_after_merging[i];
    }

    // Update the block size after merging in the output tape of that particular pass with the value computed above
    number_of_elements_in_block_after_merging[(output_index % TOTAL_TAPES)] = elements_after_merging;

    output_index--;

    // Wrapping the output index, because output index can't be negative, and we are going in circular fashion.
    output_index = (output_index + TOTAL_TAPES) % TOTAL_TAPES;

    // Since one tape is completely read and there are no more further elements in that tape we set the block size = 0 of that particular tape.
    number_of_elements_in_block_after_merging[(output_index % TOTAL_TAPES)] = 0;

    curr_pass--;
  }

  // Close all the tapes if they are open
  for (int k = 0; k < TOTAL_TAPES; k++)
  {
    if (in[k].is_open())
      in[k].close();
  }

  // Wrapping the output index, because output index can't be negative, and we are going in circular fashion.
  return (output_index + TOTAL_TAPES + 1) % TOTAL_TAPES;
}

#pragma endregion PolyPhaseMerge

#pragma region External Sort Driver Function
/*------------------------------------------------------------------------------
 The driver external sort function function
   input_file   - the file name that contains the inputs
   output_file  - the file name that contains the final sorted elements
   reverse      - if set true, sort in descending order; otherwise in ascending order
------------------------------------------------------------------------------*/
void ExternalSort(const std::string &input_file, const std::string &output_file, bool reverse = false)
{
  std::vector<std::string> ext_arrays = createTapes();

  std::vector<int> res = loadSortInputBlocks(input_file, ext_arrays, reverse);

  // Number of pass required to load input
  int cnt_pass = res[0];
  // Number of elements in the input file
  int numberofElements = res[1];
  // Number of blocks that are placed in the tapes
  int numberofBlocks = res[2];

  int final_output_index = polyPhaseMerge(cnt_pass, ext_arrays, reverse);

  writeToOutputFile(ext_arrays[final_output_index], output_file, numberofElements, numberofBlocks, reverse);

  deleteTapes();
}

#pragma endregion External Sort Driver Function

#endif
