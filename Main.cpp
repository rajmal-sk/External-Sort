#include "ExternalSort.h"
#include <stdexcept>
int main(int argc, char *argv[])
{
  if (argc < 4)
  {
    std::cout << "You must call your program using:\n";
    std::cout << "  ./YOUR_PROGRAM ${input} ${result1} ${result2}\n";
    std::cout << "    ${input} is the input file name.\n";
    std::cout << "    ${result1} and ${result2} are the file names where you want to save the result.\n";
    std::cout << "    ${result1} should save the non-descending sorted array.\n";
    std::cout << "    ${result2} should save the non-ascending sorted array.\n";
    return 0;
  }

  ExternalSort(argv[1], argv[2]);
  std::cout << "The input array is sorted in ascending order:" << (bool)IsSorted(argv[2]) << std::endl;
  ExternalSort(argv[1], argv[3], true);
  std::cout << "The input array is sorted in ascending  order:" << (bool)IsSorted(argv[3]) << std::endl;

  return 0;
}
