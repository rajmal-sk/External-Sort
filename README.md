# External Sorting Algorithm

1. Eight instances are included.
   input*[1-8].txt: the input file;
   
   output*[1-8]\_[1-2].txt: the correct output file,
   
   the \_1 suffix indicates the result of ascending sorting,
   
   the \_2 suffix indicates the result of descending sorting;

2. Compile MainTest to test data files: `g++ -std=c++11 Main.cpp -o externalsort`

3. You are encouraged to check time, memory:
   /usr/bin/time -v -o ${LOG_FILE} Lab2 ${input_[1-8].txt} ${result_[1-8]_1.txt} ${result_[1-8]_2.txt}
   
   for i in {1..8}; do /usr/bin/time -v -o "./GeneratedLogs/Log_${i}.txt" ./Lab2 "./Inputs/input*${i}.txt" "./GeneratedResults/result*${i}_1.txt" "./GeneratedResults/result_${i}\_2.txt"; done
