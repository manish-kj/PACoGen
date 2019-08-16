File description:

1. posit_mult.v 		: Top-module which takes N (posit word size) and es (posit exponent size). It also includes all the sub-modules.

Below are the files for test-module for posit mult with N=8, ES=4 (User can test for other options). It is an all exhaustive test for 8-bit operands (excluding Infinity multiplications to avoid comparision with julia pachage interupts for them).

2. posit_mult_8bit_tb.v		: Test-bench module. 	
3. posit_mult_8bit.sh		: A ModelSim bash script to invoke and run modelsim simulator to run the test-bench.
4. Pin1_8bit.txt		: Input-1 8-bit (Infinity multiplications are removed to avoid interupt from corresponding julia result comparisions)	 
5. Pin2_8bit.txt 		: Input-2 8-bit (Infinity multiplications are removed to avoid interupt from corresponding julia result comparisions)
6. Pout_8bit_ES4.txt		: Pre-stored posit multiplication results for comparison purpose. 

**. error_8bit.txt		: File will be generated during simulation which contains the difference of result produce by the Verilog module with pre-stored posit addition results. 


7. julia_posit8_mult.sh	: This is a bash shell script for posit addition using julia posit package. It is currently using 8-bit inputs.  Julia posit package can be downloaded from https://github.com/interplanetary-robot/SigmoidNumbers
