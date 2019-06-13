File description:
1. posit_add.v 			: Top-module which takes N (posit word size) and es (posit exponent size). It also contains all the required sub-module. 

Below are the files for test-module for posit adder with N=8, ES=4 (User can test for other options). It is an all exhaustive test for 8-bit operands.

2. posit_add_8bit_tb.v		: Test-bench module. 	
3. posit_add_8bit.sh		: A ModelSim bash script to invoke and run modelsim simulator to run the test-bench.
4. Pin1_8bit.txt		: Input-1 8-bit 	 
5. Pin2_8bit.txt 		: Input-2 8-bit
6. Pout_8bit_ES4.txt		: Pre-stored posit addition results for comparison purpose. 

**. error_8bit.txt		: File will be generated during simulation which contains the difference of result produce by the Verilog module with pre-stored posit addition results. 


7. julia_posit8_add.sh		: This is a bash shell script for posit addition using julia posit package. It is currently using 8-bit inputs. Julia posit package can be downloaded from https://github.com/interplanetary-robot/SigmoidNumbers
