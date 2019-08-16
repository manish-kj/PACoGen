The posit divison generator is based on Newton-Raphson method of division.


File description:

1. posit_div.v 			: Top-module which takes N (posit word size) and es (posit exponent size). It also includes all the sub-modules.

Below are the files for test-module for posit div with N=8, ES=1 (User can test for other options). It is an all exhaustive test for 8-bit operands (excluding Infinity diviplications to avoid comparision with julia pachage interupts for them).

2. posit_div_8bit_tb.v		: Test-bench module for 8-bit. 	
3. posit_div_8bit.sh		: A ModelSim bash script to invoke and run modelsim simulator to run the test-bench.
4. Pin1_8bit.txt		: Input-1 8-bit (Infinity diviplications are removed to avoid interupt from corresponding julia result comparisions)	 
5. Pin2_8bit.txt 		: Input-2 8-bit (Infinity diviplications are removed to avoid interupt from corresponding julia result comparisions)
6. Pout_8bit_ES1.txt		: Pre-stored posit diviplication results for comparison purpose. 
7. error_8bit.txt		: File will be generated during simulation which contains the difference of result produce by the Verilog module with pre-stored posit addition results. 
8. julia_posit8_div.sh	: This is a bash shell script for posit addition using julia posit package. It is currently using 8-bit inputs.  Julia posit package can be downloaded from https://github.com/interplanetary-robot/SigmoidNumbers


Below are the files for test-module for posit div with few other combinations: (N=32, ES=1) (N=32, ES=5), (N=16, ES=3).
9. posit_div_tb.v		: A general test-bench module. Edit it for required combination, currently Set for N=32, ES=5 and NR_Iter=2.
10. posit_div.sh		: A ModelSim bash script to invoke and run modelsim simulator to run the test-bench.
11. Pout_16bit_ES3.txt  	: (N=16, ES=3) Pre-stored posit diviplication results for comparison purpose. 
12. Pout_32bit_ES1.txt   	: (N=32, ES=1) Pre-stored posit diviplication results for comparison purpose. 
13. Pout_32bit_ES5.txt    	: (N=32, ES=5) Pre-stored posit diviplication results for comparison purpose. 
14. input32bit.txt		: Random input set.
