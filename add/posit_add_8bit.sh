vlib work

#All the verilog modules
vlog "posit_add_8bit_tb.v" 
vlog "posit_add.v"

vsim -t ps work.posit_add_8bit_tb_v
view wave
#add wave *
run -all
