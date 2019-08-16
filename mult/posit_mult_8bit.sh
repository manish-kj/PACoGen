vlib work

vlog "posit_mult.v" 
vlog "posit_mult_8bit_tb.v"

vsim -t ps work.posit_mult_8bit_tb_v
view wave
#add wave *
run -all
