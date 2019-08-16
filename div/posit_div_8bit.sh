vlib work

vlog "posit_div.v" 
vlog "posit_div_8bit_tb.v"

vsim -t ps work.posit_div_8bit_tb_v
view wave
#add wave *
run -all
