#!/bin/sh -x

# mkdir test

convert -density 150 bayes_direct_logconc.pdf -quality 90 bayes_direct_logconc.png

# shell script doesn't work, due to some technicalities
# see https://askubuntu.com/questions/152607/cant-make-a-file-executable


convert -density 600 bayes_direct.pdf                -quality 100 bayes_direct.png
convert -density 600 bayes_direct_par.pdf            -quality 100 bayes_direct_par.png
convert -density 600 bayes_direct_logconc.pdf        -quality 100 bayes_direct_logconc.png
convert -density 600 bayes_indirect.pdf              -quality 100 bayes_indirect.png
convert -density 600 bayes_indirect_par_even.pdf     -quality 100 bayes_indirect_par_even.png
convert -density 600 bayes_indirect_par_odd.pdf      -quality 100 bayes_indirect_par_odd.png
convert -density 600 bayes_indirect_logconc_even.pdf -quality 100 bayes_indirect_logconc_even.png
convert -density 600 bayes_indirect_logconc_odd.pdf  -quality 100 bayes_indirect_logconc_odd.png

