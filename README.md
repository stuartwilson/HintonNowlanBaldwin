# HintonNowlanBaldwin

Hinton and Nowlan's 1987 demonstration of the Baldwin Effect, and derivative code.

This code is based on the model described in the paper:

Hinton GE, Nowlan SJ (1981) How Learning Can Guide Evolution. Complex Systems. 1(3):495–502.

To run, clone morphologica into this directory:

git clone https://github.com/ABRG-Models/morphologica.git

Then: 

'''
mkdir build
cd build
cmake ..
make
cd ..
'''
Then to run:

./build/hinton config.json logs

(optionally append a seed for the random number generator)





