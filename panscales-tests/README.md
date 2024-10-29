# panscales-tests

This is a repository of analyses that run with the PanScales showers.

## How to

1. Download PanScales 

2. Make a note of the location of the .../panscales-main/2020-eeshower directory

3. Build with the following commands (putting in the correct location for the
   PanScales directory):

        cmake -DPARENT_DIR=/Users/albasotoontoso/Work/panscales-main/2020-eeshower .
        make -j

    Note that one can also simply run the scripts/build.py script e.g. with

        ../panscales-main/2020-eeshower/scripts/build.py

4. To run with PanScales + Pythia, one needs to compile with ``WITH_PYTHIA", e.g.
   with

        cmake [...] -DWITH_PYTHIA=ON -DPYTHIA_DIR=../panscales-main/2020-eeshower/pythia-interface/pythia

   or wherever your Pythia libraries are located. One can also just add these options to
   the scripts/build.py call with ``--cmake-options='-DWITH_PYTHIA=ON [...]'"


