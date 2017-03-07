#!/bin/bash

cp -RL R/snptree packagedir 
rm packagedir/snptree/src/*.o packagedir/snptree/src/makelinks packagedir/snptree/src/*.so
rm packagedir/snptree/README.md 
rm packagedir/snptree/src/tnt/gzstream.tgz
