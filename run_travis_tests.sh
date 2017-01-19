#!/usr/bin/env bash

#nosetests -v tests/test_style.py &> test_style.out &

# run nosetests
if [ $TEST_STYLE ] ; then
    nosetests -v tests/test_style.py ;
fi
