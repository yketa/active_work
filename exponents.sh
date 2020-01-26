#! /bin/bash

# This bash shell script enables one to use functions
# active_work.exponents.float_to_letters and
# active_work.exponents.letters_to_float from the shell.
#
# Information about litteral notation of floats used in this project can be
# found in active_work/exponents.py.

letters_to_float(){
	# Converts litteral expression to float expression.
	/usr/bin/python3 -c "from active_work.exponents import letters_to_float; print(letters_to_float('$1'))"
}
export -f letters_to_float

float_to_letters(){
	# Converts float expression to litteral expression.
	/usr/bin/python3 -c "from active_work.exponents import float_to_letters; print(float_to_letters($1))"
}
export -f float_to_letters
