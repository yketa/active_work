#! /bin/bash

# This bash shell script enables one to use functions
# active_work.exponents.float_to_letters and
# active_work.exponents.letters_to_float from the shell.

# Execute with `source exponents.sh` or `. exponents.sh` in order to correctly
# export functions to the current shell.

# Information about litteral notation of floats used in this project can be
# found in active_work/exponents.py.

PYTHON=${PYTHON-python} # Python executable
# NOTE: Using a weird syntax for function declarations in order to enforce this
#       choice of executable every time this script is executed.

eval "$(cat <<EOF
letters_to_float() {
  # Converts litteral expression to float expression.
  $PYTHON -c "from active_work.exponents import letters_to_float; print(letters_to_float('\$1'))"
}
EOF
)"
export -f letters_to_float

eval "$(cat <<EOF
float_to_letters() {
  # Converts float expression to litteral expression.
  $PYTHON -c "from active_work.exponents import float_to_letters; print(float_to_letters(\$1))"
}
EOF
)"
export -f float_to_letters
