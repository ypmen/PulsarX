# X_PYTHON_MODULE([PYTHON_MODULE], [ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
# -----------------------------------------------------
AC_DEFUN([X_PYTHON_MODULE],
[
  AC_MSG_CHECKING(Python module $1 installation)
  echo "import $1" | $PYTHON_BIN -
  if test $? -ne 0 ; then
    [$3]
  else
    [$2]
  fi
]) 
