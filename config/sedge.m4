# BEAR_LIB_SEDGE([ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
# -----------------------------------------------------
AC_DEFUN([BEAR_LIB_SEDGE],
[
  AC_PROVIDE([BEAR_LIB_SEDGE])
  AC_MSG_CHECKING([SEDGE installation])
  AC_LANG_PUSH(C++)
  AC_TRY_LINK([#include "se_cmd.h"], [Se_CMD_Key_Num;],
                  have_sedge=yes, have_sedge=no)
  AC_LANG_POP(C++)
  if test "$have_sedge" = "yes"; then
    AC_DEFINE([HAVE_SEDGE],[1],[Define if the SEDGE library is present])
    LIBS="-lsedge $LIBS"
    [$1]
  else
    [$2]
  fi

  AM_CONDITIONAL(HAVE_SEDGE,[test "$have_sedge" = "yes"])

])
