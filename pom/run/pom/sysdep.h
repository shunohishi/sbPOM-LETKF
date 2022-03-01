/*
      Creator: Takashi Kagimoto (kagimoto@jamstec.go.jp)
      Create:  2004/May/07
      Modify: Sergey Varlamov: limit __GNUC__ for F77 only,
              normal gnuc in 2009 is fortran90-enabled
*/
#if defined __GNUC__ && defined _FORTRAN77
  /* if using GNU FORTRAN77 */
# define NOT_HAVE_LEN_TRIM
# define NOT_HAVE_TRIM
# define NOT_HAVE_ADJUSTL
#endif
#if defined __sgi && defined _LANGUAGE_FORTRAN77
  /* if using SGI MIPS Pro FORTRAN77 */
# define NOT_HAVE_LEN_TRIM
# define NOT_HAVE_TRIM
# define NOT_HAVE_ADJUSTL
#endif
#if defined __sun && defined __SUNPRO_F77
  /* if using Sun Workshop Fortran 77 compiler */
# define NOT_HAVE_LEN_TRIM
# define NOT_HAVE_TRIM
# define NOT_HAVE_ADJUSTL
#endif
