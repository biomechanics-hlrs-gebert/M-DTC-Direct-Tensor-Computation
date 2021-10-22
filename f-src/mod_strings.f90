!==============================================================================
Module strings
  
  implicit none

contains

  !----------------------------------------------------------------------------
  Function char_to_str(char_arr) result(str)

    character, dimension(:) , intent(in) :: char_arr

    character(len=size(char_arr))        :: str

    Integer :: ii
    
    str = ''

    Do ii = 1, size(char_arr)
       str(ii:ii) = char_arr(ii)
    End Do

  End Function char_to_str

  !----------------------------------------------------------------------------
  Function char_to_linestr(char_arr) result(str)

    character, dimension(:) , intent(in) :: char_arr

    character(len=size(char_arr))        :: str

    Integer :: ii
    
    str = ''

    Do ii = 1, size(char_arr)
       if (char_arr(ii) == Char(10)) exit
       str(ii:ii) = char_arr(ii)
    End Do

  End Function char_to_linestr

  !----------------------------------------------------------------------------
  Function str_to_char(str) result(char_arr)

    character(len=*)  , intent(in) :: str

    character, dimension(len(str)) :: char_arr

    Integer :: ii
    
    char_arr = ''

    Do ii = 1, size(char_arr)
       char_arr(ii) = str(ii:ii)
    End Do

  End Function str_to_char

end Module strings
