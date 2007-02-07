#!/bin/tcsh

set num_args = ${#argv}

set myLang = "$LANG"
set tempLang = C
setenv LANG $tempLang

if ( $num_args == 3 ) then
    set input  = $argv[1]
    set bkmrk  = $argv[2]
    set output = $argv[3]

    ./bookmarkPdf.pl $input $bkmrk $output
else 
    echo "Usage: addBookmarks.csh <in.pdf> <bkmrk.txt> <out.pdf>"
    echo " "
    echo "Input PDF <in.pdf> and bookmarks text file <bkmrk.txt> must already exist,"
    echo "and results are output to bookmarked file <out.pdf>"
    echo "If <out.pdf> = <in.pdf>, the input file WILL be overwritten."
endif

setenv LANG $myLang

