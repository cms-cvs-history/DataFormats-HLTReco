#!/bin/tcsh

set num_args = ${#argv}

set myLang = "$LANG"
set tempLang = C
setenv LANG $tempLang

if ( $num_args == 3 ) then
    set input  = $argv[1]
    set bkmrk  = $argv[3]
    set output = $argv[2]

    $CMSSW_RELEASE_BASE/src/DataFormats/HLTReco/test/bookmarkPdf.pl $input $bkmrk $output
else if ( $num_args == 2 ) then 
    set input  = $argv[1]
    set base   = `echo $input | cut -d '.' -f 1`
    set bkmrk  = "$base-bookmark.txt"
    set output = $argv[2]

    $CMSSW_RELEASE_BASE/src/DataFormats/HLTReco/test/bookmarkPdf.pl $input $bkmrk $output

else if ( $num_args == 1 ) then 
    set input  = $argv[1]
    set base   = `echo $input | cut -d '.' -f 1`
    set bkmrk  = "$base-bookmark.txt"
    set output = "$base-bm.pdf"

    $CMSSW_RELEASE_BASE/src/DataFormats/HLTReco/test/bookmarkPdf.pl $input $bkmrk $output
else 
    echo "Usage: addBookmarks.csh <in.pdf> (<out.pdf>) (<bkmrk.txt>)"
    echo " "
    echo "Input PDF <in.pdf> and bookmarks text file <bkmrk.txt> must already exist,"
    echo "and results are output to bookmarked file <out.pdf>"
    echo "If <out.pdf> = <in.pdf>, the input file WILL be overwritten."
    echo "NOTE: <out.pdf> and <bkmrk.txt> are optional.  Defaults are [in]-bm.pdf and [in]-bookmark.txt, respectively"
endif

setenv LANG $myLang

