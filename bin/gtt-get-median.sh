#!/usr/bin/env bash

sort -n $1 | awk '
  BEGIN {
    counts = 0
    sum = 0
  }
  {
    values[counts++] = $1
    sum += $1
  }
  END {
    if ( NR == 2 ) {
      median = sum/2
    } else if ( (NR % 2) == 1 ) {
      median = values[ int(counts/2) ]
    } else {
      median = ( values[counts/2] + a[c/2-1] ) / 2
    }
    print median;
  }
'