#!/usr/bin/env bash

echo "Reading CSVs from $1"
echo " "
sleep 2

IFS='
'

for F in $1*.csv; do
    if [ `file -I "$F" | awk '{print $3;}'` = "charset=iso-8859-1" ]; then
      name=$(basename $F)
      # echo $name
      echo "Converting $name from iso-8859-1 to utf8"
      iconv -f iso-8859-1 -t utf-8 "$F" > "$2/$name"
    elif [ `file -I "$F" | awk '{print $3;}'` = "charset=unknown-8bit" ]; then
      name=$(basename $F)
      echo "Converting $name from unknown-8bit to utf8"
      iconv -f iso-8859-1 -t utf-8 "$F" > "$2/$name"
    else
      name=$(basename $F)
      echo "Copying $name, which is already utf8 encoded"
      cp -f "$F" "$2$name"
    fi
done
