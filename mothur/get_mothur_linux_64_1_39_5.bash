#!/usr/bin/env bash

wget https://github.com/mothur/mothur/releases/download/v1.39.5/Mothur.linux_64.zip
unzip Mothur.linux_64.zip
rm -R __MACOSX/ Mothur.linux_64.zip #tidy up folder
mv mothur code/



