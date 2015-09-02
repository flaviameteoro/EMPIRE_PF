#!/bin/bash
set -o nounset

tmp_file=tmp_instructions.sh

cat <<EOF > $tmp_file
#!/bin/bash
set -o nounset
set -o errexit
EOF

cat ../examples/lorenz96/instructions.txt >> $tmp_file

chmod u+x $tmp_file

cd ..
./tests/$tmp_file
echo ''
echo examples/lorenz96/instructions.txt returned with error code: $?
echo 

if [[ -f model_specific_backup.f90 ]]; then
    mv model_specific_backup.f90 model_specific.f90
fi

cd tests
rm $tmp_file
