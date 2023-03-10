for i in $(seq 40 5 140); do

# create directories for each calculation
mkdir test_42atom_z$i

(python3 program.py <<EOF 
cc.cif
2
1
$i
1
3 
EOF
) | tail -n 77 > test_42atom_z$i/test_42atom_z$i.scf.in 

sed "s/test_42atom/test_42atom_z$i/g" test_42atom.sub >  ./test_42atom_z$i/test_42atom_z$i.sub


done